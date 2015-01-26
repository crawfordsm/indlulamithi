# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tools associated with data reduction and image handling
"""

import os
import numpy as np
from astropy import units as u
import ccdproc

from astropy.extern import six


default_headers = ['FILENAME', 'OBJECT', 'IMAGETYP', 'TIME-OBS', 'RA', 'DEC', 'EXPTIME']

def ccd_process(ccd, oscan=None, trim=None, error=False, masterbias=None,
                bad_pixel_mask=None, gain=None, rdnoise=None,
                oscan_median=True, oscan_model=None):
    """Perform basic processing on ccd data.
    
       The following steps can be included:
        * overscan correction
        * trimming of the image
        * create edeviation frame
        * gain correction
        * add a mask to the data
        * subtraction of master bias

       The task returns a processed `ccdproc.CCDData` object.
    Parameters
    ----------
    ccd: `ccdproc.CCDData`
        Frame to be reduced

    oscan: None, str, or, `~ccdproc.ccddata.CCDData`
        For no overscan correction, set to None.   Otherwise proivde a region
        of `ccd` from which the overscan is extracted, using the FITS
        conventions for index order and index start, or a
        slice from `ccd` that contains the overscan.

    trim: None or str
        For no trim correction, set to None.   Otherwise proivde a region
        of `ccd` from which the image should be trimmed, using the FITS
        conventions for index order and index start.

    error: boolean
        If True, create an uncertainty array for ccd

    masterbias: None, `~numpy.ndarray`,  or `~ccdproc.CCDData`
        A materbias frame to be subtracted from ccd.

    bad_pixel_mask: None or `~numpy.ndarray`
        A bad pixel mask for the data. The bad pixel mask should be in given
        such that bad pixels havea value of 1 and good pixels a value of 0.

    gain: None or `~astropy.Quantity`
        Gain value to multiple the image by to convert to electrons

    rdnoise: None or `~astropy.Quantity`
        Read noise for the observations.  The read noise should be in
        `~astropy.units.electron`


    oscan_median :  bool, optional
        If true, takes the median of each line.  Otherwise, uses the mean

    oscan_model :  `~astropy.modeling.Model`, optional
        Model to fit to the data.  If None, returns the values calculated
        by the median or the mean.

    Returns
    -------
    ccd: `ccdproc.CCDData`
        Reduded ccd

    """
    # make a copy of the object
    nccd = ccd.copy()

    # apply the overscan correction
    if isinstance(oscan, ccdproc.CCDData):
        nccd = ccdproc.subtract_overscan(nccd, overscan=oscan,
                                         median=oscan_median,
                                         model=oscan_model)
    elif isinstance(oscan, six.string_types):
        nccd = ccdproc.subtract_overscan(nccd, fits_section=oscan,
                                         median=oscan_median,
                                         model=oscan_model)
    elif oscan is None:
        pass
    else:
        raise TypeError('oscan is not None, a string, or CCDData object')

    # apply the trim correction
    if isinstance(trim, six.string_types):
        nccd = ccdproc.trim_image(nccd, fits_section=trim)
    elif trim is None:
        pass
    else:
        raise TypeError('trim is not None or a string')

    # create the error frame
    if error and gain is not None and rdnoise is not None:
        nccd = ccdproc.create_deviation(nccd, gain=gain, rdnoise=rdnoise)
    elif error and (gain is None or rdnoise is None):
        raise ValueError(
            'gain and rdnoise must be specified to create error frame')

    # apply the bad pixel mask
    if isinstance(bad_pixel_mask, np.ndarray):
        nccd.mask = bad_pixel_mask
    elif bad_pixel_mask is None:
        pass
    else:
        raise TypeError('bad_pixel_mask is not None or numpy.ndarray')

    # apply the gain correction
    if isinstance(gain, u.quantity.Quantity):
        nccd = ccdproc.gain_correct(nccd, gain)
    elif gain is None:
        pass
    else:
        raise TypeError('gain is not None or astropy.Quantity')

    # test subtracting the master bias
    if isinstance(masterbias, ccdproc.CCDData):
        nccd = nccd.subtract(masterbias)
    elif isinstance(masterbias, np.ndarray):
        nccd.data = nccd.data - masterbias
    elif masterbias is None:
        pass
    else:
        raise TypeError(
            'masterbias is not None, numpy.ndarray,  or a CCDData object')

    return nccd




def create_observing_dict(image_list, header_list=default_headers):
    """Create a dictionary with important information about the
       observations
    
    Parameters
    ----------
    image_list: list
        List of images for create the image list

    Returns
    -------
    obs_dict: dict
        Dictionary of information about each image
    """
    obs_dict={}
    for img in image_list:
        ccd = ccdproc.CCDData.read(img, unit=u.adu)
        hdr_list=[img]
        for hdr in header_list[1:]:
            hdr_list.append(ccd.header[hdr])
        obs_dict[os.path.basename(img)] = hdr_list
    return obs_dict

def make_camera_flat(image_list, outfile=None):
    """Combine together the imaging flats to create the camera
        flat

    Parameters
    ----------
    image_list: list
        List of images for create the image list

    outfile: None, str
        If a name of a file is provided, the camera flat will be 
        written out to that file
    
    Returns 
    -------
    cam_flat: `~ccdproc.CCDData`
        Median combined flat field image

    """
    cameraflat_list = []
    for img in image_list:
        ccd = giraffe_reduce(img)
        cameraflat_list.append(ccd)

    #for k in obs_dict: if obs_dict[k][0] == 'CAMERA': cameraflat_list.append(k)
    cb = ccdproc.Combiner(cameraflat_list)
    cam_flat = cb.median_combine(median_func=np.median)
    if outfile is not None:
        cam_flat.write(outfile, clobber=True)
  
def giraffe_reduce(image, camera_flat = None, outfile=None):
    """Basic ccd processing for Giraffe Data

    Parameters
    ----------
    image: str
        Name of image to be processed

    camera_flat: None, str
        If a name of a file is provided, the camera flat will be 
        written out to that file

    outfile: None, str
        If a name of a file is provided, the camera flat will be 
        written out to that file
    
    Returns 
    -------
    ccd: `~ccdproc.CCDData`
        Reduced ccd frame

    """
    ccd = ccdproc.CCDData.read(image, unit=u.electron)
    oscan = ccd.header['BIASSEC']
    trim = ccd.header['TRIMSEC']
    ccd = ccd_process(ccd, oscan=oscan, trim=trim)

    if camera_flat is not None:
        cam_flat = ccdproc.CCDData.read(camera_flat, unit=u.electron)
        ccd.data = ccd.data / cam_flat.data * cam_flat.data.mean()

    if outfile is not None:
        ccd.write(outfile, clobber=True)

    return ccd

