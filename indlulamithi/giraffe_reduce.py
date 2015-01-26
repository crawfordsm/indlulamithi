# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Tools associated with data reduction and image handling
"""

import os
import numpy as np
from astropy import units as u
import ccdproc

from pyhrs import ccd_process


default_headers = ['FILENAME', 'OBJECT', 'IMAGETYP', 'TIME-OBS', 'RA', 'DEC', 'EXPTIME']

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

