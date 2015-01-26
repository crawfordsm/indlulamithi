import sys
import os

import numpy as np
from scipy import ndimage as nd

from astropy.io import fits
from astropy import stats

import scipy.optimize as optimize


def minflattenimage(data, size=10):
    """Smooth the image and flatten it using the minimum value in the image

    Parameters
    ----------
    data: ndarray
       image to flatten

    size: int
       smoothing size for image

    Returns
    -------
    data: ndarray
        flattenned image
    """
    # flatten image in y direction
    m = np.median(data, axis=1)
    m = nd.minimum_filter(m, size=size)
    m.shape = (len(m), 1)
    data = data / m

    # flatten image in x direction
    m = np.median(data, axis=0)
    m = nd.minimum_filter(m, size=size)
    data = data / m
    return data


def calc_coef(data, xc, yc):
    """Given a position of an order,
       determine the equations that defines
       its position in the image
    """
    yc = int(yc)
    cutout = data.copy()
    obj, sci_num = nd.label(cutout)
    cutout[obj != obj[yc, xc]] = 0
    y, x = np.where(cutout > 0)
    coef = np.polyfit(x, y, 2)
    return cutout, coef


def make_orders(data, xc=680, limit=1.5, image_size=10, order_size=2, outfile=None):
    """Determine coefficients that describe all of the orders in the image

       Parameters
       ----------

       data: ndarray
           image array with orders in the image

       xc: int
           Column to extract orders from

       limit: float
           Limit for select orders in flattened data

       image_size: int
           Size for minimum filtering of images

       order_size: int
           Size for minimum filtering of orders

       Returns
       -------
       order_dict: dict
           Dictionary with the key representing the y-position of the 
           order at xc and containing a list of coefficients describing
           the shape of the order

    """
    # flatten the data
    data = minflattenimage(data, image_size)

    # create a rough image of just the location of the orders
    mask = (data < limit)
    data[mask] = 0

    # clean up the orders and caculate
    # starting position for each order
    n = nd.minimum_filter(data[:, xc], size=order_size)
    o, num = nd.label(n)
    pos = nd.center_of_mass(n, o, range(1, num))
    pos = np.array(pos)

    # determine the shape of the orders
    order_dict = {}
    for yc in pos:
        yc = yc[0]
        cutout, coef = calc_coef(data, xc, yc)
        order_dict[yc] = coef

    if outfile is not None:
        keys = sorted(order_dict.keys())
        fout = open(outfile, 'w')
        for i in keys:
            coef = order_dict[i]
            output = '%i ' % i
            output += ' '.join(['%e' % x for x in coef])
            if i > 0:
                fout.write(output + '\n')
        fout.close()

    return order_dict

