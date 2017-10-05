#  
# Copyright (C) 2012-2020 Euclid Science Ground Segment      
#    
# This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General    
# Public License as published by the Free Software Foundation; either version 3.0 of the License, or (at your option)    
# any later version.    
#    
# This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied    
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more    
# details.    
#    
# You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to    
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA    
#

""" 
File: sky_image_plot.py
Created on: Sep 3, 2017
Author: Malte Tewes

This file is a standalone Euclid-agnostic module to visualize and plot sky images with matplotlib.
To visualize Euclid objects, use the wrappers in she_image_checkplot.py

""" 
from __future__ import division, print_function
from future_builtins import *

import os
import numpy as np
import astropy.io.fits

import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colors

import logging
logger = logging.getLogger(__name__)


class SkyImage(object):
    """Holds a pixel array
    
    """
    
    def __init__(self, data):
        """
        
        """
               
        self.data = data
        self.z1 = np.min(data)
        self.z2 = np.max(data)
        
        self.extent = (0, self.data.shape[0], 0, self.data.shape[1])
        
        logger.debug("Created {}".format(str(self)))

    def __str__(self):
        return "SkyImage{}".format(self.data.shape)
    

    def set_auto_z_scale(full_sample_limit = 10000, nsig=3.0):
        """Automatic z-scale determination"""
   
        #if a.size > full_sample_limit :
        #    selectionindices = np.linspace(0, a.size-1, full_sample_limit).astype(np.int)
        #    statsel = calcarray[selectionindices]
        #else :
    
        stata = self.data[:]
        stata.shape = stata.size # Flattening the array

        std = stdmad(stata)
        med = np.median(stata)
    
        nearskypixvals = stata[np.logical_and(stata > med - 2*std, stata < med + 2*std)]
        if len(nearskypixvals) > 0:
            skylevel = np.median(nearskypixvals)
        else:
            skylevel = med

        self.z1 = skylevel - nsig * std
        
        #sortedstata = np.sort(stata)
        #z2 = sortedstata[round(0.9995 * stata.size)]
        
        self.z2 = np.max(stata)
        
        logger.debug("Set zscale from {} to {}".format(self.z1, self.z2))



def draw_sky_image(si, ax, **kwargs):
    """Use imshow to draw a SkyImage to some axes
    
    """
    imshow_kwargs = {"aspect":"equal", "origin":"lower", "interpolation":"none", "cmap":matplotlib.cm.get_cmap('Greys_r')}
    imshow_kwargs.update(kwargs)
    
    return ax.imshow(si.data, vmin=si.z1, vmax=si.z2, extent=si.extent, **imshow_kwargs)



def draw_mask(si, ax, **kwargs):
    """Uses imshow to draw a binary mask to some axes
    
    """

    mask_cmap = matplotlib.colors.ListedColormap([(1.0, 1.0, 1.0, 0.0), (1.0, 0.0, 0.0, 0.6)])
    mask_bounds=[-1,0.5,1]
    mask_norm = matplotlib.colors.BoundaryNorm(mask_bounds, mask_cmap.N)

    imshow_kwargs = {"aspect":"equal", "origin":"lower", "interpolation":"none", "alpha":0.5}
    imshow_kwargs.update(kwargs)
  
    return ax.imshow(si.data, vmin=0, vmax=1, extent=si.extent, cmap=mask_cmap, norm=mask_norm, **imshow_kwargs)
    

def plot(img_array, mask_array=None, ax=None, filepath=None, figsize=None, scale=1):
    """Wraps some of the functionality into a single convenient function 
    
    """
    
    if ax is None:
        makefig = True
    else:
        makefig = False
        
    if makefig:
        
        dpi = 100.0 
        figsize = float(scale) * np.array(img_array.shape)/dpi
        
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
    
    img_si = SkyImage(img_array)
    draw_sky_image(img_si, ax)
    
    if mask_array is not None:
        mask_si = SkyImage(mask_array)
        draw_mask(mask_si, ax)
        
        
    if makefig:
        logger.info("Writing image to '{}'".format(filepath))
        fig.savefig(filepath, bbox_inches='tight')
        plt.close(fig)
 




# Some utility functions
    
def stdmad(a):
    """MAD rescaled to std of normally distributed data"""
    med = np.median(a)
    return 1.4826*np.median(np.abs(a - med))
    
    
def read_fits(filepath):
    """Reads simple 1-hdu FITS file into a numpy arrays
    
    A transposition is done, so that the indexes [x,y] of the numpy array follow the orientation of x and y in DS9
    and SExtractor.
    
    Parameters
    ----------
    filepath : str
        Filepath to read the array from
    
    """
    a = astropy.io.fits.getdata(filepath).transpose()
    logger.info("Read FITS images %s from file %s" % (a.shape, filepath))
    return a
    

def write_fits(a, filepath, clobber=True):
    """Writes a simple 2D numpy array into a FITS file
    
    As for read_fits, a transposition is applied to conserve the orientation.
    
    Parameters
    ----------
    a : array
    filepath : str
    clobber : bool
        Set this to False if an existing file should not be overwritten
    
    """
    if os.path.exists(filepath) and clobber:
        logger.info("File %s exists, I will overwrite it!" % (filepath))

    astropy.io.fits.writeto(filepath, a.transpose(), clobber=clobber)
    logger.info("Wrote %s array into %s" % (a.shape, filepath))
 
        
        
