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
#from __future__ import division, print_function
#from future_builtins import *

import os
import numpy as np
import astropy.io.fits

import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colors
import matplotlib.patches

import logging
logger = logging.getLogger(__name__)


class SkyImage(object):
    """Holds a pixel array
    
    """
    
    def __init__(self, data, z1=None, z2=None):
        """
        
        """
               
        self.data = data
        self.extent = (0, self.data.shape[0], 0, self.data.shape[1])
        self.set_z(z1, z2)
        
        logger.info("Created {}".format(str(self)))

    def __str__(self):
        return "SkyImage{}".format(self.data.shape)
    

    def set_z(self, z1, z2):
        
        if z1 is None:
            self.z1 = np.min(data)
            logger.info("Set z1 to minimum value: {}".format(self.z1))
        else:
            self.z1 = z1
        if z2 is None:
            self.z2 = np.max(data)
            logger.info("Set z2 to maximum value: {}".format(self.z2))
        else:
            self.z2 = z2
        
        autolist = ["Auto", "auto"]
        if (self.z1 in autolist) or (self.z2 in autolist):
            self.set_auto_z_scale()
            

    def set_auto_z_scale(self, full_sample_limit = 10000, nsig=5.0):
        """Automatic z-scale determination"""
   
        #if a.size > full_sample_limit :
        #    selectionindices = np.linspace(0, a.size-1, full_sample_limit).astype(np.int)
        #    statsel = calcarray[selectionindices]
        #else :
    
        stata = self.data[:]
        #stata.shape = (stata.size,) # Flattening the array

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
        
        logger.info("Set automatic zscale from {} to {}".format(self.z1, self.z2))



def draw_sky_image(ax, si, **kwargs):
    """Use imshow to draw a SkyImage to some axes
    
    """
    imshow_kwargs = {"aspect":"equal", "origin":"lower", "interpolation":"none", "cmap":matplotlib.cm.get_cmap('Greys_r')}
    imshow_kwargs.update(kwargs)
    
    return ax.imshow(si.data, vmin=si.z1, vmax=si.z2, extent=si.extent, **imshow_kwargs)



def draw_mask(ax, si, **kwargs):
    """Uses imshow to draw a binary mask to some axes
    
    """

    mask_cmap = matplotlib.colors.ListedColormap([(1.0, 1.0, 1.0, 0.0), (1.0, 0.0, 0.0, 0.6)])
    mask_bounds=[-1,0.5,1]
    mask_norm = matplotlib.colors.BoundaryNorm(mask_bounds, mask_cmap.N)

    imshow_kwargs = {"aspect":"equal", "origin":"lower", "interpolation":"none", "alpha":0.5}
    imshow_kwargs.update(kwargs)
  
    return ax.imshow(si.data, vmin=0, vmax=1, extent=si.extent, cmap=mask_cmap, norm=mask_norm, **imshow_kwargs)
    

def draw_ellipse(ax, x, y, a=5, b=None, angle=None, **kwargs):
    """Draws an ellipse patch on the axes
    
    """
    if b is None:
        b = a
    if angle is None:
        angle = 0
    
    ellipse_kwargs = {"fill":False, "edgecolor":"red", "linewidth":2, "clip_box":ax.bbox}
    ellipse_kwargs.update(kwargs)
    
    ellipse = matplotlib.patches.Ellipse((x, y), width=2*a, height=2*b, angle=angle, **ellipse_kwargs)
    
    return ax.add_artist(ellipse)

def draw_g_ellipse(ax, x, y, g1, g2, sigma, **kwargs):
    """Draws an ellipse defined by the "reduced shear" following GalSim nomenclature.
    
    Parameters
    ----------
    sigma : float
        sigma is defined as (a+b)/2
    
    """
    angle = 0.5*np.arctan2(g2, g1) * 180.0/np.pi
    g = np.hypot(g1, g2)
    a = (g + 1.0) * sigma
    b = (1.0 - g) * sigma
    return draw_ellipse(ax, x, y, a, b, angle, **kwargs)


def draw_g_ellipses(ax, cat, x="x", y="y", g1="g1", g2="g2", sigma="sigma", **kwargs):
    """Draws ellipses from a catalog (that is an astropy table or a list of dicts) of sources
    
    Parameters
    ----------
    cat : astropy table or list of dicts
        The "catalog", must be iterable like a list, with each element behaving like a dict.
    x : str
        The key of the field (or the column name) containing the x position
    y
        idem
    g1 : float
         The g1 component of the reduced shear defining the ellipse to draw
    g2
        idem
    sigma : float
        defined as (a+b)/2 of the ellipses to draw
    
    """
    for row in cat:
        # We skip silently any masked positions
        if getattr(row[x], "mask", False) or getattr(row[y], "mask", False):
            continue
        
        draw_g_ellipse(ax, row[x], row[y], row[g1], row[g2], row[sigma], **kwargs)
    

def annotate(ax, cat, x="x", y="y", text="Hello", **kwargs):
    """Annotates the positions (x, y) from a catalog
    
    """
    
    annotate_kwargs = {"horizontalalignment":"left", "verticalalignment":"top", "color":"red",
                       "xytext":(0, 0), "textcoords":'offset points'}
    annotate_kwargs.update(**kwargs)
    
    for row in cat:
        
        # We skip silently any masked positions
        if getattr(row[x], "mask", False) or getattr(row[y], "mask", False):
            continue
        
        rowtext = text.format(row=row)
        ax.annotate(rowtext,
            xy=(row[x], row[y]),
            **annotate_kwargs
            )


class SimpleFigure(object):
    """A simple plot made from scratch. If you have only one image, and want one matplotlib figure, this should do it.
    
    """
    
    def __init__(self, img_array, z1=None, z2=None, scale=1):
        """

        Parameters
        ----------
        img_array : numpy array
            A 2D numpy array containing the image

        scale : float
            A scaling for the display of the image


        """
        
        self.si = SkyImage(img_array, z1, z2)
        
        self.dpi = 72
        self.figsize = float(scale) * np.array(img_array.shape)/self.dpi
        
        self.fig = plt.figure(figsize=self.figsize)
        self.ax = self.fig.add_subplot(111)
        
        draw_sky_image(self.ax, self.si)
        

    def __str__(self):
        return "SimpleFigure({})".format(str(self.si))
            
    def draw_g_ellipses(self, cat, **kwargs):
        draw_g_ellipses(self.ax, cat, **kwargs)
    
    def annotate(self, cat, **kwargs):
        annotate(self.ax, cat, **kwargs)
    
    def show(self):
        """Update this once we settle on a minimum matplotlib version...
        
        """
        logger.info("Showing {}...".format(str(self)))
        plt.show()
        
    def save_to_file(self, filepath):
        logger.info("Saving {} to '{}'...".format(str(self), filepath))
        self.fig.savefig(filepath, bbox_inches='tight')
   

# 
# def plot(img_array, z1=None, z2=None, mask_array=None, ax=None, filepath=None, figsize=None, scale=1):
#     """Wraps some of the functionality into a single convenient function 
#     
#     Parameters
#     ----------
#     img_array : numpy array
#         A 2D numpy array containing the image
#     mask_array : numpy array
#         A 2D numpy array containing a boolean mask (1 means True means masked)
#     ax : Matplotlib Axes
#         An Axes object on which to plot the image.
#         If given, the keywords related to the figure creation are ignored
#     filepath : string
#         Path to a file in which to save the figure. If Neither ax nor fielpath is specified, the figure is shown
#         interactively.
#     scale : float
#         A scaling for the display of the image
#     
#     """
#     
#     if ax is None:
#         makefig = True
#     else:
#         makefig = False
#         
#     if makefig:
#         
#         dpi = 72.0 
#         figsize = float(scale) * np.array(img_array.shape)/dpi
#         
#         fig = plt.figure(figsize=figsize)
#         ax = fig.add_subplot(111)
#     
#     img_si = SkyImage(img_array)
#     img_si.set_auto_z_scale()
#     draw_sky_image(ax, img_si)
#     
#     if mask_array is not None:
#         mask_si = SkyImage(mask_array)
#         draw_mask(ax, mask_si)
#     
#     draw_ellipse(ax, 40, 60, a=5, b=10)
#     
#     sigma = 0.1
#     draw_g_ellipse(ax, -0.5, 0.0, -0.5, 0.0, sigma)
#     draw_g_ellipse(ax, 0.5, 0.0, 0.5, 0.0, sigma)
#     draw_g_ellipse(ax, 0.0, 0.5, 0.0, 0.5, sigma)
#     draw_g_ellipse(ax, 0.0, -0.5, 0.0, -0.5, sigma)
#     
#     
#            
#     if makefig:
#         if filepath is None:
#             plt.show()
#         else:
#             logger.info("Writing image to '{}'".format(filepath))
#             fig.savefig(filepath, bbox_inches='tight')
#         plt.close(fig)
#  




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
 
        
        
