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
    
    Todo : use properties for z1 z2 !
    """
    
    def __init__(self, data, z1=None, z2=None):
        """
        
        """
               
        self.data = data
        self.extent = get_extent(self.data)
        self.z1 = z1
        self.z2 = z2
        logger.debug("SkyImage initialization with z1={} and z2={}".format(z1, z2))
        self.set_z(z1, z2)
        
        logger.info("Created {}".format(str(self)))

    @property
    def shape(self): # Just a shortcut
        """The shape (width, height) of the image"""
        return self.data.shape
    
    def __str__(self):
        return "SkyImage{}[{}:{}]".format(self.shape, self.z1, self.z2)
    

    def set_z(self, z1, z2):
        
        if z1 is None:
            self.z1 = np.min(self.data)
            logger.info("Set z1 to minimum value: {}".format(self.z1))
        else:
            self.z1 = z1
        if z2 is None:
            self.z2 = np.max(self.data)
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
    # "origin":"lower" as well as the tranpose() within the imshow arguments both combined give the right orientation
    imshow_kwargs = {"aspect":"equal", "origin":"lower", "interpolation":"none", "cmap":matplotlib.cm.get_cmap('Greys_r')}
    imshow_kwargs.update(kwargs)
    
    return ax.imshow(si.data.transpose(), vmin=si.z1, vmax=si.z2, extent=si.extent, **imshow_kwargs)



def draw_mask(ax, si, **kwargs):
    """Uses imshow to draw a binary mask to some axes
    
    Parameters
    si : SkyImage
        The mask, can also be a plain 2D numpy array
    
    """
        
    mask_cmap = matplotlib.colors.ListedColormap([(1.0, 1.0, 1.0, 0.0), (1.0, 0.0, 0.0, 0.6)])
    mask_bounds=[-1,0.5,1]
    mask_norm = matplotlib.colors.BoundaryNorm(mask_bounds, mask_cmap.N)

    imshow_kwargs = {"aspect":"equal", "origin":"lower", "interpolation":"none", "alpha":0.5}
    imshow_kwargs.update(kwargs)
    
    if isinstance(si, SkyImage):
        return ax.imshow(si.data.tranpose(), vmin=0, vmax=1, extent=si.extent, cmap=mask_cmap, norm=mask_norm, **imshow_kwargs)
    else: # We can also work with simple numpy arrays
        return ax.imshow(si.transpose(), vmin=0, vmax=1, extent=get_extent(si), cmap=mask_cmap, norm=mask_norm, **imshow_kwargs)
    

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
        
        self.img_array = img_array
        self.z1 = z1
        self.z2 = z2
        
        self.dpi = 72
        self.figsize = float(scale) * np.array(img_array.shape)/self.dpi
        
        self.fig = plt.figure(figsize=self.figsize)
        self.ax = self.fig.add_subplot(111)
        
        self.has_been_drawn = False
        
        

    def __str__(self):
        return "SimpleFigure"#({})".format(str(self.si))
    
    def draw(self, si=None):
        """Draw the image pixels on the axes.
        
        Usually you leave si to None, in which case a new SkyImage is built from what was passed to init.
        """
        if si is None:
            si = SkyImage(self.img_array, self.z1, self.z2)
        draw_sky_image(self.ax, si)
        self.has_been_drawn = True
    
    def check_drawn(self):
        if not self.has_been_drawn:
            self.draw()
            #logger.warning("The SimpleFigure has not been drawn, you probably want to call draw() before showing or saving it!")
        
    def draw_g_ellipses(self, cat, **kwargs):
        draw_g_ellipses(self.ax, cat, **kwargs)
    
    def annotate(self, cat, **kwargs):
        annotate(self.ax, cat, **kwargs)
    
    def show(self):
        """Update this once we settle on a minimum matplotlib version...
        
        """
        self.check_drawn()
        logger.info("Showing {}...".format(str(self)))
        plt.show()
        
        
    def save_to_file(self, filepath):
        self.check_drawn()
        logger.info("Saving {} to '{}'...".format(str(self), filepath))
        self.fig.savefig(filepath, bbox_inches='tight')
   




# Some utility functions

def get_extent(a):
    """Defines the extent with which to plot an array a (we use the numpy convention)
    
    """
    return (0, a.shape[0], 0, a.shape[1])

    
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
 
        
        
