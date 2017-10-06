# The following two lines let this script find f2n, even if f2n has not yet been installed.
import sys, os
from galsim.wfirst import jitter_rms
sys.path.insert(0, os.path.abspath('../'))

# Once f2n is installed, you can skip those and start from here.

import logging
logging.basicConfig(format='%(levelname)s: %(name)s(%(funcName)s): %(message)s', level=logging.DEBUG)

import f2n
import numpy as np
import astropy.table
import galsim

# Preparing a catalog with very simple random galaxy parameters, on a (nstamps_x, nstamps_y) grid
stampsize = 100.0
nstamps_x = 5
nstamps_y = 5
catalog = []
for ix in range(nstamps_x):
    for iy in range(nstamps_y):
        row = {}
        row["ix"] = ix
        row["iy"] = iy
        row["tru_x"] = stampsize/2.0 + ix*stampsize
        row["tru_y"] = stampsize/2.0 + iy*stampsize
        row["tru_flux"] = np.random.uniform(1000.0, 5000.0)
        row["tru_sigma"] = np.random.uniform(2.0, 8.0)
        row["tru_g1"] = np.random.uniform(-0.3, 0.3)
        row["tru_g2"] = np.random.uniform(-0.3, 0.3)
        catalog.append(row)
catalog = astropy.table.Table(catalog, masked=True)
# We add masked empty columns for the HSM measurements:
for colname in ["hsm_flux", "hsm_x", "hsm_y", "hsm_g1", "hsm_g2", "hsm_sigma", "hsm_rho4"]:
    catalog.add_column(astropy.table.MaskedColumn(name=colname, mask=True, dtype=float, length=len(catalog)))
    
#print(catalog)
#exit()

# Drawing an image of stamps and running HSM

# Preparing a galsim image to contain the array of stamps
gal_image = galsim.ImageF(stampsize * nstamps_x , stampsize * nstamps_y)
gal_image.scale = 1.0

# And loop through the catalog
for row in catalog:
                
    # We will draw this galaxy in a postage stamp, but first we need the bounds of this stamp.
    ix = int(row["ix"])
    iy = int(row["iy"])
    bounds = galsim.BoundsI(ix*stampsize+1 , (ix+1)*stampsize, iy*stampsize+1 , (iy+1)*stampsize) # Default Galsim convention, index starts at 1.
    gal_stamp = gal_image[bounds]
    
    # A simple Gaussian profile
    gal = galsim.Gaussian(sigma=float(row["tru_sigma"]), flux=float(row["tru_flux"]))
    # We make this profile elliptical
    gal = gal.shear(g1=row["tru_g1"], g2=row["tru_g2"]) # This adds the ellipticity to the galaxy
    # Randomize its position with respect to the stamp center
    jitter_range = 10.0
    gal = gal.shift(np.random.uniform(-jitter_range,jitter_range), np.random.uniform(-jitter_range,jitter_range))
    
    gal.drawImage(gal_stamp)
    gal_stamp.addNoise(galsim.GaussianNoise(sigma=1.0))
    
    # We measure this stamp with HSM
    try:
        res = galsim.hsm.FindAdaptiveMom(gal_stamp, guess_sig=15.0)
    except:
        continue
    row["hsm_flux"] = res.moments_amp
    row["hsm_g1"] = res.observed_shape.g1
    row["hsm_g2"] = res.observed_shape.g2
    row["hsm_x"] = res.moments_centroid.x -0.5
    row["hsm_y"] = res.moments_centroid.y -0.5
    row["hsm_sigma"] = res.moments_sigma
    row["hsm_rho4"] = res.moments_rho4
   
#gal_image.write("gal_image.fits")
#print(catalog)

# And we visualize the grid and the measurements:
sf = f2n.SimpleFigure(gal_image.array, z1=-3.0, z2=10.0, scale=2) # There is no "tranpose" here, as we've defined the stamps with our own convention, and consistently drawn and measured them with galsim.
sf.draw_g_ellipses(catalog, x="hsm_x", y="hsm_y", g1="hsm_g1", g2="hsm_g2", sigma="hsm_sigma", edgecolor="red")
sf.annotate(catalog, x="hsm_x", y="hsm_y", text="g1 = {row[hsm_g1]:.2f}", color="white", xytext=(20, -10))
sf.annotate(catalog, x="hsm_x", y="hsm_y", text="g2 = {row[hsm_g2]:.2f}", color="white", xytext=(20, -25))


sf.show()



