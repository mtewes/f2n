# The following two lines let this script find f2n, even if f2n has not yet been installed.
import sys, os
sys.path.insert(0, os.path.abspath('../'))

# Once f2n is installed, you can skip those and start from here.

import logging
logging.basicConfig(format='%(levelname)s: %(name)s(%(funcName)s): %(message)s', level=logging.INFO)

import f2n

# Read an image as a numpy array
image_array = f2n.read_fits("example.fits")
# And generate some tiny catalog (could also be an astropy table)
catalog = [
        {"x":112.0, "y":100.5},
        {"x":187.7, "y":238.0},
    ]

# Demo of some features using the convenient SimpleFigure class:
sf = f2n.SimpleFigure(image_array, z1="auto", z2="auto", scale=3)
sf.draw()
sf.draw_g_ellipses(catalog) # Further kwargs are passed to the matplotlib Ellipse

sf.show()
#sf.save_to_file("test.png")

