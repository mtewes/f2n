# The following two lines let this script find f2n, even if f2n has not yet been installed.
import sys, os
sys.path.insert(0, os.path.abspath('../'))

# Once f2n is installed, you can skip those and start from here.

import logging
logging.basicConfig(format='%(levelname)s: %(name)s(%(funcName)s): %(message)s', level=logging.INFO)

import f2n

# Generate a tiny catalog (could also be an astropy table)
catalog = [
        {"x":112.0, "y":100.5, "g1":0.0, "g2":0.0, "sigma":30.0},
        {"x":187.7, "y":238.0, "g1":-0.4, "g2":0.0, "sigma":5.0},
        {"x":0.7, "y":20.0, "g1":0.1, "g2":0.2, "sigma":10.0},
    ]

# Demo of using SimpleFigure without any image:

sf = f2n.SimpleFigure()
sf.draw_g_ellipses(catalog, edgecolor="red") # Further kwargs are passed to the matplotlib Ellipse

sf.show()


