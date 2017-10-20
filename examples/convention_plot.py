# The following two lines let this script find f2n, even if f2n has not yet been installed.
import sys, os
sys.path.insert(0, os.path.abspath('../'))

# Once f2n is installed, you can skip those and start from here.

import logging
logging.basicConfig(format='%(levelname)s: %(name)s(%(funcName)s): %(message)s', level=logging.DEBUG)

import f2n
import numpy as np

# Read an image as a numpy array
image_array = np.random.randn(5, 5)



# And generate some tiny catalog (could also be an astropy table)
catalog = [
        {"x":1.5, "y":2.5, "g1":-0.3, "g2":0.0, "sigma":0.5},
        {"x":3.5, "y":2.5, "g1":0.3, "g2":0.0, "sigma":0.5},
        {"x":2.5, "y":1.5, "g1":0.0, "g2":-0.3, "sigma":0.5},
        {"x":2.5, "y":3.5, "g1":0.0, "g2":0.3, "sigma":0.5},    
    ]

# Demo of some features using the convenient SimpleFigure class:
sf = f2n.SimpleFigure(image_array, z1=-3.0, z2=6.0, scale=100)
sf.draw()
sf.draw_g_ellipses(catalog, edgecolor="red") # Further kwargs are passed to the matplotlib Ellipse
sf.annotate(catalog, text="g1 = {row[g1]:+.1f}", color="white", xytext=(20, -10))
sf.annotate(catalog, text="g2 = {row[g2]:+.1f}", color="white", xytext=(20, -25))

sf.ax.set_xlabel("x, first index")
sf.ax.set_ylabel("y, second index")

for item in ([sf.ax.xaxis.label, sf.ax.yaxis.label] + sf.ax.get_xticklabels() + sf.ax.get_yticklabels()):
    item.set_fontsize(20)
sf.show()
#sf.save_to_file("test.png")

