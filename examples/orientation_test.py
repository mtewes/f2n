# The following two lines let this script find f2n, even if f2n has not yet been installed.
import sys, os
sys.path.insert(0, os.path.abspath('../'))

# Once f2n is installed, you can skip those and start from here.

import logging
logging.basicConfig(format='%(levelname)s: %(name)s(%(funcName)s): %(message)s', level=logging.DEBUG)

import numpy as np
import f2n
import matplotlib.pyplot as plt


array = np.array([[00,01,02,03,04,05], [10,11,12,13,14,15], [20,21,22,23,24,25], [30,31,32,33,34,35]])
# This image looks like (values give xy coords...)
# 05 15 25 35
# 04 14 24 34
# 03 13 23 33
# 02 12 22 32
# 01 11 21 31
# 00 10 20 30
print(array.shape)

#mask = np.zeros(array.shape, dtype=bool)
#mask[10:15,:]=True

 
# Without using SimpleFigure:
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
si = f2n.SkyImage(array)
f2n.draw_sky_image(ax, si)
cat = [{"x":1.0, "y":1.0, "text":"bottom left"}, {"x":1.0, "y":5.0, "text":"top left"}, {"x":3.0, "y":1.0, "text":"bottom right"}]
f2n.annotate(ax, cat, text="{row[text]}", ha="center")
plt.show()


"""
# Demo of some features using the convenient SimpleFigure class:
sf = f2n.SimpleFigure(array, scale=100)
sf.draw()
sf.annotate(cat, text="{row[text]}", color="green", ha="center")
sf.show()
#sf.save_to_file("test.png")
"""
