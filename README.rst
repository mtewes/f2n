f2n v2: FITS to net, in python 3
================================

The purpose of the tiny ``f2n`` module is to turn FITS images (or numpy arrays) into (potentially annotated) web-friendly images (typically png, but also svg or pdf if you want vector-graphics) in a modular and flexible way.

What a strange idea to use a software like SExtractor without f2n! (Me) 

The present code is an attempt, under construction, to replace some of the (ridiculously useful) functionality of my previous `f2n.py module <https://obswww.unige.ch/~tewes/f2n_dot_py/>`_, porting it to a recent python environment with astropy and matplotlib. This new f2n will have a different and improved API, and should also have a much better coding style and bring new features.


.. code-block:: python 
	
	import f2n
	

... but also allows for a more sophisticated use.


The **demos** in the ``examples`` directory can be run without even installing f2n (just download them an run the scripts), and provide a quick overview.
Well, they will do so once they exist... 

Installation (in short)
-----------------------

.. code-block:: bash
	
	python setup.py install --user
	

Documentation
-------------

To learn more about how to install and use ``f2n``, proceed to its `documentation <http://f2n.readthedocs.org>`_.





