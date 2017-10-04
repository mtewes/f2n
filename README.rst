f2n v2: FITS to net, in python 3
================================

This is an attempt, under construction, to replicate the (ridiculously useful) functionality of my `f2n.py v1.4 <https://obswww.unige.ch/~tewes/f2n_dot_py/>`_ in a recent python environment, namely with astropy and matplotlib, and a much better coding style.

The purpose of the tiny ``f2n`` module is to turn FITS images (or numpy arrays) into (potentially annotated) web-friendly images (typically png or pdf if you want vector-graphics annotations) in a modular and flexible way.


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





