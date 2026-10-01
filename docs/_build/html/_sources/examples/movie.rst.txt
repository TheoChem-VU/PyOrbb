Orbital Evolution Along a Reaction Pathway
==========================================

This example showcases how PyOrbb can be used to generate a movie showing the evolution of an orbital interaction mechanism. In this case we show you how to draw the overlap of two SFOs along a reaction coordinate. We use the same system as in the `Orbital Interactions Along a Reaction Path <tracking.html>`_ example.


.. video:: orbitals.mp4
  :autoplay:
  :align: center
  :width: 100%
  :loop:
  :caption: Output video showing the evolution of the orbital overlap, for the two main interactions in the oxidative addition reaction of palladium into the C-H bond of methane.

This example uses the :meth:`~pyfmo.orbitals.objects.Orbital.cube_file` method to obtain the cube-file of a chosen SFO. By multiplying the values of this cube-file with another one we obtain the spatial overlap of the interaction between the two SFOs. Drawing this cube-file offers insight into the interaction itself. Specifically, we are able to see where in space the interaction occurs.


**Download** :download:`OxAdd_rkfs.zip <../../examples/tracking/OxAdd_rkfs.zip>`

**Download** :download:`movie.py <../../examples/movie/movie.py>`

.. literalinclude:: ../../examples/movie/movie.py
  :language: python

