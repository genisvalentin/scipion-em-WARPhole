=======================
Scipion WARPhole plugin
=======================

This is a `Scipion <http://scipion.i2pc.es/>`_ plugin that allows the import of cryo-em datasets preprocessed in WARP `(Tegunov and Cramer, 2019) <https://pubmed.ncbi.nlm.nih.gov/31591575/>`_ in streaming mode.


Installation requirements
-------------------------

WARPhole requires the following plugins to be installed.

- `scipion-em-facilities <https://github.com/scipion-em/scipion-em-facilities>`_
- `scipion-em-relion <https://github.com/scipion-em/scipion-em-relion>`_

Installing WARPhole
-------------------

**Clone it:**

.. code-block::

    git clone https://github.com/genisvalentin/scipion-em-WARPhole.git scipion-em-WARPhole

**Install it**

.. code-block::

    scipion3 installp -p /path/to/scipion-em-WARPhole --devel

List of protocols
-----------------

- Import WARP particles (Streaming).
This is the main protocol of WARPhole. It reads a star file writen by WARP and imports the particles, movies, micrographs and CTFs. It can also import aligned movies, provided that the star files for "RELION bayesian polishing" are exported by WARP at the end of the data collection.

- Cumulative 2D classification streamer
This protocol extends the EMFacilites 2D streamer protocol, which packages the latest extracted particles in sets of a given size and launches 2D classification jobs. It provides the additional option "cumulative". This means, it creates particle sets that keep growing in particle number to include not only the latest extracted particles but also the older ones.

- Cumulative 3D classification streamer
Same as the "Cumulative 2D classification streamer", but for 3D classifications.

- Copy/Recover a particle set from scratch
This protocol can copy particle stack files into a scratch drive in sreaming. In this way, new paticles are directly move to the scratch drive and are ready to be used by subsequent jobs.
