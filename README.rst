=======================
Scipion WARPhole plugin
=======================

This is a `Scipion <http://scipion.i2pc.es/>`_ plugin that allows the import of cryo-em data preprocessed in WARP `(Tegunov and Cramer, 2019) <https://pubmed.ncbi.nlm.nih.gov/31591575/>`_ in streaming mode.


Installation requirements
-------------------------

WARPhole requires the following plugins to be installed.

- `scipion-em-facilities <https://github.com/scipion-em/scipion-em-facilities>`_
- `scipion-em-relion <https://github.com/scipion-em/scipion-em-relion>`_

Installing WARPhole
-------------------

**Clone it:**

.. code-block::

    git clone https://github.com/genisvalentin/scipion-em-WARPhole.git scipion-em-WARPohole

**Install it**

.. code-block::

    scipion3 installp -p /path/to/scipion-em-WARPhole --devel

List of protocols
-----------------

- Import WARP particles (Streaming)
- Cumulative 2D classification streamer
