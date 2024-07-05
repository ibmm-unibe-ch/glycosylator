.. biobuild documentation master file, created by
   sphinx-quickstart on Tue Jun 13 14:40:03 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Installation
------------

Glycosylator can be directly installed via `pip` from the python package index:

.. code-block:: bash

   pip install glycosylator


Optional dependencies
---------------------

Glycosylator relies on `BuildAMol <https://biobuild.readthedocs.io>`_ for its core functionality. Biobuild itself has some optional but highly useful dependencies
such as `RDKit <https://www.rdkit.org/>`_ and `Openbabel <http://openbabel.org/wiki/Main_Page>`_. If you want to install these libraries you can use the following commands:

.. code-block:: bash

   pip install glycosylator[rdkit]
   pip install glycosylator[openbabel]

Glycosylator can use `GlycoWork <https://github.com/BojarLab/glycowork>`_ to draw glycans in SNFG format (Matplotlib-based native SNFG drawing is also available). To install `GlycoWork` together with `Glycosylator` use:

.. code-block:: bash

   # this is also included in full 
   pip install glycosylator[glycowork]

We recommend to install at least `glycosylator[rdkit]` if you want to use most features of `Glycosylator`. Or install everything with:

.. code-block:: bash

   pip install glycosylator[full]



Updating Glycosylator
-----------------------

.. warning:: 

   When updating Glycosylator to a newer version, make sure to export any custom default `PDBECompounds` or `CHARMMTopology` settings
   before updating, as the update will overwrite the pickle files where defaults are stored. You can re-set the custom settings after the update.

   In order to export the custom data it is recommended not to pickle the objects (as a pickled object may not be compatbile with a newer version of Biobuild).
   Instead, you can export the data as a JSON or XML file using the `to_json` or `to_xml` methods of the `PDBECompounds` and `CHARMMTopology` classes or the `export_custom_resources` function:

   .. code-block:: python

      import glycosylator as gl

      gl.export_custom_resources('my_settings')

      # ... update glycosylator

      # after updating glycosylator
      gl.import_custom_resources('my_settings')

