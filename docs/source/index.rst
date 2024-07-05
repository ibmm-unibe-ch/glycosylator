.. <gallery>

.. raw:: html

    <div class="slideshow-container">
        
    <div class="mySlides fade">
        <img src="_static/gallery/protglyco4.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/protglyco2.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/protglyco3.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/protglyco1.png" style="width:100%">
    </div>

    <div class="mySlides fade">
        <img src="_static/gallery/glycoshield1.png" style="width:100%">
    </div>

    </div>

.. raw:: html

    <style>
    .slideshow-container {
      max-width: 1000px;
      position: relative;
      margin: auto;
    }

   img {
      background-color: transparent !important;
   }

    .mySlides {
      display: none;
    }

    .fade {
      animation: fade 5000ms infinite;
    }

    @keyframes fade {
      from {opacity: .7} 
      to {opacity: 1}
    }
    </style>

.. raw:: html

    <script>
    var slideIndex = 0;
    showSlides();

    function showSlides() {
      var i;
      var slides = document.getElementsByClassName("mySlides");
      for (i = 0; i < slides.length; i++) {
        slides[i].style.display = "none";  
      }
      slideIndex++;
      if (slideIndex > slides.length) {slideIndex = 1}    
      slides[slideIndex-1].style.display = "block";  
      setTimeout(showSlides, 5000); // Change image every N seconds
    }
    </script>

.. <gallery>
.. Glycosylator documentation master file, created by
   sphinx-quickstart on Wed Dec 21 15:02:18 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. image:: _resources/logo_large.png
   :width: 80%
   :align: center
   :alt: glycosylator logo




Welcome to Glycosylator's documentation!
========================================

Glycosylator is a Python package to create glycan structures and glycosylate biomolecules. It is designed to be a flexible and extensible, intended to be used in conjunction with other tools for the analysis of glycan structures.
Using Glycosylator you can not only quickly generate 3D glycan models from IUPAC strings, but you can also modify existing glycans, or generate completely new or non-canonical glycans from scratch!

.. note::

   This documentation has content on `Glycosylator`. Glycosylator is built on top of `BuildAmol` so a majority of the documentation is shared between the two.
   Be sure to check out the `BuildAMol documentation <https://biobuild.readthedocs.io/en/latest/>`_ for more information on the underlying structure (and more general features) of Glycosylator.

.. .. admonition:: Build glycans 
      
..    .. image:: _resources/glycan_light.gif
..       :width: 100%
..       :align: center
..       :alt: a glycan made with glycosylator
..       :class: only-light

..    .. image:: _resources/glycan_dark.gif
..       :width: 100%
..       :align: center
..       :alt: a glycan made with glycosylator
..       :class: only-dark


.. .. admonition:: Glycosylate proteins
      

..    .. image:: _resources/prot_light.gif
..       :width: 100%
..       :align: center
..       :alt: a protein glycosylated with glycosylator
..       :class: only-light

..    .. image:: _resources/prot_dark.gif
..       :width: 100%
..       :align: center
..       :alt: a protein glycosylated with glycosylator
..       :class: only-dark

.. .. admonition:: Glycosylate membranes

..    .. image:: _resources/mem_light.gif
..       :width: 100%
..       :align: center
..       :alt: a membrane glycosylated with glycosylator
..       :class: only-light

..    .. image:: _resources/mem_dark.gif
..       :width: 100%
..       :align: center
..       :alt: a membrane glycosylated with glycosylator
..       :class: only-dark


.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :hidden:

   installation
   tutorials
   modules