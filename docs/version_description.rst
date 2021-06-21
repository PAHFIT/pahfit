###################
Version Description
###################

v1.0 - v1.4
------------

The original IDL version of PAHFIT can be found in: `http://tir.astro.utoledo.edu/jdsmith/research/pahfit.php <http://tir.astro.utoledo.edu/jdsmith/research/pahfit.php>`_. For more details and background of ``PAHFIT``, see the `background <https://pahfit.readthedocs.io/en/latest/background.html>`_ page.

v2.0
------------

A direct translation of python PAHFIT from the IDL version based on the astropy.modeling framework. Also, the modular approach has been built-in for future contributions, e.g., science/instrument packs. 

v2.1
------------
PAHFIT v2.1 is an enhanced version of v2.0 that extends the fitting range down to the near IR at 2.5 µm (total wavelength coverage: 2.5--38 µm). This version is a direct translation of IDL PAHFIT that was discussed in `Lai et al. 2020, ApJ, 905, 55 <https://iopscience.iop.org/article/10.3847/1538-4357/abc002/pdf>`_, in which the model was trained by a sample of 113 bright PAH galaxies drawn from the AKARI-Spitzer Extragalactic Spectral Survey (ASESS). 

The summary of the newly added components and the instruction of running the AKARI-Spitzer combined spectrum can be found in: :ref:`PAHFIT v2.1 summary <summary_PAHFITv21>`.

v2.2
------------
Adaptation of v2.0 to model ISO data.

v2.3
------------
Adaptation of 2.1/2.2 to JWST simulated data

v2.4
------------
Adaptation of 2.1/2.2 to JWST based on ERS data

v2.5
------------
Final ERS delivery