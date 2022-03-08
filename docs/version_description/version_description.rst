###################
Version Description
###################

v1.0 - v1.4
------------

The original IDL version of PAHFIT can be found in: `http://tir.astro.utoledo.edu/jdsmith/research/pahfit.php <http://tir.astro.utoledo.edu/jdsmith/research/pahfit.php>`_. For more details and background of ``PAHFIT``, see the `background <https://pahfit.readthedocs.io/en/latest/background.html>`_ page.

v2.0
------------

Adopting the astropy.modeling framework, this version of python PAHFIT matches the capability of the IDL version. Also, the modular approach has been built-in for future contributions, e.g., science/instrument packs. 

v2.1
------------
PAHFIT v2.1 is an enhanced version of v2.0 that extends the fitting range down to the near IR at 2.5 µm (total wavelength coverage: 2.5--38 µm). The detail of this version is discussed in `Lai et al. 2020, ApJ, 905, 55 <https://iopscience.iop.org/article/10.3847/1538-4357/abc002/pdf>`_, in which the model was trained by a sample of 113 bright PAH galaxies drawn from the AKARI-Spitzer Extragalactic Spectral Survey (ASESS). 

The summary of the newly added components and the instruction of running the AKARI-Spitzer combined spectrum can be found in: `PAHFIT v2.1 notes <https://github.com/PAHFIT/pahfit/blob/master/docs/version_description/PAHFIT_v2.1.rst>`_.

v2.2
------------
Enhancement to support Infrared Space Observatory spectra

v2.3
------------
Support the JWST simulated data

v2.4
------------
Support the JWST observations
