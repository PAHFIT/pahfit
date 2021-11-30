######
PAHFIT
######

``pahfit`` is a python package providing a decomposition model
for astronomical infrared spectra, focusing on dust and gas emission
features from the interstellar medium.

The original versions of PAHFIT (v1.x) were written in IDL and
focused mainly on Spitzer/IRS spectroscopic observations.
This python-based versions (>=v2.0) will provide expanded capabilities
including more instrument (e.g., AKARI and JWST)
and a more flexible modeling framework suitable for modeling a wider range
of astrophysical sources.

For details for the IDL version of PAHFIT see
`Smith, J.D.T., Draine B.T., et al., 2007, ApJ, 656, 770 <http://tir.astro.utoledo.edu/jdsmith/research/pahfit.php>`_.

This package is potentially an
`astropy affiliated package <http://www.astropy.org/affiliated/>`_
and uses the
`astropy.modeling <http://docs.astropy.org/en/stable/modeling/>`_
framework.

User Documentation
==================

.. toctree::
   :maxdepth: 2

   Background <background.rst>
   Fitting a spectrum <fit_spectrum.rst>
   Plotting a fit <plot_result.rst>
   File formats <file_formats.rst>
   IDL versus Python PAHFIT implementations <idl_vs_python/python_vs_IDL_fit.rst>

Installation
============

.. toctree::
  :maxdepth: 2

  How to install <install.rst>

Repository
==========

GitHub: `pahfit <https://github.com/PAHFIT/pahfit>`_

Quick Start
===========

Text, plots, and (potentially) code to quickly get a user fitting a spectrum.

Reporting Issues
================

If you have found a bug in ``pahfit`` please report it by creating a
new issue on the ``pahfit`` `GitHub issue tracker
<https://github.com/PAHFIT/pahfit/issues>`_.

Please include an example that demonstrates the issue sufficiently so that the
developers can reproduce and fix the problem.

Contributing
============

Like the `Astropy`_ project, ``pahfit`` is made both by and for its
users.  We accept contributions at all levels, spanning the gamut from fixing a
typo in the documentation to developing a major new feature. We welcome
contributors who will abide by the `Python Software Foundation Code of Conduct
<https://www.python.org/psf/conduct/>`_.

``pahfit`` follows the same workflow and coding guidelines as
`Astropy`_.  The following pages will help you get started with contributing
fixes, code, or documentation (no git or GitHub experience necessary):

* `How to make a code contribution <https://docs.astropy.org/en/latest/development/workflow/development_workflow.html>`_

* `Coding Guidelines <https://docs.astropy.org/en/latest/development/codeguide.html>`_

* `Developer Documentation <https://docs.astropy.org/en/latest/#developer-documentation>`_

For the complete list of contributors please see the `pahfit
contributors page on Github
<https://github.com/karllark/pahfit/graphs/contributors>`_.

.. _reference_API:

Reference API
=============

Base class for PAHFIT

.. automodapi:: pahfit.base

Component models not provided by astropy.models.

.. automodapi:: pahfit.component_models
