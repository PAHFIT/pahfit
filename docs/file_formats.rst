############
File Formats
############

Model Packs
===========

The model packs are given in a table with a row for each model component.
For each fit parameter, there are columns giving
the `Name` (e.g., BB3, DF10, H2 S(3), [OIV]),
the functional `Form` (e.g., BlackBody1D, Drude1D, Gaussian1D),
and multiple columns for the fit parameters.
For each fit parameter, the columns are for the actual value,
the bounds (i.e., min and max), and if the parameter is fixed or not.
For the `BlackBody1D` models, the fit parameters are temperature (`temp`)
and  "intrinsic" or "unattenuated" amplitude (`amp`).
For the `Drude1D` and `Gaussian1D` models, the fit parameters are
"intrinsic" or "unattenuated" amplitude (`amp`), center wavelength (`x_0`), and
full width at half maximum (`fwhm`).
For the dust attenuation model (`S07_att`), the fit parameters are
amplitude (`amp`).

The astropy version of `Drude1D
<https://docs.astropy.org/en/stable/modeling/physical_models.html#drude1d>`_
and `Gaussian1D
<https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.Gaussian1D.html#astropy.modeling.functional_models.Gaussian1D>`_
are used. See :ref:`reference API <reference_API>` for a
description of BlackBody1D and the dust attenuation model.

The model packs are saved as IPAC tables as these provide readable (ASCII)
files that include unit information.

See :ref:`example fit output <example_fit_output>` for an example of
a model pack/output file.

.. note::
   IDL PAHFIT uses fractional FWHM as input/output for the Drude profiles, python
   PAHFIT uses FWHM as input/output for the Drude profiles.

Fit Outputs
===========

Same as model packs.  This enables the fit output file to be used as a
model pack file, potentially useful for fitting multiple spectra that are
similar (e.g., JWST spectral cubes).
