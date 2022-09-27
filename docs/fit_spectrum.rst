.. _fit_spectrum:

#####################
How to fit a spectrum
#####################

An observed spectrum can be fit and plotted using an
interactive python session such as notebooks, using a custom script written
by the user, or using the provided convernience script.

Fitting
=======

Command line
------------

A spectrum can be fit from the command line using the ``run_pahfit`` script.
To run this script, the observed spectrum to be fit and a model pack file
are specified.

The following command fits the Spitzer IRS SL/LL spectrum of the M101
nucleus. This file can be found in the source code's ``pahfit/data`` directory
and is taken from
`Gordon et al. 2008, ApJ, 682, 336 <https://ui.adsabs.harvard.edu/abs/2008ApJ...682..336G/abstract>`_

.. code-block:: console

  $ run_pahfit M101_Nucleus_irs.ipac classic.yaml 'spitzer.irs.*.[12]'

The first argument gives the file with the observed spectrum.
The second argument gives the file with the science pack file to be used.
The third argument is a string representing the instrument pack to be used.
In this case, a model that combines 4 segments from Spitzer IRS: LL1, LL2, SL1, SL2.
See file contents of each instrument pack <file_formats.rst> for more details.

The ``run_pahfit`` script first looks in the current directory for the
model pack file, then looks in the installed repository location for the
file in the packs/ subdirectory.  Thus, a default pack can be overridden
by copying it into the current directory and modifying it.

By default, the ``run_pahfit`` script will save the output as an ASCII ECSV
table. In addition, the plot of the observed spectrum with the best fit model
(including components) is saved as a pdf file. The format of the output and
plot files can be changed through command line arguments.

Help for the possible command line options for the ``run_pahfit`` script
can be seen by:

.. code-block:: console

  $ run_pahfit --help

.. _example_fit_output:

Python notebook and scripts
---------------------------

More advanced users, such as those developing and testing their own science and
instrument packs, are encouraged to use the Model class directly. This class
acts as the main Python API.

The first step is to set up an instance of the model class. All that is needed,
is a science pack and the instrument configuration the model should represent.
See File Formats <file_formats.rst>.

.. code-block:: python

   from pahfit.model import Model
   model = Model.from_yaml('classic.yaml', 'spitzer.irs.*.[12]')

This object will then be used to perform the initial guess (if needed), the
fitting, and rudimentary analysis of the fit results.

For guessing and fitting, observational data needs to be provided in the form
of a ``Spectrum1D`` object. Any input format supported by the specutils package
will work, as long as it has assigned units and uncertainties.

The code below illustrates a minimal workflow.

.. code-block:: python

   # load in an example spectrum
   from specutils import Spectrum1D
   spec = Spectrum1D.read(file_name)
   # perform initial guess
   model.guess(spec)
   # perform fit
   model.fit(spec)
   # display plot
   model.plot(spec)
   # save results
   model.save('my_output.ecsv')

It is not required that the same spectrum object is passed every time. E.g. it
is perfectly reasonable to make the plot with another spectrum than the one
used for the fit, for comparative purposes.

The more advanced user can experiment with different initial conditions by
directly accessing the features table stored in the model object. As long as no
columns are added or removed, the user can freely edit the parameters of this
table or remove rows. Each row corresponds to one feature. Any edits to this
table will immediately reflected when the next call to guess, fit, or plot is
made.

.. code-block:: python

   model.features.loc['H2_S(7)'].power[0] = new_value

In the above example, the strength of a line is adjusted manually. The index
``[0]`` at the end accesses the value, while ``[1]`` and ``[2]`` access the
lower and upper bounds. See File Formats <file_formats.rst>. It is also
possible to use a saved model as the initial guess for another fit. This can be
done by loading the saved model as follows

..code-black:: python

  model = model.from_saved('my_ouput.ecsv')


Example Output
==============

The parameters of the best fit model and the constraint setup for each PAHFIT
model component are saved in the output file by using
``model.save(file_name)``. For the example above, the output file in ECSV table
format is below. A value of `null` means that parameter not used by that
component, and each 3-tuple represents [value, min, max]. The min/max bounds
are not changed during the fitting, but they are saved for reference. Bounds
set to `null` indicate that the parameter was fixed during the fit. Any extra
needed to reload the model from this file, is stored in the ECSV metadata.

::

   # %ECSV 1.0
   # ---
   # datatype:
   # - {name: name, datatype: string}
   # - {name: group, datatype: string}
   # - {name: kind, datatype: string}
   # - {name: temperature, unit: K, datatype: string, format: 0.4g, subtype: 'int64[3]'}
   # - {name: tau, datatype: string, format: 0.4g, subtype: 'float64[3]'}
   # - {name: wavelength, unit: um, datatype: string, format: 0.4g, subtype: 'float64[3]'}
   # - {name: power, datatype: string, format: 0.4g, subtype: 'float64[3]'}
   # - {name: fwhm, unit: um, datatype: string, format: 0.4g, subtype: 'float64[3]'}
   # - {name: model, datatype: string}
   # - {name: geometry, datatype: string}
   # meta: !!omap
   # - {redshift: 0}
   # - {instrumentname: 'spitzer.irs.*.[12]'}
   # schema: astropy-2.0
   name group kind temperature tau wavelength power fwhm model geometry
   starlight _none_ starlight [5000,null,null] [0.0,0.0,Infinity] [null,null,null] [null,null,null] [null,null,null] "" ""
   dust_cont00 dust_cont dust_continuum [300,null,null] [2.6442804784925465e-08,0.0,Infinity] [null,null,null] [null,null,null] [null,null,null] "" ""
   dust_cont01 dust_cont dust_continuum [200,null,null] [3.5775925271387995e-08,0.0,Infinity] [null,null,null] [null,null,null] [null,null,null] "" ""
   dust_cont02 dust_cont dust_continuum [135,null,null] [0.0,0.0,Infinity] [null,null,null] [null,null,null] [null,null,null] "" ""
   dust_cont03 dust_cont dust_continuum [90,null,null] [6.357043327056842e-05,0.0,Infinity] [null,null,null] [null,null,null] [null,null,null] "" ""
   dust_cont04 dust_cont dust_continuum [65,null,null] [0.00112886154478198,0.0,Infinity] [null,null,null] [null,null,null] [null,null,null] "" ""
   dust_cont05 dust_cont dust_continuum [50,null,null] [0.0003707028100634872,0.0,Infinity] [null,null,null] [null,null,null] [null,null,null] "" ""
   dust_cont06 dust_cont dust_continuum [40,null,null] [0.0,0.0,Infinity] [null,null,null] [null,null,null] [null,null,null] "" ""
   dust_cont07 dust_cont dust_continuum [35,null,null] [0.20590716546482066,0.0,Infinity] [null,null,null] [null,null,null] [null,null,null] "" ""
   H2_S(7) H2_lines line [null,null,null] [null,null,null] [5.5115,null,null] [0.0,0.0,Infinity] [null,null,null] "" ""
   H2_S(6) H2_lines line [null,null,null] [null,null,null] [6.1088,null,null] [0.0,0.0,Infinity] [null,null,null] "" ""
   H2_S(5) H2_lines line [null,null,null] [null,null,null] [6.9091,null,null] [1.7606418031978688,0.0,Infinity] [null,null,null] "" ""
   H2_S(4) H2_lines line [null,null,null] [null,null,null] [8.0258,null,null] [2.5661625021977463,0.0,Infinity] [null,null,null] "" ""
   H2_S(3) H2_lines line [null,null,null] [null,null,null] [9.6649,null,null] [3.9004292961793996,0.0,Infinity] [null,null,null] "" ""
   H2_S(2) H2_lines line [null,null,null] [null,null,null] [12.2785,null,null] [0.5683901266992609,0.0,Infinity] [null,null,null] "" ""
   H2_S(1) H2_lines line [null,null,null] [null,null,null] [17.0346,null,null] [8.747111041383896,0.0,Infinity] [null,null,null] "" ""
   H2_S(0) H2_lines line [null,null,null] [null,null,null] [28.2207,null,null] [13.538292045919352,0.0,Infinity] [null,null,null] "" ""
   [ArII] ionic_lines line [null,null,null] [null,null,null] [6.985274,null,null] [8.36886444129017,0.0,Infinity] [null,null,null] "" ""
   [ArIII] ionic_lines line [null,null,null] [null,null,null] [8.99138,null,null] [0.0,0.0,Infinity] [null,null,null] "" ""
   [SIV] ionic_lines line [null,null,null] [null,null,null] [10.5105,null,null] [0.0,0.0,Infinity] [null,null,null] "" ""
   [NeII] ionic_lines line [null,null,null] [null,null,null] [12.813,null,null] [33.19571224676041,0.0,Infinity] [null,null,null] "" ""
   [NeIII] ionic_lines line [null,null,null] [null,null,null] [15.555,null,null] [5.711816416887488,0.0,Infinity] [null,null,null] "" ""
   [SIII]_18 ionic_lines line [null,null,null] [null,null,null] [18.713,null,null] [31.455218561451822,0.0,Infinity] [null,null,null] "" ""
   [OIV] ionic_lines line [null,null,null] [null,null,null] [25.91,null,null] [8.186354417898814,0.0,Infinity] [null,null,null] "" ""
   [FeII]_26 ionic_lines line [null,null,null] [null,null,null] [25.989,null,null] [10.389897053987141,0.0,Infinity] [null,null,null] "" ""
   [SIII]_33 ionic_lines line [null,null,null] [null,null,null] [33.48,null,null] [158.62512109328276,0.0,Infinity] [null,null,null] "" ""
   [SiII] ionic_lines line [null,null,null] [null,null,null] [34.8152,null,null] [171.562392744096,0.0,Infinity] [null,null,null] "" ""
   [FeII]_35 ionic_lines line [null,null,null] [null,null,null] [35.349,null,null] [19.53628016008112,0.0,Infinity] [null,null,null] "" ""
   PAH_5.3 _none_ dust_feature [null,null,null] [null,null,null] [5.27,null,null] [3.210029634102234,0.0,Infinity] [0.17918,null,null] "" ""
   PAH_5.7 _none_ dust_feature [null,null,null] [null,null,null] [5.7,null,null] [2.630206267317981,0.0,Infinity] [0.1995,null,null] "" ""
   PAH_6.2 _none_ dust_feature [null,null,null] [null,null,null] [6.22,null,null] [33.55498609501776,0.0,Infinity] [0.1866,null,null] "" ""
   PAH_6.7 _none_ dust_feature [null,null,null] [null,null,null] [6.69,null,null] [4.0459385356533115,0.0,Infinity] [0.4683,null,null] "" ""
   PAH_7.7a PAH_7.7_cmp dust_feature [null,null,null] [null,null,null] [7.42,null,null] [9.189964220407983,0.0,Infinity] [0.93492,null,null] "" ""
   PAH_7.7b PAH_7.7_cmp dust_feature [null,null,null] [null,null,null] [7.6,null,null] [33.715242580350626,0.0,Infinity] [0.3344,null,null] "" ""
   PAH_7.7c PAH_7.7_cmp dust_feature [null,null,null] [null,null,null] [7.85,null,null] [32.68649055207141,0.0,Infinity] [0.41605,null,null] "" ""
   PAH_8.3 _none_ dust_feature [null,null,null] [null,null,null] [8.33,null,null] [8.1592792362296,0.0,Infinity] [0.4165,null,null] "" ""
   PAH_8.6 _none_ dust_feature [null,null,null] [null,null,null] [8.61,null,null] [26.67676132226628,0.0,Infinity] [0.33579,null,null] "" ""
   PAH_10.7 _none_ dust_feature [null,null,null] [null,null,null] [10.68,null,null] [1.6060415304828528,0.0,Infinity] [0.2136,null,null] "" ""
   PAH_11.3a PAH_11.3_cmp dust_feature [null,null,null] [null,null,null] [11.23,null,null] [30.13327437546274,0.0,Infinity] [0.13476,null,null] "" ""
   PAH_11.3b PAH_11.3_cmp dust_feature [null,null,null] [null,null,null] [11.33,null,null] [34.508405309456315,0.0,Infinity] [0.36256,null,null] "" ""
   PAH_12 _none_ dust_feature [null,null,null] [null,null,null] [11.99,null,null] [8.77172976153565,0.0,Infinity] [0.53955,null,null] "" ""
   PAH_12.6a PAH_12.6_cmp dust_feature [null,null,null] [null,null,null] [12.62,null,null] [21.309378952957225,0.0,Infinity] [0.53004,null,null] "" ""
   PAH_12.6b PAH_12.6_cmp dust_feature [null,null,null] [null,null,null] [12.69,null,null] [5.941460629029009,0.0,Infinity] [0.16497,null,null] "" ""
   PAH_13.48 _none_ dust_feature [null,null,null] [null,null,null] [13.48,null,null] [8.092922236753854,0.0,Infinity] [0.5392,null,null] "" ""
   PAH_14.04 _none_ dust_feature [null,null,null] [null,null,null] [14.04,null,null] [1.4771257015271373,0.0,Infinity] [0.22464,null,null] "" ""
   PAH_14.19 _none_ dust_feature [null,null,null] [null,null,null] [14.19,null,null] [7.4556013813599,0.0,Infinity] [0.35475,null,null] "" ""
   PAH_15.9 _none_ dust_feature [null,null,null] [null,null,null] [15.9,null,null] [0.0,0.0,Infinity] [0.318,null,null] "" ""
   PAH_17a PAH_17_cmp dust_feature [null,null,null] [null,null,null] [16.45,null,null] [16.34170791195906,0.0,Infinity] [0.2303,null,null] "" ""
   PAH_17b PAH_17_cmp dust_feature [null,null,null] [null,null,null] [17.04,null,null] [23.255105779785872,0.0,Infinity] [1.1076,null,null] "" ""
   PAH_17c PAH_17_cmp dust_feature [null,null,null] [null,null,null] [17.37,null,null] [9.414183385885616,0.0,Infinity] [0.2085,null,null] "" ""
   PAH_17d PAH_17_cmp dust_feature [null,null,null] [null,null,null] [17.87,null,null] [3.1748019196040747,0.0,Infinity] [0.28592,null,null] "" ""
   PAH_33.1 _none_ dust_feature [null,null,null] [null,null,null] [33.1,null,null] [30.124697454332953,0.0,Infinity] [1.655,null,null] "" ""
   silicate _none_ attenuation [null,null,null] [0.40038366018969024,0.0,Infinity] [null,null,null] [null,null,null] [null,null,null] S07_attenuation mixed
