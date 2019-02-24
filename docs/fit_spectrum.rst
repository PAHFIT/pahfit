#####################
How to fit a spectrum
#####################

Fitting
=======

A spectrum can be fit from the commandline using the `run_pahfit` script.
To run this script, the observed spectrum to be fit and a model pack file
are specified.

For example, to fit the Spitzer IRS SL/LL spectrum of the M101 nucleus
(included in the data directory), the command is:

.. code:: shell

  $ run_pahfit scipack_ExGal_SpitzerIRSSLLL.ipac

The `run_pahfit` script first looks in the current directory for the
model pack file, then looks in the installed repository location for the
file in the packs/ subdirectory.

Plotting
========

The saved results from a PAHFIT run can be plotted using the `plot_pahfit`
script.

.. code:: shell

  $ plot_pahfit pahfit_results_file

From commandline.

From code.
