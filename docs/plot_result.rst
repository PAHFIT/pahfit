####################
Plotting fit results
####################

Plotting
========

The saved results from a PAHFIT run can be plotted using the ``plot_pahfit``
script.  By default, the plot is displayed.  To save the plot instead,
use the "--savefig" command line option.

.. code-block:: console

  $ plot_pahfit M101_Nucleus_irs.ipac M101_Nucleus_irs_output.ipac

The first argument gives the file with the observed spectrum.
The second argument gives the file with the output of the ``run_pahfit``
script.

Help for the possible command line options for the ``plot_pahfit`` script
can be seen by:

.. code-block:: console

  $ plot_pahfit --help

Example Plot
============

The fit to the M101 Nucleus spectrum gives the observed spectrum as
open squares, the full PAHFIT model as a green line, the dust continua
as red lines, the sum of the continua as a black line,
the dust features as blue lines, the atomic features
as purple lines, and the attenuation as a dashed, black line.

.. image:: M101_Nucleus_irs.png
