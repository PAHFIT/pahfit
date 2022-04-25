.. _fit_spectrum:

#####################
How to fit a spectrum
#####################

An observed spectrum can be fit and plotted using command line scripts.

Fitting
=======

A spectrum can be fit from the command line using the ``run_pahfit`` script.
To run this script, the observed spectrum to be fit and a model pack file
are specified.

The following command fits the Spitzer IRS SL/LL spectrum of the M101 
nucleus. This file can be found in the source code's ``pahfit/data`` directory
and is taken from
`Gordon et al. 2008, ApJ, 682, 336 <https://ui.adsabs.harvard.edu/abs/2008ApJ...682..336G/abstract>`_

.. code-block:: console

  $ run_pahfit M101_Nucleus_irs.ipac scipack_ExGal_SpitzerIRSSLLL.ipac

The first argument gives the file with the observed spectrum.
The second argument gives the file with the model pack file to be used.

The ``run_pahfit`` script first looks in the current directory for the
model pack file, then looks in the installed repository location for the
file in the packs/ subdirectory.  Thus, a default pack can be overridden
by copying it into the current directory and modifying it.

By default, the ``run_pahfit`` script will save the output as an ASCII
IPAC table.  In addition, the plot of the observed spectrum with the
best fit model (including components) is saved as a pdf file.  The format
of the output and plot files can be changed through command line arguments.

Help for the possible command line options for the ``run_pahfit`` script
can be seen by:

.. code-block:: console

  $ run_pahfit --help

.. _example_fit_output:

Example Output
==============

The best fit model as well as the constraints for each PAHFIT model component
are saved in the output file.  For the example above, the output file in IPAC
table format is below.  A value of `null` means that parameter not used
by that component and a value of `nan` means that there is no constraint
for that specified type of constraint.

.. note:: Note that the output format is identical to the
  model pack file format, allowing the output file to be used as a model
  pack file.

::

     Name            Form   temp temp_min temp_max temp_fixed                   amp amp_min amp_max amp_fixed                x_0  x_0_min  x_0_max x_0_fixed                fwhm            fwhm_min fwhm_max fwhm_fixed               strength strength_unc                  eqw range_min range_max

       BB0     BlackBody1D 5000.0      0.0      nan       True                   0.0     0.0     nan     False               null     null     null      null                null                null     null       null                   null         null                 null      null      null 
       BB1     BlackBody1D  300.0      0.0      nan       True 3.317302580821156e-08     0.0     nan     False               null     null     null      null                null                null     null       null                   null         null                 null      null      null 
       BB2     BlackBody1D  200.0      0.0      nan       True                   0.0     0.0     nan     False               null     null     null      null                null                null     null       null                   null         null                 null      null      null 
       BB3     BlackBody1D  135.0      0.0      nan       True                   0.0     0.0     nan     False               null     null     null      null                null                null     null       null                   null         null                 null      null      null 
       BB4     BlackBody1D   90.0      0.0      nan       True 6.662041133364268e-05     0.0     nan     False               null     null     null      null                null                null     null       null                   null         null                 null      null      null 
       BB5     BlackBody1D   65.0      0.0      nan       True 0.0011513938444125192     0.0     nan     False               null     null     null      null                null                null     null       null                   null         null                 null      null      null 
       BB6     BlackBody1D   50.0      0.0      nan       True                   0.0     0.0     nan     False               null     null     null      null                null                null     null       null                   null         null                 null      null      null 
       BB7     BlackBody1D   40.0      0.0      nan       True                   0.0     0.0     nan     False               null     null     null      null                null                null     null       null                   null         null                 null      null      null 
       BB8     BlackBody1D   35.0      0.0      nan       True    0.2147129505090326     0.0     nan     False               null     null     null      null                null                null     null       null                   null         null                 null      null      null 
       DF0         Drude1D   null     null     null       null    2.7460367541124775     0.0     nan     False               5.27     5.17     5.37      True             0.17918            0.161262 0.197098       True  8.342865218986589e-14          nan  0.22034896567020126      null      null 
       DF1         Drude1D   null     null     null       null    1.5300578571478651     0.0     nan     False                5.7      5.6      5.8      True              0.1995             0.17955  0.21945       True  4.424269992286976e-14          nan  0.10132412154615557      null      null 
       DF2         Drude1D   null     null     null       null      33.2469727649945     0.0     nan     False               6.22     6.12     6.32      True              0.1866             0.16794  0.20526       True  7.551331875850883e-13          nan   1.5615479694479546      null      null 
       DF3         Drude1D   null     null     null       null    3.3848579756129307     0.0     nan     False               6.69     6.59     6.79      True              0.4683             0.42147  0.51513       True  1.667834499856808e-13          nan  0.32973108633656245      null      null 
       DF4         Drude1D   null     null     null       null     8.035441752435286     0.0     nan     False               7.42     7.32     7.52      True             0.93492            0.841428 1.028412       True   6.42564868202564e-13          nan   1.3055841719956536      null      null 
       DF5         Drude1D   null     null     null       null     34.89960508414868     0.0     nan     False                7.6      7.5      7.7      True              0.3344             0.30096  0.36784       True  9.514811236412942e-13          nan   1.9708773082671995      null      null 
       DF6         Drude1D   null     null     null       null     32.83073102824991     0.0     nan     False               7.85     7.75     7.95      True             0.41605            0.374445 0.457655       True 1.0438241629728317e-12          nan   2.2041766010245354      null      null 
       DF7         Drude1D   null     null     null       null    7.2789822706317935     0.0     nan     False               8.33     8.23     8.43      True              0.4165             0.37485  0.45815       True 2.0574829475523607e-13          nan  0.46395043920103635      null      null 
       DF8         Drude1D   null     null     null       null     27.80548523342483     0.0     nan     False               8.61     8.51     8.71      True             0.33579            0.302211 0.369369       True  5.931062682198694e-13          nan    1.387283400311155      null      null 
       DF9         Drude1D   null     null     null       null    1.9177337009185884     0.0     nan     False              10.68    10.58    10.78      True              0.2136             0.19224  0.23496       True 1.6911713911125032e-14          nan 0.056815991658161265      null      null 
      DF10         Drude1D   null     null     null       null    32.348654451636286     0.0     nan     False              11.23    11.13    11.33      True             0.13476            0.121284 0.148236       True 1.6277896792027744e-13          nan   0.5989509984076988      null      null 
      DF11         Drude1D   null     null     null       null     35.85726139126497     0.0     nan     False              11.33    11.23    11.43      True             0.36256            0.326304 0.398816       True  4.769114677426972e-13          nan   1.7748132266199077      null      null 
      DF12         Drude1D   null     null     null       null     8.837030850837156     0.0     nan     False              11.99    11.89    12.09      True             0.53955            0.485595 0.593505       True 1.5618534519685663e-13          nan    0.628701008444088      null      null 
      DF13         Drude1D   null     null     null       null     20.36707826506126     0.0     nan     False              12.62    12.52    12.72      True             0.53004            0.477036 0.583044       True  3.191973283697224e-13          nan   1.3596504824717706      null      null 
      DF14         Drude1D   null     null     null       null    6.3919812144749875     0.0     nan     False              12.69    12.59    12.79      True             0.16497            0.148473 0.181467       True  3.083598318106959e-14          nan  0.13371451063208156      null      null 
      DF15         Drude1D   null     null     null       null     7.926079427678649     0.0     nan     False              13.48    13.38    13.58      True              0.5392             0.48528  0.59312       True 1.1075646837211678e-13          nan   0.4880337819231489      null      null 
      DF16         Drude1D   null     null     null       null    1.4569850675310168     0.0     nan     False              14.04    13.94    14.14      True             0.22464            0.202176 0.247104       True   7.81895215384728e-15          nan 0.034648794112630846      null      null 
      DF17         Drude1D   null     null     null       null     7.465101636001993     0.0     nan     False              14.19    14.09    14.29      True             0.35475            0.319275 0.390225       True   6.19346828610599e-14          nan   0.2723594926464789      null      null 
      DF18         Drude1D   null     null     null       null                   0.0     0.0     nan     False               15.9     15.8     16.0      True               0.318              0.2862   0.3498       True                    0.0          nan                  0.0      null      null 
      DF19         Drude1D   null     null     null       null    16.347471641301407     0.0     nan     False              16.45    16.35    16.55      True              0.2303             0.20727  0.25333       True  6.551689484824348e-14          nan   0.2389033406506008      null      null 
      DF20         Drude1D   null     null     null       null    25.584109293688915     0.0     nan     False              17.04    16.94    17.14      True              1.1076             0.99684  1.21836       True  4.595731502350026e-13          nan    1.445339582606487      null      null 
      DF21         Drude1D   null     null     null       null     9.494902722179587     0.0     nan     False             17.375   17.275   17.475      True              0.2085             0.18765  0.22935       True    3.0880728162438e-14          nan  0.10110748258068253      null      null 
      DF22         Drude1D   null     null     null       null     3.964662252087014     0.0     nan     False              17.87    17.77    17.97      True             0.28592            0.257328 0.314512       True  1.671637888079228e-14          nan  0.05171600901384207      null      null 
      DF23         Drude1D   null     null     null       null     4.685189400059929     0.0     nan     False              18.92    18.82    19.02      True             0.35948            0.323532 0.395428       True  2.215645121229816e-14          nan 0.061009285333322755      null      null 
      DF24         Drude1D   null     null     null       null    22.781357965422234     0.0     nan     False               33.1     33.0     33.2      True               1.655              1.4895   1.8205       True  1.620549117711164e-13          nan   0.2118931863621897      null      null 
   H2 S(7)      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False  5.526811958245635   5.4615   5.5615     False              0.0265              0.0265   0.0795      False                    0.0          nan                  0.0      null      null 
   H2 S(6)      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False  6.113472648913736   6.0588   6.1588     False              0.0265              0.0265   0.0795      False                    0.0          nan                  0.0      null      null 
   H2 S(5)      Gaussian1D   null     null     null       null     4.663185944000933     0.0     nan     False  6.901018306044247   6.8591   6.9591     False              0.0265              0.0265   0.0795      False  8.279840560439328e-15          nan 0.017291542192613564      null      null 
   H2 S(4)      Gaussian1D   null     null     null       null      2.82012745235727     0.0     nan     False             8.0758   7.9758   8.0758     False  0.1346978300964964                0.05     0.15      False  1.858565843530547e-14          nan 0.042308594049899076      null      null 
   H2 S(3)      Gaussian1D   null     null     null       null     4.904783121362579     0.0     nan     False  9.711548067976954   9.6149   9.7149     False                0.15                0.05     0.15      False 2.4891660280966613e-14          nan  0.07407542335335081      null      null 
   H2 S(2)      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False            12.3285  12.2285  12.3285     False                0.15                0.05     0.15      False                    0.0          nan                  0.0      null      null 
   H2 S(1)      Gaussian1D   null     null     null       null     39.89580207188739     0.0     nan     False 17.000279628369974  16.9846  17.0846     False                0.07                0.07     0.21      False  3.083421075027653e-14          nan  0.11126351021423353      null      null 
   H2 S(0)      Gaussian1D   null     null     null       null    12.030485071120175     0.0     nan     False            28.1707  28.1707  28.2707     False                0.51 0.17000000000000004     0.51      False  2.467048928617314e-14          nan   0.0396960989965039      null      null 
    [ArII]      Gaussian1D   null     null     null       null    23.433617975392888     0.0     nan     False  6.986653221859022 6.935274 7.035274     False 0.03232551854094294              0.0265   0.0795      False    4.9518353490638e-14          nan  0.10352016750769705      null      null 
   [ArIII]      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False  8.986273164753749  8.94138  9.04138     False 0.10655551354809757                0.05     0.15      False                    0.0          nan                  0.0      null      null 
     [SIV]      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False  10.51899527488704  10.4605  10.5605     False                0.05                0.05     0.15      False                    0.0          nan                  0.0      null      null 
    [NeII]      Gaussian1D   null     null     null       null    31.066497094306175     0.0     nan     False 12.829135846841778   12.763   12.863     False                0.15                0.05     0.15      False  9.034590322165919e-14          nan  0.41775387519580875      null      null 
   [NeIII]      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False             15.505   15.505   15.605     False                0.21                0.07     0.21      False                    0.0          nan                  0.0      null      null 
 [SIII] 18      Gaussian1D   null     null     null       null    32.689858803743675     0.0     nan     False 18.732947618916615   18.663   18.763     False 0.15872885615864357                0.07     0.21      False  4.718200795811459e-14          nan  0.13941146302452906      null      null 
     [OIV]      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False              25.96    25.86    25.96     False                0.51 0.17000000000000004     0.51      False                    0.0          nan                  0.0      null      null 
    [FeII]      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False             25.939   25.939   26.039     False                0.51 0.17000000000000004     0.51      False                    0.0          nan                  0.0      null      null 
 [SIII] 33      Gaussian1D   null     null     null       null     141.3020935175383     0.0     nan     False              33.53    33.43    33.53     False                0.51 0.17000000000000004     0.51      False 2.0453676899151097e-13          nan   0.2748590976423317      null      null 
    [SiII]      Gaussian1D   null     null     null       null     306.3268736329564     0.0     nan     False            34.8652  34.7652  34.8652     False  0.2675436953319326 0.17000000000000004     0.51      False  2.151370572328281e-13          nan   0.2769147260816285      null      null 
   S07_att S07_attenuation   null     null     null       null     0.641957924103409     0.0    10.0     False               null     null     null      null                null                null     null       null                   null         null                 null      null      null 
    PAH_62            null   null     null     null       null                  null    null    null      null               null     null     null      null                null                null     null       null  7.551331875850883e-13          nan   1.5615479694479546       6.2       6.3 
  PAH_77_C            null   null     null     null       null                  null    null    null      null               null     null     null      null                null                null     null       null 2.6378701548166897e-12          nan    5.480638081287388       7.3       7.9 
    PAH_83            null   null     null     null       null                  null    null    null      null               null     null     null      null                null                null     null       null 2.0574829475523607e-13          nan  0.46395043920103635       8.3       8.4 
    PAH_86            null   null     null     null       null                  null    null    null      null               null     null     null      null                null                null     null       null  5.931062682198694e-13          nan    1.387283400311155       8.6       8.7 
 PAH_113_C            null   null     null     null       null                  null    null    null      null               null     null     null      null                null                null     null       null  6.396904356629747e-13          nan   2.3737642250276068      11.2      11.4 
   PAH_120            null   null     null     null       null                  null    null    null      null               null     null     null      null                null                null     null       null 1.5618534519685663e-13          nan    0.628701008444088      11.9      12.1 
 PAH_126_C            null   null     null     null       null                  null    null    null      null               null     null     null      null                null                null     null       null   3.50033311550792e-13          nan   1.4933649931038522      12.6      12.7 
   PAH_136            null   null     null     null       null                  null    null    null      null               null     null     null      null                null                null     null       null 1.1075646837211678e-13          nan   0.4880337819231489      13.4      13.6 
   PAH_142            null   null     null     null       null                  null    null    null      null               null     null     null      null                null                null     null       null   6.19346828610599e-14          nan   0.2723594926464789      14.1      14.2 
   PAH_164            null   null     null     null       null                  null    null    null      null               null     null     null      null                null                null     null       null  6.551689484824348e-14          nan   0.2389033406506008      16.4      16.5 
  PAH_17_C            null   null     null     null       null                  null    null    null      null               null     null     null      null                null                null     null       null  5.726871521264763e-13          nan   1.8370664148516123      16.4      17.9 
   PAH_174            null   null     null     null       null                  null    null    null      null               null     null     null      null                null                null     null       null    3.0880728162438e-14          nan  0.10110748258068253     17.35     17.45