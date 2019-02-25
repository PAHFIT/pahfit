#####################
How to fit a spectrum
#####################

An observed spectrum can be fit and plotted using command line scripts.

Fitting
=======

A spectrum can be fit from the command line using the ``run_pahfit`` script.
To run this script, the observed spectrum to be fit and a model pack file
are specified.

For example, to fit the Spitzer IRS SL/LL spectrum of the M101 nucleus
(included in the data directory, from Gordon et al. 2008, ApJ, 682, 336),
the command is:

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

Example Output
==============

The best fit model as well as the constraints for each PAHFIT model component
are saved in the output file.  For the example above, the output file in IPAC
table format is below.  A value of 'null' means that parameter not used
by that component and a value of 'nan' means that there is no constraint
for that specified type of constraint.

.. note:: Note that the output format is identical to the
  model pack file format, allowing the output file to be used as a model
  pack file.

::

  |     Name|           Form|  temp|temp_min|temp_max|temp_fixed|                  amp|amp_min|amp_max|amp_fixed|               x_0|           x_0_min|           x_0_max|x_0_fixed|                fwhm|           fwhm_min|           fwhm_max|fwhm_fixed|
  |     char|           char|double|  double|  double|      char|               double| double| double|     char|            double|            double|            double|     char|              double|             double|             double|      char|
  |         |               |      |        |        |          |                     |       |       |         |                  |                  |                  |         |                    |                   |                   |          |
  |     null|           null|  null|    null|    null|      null|                 null|   null|   null|     null|              null|              null|              null|     null|                null|               null|               null|      null|
         BB0     BlackBody1D 5000.0      0.0      nan       True                   0.0     0.0     nan     False               null               null               null      null                 null                null                null       null
         BB1     BlackBody1D  300.0      0.0      nan       True 3.153466447117577e-08     0.0     nan     False               null               null               null      null                 null                null                null       null
         BB2     BlackBody1D  200.0      0.0      nan       True                   0.0     0.0     nan     False               null               null               null      null                 null                null                null       null
         BB3     BlackBody1D  135.0      0.0      nan       True                   0.0     0.0     nan     False               null               null               null      null                 null                null                null       null
         BB4     BlackBody1D   90.0      0.0      nan       True 6.145012806725408e-05     0.0     nan     False               null               null               null      null                 null                null                null       null
         BB5     BlackBody1D   65.0      0.0      nan       True 0.0011702631947165707     0.0     nan     False               null               null               null      null                 null                null                null       null
         BB6     BlackBody1D   50.0      0.0      nan       True                   0.0     0.0     nan     False               null               null               null      null                 null                null                null       null
         BB7     BlackBody1D   40.0      0.0      nan       True                   0.0     0.0     nan     False               null               null               null      null                 null                null                null       null
         BB8     BlackBody1D   35.0      0.0      nan       True   0.21881606299706477     0.0     nan     False               null               null               null      null                 null                null                null       null
         DF0         Drude1D   null     null     null       null     2.645484365052569     0.0     nan     False  5.316714804508173               5.17  5.369999999999999     False  0.19709800000000002 0.16126200000000002 0.19709800000000002      False
         DF1         Drude1D   null     null     null       null    1.5018832810342535     0.0     nan     False                5.8 5.6000000000000005                5.8     False  0.21945000000000006 0.17955000000000004 0.21945000000000006      False
         DF2         Drude1D   null     null     null       null     33.56968008196337     0.0     nan     False  6.235101293291608               6.12  6.319999999999999     False  0.18628650328336455             0.16794             0.20526      False
         DF3         Drude1D   null     null     null       null      5.51067154143872     0.0     nan     False               6.79  6.590000000000001               6.79     False   0.5151300000000001 0.42147000000000007  0.5151300000000001      False
         DF4         Drude1D   null     null     null       null                   0.0     0.0     nan     False 7.3696145170208025               7.32               7.52     False   0.9058666251830154            0.841428            1.028412      False
         DF5         Drude1D   null     null     null       null     37.61600266521913     0.0     nan     False  7.572465444847651                7.5  7.699999999999999     False              0.36784             0.30096             0.36784      False
         DF6         Drude1D   null     null     null       null     36.94180562710299     0.0     nan     False   7.84267135535115               7.75  7.949999999999999     False  0.45765500000000003            0.374445 0.45765500000000003      False
         DF7         Drude1D   null     null     null       null    12.647020107750386     0.0     nan     False               8.43               8.23               8.43     False  0.45815000000000006             0.37485 0.45815000000000006      False
         DF8         Drude1D   null     null     null       null    19.984613014982852     0.0     nan     False  8.654706980092433               8.51  8.709999999999999     False             0.302211            0.302211            0.369369      False
         DF9         Drude1D   null     null     null       null                   0.0     0.0     nan     False 10.682928510528653              10.58              10.78     False   0.2055593981565941 0.19224000000000002 0.23496000000000003      False
        DF10         Drude1D   null     null     null       null    25.219257423552218     0.0     nan     False  11.22706845915868              11.13              11.33     False  0.14823600000000003 0.12128400000000002 0.14823600000000003      False
        DF11         Drude1D   null     null     null       null     37.09448152385599     0.0     nan     False 11.318084026036573              11.23              11.43     False  0.33194435265422223            0.326304            0.398816      False
        DF12         Drude1D   null     null     null       null     8.657919960775532     0.0     nan     False 11.946445966680784              11.89              12.09     False   0.5935050000000001            0.485595  0.5935050000000001      False
        DF13         Drude1D   null     null     null       null    29.054037330105995     0.0     nan     False 12.719999999999999              12.52 12.719999999999999     False             0.583044 0.47703599999999996            0.583044      False
        DF14         Drude1D   null     null     null       null                   0.0     0.0     nan     False 12.650226698666675              12.59              12.79     False  0.17113675997322458            0.148473            0.181467      False
        DF15         Drude1D   null     null     null       null    5.0908679304908535     0.0     nan     False 13.524030606924851              13.38              13.58     False  0.48528000000000004 0.48528000000000004  0.5931200000000001      False
        DF16         Drude1D   null     null     null       null                   0.0     0.0     nan     False              13.94              13.94 14.139999999999999     False             0.202176            0.202176            0.247104      False
        DF17         Drude1D   null     null     null       null    7.9287558654240105     0.0     nan     False              14.09              14.09              14.29     False  0.39022500000000004 0.31927500000000003 0.39022500000000004      False
        DF18         Drude1D   null     null     null       null                   0.0     0.0     nan     False               15.8               15.8               16.0     False               0.2862              0.2862 0.34980000000000006      False
        DF19         Drude1D   null     null     null       null    13.963254771554002     0.0     nan     False   16.4393289530547 16.349999999999998              16.55     False              0.25333             0.20727             0.25333      False
        DF20         Drude1D   null     null     null       null     20.78967529428516     0.0     nan     False 16.939999999999998 16.939999999999998              17.14     False              1.21836             0.99684             1.21836      False
        DF21         Drude1D   null     null     null       null    15.515404053708078     0.0     nan     False  17.35198296916865             17.275             17.475     False              0.22935 0.18764999999999998             0.22935      False
        DF22         Drude1D   null     null     null       null     4.879435979120135     0.0     nan     False  17.87475889386564              17.77 17.970000000000002     False             0.257328            0.257328            0.314512      False
        DF23         Drude1D   null     null     null       null    12.406040485240478     0.0     nan     False              18.82              18.82 19.020000000000003     False  0.39542800000000006 0.32353200000000004 0.39542800000000006      False
        DF24         Drude1D   null     null     null       null     32.11195065188575     0.0     nan     False               33.2               33.0               33.2     False   1.8205000000000005  1.4895000000000003  1.8205000000000005      False
     H2 S(7)      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False 5.5118264894552995             5.4615             5.5615     False 0.052819359225165605              0.0265              0.0795      False
     H2 S(6)      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False 6.1087495824290565             6.0588  6.158799999999999     False  0.05172457958546832              0.0265              0.0795      False
     H2 S(5)      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False  6.896685524267087             6.8591  6.959099999999999     False               0.0265              0.0265              0.0795      False
     H2 S(4)      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False  8.075800000000001 7.9758000000000004  8.075800000000001     False                 0.05                0.05 0.15000000000000002      False
     H2 S(3)      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False  9.666206502978518  9.614899999999999             9.7149     False  0.10335287021989095                0.05 0.15000000000000002      False
     H2 S(2)      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False 12.280091248181531 12.228499999999999            12.3285     False  0.10113180144096005                0.05 0.15000000000000002      False
     H2 S(1)      Gaussian1D   null     null     null       null    19.196343347494487     0.0     nan     False 17.004342248122136            16.9846 17.084600000000002     False  0.11160915274662513                0.07 0.21000000000000002      False
     H2 S(0)      Gaussian1D   null     null     null       null    13.143742969980572     0.0     nan     False            28.1707            28.1707            28.2707     False                 0.51 0.17000000000000004                0.51      False
      [ArII]      Gaussian1D   null     null     null       null     7.793199007434726     0.0     nan     False  6.992974398827724  6.935274000000001           7.035274     False               0.0795              0.0265              0.0795      False
     [ArIII]      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False   8.99110909043175  8.941379999999999            9.04138     False  0.09929646546844667                0.05 0.15000000000000002      False
       [SIV]      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False 10.512944650743767            10.4605 10.560500000000001     False  0.10962494230995277                0.05 0.15000000000000002      False
      [NeII]      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False             12.763             12.763 12.863000000000001     False  0.15000000000000002                0.05 0.15000000000000002      False
     [NeIII]      Gaussian1D   null     null     null       null     5.086789811762211     0.0     nan     False 15.527858177262358 15.504999999999999             15.605     False  0.21000000000000002                0.07 0.21000000000000002      False
   [SIII] 18      Gaussian1D   null     null     null       null    51.866423066964146     0.0     nan     False             18.663             18.663             18.763     False                 0.07                0.07 0.21000000000000002      False
       [OIV]      Gaussian1D   null     null     null       null                   0.0     0.0     nan     False 25.947859361764145              25.86              25.96     False  0.34401839709862964 0.17000000000000004                0.51      False
      [FeII]      Gaussian1D   null     null     null       null    15.777991928105198     0.0     nan     False 25.975983590567616             25.939             26.039     False                 0.51 0.17000000000000004                0.51      False
   [SIII] 33      Gaussian1D   null     null     null       null     160.8900405881582     0.0     nan     False 33.529999999999994              33.43 33.529999999999994     False   0.3744795023466675 0.17000000000000004                0.51      False
      [SiII]      Gaussian1D   null     null     null       null     304.2214561301282     0.0     nan     False 34.865199999999994            34.7652 34.865199999999994     False  0.17000000000000004 0.17000000000000004                0.51      False
     S07_att S07_attenuation   null     null     null       null   0.41836574891325695     0.0    10.0     False               null               null               null      null                 null                null                null       null
