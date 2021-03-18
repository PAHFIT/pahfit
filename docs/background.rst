##########
Background
##########

``PAHFIT`` was developed in 2004-2006 as part of the Spitzer Infrared Nearby Galaxy Survey (SINGS) to aid in the analysis of Spitzer/InfraRed Spectrograph (IRS) resolved spectral mapping of nearby galaxies, emphasizing the two low-resolution modules of IRS (Short-Low 5-15 microns and Long-Low 15-38 microns).  The model itself was designed to be as simple as possible, while still accommodating the physical variation in extragalactic MIR emission spectra.  The model combines simplified starlight, dust continuum (in the form of multi-temperature modified blackbodies), resolved dust emission features (e.g. PAH features, some comprised of constrained sub-features), unresolved emission lines from ions and molecules, and screen or mixed attenuation based on a pre-specified extinction curve. 

In practice, much of the effort of implementing PAHFIT was in the careful tuning of PAH feature sub-components, including their central positions and widths.  This tuning was achieved by designating 25 average, stitched 5-38 spectra extracted of ~kiloparsec scales in the centers of bright nearby galaxies.  Since it was created, PAHFIT's use has expanded and been adapted by other projects to handle high-resolution (R~600) IRS spectra, extended wavlengths from Akari (2.5-5microns), silicate emission in AGN systems, etc. With the immiment arrival of MIR spectra with high spectral resolution and sensitivity from JWST, the PAHFIT collaboration is working to bring the modeling philosophy into the JWST era.

Philosophy
------------

There are many methods to decompose spectra. PAHFIT implements a simple approach with the following principles:

- *Utilize as few parameters as possible to fit individual spectral components*.  For example at longer MIR wavelengths, starlight can be simply modeled as a T=5000K blackbody.  
- *Use physically-motivated model elements* driven by our best understanding of the emission source from Galactic and extra-galactic studies.
- *Avoid the use of observed or average "templates"*.  There are many template fitters that have different strengths, but they rely critically on the suite of templates employed.  PAHFIT is flexible enough to model atypical or unusual spectra.
- *Employ highly constrained fits*.  For example, sub-features are "mostly fixed" in position and width, to reduce the likelihood of apparently high quality fits yielding non-sensible results.

As PAHFIT's use expands to more data sets beyond Spitzer/IRS, we have recognized the need for broadening this approach, while still maintaining the core model philosophy.  We have introduced *science packs* which make different trade-offs between how constrained vs. how flexible the model will be, testing these carefully against observations.  And we have made it possible for PAHFIT users to easily explore this trade-off space themselves.  


