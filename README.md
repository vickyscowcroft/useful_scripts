# useful_scripts
Useful little scripts that I use a lot
---

* ads_importer.py: This script grabs all the bibliographic info from ads and puts it into a latex bibliography file.
Reference all your citeations (\citet etc.) as the ads bibcode and it will match by that.
Usage: ads_importer.py latexfile
oputputs the file you input to bibliography{}. Note that the latexfile input has no extension!!
* annotated_lightcurve.py: Makes a lightcurve of a star and annotates it. Good script for making pretty lightcurves. Examples of filling between lines, annotating a plot, and various other useful matplotlib stuff.
* inverse_wavelength_mu.py: Makes a plot of distance modulus vs inverse wavelength. Does the reddening law fit using CCM and Indebetouw laws at the corresponding wavelengths, and adds a second plot below showing residuals.
* make_smc_image.py: Example of a finder chart. Main image is the main region of the SMC with Cepheid coordinates labelled. Inset is the whole galaxy extended to show the wing, with all Cepheids shown. Uses aplpy and montage. Required files are in finder_chart_example folder.
* reddening_laws.py: CCM and Indebetouw reddning law fitters.
