# Description 

This software simulates astronomical survey detections of tidal disruptions of stars by super-massive black holes. It begins with the synthetic galaxy catalogue described in van Velzen 2008 (https://arxiv.org/abs/1707.03458). The stellar disruption rate in each galaxy is estimated based on Stone & Metzger 2016 (https://arxiv.org/abs/1410.7772), or other details that the user may specify. Based on these rates, and the present-day stellar mass function in the galaxy, disruptions are randomly sampled, and the properties of the resulting flares are sampled based on empirical distributions. The code also accounts for obscuration by dust in the host galaxy. Finally, the survey selection effects are applied. The detectable simulated flares are stored in a database, allowing histograms of their properties to be created.

# Getting Started

The code requires the use of the GNU scientific library and HDF5. The galaxy catalogue described in https://arxiv.org/abs/1707.03458 must be downloaded and properly formatted as an HDF5 file (this can be provided upon request to roth14@llnl.gov). To change default astronomical parameters, modify the constructors of the galaxy, disruption, and survey classes. The main routine which creates ntuples of flare and galaxy data is contained in main_bin_sjoert_catalogue.cpp. Routines for binning these ntuples are found in the main_make_galaxy_histograms files. Examples of how to view the histograms are provided in the ipython notebooks.

# Contributing

Please fork the repository and then submit pull requests.

# Release

LLNL-CODE- 816642
Title: tdegb - Tidal Disruption Event Galaxy Binner, Version: 1.0
Author(s) Nathaniel J. Roth








