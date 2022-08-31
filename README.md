# risk-averse-R-numbers

We present Matlab code to reproduce all analyses and figures from the paper "Risk averse reproduction numbers improve resurgence detection" by Kris V Parag and Uri Obolski. This paper is available in preprint form at:

We also present a function in R which users can easily modify to apply to custom datasets. This can be found in the folder R version and includes an update of EpiFilter to compute E. The main function empiricalEandR.R is setup to reproduce part of a COVID-19 case study in Israel.

This work focuses on designing new consensus statistics for describing transmissibility at larger scales where homogeneous mixing assumptions are likely invalid. We assume that incidence curves at the local scale (e.g., regions composing a country) are available and that homogeneous mixing is valid at this scale. We estimate local reproduction numbers from each incidence curve i.e., Rj for local group j with 1 <= j <= p groups. We then combine these to derive three statistics: R, D and E. These measure transmissibility at te global scale (e.g., countrywide).

- R is the commonly computed effective reproduction number
- D is the mean of the local reproduction numbers
- E is the contraharmonic mean of the local reproduction numbers

D and E are derived as solutions to optimal experimental design problems as described in the paper. E is especially our focus as it offers improved sensitivity to resurgence events and has important properties that might make it helpful for informing decision-making. We propose E as an alternative to R when the dynamics of local groups are heterogeneous.

System Requirements
- Should work with Matlab v2020a and above. Tested on macOS v11.6.4.
- The implementation in R has been tested on R version 3.6.2
- Slight dependence on the linespecer package (license included in main folder).
- All data for case study as simple csv files. These represent Israel COVID-19 data.
- All local reproduction number estimates use the EpiFilter package: https://github.com/kpzoo/EpiFilter.

Instructions and installation
- Should work with any standard Matlab installation and generate figures/results.
- Run FigX.m where X is the figure in the manuscript to be reproduced.
- File Fig3and4.m can be easily modified for different epidemic simulations over many groups.
- Code is self contained so no external installations required.
- Run times of all scripts are of the order of minutes or faster.

