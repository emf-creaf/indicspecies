-------------------------------
 NEWS for R Package "indicspecies"
-------------------------------

# Version 1.8.0
- Negative associations can now be selected in multipatt, when func = "r" or "r.g"

# Version 1.7.15
- Update of indicator species analysis vignette (negative associations)

# Version 1.7.14
- Bug correction in pruneindicators and print.indicators
- Bug correction in nichecentroid and plotniche for distance matrix with equal resources
- Translation of documentation to Roxygen

# Version 1.7.13
- Removed dependency from packages 'sp' and 'rgeos'

# Version 1.7.12
- Fixing error introduced unadvertedly

# Version 1.7.11
- Code modification in 'multipatt' to avoid large matrix multiplication

# Version 1.7.10
- Bug correction in 'multipatt' for single-valued restcomb
- Vignettes rewritten in Rmarkdown

# Version 1.7.9
- Improvements in 'multipatt', 'signassoc' and 'indicators' functions to accept custom permutation designs, by Noah Dell.

# Version 1.7.8
- New algorithm for 'combinespecies', similar to the one used for 'indicators'
- New option 'min.order' in function 'combinespecies'
- Permutation tests are now available in function 'indicators'
- Options 'nboot' and 'alpha' of functions 'strassoc' and 'indicators' renamed to 'nboot.ci' and 'alpha.ci'

# Version 1.7.7
- New function 'predict.indicators.cv' for cross validation of predictive value
- New option 'min.order' in function 'indicators'
