# fitMaxnet

An ensemble of functions to fit Environmental Niche Models (ENMs) using the *maxnet* implementation of the *MaxEnt* method.

## Introduction

The package *fitMaxnet* was developed to implement efficient processes for fitting Ecological Niche Models (ENMs, also often called Species Distribution Models, SDMs) using the _**R**_-package *maxnet*. ENMs are designed to describe the relationship between occurrence and environmental variables so that locations without observed occurrences can be given a score indicating the likelihood that conditions are suitable for individuals of a given species to occur there.

*maxnet* is an open-source implementation of the original Java application *MaxEnt* which fitted ENMs using the machine learning modelling method called "Maximum Entropy". Technical information on the application of the MaxEnt method to modelling presence-only species occurrence data is found in Phillips et al. (2006, 2008). For a description of the workings of MaxEnt focussed on practical applications for ecology see Elith et al. (2011).

Package *maxnet* uses the _**R**_-package *glmnet* to fit models. This is possible because it was shown that the MaxEnt approach applied to presence-only data is equivalent mathematically to fitting linear models using the *glment* method (Fithian and Hastie 2013). As a form of spatial point pattern, species presence records can also be modelled as inhomogenous Poisson spatial point patterns (Renner et al. 2013, 2015) yielding models equivalent to those produced by *glmnet*.

## What does *fitMaxnet* provide?

Functions within *fitMaxnet* are designed to provide efficient processing for a selection of critical steps in the production of ENMs. There many novel and potentially useful developments at each of the major steps towards fitting an ENM. *fitMaxnet* does not include many of them. Instead, I have implemented a select few options which personal experience indicates are robust, reliable and can be implemented with computational efficiency. This last point is significant because a major motivation for developing this package was to produce reasonable models for large assemblages of species in reasonable time.

However, there are a few experimental functions included in the package. These included: *boyce_index*, *computeMESS*, *corrAnalysis*, *energyStats*, and *pocplot*.

Functions are provided to:

* prepare covariate data into an efficient data structure for model fitting (*prepData*, *removeDuplicates*, *sampleBackground*) and projection (*prepProjData*)
* constrain background point selection (*bufferPoints*)
* account for spatial sampling bias in occurrence data using optimised spatial thinning (*occThin*)
* implement k-fold cross-validation (*kFoldSet*)
* fit models (*fit_maxnet*)
* project models (*projectMaxnet*)
* mask regions of model extrapolation beyond the range of environmental predictors used in model fitting (*maskExtrapolation*), and
* gather and present model performance and variable importance information (*feature_summary*, *varImportance*, *variable_interaction_heatmap*).
 
Unlike the original *MaxEnt* Java implementation, and most other wrappers around _**R**_-package *maxnet*, the functions in this package can work with projected predictor variables (e.g. equal area projections).


## Installing *fitMaxnet*

Before installing the package, you should ensure that the following _**R**_ packages are installed in your local _**R**_ environment:

* maxnet
* ggplot2
* dplyr
* terra
* geodist
* sf
* ggcorrplot
* energy 
* ggpubr
* dismo
* PRROC
* Rfast
* ggsci
* stringr

To install the package, please first install the _**R**_ package *remotes*.

Then issue the following command at the _**R**_ command prompt:

> remotes::install_github("peterbat1/fitMaxnet")

## Potential workflows using *fitMaxnet*

To be completed...

## References

Elith, J., S. J. Phillips, T. Hastie, M. Dudík, Y. E. Chee, and C. J. Yates. 2011. A statistical explanation of MaxEnt for ecologists. Diversity and Distributions 17:43–57.

Fithian, W., and T. Hastie. 2013. Finite-sample equivalence in statistical models for presence-only data. The Annals of Applied Statistics 7:1917–1939.

Phillips, S. J., R. P. Anderson, and R. E. Schapire. 2006. Maximum entropy modeling of species geographic distributions. Ecological Modelling 190:231–259.

Phillips, S. J., and M. Dudík. 2008. Modeling of species distributions with Maxent: new extensions and a comprehensive evaluation. Ecography 31:161–175.

Renner, I. W., J. Elith, A. Baddeley, W. Fithian, T. Hastie, S. J. Phillips, G. Popovic, and D. I. Warton. 2015. Point process models for presence-only analysis. Methods in Ecology and Evolution 6:366–379.

Renner, I. W., and D. I. Warton. 2013. Equivalence of MAXENT and Poisson Point Process Models for Species Distribution Modeling in Ecology: Equivalence of MAXENT and Poisson Point Process Models. Biometrics 69:274–281.


