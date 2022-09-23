# resampleGLM

Methods for fitting generalised linear models using resampling to reduce Type 1 error in R.

A demo R script - example/mock_regression.r - is provided to demonstrate dCRT and bootstrap methods. It works on a mock dataset which has an unobserved independent variable which causes problems for a negative binomial model. qqplot.pdf shows that the resampling process reduces the type I error observed with a normal GLM fit to within expected uncertainties. 

---

Three resampling methods are currently implemented.

### distilled CRT

Resamples the target independent variable as a function of other covariates (as a linear or logistic model). z-scores are calculated with the score test method for each resample to generate a null distribution. The z-score of the original target variable relative to the null distribution is used to evaluate a p-value.

This is the method used by https://github.com/Katsevich-Lab/sceptre which provided the inspiration for this code and which several aspects of the code are based on.

### Permutation

When the target independent variable is independent of all covariates, resampling without replacement is performed to generate the null fits. The z-score analysis then proceeds as above.

### Bootstrap

Resample with replacement to determine the uncertainty distribution around the fitted parameters.

-----

For all resampling methods, a normal distribution is currently used as the null distribution.
