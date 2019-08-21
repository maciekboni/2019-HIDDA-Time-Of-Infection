# 2019-HIDDA-Time-Of-Infection

[[ WORK IN PROGRESS ]]

C++ code for inferring the time of infection from serological data, as presented in 

   Inferring the time of infection from serological data.
   Boni MF, Mølbak K, Krogfelt KA.
   In L. Held, N. Hens, P O'Neill, J. Wallinga, editors, Handbook of Infectious Disease Data Analysis, 2019.

STEP 1:

Compile the code and run

`./mlesearch_and_profile`
   
from the command line.  This will give you the MLE estimates for the slope (antibody waning rate) as well as the estimates of mean and variance of peak antibody titer.  Note, that the search is not inferring the time of infection at this moment.

If you run

   `./mlesearch_and_profile -mlesearch_and_profile`
   
the slope parameter will be profiled and a confidence interval will be reported.
