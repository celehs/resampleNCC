# Resampling Procedures for Making Inference under Nested Case-control Studies

## Overview 

This repository contains [R codes](R/FUN.R) that were used in the paper by Cai & Zheng (2013). The `FNCC.WGT.FUN` and `PtbNCC.WGT.FUN` functions can be used to calculate and perturb the sampling weights with possible matching. With these (perturbed) sampling weights, a large number of inference procedures can be easily developed for analyzing data from nested case-control studies. For example, one can use the `NCC.Cox.FUN.Regularize` function to estimate the regression coefficients from the Cox proportional hazards model. 

Click [HERE](example.R) to see a simulated data example.  

**NOTE**: An R package is under development and will be released in the near future. 

## References

T Cai, Y Zheng. Resampling Procedures for Making Inference under Nested Case-control Studies. _JASA_, 2013, 108(504):1532-1544. <https://doi.org/10.1080/01621459.2013.856715>
