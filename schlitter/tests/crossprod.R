#! /usr/bin/R
## example cross product
## input matrix
tmp = matrix(NA, nrow = 6, ncol = 4);
tmp[1, ] = c(10, 11, 12, 13);
tmp[2, ] = c(14, 15, 13, 17);
tmp[3, ] = c(15, 16, 16, 18);
tmp[4, ] = c(13, 12, 12, 19);
tmp[5, ] = c(12, 13, 11, 17);
tmp[6, ] = c(10, 11, 10, 14);

## scale columns (just for comparison)
tmp.c = scale(tmp, TRUE, FALSE);

## matrix cross product
crossprod(tmp.c);

## covariance
cov(tmp.c);

