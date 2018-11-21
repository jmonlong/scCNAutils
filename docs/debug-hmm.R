library(devtools)
library(dplyr)
library(ggplot2)

load_all()

## One sample at a time
P = 100
SD = .1
ge = tibble(chr=1, start=1:P, end=1:P, samp1=rnorm(P, 0, SD))
## Add gain
L = 10
ge$samp1[50:(50+L-1)] = rnorm(L, log(1.3), SD)
## Add another chr
ge = rbind(ge, tibble(chr=2, start=1:P, end=1:P, samp1=rnorm(P, log(.8), SD)))

res = cnaHMM(ge, .01)
ggplot(res, aes(x=start, y=mean, color=CN)) + geom_point() + theme_bw() +
  facet_grid(.~chr) + scale_colour_brewer(palette='Set1')


## Multiple samples at a time
N = 10
P = 100
SD = .1
mat = matrix(rnorm(N*P, 0, SD), P)
## Consistent gain
L = 10
mat[50:(50+L-1),] = matrix(rnorm(L*N, log(1.3), SD), L)
## Loss in a few samples only
NN = 4
mat[10:(10+L-1),1:NN] = matrix(rnorm(L*NN, log(.5), SD), L)
## 
colnames(mat) = paste0('samp', 1:N)
ge = cbind(tibble(chr=1, start=1:P, end=1:P), mat)
## Add another chr
mat = matrix(rnorm(N*P, log(.8), SD), P)
colnames(mat) = paste0('samp', 1:N)
ge = rbind(ge, cbind(tibble(chr=2, start=1:P, end=1:P), mat))

res = cnaHMM(ge, .01)
ggplot(res, aes(x=start, y=mean, color=CN)) + geom_point() + theme_bw() +
  facet_grid(.~chr) + scale_colour_brewer(palette='Set1')


## It seems to work fine on simulated data
