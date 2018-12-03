#!/usr/bin/env python

"""
Created on May 17th, 2018 by Chen Yang
This script defines Poisson-Geometric distribution and Weibull-Geometric distribution
"""

import numpy as np
import math
from scipy.stats import rv_discrete, poisson, geom


# Scipy geometric starts with x = 1


class poisgeom_gen(rv_discrete):
    # Poisson-Geometric distribution
    def _pmf(self, x, l, p, w):
        return w * poisson.pmf(x, l) + (1 - w) * geom.pmf(x, p, loc=-1)


class weigeom_gen(rv_discrete):
    # Weibull-Geometric distribution, Geometric start from 0
    def _cdf(self, x, l, k, p, w):
        wei_cdf = 1 - np.exp(-1 * np.power(x / l, k))
        return w * wei_cdf + (1 - w) * geom.cdf(x, p, loc=-1)

    def _pmf(self, x, l, k, p, w):
        return self.cdf(x, l, k, p, w) - self.cdf(x-1, l, k, p, w)


class weigeom2_gen(rv_discrete):
    # Weibull-Geometric distribution, Geometric start from 1
    def _cdf(self, x, l, k, p, w):
        wei_cdf = 1 - np.exp(-1 * np.power(x / l, k))
        return w * wei_cdf + (1 - w) * geom.cdf(x, p)

    def _pmf(self, x, l, k, p, w):
        return self.cdf(x, l, k, p, w) - self.cdf(x-1, l, k, p, w)


def pois_geom(lam, prob, weight):
    # Draw a random number from Poisson-Geometric distribution
    # Faster to use numpy random than using Scipy rvs
    tmp_rand = np.random.random()
    if tmp_rand < weight:
        value = np.random.poisson(lam) + 1
    else:
        value = np.random.geometric(prob)
    return value


def wei_geom(lam, k, prob, weight):
    # Draw a random number from Weibull-Geometric distribution
    tmp_rand = np.random.random()
    if tmp_rand < weight:
        value = int(round(math.ceil(lam * np.random.weibull(k))))
    else:
        value = np.random.geometric(prob) - 1

    if value == 0:
        value = 1

    return value