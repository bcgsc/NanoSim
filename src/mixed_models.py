#!/usr/bin/env python
"""
Created on Apr 28th by Chen Yang

This script is used to generate random numbers following certain mixed distribution models
"""

import numpy as np
import math

# numpy.random.geometric generate positive integers, starting from 1
# the rgeom in R generate values starting from 0


def pois_geom(lam, prob, weight):
    tmp_rand = np.random.random()
    if tmp_rand < weight:
        value = np.random.poisson(lam) + 1
    else:
        value = np.random.geometric(prob)
    return value


def wei_geom(lam, k, prob, weight):
    tmp_rand = np.random.random()
    if tmp_rand < weight:
        value = int(round(math.ceil(lam * np.random.weibull(k))))
    else:
        value = np.random.geometric(prob) - 1

    if value == 0:
        value = 1

    return value