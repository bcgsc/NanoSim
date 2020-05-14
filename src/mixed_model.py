#!/usr/bin/env python

"""
Created on May 17th, 2018 by Chen Yang
This script defines Poisson-Geometric distribution and Weibull-Geometric distribution
"""

import numpy as np
from math import ceil, exp
from scipy.stats import rv_discrete, poisson, geom, lognorm


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


class trunc_lognorm_gen(rv_discrete):
    def _cdf(self, x, s, m):
        lncdf_x = lognorm.cdf(x, s, scale=m)
        lncdf_a = lognorm.cdf(self.a, s, scale=m)
        lncdf_b = lognorm.cdf(self.b, s, scale=m)
        return (lncdf_x - lncdf_a) / (lncdf_b - lncdf_a)

    def _ppf(self, q, s, m):
        lncdf_a = lognorm.cdf(self.a, s, scale=m)
        lncdf_b = lognorm.cdf(self.b, s, scale=m)
        ln_q = q * (lncdf_b - lncdf_a) + lncdf_a
        return lognorm.ppf(ln_q, s, scale=m)


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
        value = int(round(ceil(lam * np.random.weibull(k))))
    else:
        value = np.random.geometric(prob) - 1

    if value == 0:
        value = 1

    return value


def trunc_lognorm_rvs(error_type, read_type, basecaller, n):
    if basecaller == "albacore":
        a = 1
        b = 28
        if read_type == "DNA":
            if error_type == "match":
                mean = 2.7418286
                sd = 0.7578693
            elif error_type == "mis":
                mean = 1.6597215
                sd = 0.6814804
            elif error_type == "ins":
                mean = 1.9016147
                sd = 0.6842999
            elif error_type == "ht":
                mean = 2.3739153
                sd = 0.9635895
            else:  # unaligned
                mean = 2.5484921
                sd = 0.7742894
        elif read_type == "dRNA":
            if error_type == "match":
                mean = 2.236641
                sd = 0.434045
            elif error_type == "mis":
                mean = 1.8138169
                sd = 0.4535039
            elif error_type == "ins":
                mean = 1.9322685
                sd = 0.4668444
            elif error_type == "ht":
                mean = 2.0166876
                sd = 0.5714308
            else:  # unaligned
                mean = 2.1371272
                sd = 0.4763441
        else:  # cDNA
            if error_type == "match":
                mean = 2.6003978
                sd = 0.7181057
            elif error_type == "mis":
                mean = 1.6380338
                sd = 0.6695235
            elif error_type == "ins":
                mean = 1.8462438
                sd = 0.6661691
            elif error_type == "ht":
                mean = 2.510699
                sd = 1.082626
            else:  # unaligned
                mean = 2.6004634
                sd = 0.8526468
    elif basecaller == "guppy":
        a = 1
        b = 31
        if read_type == "DNA":
            if error_type == "match":
                mean = 2.9863022
                sd = 0.9493498
            elif error_type == "mis":
                mean = 1.6184245
                sd = 0.7585733
            elif error_type == "ins":
                mean = 1.8852560
                sd = 0.7623103
            elif error_type == "ht":
                mean = 1.995397
                sd = 1.008650
            else:  # unaligned
                mean = 1.2626728
                sd = 0.9012829
        elif read_type == "dRNA":
            if error_type == "match":
                mean = 2.236641
                sd = 0.434045
            elif error_type == "mis":
                mean = 1.8138169
                sd = 0.4535039
            elif error_type == "ins":
                mean = 1.9322685
                sd = 0.4668444
            elif error_type == "ht":
                mean = 2.0166876
                sd = 0.5714308
            else:  # unaligned
                mean = 2.1371272
                sd = 0.4763441
        else:  # cDNA
            if error_type == "match":
                mean = 2.7500148
                sd = 0.9195383
            elif error_type == "mis":
                mean = 1.5543628
                sd = 0.7601223
            elif error_type == "ins":
                mean = 1.765634
                sd = 0.777587
            elif error_type == "ht":
                mean = 2.001173
                sd = 1.008647
            else:  # unaligned
                mean = 1.2635415
                sd = 0.9008419

    m = exp(mean)
    truncln_gen = trunc_lognorm_gen(name="truncln", a=a, b=b, shapes="s, m")
    return truncln_gen.rvs(s=sd, m=m, size=n)
