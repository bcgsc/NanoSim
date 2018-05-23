#!/usr/bin/env python

"""
Created on May 17th, 2018
@author: Chen Yang
This script fits errors to statistical models
"""


from __future__ import print_function
from __future__ import with_statement
import sys
import numpy as np
import multiprocessing as mp
from time import strftime
from math import sqrt, ceil
from scipy.optimize import minimize
from itertools import chain
from mixed_model import *

poisgeom = poisgeom_gen(name="poisgeom")
weigeom = weigeom_gen(name="weigeom", a=1)


def read_histogram(f, error):
    hist = {}
    with open(f, 'r') as f:
        f.readline()
        for line in f:
            info = line.strip().split()
            hist[int(info[0])] = int(info[1])
    if error == 'mis':
        os = [[key - 1] * value for key, value in hist.items()]
        os = list(chain.from_iterable(os))
        pmf, bins = np.histogram(os, bins=max(os), density=True)
    elif error == 'indel':
        os = [[key] * value for key, value in hist.items()]
        os = list(chain.from_iterable(os))
        pmf, bins = np.histogram(os, bins=max(os) - 1, density=True)

    cdf = np.cumsum(pmf)

    return os, cdf


def mis_ll(params):
    global mis_os, mis_cdf
    # negative log likelihood of all observations
    l = params[0]
    p = params[1]
    w = params[2]

    diff = max(abs(poisgeom.cdf(range(0, len(mis_cdf)), l, p, w) - mis_cdf))
    return diff


def mis_fit(init):
    global mis_cdf
    # bnds = ((0.001, None), (0.001, 1), (0.001, 0.999))
    result = minimize(mis_ll, init, method='Nelder-Mead')
    diff = max(abs(poisgeom.cdf(range(0, len(mis_cdf)), result.x[0], result.x[1], result.x[2]) -
                   mis_cdf))
    return init, result.x, diff


def ins_ll(params):
    global ins_cdf
    # negative log likelihood of all observations
    l = params[0]
    k = params[1]
    p = params[2]
    w = params[3]

    diff = max(abs(weigeom.cdf(range(1, len(ins_cdf) + 1), l, k, p, w) - ins_cdf))
    return diff


def ins_fit(init):
    global ins_cdf
    result = minimize(ins_ll, init, method='Nelder-Mead')
    diff = max(abs(weigeom.cdf(range(1, len(ins_cdf) + 1), result.x[0], result.x[1], result.x[2], result.x[3]) -
                   ins_cdf))
    return init, result.x, diff


def del_ll(params):
    global del_cdf
    # negative log likelihood of all observations
    l = params[0]
    k = params[1]
    p = params[2]
    w = params[3]

    diff = max(abs(weigeom.cdf(range(1, len(del_cdf) + 1), l, k, p, w) - del_cdf))
    return diff


def del_fit(init):
    global del_cdf
    result = minimize(del_ll, init, method='Nelder-Mead')
    diff = max(abs(weigeom.cdf(range(1, len(del_cdf) + 1), result.x[0], result.x[1], result.x[2], result.x[3]) -
                   del_cdf))
    return init, result.x, diff


def model_fitting(prefix, threads):
    out = open(prefix + '_model_profile', 'w')
    out.write("Type\tlambda\tk\tprob\tweight\n")
    global mis_os, mis_cdf, ins_os, ins_cdf, del_os, del_cdf
    mis_os, mis_cdf = read_histogram(prefix + '_mis.hist', 'mis')
    ins_os, ins_cdf = read_histogram(prefix + '_ins.hist', 'indel')
    del_os, del_cdf = read_histogram(prefix + '_del.hist', 'indel')

    # Fit mismatches to Poisson-Geometric distribution
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Mismatch fitting start\n")
    pool = mp.Pool(threads)
    TASKS = [(l, p, w) for l in np.arange(0.1, 0.9, 0.1) for p in np.arange(0.1, 0.9, 0.1)
             for w in np.arange(0.1, 0.9, 0.1)]
    apply_async_res = [pool.apply_async(mis_fit, (t,)) for t in TASKS]
    pool.close()
    pool.join()

    results = [x.get() for x in apply_async_res]
    results.sort(key=lambda x: x[2])
    precision = 1.36/sqrt(len(mis_os))
    for res in results:
        params = res[1]
        if params[0] <= 0 or params[1] <= 0 or params[1] >= 1 or params[2] <= 0 or params[2] >= 1:
            continue
        else:
            diff = res[2]
            mis_params = params
            out.write("mismatch\t" + str(mis_params[0]) + '\t0\t' + str(mis_params[1]) + '\t' +
                      str(mis_params[2]) + '\n')
            if diff <= precision:
                print("Mismatch parameters: ", mis_params, diff)
                break
            else:
                print("WARNING! Mismatch parameters may not be optimal!\n", params, diff)
                break
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Mismatch fitting done\n")

    # Fit insertions to Weibull-Geometric distribution
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Insertion fitting start\n")
    pool = mp.Pool(threads)
    TASKS = [(l, p, k, w) for l in np.arange(0.1, 1.3, 0.1) for p in np.arange(0.1, 1.3, 0.1)
             for k in np.arange(0.1, 0.9, 0.1) for w in np.arange(0.1, 0.9, 0.1)]
    apply_async_res = [pool.apply_async(ins_fit, (t,)) for t in TASKS]
    pool.close()
    pool.join()

    precision = 1.36/sqrt(len(ins_os))
    results = [x.get() for x in apply_async_res]
    results.sort(key=lambda x: x[2])
    for res in results:
        params = res[1]
        if params[0] <= 0 or params[1] <= 0 or \
            params[2] <= 0 or params[2] >= 1 or params[3] <= 0 or params[3] >= 1:
            continue
        else:
            diff = res[2]
            ins_params = params
            out.write("insertion\t" + str(ins_params[0]) + '\t' + str(ins_params[1]) + '\t' +
                      str(ins_params[2]) + '\t' + str(ins_params[3]) + '\n')
            if diff <= precision:
                print("Insertion parameters: ", ins_params, diff)
                break
            else:
                print("WARNING! Insertion parameters may not be optimal!\n", params, diff)
                break
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Insertion fitting done\n")

    # Fit deletions to Weibull-Geometric distribution
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Deletion fitting start\n")
    pool = mp.Pool(threads)
    TASKS = [(l, p, k, w) for l in np.arange(0.1, 1.3, 0.1) for p in np.arange(0.1, 1.3, 0.1)
             for k in np.arange(0.1, 0.9, 0.1) for w in np.arange(0.1, 0.9, 0.1)]
    apply_async_res = [pool.apply_async(del_fit, (t,)) for t in TASKS]

    pool.close()
    pool.join()

    precision = 1.36/sqrt(len(del_os))
    results = [x.get() for x in apply_async_res]
    results.sort(key=lambda x: x[2])
    for res in results:
        params = res[1]
        if params[0] <= 0 or params[1] <= 0 or \
           params[2] <= 0 or params[2] >= 1 or params[3] <= 0 or params[3] >= 1:
            continue
        else:
            diff = res[2]
            del_params = params
            out.write("deletion\t" + str(del_params[0]) + '\t' + str(del_params[1]) + '\t' +
                      str(del_params[2]) + '\t' + str(del_params[3]) + '\n')
            if diff <= precision:
                print("Deletion parameters: ", del_params, diff)
                break
            else:
                print("WARNING! Deletion parameters may not be optimal!\n", params, diff)
                break
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Deletion fitting done\n")

    out.close()
