#!/usr/bin/env python3

import numpy as np
from math import log10
import json
import itertools
from copy import copy


X = np.logspace(log10(1e-6), log10(1), 32)
Xname = "fee_gamma"
Y = np.logspace(log10(1e-6), log10(0.199), 32)
Yname = "gamma"

other_params = dict(
    D=1,
    adjustment_step=1e-7,
    fee_gamma=0.006,
    ma_half_time=600,
    mid_fee=0.002,
    out_fee=0.01,
    gas_fee=0,
    n=3,
    log=0,
    allowed_extra_profit=1e-10,
    ext_fee=0.0003,
    gamma=1e-5,
    boost_rate=0.0,
    A=0.5)

config = {
    'configuration': [],
    'datafile': [
        "cryptogold-eth_btc_1min",
        "cryptogold-xau_btc_1min",
        "cryptogold-xau_eth_1min"],
    'debug': 0}

for x, y in itertools.product(X, Y):
    params = copy(other_params)
    params[Xname] = x
    params[Yname] = y
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
