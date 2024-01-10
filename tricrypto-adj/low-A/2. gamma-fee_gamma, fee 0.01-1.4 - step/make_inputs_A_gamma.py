#!/usr/bin/env python3

import numpy as np
from math import log10
import json
import itertools
from copy import copy


X = np.logspace(log10(4e-6), log10(4e-2), 32)
Xname = "fee_gamma"
Y = np.logspace(log10(1e-5), log10(4e-4), 32)
Yname = "gamma"

other_params = dict(
    D=40e6,
    adjustment_step=5e-4,
    fee_gamma=0.006,
    ma_half_time=600,
    mid_fee=0.0001,
    out_fee=0.014,
    gas_fee=10,
    n=3,
    log=0,
    allowed_extra_profit=1e-10,
    ext_fee=0.0003,
    gamma=1e-5,
    A=2)

config = {
    'configuration': [],
    'datafile': [
        "btcusdt-1m",
        "ethusdt-1m",
        "ethbtc-1m"],
    'debug': 0}

for x, y in itertools.product(X, Y):
    params = copy(other_params)
    params[Xname] = x
    params[Yname] = y
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
