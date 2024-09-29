#!/usr/bin/env python3

import numpy as np
from math import log10
import json
from copy import copy


X = np.logspace(log10(1e-4), log10(0.05), 128)
Xname = "mid_fee"

other_params = dict(
    D=40e3,
    adjustment_step=1e-7,
    fee_gamma=0.00022,
    ma_half_time=600,
    mid_fee=0.0004,
    out_fee=0.002,
    gas_fee=0,
    n=2,
    log=0,
    allowed_extra_profit=1e-10,
    ext_fee=0.0001,
    gamma=0.0022,
    boost_rate=0.1,
    A=40)

config = {
    'configuration': [],
    'datafile': ["eurusd-1m"],
    'debug': 0}

for x in X:
    params = copy(other_params)
    params[Xname] = x
    params["out_fee"] = x
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
