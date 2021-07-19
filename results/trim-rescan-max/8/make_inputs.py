import numpy as np
from math import log10
import json
import itertools
from copy import copy


X = np.logspace(log10(1e-3), log10(8e-2), 32)
Xname = "out_fee"
Y = np.logspace(log10(1e-5), log10(1e-1), 32)
Yname = "fee_gamma"

other_params = dict(
    D=3e8,
    adjustment_step=0.00049,
    fee_gamma=2.02e-3,
    ma_half_time=600,
    mid_fee=1.25e-3,
    out_fee=8e-3,
    n=3,
    log=0,
    price_threshold=0.00049,
    gamma=1.545e-4,
    ext_fee=3e-4,
    A=0.6561)

config = {
    'configuration': [],
    'datafile': [
        'btcusdt',
        'ethusdt',
        'ethbtc'],
    'debug': 0}

for x, y in itertools.product(X, Y):
    params = copy(other_params)
    params[Xname] = x
    params[Yname] = y
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
