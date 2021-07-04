import numpy as np
from math import log10
import json
import itertools
from copy import copy


X = np.logspace(log10(0.1), log10(1000), 32)
Xname = "A"
Y = np.logspace(log10(1e-7), log10(0.1), 32)
Yname = "gamma"

other_params = dict(
    D=5e8,
    adjustment_step=0.0015,
    fee_gamma=0.01,
    ma_half_time=600,
    mid_fee=4e-4,
    out_fee=4e-3,
    n=3,
    log=0,
    price_threshold=0.0028,
    gamma=2.5e-4,
    ext_fee=5e-4,
    A=5)

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
