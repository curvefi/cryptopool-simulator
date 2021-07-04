import numpy as np
from math import log10
import json
import itertools
from copy import copy


X = np.logspace(log10(2e-4), log10(1e-2), 16)
Xname = "mid_fee"
Y = np.logspace(log10(1e-3), log10(10e-2), 16)
Yname = "out_fee"

other_params = dict(
    D=3e8,
    adjustment_step=0.0015,
    fee_gamma=1.25e-3,
    ma_half_time=600,
    mid_fee=4e-4,
    out_fee=4e-3,
    n=3,
    log=0,
    price_threshold=0.0028,
    gamma=0.002,
    ext_fee=3e-4,
    A=0.254)

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
