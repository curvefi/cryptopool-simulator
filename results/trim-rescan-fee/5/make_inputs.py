import numpy as np
from math import log10
import json
import itertools
from copy import copy


# X = np.logspace(log10(0.25), log10(1.18), 16)
X = np.logspace(log10(0.217), log10(0.553), 16)
Xname = "A"
# Y = np.logspace(log10(9e-6), log10(0.0003), 16)
Y = np.logspace(log10(1.64e-4), log10(4.44e-3), 16)
Yname = "gamma"

other_params = dict(
    D=3e8,
    adjustment_step=0.0015,
    fee_gamma=0.00125,
    ma_half_time=600,
    mid_fee=7.5e-4,
    out_fee=5.6e-3,
    n=3,
    log=0,
    price_threshold=0.0028,
    gamma=0.00117,
    ext_fee=3e-4,
    A=0.295)

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
