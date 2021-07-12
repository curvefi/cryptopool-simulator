import numpy as np
from math import log10
import json
import itertools
from copy import copy


# X = np.logspace(log10(0.25), log10(1.18), 16)
X = np.logspace(log10(2e-4), log10(1e-2), 32)
# X = np.logspace(log10(1), log10(20), 32)
Xname = "mid_fee"
# Y = np.logspace(log10(9e-6), log10(0.0003), 16)
Y = np.logspace(log10(2e-4), log10(5e-2), 32)
# Y = np.logspace(log10(5e-6), log10(2e-4), 32)
Yname = "out_fee"

other_params = dict(
    D=3e8,
    adjustment_step=0.00049,
    fee_gamma=1.53e-3,
    ma_half_time=600,
    mid_fee=8e-4,
    out_fee=7.6e-3,
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
