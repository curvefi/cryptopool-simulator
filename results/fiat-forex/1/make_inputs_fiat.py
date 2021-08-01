import numpy as np
from math import log10
import json
import itertools
from copy import copy


# X = np.logspace(log10(0.25), log10(1.18), 16)
X = np.logspace(log10(0.1e-6), log10(1e-4), 32)
Xname = "mid_fee"
# Y = np.logspace(log10(9e-6), log10(0.0003), 16)
Y = np.logspace(log10(5e-5), log10(2e-3), 32)
Yname = "out_fee"

other_params = dict(
    D=10e6,
    adjustment_step=0.0005,
    fee_gamma=2e-8,
    ma_half_time=600,
    mid_fee=5e-4,
    out_fee=5e-4,
    n=2,
    log=0,
    price_threshold=0.0005,
    gamma=1.4e-8,
    ext_fee=0.0,
    A=100.0)

config = {
    'configuration': [],
    'datafile': [
        'eurusd'],
    'debug': 0}

for x, y in itertools.product(X, Y):
    params = copy(other_params)
    params[Xname] = x
    params[Yname] = y
    if Xname == 'mid_fee':
        other_params['price_threshold'] = x / 2
        other_params['adjustment_step'] = x / 2
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
