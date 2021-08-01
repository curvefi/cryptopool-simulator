import numpy as np
from math import log10
import json
import itertools
from copy import copy


X = np.logspace(log10(100), log10(1e4), 32)
Xname = "A"
Y = np.logspace(log10(1e-9), log10(1e-5), 64)
Yname = "gamma"

other_params = dict(
    D=10e6,
    adjustment_step=.55e-5,
    fee_gamma=7e-6,
    ma_half_time=600,
    mid_fee=1.1e-5,
    out_fee=3.8e-4,
    n=2,
    log=0,
    price_threshold=.55e-5,
    gamma=6e-7,
    ext_fee=2e-4,
    A=4500.0)

config = {
    'configuration': [],
    'datafile': [
        'eurusdt'],
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
