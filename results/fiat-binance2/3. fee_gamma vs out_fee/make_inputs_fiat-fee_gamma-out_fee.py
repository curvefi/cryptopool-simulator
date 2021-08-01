import numpy as np
from math import log10
import json
import itertools
from copy import copy


X = np.logspace(log10(1e-5), log10(1e-1), 32)
Xname = "fee_gamma"
Y = np.logspace(log10(1e-4), log10(2e-3), 32)
Yname = "out_fee"

other_params = dict(
    D=1e6,
    adjustment_step=0.55e-5,
    allowed_extra_profit=1e-8,
    fee_gamma=8e-4,
    ma_half_time=600,
    mid_fee=1e-4,
    out_fee=1e-4,
    n=2,
    log=0,
    gamma=4.869e-5,
    ext_fee=0.0,
    A=500.0)

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
