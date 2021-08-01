import numpy as np
from math import log10
import json
import itertools
from copy import copy


X = np.logspace(log10(1e-6), log10(4e-4), 32)
Xname = "mid_fee"
Y = np.logspace(log10(5e-5), log10(2e-3), 32)
Yname = "out_fee"

other_params = dict(
    D=10e6,
    adjustment_step=.55e-5,
    fee_gamma=1.73e-5,
    ma_half_time=600,
    mid_fee=1.1e-5,
    out_fee=3.8e-4,
    n=2,
    log=0,
    price_threshold=.55e-5,
    gamma=4.3e-7,
    ext_fee=0,  # 2e-4,
    A=210e3)

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
