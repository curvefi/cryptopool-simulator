import numpy as np
from math import log10
import json
from copy import copy


X = np.logspace(log10(1e-10), log10(1e-6), 192)
param = 'gamma'

other_params = dict(
    D=10e6,
    adjustment_step=1.5e-5,
    fee_gamma=2e-8,
    ma_half_time=600,
    mid_fee=3e-5,
    out_fee=6e-5,
    n=2,
    log=0,
    price_threshold=1.5e-5,
    gamma=1.4e-8,
    ext_fee=0.0,
    A=7910.0)

config = {
    'configuration': [],
    'datafile': ['eurusdt'],
    'debug': 0}

for x in X:
    params = copy(other_params)
    params[param] = x
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
