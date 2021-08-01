import numpy as np
from math import log10
import json
from copy import copy


X = np.logspace(log10(1e-8), log10(1e-1), 192)
param = 'gamma'

other_params = dict(
    D=10e6,
    adjustment_step=0.55e-5,
    fee_gamma=7e-6,
    ma_half_time=600,
    mid_fee=1.1e-5,
    out_fee=3.8e-4,
    n=2,
    log=0,
    price_threshold=0.55e-5,
    gamma=6e-7,
    ext_fee=2e-4,
    A=5000.0)

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
