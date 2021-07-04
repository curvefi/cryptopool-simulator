import numpy as np
from math import log10
import json
from copy import copy


X = np.logspace(log10(5e-6), log10(5e-4), 96)
param = 'gamma'


other_params = dict(
    D=5e8,
    adjustment_step=0.0015,
    fee_gamma=0.01,
    ma_half_time=600,
    mid_fee=1e-3,
    out_fee=1e-3,
    n=3,
    log=0,
    price_threshold=0.0028,
    gamma=0.0118,
    ext_fee=0,
    A=0.32)

config = {
    'configuration': [],
    'datafile': [
        'btcusdt',
        'ethusdt',
        'ethbtc'],
    'debug': 0}

for x in X:
    params = copy(other_params)
    params[param] = x
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
