import numpy as np
from math import log10
import json
from copy import copy


X = np.logspace(log10(1e-4), log10(1e-1), 96)
param = 'mid_fee'


other_params = dict(
    D=1e8,
    adjustment_step=0.0015,
    fee_gamma=0.0052,
    ma_half_time=400,
    mid_fee=4e-4,
    out_fee=0.016,
    n=3,
    log=0,
    price_threshold=0.0028,
    gamma=6e-5,
    ext_fee=5e-4,
    A=135)
    # gamma=1.1e-3,
    # A=5)

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
