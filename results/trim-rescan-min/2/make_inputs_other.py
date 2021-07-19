import numpy as np
from math import log10
import json
from copy import copy


X = np.logspace(log10(5e-4), log10(5e-2), 256)
param = 'gamma'


other_params = dict(
    D=3e8,
    adjustment_step=0.00049,
    fee_gamma=5e-4,
    ma_half_time=600,
    mid_fee=1.2e-3,
    out_fee=4.5e-3,
    n=3,
    log=0,
    price_threshold=0.00049,
    gamma=8e-3,
    ext_fee=3e-4,
    A=0.2)

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
