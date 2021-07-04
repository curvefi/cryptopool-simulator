import numpy as np
from math import log10
import json
from copy import copy


X = np.logspace(log10(1e-4), log10(1e-2), 96)
param = 'gamma'


other_params = dict(
    D=3e8,
    adjustment_step=0.0015,
    fee_gamma=0.01,
    ma_half_time=600,
    mid_fee=7.5e-4,
    out_fee=0.03,
    n=3,
    log=0,
    price_threshold=0.0028,
    gamma=1e-3,
    ext_fee=3e-4,
    A=0.269)

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
