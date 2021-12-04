import numpy as np
from math import log10
import json
from copy import copy


X = np.logspace(log10(1e-7), log10(1e-2), 64)
param = 'gamma'

other_params = dict(
    D=25000,
    adjustment_step=1.5e-4,
    allowed_extra_profit=1e-7,
    fee_gamma=5e-3,
    ma_half_time=600,
    mid_fee=3e-4,
    out_fee=7e-3,
    n=2,
    log=0,
    gamma=5e-6,
    ext_fee=5e-4,
    gas_fee=0.017,
    A=100)

config = {
    'configuration': [],
    'datafile': ['crveth'],
    'debug': 0}

for x in X:
    params = copy(other_params)
    params[param] = x
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
