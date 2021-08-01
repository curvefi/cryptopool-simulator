import numpy as np
from math import log10
import json
from copy import copy


X = np.logspace(log10(1e-4), log10(1e-3), 32)
param = 'gamma'

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
    gamma=1.4e-8,
    ext_fee=0.0,
    A=50.0)

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
