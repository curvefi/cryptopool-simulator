import numpy as np
from math import log10
import json
from copy import copy


X = np.logspace(log10(1e-4), log10(15e-4), 32)
param = 'mid_fee'

other_params = dict(
    D=50e6,
    adjustment_step=0.55e-5,
    allowed_extra_profit=1e-8,
    fee_gamma=1.2e-2,
    ma_half_time=600,
    mid_fee=15e-4,
    out_fee=100e-4,
    n=2,
    log=0,
    gamma=1e-4,
    ext_fee=2e-4,
    gas_fee=50,
    A=5000.0)

config = {
    'configuration': [],
    'datafile': ['eurusdt-x2'],
    'debug': 0}

for x in X:
    params = copy(other_params)
    params[param] = x
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
