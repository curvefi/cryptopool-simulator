import numpy as np
from math import log10
import json
from copy import copy


X = np.logspace(log10(1e-7), log10(1e-2), 64)
param = 'gamma'

other_params = dict(
    D=0.2e8,
    adjustment_step=1.5e-4,
    allowed_extra_profit=1e-7,
    fee_gamma=1.2e-2,
    ma_half_time=600,
    mid_fee=1e-3,
    out_fee=1e-3,
    n=2,
    log=0,
    gamma=2e-5,
    ext_fee=5e-4,
    gas_fee=70,
    A=1000)

config = {
    'configuration': [],
    'datafile': ['ethusdt-individual'],
    'debug': 0}

for x in X:
    params = copy(other_params)
    params[param] = x
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
