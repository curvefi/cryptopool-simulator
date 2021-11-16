import numpy as np
from math import log10
import json
import itertools
from copy import copy


X = np.logspace(log10(1e-4), log10(4e-3), 32)
Xname = "mid_fee"
Y = np.logspace(log10(1e-4), log10(1e-2), 32)
Yname = "out_fee"

other_params = dict(
    D=0.2e8,
    adjustment_step=1.5e-4,
    allowed_extra_profit=1e-7,
    fee_gamma=5e-3,
    ma_half_time=600,
    mid_fee=1e-3,
    out_fee=1e-3,
    n=2,
    log=0,
    gamma=5e-6,
    ext_fee=5e-4,
    gas_fee=70,
    A=1000)

config = {
    'configuration': [],
    'datafile': ['ethusdt-individual'],
    'debug': 0}

for x, y in itertools.product(X, Y):
    params = copy(other_params)
    params[Xname] = x
    params[Yname] = y
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
