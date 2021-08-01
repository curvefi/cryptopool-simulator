import numpy as np
from math import log10
import json
import itertools
from copy import copy


# X = np.logspace(log10(0.25), log10(1.18), 16)
X = np.logspace(log10(0.2), log10(200000), 32)
Xname = "A"
# Y = np.logspace(log10(9e-6), log10(0.0003), 16)
Y = np.logspace(log10(0.5e-7), log10(1e-2), 32)
Yname = "gamma"

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
    A=100.0)

config = {
    'configuration': [],
    'datafile': [
        'eurusdt'],
    'debug': 0}

for x, y in itertools.product(X, Y):
    params = copy(other_params)
    params[Xname] = x
    params[Yname] = y
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
