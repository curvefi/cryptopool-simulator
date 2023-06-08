import numpy as np
from math import log10
import json
import itertools
from copy import copy


X = np.logspace(log10(20), log10(2000), 32)
Xname = "A"
Y = np.logspace(log10(3e-5), log10(3e-3), 32)
Yname = "gamma"

other_params = dict(
    D=1000.0,
    adjustment_step=1e-7,
    fee_gamma=0.006,
    ma_half_time=600,
    mid_fee=1e-4,
    out_fee=5e-4,
    gas_fee=0,
    gas_currency=0,
    n=3,
    log=0,
    allowed_extra_profit=1e-10,
    ext_fee=0,
    gamma=3e-4,
    A=100.0)

config = {
    'configuration': [],
    'datafile': [
        "xignite-eurusd-1m",
        "xignite-xauusd-1m",
        "xaueur-1m"],
    'debug': 0}

for x, y in itertools.product(X, Y):
    params = copy(other_params)
    params[Xname] = x
    params[Yname] = y
    config['configuration'].append(params)

with open('configuration.json', 'w') as f:
    json.dump(config, f)
