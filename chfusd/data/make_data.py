#!/usr/bin/env python3

import lzma
import json
from datetime import datetime


FILENAMES = [
    'DAT_ASCII_USDCHF_M1_2019.csv.xz',
    'DAT_ASCII_USDCHF_M1_2020.csv.xz',
    'DAT_ASCII_USDCHF_M1_2021.csv.xz',
    'DAT_ASCII_USDCHF_M1_2022.csv.xz',
    'DAT_ASCII_USDCHF_M1_2023.csv.xz',
    'DAT_ASCII_USDCHF_M1_2024.csv.xz'
]

out = []


for name in FILENAMES:
    with lzma.open(name, 'rt') as f:
        for row in f.readlines():
            data = row.strip().split(';')
            row = [int(datetime.strptime(data[0], r"%Y%m%d %H%M%S").timestamp())] + [float(d) for d in data[1:-1]] + [1e6]
            out.append(row)


with open('chfusd-2019-2021.json', 'w') as f:
    json.dump(out, f)
