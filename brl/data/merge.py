#!/usr/bin/env python3
import csv
import json
from datetime import datetime

NAMES = [
    'brl_usd_1min_one_year_2023.csv',
    'brl_usd_1min_one_year_2024.csv',
    'brl_usd_1min_one_year.csv',
]
OUT_NAME = 'brlusd-1m.json'

out = []

for name in NAMES:
    with open(name, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            to_record = [
                int(datetime.fromisoformat(row['datetime']).timestamp()),
                float(row['open']),
                float(row['high']),
                float(row['low']),
                float(row['close']),
                10**7
            ]
            out.append(to_record)

with open('brlusd-1m.json', 'w') as f:
    json.dump(out, f)
