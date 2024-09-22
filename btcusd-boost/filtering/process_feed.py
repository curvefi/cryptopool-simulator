#!/usr/bin/env python3

import json
import numpy as np

name = 'btcusdt-2023-2024'
prev_number = 60
processed_data = []

with open(f'{name}.json', 'r') as f:
    data = json.load(f)

for row in data:
    r = [int(row[0])] + [float(r) for r in row[1:6]]
    if len(processed_data) > prev_number:
        thismin = min(r[1:5])
        thismax = max(r[1:5])
        prev_prices = np.array([r[1:5] for r in processed_data[-prev_number:]])
        prev_mean = np.mean(prev_prices[-1])
        max_deviation = prev_prices.max() - prev_prices.min()
        if thismin < prev_mean - max_deviation:
            print(r[0], (prev_mean - thismin) / prev_mean, max_deviation / prev_mean)
            for i in range(len(r)):
                if r[i] == thismin:
                    r[i] = prev_mean - max_deviation
        if thismax > prev_mean + max_deviation:
            print(r[0], (thismax - prev_mean) / prev_mean, max_deviation / prev_mean)
            for i in range(len(r)):
                if r[i] == thismax:
                    r[i] = prev_mean + max_deviation

    processed_data.append(r)

with open(f'{name}-processed.json', 'w') as f:
    json.dump(processed_data, f)
