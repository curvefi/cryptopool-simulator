#!/usr/bin/env python3

import json
import numpy as np

name = 'ethusdt-sample'
window_size = 10
processed_data = []

with open(f'{name}.json', 'r') as f:
    data = json.load(f)

for row in data:
    r = [int(row[0])] + [float(r) for r in row[1:6]]
    window = data[max(len(processed_data) - window_size // 2, 0): min(len(processed_data) + window_size // 2 + 1, len(data))]
    window = np.array([[float(x) for x in d[1:5]] for d in window if d[0] != r[0]])
    min_window = window.min()
    max_window = window.max()
    for i in range(1, 5):
        if r[i] > max_window:
            r[i] = max_window
        if r[i] < min_window:
            r[i] = min_window
    processed_data.append(r)

with open(f'{name}-processed.json', 'w') as f:
    json.dump(processed_data, f)
