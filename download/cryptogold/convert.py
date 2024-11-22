import csv
import json
from datetime import datetime


fnames = [
    'xau_btc_1min',
    'xau_eth_1min',
    'eth_btc_1min'
]


for fname in fnames:
    csv_name = fname + '.csv'
    json_name = fname + '.json'
    out = []

    with open(csv_name) as f:
        r = csv.reader(f)
        next(r)
        for line in r:
            t = int(datetime.fromisoformat(line[0]).timestamp() * 1000)
            prices = [float(p) for p in line[1:]]
            out.append([t, prices[0], max(prices), min(prices), prices[-1], 1000])

    with open(json_name, 'w') as f:
        json.dump(out, f)
