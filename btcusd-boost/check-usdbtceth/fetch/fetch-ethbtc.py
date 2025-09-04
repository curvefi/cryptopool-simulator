#!/usr/bin/env python3

import datetime
import time
import requests
import json

pair = 'ETHBTC'

URI_TEMPLATE = 'https://api.binance.com/api/v1/klines?symbol=%s&interval=1m&limit=500&startTime={start}&endTime={end}' % pair
data = []

begin = datetime.datetime(year=2024, month=12, day=28, hour=10)
start = datetime.datetime(year=2025, month=7, day=12, hour=10)
dt = (start - begin)
for i in range(-dt.days * 5, 0):
    d = int((start + datetime.timedelta(days=1) * i / 5).timestamp()) * 1000
    uri = URI_TEMPLATE.format(start=d, end=d + 86400 * 1000 // 5 - 1)
    resp = requests.get(uri).json()
    print(pair, datetime.datetime.fromtimestamp(d // 1000), len(resp))
    data += resp
    time.sleep(0.3)

with open(pair.lower() + '.json', 'w') as f:
    json.dump(data, f)
