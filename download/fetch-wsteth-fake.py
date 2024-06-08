#!/usr/bin/env python3

import datetime
import requests
import json
import time

pair = 'WSTETHETH'
POOL = "0xDC24316b9AE028F1497c275EB9192a3Ea0f67022"
ETH = "0xEeeeeEeeeEeEeeEeEeEeeEEEeeeeEeeeeeeeEEeE"
STETH = "0xae7ab96520DE3A18E5e111B5EaAb095312D7fE84"

URI_TEMPLATE = 'https://prices.curve.fi/v1/ohlc/ethereum/%s?main_token=%s&reference_token=%s&agg_number=1&agg_units=minute&start={start}&end={end}' % (POOL, ETH, STETH)
CHUNKS_IN_DAY = 5
YEAR = 365 * 86400
RATE = 0.03
data = []

s = requests.session()
begin = datetime.datetime(year=2021, month=1, day=10)
begin_timestamp = begin.timestamp()
start = datetime.datetime.utcnow()
dt = (start - begin)
for i in range(-dt.days * CHUNKS_IN_DAY, 0):
    d = int((start + datetime.timedelta(days=1) * i / CHUNKS_IN_DAY).timestamp())
    uri = URI_TEMPLATE.format(start=d, end=d + 86400 // CHUNKS_IN_DAY - 1)
    print(uri)
    while True:
        try:
            time.sleep(0.1)
            resp = s.get(uri).json()
            break
        except Exception as e:
            print(e)
            time.sleep(30)
    print(pair, datetime.datetime.fromtimestamp(d), len(resp['data']))
    additional_mul = [(1 + RATE) ** ((int(r['time']) - begin_timestamp) / YEAR) for r in resp['data']]
    data += [[r['time'] * 1000, r['open'] * m, r['high'] * m, r['low'] * m, r['close'] * m, 1000000]
             for r, m in zip(resp['data'], additional_mul)]

with open(pair.lower() + '.json', 'w') as f:
    json.dump(data, f)
