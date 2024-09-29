#!/usr/bin/env python3

import datetime
import time
import requests
import json

POOL = "0x4eBdF703948ddCEA3B11f675B4D1Fba9d2414A14"
crvUSD = "0xf939E0A03FB07F59A73314E73794Be0E57ac1b4E"
ETH = "0xC02aaA39b223FE8D0A0e5C4F27eAD9083C756Cc2"
CRV = "0xD533a949740bb3306d119CC777fa900bA034cd52"

pairs = {
    'ETHUSD': (crvUSD, ETH),
    'CRVUSD': (crvUSD, CRV),
    'CRVETH': (ETH, CRV)
}

CHUNKS_IN_DAY = 5
URI_TEMPLATE = 'https://prices.curve.fi/v1/ohlc/ethereum/%s?main_token={base_token}&reference_token={trade_token}&agg_number=1&agg_units=minute&start={start}&end={end}' % POOL

s = requests.session()
begin = datetime.datetime(year=2023, month=8, day=1)
start = datetime.datetime.utcnow()
dt = (start - begin)

for pair, (base_token, trade_token) in pairs.items():
    data = []
    for i in range(-dt.days * CHUNKS_IN_DAY, 0):
        d = int((start + datetime.timedelta(days=1) * i / CHUNKS_IN_DAY).timestamp())
        uri = URI_TEMPLATE.format(start=d, end=d + 86400 // CHUNKS_IN_DAY - 1, base_token=base_token, trade_token=trade_token)
        while True:
            try:
                time.sleep(0.1)
                resp = s.get(uri).json()
                break
            except Exception as e:
                print(e)
                time.sleep(30)
        print(pair, datetime.datetime.fromtimestamp(d), len(resp['data']))
        data += [[r['time'] * 1000, r['open'], r['high'], r['low'], r['close'], 1000000] for r in resp['data']]

    with open(pair.lower() + '.json', 'w') as f:
        json.dump(data, f)
