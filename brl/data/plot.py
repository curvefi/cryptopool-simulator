#!/usr/bin/python3

import json
import pylab

with open('brlusd-1m.json', 'r') as f:
    data = json.load(f)

t = [d[0] for d in data]
p = [d[1] for d in data]

pylab.plot(t[::1000], p[::1000])
pylab.show()
