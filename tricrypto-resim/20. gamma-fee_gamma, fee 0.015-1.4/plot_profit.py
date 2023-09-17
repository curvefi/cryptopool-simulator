import re
import pylab
import numpy as np
from datetime import datetime

r = re.compile(r'^t\=(\d+).*xCP\-growth\:\ \{([\.\d]+)\}.*')
shift = 100

times = []
vps = []

with open('simulation.log', 'r') as f:
    for line in f.readlines():
        m = r.match(line)
        if m:
            t, vp = m.groups()
            times.append(int(t))
            vps.append(float(vp))

vps = np.array(vps)
times = np.array(times)

apys = np.array(((vps[shift:]/vps[:-shift] - 1) / (times[shift:] - times[:-shift]) + 1) ** (365*86400) - 1) * 100

pylab.semilogy([datetime.fromtimestamp(d) for d in times[:-shift]], apys)
pylab.xlabel('t')
pylab.ylabel('APY (%)')
pylab.show()
