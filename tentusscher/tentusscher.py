import matplotlib.pyplot as plt
import myokit
import numpy as np
import time

# 한글 깨짐 방지
plt.rc('font', family='malgun gothic')
# 마이너스 깨짐 방지
plt.rcParams['axes.unicode_minus'] = False

#
# This example file uses the O'Hara-Rudy model in epicardial mode to produce
# and display a single AP.
#

# Get model from magic method
m, p, x = myokit.load('tentusscher.mmt')


# Create a pacing protocol
# bcl = 1000
# p = myokit.pacing.blocktrain(bcl, 1 , offset=20)

# Create simulation and run
s = myokit.Simulation(m, p)

# Pre-pace for a hundred beats
# s.pre(bcl*100)

# Run simulation
d = s.run(1000, log_interval=0.001)
# a = myokit.SimulationOpenCL(d['membrane.V'])

# print(a)
print(np.shape(d['membrane.V']))

# Display the result
fig, axes = plt.subplots(2, 4, figsize=(15,15))
plt.suptitle('myokit in tentusscher')
axes[0][0].plot(d['engine.time'], d['membrane.V'])
axes[0][1].plot(d['engine.time'], d['CD.Cai'])
axes[0][2].plot(d['engine.time'], d['FNC.INa'])
axes[0][3].plot(d['engine.time'], d['TOC.Ito'])
axes[1][0].plot(d['engine.time'], d['RDRC.IKr'])
axes[1][1].plot(d['engine.time'], d['SDRC.IKs'])
axes[1][2].plot(d['engine.time'], d['LCC.ICaL'])
axes[1][3].plot(d['engine.time'], d['IRKC.IK1'])
plt.show()
