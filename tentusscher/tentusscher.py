import matplotlib.pyplot as plt
import myokit
import numpy as np
import time

# 한글 깨짐 방지
plt.rc('font', family='malgun gothic')
# 마이너스 깨짐 방지
plt.rcParams['axes.unicode_minus'] = False

# Get model from magic method
m, p, x = myokit.load('tentusscher.mmt')

# Create simulation and run
s = myokit.Simulation(m, p)

# Run simulation
d = s.run(1000, log_interval=0.001)

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
