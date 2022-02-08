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
plt.suptitle('myokit in tentusscher')
plt.subplot(2, 4, 1)
plt.plot(d['engine.time'], d['membrane.V'])
plt.title('Action potential')
plt.subplot(2, 4, 2)
plt.plot(d['engine.time'], d['CD.Cai'])
plt.title('Cai')
plt.subplot(2, 4, 3)
plt.plot(d['engine.time'], d['FNC.INa'])
plt.title('INa')
plt.subplot(2, 4, 4)
plt.plot(d['engine.time'], d['TOC.Ito'])
plt.title('Ito')
plt.subplot(2, 4, 5)
plt.plot(d['engine.time'], d['RDRC.IKr'])
plt.title('IKr')
plt.subplot(2, 4, 6)
plt.plot(d['engine.time'], d['SDRC.IKs'])
plt.title('IKs')
plt.subplot(2, 4, 7)
plt.plot(d['engine.time'], d['LCC.ICaL'])
plt.title('ICaL')
plt.subplot(2, 4, 8)
plt.plot(d['engine.time'], d['IRKC.IK1'])
plt.title('iK1')
plt.show()
