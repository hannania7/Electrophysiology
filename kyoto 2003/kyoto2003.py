import matplotlib.pyplot as plt
import myokit
import numpy as np

# 한글 깨짐 방지
plt.rc('font', family='malgun gothic')
# 마이너스 깨짐 방지
plt.rcParams['axes.unicode_minus'] = False

#
# This example file uses the O'Hara-Rudy model in epicardial mode to produce
# and display a single AP.
#

# Get model from magic method
def g():
    m, p, x = myokit.load('kyoto2003.mmt')

    # Create a pacing protocol
    bcl = 400
    p = myokit.pacing.blocktrain(bcl, 2, offset=50)

    # Create simulation and run
    s = myokit.Simulation(m, p)

    # Pre-pace for a hundred beats
    s.pre(bcl*100)

    # Run simulation
    indice = m.get('membrane.Vm').indice()
    vt = 1*s.state()[indice]
    d, apds = s.run(400, log_interval=0.01, apd_variable='membrane.Vm', apd_threshold=vt)
    return d

print(np.shape(g()['membrane.Vm']))

# Display the result
plt.figure()
plt.suptitle('kyoto model in 2013')
plt.plot(g()['engine.time'], g()['membrane.Vm'])
plt.xlim([0,401])
plt.ylim([-100,60])
# plt.xticks(np.arange(0,401,200))
plt.yticks(np.arange(-100,51,50))
# plt.grid(True, color='blue', alpha=0.5, linestyle='-')
# plt.legend(('Action potential', 'Action potential in IKr = 0'))     
plt.title('Action potential')
plt.show()