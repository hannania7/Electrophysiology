#-*- coding:utf-8 -*-
import matplotlib.pyplot as plt
import myokit
import numpy as np

# 한글 깨짐 방지
plt.rc('font', family='malgun gothic')
# 마이너스 깨짐 방지
plt.rcParams['axes.unicode_minus'] = False


# Get model from magic method


def f():
    m, p, x = myokit.load('Tor_ord_endo.mmt')
    m.set_value('RTF.celltype', 0)
    types = ['Endocardial', 'Epicardial', 'Midmyocardial']
    # Create a pacing protocol
    bcl = 1000
    p = myokit.pacing.blocktrain(bcl, 1 , offset=0)

    # Create simulation and run
    s = myokit.Simulation(m, p)

    # Pre-pace for a hundred beats
    s.pre(bcl*100)

    # Run simulation
    indice = m.get('membrane.V').indice()
    vt = 0.9*s.state()[indice]
    d, apds = s.run(10 * 500, log_interval=0.01, apd_variable='membrane.V', apd_threshold=vt)
    return d

def g():
    m, p, x = myokit.load('Tor_ord_epi.mmt')
    m.set_value('RTF.celltype', 1)
    types = ['Endocardial', 'Epicardial', 'Midmyocardial']
    # Create a pacing protocol
    bcl = 1000
    p = myokit.pacing.blocktrain(bcl, 1 , offset=0)

    # Create simulation and run
    s = myokit.Simulation(m, p)

    # Pre-pace for a hundred beats
    s.pre(bcl*100)

    # Run simulation
    indice = m.get('membrane.V').indice()
    vt = 0.9*s.state()[indice]
    d, apds = s.run(10 * 500, log_interval=0.01, apd_variable='membrane.V', apd_threshold=vt)
    return d

def k():
    m, p, x = myokit.load('Tor_ord_mid.mmt')
    m.set_value('RTF.celltype', 2)
    types = ['Endocardial', 'Epicardial', 'Midmyocardial']
    # Create a pacing protocol
    bcl = 1000
    p = myokit.pacing.blocktrain(bcl, 1 , offset=0)

    # Create simulation and run
    s = myokit.Simulation(m, p)

    # Pre-pace for a hundred beats
    s.pre(bcl*100)

    # Run simulation
    indice = m.get('membrane.V').indice()
    vt = 0.9*s.state()[indice]
    d, apds = s.run(10 * 500, log_interval=0.01, apd_variable='membrane.V', apd_threshold=vt)
    return d

print(np.shape(f()['membrane.V']))

# Display the result
plt.suptitle('myokit in Tor_ord')
plt.plot(f()['engine.time'], f()['membrane.V'], 'blue' \
        ,g()['engine.time'], g()['membrane.V'], 'red' \
        ,k()['engine.time'], k()['membrane.V'], 'orange')
plt.xlim([0,5000])
plt.ylim([-100,40])
plt.legend(('Endocardial', 'Epicardial', 'Midmyocardial'))
plt.title('Action potential')
plt.show()
