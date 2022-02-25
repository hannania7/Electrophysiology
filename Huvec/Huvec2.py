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
    m, p, x = myokit.load('Huvec2.mmt')


    # Create a pacing protocol
    bcl = 1000
    p = myokit.pacing.blocktrain(bcl, 3 , offset=50)

    # Create simulation and run
    s = myokit.Simulation(m, p)

    # Pre-pace for a hundred beats
    s.pre(bcl*100)

    # Run simulation
    d = s.run(800, log_interval=0.01)
    return d

def f():
    m, p, x = myokit.load('Huvec2.mmt')
    m.set_value('CKR.I_Kr', 0)

    # Create a pacing protocol
    bcl = 1000
    p = myokit.pacing.blocktrain(bcl, 3 , offset=50)

    # Create simulation and run
    s = myokit.Simulation(m, p)

    # Pre-pace for a hundred beats
    s.pre(bcl*100)

    # Run simulation
    d = s.run(800, log_interval=0.01)
    return d


print(np.shape(g()['MP.Vm']))

# Display the result
plt.figure()
plt.suptitle('Huvec model')
plt.plot(g()['engine.time'], g()['MP.Vm'], 'black' \
        ,f()['engine.time'], f()['MP.Vm'], 'red')
plt.xlim([0,800])
plt.ylim([-100,60])
plt.xticks(np.arange(0,801,200))
plt.grid(True, color='blue', alpha=0.5, linestyle='-')
plt.legend(('Action potential', 'Action potential in IKr = 0'))     
plt.title('Action potential')
plt.show()