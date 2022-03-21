"""Contains protocols to act in silico to probe cellular mechanics."""
import os, sys
import bisect
from typing import List, Union
import random
from math import floor
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.append('../Lib')
import mod_protocols



def mutate(bounds, value, normal_denom=20):
    new_val_offset = np.random.normal(
        loc=0,
        scale=abs(bounds[0] - bounds[1]) / normal_denom)

    new_value = value + new_val_offset
    while ((new_value > bounds[1]) or
            (new_value < bounds[0])):
        new_val_offset = np.random.normal(
            loc=0,
            scale=abs(bounds[0] - bounds[1]) / normal_denom)
        new_value = value + new_val_offset

    return new_value

class VoltageClampStep():
    """A step in a voltage clamp protocol."""

    def __init__(self, voltage : float=None, duration : float=None) -> None:
        self.voltage = voltage
        self.duration = duration

    def __str__(self):
        return '|STEP: Voltage: {}, Duration: {}|'.format(self.voltage, self.duration)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return (abs(self.voltage - other.voltage) < 0.001 and
                abs(self.duration - other.duration) < 0.001)

    def set_to_random_step(self, voltage_bounds, duration_bounds):
        self.voltage = random.uniform(*voltage_bounds)
        self.duration = random.uniform(*duration_bounds)
        
    def mutate(self, voltage_bounds, duration_bounds):
        self.voltage = mutate(voltage_bounds, self.voltage)
        self.duration = mutate(duration_bounds, self.duration)


class VoltageClampRamp():
    """A step in a voltage clamp protocol."""

    def __init__(self, voltage_start : float=None, 
                       voltage_end : float=None,
                       duration : float=None) -> None:
        self.voltage_start = voltage_start
        self.voltage_end = voltage_end
        self.duration = duration

    def __str__(self):
        return '|RAMP: Voltage Start: {}, Voltage End: {}, Duration: {}|'.format(
                self.voltage_start, self.voltage_end, self.duration)
                
    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return (abs(self.voltage - other.voltage) < 0.001 and
                abs(self.duration - other.duration) < 0.001)

    def get_voltage(self, time):        
        fraction_change = time / self.duration
        voltage_change = self.voltage_end - self.voltage_start
        return self.voltage_start + fraction_change * voltage_change 

    def set_to_random_step(self, voltage_bounds, duration_bounds):
        self.voltage_start = random.uniform(*voltage_bounds)
        self.voltage_end = random.uniform(*voltage_bounds)
        self.duration=random.uniform(*duration_bounds)

    def mutate(self, voltage_bounds, duration_bounds):
        self.voltage_start = mutate(voltage_bounds, self.voltage_start)
        self.voltage_end = mutate(voltage_bounds, self.voltage_end)
        self.duration = mutate(duration_bounds, self.duration)


class VoltageClampProtocol():
    """Encapsulates state and behavior of a voltage clamp protocol."""

    HOLDING_STEP = VoltageClampStep(voltage=-87.0, duration=1.0)

    def __init__(self, steps: list=None):        
        self.steps = steps if steps else []
        
    def __str__(self):
        return ' + '.join([i.__str__() for i in self.steps])

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):

        if not isinstance(other, self.__class__):
            return False

        if len(other.steps) != len(self.steps):
            return False

        for i in range(len(other.steps)):
            if other.steps[i] != self.steps[i]:
                return False
        return True    

    def add(self, step):
        self.steps.append(step)

    def get_voltage_change_endpoints(self) -> List[float]:
        """Initializes voltage change endpoints based on the steps provided.

        For example, if the steps provided are:
            VoltageClampStep(voltage=1, duration=1),
            VoltageClampStep(voltage=2, duration=0.5),
        the voltage change points would be at 1 second and 1.5 seconds.

        Returns:
            A list of voltage change endpoints.
        """

        voltage_change_endpoints = []
        cumulative_time = 0
        for i in self.steps:
            cumulative_time += i.duration
            voltage_change_endpoints.append(cumulative_time)
        return voltage_change_endpoints

    def get_voltage_change_startpoints(self):
        voltage_change_endpoints = [0]
        cumulative_time = 0
        for i in self.steps[0:-1]:
            cumulative_time += i.duration
            voltage_change_endpoints.append(cumulative_time)
        return voltage_change_endpoints

    def get_voltage_at_time(self, time: float) -> float:
        """Gets the voltage based on provided steps for the specified time."""
        step_index = bisect.bisect_left( self.get_voltage_change_endpoints(), time )
        current_step = self.steps[step_index]
        time_into_step = time - self.get_voltage_change_startpoints()[step_index]            
        if isinstance(current_step, VoltageClampStep) or isinstance(current_step, mod_protocols.VoltageClampStep):
            return current_step.voltage
        elif isinstance(current_step, VoltageClampRamp) or isinstance(current_step, mod_protocols.VoltageClampRamp):
            return current_step.get_voltage(time_into_step)
        else:
            return current_step.get_voltage(time_into_step)

    def get_voltage_clamp_protocol(self, times=np.arange(0,1000,1)):            
        return [self.get_voltage_at_time(t) for t in times]
        
    
    def plot_voltage_clamp_protocol(self, times=None, saved_to=None, is_plotted=True, ax=None, fig=None, unit: str='ms'):
        # duration = self.get_voltage_change_endpoints()[-1]
        
        # times = np.arange(0, duration, dt)
        scale = 1000 if unit=='s' else 1         
        voltages = [self.get_voltage_at_time(t)*scale for t in times]
        
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(15, 5))
            ax.plot(times, voltages, 'k')
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            plt.xlabel("Time (ms)", fontsize=18)
            plt.ylabel("Voltages (mV)", fontsize=18)
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
        else:
            ax.plot(times, voltages, 'k')

        if saved_to:
            plt.savefig(saved_to)
            np.save("times", times)
            np.save("voltage_protocol", voltages)
            
        if is_plotted:
            plt.show()
        else:
            return fig, ax
  
    def to_csv(self, path, period=0.1):
        length_of_proto = self.get_voltage_change_endpoints()[-1]

        times = np.arange(0, length_of_proto, period)
        voltages = np.zeros(len(times))

        for i, t in enumerate(times):
            voltages[i] = self.get_voltage_at_time(t)

        vcp_pd = pd.DataFrame(
                {'Times (ms)': times, 'Voltages (mV)': voltages})

        vcp_pd.to_csv(path, index=False)

        return vcp_pd
    
    def to_csp(self, path):
        f = open(path, "w")
        f.write('<!DOCTYPE ClampProtocolML>\n')
        f.write('<Clamp-Suite-Protocol-v1.0>\n')
        f.write('<segment numSweeps="1">\n')

        step_number = 0
        for step in self.steps:
            if isinstance(step, VoltageClampStep) or isinstance(step, mod_protocols.VoltageClampStep):
                formatted_step = f'<step stepDuration="{step.duration}" holdingLevel1="{step.voltage}" stepNumber="{step_number}" stepType="0"/>\n'
            else:
                formatted_step = f'<step holdingLevel2="{step.voltage_end}" stepDuration="{step.duration}" holdingLevel1="{step.voltage_start}" stepNumber="{step_number}" stepType="1"/>\n'

            f.write(formatted_step)
            step_number += 1

        f.write('</segment>\n')
        f.write('</Clamp-Suite-Protocol-v1.0>')

        f.close()



class PacingProtocol(): # suggested by 
    def __init__(self, level=1, start=20, length=0.5, period=1000, multiplier=0, default_time_unit='ms'):
        '''
        '''
        self._pace = 1        
        self._level = level
        self._start = start
#         self._end = end
        self._length = length
        self._period = period
        self._multiplier = multiplier

        self._time_conversion = 1000.0
        if default_time_unit == 's':
            self._time_conversion = 1.0            
        else:
            self._time_conversion = 1000.0
            
    def pacing(self, t):  
        if self._start<=0 or self._length<=0:
            return 0  
           
        # stim_amplitude = protocol.stim_amplitude * 1E-3 * self._time_conversion
        stim_start = self._start * 1E-3 * self._time_conversion
        stim_duration = self._length * 1E-3 * self._time_conversion
#         stim_end = self._end * 1E-3 * self._time_conversion
        i_stim_period = self._period * 1E-3 *self._time_conversion / self._pace

        if self._time_conversion == 1:
            denom = 1E9
        else:
            denom = 1

        pace = (1.0 if t-stim_start - i_stim_period*floor((t-stim_start)/i_stim_period) <= stim_duration \
#                 and time <= stim_end \
                and t >= stim_start else 0) / denom
  
        return pace

    





class SpontaneousProtocol:
    """Encapsulates state and behavior of a single action potential protocol."""

    def __init__(self, duration=1800):
        self.duration = duration


class PacedProtocol:
    """
    Encapsulates state and behavior of a paced protocol
    
    model_name: "Paci", "Kernik", "OR"
    """
    def __init__(self,
            model_name,
            stim_end=6000,
            stim_start=10,
            pace=1,
            stim_mag=1):
        """

        """
        if (model_name == "Kernik"):
            self.stim_amplitude = 220 * stim_mag
            self.stim_duration = 5
        elif model_name == "OR":
            self.stim_amplitude = 80 * stim_mag
            self.stim_duration = 1
        elif (model_name == "Paci"):
            self.stim_amplitude = 220 * stim_mag
            self.stim_duration = 5

        self.pace = pace
        self.stim_end = stim_end
        self.stim_start = stim_start


class IrregularPacingProtocol:
    """Encapsulates state and behavior of a irregular pacing protocol.

    Attributes:
        duration: Duration of integration.
        stimulation_offsets: Each offset corresponds to the
            seconds after diastole begins that stimulation will
            occur. Cannot exceed `max_stim_interval_duration`, which is the
            time between beats when cell is pacing naturally.
    """

    # The start of a diastole must be below this voltage, in Vm.
    DIAS_THRESHOLD_VOLTAGE = -0.06

    # Set to time between naturally occurring spontaneous beats.
    _MAX_STIM_INTERVAL = 1.55

    STIM_AMPLITUDE_AMPS = 7.5e-10
    STIM_DURATION_SECS = 0.005

    def __init__(self, duration: int, stimulation_offsets: List[float]) -> None:
        self.duration = duration
        self.stimulation_offsets = stimulation_offsets
        self.all_stimulation_times = []

    @property
    def stimulation_offsets(self):
        return self._stimulation_offsets

    @stimulation_offsets.setter
    def stimulation_offsets(self, offsets):
        for i in offsets:
            if i > self._MAX_STIM_INTERVAL:
                raise ValueError(
                    'Stimulation offsets from diastolic start cannot be '
                    'greater than `self.max_stim_interval_duration` because '
                    'the cell will have started to spontaneously beat.')
        self._stimulation_offsets = offsets

    def make_offset_generator(self):
        return (i for i in self._stimulation_offsets)






class AperiodicPacingProtocol():
    """
    Encapsulates state and behavior of a paced protocol
    
    model_name: "Paci", "Kernik", "OR"
    """
    def __init__(self, model_name, duration=10000,
                 stim_starts=[442.80, 942.80, 2942.8,
                              4142.8, 4742.8, 5142.8,
                              6142.8, 6442.8, 7142.8,
                              8642.8, 9442.8]):

        if (model_name == "Kernik") or (model_name == "Paci"):
            self.stim_amplitude = 30
            self.stim_duration = 2
        elif model_name == "OR":
            self.stim_amplitude = 80
            self.stim_duration = 1

        self.duration = duration
        self.stim_starts = stim_starts





class VoltageClampSinusoid:
    """A sinusoidal step in a voltage clamp protocol."""

    def __init__(self, voltage_start=None,
                    voltage_amplitude=None,
                    voltage_frequency=None,
                    duration=None) -> None:
        self.voltage_start = voltage_start
        self.amplitude = voltage_amplitude
        self.frequency = voltage_frequency
        self.duration = duration

    def set_to_random_step(self,
                           voltage_bounds,
                           duration_bounds,
                           frq_bounds=[.005, .25],
                           amp_bounds=[5, 75]):
        v_start = -200
        amplitude = -200
        while (((v_start - abs(amplitude)) < voltage_bounds[0]) or
                ((v_start + abs(amplitude)) > voltage_bounds[1])):
            v_start = random.uniform(*voltage_bounds)
            amplitude = random.choice([-1, 1])*random.uniform(*amp_bounds)
        
        self.duration = random.uniform(*duration_bounds)
        self.frequency = random.uniform(*frq_bounds)
        self.voltage_start = v_start
        self.amplitude = amplitude

    def get_voltage(self, time):
        return np.sin(
                self.frequency*time) * self.amplitude + self.voltage_start

    def __str__(self):
        return '|SINUSOID: Voltage Start: {}, Duration: {}, Amplitude: {}, Frequency: {}|'.format(
                self.voltage_start, self.duration, self.amplitude, self.frequency)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return False
        if not isinstance(other, self.__class__):
            return False
        return (abs(self.voltage - other.voltage) < 0.001 and
                abs(self.duration - other.duration) < 0.001)

    def mutate(self, vcga_params, frq_bounds=[.005, .25]):
        v_bounds = vcga_params.config.ga_config.step_voltage_bounds
        d_bounds = vcga_params.config.ga_config.step_duration_bounds

        self.duration = mutate(d_bounds, self.duration)
        self.frequency = mutate(frq_bounds, self.frequency)

        v_start = self.voltage_start
        amplitude = self.amplitude

        while (((v_start - abs(amplitude)) < v_bounds[0]) or
                ((v_start + abs(amplitude)) > v_bounds[1])):
            v_offset = np.random.normal(
                loc=0,
                scale=abs(v_bounds[0] - v_bounds[1]) / 20)
            a_offset = np.random.normal(
                loc=0,
                scale=abs(v_bounds[0] - v_bounds[1]) / 20)
            v_start += v_offset
            amplitude += a_offset

        self.voltage_start = v_start
        self.amplitude = amplitude



PROTOCOL_TYPE = Union[
    SpontaneousProtocol,
    IrregularPacingProtocol,
    VoltageClampProtocol,
    PacedProtocol
]


