from mpmath import mp, exp, nprint
import matplotlib.pyplot as plt

__author__ = 'Kyle Vitautas Lopin'

voltages = range(-100, 110, 10)
ion_concs = {'Nai': 0.002, 'Nae': 0.002,
             'Cai': 0.002, 'Cae': 0.002,
             'Mgi': 0.001, 'Mge': 0.0}
energy_barriers = {# 'distance': [0.1667, 0.3333, 0.5, 0.6667, 0.8333],
                   'distance': [0.1, 0.3, 0.5, 0.8, 0.9],
                   'GNa': [8, -2, 8, -6, 8],
                   'GCa': [9, -8, 9, -12, 9],
                   'GMg': [9, -12, 15, -12, 9]}

Qs = [1]
Rs = [1, 1]


k0 = mp.mpf(6.1*10**12)
# k0 = mp.mpf(1.0)

q = 1/25.  # unit is e- / kT

Nai = ion_concs['Nai']
Nae = ion_concs['Nae']
Cai = ion_concs['Cai']
Cae = ion_concs['Cae']

# global Q12, R02, R13

Q12 = Qs[0]
R02 = Rs[0]
R13 = Rs[1]

d1 = energy_barriers['distance'][0]
d2 = energy_barriers['distance'][1]
d3 = energy_barriers['distance'][2]
d4 = energy_barriers['distance'][3]
d5 = energy_barriers['distance'][4]

GNa1 = energy_barriers['GNa'][0]
GNa2 = energy_barriers['GNa'][1]
GNa3 = energy_barriers['GNa'][2]
GNa4 = energy_barriers['GNa'][3]
GNa5 = energy_barriers['GNa'][4]

GCa1 = energy_barriers['GCa'][0]
GCa2 = energy_barriers['GCa'][1]
GCa3 = energy_barriers['GCa'][2]
GCa4 = energy_barriers['GCa'][3]
GCa5 = energy_barriers['GCa'][4]

# global k_0_1_Na, k_1_0_Na, k_0_1_Ca, k_1_0_Ca, k_1_2_Na, k_2_1_Na
# global k_1_2_Ca, k_2_1_Ca, k_2_3_Na, k_3_2_Na, k_2_3_Ca, k_3_2_Ca

k_0_1_Na = []
k_1_0_Na = []
k_0_1_Ca = []
k_1_0_Ca = []
k_1_2_Na = []
k_2_1_Na = []
k_1_2_Ca = []
k_2_1_Ca = []
k_2_3_Na = []
k_3_2_Na = []
k_2_3_Ca = []
k_3_2_Ca = []


for voltage in voltages:
    V = voltage
    k_0_1_Na.append(mp.mpf(Nae*k0*exp(-GNa1)*exp(1*q*-d1*V)))
    k_1_0_Na.append(mp.mpf(k0*exp(GNa2-GNa1)*exp(1*q*(d2-d1)*V)))
    k_0_1_Ca.append(mp.mpf(Cae*k0*exp(-GCa1)*exp(2*q*-d1*V)))
    k_1_0_Ca.append(mp.mpf(k0*exp(GCa2-GCa1)*exp(2*q*(d2-d1)*V)))
    k_1_2_Na.append(mp.mpf(k0*exp(GNa2-GNa3)*exp(1*q*(d2-d3)*V)))
    k_2_1_Na.append(mp.mpf(k0*exp(GNa4-GNa3)*exp(1*q*(d4-d3)*V)))
    k_1_2_Ca.append(mp.mpf(k0*exp(GCa2-GCa3)*exp(2*q*(d2-d3)*V)))
    k_2_1_Ca.append(mp.mpf(k0*exp(GCa4-GCa3)*exp(2*q*(d4-d3)*V)))
    k_2_3_Na.append(mp.mpf(k0*exp(GNa4-GNa5)*exp(1*q*(d4-d5)*V)))
    k_3_2_Na.append(mp.mpf(Nai*k0*exp(-GNa5)*exp(1*q*(1-d5)*V)))
    k_2_3_Ca.append(mp.mpf(k0*exp(GCa4-GCa5)*exp(2*q*(d4-d5)*V)))
    k_3_2_Ca.append(mp.mpf(Cai*k0*exp(-GCa5)*exp(2*q*(1-d5)*V)))

rates_figure = plt.figure()
rates_axis = rates_figure.add_subplot(111)
rates_axis.plot(voltages, k_0_1_Na,
                voltages, k_1_0_Na,
                voltages, k_0_1_Ca,
                voltages, k_1_0_Ca,
                voltages, k_1_2_Na,
                voltages, k_2_1_Na,
                voltages, k_1_2_Ca,
                voltages, k_2_1_Ca,
                voltages, k_2_3_Na,
                voltages, k_3_2_Na,
                voltages, k_2_3_Ca,
                voltages, k_3_2_Ca)
handles = ["k_0_1_Na", "k_1_0_Na", "k_0_1_Ca", "k_1_0_Ca",
           "k_1_2_Na", "k_2_1_Na", "k_1_2_Ca", "k_2_1_Ca",
           "k_2_3_Na", "k_3_2_Na", "k_2_3_Ca", "k_3_2_Ca"]
rates_axis.legend(handles)
rates_axis.set_yscale('log')


nprint(k_0_1_Na)
nprint(k_1_0_Na)
nprint(k_0_1_Ca)
nprint(k_1_0_Ca)
nprint(k_1_2_Na)
nprint(k_2_1_Na)
nprint(k_1_2_Ca)
nprint(k_2_1_Ca)
nprint(k_3_2_Na)
nprint(k_2_3_Na)
nprint(k_3_2_Ca)
nprint(k_2_3_Ca)

rates1_figure = plt.figure()
rates1_axis = rates1_figure.add_subplot(111)
rates1_axis.plot(voltages, k_0_1_Na,
                 voltages, k_1_0_Na,
                 voltages, k_0_1_Ca,
                 voltages, k_1_0_Ca)
handles1 = ["k_0_1_Na", "k_1_0_Na", "k_0_1_Ca", "k_1_0_Ca"]
rates1_axis.legend(handles1)
rates1_axis.set_yscale('log')


rates2_figure = plt.figure()
rates2_axis = rates2_figure.add_subplot(111)
rates2_axis.plot(voltages, k_1_2_Na,
                 voltages, k_2_1_Na,
                 voltages, k_1_2_Ca,
                 voltages, k_2_1_Ca)
handles2 = ["k_1_2_Na", "k_2_1_Na", "k_1_2_Ca", "k_2_1_Ca"]
rates2_axis.legend(handles2)
rates2_axis.set_yscale('log')


rates3_figure = plt.figure()
rates3_axis = rates3_figure.add_subplot(111)
rates3_axis.plot(voltages, k_2_3_Na,
                 voltages, k_3_2_Na,
                 voltages, k_2_3_Ca,
                 voltages, k_3_2_Ca)
handles3 = ["k_2_3_Na", "k_3_2_Na", "k_2_3_Ca", "k_3_2_Ca"]
rates3_axis.legend(handles3)
rates3_axis.set_yscale('log')


plt.show()