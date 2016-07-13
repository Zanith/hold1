import eyring_rate_script as test
import matplotlib.pyplot as plt
import shelve
from mpmath import nprint
import numpy as np

__author__ = 'Kyle Vitautas Lopin'


ion_conc22 = {'Nai': 0.12, 'Nae': 0.12,
              'Cai': 0.002, 'Cae': 0.002}
energy_barriers22 = {'distance': [0.167, 0.333, 0.5, 0.667, 0.8333],
                     'GNa': [8, -4, 20, -4, 8],
                     'GCa': [9, -12, 20, -12, 9]}

ion_conc32 = {'Nai': 0.145, 'Nae': 0.145,
              'Cai': 0.00000004, 'Cae': 0.002,
              'Mgi': 0.001, 'Mge': 0.0,
              'Fei': 0.0, 'Fee': 0.000001,
              'Bai': 0.1, 'Bae': 0.000001}
energy_barriers32 = {'distance': [0.095, 0.301, 0.353, 0.544, 0.999],
                     'GNa': [10.99, -2, 6.49, -2.9, 10.16, -2, 10, -2, 10],
                     'GCa': [8.12, -13.45, 0.96, -11.25, 10.64, -8, 10, -8, 10],
                     'GMg': [10.47, -9.62, 9.68, -6.44, 6.3, -8, 10, -15, 8],
                     'GFe': [10.47, -12.62, 9.68, -12.44, 6.3, -8, 10, -12, 10],
                     'GBa': [10.47, -12.62, 9.68, -12.44, 6.3, -8, 10, -8, 8]}

voltages = range(-150, 110, 10)
voltages2 = [-100, 100]
dps = range(50, 150, 10)
# print voltages
num_barriers = 3
# results_eig,
ion_conc = {'solute_1i': 0.001, 'solute_1e': 0.001}
energy_barriers = {'distance': [0.25, 0.5, 0.75], 'solute_1': [8.0, -10.0, 8.0]}

results_eig, results_svd = test.eyring_rate_algo(voltages,
                                                 ion_conc,
                                                 energy_barriers,
                                                 2,
                                                 [1],
                                                 [0.5, 0.9, 1, 0.5, 0.9, 0.5, 1, 0.9, 0.5],
                                                 10)
print 'test1'
"""
results has the fields
voltage
matrix_spec
ion_transport
fitting
current

"""
print 'test a'
result_shelf = shelve.open("Results1")
result_shelf["results eig"] = results_eig
result_shelf["results svd"] = results_svd
result_shelf.close()
ion_transport = []
currents = []
#
# Na_transport0 = []
# Na_transport1 = []
# Na_transport2 = []
#
# Ca_transport0 = []
# Ca_transport1 = []
# Ca_transport2 = []
# #
# Fe_transport0 = []
# Fe_transport1 = []
# Fe_transport2 = []
#
solute_transport = []
for result in results_svd:
    currents.append(result.current)
    print result.voltage
#     print result.ion_transport['Ca']
#     matrix_specs.append(result.matrix_spec.element_different)
#     transport_errors.append(result.fitting.transport_errors)
#     residue_errors.append(result.fitting.sae_residues)
#     steady_states.append(result.steady_state)
#     # print 'tester22: ', result.ion_transport['Fe'][0]
#     print result.ion_transport
#     Na_transport0.append(result.ion_transport['Na'][0])
#     Na_transport1.append(result.ion_transport['Na'][1])
#     Na_transport2.append(result.ion_transport['Na'][2])
# #
#     Ca_transport0.append(result.ion_transport['Ca'][0])
#     Ca_transport1.append(result.ion_transport['Ca'][1])
#     Ca_transport2.append(result.ion_transport['Ca'][2])
# print 'test4'
current_figure = plt.figure()
current_axis = current_figure.add_subplot(111)
currents2 = []
for current in currents:
    currents2.append(current[0])
# print "test||: ", len(currents2), len(voltages)
print voltages
print currents2
current_axis.plot(voltages, currents)
current_figure.suptitle("Current")
plt.ylabel("current (pA)")
plt.xlabel("voltage (mV)")
#
# Na_transport_figure = plt.figure()
# Na_transport_axis = Na_transport_figure.add_subplot(111)
# # Na_transport_axis.plot(voltages, Na_transport0, voltages, Na_transport1, voltages, Na_transport2,
# #                        voltages, Ca_transport0, voltages, Ca_transport1, voltages, Ca_transport2)
# Na_transport_axis.plot(voltages, Fe_transport0, voltages, Fe_transport1, voltages, Fe_transport2)
# Na_label_handlers = ["Fe 1", "Fe 2", "Fe 3"]
# Na_transport_axis.legend(Na_label_handlers)
# Na_transport_figure.suptitle("Fe transport")
#
# Ca_transport_figure = plt.figure()
# Ca_transport_axis = Ca_transport_figure.add_subplot(111)
# Ca_transport_axis.plot(voltages, Ca_transport0, voltages, Ca_transport1, voltages, Ca_transport2)
# Ca_transport_figure.suptitle("Ca transport")

# #  'steady states:'
# # print steady_states
# plt.xlabel("voltage (mV)")c
# print 'currents: '
# for current in currents:
#     nprint(current[1])
#
# print "Fe transport 0"
# for fe1 in Fe_transport0:
#     nprint(fe1)
#
# print "Fe transport 1"
# for fe1 in Fe_transport1:
#     nprint(fe1)
#
# print "Fe transport 2"
# for fe1 in Fe_transport2:
#     nprint(fe1)
# print 'show'
# print 'done'
#
# error_figure = plt.figure()
# error_axis = error_figure.add_subplot(111)
# error_plt = plt.plot(voltages, residue_errors)
plt.show()