__author__ = 'Kyle Vitautas Lopin'


import template_maker as tm
import shelve
import importlib


ions = ['Na', 'Ca', 'Mg', 'Fe', 'Ba', 'Cd']
charges = [1, 2, 2, 2, 2, 2]
voltages = range(-200, 220, 20)
ion_conc = {'Nai': 0.1, 'Nae': 0.000001,
            'Cai': 0.000001, 'Cae': 0.1,
            'Mgi': 0.1, 'Mge': 0.000001,
            'Fei': 0.000001, 'Fee': 0.1,
            'Bai': 0.1, 'Bae': 0.000001,
            'Cdi': 0.000001, 'Cde': 0.1}

energy_barriers = {'GNa': [8, -12, 8, -12, 8, -12, 8, -12, 8, -12, 8, -12, 8],
                   'GCa': [8, -12, 8, -12, 8, -12, 8, -12, 8, -12, 8, -12, 8],
                   'GMg': [8, -12, 8, -12, 8, -12, 8, -12, 8, -12, 8, -12, 8],
                   'GFe': [8, -12, 8, -12, 8, -12, 8, -12, 8, -12, 8, -12, 8],
                   'GBa': [8, -12, 8, -12, 8, -12, 8, -12, 8, -12, 8, -12, 8],
                   'GCd': [8, -12, 8, -12, 8, -12, 8, -12, 8, -12, 8, -12, 8]}
Qs = 30*[5]
Rs = 30*[0.5]

total_results_eig = []
total_results_svd = []
num_binding_sites = 1
num_ions = 1
run_options = (num_binding_sites, num_ions)

# make the script make an eyring rate model and save it to test_script.py
eyring_rate_script = tm.make_template('numpy', num_binding_sites, ions[:num_ions], charges, 'full QR')
with open("test_script.py", "w") as file:
    file.write(eyring_rate_script)
file.close()

# import the eyring rate model file that was just made
import test_script

# make a list of the electrical distances to update the energy_barriers dict with
distance_buffer = 2 * (num_binding_sites + 1)
distance = []
for i in range(1, distance_buffer):
    distance.append(float(i)/distance_buffer)
energy_barriers['distance'] = distance

# call the erying rate model script made earlier
results_eig, results_svd = test_script.eyring_rate_algo(voltages,
                                                        ion_conc,
                                                        energy_barriers,
                                                        num_binding_sites, Qs, Rs)
print results_eig
print 'humm'
total_results_eig.append(results_eig)
total_results_svd.append(results_svd)


total_results_shelf = shelve.open("Full_Results")
total_results_shelf["total results eig"] = total_results_eig
total_results_shelf["total resutls svd"] = total_results_svd
total_results_shelf.close()
