from numpy import exp, matrix

import np_func as helper

states = [[0, 0, 0],
          ['solute_1', 0, 0],
          ['solute_2', 0, 0],
          [0, 'solute_1', 0],
          [0, 'solute_2', 0],
          ['solute_1', 'solute_1', 0],
          ['solute_2', 'solute_1', 0],
          [0, 0, 'solute_1'],
          ['solute_1', 'solute_2', 0],
          ['solute_2', 'solute_2', 0],
          [0, 0, 'solute_2'],
          ['solute_1', 0, 'solute_1'],
          ['solute_2', 0, 'solute_1'],
          ['solute_1', 0, 'solute_2'],
          ['solute_2', 0, 'solute_2'],
          [0, 'solute_1', 'solute_1'],
          [0, 'solute_2', 'solute_1'],
          [0, 'solute_1', 'solute_2'],
          [0, 'solute_2', 'solute_2'],
          ['solute_1', 'solute_1', 'solute_1'],
          ['solute_2', 'solute_1', 'solute_1'],
          ['solute_1', 'solute_2', 'solute_1'],
          ['solute_2', 'solute_2', 'solute_1'],
          ['solute_1', 'solute_1', 'solute_2'],
          ['solute_2', 'solute_1', 'solute_2'],
          ['solute_1', 'solute_2', 'solute_2'],
          ['solute_2', 'solute_2', 'solute_2']]

ions = ['solute_1', 'solute_2']
ion_charges = {'solute_1': 0, 'solute_2': 0}
q_charge = 1.602e-19


def eyring_rate_algo(voltages, ion_conc, energy_barriers, Qs=None, Rs=1, mp_dps=15):
    
    if not Qs:
        Qs = [1]
    results_eig = []
    results_svd = []
    results_qr = []
    for voltage in voltages:
        print 'voltage: ', voltage

        # 1: make the transition matrix
        trans_matrix = eyring_rate_matrix(voltage, ion_conc, energy_barriers, Qs, Rs)

        result_eig, result_svd, result_qr = helper.solve_eyring_rate_model(voltage, trans_matrix)

        results_eig.append(result_eig)
        results_svd.append(result_svd)
        results_qr.append(result_qr)

    return results_eig, results_svd, results_qr


def convert_mp_int(_vector):
    int_vector = []
    for num in _vector:
        int_vector.append(num)
    return int_vector


def eyring_rate_matrix(voltage, ion_concs, energy_barriers, Qs, Rs):
    k0 = 6.1*10**12
    V = voltage
    q = 1/25.  # unit is e- / kT

    solute_1i = ion_concs['solute_1i']
    solute_1e = ion_concs['solute_1e']
    solute_2i = ion_concs['solute_2i']
    solute_2e = ion_concs['solute_2e']

    global Q  

    Q = Qs[0]

    d1 = energy_barriers['distance'][0]
    d2 = energy_barriers['distance'][1]
    d3 = energy_barriers['distance'][2]
    d4 = energy_barriers['distance'][3]
    d5 = energy_barriers['distance'][4]
    d6 = energy_barriers['distance'][5]
    d7 = energy_barriers['distance'][6]

    Gsolute_11 = energy_barriers['solute_1'][0]
    Gsolute_12 = energy_barriers['solute_1'][1]
    Gsolute_13 = energy_barriers['solute_1'][2]
    Gsolute_14 = energy_barriers['solute_1'][3]
    Gsolute_15 = energy_barriers['solute_1'][4]
    Gsolute_16 = energy_barriers['solute_1'][5]
    Gsolute_17 = energy_barriers['solute_1'][6]

    Gsolute_21 = energy_barriers['solute_2'][0]
    Gsolute_22 = energy_barriers['solute_2'][1]
    Gsolute_23 = energy_barriers['solute_2'][2]
    Gsolute_24 = energy_barriers['solute_2'][3]
    Gsolute_25 = energy_barriers['solute_2'][4]
    Gsolute_26 = energy_barriers['solute_2'][5]
    Gsolute_27 = energy_barriers['solute_2'][6]

    global k_0_1_solute_1, k_1_0_solute_1, k_0_1_solute_2, k_1_0_solute_2
    global k_1_2_solute_1, k_2_1_solute_1, k_1_2_solute_2, k_2_1_solute_2
    global k_2_3_solute_1, k_3_2_solute_1, k_2_3_solute_2, k_3_2_solute_2
    global k_3_4_solute_1, k_4_3_solute_1, k_3_4_solute_2, k_4_3_solute_2

    k_0_1_solute_1 = solute_1e*k0*exp(-Gsolute_11)*exp(0*q*-d1*V)
    k_1_0_solute_1 = k0*exp(Gsolute_12-Gsolute_11)*exp(0*q*(d2-d1)*V)
    k_0_1_solute_2 = solute_2e*k0*exp(-Gsolute_21)*exp(0*q*-d1*V)
    k_1_0_solute_2 = k0*exp(Gsolute_22-Gsolute_21)*exp(0*q*(d2-d1)*V)
    k_1_2_solute_1 = k0*exp(Gsolute_12-Gsolute_13)*exp(0*q*(d2-d3)*V)
    k_2_1_solute_1 = k0*exp(Gsolute_14-Gsolute_13)*exp(0*q*(d4-d3)*V)
    k_1_2_solute_2 = k0*exp(Gsolute_22-Gsolute_23)*exp(0*q*(d2-d3)*V)
    k_2_1_solute_2 = k0*exp(Gsolute_24-Gsolute_23)*exp(0*q*(d4-d3)*V)
    k_2_3_solute_1 = k0*exp(Gsolute_14-Gsolute_15)*exp(0*q*(d4-d5)*V)
    k_3_2_solute_1 = k0*exp(Gsolute_16-Gsolute_15)*exp(0*q*(d6-d5)*V)
    k_2_3_solute_2 = k0*exp(Gsolute_24-Gsolute_25)*exp(0*q*(d4-d5)*V)
    k_3_2_solute_2 = k0*exp(Gsolute_26-Gsolute_25)*exp(0*q*(d6-d5)*V)
    k_3_4_solute_1 = k0*exp(Gsolute_16-Gsolute_17)*exp(0*q*(d6-d7)*V)
    k_4_3_solute_1 = solute_1i*k0*exp(-Gsolute_17)*exp(0*q*(1-d7)*V)
    k_3_4_solute_2 = k0*exp(Gsolute_26-Gsolute_27)*exp(0*q*(d6-d7)*V)
    k_4_3_solute_2 = solute_2i*k0*exp(-Gsolute_27)*exp(0*q*(1-d7)*V)

    transition_matrix = matrix([
    [-(k_0_1_solute_1 + k_0_1_solute_2 + k_4_3_solute_1 + k_4_3_solute_2), k_1_0_solute_1, k_1_0_solute_2, 0, 0, 0, 0, k_3_4_solute_1, 0, 0, k_3_4_solute_2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [k_0_1_solute_1, -(k_1_0_solute_1 + k_1_2_solute_1 + k_4_3_solute_1 + k_4_3_solute_2), 0, k_2_1_solute_1, 0, 0, 0, 0, 0, 0, 0, k_3_4_solute_1, 0, k_3_4_solute_2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [k_0_1_solute_2, 0, -(k_1_0_solute_2 + k_1_2_solute_2 + k_4_3_solute_1 + k_4_3_solute_2), 0, k_2_1_solute_2, 0, 0, 0, 0, 0, 0, 0, k_3_4_solute_1, 0, k_3_4_solute_2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, k_1_2_solute_1, 0, -(k_2_1_solute_1 + k_0_1_solute_1 + k_0_1_solute_2 + k_2_3_solute_1 + k_4_3_solute_1 + k_4_3_solute_2), 0, Q**0 * k_1_0_solute_1, Q**0 * k_1_0_solute_2, k_3_2_solute_1, 0, 0, 0, 0, 0, 0, 0, Q**0 * k_3_4_solute_1, 0, Q**0 * k_3_4_solute_2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, k_1_2_solute_2, 0, -(k_2_1_solute_2 + k_0_1_solute_1 + k_0_1_solute_2 + k_2_3_solute_2 + k_4_3_solute_1 + k_4_3_solute_2), 0, 0, 0, Q**0 * k_1_0_solute_1, Q**0 * k_1_0_solute_2, k_3_2_solute_2, 0, 0, 0, 0, 0, Q**0 * k_3_4_solute_1, 0, Q**0 * k_3_4_solute_2, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, k_0_1_solute_1, 0, -(Q**0 * k_1_0_solute_1 + Q**0 * k_2_3_solute_1 + k_4_3_solute_1 + k_4_3_solute_2), 0, 0, 0, 0, 0, k_3_2_solute_1, 0, 0, 0, 0, 0, 0, 0, Q**0 * k_3_4_solute_1, 0, 0, 0, Q**0 * k_3_4_solute_2, 0, 0, 0],
    [0, 0, 0, k_0_1_solute_2, 0, 0, -(Q**0 * k_1_0_solute_2 + Q**0 * k_2_3_solute_1 + k_4_3_solute_1 + k_4_3_solute_2), 0, 0, 0, 0, 0, k_3_2_solute_1, 0, 0, 0, 0, 0, 0, 0, Q**0 * k_3_4_solute_1, 0, 0, 0, Q**0 * k_3_4_solute_2, 0, 0],
    [k_4_3_solute_1, 0, 0, k_2_3_solute_1, 0, 0, 0, -(k_3_4_solute_1 + k_3_2_solute_1 + k_0_1_solute_1 + k_0_1_solute_2), 0, 0, 0, k_1_0_solute_1, k_1_0_solute_2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, k_0_1_solute_1, 0, 0, 0, -(Q**0 * k_1_0_solute_1 + Q**0 * k_2_3_solute_2 + k_4_3_solute_1 + k_4_3_solute_2), 0, 0, 0, 0, k_3_2_solute_2, 0, 0, 0, 0, 0, 0, 0, Q**0 * k_3_4_solute_1, 0, 0, 0, Q**0 * k_3_4_solute_2, 0],
    [0, 0, 0, 0, k_0_1_solute_2, 0, 0, 0, 0, -(Q**0 * k_1_0_solute_2 + Q**0 * k_2_3_solute_2 + k_4_3_solute_1 + k_4_3_solute_2), 0, 0, 0, 0, k_3_2_solute_2, 0, 0, 0, 0, 0, 0, 0, Q**0 * k_3_4_solute_1, 0, 0, 0, Q**0 * k_3_4_solute_2],
    [k_4_3_solute_2, 0, 0, 0, k_2_3_solute_2, 0, 0, 0, 0, 0, -(k_3_4_solute_2 + k_3_2_solute_2 + k_0_1_solute_1 + k_0_1_solute_2), 0, 0, k_1_0_solute_1, k_1_0_solute_2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, k_4_3_solute_1, 0, 0, 0, Q**0 * k_2_3_solute_1, 0, k_0_1_solute_1, 0, 0, 0, -(k_3_4_solute_1 + k_3_2_solute_1 + k_1_0_solute_1 + k_1_2_solute_1), 0, 0, 0, Q**0 * k_2_1_solute_1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, k_4_3_solute_1, 0, 0, 0, Q**0 * k_2_3_solute_1, k_0_1_solute_2, 0, 0, 0, 0, -(k_3_4_solute_1 + k_3_2_solute_1 + k_1_0_solute_2 + k_1_2_solute_2), 0, 0, 0, Q**0 * k_2_1_solute_2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, k_4_3_solute_2, 0, 0, 0, 0, 0, 0, Q**0 * k_2_3_solute_2, 0, k_0_1_solute_1, 0, 0, -(k_3_4_solute_2 + k_3_2_solute_2 + k_1_0_solute_1 + k_1_2_solute_1), 0, 0, 0, Q**0 * k_2_1_solute_1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, k_4_3_solute_2, 0, 0, 0, 0, 0, 0, Q**0 * k_2_3_solute_2, k_0_1_solute_2, 0, 0, 0, -(k_3_4_solute_2 + k_3_2_solute_2 + k_1_0_solute_2 + k_1_2_solute_2), 0, 0, 0, Q**0 * k_2_1_solute_2, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, k_4_3_solute_1, 0, 0, 0, 0, 0, 0, 0, k_1_2_solute_1, 0, 0, 0, -(Q**0 * k_3_4_solute_1 + Q**0 * k_2_1_solute_1 + k_0_1_solute_1 + k_0_1_solute_2), 0, 0, 0, Q**0 * k_1_0_solute_1, Q**0 * k_1_0_solute_2, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, k_4_3_solute_1, 0, 0, 0, 0, 0, 0, 0, k_1_2_solute_2, 0, 0, 0, -(Q**0 * k_3_4_solute_1 + Q**0 * k_2_1_solute_2 + k_0_1_solute_1 + k_0_1_solute_2), 0, 0, 0, 0, Q**0 * k_1_0_solute_1, Q**0 * k_1_0_solute_2, 0, 0, 0, 0],
    [0, 0, 0, k_4_3_solute_2, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_1_2_solute_1, 0, 0, 0, -(Q**0 * k_3_4_solute_2 + Q**0 * k_2_1_solute_1 + k_0_1_solute_1 + k_0_1_solute_2), 0, 0, 0, 0, 0, Q**0 * k_1_0_solute_1, Q**0 * k_1_0_solute_2, 0, 0],
    [0, 0, 0, 0, k_4_3_solute_2, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_1_2_solute_2, 0, 0, 0, -(Q**0 * k_3_4_solute_2 + Q**0 * k_2_1_solute_2 + k_0_1_solute_1 + k_0_1_solute_2), 0, 0, 0, 0, 0, 0, Q**0 * k_1_0_solute_1, Q**0 * k_1_0_solute_2],
    [0, 0, 0, 0, 0, k_4_3_solute_1, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_0_1_solute_1, 0, 0, 0, -(Q**0 * k_3_4_solute_1 + Q**0 * k_1_0_solute_1), 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, k_4_3_solute_1, 0, 0, 0, 0, 0, 0, 0, 0, k_0_1_solute_2, 0, 0, 0, 0, -(Q**0 * k_3_4_solute_1 + Q**0 * k_1_0_solute_2), 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, k_4_3_solute_1, 0, 0, 0, 0, 0, 0, 0, k_0_1_solute_1, 0, 0, 0, 0, -(Q**0 * k_3_4_solute_1 + Q**0 * k_1_0_solute_1), 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, k_4_3_solute_1, 0, 0, 0, 0, 0, 0, k_0_1_solute_2, 0, 0, 0, 0, 0, -(Q**0 * k_3_4_solute_1 + Q**0 * k_1_0_solute_2), 0, 0, 0, 0],
    [0, 0, 0, 0, 0, k_4_3_solute_2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_0_1_solute_1, 0, 0, 0, 0, 0, -(Q**0 * k_3_4_solute_2 + Q**0 * k_1_0_solute_1), 0, 0, 0],
    [0, 0, 0, 0, 0, 0, k_4_3_solute_2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_0_1_solute_2, 0, 0, 0, 0, 0, 0, -(Q**0 * k_3_4_solute_2 + Q**0 * k_1_0_solute_2), 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, k_4_3_solute_2, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_0_1_solute_1, 0, 0, 0, 0, 0, 0, -(Q**0 * k_3_4_solute_2 + Q**0 * k_1_0_solute_1), 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, k_4_3_solute_2, 0, 0, 0, 0, 0, 0, 0, 0, k_0_1_solute_2, 0, 0, 0, 0, 0, 0, 0, -(Q**0 * k_3_4_solute_2 + Q**0 * k_1_0_solute_2)]
    ])

    return transition_matrix


def eyring_rate_transport(steady_state):

    inward = dict()
    inward['solute_1'] = [0, 0, 0, 0]
    inward['solute_1'][0] = (k_0_1_solute_1 * steady_state[0] +
                             k_0_1_solute_1 * steady_state[3] +
                             k_0_1_solute_1 * steady_state[4] +
                             k_0_1_solute_1 * steady_state[7] +
                             k_0_1_solute_1 * steady_state[10] +
                             k_0_1_solute_1 * steady_state[15] +
                             k_0_1_solute_1 * steady_state[16] +
                             k_0_1_solute_1 * steady_state[17] +
                             k_0_1_solute_1 * steady_state[18])
    inward['solute_1'][1] = (k_1_2_solute_1 * steady_state[1] +
                             k_1_2_solute_1 * steady_state[11] +
                             k_1_2_solute_1 * steady_state[13])
    inward['solute_1'][2] = (k_2_3_solute_1 * steady_state[3] +
                             Q**0 * k_2_3_solute_1 * steady_state[5] +
                             Q**0 * k_2_3_solute_1 * steady_state[6])
    inward['solute_1'][3] = (k_3_4_solute_1 * steady_state[7] +
                             k_3_4_solute_1 * steady_state[11] +
                             k_3_4_solute_1 * steady_state[12] +
                             Q**0 * k_3_4_solute_1 * steady_state[15] +
                             Q**0 * k_3_4_solute_1 * steady_state[16] +
                             Q**0 * k_3_4_solute_1 * steady_state[19] +
                             Q**0 * k_3_4_solute_1 * steady_state[20] +
                             Q**0 * k_3_4_solute_1 * steady_state[21] +
                             Q**0 * k_3_4_solute_1 * steady_state[22])
    inward['solute_2'] = [0, 0, 0, 0]
    inward['solute_2'][0] = (k_0_1_solute_2 * steady_state[0] +
                             k_0_1_solute_2 * steady_state[3] +
                             k_0_1_solute_2 * steady_state[4] +
                             k_0_1_solute_2 * steady_state[7] +
                             k_0_1_solute_2 * steady_state[10] +
                             k_0_1_solute_2 * steady_state[15] +
                             k_0_1_solute_2 * steady_state[16] +
                             k_0_1_solute_2 * steady_state[17] +
                             k_0_1_solute_2 * steady_state[18])
    inward['solute_2'][1] = (k_1_2_solute_2 * steady_state[2] +
                             k_1_2_solute_2 * steady_state[12] +
                             k_1_2_solute_2 * steady_state[14])
    inward['solute_2'][2] = (k_2_3_solute_2 * steady_state[4] +
                             Q**0 * k_2_3_solute_2 * steady_state[8] +
                             Q**0 * k_2_3_solute_2 * steady_state[9])
    inward['solute_2'][3] = (k_3_4_solute_2 * steady_state[10] +
                             k_3_4_solute_2 * steady_state[13] +
                             k_3_4_solute_2 * steady_state[14] +
                             Q**0 * k_3_4_solute_2 * steady_state[17] +
                             Q**0 * k_3_4_solute_2 * steady_state[18] +
                             Q**0 * k_3_4_solute_2 * steady_state[23] +
                             Q**0 * k_3_4_solute_2 * steady_state[24] +
                             Q**0 * k_3_4_solute_2 * steady_state[25] +
                             Q**0 * k_3_4_solute_2 * steady_state[26])
    outward = dict()
    outward['solute_1'] = [0, 0, 0, 0]
    outward['solute_1'][0] = (k_1_0_solute_1 * steady_state[1] +
                              Q**0 * k_1_0_solute_1 * steady_state[5] +
                              Q**0 * k_1_0_solute_1 * steady_state[8] +
                              k_1_0_solute_1 * steady_state[11] +
                              k_1_0_solute_1 * steady_state[13] +
                              Q**0 * k_1_0_solute_1 * steady_state[19] +
                              Q**0 * k_1_0_solute_1 * steady_state[21] +
                              Q**0 * k_1_0_solute_1 * steady_state[23] +
                              Q**0 * k_1_0_solute_1 * steady_state[25])
    outward['solute_1'][1] = (k_2_1_solute_1 * steady_state[3] +
                              Q**0 * k_2_1_solute_1 * steady_state[15] +
                              Q**0 * k_2_1_solute_1 * steady_state[17])
    outward['solute_1'][2] = (k_3_2_solute_1 * steady_state[7] +
                              k_3_2_solute_1 * steady_state[11] +
                              k_3_2_solute_1 * steady_state[12])
    outward['solute_1'][3] = (k_4_3_solute_1 * steady_state[0] +
                              k_4_3_solute_1 * steady_state[1] +
                              k_4_3_solute_1 * steady_state[2] +
                              k_4_3_solute_1 * steady_state[3] +
                              k_4_3_solute_1 * steady_state[4] +
                              k_4_3_solute_1 * steady_state[5] +
                              k_4_3_solute_1 * steady_state[6] +
                              k_4_3_solute_1 * steady_state[8] +
                              k_4_3_solute_1 * steady_state[9])
    outward['solute_2'] = [0, 0, 0, 0]
    outward['solute_2'][0] = (k_1_0_solute_2 * steady_state[2] +
                              Q**0 * k_1_0_solute_2 * steady_state[6] +
                              Q**0 * k_1_0_solute_2 * steady_state[9] +
                              k_1_0_solute_2 * steady_state[12] +
                              k_1_0_solute_2 * steady_state[14] +
                              Q**0 * k_1_0_solute_2 * steady_state[20] +
                              Q**0 * k_1_0_solute_2 * steady_state[22] +
                              Q**0 * k_1_0_solute_2 * steady_state[24] +
                              Q**0 * k_1_0_solute_2 * steady_state[26])
    outward['solute_2'][1] = (k_2_1_solute_2 * steady_state[4] +
                              Q**0 * k_2_1_solute_2 * steady_state[16] +
                              Q**0 * k_2_1_solute_2 * steady_state[18])
    outward['solute_2'][2] = (k_3_2_solute_2 * steady_state[10] +
                              k_3_2_solute_2 * steady_state[13] +
                              k_3_2_solute_2 * steady_state[14])
    outward['solute_2'][3] = (k_4_3_solute_2 * steady_state[0] +
                              k_4_3_solute_2 * steady_state[1] +
                              k_4_3_solute_2 * steady_state[2] +
                              k_4_3_solute_2 * steady_state[3] +
                              k_4_3_solute_2 * steady_state[4] +
                              k_4_3_solute_2 * steady_state[5] +
                              k_4_3_solute_2 * steady_state[6] +
                              k_4_3_solute_2 * steady_state[8] +
                              k_4_3_solute_2 * steady_state[9])


    transport_rate = dict()
    for ion in ions:
        transport_rate[ion] = []

        for i in range(len(inward.values()[0])):
            transport_rate[ion].append((outward[ion][i].real - inward[ion][i].real))

    return transport_rate


def current_calc(transport_rates):
    """
    Calculate the current for each barrier using the ion transport rates
    :param transport_rates:  dictionary of ions (as keys) with list of rates as values
    :return: dictionary of currents (in picoamps) caused by each ion and the 'total' current
    """
    all_currents = []
    # get number of barriers by looking at a transport value
    num_barriers = len(transport_rates[transport_rates.keys()[0]])
    for barrier in range(num_barriers):
        current = 0
        for ion in ions:
            # calcuate current as the ion charge (coulomb) times the transport rate (per second)
            current += ion_charges[ion] * q_charge * 10**12 * transport_rates[ion][barrier]  # 10e12 is to convert to picoamps
        all_currents.append(current)
    return all_currents