from numpy import exp, matrix

import numpy_helper_function as helper


states = [[0],
          ['solute_1']]

ions = ['solute_1']
ion_charges = {'solute_1': 0}
q_charge = 1.602e-19

def eyring_rate_algo(voltages, ion_conc, energy_barriers, num_barriers=None, Qs=[1], Rs=1, mp_dps=15):
    
    results_eig = []
    results_svd = []
    results_qr = []
    if not num_barriers:
        num_barriers = (len(energy_barriers) + 1) / 2
    for voltage in voltages:
        print 'voltage: ', voltage

        # 1: make the transition matrix
        trans_matrix = eyring_rate_matrix(voltage, ion_conc, energy_barriers, Qs, Rs)

        # 2: get the largest and smallest elements from the matrix to characterize the difficulty of solving null
        # space of the matrix and save it to be retrieved later
        # largest_matrix_element, smallest_matrix_element = helper.smallest_largest_elements(trans_matrix)
        # condition_number = helper.cond(trans_matrix)
        # matrix_specs = helper.MatrixSpecs(largest_matrix_element,
        #                                   smallest_matrix_element,
         #                                  condition_number)

        # 3a eig: get the steady state (ss) solution by using the eigenvector of the lowest eigenvalue
         #ss_by_eig, test_eigs_by_eig = helper.steady_state_eig(trans_matrix)

        # 4a calculate the ion transport rates
        # ion_transport_eig = eyring_rate_transport(ss_by_eig)

        # 5a for the steady state calculated by the eigenvalue method calculate the current and fitting
        # specs of the steady state
        # eig_results = helper.solve_for_steady_state(ss_by_eig,
        #                                             trans_matrix,
        #                                             num_barriers)
        # unpack the results
        # sse_eig = eig_results[0]
        # sae_eig = eig_results[1]
        # current_eig = eig_results[2]

        # 6a eig: Save the fitting results in custom class
        # fitting_specs_eig = helper.FittingMetrics(test_eigs_by_eig,
        #                                           sae_eig,
        #                                           sse_eig,
        #                                           ion_transport_eig)

        # 7a eig: save all results in custom data class
        # results_eig.append(helper.Results(voltage, matrix_specs,
        #                                   ion_transport_eig, fitting_specs_eig,
        #                                   current_eig, ss_by_eig))

        # 3b svd: try to get the ss using svd decomposition and getting the null space
        # NOTE: this method may not converge
        # ss_by_svd, test_eig_by_svd = helper.svd_func(trans_matrix)

        # if any(ss_by_svd):  # test if the svd decomposition converged before using it
            # 4b calculate the ion transport rates
            # ion_transport_svd = eyring_rate_transport(ss_by_svd)

            # 5b for the steady state calculated by the svd method calculate the ion transport, current and fitting
            # specs of the steady state
            # svd_results = helper.solve_for_steady_state(ss_by_svd, trans_matrix, num_barriers)

            # sse_svd = svd_results[0]
            # sae_svd = svd_results[1]
            # current_svd = svd_results[2]

            # 6b svd: Save the fitting results in custom class
            # fitting_specs_svd = helper.FittingMetrics(test_eig_by_svd,
            #                                           sae_svd,
            #                                           sse_svd,
            #                                           ion_transport_svd)
            # 7b: save all results in custom data class
            # results_svd.append(helper.Results(voltage, matrix_specs,
            #                                   ion_transport_svd, fitting_specs_svd,
            #                                   current_svd, ss_by_svd))
        # else:
            # print "SVD fail at voltage: ", voltage
            # results_svd.append(helper.Results(voltage, matrix_specs,
             #                                  None, None, None, None))

        # 3c qr: try to get the ss using qr decomposition and getting the null space
        # NOTE: this method may not converge
        # ss_by_qr, test_eig_by_qr = helper.qr_func(trans_matrix)

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



    Q = Qs[0]

    d1 = energy_barriers['distance'][0]
    d2 = energy_barriers['distance'][1]
    d3 = energy_barriers['distance'][2]

    Gsolute_11 = energy_barriers['solute_1'][0]
    Gsolute_12 = energy_barriers['solute_1'][1]
    Gsolute_13 = energy_barriers['solute_1'][2]

    global k_0_1_solute_1, k_1_0_solute_1, k_1_2_solute_1, k_2_1_solute_1

    k_0_1_solute_1 = solute_1e*k0*exp(-Gsolute_11)*exp(0*q*-d1*V)
    k_1_0_solute_1 = k0*exp(Gsolute_12-Gsolute_11)*exp(0*q*(d2-d1)*V)
    k_1_2_solute_1 = k0*exp(Gsolute_12-Gsolute_13)*exp(0*q*(d2-d3)*V)
    k_2_1_solute_1 = solute_1i*k0*exp(-Gsolute_13)*exp(0*q*(1-d3)*V)

    transition_matrix = matrix([
    [-(k_0_1_solute_1 + k_2_1_solute_1), k_1_0_solute_1 + k_1_2_solute_1],
    [k_0_1_solute_1 + k_2_1_solute_1, -(k_1_0_solute_1 + k_1_2_solute_1)]
    ])

    return transition_matrix


def eyring_rate_transport(steady_state):

    inward = dict()
    inward['solute_1'] = [0, 0]
    inward['solute_1'][0] = (k_0_1_solute_1 * steady_state[0])
    inward['solute_1'][1] = (k_1_2_solute_1 * steady_state[1])
    outward = dict()
    outward['solute_1'] = [0, 0]
    outward['solute_1'][0] = (k_1_0_solute_1 * steady_state[1])
    outward['solute_1'][1] = (k_2_1_solute_1 * steady_state[0])


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