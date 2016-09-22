from numpy import exp, matrix, asscalar, any
from numpy.linalg import cond

import np_func as helper

states = [[0],
          ['solute 1']]

ions = ['solute 1']
ion_charges = {'solute 1': 0}
q_charge = 1.602e-19

def eyring_rate_algo(voltages, ion_conc, energy_barriers, num_barriers=None, Qs=1, Rs=1, mp_dps=0):
    
    results_eig = []
    results_svd = []
    if not num_barriers:
        num_barriers = (len(energy_barriers) + 1) / 2
    for voltage in voltages:
        # 1: make the transition matrix
        trans_matrix = eyring_rate_matrix(voltage, ion_conc, energy_barriers, Qs, Rs)

        # 2: get the largest and smallest elements from the matrix to characterize the difficulty of solving null
        # space of the matrix and save it to be retrieved later
        largest_matrix_element, smallest_matrix_element = helper.smallest_largest_elements(trans_matrix)
        condition_number = cond(trans_matrix)
        matrix_specs = helper.MatrixSpecs(largest_matrix_element,
                                          smallest_matrix_element,
                                          condition_number)

        # 3a eig: get the steady state (ss) solution by using the eigenvector of the lowest eigenvalue
        ss_by_eig, test_eigs_by_eig = helper.steady_state_eig(trans_matrix)

        # for the steady state calculated by the eigenvalue method calculate the ion transport, current and fitting
        # specs of the steady state
        eig_results = helper.solve_for_steady_state(ss_by_eig,
                                                    trans_matrix,
                                                    num_barriers)
        # unpack the results
        sse_eig = eig_results[0]
        sae_eig = eig_results[1]
        ion_transport_eig = eig_results[2]
        current_eig = eig_results[3]
        int_ss_by_eig = eig_results[4]

        # 4a eig: compute the sum of absolute errors (sae) and sum of squared errors of the residues left by multiplying the
        # transition matrix by the steady state solution (which should be all zeros if it is the true steady state)
        # sae_eig, sse_eig = helper.characterize_solution(ss_by_eig, trans_matrix)

        # 5a eig: calculate the ion transport from the ss calculated by the lowest eigenvalue
        # ion_transport_eig = eyring_rate_transport(ss_by_eig)

        # 6a eig: calculate the current from the ions being transported for the eig and svd methods
        # current_eig = current_calc(ion_transport_eig, num_barriers)

        # 7a eig: convert the ss to ints to save in the results section
        # int_ss_by_eig = convert_mp_int(ss_by_eig)

        # 8a eig: Save the fitting results in custom class
        fitting_specs_eig = helper.FittingMetrics(test_eigs_by_eig,
                                                  sae_eig,
                                                  sse_eig,
                                                  ion_transport_eig)

        # 9a eig: save all results in custom data class
        results_eig.append(helper.Results(voltage, matrix_specs,
                                          ion_transport_eig, fitting_specs_eig,
                                          current_eig, int_ss_by_eig))

        # 3b svd: try to get the ss using svd decomposition and getting the null space
        # NOTE: this method may not converge
        ss_by_svd, test_eig_by_svd = helper.svd_func(trans_matrix)

        if any(ss_by_svd):  # test if the svd decomposition converged before using it
            # for the steady state calculated by the svd method calculate the ion transport, current and fitting
            # specs of the steady state
            svd_results = helper.solve_for_steady_state(ss_by_svd, trans_matrix, num_barriers)

            sse_svd = svd_results[0]
            sae_svd = svd_results[1]
            ion_transport_svd = svd_results[2]
            current_svd = svd_results[3]
            int_ss_by_svd = svd_results[4]

            # 4b svd: compute the sae and sse of residues of the svd ss solution
            # sae_svd, sse_svd = helper.characterize_solution(ss_by_svd, trans_matrix)
            
            # 5b svd: calculate the ion transport from the ss calculated by the svd method
            # ion_transport_svd = eyring_rate_transport(ss_by_svd)

            # 6b svd: calculate the current from the ions being transported for the svd method
            # current_svd = current_calc(ion_transport_svd, num_barriers)

            # 7b svd: convert the steady states' to ints to save in the results section
            # int_ss_by_svd = convert_mp_int(ss_by_svd)

            # 8b svd: Save the fitting results in custom class
            fitting_specs_svd = helper.FittingMetrics(test_eig_by_svd,
                                                      sae_svd,
                                                      sse_svd,
                                                      ion_transport_svd)
            results_svd.append(helper.Results(voltage, matrix_specs,
                                              ion_transport_svd, fitting_specs_svd,
                                              current_svd, int_ss_by_svd))
        else:
            print "SVD fail at voltage: ", voltage
            results_svd.append(helper.Results(voltage, matrix_specs,
                                              None, None, None, None))

    return results_eig, results_svd


def convert_mp_int(_vector):
    int_vector = []
    for num in _vector:
        int_vector.append(num)
    return int_vector


def eyring_rate_matrix(voltage, ion_concs, energy_barriers, Qs, Rs):
    k0 = 6.1*10**12
    V = voltage
    q = 1/25.  # unit is e- / kT

    solute 1i = ion_concs['solute 1i']
    solute 1e = ion_concs['solute 1e']



    R = Rs[0]

    d1 = energy_barriers['distance'][0]
    d2 = energy_barriers['distance'][1]
    d3 = energy_barriers['distance'][2]

    Gsolute 11 = energy_barriers['Gsolute 1'][0]
    Gsolute 12 = energy_barriers['Gsolute 1'][1]
    Gsolute 13 = energy_barriers['Gsolute 1'][2]

    global k_0_1_solute 1, k_1_0_solute 1, k_1_2_solute 1, k_2_1_solute 1

    k_0_1_solute 1 = solute 1e*k0*exp(-Gsolute 11)*exp(0*q*-d1*V)
    k_1_0_solute 1 = k0*exp(Gsolute 12-Gsolute 11)*exp(0*q*(d2-d1)*V)
    k_1_2_solute 1 = k0*exp(Gsolute 12-Gsolute 13)*exp(0*q*(d2-d3)*V)
    k_2_1_solute 1 = solute 1i*k0*exp(-Gsolute 13)*exp(0*q*(1-d3)*V)

    transition_matrix = matrix([
    [-(k_0_1_solute 1 + k_2_1_solute 1), k_1_0_solute 1 + k_1_2_solute 1],
    [k_0_1_solute 1 + k_2_1_solute 1, -(k_1_0_solute 1 + k_1_2_solute 1)]
    ])

    return transition_matrix

def eyring_rate_transport(steady_state):

    inward = dict()
    inward['solute 1'] = [0, 0]
    inward['solute 1'][0] = (k_0_1_solute 1 * steady_state[0])
    inward['solute 1'][1] = (k_1_2_solute 1 * steady_state[1])
    outward = dict()
    outward['solute 1'] = [0, 0]
    outward['solute 1'][0] = (k_1_0_solute 1 * steady_state[1])
    outward['solute 1'][1] = (k_2_1_solute 1 * steady_state[0])


    transport_rate = dict()
    for ion in ions:
        transport_rate[ion] = []

        for i in range(len(inward.values()[0])):
            transport_rate[ion].append(asscalar(outward[ion][i].real - inward[ion][i].real))
    return transport_rate


def current_calc(transport_rates, num_barriers):
    """
    Calculate the current for each barrier using the ion transport rates
    :param transport_rates:  dictionary of ions (as keys) with list of rates as values
    :return: dictionary of currents (in picoamps) caused by each ion and the 'total' current
    """
    all_currents = []
    for barrier in range(num_barriers):
        current = 0
        for ion in ions:
            # calcuate current as the ion charge (coulomb) times the transport rate (per second)
            current += ion_charges[ion] * q_charge * 10**12 * transport_rates[ion][barrier]  # 10e12 is to convert to picoamps
        all_currents.append(current)
    return all_currents