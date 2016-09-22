# Copyright (c) 2015-2016 Kyle Lopin (Naresuan University) <kylel@nu.ac.th>
# Licensed under the GPL

""" Helper functions to solve an eyring rate model and custom classes to save results
"""
# standard libraries
import logging
import time

import numpy as np

import eyring_rate_script as eyring_script

__author__ = 'Kyle Vitautas Lopin'

# logging setup to measure time to calculate steady states
LOGGER = logging.getLogger('timer')
HANDLER = logging.FileHandler('np_timing.txt')
LOGGER.addHandler(HANDLER)
LOGGER.setLevel(logging.INFO)


def solve_eyring_rate_model(voltage, transition_matrix):
    """
    Save all the matrix specifications, calculate the steady state of the matrix,
    NOTE: using 3 methods for testing purposes
    then save all the results.  Results saves in custom class.
    NOTE: too many local variables for testing; TODO: pick best method for final use
    :param voltage:  voltage the matrix is at, used to save conditions, TODO: take out for final use
    :param transition_matrix: numpy matrix of transition rates
    :return: custom class that saves the results, see Results class at bottom of file
    """
    # get the largest and smallest elements from the matrix to characterize the difficulty of
    # solving null space of the matrix and save it to be retrieved later
    len_matrix = len(transition_matrix)
    largest_matrix_element, smallest_matrix_element = smallest_largest_elements(transition_matrix)
    condition_number = np.linalg.cond(transition_matrix)
    matrix_specs = MatrixSpecs(largest_matrix_element,
                               smallest_matrix_element,
                               condition_number)

    # get the steady state (ss) solution by using the eigenvector of the lowest eigenvalue
    start = time.time()
    ss_by_eig, test_eigs_by_eig = steady_state_eig(transition_matrix)
    eig_time = time.time()-start
    # save all the results in a custom data class and return it
    results_eig = solve_eyring_rate_model_ss(voltage, ss_by_eig,
                                             transition_matrix, test_eigs_by_eig,
                                             matrix_specs)

    # get the steady state (ss) solution by using the svd decomposition
    start = time.time()
    ss_by_svd, test_eig_by_svd = svd_func(transition_matrix)
    svd_time = time.time()-start
    if any(ss_by_svd):  # in case the svd fails because of singularity
        results_svd = solve_eyring_rate_model_ss(voltage, ss_by_svd,
                                                 transition_matrix, test_eig_by_svd,
                                                 matrix_specs)
    else:
        results_svd = results_eig  # hack to make the program work if svd fails

    # get the steady state (ss) solution by using qr factorization
    start = time.time()
    ss_by_qr, test_eig_by_qr = qr_func(transition_matrix)
    qr_time = time.time()-start
    results_qr = solve_eyring_rate_model_ss(voltage, ss_by_qr,
                                            transition_matrix, test_eig_by_qr,
                                            matrix_specs)
    LOGGER.info('numpy times eig; svd; qr for matrix size %d: %5.10f %5.10f %5.10f',
                len_matrix, eig_time, svd_time, qr_time)
    return results_eig, results_svd, results_qr


def solve_eyring_rate_model_ss(voltage, _ss, _matrix, _test, _specs):
    """
    Take the steady state of a matrix, the matrix and calculate the transport rates of
    the solutes and current created
    :param voltage: voltage of the data point, to pass to results
    :param _ss: numpy vector of the steady state of the transition matrix
    :param _matrix: numpy matrix of transition rates
    :param _test: tuple of testing characteristics to pass through to results
    :param _specs: matrix specifications to pass to the results
    :return: custom class Results
    """
    # 3 calculate the ion transport rates
    _solute_transport = eyring_script.eyring_rate_transport(_ss)
    # transports are in 1x1 matrix form, so convert them to just scalars
    solute_transport = dict()
    for solute in _solute_transport:
        solute_transport[solute] = [np.asscalar(x) for x in _solute_transport[solute]]

    # characterize how well the steady state is by getting the residues of the steady
    # state times the transition matrix
    sum_squared_errors, sum_absolute_errors = characterize_solution(_ss, _matrix)

    # 4 calculate the current from the solute transport rates
    current = eyring_script.current_calc(solute_transport)

    # save the fitting results in a custom class
    fitting_specs_eig = FittingMetrics(_test, sum_absolute_errors,
                                       sum_squared_errors, solute_transport)

    return Results(voltage, _specs, solute_transport, fitting_specs_eig,
                   current, _ss)


def svd_func(_matrix):
    """
    Calculate the steady state of a matrix using Singular Value Decomposition
    :param _matrix: matrix to solve steady state for
    :return: vector of steady state and 2 smallest singular values
    """
    try:
        # pylint: disable=invalid-name
        U, S, V = np.linalg.svd(_matrix)  # ignore python naming convention for math convention
        smallest_singular_values = (S[len(S)-1], S[len(S)-2])
        steady_state_raw = V[len(V)-1, :]
        steady_state = steady_state_raw / np.sum(steady_state_raw)
    except np.linalg.linalg.LinAlgError as error:  # catch SVD not converging error,
        print "SVD error", error
        return [False], 0
    return steady_state.T, smallest_singular_values


def steady_state_eig(_matrix):
    """
    Find the steady state of a matrix by calculating the eigenvalues
    and eignevectors and using the eigenvector of the lowest eigenvalue
    :param _matrix: matrix to solve steady state for
    :return: vector of steady state and 2 smallest eigenvalues
    """
    values, vectors = np.linalg.eig(_matrix)
    real_values = []
    for num in values:
        # numerical errors cause the values to become complex for large matrices
        real_values.append(num.real)
        # want the largest numbers because they should all be negative
    largest_eig_values = two_largest(real_values)
    _index = real_values.index(max(real_values))
    steady_state_raw = vectors[:, _index].real
    steady_state = steady_state_raw / np.sum(steady_state_raw)
    return steady_state, largest_eig_values


def qr_func(_matrix):
    """
    Find the steady state of a matrix by qr factorization
    :param _matrix: matrix to solve steady state for
    :return: vector of steady state and 2 smallest upper-triangle diagonal values
    """
    # pylint: disable=invalid-name
    q, r = np.linalg.qr(_matrix.T)  # ignore python naming convention for math convention
    return q[:, -1]/np.sum(q[:, -1]), (r[-1, -1], r[-2, -2])


def characterize_solution(steady_state, _matrix):
    """
    Tell how well a steady state vector is in the null space of a matrix by
    finding the sum of absolute and sum of squared errors
    :param steady_state: steady state of a transition matrix (null space)
    :param _matrix: matrix to test
    :return: sum of absoluteerrors, sum of squared errors
    """
    # residues should be all zeros if fit perfectly, i.e. steady state is the null space
    residues = _matrix * steady_state
    sae_residues = np.linalg.norm(residues, 1)
    sse_residues = np.linalg.norm(residues, 2)
    return sae_residues, sse_residues


def two_largest(num_list):
    """
    Find the two largest numbers, for this module its use to find the 2 least negative numbers
    :param num_list: numbers to sort
    :return: 2 least negative numbers
    """
    return np.sort(num_list)[-2:]


def smallest_largest_elements(_matrix):
    """
    get the smallest and largest non-zero elements of a matrix
    :param _matrix:
    :return:
    """
    abs_matrix = np.fabs(_matrix)
    smallest_element = np.amin(abs_matrix[np.nonzero(abs_matrix)])
    largest_element = np.amax(abs_matrix[np.nonzero(abs_matrix)])
    return largest_element, smallest_element


# bad style below for numerical testing, set to be fixed after testing
class MatrixSpecs(object):
    """
    Store the specification of a matrix
    """
    def __init__(self, _largest, _smallest, _condition_number):
        self.largest_element = _largest
        self.smallest_element = _smallest
        self.element_different = _largest - _smallest
        self.condition_number = _condition_number


class FittingMetrics(object):
    """
    Store metrics about how well the calculated steady state solution to
    a transition matrix is actually in the null space
    """
    def __init__(self, _smallest, _sae, _sse, transport_rates):
        self.smallest_eig = _smallest[1]
        self.second_smallest_eig = _smallest[0]
        self.sae_residues = _sae
        self.sse_residues = _sse
        self.transport_errors = self.set_transport_errors(transport_rates)

    @staticmethod
    def set_transport_errors(transport_rates):
        """
        Calculate the errors in the transport rates, i.e. what is the difference in the
        transport rates between different barriers
        :param transport_rates: dict of transport rates, keys are solutes names
        and the values are lists of transport rate for each barrier
        :return: dict of errors, keys are solutes and values are errors
        """
        _errors = dict()
        for key in transport_rates:
            _errors[key] = max(transport_rates[key]) - min(transport_rates[key])
        return _errors


class Results(object):
    """
    Custom class to store the results, currently stores all ion transported over each barrier
    """
    def __init__(self, _voltage, _matrix_specs, _ion_transport, _fitting_specs, _current, _ss):
        self.voltage = _voltage
        self.matrix_spec = _matrix_specs
        self.ion_transport = _ion_transport
        self.fitting = _fitting_specs
        self.current = _current
        self.steady_state = _ss
