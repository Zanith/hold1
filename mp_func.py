# Copyright (c) 2015-2016 Kyle Lopin (Naresuan University) <kylel@nu.ac.th>
# Licensed under the GPL

""" Helper functions to solve an eyring rate model and custom classes to save results
"""
# standard libraries
import logging
import time

from mpmath import mp, eig, fsum, fabs, norm, qr, mnorm

import eyring_rate_script as eyring_script

__author__ = 'Kyle Vitautas Lopin'

LOGGER = logging.getLogger('timer')
HANDLER = logging.FileHandler('timing_logger_mp_3x4extreme1_15dps.txt')
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
    len_matrix = len(transition_matrix)
    # get the largest and smallest elements from the matrix to characterize the difficulty of
    # solving null space of the matrix and save it to be retrieved later
    largest_matrix_element, smallest_matrix_element = smallest_largest_elements(transition_matrix)
    condition_number = mnorm(transition_matrix, 1)
    matrix_specs = MatrixSpecs(largest_matrix_element,
                               smallest_matrix_element,
                               condition_number)

    # get the steady state (ss) solution by using the eigenvector of the lowest eigenvalue
    start = time.time()
    ss_by_eig, test_eigs_by_eig = steady_state_eig(transition_matrix)
    eig_time = time.time()-start
    print 'eig: ', eig_time
    # save all the results in a custom data class and return it
    results_eig = solve_eyring_rate_model_ss(voltage, ss_by_eig,
                                             transition_matrix, test_eigs_by_eig,
                                             matrix_specs)
    # get the steady state (ss) solution by using the svd decomposition
    start = time.time()
    ss_by_svd, test_eig_by_svd = svd_func(transition_matrix)
    svd_time = time.time()-start
    print 'svd: ', svd_time
    if any(ss_by_svd):  # incase the svd fails because of singularity
        results_svd = solve_eyring_rate_model_ss(voltage, ss_by_svd,
                                                 transition_matrix, test_eig_by_svd,
                                                 matrix_specs)
    else:
        pass
        # results_svd = results_eig  # hack to make the program work
    # get the steady state (ss) solution by using qr factorization
    start = time.time()
    ss_by_qr, test_eig_by_qr = qr_func(transition_matrix)
    qr_time = time.time()-start
    print 'qr: ', qr_time
    results_qr = solve_eyring_rate_model_ss(voltage, ss_by_qr,
                                            transition_matrix, test_eig_by_qr,
                                            matrix_specs)
    LOGGER.info('mpmath times eig; svd; qr for matrix size %d: %5.10f %5.10f %5.10f',
                len_matrix, eig_time, svd_time, qr_time)
    print eig_time
    return results_eig, results_svd, results_qr


def solve_eyring_rate_model_ss(voltage, _ss, _matrix, _test, _specs):
    """
    Take the steady state of a matrix, the matrix and calculate the transport rates of
    the solutes and current created
    :param voltage: voltage of the data point, to pass to results
    :param _ss: vector of the steady state of the transition matrix
    :param _matrix: matrix of transition rates
    :param _test: tuple of testing characteristics to pass through to results
    :param _specs: matrix specifications to pass to the results
    :return: custom class Results
    """
    # 3 calculate the ion transport rates
    solute_transport = eyring_script.eyring_rate_transport(_ss)
    # characterize how well the steady state is by getting the residues of the steady state times
    # the transition matrix
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
    U, S, V = mp.svd(_matrix)
    smallest_eig_values = (S[len(S)-1], S[len(S)-2])
    steady_state_raw = V[len(V)-1, :]
    steady_state = steady_state_raw / fsum(steady_state_raw)
    return steady_state.T, smallest_eig_values


def steady_state_eig(_matrix):
    """
    Find the steady state of a matrix by calculating the eigenvalues
    and eignevectors and using the eigenvector of the lowest eigenvalue
    :param _matrix: matrix to solve steady state for
    :return: vector of steady state and 2 smallest eigenvalues
    """
    values, vectors = eig(_matrix)
    real_values = []
    for num in values:
        real_values.append(num.real)
    # want the largest numbers because they should all be negative
    largest_eig_values = two_largest(real_values)
    _index = real_values.index(max(real_values))
    steady_state_raw = vectors[:, _index]
    steady_state = steady_state_raw / fsum(steady_state_raw)
    return steady_state, largest_eig_values


def qr_func(_matrix):
    """
    Find the steady state of a matrix by qr factorization
    :param _matrix: matrix to solve steady state for
    :return: vector of steady state and 2 smallest upper-triangle diagonal values
    """
    Q, R = qr(_matrix.T)
    steady_state = Q[:, Q.cols-1]
    return steady_state/fsum(steady_state), (R[R.rows-1, R.cols-1], R[R.rows-2, R.cols-2])


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
    sae_residues = norm(residues, 1)
    sse_residues = norm(residues, 2)
    return sae_residues, sse_residues


def two_largest(num_list):
    """
    Find the two largest numbers, for this module its use to find the 2 least negative numbers
    :param num_list: numbers to sort
    :return: 2 least negative numbers
    """
    largest_num = mp.mpf(-1e100)
    second_largest_num = mp.mpf(-1e99)
    for num in num_list:
        if num > largest_num:
            second_largest_num = largest_num
            largest_num = num
        if largest_num < num < second_largest_num:
            second_largest_num = num
    return (largest_num, second_largest_num)


def smallest_largest_elements(_matrix):
    """
    get the smallest and largest non-zero elements of a matrix in absolute terms
    :param _matrix:  matrix to search
    :return: largest number, smallest number
    """
    _length = len(_matrix) - 1
    smallest_element = mp.mpf(1e100)  # initialize to very large number
    largest_element = mp.mpf(0)
    for i in range(_length):
        for j in range(_length):
            num = _matrix[i, j]
            if num != 0:
                abs_num = fabs(num)
                if abs_num < smallest_element:
                    smallest_element = abs_num
                if abs_num > largest_element:
                    largest_element = abs_num
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
