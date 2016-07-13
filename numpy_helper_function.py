# from mpmath import mp, eig, fsum, fabs, matrix, inverse, norm, nprint, fp
import numpy as np
import eyring_rate_script as test

__author__ = 'Kyle Vitautas Lopin'


def solve_for_steady_state(_ss, _trans_matrix, _num_barriers=None):
    # 1: compute the sum of absolute errors (sae) and sum of squared errors of the residues left by multiplying the
    # transition matrix by the steady state solution (which should be all zeros if it is the true steady state)
    sae, sse = characterize_solution(_ss, _trans_matrix)

    # 2: calculate the ion transport from the ss calculated by the lowest eigenvalue
    ion_transport = test.eyring_rate_transport(_ss)

    # 3: calculate the current from the ions being transported for the eig and svd methods
    current = test.current_calc(ion_transport, _num_barriers)

    # 4: convert the ss to ints to save in the results section
    # TODO: this is broken?
    # int_ss = convert_to_int(_ss)
    # print 'check: ', int_ss
    # print 'why?: ', _ss
    # int_ss = 0

    # 5: Save the fitting results in custom class
    # fitting_specs = FittingMetrics(test_eig, sae, sse, ion_transport)

    return sse, sae, ion_transport, current, _ss


def svd_func(_matrix):
    try:
        U, S, V = np.linalg.svd(_matrix)
        smallest_eig_values = (S[len(S)-1], S[len(S)-2])
        ss_raw = V[len(V)-1, :]
        ss = ss_raw / np.sum(ss_raw)
    except Exception as e:
        print "SVD error", e
        return None, 0
    return ss.T, smallest_eig_values


def steady_state_eig(_matrix):
    """
    Find the steady state of a transition matrix by calculating the eigenvalues and eignevectors and using
    the eigenvector of the lowest eigenvalue
    :param _matrix:
    :return: the steady state distribution and the two
    """
    values, vectors = np.linalg.eig(_matrix)
    real_values = []
    for num in values:
        real_values.append(num.real)
    largest_eig_values = two_largest(real_values)  # want the largest numbers because they should all be negative
    _index = real_values.index(max(real_values))
    steady_state_raw = vectors[:, _index].real
    steady_state = steady_state_raw / np.sum(steady_state_raw)
    return steady_state, largest_eig_values


def characterize_solution(ss, _matrix):
    residues = residues_calc(ss, _matrix)
    sae_residues = np.linalg.norm(residues, 1)
    sse_residues = np.linalg.norm(residues, 2)
    return sae_residues, sse_residues


def residues_calc(null_vector, _matrix):
    return _matrix * null_vector


def two_largest(num_list):
    """
    Find the two largest numbers, for this module its use to find the 2 least negative numbers
    :param num_list: numbers to sort
    :return: 2 least negative numbers
    """
    return np.sort(num_list)[-2:]


def smallest_largest_elements(_matrix):
    l = len(_matrix) - 1
    smallest_element = np.amin(_matrix)
    largest_element = np.amax(_matrix)
    for i in range(l):
        for j in range(l):
            num = _matrix[i, j]
            if num != 0:
                abs_num = np.fabs(num)
                if abs_num < smallest_element:
                    smallest_element = abs_num
                if abs_num > largest_element:
                    largest_element = abs_num
    return largest_element, smallest_element


def convert_to_int(_vector):
    int_vector = []
    print 'hum: ', _vector
    for num in _vector:
        int_vector.append(num)
    print 'hum2 ', int_vector
    return int_vector


class MatrixSpecs(object):
    def __init__(self, _largest, _smallest, _condition_number):
        self.largest_element = _largest
        self.smallest_element = _smallest
        self.element_different = _largest - _smallest
        self.condition_number = _condition_number


class FittingMetrics(object):
    def __init__(self, _smallest, _sae, _sse, transport_rates):
        self.smallest_eig = _smallest[0]
        self.second_smallest_eig = _smallest[1]
        self.sae_residues = _sae
        self.sse_residues = _sse
        self.transport_errors = self.set_transport_errors(transport_rates)

    def set_transport_errors(self, transport_rates):
        _errors = dict()
        for key in transport_rates:
            _errors[key] = max(transport_rates[key]) - min(transport_rates[key])
        return _errors


class Results(object):
    def __init__(self, _voltage, _matrix_specs, _ion_transport, _fitting_specs, _current, _ss):
        self.voltage = _voltage
        self.matrix_spec = _matrix_specs
        self.ion_transport = _ion_transport
        self.fitting = _fitting_specs
        self.current = _current
        self.steady_state = _ss
