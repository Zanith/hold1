import numpy as np

import eyring_rate_script as eyring_script


__author__ = 'Kyle Vitautas Lopin'


def solve_eyring_rate_model(voltage, transition_matrix):
    # get the largest and smallest elements from the matrix to characterize the difficulty of solving null
    # space of the matrix and save it to be retrieved later
    largest_matrix_element, smallest_matrix_element = smallest_largest_elements(transition_matrix)
    condition_number = cond(transition_matrix)
    matrix_specs = MatrixSpecs(largest_matrix_element,
                               smallest_matrix_element,
                               condition_number)

    # 2a eig: get the steady state (ss) solution by using the eigenvector of the lowest eigenvalue
    ss_by_eig, test_eigs_by_eig = steady_state_eig(transition_matrix)

    # save all the results in a custom data class and return it
    results_eig = solve_eyring_rate_model_ss(voltage, ss_by_eig,
                                             transition_matrix, test_eigs_by_eig,
                                             matrix_specs)

    ss_by_svd, test_eig_by_svd = svd_func(transition_matrix)

    if any(ss_by_svd):  # incase the svd fails because of singularity
        results_svd = solve_eyring_rate_model_ss(voltage, ss_by_svd,
                                                 transition_matrix, test_eig_by_svd,
                                                 matrix_specs)
    else:
        results_svd = results_eig  # hack to make the program work

    ss_by_qr, test_eig_by_qr = qr_func(transition_matrix)

    results_qr = solve_eyring_rate_model_ss(voltage, ss_by_qr,
                                            transition_matrix, test_eig_by_qr,
                                            matrix_specs)

    return results_eig, results_svd, results_qr


def solve_eyring_rate_model_ss(voltage, _ss, _matrix, _test, _specs):
    # 3 calculate the ion transport rates
    _solute_transport = eyring_script.eyring_rate_transport(_ss)
    # transports are in 1x1 matrix form, so convert them to just scalars
    solute_transport = dict()
    for solute in _solute_transport:
        solute_transport[solute] = [np.asscalar(x) for x in _solute_transport[solute]]

    # characterize how well the steady state is by getting the residues of the steady state times the transition matrix
    sum_squared_errors, sum_absolute_errors = characterize_solution(_ss, _matrix)

    # 4 calculate the current from the solute transport rates
    current = eyring_script.current_calc(solute_transport)

    # save the fitting results in a custom class
    fitting_specs_eig = FittingMetrics(_test, sum_absolute_errors, sum_squared_errors, solute_transport)

    return Results(voltage, _specs, solute_transport, fitting_specs_eig,
                   current, _ss)


def svd_func(_matrix):
    try:
        U, S, V = np.linalg.svd(_matrix)
        smallest_eig_values = (S[len(S)-1], S[len(S)-2])
        ss_raw = V[len(V)-1, :]
        ss = ss_raw / np.sum(ss_raw)
    except Exception as e:
        print "SVD error", e
        return [False], 0
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
        real_values.append(num.real)  # numerical errors cause the values to become complex for large matrices
    largest_eig_values = two_largest(real_values)  # want the largest numbers because they should all be negative
    _index = real_values.index(max(real_values))
    steady_state_raw = vectors[:, _index].real
    steady_state = steady_state_raw / np.sum(steady_state_raw)
    return steady_state, largest_eig_values


def qr_func(_matrix):
    q, r = np.linalg.qr(_matrix.T)
    return q[:, -1]/np.sum(q[:, -1]), (r[-1, -1], r[-2, -2])


def cond(_matrix):
    return np.linalg.cond(_matrix)


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
    """
    get the smallest and largest non-zero elements of a matrix
    :param _matrix:
    :return:
    """
    abs_matrix = np.fabs(_matrix)
    smallest_element = np.amin(abs_matrix[np.nonzero(abs_matrix)])
    largest_element = np.amax(abs_matrix[np.nonzero(abs_matrix)])
    return largest_element, smallest_element


class MatrixSpecs(object):
    def __init__(self, _largest, _smallest, _condition_number):
        self.largest_element = _largest
        self.smallest_element = _smallest
        self.element_different = _largest - _smallest
        self.condition_number = _condition_number


class FittingMetrics(object):
    def __init__(self, _smallest, _sae, _sse, transport_rates):
        self.smallest_eig = _smallest[1]
        self.second_smallest_eig = _smallest[0]
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
