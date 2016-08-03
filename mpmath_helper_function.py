from mpmath import mp, eig, fsum, fabs, norm, qr, mnorm

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
    solute_transport = eyring_script.eyring_rate_transport(_ss)
    # characterize how well the steady state is by getting the residues of the steady state times the transition matrix
    sum_squared_errors, sum_absolute_errors = characterize_solution(_ss, _matrix)

    # 4 calculate the current from the solute transport rates
    current = eyring_script.current_calc(solute_transport)

    # save the fitting results in a custom class
    fitting_specs_eig = FittingMetrics(_test, sum_absolute_errors, sum_squared_errors, solute_transport)

    return Results(voltage, _specs, solute_transport, fitting_specs_eig,
                   current, _ss)


def svd_func(_matrix):
    U, S, V = mp.svd(_matrix)
    smallest_eig_values = (S[len(S)-1], S[len(S)-2])
    ss_raw = V[len(V)-1, :]
    ss = ss_raw / fsum(ss_raw)
    return ss.T, smallest_eig_values


def steady_state_eig(_matrix):
    values, vectors = eig(_matrix)
    real_values = []
    for num in values:
        real_values.append(num.real)
    largest_eig_values = two_largest(real_values)  # want the largest numbers because they should all be negative
    _index = real_values.index(max(real_values))
    steady_state_raw = vectors[:, _index]
    steady_state = steady_state_raw / fsum(steady_state_raw)
    return steady_state, largest_eig_values


def qr_func(_matrix):
    Q, R = qr(_matrix.T)
    ss = Q[:, Q.cols-1]
    return ss/fsum(ss), (R[R.rows-1, R.cols-1], R[R.rows-2, R.cols-2])


def cond(_matrix):
    return mnorm(_matrix, 1)


def characterize_solution(ss, _matrix):
    residues = residues_calc(ss, _matrix)
    sae_residues = norm(residues, 1)
    sse_residues = norm(residues, 2)
    return sae_residues, sse_residues


def residues_calc(null_vector, _matrix):
    return _matrix * null_vector


def two_largest(num_list):
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
    l = len(_matrix) - 1
    smallest_element = mp.mpf(1e100)  # initialize to very large number
    largest_element = mp.mpf(0)
    for i in range(l):
        for j in range(l):
            num = _matrix[i, j]
            if num != 0:
                abs_num = fabs(num)
                if abs_num < smallest_element:
                    smallest_element = abs_num
                if abs_num > largest_element:
                    largest_element = abs_num
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
