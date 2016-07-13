from mpmath import mp, eig, fsum, fabs, matrix, inverse, norm, nprint, fp


__author__ = 'Kyle Vitautas Lopin'


def svd_func(_matrix):
    U, S, V = mp.svd(_matrix)
    smallest_eig_values = (S[len(S)-1], S[len(S)-2])
    ss_raw = V[len(V)-1, :]
    ss = ss_raw / fsum(ss_raw)
    return ss.T, smallest_eig_values


def steady_state_eig(_matrix):
    print 'test c'
    values, vectors = eig(_matrix)
    print 'test d'
    real_values = []
    for num in values:
        real_values.append(num.real)
    largest_eig_values = two_largest(real_values)  # want the largest numbers because they should all be negative
    _index = real_values.index(max(real_values))
    steady_state_raw = vectors[:, _index]
    steady_state = steady_state_raw / fsum(steady_state_raw)
    return steady_state, largest_eig_values


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
    def __init__(self, _largest, _smallest):
        self.largest_element = _largest
        self.smallest_element = _smallest
        self.element_different = _largest - _smallest


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
