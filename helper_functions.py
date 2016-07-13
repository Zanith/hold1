import numpy.linalg
import numpy as np

# np.set_printoptions(threshold=np.nan)

__author__ = 'Kyle Vitautas Lopin'


def steady_state_lowest_eig(_matrix, _min_value=1e-10):
    values, vectors = numpy.linalg.eig(_matrix)
    print "eigen values: "
    print values
    print "eigen vectors: "
    print vectors
    min_index = np.argmin(np.absolute(values))

    min_eig = vectors[:, min_index]

    print "make anything less than min_value 0 here"

    ss_eigen = np.divide(min_eig, np.sum(min_eig))

    return ss_eigen


def steady_state_null_space(_matrix, eps=1e-4):
    null_matrix = null(_matrix, eps)
    return np.divide(null_matrix, np.sum(null_matrix))


def null(_matrix, eps):
    u, s, vh = np.linalg.svd(_matrix)
    print 's:'
    print s
    print 'vh'
    print vh
    print 'test'
    print vh[-1]
    null_mask = (s <= eps)
    print 'null mask'
    print null_mask
    null_space = np.compress(null_mask, vh, axis=0)
    print "null_space"
    print null_space
    return np.transpose(null_space)


def matrix_solve(_matrix):
    _b = np.zeros((len(_matrix), 1))
    print "check2"
    print _matrix
    print _b
    ss_solved = numpy.linalg.solve(_matrix, _b)
    print ss_solved


def check_soln(_matrix, _steady_state):
    individual_errors = np.dot(_matrix, _steady_state)
    print individual_errors
    get_errors(_matrix, _steady_state, -1)
    a = .02
    step = _steady_state

    for i in range(2000):
        _steady_state = _steady_state - np.multiply(a, _steady_state)
        _steady_state = np.divide(_steady_state, np.sum(_steady_state))
        get_errors(_matrix, _steady_state, i)
    print _steady_state


def get_errors(_matrix, _steady_state, _i):
    error = np.dot(_matrix, _steady_state)
    sum_error = np.sum(np.absolute(error))

    print "i: ", _i, sum_error


def steady_state_populations(steady_state, states):
    # state_index = np.argsort(steady_state, axis=0)
    print np.ravel(steady_state)
    state_index = np.argsort(np.ravel(steady_state))  # get the indices of
    print "steady state, ", steady_state.shape
    print steady_state[:, 0]
    ordered_steady_state, indexes = np.unique(np.ravel(steady_state),
                                              return_index=True)
    ordered_states = []
    for _i in indexes:
        ordered_states.append(states[_i])
    ordered_steady_state = ordered_steady_state[::-1]
    ordered_states = ordered_states[::-1]
    print 'reversed'
    print ordered_steady_state
    print 'states'
    print ordered_states
    ordered_pairs = []
    for i in range(len(ordered_states)):
        ordered_pairs.append((ordered_steady_state[i], ordered_states[i]))

    for j in range(len(ordered_pairs)):
        print ordered_pairs[j]


def compare_two_states(vector1, vector2, tran_matrix):
    _compare_list = []
    print "comparing vectors", vector1[0][0]
    for i in range(len(vector1)):
        _compare_list.append([vector1[i][0], vector2[i][0]])
        print [vector1[i], vector2[i]]
    error1 = np.dot(tran_matrix, vector1)
    print "error1:"
    print error1
    sum_error1 = np.sum(np.absolute(error1))
    error2 = np.dot(tran_matrix, vector2)
    print "error2:"
    print error2
    sum_error2 = np.sum(np.absolute(error2))
    print "errors: ", sum_error1, sum_error2
    return sum_error1, sum_error2
