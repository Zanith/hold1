import rates_algo
import step_algo

__author__ = 'Kyle Vitautas Lopin'


def make_template(_math_type, num_binding_sites, ions, charges, q_type):
    """
    q_type options are: 'single Q', 'single QR', 'full Q', 'full QR'
    :param _math_type:
    :param num_binding_sites:
    :param ions:
    :param charges:
    :param q_type: what type of repulsion / attraction coefficents their are
    :return:
    """
    print "makeing model with ", num_binding_sites, "binding sites and ", ions, " ions"
    with open('template.txt', 'r') as f:
        template = f.read()  # get the template of the file to make
    if _math_type == 'mpmath':
        _import_statement = "from mpmath import mp, exp, matrix, nstr, chop\n" \
                            "import mpmath_helper_function as helper"
        _mp_dps_statement = "mp.dps = mp_dps"
    elif _math_type == 'numpy':
        _import_statement = "import numpy_helper_function as helper\n" \
                            "from numpy import exp, matrix, asscalar, any\n" \
                            "from numpy.linalg import cond"
        _mp_dps_statement = ""
    else:
        raise IOError("_math_type should be 'mpmath or numpy")

    template = template.replace('%%import statement%%', _import_statement)

    template = template.replace('%%mp.dps statement%%', _mp_dps_statement)

    algo_instance = step_algo.EryingRateModelMaker(num_binding_sites, ions, charges, q_type)
    rates_helper_instance = rates_algo.EryingRateMaker(_math_type, num_binding_sites, ions, charges)

    _states = algo_instance.get_states_str()
    template = template.replace('%%states%%', _states)

    ions_str = 'ions = ' + str(ions) + "\nion_charges = {'"
    for ion in ions:
        index = ions.index(ion)
        ions_str += ion + "': " + str(charges[index]) + ", '"
    template = template.replace('%%ions%%', ions_str[:-3] + '}')
    if q_type:
        _global_q = algo_instance.get_q_global_str()
        template = template.replace('%%global Q%%', _global_q)
    else:
        template = template.replace('\n%%global Q%%\n', "")

    _ion_assignment = rates_helper_instance.get_ion_assignment_str()
    template = template.replace('%%ion assignment%%', _ion_assignment)

    if q_type:
        _Q_assignment = algo_instance.get_q_str()
        template = template.replace('%%Q assignment%%', _Q_assignment)
    else:
        template = template.replace('\n%%Q assignment%%\n', "")

    _distance_assignment = rates_helper_instance.get_electrical_distance_str()
    template = template.replace('%%distance assignment%%', _distance_assignment)

    _energy_barriers = rates_helper_instance.get_energy_barrier_str()

    template = template.replace('%%energy barriers%%', _energy_barriers)

    _global_calls = rates_helper_instance.get_global_variables_str()
    template = template.replace('%%global rates%%', _global_calls)

    _ion_rates = rates_helper_instance.get_rates_str()
    template = template.replace('%%ion rates%%', _ion_rates)

    _matrix = algo_instance.get_matrix_str()
    template = template.replace('%%transition matrix%%', _matrix)

    _transport_rates = (algo_instance.get_forward_transport_rates_str()
                        + algo_instance.get_backward_transport_rates_str())
    template = template.replace('%%transport rates%%', _transport_rates)

    return template


if __name__ == '__main__':
    ions = ['Na', 'Ca', 'Mg', 'Fe', 'Ba', 'Cd']
    eyring_rate_script = make_template('numpy', 2, ions[:3], [1, 2, 2, 2, 2, 2], 'single Q')
    # make_template('numpy', 2, ['Na', 'Ca'], [1, 2], None)
    # print eyring_rate_script
    with open("test_script.py", "w") as file:
        file.write(eyring_rate_script)
