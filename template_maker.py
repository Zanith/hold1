# Copyright (c) 2015-2016 Kyle Lopin (Naresuan University) <kylel@nu.ac.th>
# Licensed under the GPL

""" function to make file that solves an eyring rate model
"""

# local files
import rates_algo
import step_algo

__author__ = 'Kyle Vitautas Lopin'


def make_template(_math_type, num_binding_sites, solutes, charges, q_type):
    """
    Make a string that can be saved as a file to solve an eyring rate model
    :param _math_type:  string of type of math package to use, 'numpy' or 'mpmath'
    :param num_binding_sites:  int of the number of binding sites in the model
    :param solutes:  list of strings of the solutes in the model
    :param charges:  list of int of the charges of the solutes
    :param q_type: what type of repulsion / attraction coefficents their are;
    options are: 'single Q', 'single QR', 'full Q', 'full QR'
    :return: string that is a python file to run an eyring rate model
    """
    with open('template.txt', 'r') as temp_file:
        template = temp_file.read()  # get the template of the file to make
    # set the correct imports for the selected math type
    if _math_type == 'mpmath':
        _import_statement = "from mpmath import mp, exp, matrix, nstr, chop\n" \
                            "import mp_func as helper"
        _mp_dps_statement = "mp.dps = mp_dps"
    elif _math_type == 'numpy':
        _import_statement = "import np_func as helper\n" \
                            "from numpy import exp, matrix\n"
        _mp_dps_statement = ""
    else:
        raise IOError("_math_type should be 'mpmath' or 'numpy'")

    # put in the import statements
    template = template.replace('%%import statement%%', _import_statement)
    # set the precision levels if mpmath is being used
    template = template.replace('%%mp.dps statement%%', _mp_dps_statement)

    # create instances of the 2 helper modules needed to make the eyring rate model
    algo_instance = step_algo.EryingRateModelMaker(num_binding_sites, solutes, charges, q_type)
    rates_helper_instance = rates_algo.EryingRateMaker(_math_type, num_binding_sites, solutes, charges)

    # put in a list of the states that are possible for the model
    _states = algo_instance.get_states_str()
    template = template.replace('%%states%%', _states)

    # make a list of the names of the solutes in the model and a separate dict of charges with
    # the solutes as keys and charge as values an put in the top of the eyring rate file
    ions_str = 'ions = ' + str(solutes) + "\nion_charges = {'"
    for ion in solutes:
        index = solutes.index(ion)
        ions_str += ion + "': " + str(charges[index]) + ", '"
    template = template.replace('%%ions%%', ions_str[:-3] + '}')
    # check if Q or R values are used and add statements to make them global if they are
    if q_type:
        _global_q = algo_instance.get_qr_global_str()
        template = template.replace('%%global Q%%', _global_q)
    else:
        template = template.replace('\n%%global Q%%\n', "")

    # the energy barriers are inputted as a dict, make a series of assignments
    # to assign the individual values to the correct dict value
    _ion_assignment = rates_helper_instance.get_ion_assignment_str()
    template = template.replace('%%ion assignment%%', _ion_assignment)

    # make statement to assign Q and R values
    if q_type:
        _Q_assignment = algo_instance.get_q_str()
        template = template.replace('%%Q assignment%%', _Q_assignment)
    else:
        template = template.replace('\n%%Q assignment%%\n', "")

    # assign the distance values to the input argument
    _distance_assignment = rates_helper_instance.get_electrical_distance_str()
    template = template.replace('%%distance assignment%%', _distance_assignment)

    # assign the energy barriers variables to the correct input argument
    _energy_barriers = rates_helper_instance.get_energy_barrier_str()
    template = template.replace('%%energy barriers%%', _energy_barriers)

    # make all the rates global so that all the function in the script can use them
    _global_calls = rates_helper_instance.get_global_variables_str()
    template = template.replace('%%global rates%%', _global_calls)

    # calculate the equations to make the rates,
    # ex. k_0_1_solute_1 = solute_1e*k0*exp(-Gsolute_11)*exp(0*q*-d1*V)
    _ion_rates = rates_helper_instance.get_rates_str()
    template = template.replace('%%ion rates%%', _ion_rates)

    # make and put in the transition matrix
    _matrix = algo_instance.get_matrix_str()
    template = template.replace('%%transition matrix%%', _matrix)

    # make the equation to calculate the rate of solutes moving over a barrier
    _transport_rates = (algo_instance.get_forward_transport_rates_str()
                        + algo_instance.get_backward_transport_rates_str())
    template = template.replace('%%transport rates%%', _transport_rates)

    return template


if __name__ == '__main__':
    _ions = ['Na', 'Ca', 'Mg', 'Fe', 'Ba', 'Cd']
    EYRING_RATE_SCRIPT = make_template('numpy', 2, _ions[:3], [1, 2, 2, 2, 2, 2], 'single Q')
    print make_template('numpy', 2, ['Na', 'Ca'], [1, 2], None)
    # print eyring_rate_script
    # with open("test_script.py", "w") as _file:
    #     _file.write(EYRING_RATE_SCRIPT)
