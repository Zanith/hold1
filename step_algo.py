__author__ = 'Kyle Vitautas Lopin'

power_symbol = '**'


class QClass(object):
    def __init__(self, _q, _index):
        self.place = str(_index + 1)
        self.charge = _q


class TransportClass(object):
    """
    A custom class to be used by EryingRateModelMaker to keep track of how to calculate
    the transport rates of ions over different barriers.
    note: each instance of a class calculates for all barriers and all ions
    but only in one direction.  Either the forward direction (eg towards the intracellular
    side) or the backward direction (eg towards the extracellular side).
    """
    def __init__(self, ions, number_of_barriers):
        """
        Initialize the class by make the proper data structure.  The basis of the class is a list
        where the ith entry is the rate for the (i+1) barrier if the first barrier is counted as 1
        Each element in the base list is a dictionary where the keys of the dictionary is the string of the
        ion being transported and the value of the dictionary is a list that contains the rates of ion
        movement over a barrier.

        To access a value go to self.list[number of the barrier - 1][name of the ion]

        example:
        [ {'Na': [('k_0_1_Na', 1), ('k_0_1_Na', 3)]},
          {'Ca': [('k_0_1_Ca', 2), ('k_0_1_Ca', 4)]}
        """
        self.list = []
        for barrier_number in range(1, number_of_barriers+1):
            self.list.append(dict())
            for ion in ions:
                self.list[barrier_number-1][ion] = []

    def update(self, barrier, ion, rate_str, state_index):
        self.list[barrier-1][ion].append((rate_str, state_index))

    def display(self):
        barrier_number = 0
        for barrier in self.list:
            barrier_number += 1
            for ion in barrier:
                print 'ion: '
                print ion
                print 'value: '
                print barrier[ion]

    def output_rate(self, barrier_number, ion, end_string=""):
        """
        Output a string that takes all the rates in the instance and makes a mathematical
        expression to calculate transport rate of in ion (ion) over a barrier (barrier_number)
        Note, each instance only contains one direction so this is only for forward, or backward
        rates
        """
        rates_states = self.list[barrier_number-1][ion]
        _str = ""
        for rates, states in rates_states:
            _str += rates + ' * steady_state[' + str(states) + '] +\n' + end_string
        return _str[:(-3-len(end_string))]


class EryingRateModelMaker(object):
    """
    Class to create a transition matrix to describe an Erying Rate model with an
    arbitrary amount of binding sites and solutes that may bind at the sites
    """
    def __init__(self, num_binding_sites, list_ions, ion_charges, q_type=None):
        """
        Take an arbitrary number and a list of solutes that can interact with the channel and create
        a transition matrix of the different states and rates that can occur
        :param num_binding_sites: integer, number of binding sites in the channel
        :param list_ions: list of strings that represent the ions
        :param ion_charges: list of charges the corresponds to the list of ions
        :param q_type: what time of Q values to use for the model,
        'single Q' puts a Q when 2 ions are adjacent
        'single QR' puts a Q or a R when 2 ions are adjacent (Q) or moving to be adjacent
        'full Q' puts a Q for each ion-ion interaction
        'full QR' puts Q's and Rs for each ion-ion interact
        :return: makes a matrix bound to self that represents the transition matrix of the form
        X' = A*X, where X is a vector of possible channel configurations and A is matrix that describes how
        the channel configurations change
        also make a vector (channel_configs) that represents the different configurations
        (states) the channel can take
        """
        # initialize the class attributes
        self.num_binding_sites = num_binding_sites
        self.list_ions = list_ions
        self.ion_charges = ion_charges
        self.q_type = q_type
        self.q_str = ""
        self.q_global_list = []
        if q_type:
            self.make_Q_assignment()
        self.channel_configs = [[0] * num_binding_sites]
        self.new_channel_configs = self.channel_configs[:]  # prevent aliasing
        self.matrix = [['0']]
        # to make the equations for the transport equation
        self.forward_transport = TransportClass(list_ions, num_binding_sites+1)
        self.backward_transport = TransportClass(list_ions, num_binding_sites+1)
        # go through all channel configurations and check if they can transition to another state
        while self.new_channel_configs:
            # get the first config in the new_channel_config list and see if that state can transition to another state
            channel_state = self.new_channel_configs.pop(0)
            # keep track of what the original state is
            original_state_index = self.channel_configs.index(channel_state)
            # check if an ion can enter the channel configuration from the extracellular side
            self.check_for_ion_entry(channel_state, original_state_index)
            # check if an ions can move forward in the pore
            self.check_for_hop(channel_state, original_state_index)
            # check if an ion can exit the channel to the intracellular side
            self.check_for_ion_exit(channel_state, original_state_index)

        self.add_diagonal_rates()
        # self.forward_transport.output_rate(1, "Na")
        # self.backward_transport.output_rate(1, "Na")

    def check_for_ion_entry(self, config, original_state_index):
        """
        Check if a channel state can change if an ion can enters it and update
        the channel_configs and matrix rates if an ion can enter
        :param config: state configuration to test
        :param original_state_index: integer representing the index where _config is in self.channel_configs
        :return: none, update self.channel_configs and self.matrix if possible
        """
        # check if the channel config is empty in the first binding site
        if config[0] == 0:
            # if there is an empty site so go through each ion and make update states and rates
            for ion in self.list_ions:
                _new_config = config[:]  # prevent aliasing
                _new_config[0] = ion  # add the ion to the existing channel state
                sites = (0, 1)  # the sites of ion movement are from the extracellular side (0) to the first
                # binding site (1)
                if _new_config not in self.channel_configs:
                    # if the state has not been seen before, add the state and rates to the channel_config and matrix
                    self.add_new_state(_new_config, sites, ion, original_state_index)
                else:
                    # else if the state has already been added, just update the rate matrix
                    new_state_index = self.channel_configs.index(_new_config)
                    self.add_rates(sites, ion, original_state_index, new_state_index)

    def check_for_hop(self, config, original_state_index):
        """
        Check if it is possible for an ion hopping to occur in the channel config, it it
        is possible call the add_hop method to perform the action
        :param config: channel configuration
        :return: none, the channel_config and matrix attributes are updated
        """
        # check if ion can hop inside the channel by going through each binding site (except the last, that
        # situation is dealt with in the check_for_ion_exit method)
        for site_index in range(self.num_binding_sites-1):
            # check if there is an ion present and a place open for it to move to
            if (config[site_index] != 0) and (config[site_index+1] == 0):
                _config = config[:]  # copy the state so it is not aliased
                # it is possible for an ion to move, call the add_hop sub routine to peform the action
                self.add_hop(_config, site_index, original_state_index)

    def check_for_ion_exit(self, config, original_state_index):
        """
        check if it is possible for an ion to exit the channel to the intracellular side (right side of the channel
        configuration list)
        :param config: state configuration to test
        :param original_state_index: integer representing the index where _config is in self.channel_configs
        :return: none, update self.channel_configs and self.matrix if possible
        """
        # check if an ion is in the binding site next to the intracelluar side and can move out of the channel
        if config[-1] != 0:
            _new_config = config[:]  # prevent aliasing by make a new deep copy
            ion = config[-1]  # get the ion that can move out so that the proper rate can be made
            _new_config[-1] = 0  # make the new configuration state where the ion has exited
            sites = (self.num_binding_sites, self.num_binding_sites+1)  # the ion moves from the last binding site
            # (indexed as num_binding_sites) to the extracellular side (indexed as num_binding_sites + 1)
            if _new_config not in self.channel_configs:
                # if the state has not been seen before, add the state and rates to the channel_config and matrix
                self.add_new_state(_new_config, sites, ion, original_state_index)
            else:
                # else if the state has already been added, just update the rate matrix
                new_state_index = self.channel_configs.index(_new_config)
                self.add_rates(sites, ion, original_state_index, new_state_index)

    def add_hop(self, _config, binding_site, original_state_index):
        """
        Take a channel configuration (_config) and move the ion in the binding_site toward the intracellular
        side, add the corresponding config to the channel_config attribute if it is not already present,
        add the rates to the matrix
        :param _config: channel configuration that can allow an ion to hop towards the intracellular side
        :param binding_site: number of the binding site that has the ion that is to move
        :param original_state_index: number that is index of where _config is in the self.channel_configs list
        :return: none, update channel_config and matrix attributes
        """
        ion = _config[binding_site]
        _config[binding_site+1] = _config[binding_site]
        _config[binding_site] = 0
        sites = (binding_site+1, binding_site+2)
        if _config not in self.channel_configs:
            # the first binding site is labeled as 1, as the extracellular
            # side is 0 so increment the binding_site to reflect this
            self.add_new_state(_config, sites, ion, original_state_index)
        else:  # add the rates in here
            new_state_index = self.channel_configs.index(_config)
            self.add_rates(sites, ion, original_state_index, new_state_index)

    def add_new_state(self, _config, sites, ion, original_state_index):
        """
        Add a new ion channel state to the list of ion channel configurations with
        a new way of distributing the ions in the channel in the binding sites
        :param _config: new ion configuration to add to channel_configs list
        :param sites: a tuple with two numbers, the first is the binding site number
        where the ion starts (or 0 if is extracellular ion entering channel),
        and the second number is the binding site the ion moves to
        :param ion:  string with the name of the ion that is moving
        :param original_state_index:  index of the place the state the new state
        was derived from
        :return:
        """
        self.channel_configs.append(_config)  # add the configuration to the list of all configs
        if _config not in self.new_channel_configs:
            self.new_channel_configs.append(_config)  # add the config to list to check for further transitions
        for i in range(len(self.matrix)):  # add another column to the transition matrix
            self.matrix[i].append('0')
        # add a new row to the bottom of the matrix
        self.matrix.append(['0'] * (len(self.matrix)+1))

        self.add_rates(sites, ion, original_state_index, len(self.channel_configs)-1)

    def add_rates(self, sites, ion, original_state_index, new_config_state_index):
        """
        Add rates to the transition matrix
        :param sites: a tuple with two numbers, the first is the binding site number
        where the ion starts (or 0 if is extracellular ion entering channel),
        and the second number is the binding site the ion moves to
        :param ion: string with the name of the ion that is moving
        :param original_state_index:  index of the place the state the new state
        :param new_config_state_index: index of the channel configuration that is being transitioned to
        :return:  none, the matrix is updated
        """
        _config = self.channel_configs[original_state_index]
        # make a string that contains the pre-exponential Q values
        if self.q_type:
            ion_charge = self.ion_charges[self.list_ions.index(ion)]
            q_str_forward, q_str_backward = self.check_for_Q_values(_config, sites[0], ion_charge)
        else:
            q_str_forward, q_str_backward = "", ""
        # make a string to describe the rate for the ion moving to the extracellular side
        backward_rate_str = q_str_backward+'k_'+str(sites[1])+'_'+str(sites[0])+'_'+ion

        if self.matrix[original_state_index][new_config_state_index] == '0':
            self.matrix[original_state_index][new_config_state_index] = backward_rate_str
        else:
            self.matrix[original_state_index][new_config_state_index] += (' + ' + backward_rate_str)
        # update the transport rates
        self.backward_transport.update(sites[1], ion, backward_rate_str, new_config_state_index)

        # make a string to describe the rate for the ion moving to the intracellular side
        forward_rate_str = q_str_forward+'k_'+str(sites[0])+'_'+str(sites[1])+'_'+ion
        if self.matrix[new_config_state_index][original_state_index] == '0':
            self.matrix[new_config_state_index][original_state_index] = forward_rate_str
        else:
            self.matrix[new_config_state_index][original_state_index] += (' + ' + forward_rate_str)
        # update the transport rates
        self.forward_transport.update(sites[1], ion, forward_rate_str, original_state_index)

    def add_diagonal_rates(self):
        """
        Add the rates along the diagonal of the transition matrix by summing all numbers in
        each column and subtracting from the diagonal element
        :return: none, update self.matrix
        """
        # for each element in the column
        for i in range(len(self.matrix)):
            # make a list of elements in the column to make the rates with at the end
            _new_rate_list = []
            # go through each row in the column
            for j in range(len(self.matrix[i])):
                # if there an element there, add it to the list
                if self.matrix[j][i] != '0':
                    _new_rate_list.append(self.matrix[j][i])
            # make a string by subtracting all the elements in the column
            _new_rate = '-(' + " + ".join(_new_rate_list) + ')'
            # and add to the diagonal element of that column
            self.matrix[i][i] = _new_rate

    def check_for_Q_values(self, _config, site_index, ion_moving_charge):
        """
        Make a string for the pre exponential Q or R values that account for the increase in the rate
        of ion movement due to electrostatic repulsion
        :param _config: configuration of the channel to find the Q values for
        :param site_index: binding site the ion is currently in
        :return: string that has the q
        """
        return_str_forward = ""  # initialize string to save Q or R values to
        return_str_backward = ""
        forward_rate = []
        backward_rate = []
        left_ions = []  # make list of ions that are in sites towards the extracellular side of the current ion
        right_ions = []  # make list of ions that are in sites towards the intracellular side of the current ion
        for left_site_index in range(0, site_index-1):
            if _config[left_site_index] != 0:  # there is an ion to the left of the hopping ion
                _ion_in_channel = _config[left_site_index]
                ion_charge_list_index = self.list_ions.index(_ion_in_channel)
                _charge = self.ion_charges[ion_charge_list_index]
                left_ions.append(QClass(_charge, left_site_index))
        for right_site_index in range(site_index+1, self.num_binding_sites):
            if _config[right_site_index] != 0:
                _ion_in_channel = _config[right_site_index]
                ion_charge_list_index = self.list_ions.index(_ion_in_channel)
                _charge = self.ion_charges[ion_charge_list_index]
                right_ions.append(QClass(_charge, right_site_index))
        for left_ion in left_ions:
            if left_ion.charge*ion_moving_charge == 1:
                power_str = ""
            else:
                power_str = power_symbol + str(left_ion.charge*ion_moving_charge)
            if 'full Q' in self.q_type:
                forward_rate.append('Q'
                                    + str(left_ion.place)
                                    + str(site_index)
                                    + power_str)
            elif 'single Q' in self.q_type and int(left_ion.place)+1 == site_index:
                forward_rate.append('Q' + power_str)

            if 'full QR' == self.q_type:
                backward_rate.append('R'
                                     + str(left_ion.place)
                                     + str(site_index+1)
                                     + power_str)
            elif 'single QR' == self.q_type and left_ion.place+2 == site_index:
                backward_rate.append('R' + power_str)
        for right_ion in right_ions:
            if right_ion.charge*ion_moving_charge == 1:
                power_str = ""
            else:
                power_str = power_symbol + str(right_ion.charge*ion_moving_charge)
            if 'full QR' == self.q_type:
                forward_rate.append('R'
                                    + str(site_index)
                                    + str(right_ion.place)
                                    + power_str)
            elif 'single QR' in self.q_type and int(right_ion.place)-1 == site_index:
                forward_rate.append('R' + power_str)

            if ('full Q' in self.q_type):
                backward_rate.append('Q'
                                     + str(site_index+1)
                                     + str(right_ion.place)
                                     + power_str)
            elif 'single Q' in self.q_type and int(right_ion.place)-2 == site_index:
                backward_rate.append('Q' + power_str)

        if forward_rate:
            return_str_forward += " * ".join(forward_rate) + " * "
        if backward_rate:
            return_str_backward += " * ".join(backward_rate) + " * "
        return return_str_forward, return_str_backward

    def make_Q_assignment(self):
        _str = ""
        q_global_list = []  # make list of Qs to put in global call
        q_list_index = 0
        r_list_index = 0
        if 'full' in self.q_type:
            num = self.num_binding_sites
            if 'Q' in self.q_type:  # the full Q values are added in
                for i in range(1, num):
                    for j in range(i+1, num+1):
                        _q_str = 'Q' + str(i) + str(j)
                        q_global_list.append(_q_str)
                        _str += '    ' + _q_str + " = Qs[" + str(q_list_index) + ']\n'
                        q_list_index += 1
            if 'R' in self.q_type:  # full R values should be added
                for i in range(num-1):
                    for j in range(i+2, num+1):
                        _r_str = 'R' + str(i) + str(j)
                        q_global_list.append(_r_str)
                        _str += '    ' + _r_str + " = Rs[" + str(r_list_index) + ']\n'
                        r_list_index += 1
                for k in range(1, num):
                    _r_str = 'R' + str(k) + str(num+1)
                    q_global_list.append(_r_str)
                    _str += '    ' + _r_str + " = Rs[" + str(r_list_index) + ']\n'
                    r_list_index += 1
        elif 'single' in self.q_type:
            if 'Q' in self.q_type:
                _str += '    Q = Qs[0]\n'
            if 'R' in self.q_type:
                _str += '    R = Rs[0]\n'
        self.q_str = _str[:-1]
        self.q_global_list = q_global_list

    def get_transport_rate(self, barrier_number, ion):
        forward_eqn = self.forward_transport.output_rate(barrier_number, ion)
        backward_eqn = self.backward_transport.output_rate(barrier_number, ion)
        return forward_eqn, backward_eqn

    def get_forward_transport_rates_str(self):
        _str = "    inward = dict()\n"
        zeros_string = '0, ' * (self.num_binding_sites + 1)  # +1 because there is one more barrier than binding site
        init_string = "'] = [" + zeros_string[:-2] + "]\n"
        for ion in self.list_ions:
            _str += "    inward['" + ion + init_string
            filler_str = " " * (21 + len(ion))

            for i in range(self.num_binding_sites+1):
                _str += "    inward['" + ion + "'][" + str(i) + '] = ('
                _str += self.forward_transport.output_rate(i+1, ion, filler_str)
                _str += ')\n'
        return _str

    def get_backward_transport_rates_str(self):
        _str = "    outward = dict()\n"
        zeros_string = '0, ' * (self.num_binding_sites + 1)  # +1 because there is one more barrier than binding site
        init_string = "'] = [" + zeros_string[:-2] + "]\n"
        for ion in self.list_ions:
            _str += "    outward['" + ion + init_string
            filler_str = " " * (22 + len(ion))
            for i in range(self.num_binding_sites+1):
                _str += "    outward['" + ion + "'][" + str(i) + '] = ('
                _str += self.backward_transport.output_rate(i+1, ion, filler_str)
                _str += ')\n'
        return _str

    def get_transition_matrix(self):
        return self.matrix

    def get_states_vector(self):
        return self.channel_configs

    def get_states_str(self):
        _str = "states = ["
        for state in self.channel_configs:
            _str += str(state) + ',\n          '
        _str = _str[:-12] + ']'
        return _str

    def get_matrix_str(self):
        output_str = 'matrix([\n    '
        for row in self.matrix:
            output_str += str(row) + ',\n    '
        output_str = output_str.translate(None, "'")
        output_str = output_str[:-6] + '\n    ])'
        return output_str

    def get_number_states(self):
        return len(self.channel_configs)

    def get_q_str(self):
        return self.q_str

    def get_q_global_str(self):
        if self.num_binding_sites < 2:  # if there is only 1 binding site, there are not Qs
            return ""
        if 'full Q' in self.q_type:
            _str = "    global "
            _str += ', '.join(self.q_global_list)
        elif 'single Q' in self.q_type:
            _str = "    global Q  "
        return _str

    def print_matrix(self):
        print 'matrix:'
        for row in self.matrix:
            print row

    def make_chain(self):
        """
        Make a dictionary of the transition states and rates that can change each state
        the values of chain is another dictionary with each state and rate than can transition to
        the key of chain

        chain = { '0' : { 'Na' : rate_of_Na_exiting },
                        { 'Ca' : rate_of_Ca_exiting } }
        :return: chain of the states and rates that can reach each state
        """
        chain = {}
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[0])):
                if self.matrix[i][j] != 0:
                    # add states and rate to a chain
                    if tuple(self.channel_configs[i]) in chain:
                        chain[tuple(self.channel_configs[i])][tuple(self.channel_configs[j])] = self.matrix[i][j]
                    else:
                        chain[tuple(self.channel_configs[i])] = {tuple(self.channel_configs[j]): self.matrix[i][j]}
                    # chain[(tuple(self.channel_configs[i]), tuple(self.channel_configs[j]))] = self.matrix[i][j]
        return chain

    def __str__(self):
        """
        represent the states and matrix in a reader friendly version
        :return:
        """
        states_str = "states = \n["
        for state_row in self.channel_configs:
            states_str += str(state_row) + ',\n'
        states_str = states_str[:-2] + ']\n'

        matrix_str = "matrix([\n"
        for row in self.matrix:
            matrix_str += '['
            matrix_str += ', '.join(row)
            matrix_str += '],\n'
        matrix_str = matrix_str[:-2] + '\n])'
        return states_str + matrix_str.translate(None, "'")


def print_chain(_chain, _configs):
    for config in _configs:
        for key, value in _chain[tuple(config)].iteritems():
            print key, value
    print len(_chain)


if __name__ == "__main__":
    instance = EryingRateModelMaker(2, ['Na', 'Ca'], [1, 2, 2], 'full QR')
    # print instance.get_states_str()
    # print instance.get_matrix_str()
    if False:
        forward_trans_eqn, backward_trans_eqn = instance.get_transport_rate(1, 'Ca')
        print forward_trans_eqn
        print 'split'
        print backward_trans_eqn
        print 'break'
        forward_trans_eqn, backward_trans_eqn = instance.get_transport_rate(1, 'Na')
        print forward_trans_eqn
        print 'split'
        print backward_trans_eqn
        print 'break1'
        forward_trans_eqn, backward_trans_eqn = instance.get_transport_rate(0, 'Na')
        print forward_trans_eqn
        print 'split'
        print backward_trans_eqn

    hold_f = instance.get_forward_transport_rates_str()
    hold_b = instance.get_backward_transport_rates_str()
    print instance.get_q_global_str()
    print ""
    print hold_f[:-1]
    print ""
    print hold_b[:-1]
    # forward_trans_eqn, backward_trans_eqn = instance.get_transport_rate(2, 'Na')
    # print forward_trans_eqn
