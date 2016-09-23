# Copyright (c) 2015-2016 Kyle Lopin (Naresuan University) <kylel@nu.ac.th>
# Licensed under the GPL
""" Classes and methods to use to make an Eyring Rate Model
"""

__author__ = 'Kyle Vitautas Lopin'


class EryingRateMaker(object):
    """
    Class to make rates to calculate transition rates of an erying rate model
    """
    def __init__(self, _math_package, num_binding_sites, list_ions, ion_charges):
        # make the rates of ions moving forward (towards the intracellular side)
        # and backwards separately
        forward_rates = []
        backward_rates = []

        # initialize and assign variables
        self.transport_rates = dict()
        self.ion_charges = ion_charges
        self.list_ions = list_ions
        self.num_binding_sites = num_binding_sites
        self.rates = []

        # make all the rates that need to be assigned, i.e. k_0_1_
        for i in range(num_binding_sites+1):
            forward_rates.append('k_'+str(i)+'_'+str(i+1)+'_')
            backward_rates.append('k_'+str(i+1)+'_'+str(i)+'_')

        # add ion assignment to the rates, i.e. k_0_1_+'Na'
        for k in range(len(forward_rates)):
            for ion in list_ions:
                self.rates.append(forward_rates[k]+ion)
                self.rates.append(backward_rates[k]+ion)

        # make a copy of the rates to be used to make a global statement in the script
        self.global_rate_names = self.rates[:]

        # make the right side of the rate assignment equations,
        # ie. = mp.mpf(Nae*k0*exp(-GNa1)*exp(1*q*-d1*V))
        self.make_right_side_rate_equation(_math_package)
        self._make_transport_rates()
        self.ion_assign_str = ""
        self._make_ion_assignment(list_ions)
        self.electrical_distance_str = ""
        self._make_elec_dist_str(num_binding_sites)
        self.energy_barrier_str = ""
        self._make_energy_barrier_str(num_binding_sites)

    def make_right_side_rate_equation(self, math_package):
        """
        Make the right side of the ion assignment equations

        ie. = mp.mpf(Nae*k0*exp(-GNa1)*exp(1*q*-d1*V)) for k_0_1_Na

        :param math_package: check if mpmath is being used
        :return: bind all results to self.rates
        """
        # check if mpmath is being used
        if math_package == 'mpmath':
            _start = "mp.mpf("
            _end = "))"
        elif math_package == 'numpy':
            _start = ""
            _end = ")"
        # go through each rate int  self.rates
        for i in range(len(self.rates)):
            # get the binding sides and the ion that is moving, k_0_1_Na => k 0 1 Na
            # and assign them proper place
            rate_elements = self.rates[i].split('_', 3)
            ion = rate_elements[3]
            rate_elements[1] = int(rate_elements[1])  # convert str to int
            rate_elements[2] = int(rate_elements[2])  # convert str to int

            # get the charge of the ion
            charge_index = self.list_ions.index(ion)
            ion_charge = self.ion_charges[charge_index]
            rate_start = ""  # initialize a string to put at the beginning
            # this rate is an ion moving into the channel from extracellular side
            if rate_elements[1] == 0:
                energy_barrier_str = '-G'+ion+'1'
                electrical_distance_str = '-d1'
                rate_start += (ion+"e*")
            elif (rate_elements[1] == self.num_binding_sites+1 and  # the ion is from the last
                  rate_elements[2] == self.num_binding_sites):  # site to the intracellular side
                str_num = str(2*self.num_binding_sites+1)
                energy_barrier_str = '-G'+ion+str_num
                electrical_distance_str = '(1-d'+str_num+')'
                rate_start += (ion+"i*")
            else:
                first_barrier_num = 2*rate_elements[1]
                second_barrier_num = str(first_barrier_num - (rate_elements[1]-rate_elements[2]))
                first_barrier_num = str(first_barrier_num)
                energy_barrier_str = 'G'+ion+first_barrier_num+'-G'+ion+second_barrier_num
                electrical_distance_str = '(d'+first_barrier_num+'-d'+second_barrier_num+')'

            self.rates[i] += ' = ' + _start + rate_start \
                             + 'k0*exp(' + energy_barrier_str + _end \
                             + '*exp(' + str(ion_charge) + '*q*' \
                             + electrical_distance_str + '*V)'

    def _make_transport_rates(self):
        """
        make the rates of how the ions are transported over a barrier
        """
        for ion in self.list_ions:
            self.transport_rates[ion] = []
            for i in range(self.num_binding_sites+1):
                _str = 'k_'+str(i)+'_'+str(i+1)
                self.transport_rates[ion].append([])

    def get_rates_str(self):
        """
        Export a string that assign ion movement rates
        (ie. k_0_1_Na = Nai*ko*exp(-GNa1)*exp(1*q*-d1*V))
        :return:
        """
        return '    ' + '\n    '.join(self.rates)

    def get_global_variables_str(self):
        """
        Get a string to make the rates global
        :return:
        """
        line_cutoff = 70  # make a new line if the current one gets too long
        _str = "    global "
        for rate in self.global_rate_names:
            extra_ending = False
            _str += rate + ', '
            if (_str.find('\n', -line_cutoff) == -1) \
               and (len(_str) > line_cutoff):
                _str = _str[:-2] + '\n    global '
                extra_ending = True
        if extra_ending:
            _str = _str[:-10]

        return _str[:-2]

    def _make_ion_assignment(self, solutes):
        """
        Make a string to assign the concentrations from the input variable ion_conc in the
        Eyring rate script
        :param solutes: list of solutes
        """
        _str = ""
        for ion in solutes:
            _str += '    ' + ion + 'i' + " = ion_concs['" + ion + "i']\n"
            _str += '    ' + ion + 'e' + " = ion_concs['" + ion + "e']\n"
        self.ion_assign_str = _str[:-1]

    def get_ion_assignment_str(self):
        """
        Get the string representing how the solutes should be assigned
        :return: string
        """
        return self.ion_assign_str

    def _make_elec_dist_str(self, num):
        """
        Make a string to assign the electrical distance variables to the inputted dict
        energy_barriers in the Eyring Rate script
        :param num: number of binding sites
        """
        _str = ""
        for i in range(2*num+1):
            _str += '    d' + str(i+1) \
                    + " = energy_barriers['distance'][" \
                    + str(i) + ']\n'
        self.electrical_distance_str = _str[:-1]

    def get_electrical_distance_str(self):
        """
        get the string assigning the electrical distance variables
        :return: string
        """
        return self.electrical_distance_str

    def _make_energy_barrier_str(self, num):
        """
        Make a string to assign the individual energy levels from the inputted dict
        energy_barriers in the Eyring Rate script
        :param num: number of binding sites
        """
        _str = ""
        for ion in self.list_ions:
            for i in range(2*num+1):
                _str += '    G' + ion + str(i+1) + " = energy_barriers['" \
                        + ion + "'][" + str(i) + "]\n"
            _str += '\n'
        self.energy_barrier_str = _str[:-2]

    def get_energy_barrier_str(self):
        """
        Get the string assigning the energy barrier variables
        :return: string
        """
        return self.energy_barrier_str

    def __str__(self):
        _str = "\n".join(self.rates)
        return _str


if __name__ == "__main__":
    RATES = EryingRateMaker('numpy', 2, ['Na', 'Ca'], [1, 2])
    print RATES.get_rates_str()
    print RATES.get_global_variables_str()
