import Tkinter as tk


__author__ = 'Kyle Vitautas Lopin'

"""
Notes: store solute
"""

conc_options = ['M', 'mM', 'uM', 'nM', 'pM', 'fM']


class SimulationWindow(tk.Frame):
    """
    Make a frame that allows the user to input a voltage range and ion concentrations
    Use get_run_simulation_settings() to get the values the user entered
    get_settings() is used to redisplay previous values into a new display by putting the values into _settings
    """
    def __init__(self, _master, solutes, _settings=None, *args):
        """
        make new frame by checking if user entered previous values and calling the 2 subroutines:
        make_voltage_selection() and make_solute_concentration()
        :param _master: tkinter frame this frame is embedded in
        :param solutes: list of solutes to make widgets for to allow user to enter
        intra and extracellular concentrations for
        :param _settings: dict of previously used values to initialize the concentrations with
        :param args: just incase
        :return:
        """
        tk.Frame.__init__(self, master=_master, bd=5, relief='raised')
        self.voltage_setting = [-150, 100, 10]  # loading from settings not implemented yet
        if _settings:  # if settings were inputted load those, else make a blank dict to put them in
            self.solute_conc_vars_settings = _settings  # to set the tk variables for the user selected concentrations
        else:
            self.solute_conc_vars_settings = dict()
        self.solute_conc_vars = dict()  # to save tk variables in, values are lists

        # make a section for the user to select the voltage range
        voltage_parameter_frame = tk.Frame(self)
        self.voltage_range = self.make_voltage_selection(voltage_parameter_frame)
        voltage_parameter_frame.pack(side='left')

        # make a section for the user to set the solute concentrations
        self.make_solute_concentration(voltage_parameter_frame, solutes)

    def get_settings(self):
        """
        Get settings user has entered, the is used to redisplay a new sim frame instance
        with the previous' frames values
        :return:  voltages user choose (ex. [-150, 100, 10]) and ion concentration setting values
        (ex. {'Nai': [1.0, 'mM'], 'Nae': [120.0, 'mM']

        """
        # voltages = []
        # for voltage in self.voltage_range:
        #     voltages.append(voltage.get())
        voltages  = [x.get() for x in self.voltage_range]  # get voltages
        ion_concs = dict()
        for key in self.solute_conc_vars:  # get ion concentrations
            ion_concs[key] = []
            ion_concs[key].append(self.solute_conc_vars[key][0].get())
            ion_concs[key].append(self.solute_conc_vars[key][1].get())
        return voltages, ion_concs

    def get_run_simulation_settings(self):
        """
        Send the information needed for the simulation that the user entered
        :return:  list of voltage ranges
        dict of numerical values of ion concentrations used in units of M
        dict of strings with the concentrations the user entered, i.e. '1 nM'
        """
        # get the voltages the user entered
        voltages = []
        for voltage in self.voltage_range:
            voltages.append(voltage.get())
        # make list of the voltage range the user entered
        voltage_range = range(voltages[0], voltages[1]+voltages[2], voltages[2])
        # get the ion concentrations the user enetered
        ion_concs = dict()  # to put the numerical values in Molar, i.e. 0.001 M = 1 mM
        conc_label = dict()  # make string labels to print to the user, i.e. '1mM'
        for key in self.solute_conc_vars:
            _power = 10**(-3*conc_options.index(self.solute_conc_vars[key][1].get()))  # convert mM to 10**-3 etc.
            ion_concs[key] = self.solute_conc_vars[key][0].get() * _power
            conc_label[key] = str(self.solute_conc_vars[key][0].get()) + ' ' + self.solute_conc_vars[key][1].get()
        return voltage_range, ion_concs, conc_label

    def make_voltage_selection(self, _frame):
        """
        enter spinboxes that will allow the user to select the voltage range to simulate
        TODO: clean this up with a subroutine
        :param _frame: the tkinter frame to put the spinboxes in
        :return:  3 instances of tkinter IntVars
        """
        tk.Label(_frame, text="voltage range:").pack(side='top')
        voltage_widget_frame = tk.Frame(_frame)
        voltage_widget_frame.pack(side='top')
        tk.Label(voltage_widget_frame, text="low:").pack(side='left')
        start_voltage = tk.IntVar()
        start_voltage.set(self.voltage_setting[0])
        start_voltage_spinbox = tk.Spinbox(voltage_widget_frame,
                                           textvariable=start_voltage,
                                           from_=-250, to=250, increment=10, width=5)
        start_voltage_spinbox.pack(side='left')
        tk.Label(voltage_widget_frame, text="high:").pack(side='left')
        end_voltage = tk.IntVar()
        end_voltage.set(self.voltage_setting[1])
        end_voltage_spinbox = tk.Spinbox(voltage_widget_frame,
                                         textvariable=end_voltage,
                                         from_=-250, to=250, increment=10, width=5)
        end_voltage_spinbox.pack(side='left')
        tk.Label(voltage_widget_frame, text="increment:").pack(side='left')
        increment_voltage = tk.IntVar()
        increment_voltage.set(self.voltage_setting[2])
        increment_voltage_spinbox = tk.Spinbox(voltage_widget_frame,
                                               textvariable=increment_voltage,
                                               from_=1, to=50, width=5)
        increment_voltage_spinbox.pack(side='left')
        return [start_voltage, end_voltage, increment_voltage]

    def make_solute_concentration(self, _frame, solutes):
        """
        Make a space where the user can enter the intracellular and extracellular concentrations to use
        :param _frame:  frame where the entry and spinboxes will be placed
        :param solutes:  list of the solutes to make entries for
        :return:  nothing, bind the instances to self.solute_conc_vars
        """
        for solute in solutes:
            # for each solute make a frame to put a entry and spinbox for the
            # intracellular and extracellular concentrations
            solute_frame_instance = tk.Frame(_frame, highlightbackground='light slate gray', highlightthickness=2)
            solute_frame_instance.pack(side='top')

            tk.Label(solute_frame_instance, text=solute+" concentrations:").pack(side='top')
            # use a subroutine make the options for the intra and extracellular choices
            # this routine will return a list of an intvar and a StringVar to represent the concentration and units
            self.solute_conc_vars[solute+'i'] = self.make_single_solute_conc(solute_frame_instance,
                                                                             solute, "intra")
            self.solute_conc_vars[solute+'e'] = self.make_single_solute_conc(solute_frame_instance,
                                                                             solute, "extra")

    def make_single_solute_conc(self, _frame, solute, _type):
        """
        Subroutine to make a spinbox and an optionmenu for the user to select the
        amount and units of the concentrations to use
        :param _frame:  frame to put the widgets in
        :param solute:  the name of the solute the user is entering for
        :param _type:  'intra' or 'extra' to tell which side is being used
        :return:  list of the intVar and StringVar the user has to put in the
        amount (IntVar) and the unit (StringVar) the user selected
        """
        # make a frame with an outline to put the widgets in
        intra_frame_instance = tk.Frame(_frame, highlightbackground='LightBlue3', highlightthickness=2)
        intra_frame_instance.pack(side='left')
        # make a label above the widgets
        tk.Label(intra_frame_instance, text=_type+"cellular concentration:").pack(side='top')
        # make frame to put the widgets in
        spinbox_instance_frame = tk.Frame(intra_frame_instance)
        spinbox_instance_frame.pack(side='top')
        # make tkinter variable instances
        intra_var_instance = tk.DoubleVar()
        conc_option_var = tk.StringVar()
        # _type[0] is either 'i' or 'e'
        if solute+_type[0] in self.solute_conc_vars_settings:  # check if the solute has been used before and if
            # so then enter the previous values into the new instances
            previously_entered_values = self.solute_conc_vars_settings[solute+_type[0]]
            previous_conc = previously_entered_values[0]
            previous_order = previously_entered_values[1]
            intra_var_instance.set(previous_conc)
            _option_index = conc_options.index(previous_order)
            conc_option_var.set(conc_options[_option_index])
        else:  # else just set the concentration to 1 mM
            intra_var_instance.set(1)
            conc_option_var.set(conc_options[1])
        # make spinbox the user to set the concentration with
        spinbox_instance_intra = tk.Spinbox(intra_frame_instance,
                                            textvariable=intra_var_instance,
                                            from_=0, to=999, width=7)
        spinbox_instance_intra.pack(side='left')

        # make option menu for use to select concentration units, i.e. mM, uM, nM etc.
        conc_option_menu = tk.OptionMenu(intra_frame_instance, conc_option_var,
                                         *conc_options)
        conc_option_menu.config(width=7)
        conc_option_menu.pack(side='left')
        return [intra_var_instance, conc_option_var]


if __name__ == '__main__':
    app = tk.Tk()
    solutes1 = ['Na', 'Ca', 'Mg']
    apptop = SimulationWindow(app, solutes1, None)
    app.mainloop()