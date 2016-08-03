__author__ = 'Kyle Vitautas Lopin'

import Tkinter as tk
import tkFileDialog

import gui_parameter_frame as param_frame
import gui_simulation_frame as sim_frame
import template_maker as script_maker
import eyring_rate_script
import result_toplevel


max_charge = 5
FILEOPTIONS = {'defaultextension': '*.ert_model'}


class EyringGUI(tk.Tk):
    """
    A graphical user interface to make a multiple binding site Eyring rate model and to view the
    resulting transport rates
    """
    def __init__(self, parent=None):
        """
        Make the initiale frame for the user to select the number of binding sites and the solutes to model.
        The user also enters the charge of each solute and selects if repulsion factor should be used.  The user
        can also select to use the mpmath package to solve the problem with
        """
        tk.Tk.__init__(self, parent)
        # initialize variables to be used later
        self.Q_value = None  # not sure why but these need to be passed this way
        self.R_value = None
        self.parameter_frame = None
        self.simulation_frame = None
        self.run_button = None

        self.top_frame = tk.Frame(self)
        self.top_frame.pack(side='top')
        self.bottom_button_frame = tk.Frame(self)
        self.bottom_button_frame.pack(side='top')
        # frame to pack run simulation button in after initialization
        self.run_button_frame = tk.Frame(self.bottom_button_frame)
        self.run_button_frame.pack(side='top')

        self.initial_frame = InitialFrame(self.top_frame)  # hacked it to self to make menu work

        # make the Menu bar
        self.make_menubar()

        self.hold_parameter_frame = tk.Frame(self.top_frame)

        self.initial_frame.pack(side='left', anchor='n')
        make_params_button = tk.Button(self.initial_frame, text="Set model parameters",
                                       command=lambda: self.process_parameters(self.initial_frame))
        make_params_button.pack(side='bottom')

    def process_parameters(self, initial_settings, _saved_settings=None):
        """
        Process information to properly set the initial settings
        """
        if initial_settings.R_value_used:
            R_str = "R"
        else:
            R_str = ""

        num_binding_sites = initial_settings.num_binding_sites.get()

        solutes, charges = self.get_settings(initial_settings)

        print 'making script with num binding sties= ', num_binding_sites
        QR_str = "single Q" + R_str

        # make the script for the Eyring rate model
        print 'mpmath used = ', initial_settings.mpmath_used
        if initial_settings.mpmath_used:
            math_package = 'mpmath'
        else:
            math_package = 'numpy'
        print 'math package: ', math_package
        _eyring_rate_script = script_maker.make_template(math_package, num_binding_sites,
                                                         solutes, charges, QR_str)
        with open("eyring_rate_script.py", "w") as _file:
            _file.write(_eyring_rate_script)
        _file.close()
        reload(eyring_rate_script)

        self.set_parameters(solutes, initial_settings.charges,
                            initial_settings.num_binding_sites.get(), _saved_settings)

    def set_parameters(self, solute_text, charges, num_binding_sites, _saved_settings=None):
        """
        set all the parameters for the parameters frame and the simulation frame
        :param solute_text:  list of solutes the model is made with
        :param charges: list of the charges each solute has
        :param num_binding_sites: number of binding sites the model has
        :param _saved_settings:  dict of settings to pass to the simulation and params frame, used by the load function
        should have a 'param' key with setting for param frame,
        e.g.{'distance': [0.25, 0.5, 0.75], 'Na': [8.0, -10.0, 8.0], 'Mg': [8.0, -10.0, 8.0]}
        and a 'sim' key with setting for the simulation frame, e.g.
        {'sim': {'Mge': ['0.0', 'mM'], 'Nae': ['145.0', 'mM'], 'Mgi': ['1.0', 'mM'], 'Nai': ['145.0', 'mM']}
        """
        print 'settings in set: ', _saved_settings
        # check if the user has inputted any saved parameters, or if there is are already parameters in the frame
        saved_params = None
        if _saved_settings:
            saved_params = _saved_settings['param']  # load the users parameters
        elif self.parameter_frame:
            # take the old parameters of the frame and pass them to the new frame to make it with
            saved_params = self.parameter_frame.energy_barriers
        if self.parameter_frame:  # if the there is already a parameter frame, destroy it so it can be replaced
            self.parameter_frame.destroy()
        # have a frame to hold the position of the parameter frome so it doesn't move
        self.hold_parameter_frame.pack(side='left')

        # make a frame for the user to make the energy profiles to use
        self.parameter_frame = param_frame.ParameterSelectionFrame(self.hold_parameter_frame,
                                                                   solute_text, charges,
                                                                   num_binding_sites,
                                                                   self.Q_value, self.R_value,
                                                                   saved_params)
        self.parameter_frame.pack(side='left', expand=True, fill=tk.X)

        # make a frame for the user to select the voltages to use and the ion concentrations to use
        _saved_sim_params = None
        # check if the user has inputted any saved parameters, or if there is are already parameters in the frame
        if _saved_settings:
            _saved_sim_params = _saved_settings['sim']
        elif self.simulation_frame:
            # take the old parameters of the frame and pass them to the new frame to make it with
            _saved_sim_params = self.simulation_frame.get_settings()
        if self.simulation_frame:  # if the there is already a simulation frame, destroy it so it can be replaced
            self.simulation_frame.destroy()

        # make the frame
        self.simulation_frame = sim_frame.SimulationWindow(self.top_frame, solute_text, _saved_sim_params)
        self.simulation_frame.pack(side='left')

        # add the run button the first time, afterwards change the config, not sure why this is needed but
        # num_binding_sites will not be correct if the else statement is not included
        if not self.run_button:
            self.run_button = tk.Button(self.run_button_frame, text="Run Simulation full",
                                        command=lambda: self.run_simulation(num_binding_sites))
            self.run_button.pack()
        else:
            self.run_button.config(command=lambda: self.run_simulation(num_binding_sites))  # hackish

        # change the size of the application depending on how many binding sites and solutes there are
        num_solutes = len(solute_text)
        _heigth = int(69 * num_solutes + 450)
        _width = int(88 * num_binding_sites + 980)
        self.geometry("%dx%d" % (_width, _heigth))

    def run_simulation(self, num_sites):
        voltages, concentrations, conc_labels = self.simulation_frame.get_run_simulation_settings()
        if not self.Q_value:
            self.Q_value = [1]  # hack
        barriers = self.parameter_frame.energy_barriers
        """
        results have the class Results found in numpy_helper_functions and has the attributes
        voltage, matrix_specs, ion_transport self.fitting,  current and steady_state
        """
        results_eig, results_svd, results_qr = eyring_rate_script.eyring_rate_algo(voltages, concentrations,
                                                                                   barriers,
                                                                                   num_barriers=num_sites,
                                                                                   Qs=self.Q_value, Rs=self.R_value)

        solutes = []
        for solute in results_eig[0].ion_transport:
            solutes.append(solute)

        result_toplevel.MultiPlotWindows(self, voltages, barriers, results_eig, solutes, conc_labels, "Eig results")

        result_toplevel.MultiPlotWindows(self, voltages, barriers, results_svd, solutes, conc_labels, "SVD results")

        result_toplevel.MultiPlotWindows(self, voltages, barriers, results_qr, solutes, conc_labels, "QR results")

    def get_settings(self, initial_settings):
        if initial_settings.Q_value_used:
            self.Q_value = [float(initial_settings.Q_entry.get())]  # hack
        else:
            self.Q_value = None
        if initial_settings.R_value_used:
            self.R_value = initial_settings.R_entry.get()
        else:
            self.R_value = None

        # make a list of the solute texts by getting them from the tkinter varstring
        solutes = []
        for _solute in initial_settings.solutes:
            solute = _solute.get()
            solutes.append(solute.replace(" ", "_"))  # replace spaces because they can not be used for variable names
        # solutes = [i.get().replace(" ", "_") for i in initial_settings.solutes]

        # make a list of the charges used, the fancy way
        charges = [i.get() for i in initial_settings.charges]

        return solutes, charges

    def make_menubar(self):
        menubar = tk.Menu(self)
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Open", command=self.open_file)
        filemenu.add_command(label="Save", command=self.save_settings)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self.quit)
        menubar.add_cascade(label="File", menu=filemenu)
        self.config(menu=menubar)

    def open_file(self):
        filename = tkFileDialog.askopenfilename()
        if not filename:  # user cancelled so just return
            return
        with open(filename, 'r') as model_file:
            model = model_file.read()
            num_sites = self.find_between(model, 'num sites: ', '\n')
            self.initial_frame.num_binding_sites.set(num_sites)
            print 'num sites: ', num_sites
            solutes = self.find_between(model, '\nsolutes: ', ' \n').split(' ')
            num_solutes = len(solutes)
            print 'num solutes: ', num_solutes, solutes
            # self.initial_frame.solutes = solutes
            self.initial_frame.remove_all_solute_entries()
            self.initial_frame.num_solutes_entries = num_solutes
            self.initial_frame.num_solutes.set(num_solutes)
            self.initial_frame.solute_num_event(num_solutes)


            charges_str = self.find_between(model, '\ncharges: ', ' \n').split(' ')
            charges = [int(x) for x in charges_str]
            for i, solute in enumerate(solutes):
                self.initial_frame.add_solute_entry(solute, charges[i])
            print 'charges: ', charges
            if 'no Q factor' in model:
                self.initial_frame.Q_value_used = False
                self.initial_frame.Q_on_off.set(0)
                self.initial_frame.Q_check_call()
            else:
                self.initial_frame.Q_value_used = True
                self.initial_frame.Q_on_off.set(1)
                self.initial_frame.Q_check_call()
                self.initial_frame.Q_value = float(self.find_between(model, 'Q factor: ', '\n'))
        if 'energy barriers:' in model:
            # the model has a parameter and simulation frame settings
            _parameter_settings_str = self.find_between(model, '\nenergy barriers: ', '\nvoltages:')
            print 'param setting str ', _parameter_settings_str
            _energy_barriers = dict()
            for _setting_str in _parameter_settings_str.split('\n'):
                _key, _values = _setting_str.split(':')
                _energy_barriers[_key] = [float(x) for x in _values.split(', ')]
            print 'energy build: ', _energy_barriers



            _simulation_settings_str = model.split('solute concentrations:\n')[1]
            print 'sim settings: ', _simulation_settings_str
            _sim_settings = dict()
            for _conc in _simulation_settings_str.split('\n'):
                if _conc:
                    print 'conc: ', _conc
                    _key, _value = _conc.split(': ')
                    _sim_settings[_key] = [_value[:-2], _value[-2:]]
            print 'sim settings: done', _sim_settings
            _settings = dict()
            _settings['param'] = _energy_barriers
            _settings['sim'] = _sim_settings
            self.process_parameters(self.initial_frame, _settings)

        # self.initial_frame.update()

    def find_between(self, target, _start, _end):
        """
        return the substring between the strings _start and _end
        """
        return target.split(_start)[1].split(_end)[0]


    def save_settings(self):
        """
        Save the settings of the model into a file so the user can load them later
        """
        solutes, charges = self.get_settings(self.initial_frame)
        if self.initial_frame.Q_value_used:
            _Qstr = 'Q factor: ' + self.initial_frame.Q_value.get()
        else:
            _Qstr = 'no Q factor'

        filename = tkFileDialog.asksaveasfilename(**FILEOPTIONS)
        if not filename:
            return
        with open(filename, 'w') as model_file:
            model_file.write('num sites: %d\n' % self.initial_frame.num_binding_sites.get())
            model_file.write('num solutes: %d\n' % self.initial_frame.num_solutes_entries)
            model_file.write('solutes: ')
            for solute in solutes:
                model_file.write('%s ' % solute)
            model_file.write('\ncharges: ')
            for charge in charges:
                model_file.write('%d ' % charge)
            model_file.write('\n%s\n' % _Qstr)
            if self.parameter_frame:
                model_file.write('energy barriers: ')
                for k, v in self.parameter_frame.energy_barriers.items():
                    v_str = ', '.join(map(str, v))
                    model_file.write('%s: %s\n' % (k, v_str))
                voltages, concentrations = self.simulation_frame.get_settings()
                voltage_str = ', '.join(map(str, voltages))
                model_file.write('voltages: %s\nsolute concentrations:\n' % voltage_str)
                for k, v in concentrations.items():
                    conc_str = str(v[0]) + v[1]
                    model_file.write('%s: %s\n' % (k, conc_str))


class InitialFrame(tk.Frame):
    """
    Frame for the user select parameters to run an Eyring rate model with by selecting the number
    of binding sites, the number of solutes, the names of the solutes, and if there is a Q and or an R factor
    public attributes::
    solutes: list of solutes the user wants to use
    charges: list of charges of the solutes, indexed to be the same order as solutes
    num_binding_sites: tk.IntVar :number of binding sites the user has selected
    num_solute_entries: int: number to keep track of how many entries the user wants
    Q_on_off: tk.IntVar: to check if the user wants to use a Q factor
    """
    def __init__(self, master=None):
        tk.Frame.__init__(self, master=master)
        self._solute_entry_boxes = []  # list of the tk.Entry for the user to entry solutes' names
        self._charge_spinboxes = []  # list of tk.Spinboxes for the user to select the solutes' charges
        self.solutes = []  # list to store the string of names of the solutes the user inputs
        self.charges = []  # list to store the valence charge of the solutes
        self.Q_value_used = False  # to tell if the user wants to use a Q value
        self.R_value_used = False  # to tell if the user wants to use a R value
        self.Q_entry = None  # initialize
        self.R_entry = None  # initialize
        self.mpmath_used = False

        # make region for the user to select the number of binding sites and bind the number to self.num_binding_sites
        self.num_binding_sites = tk.IntVar()
        tk.Label(master=self, text="Number of binding sites").pack(side='top', pady=5)
        tk.Spinbox(self, from_=1, to=6,
                   textvariable=self.num_binding_sites,
                   width=5).pack(side='top', pady=1)

        # make a place wher ethe user can choose the number of solutes to include in the model and bind the
        # number to solute_num_event to be called when the user clicks on the spinbox
        self.num_solutes_entries = 1
        self.num_solutes = tk.IntVar()
        tk.Label(master=self, text="Numebr of solutes").pack(side='top', pady=5)
        solute_spinbox = tk.Spinbox(self, from_=1, to=6,
                                    textvariable=self.num_solutes,
                                    command=lambda: self.solute_num_event(self.num_solutes.get()),
                                    width=5)
        solute_spinbox.pack(side='top', pady=1)

        solute_charges_frame = tk.Frame(self)
        solute_charges_frame.pack(side='top')
        self._solute_frame = tk.Frame(solute_charges_frame)
        self._solute_frame.pack(side='left')
        self._charges_frame = tk.Frame(solute_charges_frame)
        self._charges_frame.pack(side='left')

        self.init_solute_frame()
        self.init_charge_frame()

        # make options for the use to select if a Q factor should be used and what the value is
        self.Q_on_off = tk.IntVar()
        self.Q_value = tk.StringVar()
        self.R_on_off = tk.IntVar()
        self.R_value = tk.StringVar()
        self.make_Q_R_boxes()

        # allow the user to use the mpmath package
        self.mpmath_used_var = tk.IntVar()
        tk.Checkbutton(self, text="Use mpmath?",
                       variable=self.mpmath_used_var,
                       command=self.set_math_package).pack(side='bottom')

    def make_Q_R_boxes(self):
        QR_frame = tk.Frame(self)
        QR_frame.pack(side='top')

        Q_checkbox = tk.Checkbutton(QR_frame, text="Use Q factor?",
                                    variable=self.Q_on_off,
                                    command=self.Q_check_call)

        Q_checkbox.pack(side='top')
        self.Q_value.set(5.0)
        self.Q_entry = tk.Entry(QR_frame, state='disabled', width=10,
                                textvariable=self.Q_value, disabledbackground='grey')
        self.Q_entry.pack(side='top')

        R_checkbox = tk.Checkbutton(QR_frame, text="Use R factor?",
                                    variable=self.R_on_off,
                                    command=self.R_check_call)
        R_checkbox.pack(side='top')
        self.R_value.set(0.5)
        self.R_entry = tk.Entry(QR_frame, state='disabled', width=10,
                                textvariable=self.R_value, disabledbackground='grey')
        self.R_entry.pack(side='top')

    def Q_check_call(self):
        if self.Q_on_off.get():
            self.Q_entry.config(state='normal')
            self.Q_value_used = True
        else:
            self.Q_entry.config(state='disabled')
            self.Q_value_used = False

    def R_check_call(self):
        if self.R_on_off.get():
            self.R_entry.config(state='normal')
            self.R_value_used = True
        else:
            self.R_entry.config(state='disabled')
            self.R_value_used = False

    def solute_num_event(self, num):
        if num > self.num_solutes_entries:
            self.add_solute_entry()
        elif num < self.num_solutes_entries:
            self.remove_solute_entry()

    def init_charge_frame(self):
        tk.Label(self._charges_frame, text="Charge").pack(side='top')
        self.charges.append(tk.IntVar())
        self._charge_spinboxes.append(tk.Spinbox(self._charges_frame, from_=-max_charge, to=max_charge,
                                                 textvariable=self.charges[0], width=3))
        self._charge_spinboxes[0].pack(side='top')

    def init_solute_frame(self):
        """
        initialize the solute frame by making 1 solute entry
        :return:
        """
        tk.Label(self._solute_frame, text="Solutes").pack(side='top')
        self.solutes.append(tk.StringVar())
        self.solutes[0].set("solute 1")
        self._solute_entry_boxes.append(tk.Entry(self._solute_frame, textvariable=self.solutes[0]))
        self._solute_entry_boxes[0].pack(side='top')

    def add_solute_entry(self, _solute=None, _charge=None):
        """
        add a new entry box for the user to type in the name of a new solute and charge spinbox
        :return:
        """
        # make new solute box
        self.solutes.append(tk.StringVar())
        if not _solute:
            self.solutes[-1].set("solute "+str(self.num_solutes_entries+1))
        else:
            self.solutes[-1].set(_solute)
        self._solute_entry_boxes.append(tk.Entry(self._solute_frame, textvariable=self.solutes[-1]))
        self._solute_entry_boxes[-1].pack(side='top')

        # add new charges spin box
        self.charges.append(tk.IntVar())
        if _charge:
            self.charges[-1].set(_charge)
        self._charge_spinboxes.append(tk.Spinbox(self._charges_frame, from_=-max_charge, to=max_charge,
                                                 textvariable=self.charges[-1], width=3))
        self._charge_spinboxes[-1].pack(side='top')
        # update number of solute entries and update the frame
        self.num_solutes_entries += 1
        self.update()

    def remove_all_solute_entries(self):
        while self.solutes:
            self.remove_solute_entry()

    def remove_solute_entry(self):
        """
        When the user lowers the number of solutes, remove the last solute entered
        :return:
        """
        # remove solute entry and entry box
        self.solutes.pop()
        item = self._solute_entry_boxes.pop()
        item.destroy()
        # remove charge box
        self.charges.pop()
        item = self._charge_spinboxes.pop()
        item.destroy()
        self.num_solutes_entries -= 1
        self.pack()

    def set_math_package(self):
        self.mpmath_used = self.mpmath_used_var.get()


if __name__ == '__main__':
    app = EyringGUI()
    app.title("Eyring Rate modeler")
    app.geometry("200x450")
    app.mainloop()