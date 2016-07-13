__author__ = 'Kyle Vitautas Lopin'

import Tkinter as tk
import gui_parameter_frame as param_frame
import gui_simulation_frame as sim_frame
import template_maker as script_maker
import eyring_rate_script
import result_toplevel

max_charge = 5


class EyringGUI(tk.Tk):
    def __init__(self, parent=None):
        tk.Tk.__init__(self, parent)
        self.top_frame = tk.Frame(self)
        self.top_frame.pack(side='top')
        self.bottom_button_frame = tk.Frame(self)
        self.bottom_button_frame.pack(side='top')
        self.run_button_frame = tk.Frame(self.bottom_button_frame)
        self.run_button_frame.pack(side='top')

        initial_frame = InitialFrame(self.top_frame)
        self.hold_parameter_frame = tk.Frame(self.top_frame)
        self.parameter_frame = None
        self.simulation_frame = None
        self.run_button = None
        initial_frame.pack(side='left', anchor='n')
        make_params_button = tk.Button(initial_frame, text="Set model parameters",
                                       command=lambda: self.process_parameters(initial_frame))
        make_params_button.pack(side='bottom')

    def process_parameters(self, initial_settings):
        """
        Process information to properly set the initial settings
        """
        if initial_settings.Q_value_used:
            Q_value = initial_settings.Q_entry.get()
        else:
            Q_value = None
        if initial_settings.R_value_used:
            R_value = initial_settings.R_entry.get()
            R_str = "R"
        else:
            R_value = None
            R_str = ""

        # make a list of the solute texts by getting them from the tkinter varstring
        solutes = []
        for _solute in initial_settings.solutes:
            solute = _solute.get()
            solutes.append(solute.replace(" ", "_"))
        # solutes = [i.get().replace(" ", "_") for i in initial_settings.solutes]

        # make a list of the charges used
        charges = [i.get() for i in initial_settings.charges]
        num_binding_sites = initial_settings.num_binding_sites.get()

        QR_str = "single Q" + R_str
        print "make script with variables:"
        print "num binding sites: ", num_binding_sites
        print 'solutes: ', solutes
        print "charges: ", charges
        print "QR string: ", QR_str

        # make the script for the Eyring rate model

        _eyring_rate_script = script_maker.make_template("numpy", num_binding_sites,
                                                         solutes, charges, QR_str)
        with open("eyring_rate_script.py", "w") as _file:
            _file.write(_eyring_rate_script)
        _file.close()
        reload(eyring_rate_script)

        self.set_parameters(solutes,
                            initial_settings.charges,
                            initial_settings.num_binding_sites.get(),
                            Q_value, R_value)

    def set_parameters(self, solute_text, charges, num_binding_sites, Q_value, R_value):
        saved_params = None
        if self.parameter_frame:
            # take the old parameters of the frame and pass them to the new frame to make it with
            saved_params = self.parameter_frame.energy_barriers
            self.parameter_frame.destroy()
        self.hold_parameter_frame.pack(side='left')

        # make a frame for the user to make the energy profiles to use
        self.parameter_frame = param_frame.ParameterSelectionFrame(self.hold_parameter_frame,
                                                                   solute_text, charges,
                                                                   num_binding_sites,
                                                                   Q_value, R_value,
                                                                   saved_params)
        self.parameter_frame.pack(side='left', expand=True, fill=tk.X)

        # make a frame for the user to select the voltages to use and the ion concentrations to use
        _saved_sim_params = None
        if self.simulation_frame:
            _saved_sim_params = self.simulation_frame.get_settings()
            self.simulation_frame.destroy()

        self.simulation_frame = sim_frame.SimulationWindow(self.top_frame, solute_text, _saved_sim_params)
        self.simulation_frame.pack(side='left')

        if not self.run_button:
            self.run_button = tk.Button(self.run_button_frame, text="Run Simulation full",
                                        command=lambda: self.run_simulation(num_binding_sites, Q_value, R_value))
            self.run_button.pack()

        num_solutes = len(solute_text)
        _heigth = int(69 * num_solutes + 450)
        _width = int(88 * num_binding_sites + 980)
        self.geometry("%dx%d" % (_width, _heigth))

    def run_simulation(self, num_sites, Q, R):
        print "run simulation"
        voltages, concentrations, conc_labels = self.simulation_frame.get_run_simulation_settings()
        print 'voltages: ', voltages
        print "concentrations: ", concentrations
        print "energy barriers: ", self.parameter_frame.energy_barriers
        print "num sites: ", num_sites

        if not Q:
            Q = [1]  # hack
        print 'Qs: ', Q
        print 'Rs: ', R
        barriers = self.parameter_frame.energy_barriers
        # results have have the class
        results_eig, results_svd = eyring_rate_script.eyring_rate_algo(voltages, concentrations,
                                                                       barriers,
                                                                       num_barriers=num_sites,
                                                                       Qs=Q, Rs=R)
        ions_transported = dict()
        transport_errors = dict()
        current = []
        solutes = []
        int_ss = []
        for solute in results_eig[0].ion_transport:
            solutes.append(solute)
            ions_transported[solute] = []
            transport_errors[solute] = []
        for result in results_eig:
            int_ss.append(result.steady_state)
            current.append(result.current[0])
            for solute in result.ion_transport:
                # get the ions transported over the 2nd barrier
                ions_transported[solute].append(result.ion_transport[solute][1])
                # get errors calculated for the ions transported over the different barriers
                transport_errors[solute].append(max(result.ion_transport[solute])-min(result.ion_transport[solute]))
        print 'ions transported'
        print ions_transported
        print 'int ss: ', int_ss
        print transport_errors
        print current
        print solutes, len(solutes)
        if len(solutes) < 2:
            size = (2, 2)
        elif len(solutes) > 3:
            size = (3, 3)
        else:
            size = (2, 3)
        print 'size: ', size
        result_toplevel.MultiPlotWindows(self, voltages, current, ions_transported,
                                         barriers, conc_labels, solutes, size)


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

        # make region for the user to select the number of binding sites and bind the number to self.num_binding_sites
        self.num_binding_sites = tk.IntVar()
        tk.Label(master=self, text="Number of binding sites").pack(side='top', pady=5)
        tk.Spinbox(self, from_=1, to=6,
                   textvariable=self.num_binding_sites,
                   width=5).pack(side='top', pady=1)

        # make a place wher ethe user can choose the number of solutes to include in the model and bind the
        # number to solute_num_event to be called when the user clicks on the spinbox
        self.num_solutes_entries = 1
        num_solutes = tk.IntVar()
        tk.Label(master=self, text="Numebr of solutes").pack(side='top', pady=5)
        solute_spinbox = tk.Spinbox(self, from_=1, to=6,
                                    textvariable=num_solutes,
                                    command=lambda: self.solute_num_event(num_solutes.get()),
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

    def add_solute_entry(self):
        """
        add a new entry box for the user to type in the name of a new solute and charge spinbox
        :return:
        """
        # make new solute box
        self.solutes.append(tk.StringVar())
        self.solutes[-1].set("solute "+str(self.num_solutes_entries+1))
        self._solute_entry_boxes.append(tk.Entry(self._solute_frame, textvariable=self.solutes[-1]))
        self._solute_entry_boxes[-1].pack(side='top')

        # add new charges spin box
        self.charges.append(tk.IntVar())
        self._charge_spinboxes.append(tk.Spinbox(self._charges_frame, from_=-max_charge, to=max_charge,
                                                 textvariable=self.charges[-1], width=3))
        self._charge_spinboxes[-1].pack(side='top')
        # update number of solute entries and update the frame
        self.num_solutes_entries += 1
        self.update()

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


if __name__ == '__main__':
    app = EyringGUI()
    app.title("Eyring Rate modeler")
    app.geometry("200x450")
    app.mainloop()