import Tkinter as tk
import gui_simulation_frame as sim_frame
import pyplot_to_tkinter as pytink
import template_maker as tm

__author__ = 'Kyle Vitautas Lopin'


voltages = range(-150, 110, 10)

class ParameterSelectionFrame(tk.Frame):
    """
    Class to make and display a energy barrier and binding site profile
    Public attributes:
    ::self.energy_barriers: dictionary with the keys of "energy
    """
    def __init__(self, _master, solutes, _charges, num_binding_sites, Q_value, R_value, _settings=None):
        tk.Frame.__init__(self, master=_master, bd=5, relief='raised')
        self.done_initialize = False
        charges = []
        if not _settings:
            _settings = dict()
        # process the _settings variable
        if ('distance' not in _settings) or (len(_settings['distance']) != 2*num_binding_sites+1):
            _settings['distance'] = None
        self.parameters = _settings
        self.solutes = solutes
        for _charge in _charges:
            charges.append(_charge.get())
        self.distances = []
        self.energies_var = dict()
        self.energy_barriers = dict()
        self.parameter_plot = self.pyplot_embed(_master)
        self.parameter_plot.pack(side='top')

        self.make_Q_value_display(Q_value, R_value)

        self.make_parameters_frame(num_binding_sites, solutes)
        buttons_frame = tk.Frame(self)
        buttons_frame.pack(side='bottom')
        self.done_initialize = True
        self.update_energies(solutes)

    def make_Q_value_display(self, Q_values, R_values):
        if Q_values:
            Q_str = "Q value = " + str(Q_values)
        else:
            Q_str = "no Q value"
        if R_values:
            R_str = "R value = " + str(R_values)
        else:
            R_str = "no R value"
        tk.Label(self, text=Q_str+"     "+R_str, bd=2, relief='groove').pack(side='top', expand=True, fill=tk.X)

    def update_energies(self, solutes_list):
        """
        Check if the user has input sensible values and update the energy profile plot
        :param solutes_list:
        :return:
        """
        if self.done_initialize:
            # convert the distances into a list of floats
            distances_float = []
            for distance in self.distances:
                distances_float.append(distance.get())
                if len(distances_float) > 1:  # check if distances are in order
                    if distances_float[-1] < distances_float[-2]:
                        self.distances[-2].set(distances_float[-1])
            self.energy_barriers['distance'] = distances_float
            # convert the energy barriers into a list of floats
            for solute in self.energies_var:
                self.energy_barriers[solute] = []
                for energy in self.energies_var[solute]:
                    try:
                        self.energy_barriers[solute].append(energy.get())
                    except:
                        pass  # ignore if user put in something weird and keep old value

            # check that the inputted values make sense

            # self.run_button.config(state="active")
            self.parameter_plot.draw_barriers(self.energy_barriers, solutes_list)

    def run_simulation_settings(self, num_sites, solutes, charges, Q, R):
        print 'run simulation'
        print 'num binding sites: ', num_sites
        print 'solutes: ', solutes
        print "charges: ", charges
        print "making script"
        if Q and R:
            QR_str = "single QR"
        elif Q:
            QR_str = "single Q"
        else:
            QR_str = "single R"
        eyring_rate_script = tm.make_template('numpy', num_sites, solutes, charges, QR_str)
        with open("gui_eyring_rate_script.py", "w") as _file:
            _file.write(eyring_rate_script)
        _file.close()

        # import the eyring rate model file that was just made

        #
        #
        # results_eig, results_svd = script.eyring_rate_algo(voltages,
        #                                                    ion_conc,
        #                                                    energy_barriers,
        #                                                    num_binding_sites, Qs, Rs)
        return self.num_binding_sites

    def parameter_changed(self, *args):
        """
        disable the run button if the user has changed a parameter so that the user has to
        update the energy profile before running the simulation
        :param args: needed so that the method can be bound to an event trace
        :return:
        """
        # check to make sure the full app has been made first, not sure why this is needed
        self.update_energies(self.solutes)

    def pyplot_embed(self, _master):
        _embedded_pyplot = pytink.PyplotEmbed(self)
        return _embedded_pyplot

    def make_parameters_frame(self, num_sites, solutes):
        frame = tk.Frame(self)  # make a frame to put all the parameter selection widgets in
        barrier_distance_frame = tk.Frame(frame)  # frame to put all the electrical distances widgets in
        # make all the electrical distance widgets in a seperate method and pack the frame
        self.make_electrical_distances(barrier_distance_frame, num_sites)
        barrier_distance_frame.pack(side='top', fill=tk.X, expand=True, padx=15)

        # make a frame to select all the barrier and binding site energies, call a method to make it, then pack it
        barrier_frame = tk.Frame(frame)
        self.make_barrier_params(barrier_frame, num_sites, solutes)
        barrier_frame.pack(side='top', fill=tk.X, expand=True, padx=15)

        frame.pack(side='top', fill=tk.X, expand=True)

    def make_electrical_distances(self, _frame, num_sites):
        """
        Make a a frame with a spinbox for each electrical distance
        :param _frame: parent frame to put everything in
        :return:
        """
        num_barriers = num_sites + 1
        tk.Label(master=_frame, text="electrical distances of barriers:").pack(side='top', anchor='w')
        distance_padding = 2*num_sites + 2
        barrier_frame = tk.Frame(master=_frame)
        for barrier_index in range(num_barriers):
            tk.Label(master=barrier_frame, text='d'+str(2*barrier_index+1)).pack(side='left',
                                                                                 fill=tk.X, expand=True)
            distance_instance = tk.DoubleVar()

            if self.parameters['distance']:
                distance_instance.set(self.parameters['distance'][2*barrier_index])
            else:
                distance_instance.set((2.*barrier_index+1)/distance_padding)
            # if the user writes a number into the spinbox, disable run button till user updates the profile
            distance_instance.trace("w", self.parameter_changed)
            spin_box_instance = tk.Spinbox(barrier_frame, from_=0, to=1, increment=0.01,
                                           width=6, textvariable=distance_instance)

            spin_box_instance.pack(side='left')

            self.distances.append(distance_instance)
            self.distances.append(-1)  # place holder to fill in with binding site distances later

        self.distances.pop()  # delete last instance because it is not used

        tk.Label(master=barrier_frame, text=' ').pack(side='left', fill=tk.X, expand=True)
        barrier_frame.pack(side='top', fill=tk.X, expand=True, padx=15)

        tk.Label(master=_frame, text="electrical distances of binding sites:").pack(side='top', anchor='w')
        site_distance_frame = tk.Frame(_frame)
        for site_index in range(num_sites):
            tk.Label(master=site_distance_frame, text='d'+str(2*site_index+2)).pack(side='left', fill=tk.X, expand=True)
            distance_instance = tk.DoubleVar()
            if self.parameters["distance"]:
                distance_instance.set(self.parameters["distance"][2*site_index+1])
            else:
                distance_instance.set((2.*site_index+2)/distance_padding)
            # if the user writes a number into the spinbox, disable run button till user updates the profile
            distance_instance.trace("w", self.parameter_changed)
            spin_box_instance = tk.Spinbox(site_distance_frame, from_=0, to=1, increment=0.01,
                                           width=6, textvariable=distance_instance)

            spin_box_instance.pack(side='left')
            # self.electrical_distance_spinboxes[2*site_index+1] = spin_box_instance
            self.distances[2*site_index+1] = distance_instance
        tk.Label(master=site_distance_frame, text=' ').pack(side='left', fill=tk.X, expand=True)
        site_distance_frame.pack(side='top', fill=tk.X, expand=True, padx=15)

    def make_barrier_params(self, _frame, num_sites, _solutes):
        num_barriers = num_sites + 1
        for solute in _solutes:
            self.energies_var[solute] = []
            solute_frame = tk.Frame(_frame, bd=2, relief='ridge')
            tk.Label(master=solute_frame, text='Barrier energies for '+solute).pack(side='top', anchor='w')
            barrier_frame = tk.Frame(solute_frame)
            for barrier_index in range(num_barriers):
                tk.Label(master=barrier_frame, text='barrier '+str(2*barrier_index+1)).pack(side='left',
                                                                                            fill=tk.X, expand=True)
                self.energies_var[solute].append(tk.DoubleVar())
                if (solute in self.parameters) and len(self.parameters[solute]) >= 2*barrier_index:
                    self.energies_var[solute][-1].set(self.parameters[solute][2*barrier_index])
                else:
                    self.energies_var[solute][-1].set(8)

                spin_box_instance = tk.Spinbox(barrier_frame, from_=-20, to=20, increment=0.1,
                                               width=6, textvariable=self.energies_var[solute][-1])

                spin_box_instance.pack(side='left')
                self.energies_var[solute][-1].trace("w", self.parameter_changed)
                self.energies_var[solute].append(0)  # place holder to fill in with binding energy
            self.energies_var[solute].pop()  # remove last item as it wont be filled
            barrier_frame.pack(side='top', fill=tk.X, expand=True)

            # make frame tp put in binding site energies
            site_frame = tk.Frame(solute_frame)
            for site_index in range(num_sites):
                tk.Label(master=site_frame, text="site"+str(2*site_index+2)).pack(side='left', fill=tk.X, expand=True)

                self.energies_var[solute][2*site_index+1] = tk.DoubleVar()
                if (solute in self.parameters) and len(self.parameters[solute]) >= 2*barrier_index+1:
                    self.energies_var[solute][2*site_index+1].set(self.parameters[solute][2*site_index+1])
                else:
                    self.energies_var[solute][2*site_index+1].set(-10)

                spin_box_instance = tk.Spinbox(site_frame, from_=-20, to=5, increment=0.1,
                                               width=6, textvariable=self.energies_var[solute][2*site_index+1])

                spin_box_instance.pack(side='left')
                self.energies_var[solute][2*site_index+1].trace("w", self.parameter_changed)

            # add spacer to look more centered
            tk.Label(site_frame, text="  ").pack(side='left', fill=tk.X, expand=True)
            site_frame.pack(side='left', fill=tk.X, expand=True)

            solute_frame.pack(side='top', fill=tk.X, expand=True)
