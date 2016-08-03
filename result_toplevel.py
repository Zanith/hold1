import Tkinter as tk
from tkFileDialog import asksaveasfilename
import csv

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg


__author__ = 'Kyle Vitautas Lopin'

colors = ['k', 'r', 'b', 'g', 'm', 'c']


class MultiPlotWindows(tk.Toplevel):
    def __init__(self, master, voltage, energy_profile, results, solutes, conc_label, _title):
        """
        Make a tkinter toplevel that displays the 1) energy profiles, and 2) concentrations used for the simulation with
        3) the current per channel and 4) each solutes transport rates calculated in matplotlib subplots
        :param master:  root tk.Tk the toplevel is made from
        :param voltage: list of voltages to display
        :param energy_profile:  dict of energy profiles used in the calculation.  The keys are the solutes used plus a
        'distance'  entry, the values are a list of the energy potentials of the barriers and binding sites except for
        the 'distance' value which is a list of the electrical distances used
        :param results: claa Results found in numpy_helper_function
        :param solutes:  list of the strings of the solutes used
        :param conc_label: dict of concentrations used, the keys are the solutes + a 'e' or 'i' to represent if the
        concentration is the extra or intracellular concentration, the values are strings of the user
        entered concentrations
        :return:
        """
        tk.Toplevel.__init__(self, master=master)
        self.title(_title)
        # bind variables to self so they don't have to be passed around as much
        self.voltages = voltage
        self.energy_profile = energy_profile
        self.results = results
        self.conc_label = conc_label
        self.conc_label = conc_label
        self.solutes = solutes
        # make variables to hold measurements
        # dictionary to hold transport rates, keys are solutes and values are the a list
        # transport rates over each barrier
        self.transport = dict()
        transport_errors = dict()  # keep any differences between different transport rates over different barriers
        self.current = []  # calculated current
        ss = []  # steady state solution of states
        for solute in solutes:
            self.transport[solute] = []
            transport_errors[solute] = []
        for result in results:
            ss.append(result.steady_state)
            self.current.append(result.current[1])  # use the transport over the second barrier
            for solute in solutes:
                # get the ions transported over the 2nd barrier
                self.transport[solute].append(result.ion_transport[solute][1])
                # get errors calculated for the ions transported over the different barriers
                transport_errors[solute].append(max(result.ion_transport[solute])-min(result.ion_transport[solute]))
        if len(solutes) < 2:
            self._size = (2, 2)
        elif len(solutes) > 3:
            self._size = (3, 3)
        else:
            self._size = (2, 3)

        # set the geometry based on how many plots are to be made
        geometry_size = "%dx%d" % (self._size[1]*300, self._size[0]*300)
        self.geometry(geometry_size)
        # make matplotlib figure and get the tkinter canvas backend for the figure
        self.fig = plt.figure()
        self.fig.set_facecolor('white')  # TODO: think this looks good on mac but not windows
        canvas = FigureCanvasTkAgg(self.fig, self)
        canvas._tkcanvas.config(highlightthickness=0)  # get rid of outline

        # make plot of energy profiles
        energy_profile_plot = self.fig.add_subplot(self._size[0], self._size[1], 1)
        energy_plot(energy_profile_plot, energy_profile)

        # make graph to show the concentration on either side of the membrane
        self.display_concentrations(conc_label)

        # make current graph
        self.plot_routine(self.current, "current (pA)", "Total Current per Channel", 3, sci_format=True)

        # make transport plots for each solute transported
        _index = 4
        for solute in solutes:
            self.plot_routine(self.transport[solute], "ions trasnported per second",
                              solute+" transport", _index, sci_format=True)
            _index += 1

        self.fig.tight_layout()
        canvas.draw()
        canvas.get_tk_widget().pack(expand=True, fill=tk.BOTH)

        # make toolbar
        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()

        # make button to save data
        tk.Button(self, text='Save data', command=self.save_data).pack(side='bottom')

        # make button to make a custom save data function
        tk.Button(self, text='Save custom_data', command=self.save_custom_data).pack(side='bottom')

        # make button to make a custom save data function
        tk.Button(self, text='Save custom_data2', command=self.save_custom_data2).pack(side='bottom')

    def display_concentrations(self, conc_label):
        """
        Make a plot to display the concentrations of the solutes used for the simulation
        Spacing is mostly done by magic numbers
        :param conc_label:  dictionary of the concentrations used
        :return:
        """
        concentration_plot = self.fig.add_subplot(self._size[0], self._size[1], 2)
        concentration_plot.axis('off')
        concentration_plot.plot([15, 15], [0, 15], 'k', linewidth=3.0)
        concentration_plot.set_xlim([0, 30])
        y_spacing = 15
        concentration_plot.text(16.5, y_spacing, r'Intracellular')
        concentration_plot.text(0.5, y_spacing, r'Extracellular')
        for solute in self.solutes:
            y_spacing -= 3
            concentration_plot.text(16, y_spacing, '%s:\n    %s' % (solute, conc_label[solute+'i']))
            concentration_plot.text(0, y_spacing, '%s:\n    %s' % (solute, conc_label[solute+'e']))

    def plot_routine(self, y, _ylabel, _title, plot_index, sci_format=False):
        """
        routine to plot currents or ion transport rates
        :param y: list of y axis data
        :param _ylabel: string to label the y axis
        :param _title: title to give the plot
        :param plot_index: what position to put the plot in, in the overall figure
        :param sci_format: turn on or off scientific plotting which scales the y axis
        :return:
        """
        _plot = self.fig.add_subplot(self._size[0], self._size[1], plot_index)
        _plot.plot(self.voltages, y)
        _plot.set_xlabel("voltage (mV)")
        _plot.set_ylabel(_ylabel)
        _plot.set_title(_title, y=1.08)  # move the title up to get it out of the way of the science formatting thing
        if sci_format:
            _plot.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    def save_data(self):
        """
        Save the data with the columns in the format:
        voltage  current  solute_transport1  solute_transport2
        to a comma separated values file
        :return:
        """
        # get filename from user and open file
        filename = asksaveasfilename(defaultextension=".csv")
        with open(filename, 'wb') as csvfile:
            writer = csv.writer(csvfile, dialect='excel')  # write file as excel dialect
            # make a header with simulation details as the first line
            header_row = [str(self.energy_profile), str(self.conc_label)]
            writer.writerow(header_row)
            # make a row to label what each column will be
            _row = ['voltage', 'current']
            for solute in self.solutes:
                _row.append(solute)
            writer.writerow(_row)
            # make columns with the voltage, current and transport rates
            for i, voltage in enumerate(self.voltages):
                _row = [voltage, self.current[i]]
                for solute in self.solutes:
                    _row.append(self.transport[solute][i])
                writer.writerow(_row)  #write the row
        # file automatically closes with with statement

    def save_custom_data(self):
        """
        Make a list to test how well the programming is calculating the currents and tranport rate
        columns are: voltage:largest element:smallest_element:element difference:element ratio:
        condition number:transport difference solute 1:transport difference solute 2: ... :
        difference in current :
        smallest eig : second smallest eig : absolute errors in ss : sum squared errors in ss
        """
        voltages = self.voltages
        # make predictor variables
        largest_element = []
        smallest_element = []
        element_difference = []
        ratio_of_elements = []
        condition_number = []

        # make variables of how well the calculation was performed
        smallest_eig = []
        second_smallest_eig = []
        sae = []
        sse = []

        transport_errors = dict()
        current_errors = []

        for solute in self.solutes:

            # self.transport[solute] = []
            transport_errors[solute] = []

        for result in self.results:
            largest_element.append(result.matrix_spec.largest_element)
            smallest_element.append(result.matrix_spec.smallest_element)
            element_difference.append(largest_element[-1]-smallest_element[-1])
            ratio_of_elements.append(largest_element[-1]/smallest_element[-1])
            condition_number.append(result.matrix_spec.condition_number)

            smallest_eig.append(result.fitting.smallest_eig)
            second_smallest_eig.append(result.fitting.second_smallest_eig)
            sae.append(result.fitting.sae_residues)
            sse.append(result.fitting.sse_residues)

            current_errors.append(max(result.current)-min(result.current))

            for solute in self.solutes:
                # get errors calculated for the ions transported over the different barriers
                print 'voltage: ', result.voltage
                print 'solute: ', solute
                print 'transport: ', result.ion_transport[solute]
                transport_errors[solute].append(max(result.ion_transport[solute])-min(result.ion_transport[solute]))

        # get filename from user and open file
        filename = asksaveasfilename(defaultextension=".csv")
        if not filename:
            return  # exit if save was cancelled
        with open(filename, 'wb') as csvfile:
            writer = csv.writer(csvfile, dialect='excel')  # write file as excel dialect
            # make a header with simulation details as the first line
            header_row = [str(self.energy_profile), str(self.conc_label)]
            writer.writerow(header_row)
            _row = ['voltage', 'largest element', 'smallest element', 'element difference', 'element ratio difference',
                    'condition number']
            for solute in self.solutes:
                _row.append('transport error '+solute)
            _row.extend(['current errors', 'smallest eigenvalue', 'second smallest eigenvalue',
                         'sum absolute errors', 'sum squared errors'])
            writer.writerow(_row)

            for i, voltage in enumerate(voltages):
                _row = [voltage, largest_element[i], smallest_element[i], element_difference[i], ratio_of_elements[i],
                        condition_number[i]]
                for solute in self.solutes:
                    _row.append(transport_errors[solute][i])
                _row.extend([current_errors[i], smallest_eig[i], second_smallest_eig[i], sae[i], sse[i]])
                writer.writerow(_row)

    def save_custom_data2(self):
        """
        Class to save the current and solute transport over each barrier
        """
        # ask the use what to call the data file
        filename = asksaveasfilename(defaultextension=".csv")
        if not filename:
            return  # exit if save was cancelled
        with open(filename, 'wb') as csvfile:
            writer = csv.writer(csvfile, dialect='excel')  # write file as excel dialect

            # make a header with simulation details as the first line
            header_row = [str(self.energy_profile), str(self.conc_label)]
            writer.writerow(header_row)
            # make column headers with details of the columns
            _row = ['voltage']
            _i = len(self.results[0].current)
            for i in range(_i):
                _row.append('current, barrier %d' % i)
            for solute in self.solutes:
                for i in range(_i):
                    _row.append('%s transport over barrier %d' % (solute, i))
            writer.writerow(_row)

            # write the results for each variable and write the row
            for result in self.results:
                _row = [result.voltage]
                _row.extend(result.current)
                for solute in self.solutes:
                    _row.extend(result.ion_transport[solute])
                writer.writerow(_row)
        # file automatically closes with with statement


def energy_plot(plot_axis, profiles):
    """
    Make a graph of energy profiles
    :param plot_axis: the pyplot figure to graph in
    :param profiles: dictionary of energy profiles to graph, must have a key 'distance' with the electrical
    distances to graph on the x-axis and the other values are lists with the energy barrier values
    :return:
    """

    x = [-0.2, 0] + profiles['distance'] + [1, 1.2]
    i = 0

    for solute in profiles:
        if solute != 'distance':
            y = [0, 0] + profiles[solute] + [0, 0]
            plot_axis.plot(x, y, colors[i], label=solute)
            i += 1
    # TODO: put concentrations in energy barrier graph?
    # ymin, ymax = plot_axis.get_ylim()
    # y_spread = float(ymax) - ymin
    # top_y_place = ymin + y_spread/2.
    # y_difference = y_spread / 14.
    # plot_axis.text(-0.15, top_y_place, r'Extracellular', fontsize=10)
    # plot_axis.text(0.85, top_y_place, r'Intracellular', fontsize=10)
    plot_axis.set_xlabel("electrical distance")
    plot_axis.set_ylabel("energy (kT)")
    plot_axis.legend(prop={'size': 10}, loc='upper right')
    plot_axis.set_title("Energy Profiles")


if __name__ == '__main__':
    # to test the toplevel pass in these test values to run it
    app = tk.Tk()
    voltage = range(-150, 110, 10)
    current = range(-15, 11, 1)
    transported = {'solute_1': [-378450.7354495374, -338463.9353080591, -301864.595992801, -268286.4190144333, -237393.3426917416, -208876.17873336677, -182449.51777746368, -157848.87291887807, -134828.03263530295, -113156.59661960867, -92617.66985612901, -73005.69186244335, -54124.37937097749, -35784.761860088474, -17803.290273569044, -1.4551915228366852e-11, 17803.290273569044, 35784.761860088474, 54124.379370977476, 73005.69186244336, 92617.66985612901, 113156.59661960864, 134828.03263530298, 157848.87291887804, 182449.5177774637, 208876.1787333668]}
    profile = {'distance': [0.25, 0.5, 0.75], 'solute_1': [8.0, -10.0, 8.0]}
    conc_label = {'solute_1i': '120.0 mM', 'solute_1e': '1.0 M', 'solute 2i': '120.0 mM', 'solute 2e': '1.0 M'}
    concs = {'solute_1i': 0.001, 'solute_1e': 0.001}
    size = (2, 2)
    MultiPlotWindows(app, voltage, current, transported, profile, conc_label)
    app.mainloop()
