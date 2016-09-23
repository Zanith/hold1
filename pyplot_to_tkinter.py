# Copyright (c) 2015-2016 Kyle Lopin (Naresuan University) <kylel@nu.ac.th>
# Licensed under the GPL
""" Class to make frame with a pyplot embedded in it, and will draw energy barriers
in the pyplot when draw_barriers is called
"""
# standard libraries
import Tkinter as tk

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

__author__ = 'Kyle Vitautas Lopin'


COLORS = ['k', 'r', 'b', 'g', 'm', 'c']


class PyplotEmbed(tk.Frame):
    """
    Class with a pyplot embedded in a tkinter frame
    """
    def __init__(self, _master, _size=(6.2, 2.5)):

        tk.Frame.__init__(self, master=_master)

        self.extra_text = None  # initialize the text options variables
        self.intra_text = None

        self.figure_bed = plt.figure(figsize=_size)
        self.axis = self.figure_bed.add_subplot(111)

        box = self.axis.get_position()
        self.axis.set_position([box.x0, box.y0, box.width*0.8, box.height])

        # self.axis.set_axis_bgcolor('red')
        self.figure_bed.set_facecolor('white')
        self.canvas = FigureCanvasTkAgg(self.figure_bed, master=self)
        self.canvas._tkcanvas.config(highlightthickness=0)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side='top')

    def draw_barriers(self, energy_barriers, solutes):
        """
        Draw the energy profile
        :param energy_barriers: dict with keys of 'distance' and solutes the energy
        barriers are for
        :param solutes: list of solutes
        :return: display the pyplot
        """
        while self.axis.lines:  # if lines are already plotted, go through and delete them all
            self.axis.lines.pop()
        x = [-0.2, 0] + energy_barriers['distance'] + [1, 1.2]  # pad the ends of the plot
        for i, solute in enumerate(solutes):
            y = [0, 0] + energy_barriers[solute] + [0, 0]  # pad the ends of the plot
            if len(x) == len(y):  # incase the user inputted something weird just ignore it
                self.axis.plot(x, y, COLORS[i], label=solute)  # draw the plot
        # move the legend off to the side
        self.axis.legend(prop={'size': 11}, loc='center left', bbox_to_anchor=(1, 0.5))
        ymin, ymax = self.axis.get_ylim()
        y_spacing = ymin + (ymax-ymin)/7.
        if self.extra_text:  # remove previous labels in case they have moved
            self.extra_text.remove()
            self.intra_text.remove()

        self.extra_text = self.axis.text(-0.15, y_spacing, r'Extracellular',
                                         color='#323280', fontsize=12)
        self.intra_text = self.axis.text(0.83, y_spacing, r'Intracellular',
                                         color='#323280', fontsize=12)
        # rescale the axis
        self.axis.relim()
        self.axis.autoscale()
        self.canvas.show()  # redraw the canvas to update it
