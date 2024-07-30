#Code written by Jacob St. Martin

import tkinter

import numpy as np
from matplotlib.patches import Rectangle, Circle
from tkinter.ttk import Notebook
from tkinter import filedialog as fd
from timeit import default_timer as timer
from datetime import timedelta
from PIL import Image
import scipy
import PIL
from tkinter import *
import matplotlib
from PIL import ImageOps
import csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,NavigationToolbar2Tk)

matplotlib.use('TkAgg')


class CustomToolbar(NavigationToolbar2Tk):
    '''
    This class initiates the toolbar on the graph. I have added the save figure function to allow for saving.
    '''
    def save_Figure(self):
        filename = fd.asksaveasfilename(initialdir="/", title="Select file", filetypes=(
            ('png files', '*.png'), ("jpeg files", "*.jpg"), ("all files", "*.*")))
        plt.savefig(filename)

    def __init__(self, canvas_, parent_):
        self.toolitems = (
            ('Home', 'Reset original view', 'home', 'home'),
            ('Back', 'Back to previous view', 'back', 'back'),
            ('Forward', 'Forward to next view', 'forward', 'forward'),
            ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
            ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
            ('save_fig', 'Save Figure', 'Filesave', 'save_Figure')
        )
        NavigationToolbar2Tk.__init__(self, canvas_, parent_)


class interface(Frame):

    def __init__(self):
        super().__init__()

        self.initUI()     # Calls the initUI function which sets up the program
        self.mainloop()   # Starts tkinter

    def initUI(self):
        '''
        This initiates the GUI and sets default values. The GUI was built using tkinter libraries. Any questions about the
        use of tkinter here should be well addressed in the documentation for the library. In order to be able to interact
        with the figure in the ways that have been implemented here, matplotlib Figure must be used instead of pyplot.
        Although the Figure class is less straightforward, it is much more customizable and the library FigureCanvasTkAgg
        allows for the integration of Figures and tkinter.
        :return:
        '''

        #Initializes variables/class attributes that need to have default values
        self.coordinate = (0, 0)   # self.coordinate represents the current coordinate selected on the image, sets  (0,0) as the initial value
        self.voltage = 0           # self.voltage represents the current voltage selected on the image, sets 0 as the initial value
        self.rectangle = False     # Used to track whether there is an area selection, starts as false
        self.circ1 = None          # initializes the selected point indicator in subplot 1
        self.circ3 = None
        self.colorbar_init = False # To see if the colorbar has been initialized
        self.initialized = False   # this is a way to track whether an STS file has been loaded yet
        self.unclicked = True      # A way to check if there has been a click event on the figure yet
        self.rect1 = None          # initializes the area select indicator for subplot 1
        self.color_min = 0         # initializes a min value for the colorbar
        self.color_max = 1         # initializes a max value for the colorbar
        self.coords_to_lines = {}  # Initializes a dictionary that maps coordinates to the IV and didv curves for that point
        self.v_to_xyi = {}         # Initializes the dictionary that maps a voltage to the x,y, and i values for that voltage slice
        self.fig = plt.figure()    # initializes the matplotlib figure

        # adds the 4 subplots to the figure
        self.ax1 = self.fig.add_subplot(2, 2, 1)
        self.ax2 = self.fig.add_subplot(2, 2, 2)
        self.ax3 = self.fig.add_subplot(2, 2, 3)
        self.ax4 = self.fig.add_subplot(2, 2, 4)

        # formatting for the graphs
        self.ax2.grid()
        self.ax4.grid()
        self.ax1.title.set_text('Topography Image')
        self.ax2.title.set_text('IV Plot at (x,y) = ' + str(self.coordinate))
        self.ax2.set_ylabel('Current (nA)')
        self.ax3.title.set_text('Normalized DIdV Image at V = ' + str(self.voltage))
        self.ax4.title.set_text('Norm dIdV Plot at (x,y) = ' + str(self.coordinate))
        self.ax4.set_xlabel('Voltage (V)')
        self.ax4.set_ylabel('dIdv (nA/V)')

        self.ax1.format_coord = lambda x, y: format(x, '1.4f') + ', ' + format(y, '1.4f')
        self.ax3.format_coord = lambda x, y: format(x, '1.4f') + ', ' + format(y, '1.4f')

        self.ax1.set_aspect(1)
        self.ax3.set_aspect(1)

        # this starts sensing for mouse clicks over the figure
        self.fig.canvas.mpl_connect('button_press_event', self.click_on_image)
        self.fig.canvas.mpl_connect('button_release_event', self.click_on_image)

        # This is where most of the GUI magic happens
        # Creates the frame for the widget to go in, then packs the widget into the frame
        self.options_frame = Frame(self.master)
        self.options_frame.pack(side=LEFT)

        self.tabs = Notebook(self.master)                     # Makes the tabs at the top of the figure
        self.tabs.pack(side=LEFT, expand=True, fill='both')

        self.canvas_frame = Frame(self.tabs, bg='')                    # Makes the frame for the figure to go in
        self.canvas_frame.pack(side=LEFT, expand=True, fill='both')
        self.tabs.add(self.canvas_frame)

        self.canvas = FigureCanvasTkAgg(self.fig, self.canvas_frame)   # Converts figure to a canvas widget so it can be used by tkinter
        self.toolbar = CustomToolbar(self.canvas, self.canvas_frame)   # replaces original toolbar with a custom one, don't remember why this is important
        self.canvas.get_tk_widget().pack(side=LEFT, expand=True, fill='both')

        self.tab2_frame = Frame(self.tabs)
        self.tab2_frame.pack(side=LEFT, expand=True, fill='both')
        self.tabs.add(self.tab2_frame)

        '''
        Creates open file button
        '''

        open_file_frame = Frame(self.options_frame)
        open_file_frame.pack(expand=True, fill='x')
        self.open_file_button = Button(open_file_frame, text='Open STS File', command=self.select_sts_file)
        self.open_file_button.pack()

        '''
        Creates open file button
        '''

        open_file_frame = Frame(self.options_frame)
        open_file_frame.pack(expand=True, fill='x')
        self.open_file_button = Button(open_file_frame, text='Open Topography Image', command=self.select_image_file)
        self.open_file_button.pack()

        '''
        Creates the data type box
        '''

        data_type_frame = Frame(self.options_frame)
        data_type_frame.pack(expand=True, fill='x')
        self.data_type_label = Label(data_type_frame, text='Data Source')
        self.data_type_label.pack(expand=True, fill='x')
        self.data_type_box = Entry(data_type_frame)
        self.data_type_box.pack(side=LEFT, fill='x', expand=True)

        '''
        Creates the Coordinate Box
        '''

        coord_frame = Frame(self.options_frame)
        coord_frame.pack(expand=True, fill='x')
        coord_label = Label(coord_frame, text='Coordinate')
        coord_label.pack(expand=True, fill='x')
        self.coord_box = Entry(coord_frame)
        self.coord_box.pack(side=LEFT, expand=True, fill='x')
        self.coord_box.insert(0, '(0,0)')

        '''
        Creates the color bar limits boxes
        '''

        colorbar_lim_frame = Frame(self.options_frame)
        colorbar_lim_frame.pack(expand=True, fill='x')
        colorbar_label = Label(colorbar_lim_frame, text='Colorbar Scale')
        colorbar_label.pack(expand=True, fill='x')
        colorbar_lower_label = Label(colorbar_lim_frame, text='Lower Limit:')
        colorbar_lower_label.pack(expand=True, fill='x', side=LEFT)
        self.colorbar_lower_box = Entry(colorbar_lim_frame, width=4)
        self.colorbar_lower_box.pack(side=LEFT, expand=True, fill='x')
        colorbar_upper_label = Label(colorbar_lim_frame, text='Upper Limit:')
        colorbar_upper_label.pack(side=LEFT, expand=True, fill='x')
        self.colorbar_higher_box = Entry(colorbar_lim_frame, width=4)
        self.colorbar_higher_box.pack(side=RIGHT, expand=True, fill='x')
        self.colorbar_higher_box.insert(0, '1')
        self.colorbar_lower_box.insert(0, '0')

        average_area_frame = Frame(self.options_frame)
        average_area_frame.pack(expand=True, fill='x')
        self.avg_area_state = tkinter.IntVar()
        self.average_area_check = Checkbutton(average_area_frame, command=self.average_area, text='Average Area' ,
                                              variable= self.avg_area_state)
        self.average_area_check.pack(side=LEFT, fill='x', expand = True)

        self.average_init = 0

        def close(self):
            '''
            closes the tkinter window
            :param self:
            :return:
            '''
            self.master.destroy()
            self.master.quit()

        self.master.protocol("WM_DELETE_WINDOW", lambda: close(self))

    def initialize(self, filepath):
        '''
        This function opens the selected sts file and calls a number of other functions that prepare
        that process the data so that it will run faster. This function will also print the time taken to complete each step
        so that you can check to see if a particular function is taking a long time. If this fails it is typically an
        issue with one of the functions it calls.

        :param filepath:
        :return:
        '''
        raw = self.get_data(filepath)                                                # gets the raw data
        self.data_type = self.data_type_box.get()                                    # this gets the value from the data type box. If the data is in a different format you will need to make a new parse function AND add it as an option here
        start = timer()
        if self.data_type == 'ORNL':
            data = self.parse_ORNL_data(raw)
        else:
            data = self.parse_data(raw)
        end = timer()
        print('Time to Parse: ', timedelta(seconds=end - start))
        start = timer()
        self.seperate_voltages(data)
        end = timer()
        self.init_lines()
        print('Time to seperate_voltages: ', timedelta(seconds=end - start))
        start = timer()
        self.generate_lines()
        end = timer()
        print('Time to Generate Lines: ', timedelta(seconds=end - start))

    def select_sts_file(self):
        '''
        Opens filedialog and allows the user to select an sts file. This program is designed for txt files. Also
        initializes some of the boxes and fills them with default values.
        :return:
        '''
        filetypes = (('text files', '*.txt'), ('All files', '*.*'))                                  # specifies the allowable filetypes
        filename = fd.askopenfilename(title='Open a file', initialdir='/', filetypes=filetypes)      # opens filedialog box

        if self.initialized == False:                                                                # checks to see if a file has already been initialized
            self.initialize(filename)                                                                # runs initialize
            self.plot_didv()                                                                         # plots the didv curve for (0,0)
            if self.data_type != 'ORNL':                                                             # plots IV curve as long as it isn't lock in STS data (ORNL is)
                self.plot_IV()
            self.didv_image_plot()                                                                   # plots the didv image at the initial voltage

            '''
            Creates the voltage box and apply button
            '''

            self.volt_frame = Frame(self.options_frame)
            self.volt_frame.pack(expand=True, fill='x')
            volt_label = Label(self.volt_frame, text='Voltage')
            volt_label.pack(expand=True, fill='x')

            self.volt_min = min(self.v_to_xyi.keys())
            self.volt_max = max(self.v_to_xyi.keys())
            self.volt_step = round((self.volt_max - self.volt_min) / (len(self.v_to_xyi.keys()) - 1), 7)
            self.volt_slider = Scale(self.volt_frame, from_=self.volt_min, to_=self.volt_max, orient=HORIZONTAL,
                                     resolution=self.volt_step)
            self.volt_slider.pack(side=LEFT, expand=True, fill='x')

            apply_button_frame = Frame(self.options_frame)
            apply_button_frame.pack()

            self.apply_button = Button(apply_button_frame, text='Apply', command=self.apply)
            self.apply_button.grid(row=0, column=0)
        else:                                                                                                     # my intent is for this to clear the old data and open the new data if the user selects a new file
            self.initialize(filename)
            self.plot_didv()
        self.ax4.set_xlim(self.volt_min, self.volt_max)                                                           # sets the x axis of the two line plots
        #self.ax4.set_ylim(self.didv_norm_min, self.didv_norm_max)
        self.ax2.set_xlim(self.volt_min, self.volt_max)
        #self.ax2.set_ylim(self.i_min, self.i_max)

    def select_image_file(self):

        '''
        Allows the user to open an STM image that corresponds with the sts grid. This is especially helpful if you are trying to look at the sts spectra on specific
        features in the image. If you load in the image you can click on the feature and it will pull up the spectra for that point. This is intended to work with png files
        only.
        :return:
        '''

        filetypes = (('Image files', '*.png'), ('All files', '*.*'))                                   # specifies the filetypes that appear in the filedialog
        filename = fd.askopenfilename(title='Open a file', initialdir='/', filetypes=filetypes)        # opens the filedialog
        old_im = PIL.Image.open(filename)                                                              # Creates a PIL image using the file selected
        old_im = ImageOps.mirror(old_im)                                                               # for whatever reason, the images load in mirrored and rotated 270 degrees so these operations fix that
        old_im = old_im.rotate(270)
        '''
        I believe that this next section was necessary because I wanted the coordinates in the image to match 
        those of the sts grid. In order to do that I had to get the colorvalue from each pixel and recreate the image 
        as a matplotlib colormesh. This also allows for the color palette to be changed. 
        '''
        w,h = old_im.size                                                                              # gets the width and height of the image in pixels
        x = np.linspace(0,self.image_size,num=w)                                                       # creates a numpy array with values from 0 to the width
        c_grid = np.zeros((w,h))                                                                       # creates a numpy array in a grid that is the size of the image.
        for i in range(w):                                                                             # loops across the rows
            for j in range(h):                                                                         # loops through the columns
                c = old_im.getpixel((i,j))                                                             # gets the color value from the corresponding spot in the old image
                c_grid[i][j] = self.grayscale(c)                                                 # converts it to a black and white image
        self.ax1.pcolormesh(x, x, c_grid, cmap='YlOrBr_r', shading='nearest', picker=True)               # plots the image in ax1
        self.canvas.draw()                                                                             # updates the canvas

    def get_data(self, filepath):

        """
        Opens the file at the location specified and converts it to a raw string. This will always
        need to be parsed later. It is necessary to include the .txt file extension in order for it
        to open properly.
        """

        with open(filepath) as raw_txt:
            raw = raw_txt.readlines()
        return raw

    def parse_data(self, raw):

        """
        Parses the data into a usable format. Depending on the source of this data this may or may not work.
        This specific parser is designed for use with our lab's STS data. If this parser does not work, create a
        new one (do not change this one) after this function with the name parse_XXXXX_data() with XXXXX being
        source of the data. This also creates a dictionary, coords_to_iv which maps a coordinate to its corresponding v
        and i values.
        :param raw: This parameter is intended to be the result of the previous function, get_data(), raw is in this case
        a list of strings where each row is tab delimited and each row terminates in a \n
        """

        str_data = raw[4:]  # removes header from data, you may need to change this value depending on the number of header lines.
        x_data = []  #initializing the lists that we will sort the data into
        y_data = []
        v_data = []
        i_data = []

        for row in str_data:                                          # Loops through each row of the raw string
            x_str, y_str, v_str, i_str = row.strip('\n').split('\t')  # splits the values in the row
            x_val = float(x_str)  # converts values into the correct variable type
            y_val = float(y_str)
            x_data.append(x_val)  # accumulates x and y lists
            y_data.append(y_val)
            v_val = float(v_str)
            i_val = float(i_str)
            v_data.append(v_val)
            i_data.append(i_val)
            coord = (x_val, y_val)
            if coord not in self.coords_to_lines.keys():                  # initializes a dict that maps the coordinate to the iv curves for that coord
                self.coords_to_lines[coord] = [[v_val], [i_val], [], []]  # the if/else isto ensure that there is only one key per coordinate and that it maps to all IV values
            else:
                self.coords_to_lines[coord][0].append(v_val)              # if there is already a key for that coord it appends the v and i values to the lists for that coord
                self.coords_to_lines[coord][1].append(i_val)

        self.i_max = max(i_data)
        self.i_min = min(i_data)
        largest_didv = -1000
        smallest_didv = 1000
        largest_didv_norm = -10000
        smallest_didv_norm = 10000

        for coord in self.coords_to_lines.keys():   # for each coordinate, do a two point derivative approximation to get didv
            x = coord[0]
            y = coord[1]
            didv_vals, didv_norm_vals = self.tp_der_approx_at_point(x, y)  # does a two point derivative approximation
            didv_max = max(didv_vals)
            didv_min = min(didv_vals)
            didv_norm_max = max(didv_norm_vals)
            didv_norm_min = min(didv_norm_vals)
            if didv_norm_max > largest_didv_norm:  # I am tracking the min and max this way because you cant do max(___) of a dictionary
                largest_didv_norm = didv_norm_max
            if didv_norm_min < smallest_didv_norm:
                smallest_didv_norm = didv_norm_min
            if didv_max > largest_didv:
                largest_didv = didv_max
            if didv_min < smallest_didv:
                smallest_didv = didv_min
            self.coords_to_lines[coord][2] = didv_vals
            self.coords_to_lines[coord][3] = didv_norm_vals

        self.didv_norm_min = smallest_didv_norm
        self.didv_norm_max = largest_didv_norm
        self.didv_max = largest_didv
        self.didv_min = smallest_didv
        self.coord_step = x_data[1] - x_data[0]  # records the step size of the piezo
        data = [x_data, y_data, v_data, i_data]
        self.voltage = v_data[0]  # sets the current voltage to be the first voltage in the file
        self.volt_list = self.coords_to_lines[self.coordinate][0]  # records a list of voltage values
        return data

    def parse_ORNL_data(self, raw):

        """
        Parses the data into a usable format. Depending on the source of this data this may or may not work.
        This specific parser is designed for use with our lab's STS data. If this parser does not work, create a
        new one (do not change this one) after this function with the name parse_XXXXX_data() with XXXXX being
        source of the data. This also creates a dictionary, coords_to_iv which maps a coordinate to its corresponding z
        and f values. This is different from the previous parser because the ORNL STS uses the lock in technique which measures didv rather than IV. It also provides
        the x and y values in meters not nm so it converts to make the numbers easier to deal with.
        :param raw: This parameter is intended to be the result of the previous function, get_data(), raw is in this case
        a list of strings where each row is tab delimited and each row terminates in a \n
        """
        str_data = raw[4:]
        x_data = []
        y_data = []
        v_data = []
        didv_data = []
        largest_didv = -1000
        smallest_didv = 1000
        coords = []
        for row in str_data:
            x_str, y_str, v_str, didv_str = row.strip('\n').split('\t')
            x_val =10**9*float(x_str)                                        # These values are written given in meters not nanometers so here you convert to nm to make it easier to work with
            y_val = 10**9*float(y_str)
            #Decimal(10**9*float(y_str)).quantize(Decimal('.0001'), rounding=ROUND_UP)
            coord = (x_val, y_val)
            coords.append(coord)
            x_data.append(x_val)
            y_data.append(y_val)
            v_val = float(v_str)
            didv_val = float(didv_str)
            # didv_norm = didv_val / ((i_data[i] / v_data[i]) + shift)
            if didv_val < smallest_didv:
                smallest_didv = didv_val
            if didv_val > largest_didv:
                largest_didv = didv_val
            v_data.append(v_val)
            didv_data.append(didv_val)
            if coord not in self.coords_to_lines.keys():
                self.coords_to_lines[coord] = [[v_val], [didv_val], [didv_val], [didv_val]]
            else:
                self.coords_to_lines[coord][0].append(v_val)
                self.coords_to_lines[coord][1].append(didv_val)
                self.coords_to_lines[coord][2].append(didv_val)
                self.coords_to_lines[coord][3].append(didv_val)  # this should be changed once I figure out normalization for ORNL data
        self.didv_max = largest_didv
        self.didv_norm_max = largest_didv
        self.didv_min = smallest_didv
        self.didv_norm_min = smallest_didv
        self.coord_step = x_data[1] - x_data[0]
        data = [x_data, y_data, v_data, didv_data]
        self.voltage = v_data[0]
        self.coords = coords
        return data

    def apply(self):
        '''
        Calls functions to update the graphs when the apply button is hit
        :return:
        '''
        self.set_voltage()
        self.didv_image_plot()
        self.set_coord()


    def init_lines(self):
        '''
        Initiates the numpy arrays that will be used to store IV and didv curves
        :return:
        '''
        len_1d = len(self.one_dim_coords_x)
        self.iv_lines = np.zeros((len_1d, len_1d, self.v_len + 1, 2), dtype='float16')
        self.didv_lines = np.zeros((len_1d, len_1d, self.v_len + 1, 2), dtype='float16')

    def generate_lines(self):
        '''
        Generates the IV and didv curves and stores them in arrays. In order to allow for the quick plotting of many lines at once
        for the area selection, this required a somewhat complex data structure. self.iv_lines and self.didv_lines are 4 dimensional arrays. The
        first two dimensions represent the x and y coordinates of the STS. at each point in that grid there is an array that is the length of the
        number of total points in one iv curve (or didv curve). At each spot in that array there is an array with an x coordinate and a y coordinate.
        Essentially, the line is stores as an array of coordinates rather than an x list and y list. This dramatically reduces the plotting time for a
        collection of lines and is utilized by matplotlib's Collection class.
        :return:
        '''
        for i in range(len(self.one_dim_coords_x)):
            x = self.one_dim_coords_x[i]
            for j in range(len(self.one_dim_coords_x)):                                   # loops through the x and y coordinates of the grid
                y = self.one_dim_coords_y[j]
                coord = (x, y)
                v_vals, i_vals, didv_vals, didv_norm_vals = self.coords_to_lines[coord]   # extracts the data for each point

                for k in range(len(v_vals)):                                              # finds the length of each line and loops through that number
                    v = v_vals[k]
                    i_val = i_vals[k]
                    didv_norm = didv_norm_vals[k]
                    self.iv_lines[i][j][k][0] = v                                         # stores the values in the array so that they can be plotted quickly later
                    self.iv_lines[i][j][k][1] = i_val
                    self.didv_lines[i][j][k][0] = v
                    self.didv_lines[i][j][k][1] = didv_norm


    def set_voltage(self):
        '''
        Gets the value from the voltage slider and updates the current voltage to be the voltage of the slider if they are different, then updates the didv image
        with the new voltage
        :return:
        '''
        old_volt = self.voltage
        v_slider = self.volt_slider.get()
        self.voltage = float(self.nearest_voltage(v_slider))

        if old_volt != self.voltage:
            self.didv_image_plot()

    def nearest_voltage(self, voltage):
        '''

        takes the voltage from the slider and finds the closest match in the data.

        :param voltage:
        :return:
        '''
        nearest_d = 9999
        closest = 0
        for v in self.v_to_xyi.keys():
            dist = abs(voltage - v)
            if dist < nearest_d:
                nearest_d = dist
                closest = v

        return closest

    def set_coord(self):
        '''
        Pulls the coordinate from the coordinate box and sets the current coordinate to that value.
        Plots the IV and didv curves for the new coordinate
        :return:
        '''

        new_coord_str = self.coord_box.get()
        new_coord = self.str_to_coord(new_coord_str)
        new_coord = self.nearest_coord(new_coord)
        self.coordinate = new_coord

        if self.circ1 != None:
            self.circ1.remove()
            self.circ3.remove()
            self.circ1 = None
            print(new_coord)
        self.circ1 = Circle(new_coord, radius=max(self.one_dim_coords_x) / 50, color='black')
        self.circ3 = Circle(new_coord, radius=max(self.one_dim_coords_x) / 50, color='black')
        self.ax3.add_patch(self.circ3)
        self.ax1.add_patch(self.circ1)


        if self.data_type != 'ORNL':
            self.plot_IV()                           # the lock in/ORNL does not have an IV plot so this keeps it from throwing an erro
        self.plot_didv()


    def str_to_coord(self, coord_string):
        '''
        Converts a acoordinate string formatted like (num1, num2) into a tuple with the corresponding numerical values
        This allows the user to input a strong into the coordinate box and for the program to use the coordinates as
        numerical values.
        :param coord_string:
        :return:
        '''

        no_parenths = coord_string[1:len(coord_string) - 1]
        str_vals = no_parenths.split(', ')
        x = float(str_vals[0])
        y = float(str_vals[1])
        coord = (x, y)
        return coord

    def coord_to_str(self, x, y):
        '''
        Converts a numerical coordinate (tuple )formatted like: (num1, num2) into a string
        :param x:
        :param y:
        :return:
        '''
        if x < 10**-6:
            coord_str = '(' + '{: .4f}'.format(x*10**9) + ',' + '{: .4f}'.format(y*10**9) + ')'
        else:
            coord_str = '(' + '{: .4f}'.format(x) + ',' + '{: .4f}'.format(y) + ')'
        return coord_str

    def get_dIdV(self, coord, v):
        '''
        gets the didv and normalized didv values at a particular coordinate and voltage
        :param coord:
        :param v:
        :return:
        '''

        v_data, i_data, didv_data, didv_norm_data = self.coords_to_lines[coord]
        ind = v_data.index(v)
        didv = didv_data[ind]
        didv_norm = didv_norm_data[ind]

        return [didv, didv_norm]

    def seperate_voltages(self, data):
        '''
        If you open the sts txt file you can see the way that the data is stored across the grid. At a given voltage first x is incremented from 0 to the max value,
        then y is increased, and it increments through the x values again and this process repeats to scan across the entire grid. Once the grid is complete the voltage
        is incremented and the process repeats. Although this is not the way that the STM actually measures, it is the way that the data is stored. This function splits
        the dataset at each point that the voltage changes then creates a dictionary where a given voltage maps to the x and y coordinates and the current values at that voltage
        for those coordinates.

        :param data:
        :return:
        '''

        x_data, y_data, v_data, i_data = data
        self.v_len = len(v_data)
        v_last = v_data[0]
        first_ind = 0
        first_run = True
        for v in v_data:
            if v_last != v:
                ind = v_data.index(v)
                x_vals = x_data[first_ind:ind]
                y_vals = y_data[first_ind:ind]
                didv_vals = []
                didv_norm_vals = []
                for i in range(0, len(x_vals)):
                    x = x_vals[i]
                    y = y_vals[i]
                    coord = (x, y)
                    didv, didv_norm = self.get_dIdV(coord, v)
                    didv_vals.append(didv)
                    didv_norm_vals.append(didv_norm)
                dataset = [didv_vals, didv_norm_vals]
                self.v_to_xyi[v_last] = dataset
                first_ind = ind
                if first_run == True:
                    self.x_vals = x_vals
                    self.y_vals = y_vals
                    self.one_dim_coords_x = x_vals[0:int(len(x_vals) ** .5)]      # generates an array that just has the one dimensional coordinates and not a grid
                    self.image_size = max(self.one_dim_coords_x)
                    self.one_dim_coords_y = np.unique(np.array(y_vals)).tolist()
                    first_run = False

            v_last = v

        self.v_len = len(self.v_to_xyi.keys())
        return

    def tp_der_approx_at_point(self, x, y):  # two point derivative approximation at a single point across all voltages
        v_data, i_data, trash1, trash2 = self.coords_to_lines[(x, y)]
        if min(v_data) <= 0:
            shift = -min(v_data) + .5
        else:
            shift = 0

        if self.data_type == 'ORNL':
            didv_data = i_data
            didv_norm_data = i_data
        else:
            didv_data = []
            didv_norm_data = []
            for i in range(len(v_data)):
                if i == 0:  # forward difference for the first point
                    v1 = v_data[i]
                    i1 = i_data[i]
                    i2 = i_data[i + 1]
                    v2 = v_data[i + 1]
                elif i == len(v_data) - 1:  # backward difference for the last point
                    v1 = v_data[i]
                    i1 = i_data[i]
                    i2 = i_data[i - 1]
                    v2 = v_data[i - 1]
                else:  # centered difference for all of the other points
                    v1 = v_data[i - 1]
                    i1 = i_data[i - 1]
                    i2 = i_data[i + 1]
                    v2 = v_data[i + 1]
                didv = (i2 - i1) / (v2 - v1)
                didv_norm = didv / (( i_data[i] / v_data[i]) + shift)
                if didv_norm > 10 or didv_norm < -10:
                    didv_norm = 0
                didv_norm_data.append(didv_norm)
                didv_data.append(didv)
        return [didv_data, didv_norm_data]

    def didv_image_plot(self):
        self.ax3.cla()

        volt = self.voltage
        didv_data, didv_norm_data = self.v_to_xyi[volt]
        didv_norm_ar = np.array(didv_norm_data)

        didv_norm_grid = didv_norm_ar.reshape(len(self.one_dim_coords_x),len(self.one_dim_coords_y))
        self.color_max = float(self.colorbar_higher_box.get())
        self.color_min = float(self.colorbar_lower_box.get())
        self.im3 = self.ax3.pcolormesh(self.one_dim_coords_x, self.one_dim_coords_y, didv_norm_grid, shading='nearest', picker=True, vmin=self.color_min,
                                       vmax=self.color_max, cmap= 'YlOrBr_r')

        if self.colorbar_init == False:
            self.cb3 = self.fig.colorbar(self.im3, ax=self.ax3, orientation='vertical')
            self.colorbar_init = True
        else:
            self.cb3.update_normal(self.im3)
        self.canvas.draw()

    def plot_IV(self):
        if self.rectangle == False:
            self.ax2.cla()
            self.ax2.grid()
            self.ax2.title.set_text('IV Plot at (x,y) = ' + self.coord_to_str(self.coordinate[0], self.coordinate[1]))
            self.ax2.set_ylabel('Current (nA)')
        v_data, i_data, didv_vals, didv_norm_vals = self.coords_to_lines[self.coordinate]
        self.ax2.plot(v_data, i_data, c='#1f77b4')
        if self.rectangle == False:
            self.canvas.draw()
        try:
            self.ax4.set_xlim(self.volt_min, self.volt_max)
            #self.ax4.set_ylim(self.didv_norm_min, self.didv_norm_max)
            self.ax2.set_xlim(self.volt_min, self.volt_max)
            #self.ax2.set_ylim(self.i_min, self.i_max)
        except:
            pass

    def sort_by_position(self):
        with open('sorted_by_pos.csv', 'w') as csvfile:
            datawriter = csv.writer(csvfile, delimiter=',')
            for coord in self.coords_to_lines.keys():
                x, y = coord
                v_data, i_data = self.coords_to_lines[coord]
                for i in range(len(v_data)):
                    v_val = v_data[i]
                    i_val = i_data[i]
                    line = [x, y, v_val, i_val]
                    datawriter.writerow(line)

    def scatter_IV(self):
        raw_data = self.get_data()
        x_data, y_data, v_data, i_data = self.parse_ORNL_data(raw_data)
        plt.scatter(v_data, i_data, c='b')
        plt.grid()
        plt.xlabel('Voltage')
        plt.ylabel('Current')
        plt.show()

    def plot_didv(self):
        if self.rectangle == False:
            self.ax4.cla()
            self.ax4.grid()
            self.ax4.title.set_text('Norm dIdV Plot at (x,y) = ' + self.coord_to_str(self.coordinate[0], self.coordinate[1]))
            self.ax4.set_xlabel('Voltage (V)')
            self.ax4.set_ylabel('dIdv (nA/V)')
        try:
            self.ax4.set_xlim(self.volt_min, self.volt_max)
            #self.ax4.set_ylim(self.didv_norm_min, self.didv_norm_max)
            self.ax2.set_xlim(self.volt_min, self.volt_max)
            #self.ax2.set_ylim(self.i_min, self.i_max)
        except:
            pass
        v_data, i_data, didv_vals, didv_norm_vals = self.coords_to_lines[self.coordinate]
        self.ax4.plot(v_data, didv_vals, c='#1f77b4')
        #self.band_gap()
        if self.rectangle == False:
            self.canvas.draw()

    def nearest_coord(self, c):
        x = c[0]
        y = c[1]
        closest_x = 0
        best_xdist = 9999999
        closest_y = 0
        best_ydist = 9999999
        for i in range(len(self.one_dim_coords_y)):
            x2 = self.one_dim_coords_x[i]
            y2 = self.one_dim_coords_y[i]
            x_dist = abs(x2 - x)
            if x_dist < best_xdist:
                best_xdist = x_dist
                closest_x = x2
            y_dist = abs(y2 - y)
            if y_dist < best_ydist:
                best_ydist = y_dist
                closest_y = y2

        return (closest_x,closest_y)

    def reset_ax(self, ax):

        if ax == self.ax2:
            if self.data_type != 'ORNL':
                self.ax2.cla()
                self.ax2.grid()
                self.ax2.title.set_text(
                    'IV Plot at (x,y) = ' + self.coord_to_str(self.coordinate[0], self.coordinate[1]))
                self.ax2.set_ylabel('Current (nA)')

        if ax == self.ax4:
            self.ax4.cla()
            self.ax4.grid()
            self.ax4.title.set_text(
                'dIdV Plot at (x,y) = ' + self.coord_to_str(self.coordinate[0], self.coordinate[1]))
            self.ax4.set_xlabel('Voltage (V)')
            self.ax4.set_ylabel('dIdv (nA/V)')

    def average_area(self):
        if self.avg_area_state.get() == 1:
            i_sums = np.zeros(len(self.volt_list))
            didv_sums = np.zeros(len(self.volt_list))
            num_cells = 0
            for row in self.line_data:
                for cell in row:
                    ivline, didv_line = cell
                    num_cells += 1
                    for i in range(len(ivline)):
                        didv_val = didv_line[i][1]
                        i_val = ivline[i][1]
                        i_sums[i] += i_val
                        didv_sums[i] += didv_val
            i_avg = i_sums/num_cells
            didv_avg = didv_sums/num_cells
            self.average_i_line = self.ax2.plot(self.volt_list,i_avg, c='r')
            self.average_didv_line = self.ax4.plot(self.volt_list, didv_avg, c='r')
            self.canvas.draw()
            self.average_init = 1

        if self.avg_area_state.get() == 0 and self.average_init == 1:
            pop1 = self.average_i_line.pop()
            pop2 = self.average_didv_line.pop()
            pop1.remove()
            pop2.remove()
            self.canvas.draw()

    def click_on_image(self, event):
        if event.inaxes == self.ax1 or event.inaxes == self.ax3:
            if event.name == 'button_press_event':
                self.last_click = event
            if event.name == 'button_release_event':
                self.last_release = event

            if self.unclicked == True:
                xdata = event.xdata
                ydata = event.ydata
                click_point = (xdata, ydata)
                image_point = self.nearest_coord(click_point)
                self.coordinate = image_point
                self.coord_box.delete(0, END)
                self.coord_box.insert(0, self.coord_to_str(image_point[0], image_point[1]))
                if self.data_type != 'ORNL':
                    self.plot_IV()
                self.plot_didv()
                self.unclicked = False
            else:
                if self.avg_area_state.get() == 1 and self.average_init == 1:
                    self.average_area_check.toggle()
                click_coord = (self.last_click.xdata, self.last_click.ydata)
                release_coord = (self.last_release.xdata, self.last_release.ydata)
                if click_coord != release_coord and event.name != 'button_press_event':
                    self.reset_ax(self.ax2)
                    self.reset_ax(self.ax4)
                    self.rectangle = True
                    nearest_to_click = self.nearest_coord(click_coord)
                    nearest_to_release = self.nearest_coord(release_coord)
                    x_bounds = [nearest_to_click[0], nearest_to_release[0]]
                    y_bounds = [nearest_to_click[1], nearest_to_release[1]]
                    x_min = self.one_dim_coords_x.index(min(x_bounds))
                    x_max = self.one_dim_coords_x.index(max(x_bounds))
                    y_min = self.one_dim_coords_y.index(min(y_bounds))
                    y_max = self.one_dim_coords_y.index(max(y_bounds))
                    if self.circ1 != None:
                        self.circ1.remove()
                        self.circ3.remove()
                        self.circ1 = None
                    if self.rect1 != None:
                        self.rect1.remove()
                        self.rect3.remove()
                    rect_c = (min(x_bounds), min(y_bounds))
                    rect_w = max(x_bounds) - min(x_bounds)
                    rect_h = max(y_bounds) - min(y_bounds)
                    self.rect1 = Rectangle(rect_c, rect_w, rect_h, edgecolor= (0,0,0,1),
                                          linewidth=1, facecolor= (0,0,0,.25))
                    self.rect3 = Rectangle(rect_c, rect_w, rect_h, edgecolor=(0, 0, 0, 1),
                                           linewidth=1, facecolor=(0, 0, 0, .25))
                    self.ax3.add_patch(self.rect3)
                    self.ax1.add_patch(self.rect1)
                    x_vals = np.arange(x_min, x_max)
                    y_vals = np.arange(y_min, y_max)
                    ixgrid = np.ix_(x_vals, y_vals)
                    iv_line_grid = self.iv_lines[ixgrid]
                    didv_line_grid = self.didv_lines[ixgrid]
                    iv_line_data = []
                    didv_line_data = []
                    for row in iv_line_grid:
                        for line in row:
                            iv_line_data.append(line)
                    for row in didv_line_grid:
                            for line in row:
                                didv_line_data.append(line)

                    iv_plot = matplotlib.collections.LineCollection(iv_line_data)
                    didv_plot = matplotlib.collections.LineCollection(didv_line_data)
                    self.ax4.add_collection(didv_plot)
                    if self.data_type != 'ORNL':
                        self.ax2.add_collection(iv_plot)
                    self.ax4.autoscale()
                    self.ax2.autoscale()

                    self.canvas.draw()
                    self.rectangle = False


                if click_coord == release_coord:
                    if self.rect1 != None:
                        self.rect1.remove()
                        self.rect3.remove()
                        self.rect1 = None
                    if self.circ1 != None:
                        self.circ1.remove()
                        self.circ3.remove()
                    xdata = event.xdata
                    ydata = event.ydata
                    click_point = (xdata, ydata)
                    image_point = self.nearest_coord(click_point)
                    self.circ1 = Circle(image_point, radius=max(self.one_dim_coords_x)/50, color='black')
                    self.circ3 = Circle(image_point, radius=max(self.one_dim_coords_x)/50, color='black')
                    self.ax1.add_patch(self.circ1)
                    self.ax3.add_patch(self.circ3)
                    self.coordinate = image_point
                    self.coord_box.delete(0, END)
                    self.coord_box.insert(0, self.coord_to_str(image_point[0], image_point[1]))
                    if self.data_type != 'ORNL':
                        self.plot_IV()
                    self.plot_didv()

    def grayscale(self, color):
        r = color[0]
        g = color[1]
        b = color[2]
        grayscale = int(r*.299+g*.587+b*.114)
        return grayscale

    def band_gap(self):
        '''
        current thoughts on the process here:
        split the values into two halves, take the maximum values of each half then draw a line between that value
        and the point where the background line intersects the didv line on each side
        :return:
        '''
        v_data, i_data, didv_vals, didv_norm_vals = self.coords_to_lines[self.coordinate]
        max_didv = max(didv_norm_vals)
        left_side = didv_norm_vals[0:len(didv_norm_vals)//2]
        left_arr = np.array(left_side)
        right_side = didv_norm_vals[len(didv_norm_vals)//2:]
        right_arr = np.array(right_side)
        left_max = max(left_side)
        right_max = max(right_side)
        background_arr_left = np.array([max_didv * .08] * len(left_side))
        background_arr_right = np.array([max_didv * .08] * len(right_side))
        idx_left = np.argwhere(np.diff(np.sign(left_arr - background_arr_left))).flatten()
        idx_right = np.argwhere(np.diff(np.sign(right_arr - background_arr_right))).flatten()
        left_vals = [left_max,max_didv*.08]
        cond_band_min = v_data[idx_left[-1]]
        left_x = [v_data[left_side.index(left_max)],cond_band_min]
        right_vals = [max_didv*.08,right_max]
        val_band_max = v_data[idx_right[0]+len(didv_norm_vals)//2]
        right_x = [val_band_max, v_data[right_side.index(right_max)+len(didv_norm_vals)//2]]
        background_line = [max_didv*.08]*2
        background_line_x = [v_data[0], v_data[-1]]
        self.ax4.plot(background_line_x, background_line, c='red')
        self.ax4.plot(left_x, left_vals, c='red')
        self.ax4.plot(right_x, right_vals, c='red')
        print('Valence Band Maximum:',val_band_max)
        print('Conduction Band Minimum:', cond_band_min)
        print('Band Gap:', cond_band_min-val_band_max)





if __name__ == "__main__":
    interface()
