#Author: Charles Giese


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mesa_reader as mr
import os

plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['font.size'] = 10
plt.rcParams['figure.titlesize'] = 18
plt.rcParams['axes.labelsize'] = 12






class history_plot():
    """Class for plotting outputs from MESA 1-D stellar evolution code"""

    def __init__(self, mesa_log_dir = os.getcwd()+'/LOGS', mesa_profile = None):
        self.mesa_log_dir = mesa_log_dir


    def load_history(self):
        """Load History data"""

        history_data = mr.MesaData(self.mesa_log_dir+'/history.data')
        return history_data

    def plot_hr_diagram(self):

        history_data = self.load_history()
        star_mass = round(history_data.star_mass[0], 1)

        L = history_data.log_L
        T_eff = history_data.log_Teff

        plt.plot(T_eff, L, c = 'k')
        plt.xlim(plt.xlim()[::-1])
        plt.xlabel('log(T)')
        plt.ylabel('log(L)')
        plt.title(' %s Solar Mass HR Diagram' % star_mass)

    def plot_T_rho(self):
        """This function creates an evolutionary track in the density-temperature plane"""

        history_data = self.load_history()
        star_mass = round(history_data.star_mass[0], 1)

        T_c = history_data.log_center_T
        rho_c = history_data.log_center_Rho

        plt.plot(rho_c, T_c, c='k')
        plt.xlabel('rho_C')
        plt.ylabel('T_C')
        plt.title(' %s Solar Masses' %star_mass)
