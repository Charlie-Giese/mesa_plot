#Author: Charles Giese


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mesa_reader as mr
import os
import pandas as pd

plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['font.size'] = 10
plt.rcParams['figure.titlesize'] = 18


class profile_plot():
    """Class for plotting outputs from MESA 1-D stellar evolution code"""

    def __init__(self, mesa_log_dir = os.getcwd()+'/LOGS', mesa_profile = None):
        self.mesa_log_dir = mesa_log_dir
        self.mesa_profile = mesa_profile


    def load_profile(self):
        """Load requisite profile"""

        log_dir = mr.MesaLogDir(self.mesa_log_dir)
        df = pd.read_table(self.mesa_log_dir+"/profile"+str(self.mesa_profile)+'.data', sep="\s+", header = 4)
        return df

    def load_metadata(self):

        df = pd.read_table(self.mesa_log_dir+'/profile'+str(self.mesa_profile)+'.data', sep="\s+", header = 1, nrows = 1)
        return df


    def plot_massfrac(self):
        """Plot profile of mass fraction inside star"""

        df = self.load_profile()
        mass = df['mass']
        h1 = np.log(df['h1'])
        he3 = np.log(df['he3'])
        he4 = np.log(df['he4'])
        c12 = np.log(df['he4'])
        n14 = np.log(df['he4'])
        o16 = np.log(df['he4'])


        labels = ['h1', 'he3', 'he4', 'c12', 'n14', 'o16']
        i=0
        for frac in [h1, he3, he4, c12, n14, o16]:
            #if np.min(frac) >= -100:
            plt.plot(mass, frac, label = labels[i])
            plt.ylim(-5,0)
            plt.legend()
            i+=1

    def rho_T_profile(self):

        df = self.load_profile()
        logT = df['logT']
        logRho = df['logRho']
        logL = df['logL']
        opacity = df['opacity']
        radius = df['radius']

        fig, axs = plt.subplots(2,1)
        axs[0].plot(radius, logT, label = 'Log(T)', c='r')
        axs[0].plot(radius, opacity, label = 'Opacity', c='k')

        axs[1].plot(radius, logL, label = 'Log(L)', c='k')
        axs[1].plot(radius, logRho, label = 'Log(\rho)', c='r')


    def P_rho(self):
        df = self.load_profile()

        pressure = np.exp(df['logP'])
        density = np.exp(df['logRho'])

        plt.plot(density, pressure, c='k')
        plt.xlabel('Density')
        plt.ylabel('Pressure')

    def mass_temp(self):
        df = self.load_profile()

        mass = df['mass']
        T = np.exp(df['logT'])
        plt.plot(mass, T, c='k')
        plt.xlabel('Mass')
        plt.ylabel('T')


    def opacity(self):
        df = self.load_profile()

        opacity = df['opacity']
        mass = df['mass']

        plt.plot(mass, opacity, c='k')
        plt.xlabel('mass')
        plt.ylabel('opacity')

    def lum_mass(self):
        df = self.load_profile()

        lum = np.exp(df['logL'])
        mass = df['mass']

        plt.plot(mass, lum, c='k')
        plt.xlabel('mass')
        plt.ylabel('l')
