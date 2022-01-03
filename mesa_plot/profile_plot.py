#Author: Charles Giese


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mesa_reader as mr
import os
import pandas as pd




class profile_plot():
    """Class for plotting outputs from MESA 1-D stellar evolution code"""

    def __init__(self, mesa_log_dir = os.getcwd(), mesa_profile = None):
        self.mesa_log_dir = mesa_log_dir
        self.mesa_profile = mesa_profile


    def load_profile(self):
        """Load requisite profile"""

        df = pd.read_table("profile"+str(self.mesa_profile)+'.data', sep="\s+", header = 4)
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
