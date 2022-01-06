#Author: Charles Giese


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mesa_reader as mr
import os
import pandas as pd
import astropy.constant as c

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
        plt.xlabel('log(T) (K)')
        plt.ylabel('log(L) (L/L_solar)')
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

    def energy_production(self):

        history_data = self.load_history()

        LH = history_data.log_LH
        LHe = history_data.log_LHe

        pp = history_data.pp
        cno = history_data.cno
        tri_alfa = history_data.tri_alfa

        age = history_data.star_age

        fig, axs = plt.subplots(2,1)

        for ax in axs:
            #   ax.set_xscale('log')
            ax.set_yscale('log')

        axs[0].plot(age, LH, label = 'H Luminosity', c='k')
        axs[0].plot(age, LHe, label = 'He Luminosity', c='r')

        axs[0].legend()
        axs[0].set_title('Energy Production$')
        axs[0].set_ylabel('$Log(L/L_{\odot})$')


        axs[1].plot(age, pp, label = 'PP Cycle', c='k')
        axs[1].plot(age, cno, label = 'CNO Cycle', c='r')
        axs[1].plot(age, tri_alfa, label = 'Triple Alpha', c='b')
        axs[1].legend()
        axs[1].set_ylabel('$Log(L/L_{\odot})$')
        axs[1].set_xlabel('Stellar Age (Years)')

    def abundances(self):

        h = self.load_history()

        center_h1 = h.center_h1
        center_he4 = h.center_he4
        center_c12 = h.center_c12
        center_n14 = h.center_n14
        center_o16 = h.center_o16
        center_fe56 = h.center_fe56

        surface_h1 = h.surface_h1
        surface_he4 = h.surface_he4
        surface_c12 = h.surface_c12
        surface_n14 = h.surface_n14
        surface_o16 = h.surface_o16

        star_age = h.star_age

        fig = plt.figure()
        axc = fig.add_subplot(121)
        axs = fig.add_subplot(122)
        axc.set_xscale('log')
        axs.set_xscale('log')
        axc.set_yscale('log')
        axs.set_yscale('log')
        axc.set_xlabel('Stellar Age (Years)')
        axs.set_xlabel('Stellar Age (Years)')
        axc.set_ylabel('Abundance')

        axc.plot(star_age, center_h1, label = 'H1, Centre')
        axc.plot(star_age, center_he4, label = 'He4, Centre')
        axc.plot(star_age, center_c12, label = 'C12, Centre')
        axc.plot(star_age, center_n14, label = 'N14, Centre')
        axc.plot(star_age, center_o16, label = 'O16, Centre')
        axc.plot(star_age, center_fe56, label = 'Fe56, Centre')
        axc.set_ylim(1e-5, 1)

        axs.plot(star_age, surface_h1, label = 'H1, Surface')
        axs.plot(star_age, surface_he4, label = 'He4, Surface')
        axs.plot(star_age, surface_c12, label = 'C12, Surface')
        axs.plot(star_age, surface_n14, label = 'N14, Surface')
        axs.plot(star_age, surface_o16, label = 'O16, Surface')


        axc.legend()
        axs.legend()


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

    def energy_transport():

        df = self.load_profile()

        mass = df['mass']
        conv_vel = df['conv_vel']
        opacity = df['opacity']
        pressure = 10**df['logP']
        temp = 10**df['logT']

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(mass, conv_vel)
        ax.set_xscale('Mass (m/M)')
        
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
        plt.xlabel('Mass (m/M)')
        plt.ylabel('opacity')

    def lum_mass(self):
        df = self.load_profile()

        lum = np.exp(df['logL'])
        mass = df['mass']

        plt.plot(mass, lum, c='k')
        plt.xlabel('Mass (m/M)')
        plt.ylabel('l')

class ZAMS_TAMS():

    def __init__(self, mesa_log_dir = os.getcwd()+'/LOGS'):
        self.mesa_log_dir = mesa_log_dir

    def calc_points(self):

        import os

        path, dirs, files = next(os.walk(self.mesa_log_dir))
        file_count = len(files) -3

        #h_power = np.zeros(shape=(file_count))
        #power_nuc = np.zeros_like(h_power)
        #star_age = np.zeros_like(h_power)
        #r = np.zeros_like(h_power)

        h = history_plot()
        data = h.load_history()
        central_H1 = data.center_h1
        pp = data.pp
        cno = data.cno
        star_age = data.star_age
        model_number  = data.model_number


        """
        for i in range(1, file_count):
            profile = profile_plot(mesa_profile=i)
            meta_data = profile.load_metadata()
            data = profile.load_profile()
            h_power[i-1] = meta_data['power_h_burn'].values
            star_age[i-1] = meta_data['star_age'].values
            power_nuc[i-1] = meta_data['power_nuc_burn'].values
            r[i-1] = (data['radius'].values)[0]


        indices = np.array([0,0])
        for p in range(0, file_count-1):
            if np.abs(h_power[p] - h_power[p+1]) > h_power[p+1]/10 and (h_power[p+1] - h_power[p]) > 0:
                print('ZAMS is profile %i' %p)
                indices[0] = p
                break
        for p in range(0, file_count-1):
            if np.abs(h_power[p] - h_power[p+1]) > h_power[p+1]/10 and (h_power[p+1] - h_power[p]) < 0:
                print('TAMS is profile %i' %p)
                indices[1] = p
                break
        """


        fig = plt.figure()
        ax = fig.add_subplot(211)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.plot(model_number, central_H1, c='k')
        ax.set_ylabel('Central H1 Abundance')
        ax.set_xlabel('Model Number')


        ax2 = fig.add_subplot(212)
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.plot(model_number, pp, label = 'PP Chain')
        ax2.plot(model_number, cno, label = 'CNO Cycle')
        ax2.set_xlabel('Model Number')
        ax2.legend()
