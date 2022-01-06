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

    def energy_production(self):

        history_data = self.load_history()

        LH = history_data.log_LHe
        LHe = history_data.log_LHe
        LZ = history_data.log_LHe

        pp = history_data.pp
        cno = history_data.cno
        tri_alfa = history_data.tri_alfa

        age = history_data.star_age

        fig, axs = plt.subplots(2,1)
        axs[0].plot(age, LH, label = 'H Luminosity', c='k')
        axs[0].plot(age, LHe, label = 'He Luminosity', c='r')
        axs[0].plot(age, LZ, label = 'Z Luminosity', c='b')
        axs[0].legend()
        axs[0].set_title('Energy Production for 0.4 $M_{\odot}$')
        axs[0].set_ylabel('$Log(L/L_{\odot})$')

        axs[1].plot(age, pp, label = 'PP Cycle', c='k')
        axs[1].plot(age, cno, label = 'CNO Cycle', c='r')
        axs[1].plot(age, tri_alfa, label = 'Triple Alpha', c='b')
        axs[1].legend()
        axs[1].set_ylabel('$Log(L/L_{\odot})$')
        axs[1].set_xlabel('Stellar Age (Years)')

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

        h_power = np.zeros(shape=(file_count))
        power_nuc = np.zeros_like(h_power)
        star_age = np.zeros_like(h_power)

        for i in range(1, file_count):
            profile = profile_plot(mesa_profile=i)
            meta_data = profile.load_metadata()
            h_power[i-1] = meta_data['power_h_burn'].values
            star_age[i-1] = meta_data['star_age'].values
            power_nuc[i-1] = meta_data['power_nuc_burn'].values

        """
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
        axh = fig.add_subplot(211)
        axt = fig.add_subplot(212)
        axh.set_xscale('log')
        axh.set_yscale('log')
        axt.set_xscale('log')
        axt.set_yscale('log')
        axh.scatter(star_age, h_power, c='k')
        axt.scatter(star_age, power_nuc, c='r')
        axh.set_ylabel('Power')
        axt.set_xlabel('Stellar Age (Years)')
        axt.set_ylabel('Power')
