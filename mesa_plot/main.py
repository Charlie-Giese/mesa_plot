#Author: Charles Giese


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mesa_reader as mr
import os
import pandas as pd
import astropy.constants as c

plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['font.size'] = 14
plt.rcParams['figure.titlesize'] = 18
plt.rcParams['figure.figsize'] = [7.5,5.5]

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

        plt.plot(T_eff, L, c = 'k', label = (r'$%.1f M_{\odot}$' %star_mass))
        plt.xlim(plt.xlim()[::-1])
        plt.xlabel(r'$log(T) (K)$')
        plt.ylabel(r'$log(L/L_{\odot})$')
        plt.legend(fontsize = 14)


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
        Lgrav = history_data.log_abs_Lgrav

        pp = history_data.pp
        cno = history_data.cno
        tri_alfa = history_data.tri_alfa

        age = history_data.star_age
        model = history_data.model_number

        fig, axs = plt.subplots(2,1)

        for ax in axs:
            #ax.set_xscale('log')
            ax.set_yscale('log')

        #axs[0].plot(age, 10**LH, label = r'$log(L_H/L_{\odot})$', c='k')
        axs[0].plot(age, 10**LHe, label = r'$log(L_{He}/L_{\odot})$', c='r')
        #axs[0].plot(age, 10**Lgrav, label = r'$log(L_{grav}/L_{\odot})$', c='b')
        #axs[0].set_ylim(1e-20, 1e3)
        #axs[0].set_ylim(1e-5, 1e10)
        #axs[0].set_xlim(1.2e7, 1.4e7)

        axs[0].legend(fontsize=14)
        axs[0].set_ylabel('$L/L_{\odot}$')

        """
        axs[1].plot(age, 10**pp, label = 'PP Cycle', c='k')
        axs[1].plot(age, 10**cno, label = 'CNO Cycle', c='r')
        axs[1].plot(age, 10**tri_alfa, label = 'Triple Alpha', c='b')
        axs[1].legend(fontsize=14)
        axs[1].set_ylabel('$erg/s/g$')
        axs[1].set_xlabel('Stellar Age (Years)')
        axs[1].set_ylim(1e-6, 1e8)
        axs[1].set_xlim(1.2e7, 1.4e7)

        axs[0].plot(np.ones_like(age)*age[595-1], np.linspace(1e-10, 1e10, 3092), '--')
        axs[1].plot(np.ones_like(age)*age[595-1], np.linspace(1e-10, 1e10, 3092), '--')
        """


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
        #axc = fig.add_subplot(111)
        axs = fig.add_subplot(111)
        #axc.set_xscale('log')
        #axs.set_xscale('log')
        #axc.set_yscale('log')
        axs.set_yscale('log')
        #axc.set_xlabel('Stellar Age (Years)')
        axs.set_xlabel('Stellar Age (Years)')
        axs.set_ylabel('Abundance')

        #axc.plot(star_age, center_h1, label = 'H1, Centre')
        #axc.plot(star_age, center_he4, label = 'He4, Centre')
        #axc.plot(star_age, center_c12, label = 'C12, Centre')
        #axc.plot(star_age, center_n14, label = 'N14, Centre')
        #axc.plot(star_age, center_o16, label = 'O16, Centre')
        #axc.plot(star_age, center_fe56, label = 'Fe56, Centre')
        #axc.set_ylim(1e-5, 1)

        axs.plot(star_age, surface_h1, label = 'H1, Surface')
        axs.plot(star_age, surface_he4, label = 'He4, Surface')
        axs.plot(star_age, surface_c12, label = 'C12, Surface')
        axs.plot(star_age, surface_n14, label = 'N14, Surface')
        axs.plot(star_age, surface_o16, label = 'O16, Surface')
        axs.set_ylim(1e-5,1e1)


        #axc.legend(fontsize = 14)
        axs.legend(fontsize = 14)

    def zams_tams(self):

        h = self.load_history()
        model = h.model_number
        h1 = h.center_h1
        he4 = h.center_he4

        plt.plot(model, h1, label = 'H1', c = 'b', lw = 10)
        plt.plot(model, he4, label = 'He4', c = 'g', lw = 10)
        plt.xlabel('Model Number')
        plt.ylabel('Core Mass Fraction')


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
        dm = self.load_metadata()
        mass = df['mass'].values
        radius = df['radius'].values
        h1 = df['h1']
        he3 = df['he3']
        he4 = df['he4']
        c12 = df['he4']
        n14 = df['he4']
        o16 =   df['he4']


        labels = ['h1', 'he3', 'he4', 'c12', 'n14', 'o16']
        i=0
        for frac in [h1, he3, he4, c12, n14, o16]:
            #if np.min(frac) >= -100:
            plt.plot(mass, frac, label = labels[i])
            #plt.ylim(-5,0)
            plt.legend(fontsize = 14)
            #plt.xscale('log')
            plt.yscale('log')
            plt.xlabel(r'$M/M_{\odot}$', fontsize = 16)
            plt.ylabel(r'$Mass Fraction$', fontsize = 16)
            i+=1

    def rho_T_profile(self):

        df = self.load_profile()
        logT = df['logT']
        logRho = df['logRho']

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(logRho, logT, c='k')
        ax.set_xlabel(' log(Density)')
        ax.set_ylabel('$log(T)$')


    def energy_transport(self):

        df = self.load_profile()

        from astropy import constants as const
        import astropy.units as u

        m = df['mass'] * const.M_sun.cgs
        conv_vel = df['conv_vel']
        k = df['opacity']
        P = 10**df['logP']
        T = 10**df['logT']
        l = 10**df['logL'] * const.L_sun.cgs
        r = df['radius']
        import math

        def radiative_tempgrad(k, P, T, l, m):
            rad_const = (8*(math.pi)**5 * (const.k_B.cgs)**4)/(15 * (const.c.cgs)**3 * (const.h.cgs)**3)
            cons = 3/(const.c.cgs * const.G.cgs * 16 * math.pi * rad_const)
            rad_grad = cons * (k*l*P)/(m* T**4)
            return rad_grad


        #plt.plot(m/const.M_sun.cgs, radiative_tempgrad(k,P,T,l,m), label = r'$\nabla_{rad}$', c='k')
        #plt.plot(m/const.M_sun.cgs, np.ones_like(m)*0.4, label = r'$\nabla_{ad}$', c='r')
        plt.plot(r, radiative_tempgrad(k,P,T,l,m), label = r'$\nabla_{rad}$', c='k')
        plt.plot(r, np.ones_like(r)*0.4, label = r'$\nabla_{ad}$', c='r')
        plt.xlabel(r'$r/R_{\odot}$')
        plt.ylabel(r'$\nabla$')
        plt.legend(fontsize=14)
        plt.ylim(0,2)

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

    def mass_loss(self):

        df = self.load_profile()
        mass_loss = df['star_mdot'].values
        star_age = df['star_age'].values

        plt.plot(star_age, mass_loss, c='k')
        plt.xlabel('Stellar Age (Years)')
        plt.ylabel('Mass Loss')
        plt.xscale('log')
        plt.yscale('log')

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
