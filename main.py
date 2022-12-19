import numpy as np
import math as m
import mesa_web as mw
import matplotlib.pyplot as plt
import os


class star:
    """
    Class for a folder of data
    """

    def __init__(self, folder):
        self.folder = "./{dir}".format(dir=folder)
        historypath = "./{dir}/trimmed_history.data".format(dir=folder)
        self.nprofiles = len(os.listdir("./{dir}".format(dir=folder)))-5
        self.datadict = mw.read_history(historypath)

        self.log_L = []
        self.log_Teff = []
        self.star_age = []
        self.rho = []
        self.mass = []
        self.temperature = []
        self.pressure = []
        self.pp = []
        self.cno = []
        self.tri_alfa = []
        self.h1 = []
        self.he3 = []
        self.he4 = []
        self.he = []
        self.c12 = []
        self.o16 = []
        self.grada = []
        self.gradr = []
        self.rad = []
        self.conv = []
        self.ZAMS = []

        self.i1 = 0
        self.i2 = 0
        self.i3 = 0
        self.i4 = 0

        self.profiles = []
        # print(self.nprofiles)
        for i in range(self.nprofiles):
            self.profiles.append(mw.read_profile("{dir}/profile{index}.data".format(dir=self.folder, index=i+1)))

    @staticmethod
    def meanarr(arr, log):
        """
        Static function useful for getting the mean of an array can also log the numbers
        and then get the the mean of that


        If the array does contain negatives:
            - For positive values, it will use the log of the value
            - For negative values, it will use the negative of the log of the absolute magnitude of the value

        :param arr: Array to be passed through
        :param log: Take the log of the elements or not
        :return: output: The mean of the array with either log(elements) or elements
        """
        output = []
        if log:
            for i in arr:
                if i < 0:
                    output.append(-np.log(np.abs(i)))
                else:
                    output.append(np.log(i))
        else:
            output = arr

        return np.mean(output)

    @staticmethod
    def getindex(arr, point):
        """
        Static method used for finding the index of the entry provided or closest to it in the case of approximation
        by looking for smallest absolute difference between all the entries of the array and the point

        :param arr: Array to be searched
        :param point: Point to be searched for
        :return: index: Index of point or closest point in array
        """
        points = np.full_like(arr, point)
        temp = np.subtract(arr, points)
        # print(temp)
        # print(np.where(np.abs(temp) == np.min(np.abs(temp)))[0])
        index = int(np.where(np.abs(temp) == np.min(np.abs(temp)))[0])
        return index

    def cutdict(self, age):
        """
        Method for cutting the dictionary to relevant data

        :return: None
        """
        temp = self.datadict["star_age"]
        index = self.getindex(temp, age)
        keys = list(self.datadict.keys())[8:]
        for i in keys:
            self.datadict[i] = self.datadict[i][:index]

    def buildHR(self, show):
        """
        Function for building an H-R diagram based of class properties

        :param show: Boolean responsible for showing the plot
        :return: None
        """

        self.log_L = self.datadict["log_L"]
        self.log_Teff = self.datadict["log_Teff"]

        fig, ax = plt.subplots()
        ax.plot(self.log_Teff, self.log_L)
        ax.set_title("H-R Diagram")
        ax.set_xlabel("log(Temperature) [K]")
        ax.set_ylabel("log(Luminosity) [Lsun]")
        ax.invert_xaxis()
        if show:
            plt.show()
        try:
            fig.savefig(fname="{dir}/outputs/HRdiagram.png".format(dir=self.folder))
        except FileNotFoundError:
            os.mkdir("./{dir}/outputs".format(dir=self.folder))
            fig.savefig(fname="{dir}/outputs/HRdiagram.png".format(dir=self.folder))
        plt.close()

    def builddensity(self, show):
        """
        Build density and mass plot
        :return: None
        """

        for i in self.profiles:
            self.star_age.append(np.mean(i["star_age"], dtype=float))
            self.rho.append(np.mean(np.power(i["logRho"], 10)))
            self.mass.append(np.mean(i["mass"]))

        self.i1 = self.getindex(self.star_age, 73895630)
        self.i2 = self.getindex(self.star_age, 11800000000)
        self.i3 = self.getindex(self.star_age, 12477950000)
        self.i4 = self.getindex(self.star_age, 12477990000)

        ZAMS = [self.mass[:self.i1], self.rho[:self.i1]]
        TAMS = [self.mass[self.i1-1:self.i2], self.rho[self.i1-1:self.i2]]
        RGB = [self.mass[self.i2-1:self.i3], self.rho[self.i2-1:self.i3]]
        HeFlash = [self.mass[self.i3-1:self.i4+1], self.rho[self.i3-1:self.i4+1]]

        fig, ax = plt.subplots()
        ax.plot(ZAMS[0], ZAMS[1], label="ZAMS")
        ax.plot(TAMS[0], TAMS[1], label="TAMS")
        ax.plot(RGB[0], RGB[1], label="RGB")
        ax.plot(HeFlash[0], HeFlash[1], label="He Flash")
        ax.set_title("Rho vs Mass")
        ax.set_ylabel("Rho [g/cm^3]")
        ax.set_xlabel("Mass [M_sun]")
        ax.legend(loc="upper right")
        if show:
            plt.show()

        try:
            fig.savefig(fname="{dir}/outputs/rho_mass.png".format(dir=self.folder))
        except FileNotFoundError:
            os.mkdir("./{dir}/outputs".format(dir=self.folder))
            fig.savefig(fname="{dir}/outputs/rho_mass.png".format(dir=self.folder))
        plt.close()

    def buildtemp(self, show):
        """
        Build temp and mass plot
        :return: None
        """
        for i in self.profiles:
            self.temperature.append(np.mean(np.power(i["logT"], 10)))

        ZAMS = [self.mass[:self.i1], self.temperature[:self.i1]]
        TAMS = [self.mass[self.i1-1:self.i2], self.temperature[self.i1-1:self.i2]]
        RGB = [self.mass[self.i2-1:self.i3], self.temperature[self.i2-1:self.i3]]
        HeFlash = [self.mass[self.i3-1:self.i4], self.temperature[self.i3-1:self.i4]]

        fig, ax = plt.subplots()
        ax.plot(ZAMS[0], ZAMS[1], label="ZAMS")
        ax.plot(TAMS[0], TAMS[1], label="TAMS")
        ax.plot(RGB[0], RGB[1], label="RGB")
        ax.plot(HeFlash[0], HeFlash[1], label="He Flash")
        ax.set_title("Temperature vs Mass")
        ax.set_ylabel("Temperature [K]")
        ax.set_xlabel("Mass [M_sun]")
        ax.legend(loc="upper right")
        if show:
            plt.show()

        try:
            fig.savefig(fname="{dir}/outputs/logt_mass.png".format(dir=self.folder))
        except FileNotFoundError:
            os.mkdir("./{dir}/outputs".format(dir=self.folder))
            fig.savefig(fname="{dir}/outputs/logt_mass.png".format(dir=self.folder))
        plt.close()

    def buildpressure(self, show):
        """
        Building log pressure vs mass

        :param show: Show graph
        :return: None
        """
        for i in self.profiles:
            self.pressure.append(np.mean(i["pressure"]))

        ZAMS = [self.mass[:self.i1], self.pressure[:self.i1]]
        TAMS = [self.mass[self.i1-1:self.i2], self.pressure[self.i1-1:self.i2]]
        RGB = [self.mass[self.i2-1:self.i3], self.pressure[self.i2-1:self.i3]]
        HeFlash = [self.mass[self.i3-1:self.i4], self.pressure[self.i3-1:self.i4]]

        fig, ax = plt.subplots()
        ax.plot(ZAMS[0], ZAMS[1], label="ZAMS")
        ax.plot(TAMS[0], TAMS[1], label="TAMS")
        ax.plot(RGB[0], RGB[1], label="RGB")
        ax.plot(HeFlash[0], HeFlash[1], label="He Flash")
        ax.set_title("Pressure vs Mass")
        ax.set_ylabel("Pressure [dyn/cm^2]")
        ax.set_xlabel("Mass [M_sun]")
        ax.legend(loc="upper right")
        if show:
            plt.show()

        try:
            fig.savefig(fname="{dir}/outputs/pressure_mass.png".format(dir=self.folder))
        except FileNotFoundError:
            os.mkdir("./{dir}/outputs".format(dir=self.folder))
            fig.savefig(fname="{dir}/outputs/pressure_mass.png".format(dir=self.folder))
        plt.close()

    def buildpp(self, show):
        """
        Build pp reaction rate vs mass plots

        :return: None
        """
        for i in self.profiles:
            self.pp.append(np.mean(i["pp"]))

        ZAMS = [self.mass[:self.i1], self.pp[:self.i1]]
        TAMS = [self.mass[self.i1-1:self.i2], self.pp[self.i1-1:self.i2]]
        RGB = [self.mass[self.i2-1:self.i3], self.pp[self.i2-1:self.i3]]
        HeFlash = [self.mass[self.i3-1:self.i4], self.pp[self.i3-1:self.i4]]

        fig, ax = plt.subplots()
        ax.plot(ZAMS[0], ZAMS[1], label="ZAMS")
        ax.plot(TAMS[0], TAMS[1], label="TAMS")
        ax.plot(RGB[0], RGB[1], label="RGB")
        ax.plot(HeFlash[0], HeFlash[1], label="He Flash")
        ax.set_title("Proton-Proton reaction rate vs Mass")
        ax.set_ylabel("PP reaction rate [erg/s/g]")
        ax.set_xlabel("Mass [M_sun]")
        ax.legend(loc="upper right")
        if show:
            plt.show()

        try:
            fig.savefig(fname="{dir}/outputs/pp_mass.png".format(dir=self.folder))
        except FileNotFoundError:
            os.mkdir("./{dir}/outputs".format(dir=self.folder))
            fig.savefig(fname="{dir}/outputs/pp_mass.png".format(dir=self.folder))
        plt.close()

    def buildcno(self, show):
        """
        Build cno reaction rate vs mass plots

        :return: None
        """
        for i in self.profiles:
            self.cno.append(np.mean(i["cno"]))

        ZAMS = [self.mass[:self.i1], self.cno[:self.i1]]
        TAMS = [self.mass[self.i1-1:self.i2], self.cno[self.i1-1:self.i2]]
        RGB = [self.mass[self.i2-1:self.i3], self.cno[self.i2-1:self.i3]]
        HeFlash = [self.mass[self.i3-1:self.i4], self.cno[self.i3-1:self.i4]]

        fig, ax = plt.subplots()
        ax.plot(ZAMS[0], ZAMS[1], label="ZAMS")
        ax.plot(TAMS[0], TAMS[1], label="TAMS")
        ax.plot(RGB[0], RGB[1], label="RGB")
        ax.plot(HeFlash[0], HeFlash[1], label="He Flash")
        ax.set_title("CNO reaction rate vs Mass")
        ax.set_ylabel("CNO reaction rate [erg/s/g]")
        ax.set_xlabel("Mass [M_sun]")
        ax.legend(loc="upper right")
        if show:
            plt.show()

        try:
            fig.savefig(fname="{dir}/outputs/cno_mass.png".format(dir=self.folder))
        except FileNotFoundError:
            os.mkdir("./{dir}/outputs".format(dir=self.folder))
            fig.savefig(fname="{dir}/outputs/cno_mass.png".format(dir=self.folder))
        plt.close()

    def buildtriplealpha(self, show):
        """
        Build triple alpha reaction rate vs mass plots
        :return:
        """
        for i in self.profiles:
            self.tri_alfa.append(np.mean(i["tri_alfa"]))

        ZAMS = [self.mass[:self.i1], self.tri_alfa[:self.i1]]
        TAMS = [self.mass[self.i1-1:self.i2], self.tri_alfa[self.i1-1:self.i2]]
        RGB = [self.mass[self.i2-1:self.i3], self.tri_alfa[self.i2-1:self.i3]]
        HeFlash = [self.mass[self.i3-1:self.i4], self.tri_alfa[self.i3-1:self.i4]]

        fig, ax = plt.subplots()
        ax.plot(ZAMS[0], ZAMS[1], label="ZAMS")
        ax.plot(TAMS[0], TAMS[1], label="TAMS")
        ax.plot(RGB[0], RGB[1], label="RGB")
        ax.plot(HeFlash[0], HeFlash[1], label="He Flash")
        ax.set_title("Triple Alpha reaction rate vs Mass")
        ax.set_ylabel("Triple Alpha reaction rate")
        ax.set_xlabel("Mass [M_sun]")
        ax.legend(loc="upper right")
        if show:
            plt.show()

        try:
            fig.savefig(fname="{dir}/outputs/trialfa_mass.png".format(dir=self.folder))
        except FileNotFoundError:
            os.mkdir("./{dir}/outputs".format(dir=self.folder))
            fig.savefig(fname="{dir}/outputs/trialfa_mass.png".format(dir=self.folder))
        plt.close()

    def buildh1abundance(self, show):
        """
        Build Hydrogen vs mass plots
        :return: None
        """

        for i in self.profiles:
            self.h1.append(np.mean(i["star_mass_h1"]))

        ZAMS = [self.mass[:self.i1], self.h1[:self.i1]]
        TAMS = [self.mass[self.i1-1:self.i2], self.h1[self.i1-1:self.i2]]
        RGB = [self.mass[self.i2-1:self.i3], self.h1[self.i2-1:self.i3]]
        HeFlash = [self.mass[self.i3-1:self.i4], self.h1[self.i3-1:self.i4]]

        fig, ax = plt.subplots()
        ax.plot(ZAMS[0], ZAMS[1], label="ZAMS")
        ax.plot(TAMS[0], TAMS[1], label="TAMS")
        ax.plot(RGB[0], RGB[1], label="RGB")
        ax.plot(HeFlash[0], HeFlash[1], label="He Flash")
        ax.set_title("Hydrogen vs Mass")
        ax.set_ylabel("Hydrogen abundance [M_sun]")
        ax.set_xlabel("Mass [M_sun]")
        ax.legend(loc="upper right")
        if show:
            plt.show()

        try:
            fig.savefig(fname="{dir}/outputs/h1_mass.png".format(dir=self.folder))
        except FileNotFoundError:
            os.mkdir("./{dir}/outputs".format(dir=self.folder))
            fig.savefig(fname="{dir}/outputs/h1_mass.png".format(dir=self.folder))
        plt.close()

    def buildhe3abundance(self, show):
        """
        Build abundance plot for Helium
        :param show: Show plot
        :return: None
        """
        for i in self.profiles:
            self.he3.append(np.mean(i["star_mass_he3"]))

        ZAMS = [self.mass[:self.i1], self.he3[:self.i1]]
        TAMS = [self.mass[self.i1-1:self.i2], self.he3[self.i1-1:self.i2]]
        RGB = [self.mass[self.i2-1:self.i3], self.he3[self.i2-1:self.i3]]
        HeFlash = [self.mass[self.i3-1:self.i4], self.he3[self.i3-1:self.i4]]

        fig, ax = plt.subplots()
        ax.plot(ZAMS[0], ZAMS[1], label="ZAMS")
        ax.plot(TAMS[0], TAMS[1], label="TAMS")
        ax.plot(RGB[0], RGB[1], label="RGB")
        ax.plot(HeFlash[0], HeFlash[1], label="He Flash")
        ax.set_title("Helium-3 abundance vs Mass")
        ax.set_ylabel("Helium-3 abundance [M_sun]")
        ax.set_xlabel("Mass [M_sun]")
        ax.legend(loc="upper right")
        if show:
            plt.show()

        try:
            fig.savefig(fname="{dir}/outputs/he3_mass.png".format(dir=self.folder))
        except FileNotFoundError:
            os.mkdir("./{dir}/outputs".format(dir=self.folder))
            fig.savefig(fname="{dir}/outputs/he3_mass.png".format(dir=self.folder))
        plt.close()

    def buildhe4abundance(self, show):
        """
        Build abundance plot for Helium-4
        :param show: Show plot
        :return: None
        """
        for i in self.profiles:
            self.he4.append(np.mean(i["star_mass_he4"]))

        ZAMS = [self.mass[:self.i1], self.he4[:self.i1]]
        TAMS = [self.mass[self.i1-1:self.i2], self.he4[self.i1-1:self.i2]]
        RGB = [self.mass[self.i2-1:self.i3], self.he4[self.i2-1:self.i3]]
        HeFlash = [self.mass[self.i3-1:self.i4], self.he4[self.i3-1:self.i4]]

        fig, ax = plt.subplots()
        ax.plot(ZAMS[0], ZAMS[1], label="ZAMS")
        ax.plot(TAMS[0], TAMS[1], label="TAMS")
        ax.plot(RGB[0], RGB[1], label="RGB")
        ax.plot(HeFlash[0], HeFlash[1], label="He Flash")
        ax.set_title("Helium-4 abundance vs Mass")
        ax.set_ylabel("Helium-4 abundance [M_sun]")
        ax.set_xlabel("Mass [M_sun]")
        ax.legend(loc="upper right")
        if show:
            plt.show()

        try:
            fig.savefig(fname="{dir}/outputs/he4_mass.png".format(dir=self.folder))
        except FileNotFoundError:
            os.mkdir("./{dir}/outputs".format(dir=self.folder))
            fig.savefig(fname="{dir}/outputs/he4_mass.png".format(dir=self.folder))
        plt.close()

    def buildc12abundance(self, show):
        """
        Build abundance plot for Carbon-12
        :param show: Show plot
        :return: None
        """
        for i in self.profiles:
            self.c12.append(np.mean(i["star_mass_c12"]))

        ZAMS = [self.mass[:self.i1], self.c12[:self.i1]]
        TAMS = [self.mass[self.i1-1:self.i2], self.c12[self.i1-1:self.i2]]
        RGB = [self.mass[self.i2 - 1:self.i3], self.c12[self.i2 - 1:self.i3]]
        HeFlash = [self.mass[self.i3 - 1:self.i4], self.c12[self.i3 - 1:self.i4]]

        fig, ax = plt.subplots()
        ax.plot(ZAMS[0], ZAMS[1], label="ZAMS")
        ax.plot(TAMS[0], TAMS[1], label="TAMS")
        ax.plot(RGB[0], RGB[1], label="RGB")
        ax.plot(HeFlash[0], HeFlash[1], label="He Flash")
        ax.set_title("Carbon-12 abundance vs Mass")
        ax.set_ylabel("Carbon-12 abundance [M_sun]")
        ax.set_xlabel("Mass [M_sun]")
        ax.legend(loc="upper right")
        if show:
            plt.show()

        try:
            fig.savefig(fname="{dir}/outputs/c12_mass.png".format(dir=self.folder))
        except FileNotFoundError:
            os.mkdir("./{dir}/outputs".format(dir=self.folder))
            fig.savefig(fname="{dir}/outputs/c12_mass.png".format(dir=self.folder))
        plt.close()

    def buildo16abundance(self, show):
        """
        Build abundance plot for Carbon-12
        :param show: Show plot
        :return: None
        """
        for i in self.profiles:
            self.o16.append(np.mean(i["star_mass_o16"]))

        ZAMS = [self.mass[:self.i1], self.o16[:self.i1]]
        TAMS = [self.mass[self.i1-1:self.i2], self.o16[self.i1-1:self.i2]]
        RGB = [self.mass[self.i2 - 1:self.i3], self.o16[self.i2 - 1:self.i3]]
        HeFlash = [self.mass[self.i3 - 1:self.i4], self.o16[self.i3 - 1:self.i4]]

        fig, ax = plt.subplots()
        ax.plot(ZAMS[0], ZAMS[1], label="ZAMS")
        ax.plot(TAMS[0], TAMS[1], label="TAMS")
        ax.plot(RGB[0], RGB[1], label="RGB")
        ax.plot(HeFlash[0], HeFlash[1], label="He Flash")
        ax.set_title("Oxygen-16 abundance vs Mass")
        ax.set_ylabel("Oxygen-16 abundance [M_sun]")
        ax.set_xlabel("Mass [M_sun]")
        ax.legend(loc="upper right")
        if show:
            plt.show()

        try:
            fig.savefig(fname="{dir}/outputs/o16_mass.png".format(dir=self.folder))
        except FileNotFoundError:
            os.mkdir("./{dir}/outputs".format(dir=self.folder))
            fig.savefig(fname="{dir}/outputs/o16_mass.png".format(dir=self.folder))
        plt.close()

    def radvsconvec(self, show):
        """
        Determining which profiles are are convective or radiating and then plotting hr diagram with
        :param show:
        """

        for i in range(len(self.profiles)):
            if i < self.i4+1:
                self.grada.append(np.mean(self.profiles[i]["grada"]))
                self.gradr.append(np.mean(self.profiles[i]["gradr"]))
            else:
                break

        for i in range(len(self.gradr)):
            if self.gradr[i] < self.grada[i]:
                self.rad.append([self.mass[i], self.rho[i]])
            else:
                self.conv.append([self.mass[i], self.rho[i]])

        self.rad = np.stack(self.rad, axis=-1)
        self.conv = np.stack(self.conv, axis=-1)

        ZAMS = [self.mass[:self.i1+1], self.rho[:self.i1+1]]
        TAMS = [self.mass[self.i1:self.i2+1], self.rho[self.i1:self.i2+1]]
        RGB = [self.mass[self.i2:self.i3+1], self.rho[self.i2:self.i3+1]]
        HeFlash = [self.mass[self.i3:self.i4+1], self.rho[self.i3:self.i4+1]]

        fig2, ax2 = plt.subplots()
        ax2.plot(ZAMS[0], ZAMS[1], label="ZAMS")
        ax2.plot(TAMS[0], TAMS[1], label="TAMS")
        ax2.plot(RGB[0], RGB[1], label="RGB")
        ax2.plot(HeFlash[0], HeFlash[1], label="He Flash")
        plt.scatter(self.conv[0], self.conv[1], label="Convective", color='grey', marker='d')
        plt.scatter(self.rad[0], self.rad[1], label="Radiative", color='black', marker='s')
        ax2.set_title("Rho vs Mass \nwith Radiative and Convective sections overplotted")
        ax2.set_ylabel("Rho [g/cm^3]")
        ax2.set_xlabel("Mass [M_sun]")
        ax2.legend(loc="upper right")
        if show:
            plt.show()

        try:
            fig2.savefig(fname="{dir}/outputs/radconv.png".format(dir=self.folder))
        except FileNotFoundError:
            os.mkdir("./{dir}/outputs".format(dir=self.folder))
            fig2.savefig(fname="{dir}/outputs/radconv.png".format(dir=self.folder))
        plt.close()

    def definezams(self, age):
        """
        Define the parameters age range of the ZAMS
        :param age:
        :return:
        """
        self.ZAMS = dict.fromkeys(self.datadict.keys())
        temp = self.datadict["star_age"]
        index = self.getindex(temp, age)
        keys = list(self.datadict.keys())[8:]
        for i in keys:
            self.ZAMS[i] = self.datadict[i][:index]


def estimatelogTlogL(mass):
    """
    Create estimate of H-R diagram based on mas
    :param mass: Stellar mass in M_sun
    :return: 2-d array with list of logT and logL
    """
    massarr = mass
    logT = 0.25*np.log(np.power(massarr, 1.5)/(4*m.pi))
    logL = 3.5*np.log(massarr)
    return [logT, logL]


if __name__ == "__main__":
    allplots = False
    allhrs = False
    doestimates = True
    a = star("MESA1")
    a.definezams(age=7.389563 * 10 ** 7)

    # Plotting various quantities from one model
    if allplots:
        a.cutdict(age=7.389563*10**7)
        a.buildHR(show=True)
        a.buildHR(show=True)
        a.builddensity(show=True)
        a.buildtemp(show=True)
        a.buildpressure(show=True)
        a.buildpp(show=True)
        a.buildcno(show=True)
        a.buildtriplealpha(show=True)
        a.buildh1abundance(show=True)
        a.buildhe3abundance(show=True)
        a.buildhe4abundance(show=True)
        a.buildc12abundance(show=True)
        a.buildo16abundance(show=True)
        a.radvsconvec(show=True)

    # Plotting HR diagrams from models of various masses
    if allhrs:
        tenth = star("MESA28")
        ten = star("MESA27")
        twentyfive = star("MESA23")
        fourty = star("MESA10")
        fifty = star("MESA11")
        seventy = star("MESA12")
        eightyfive = star("MESA13")
        ninety = star("MESA14")
        hundo = star("MESA15")

        tenth.definezams(age=1.248102*10**5)
        ten.definezams(age=4.202688 * 10 ** 5)
        twentyfive.definezams(age=3.420263 * 10 ** 4)
        fourty.definezams(age=2.092475 * 10 ** 4)
        fifty.definezams(age=2.081969 * 10 ** 4)
        seventy.definezams(age=1.147982 * 10 ** 4)
        eightyfive.definezams(age=2.262014 * 10 ** 4)
        ninety.definezams(age=1.213162 * 10 ** 4)
        hundo.definezams(age=1.182317 * 10 ** 4)

        fig, ax = plt.subplots()
        ax.plot(tenth.ZAMS["log_Teff"], tenth.ZAMS["log_L"], label="0.1 M_sun")
        ax.plot(a.ZAMS["log_Teff"], a.ZAMS["log_L"], label="1 M_sun")
        ax.plot(ten.ZAMS["log_Teff"], ten.ZAMS["log_L"], label="10 M_sun")
        ax.plot(twentyfive.ZAMS["log_Teff"], twentyfive.ZAMS["log_L"], label="25 M_sun")
        ax.plot(fourty.ZAMS["log_Teff"], fourty.ZAMS["log_L"], label="40 M_sun")
        ax.plot(fifty.ZAMS["log_Teff"], fifty.ZAMS["log_L"], label="50 M_sun")
        ax.plot(seventy.ZAMS["log_Teff"], seventy.ZAMS["log_L"], label="70 M_sun")
        ax.plot(eightyfive.ZAMS["log_Teff"], eightyfive.ZAMS["log_L"], label="85 M_sun")
        ax.plot(ninety.ZAMS["log_Teff"], ninety.ZAMS["log_L"], label="90 M_sun")
        ax.plot(hundo.ZAMS["log_Teff"], hundo.ZAMS["log_L"], label="100 M_sun")
        ax.set_title("Multi H-R Diagram")
        ax.set_xlabel("log(Temperature) [K]")
        ax.set_ylabel("log(Luminosity) [Lsun]")
        ax.legend(loc="lower left")
        ax.invert_xaxis()
        plt.show()
        fig.savefig(fname="/multihr.png")
        plt.close()

        plt.close(fig)

    # Plotting estimated and actual H-R Diagrams of masses from 0.1-100 M_sun stars
    if doestimates:
        stellarmasses = [0.1*a.ZAMS["star_mass"], 1*a.ZAMS["star_mass"], 10*a.ZAMS["star_mass"], 25*a.ZAMS["star_mass"], 50*a.ZAMS["star_mass"], 70*a.ZAMS["star_mass"], 85*a.ZAMS["star_mass"], 90*a.ZAMS["star_mass"], 100*a.ZAMS["star_mass"]]
        stellarmasseslabel = [0.1, 1, 10, 25, 50, 70, 85, 90, 100]

        estimations = []
        for i in stellarmasses:
            print("Calculating estimations: i={}".format(i))
            estimations.append(estimatelogTlogL(np.array(i)))

        fig, ax = plt.subplots()
        for i in range(len(estimations)):
            print("Plotting estimations: i={}".format(i))
            ax.scatter(estimations[i][0], estimations[i][1], label="{} M_sun".format(stellarmasseslabel[i]))

        ax.set_title("H-R Diagram Estimations")
        ax.set_xlabel("log(Temperature) [K]")
        ax.set_ylabel("log(Luminosity) [Lsun]")
        ax.invert_xaxis()
        ax.legend(loc="lower left")
        plt.show()
        plt.close(fig)
