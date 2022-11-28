import numpy as np
import astropy as ap
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

    def cutdict(self):
        temp = self.datadict["star_age"]
        index1 = self.getindex(temp, 1.2477973*10**10)
        index2 = self.getindex(temp, 7.389563*10**7)
        keys = list(self.datadict.keys())[8:]
        for i in keys:
            self.datadict[i] = self.datadict[i][:index2]

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
        ax.set_xlabel("Teff [K]")
        ax.set_ylabel("Luminosity [Lsun]")
        ax.invert_xaxis()
        if show:
            plt.show()
        try:
            fig.savefig(fname="{dir}/outputs/HRdiagram.png".format(dir=self.folder))
        except FileNotFoundError:
            os.mkdir("./{dir}/outputs".format(dir=self.folder))
            fig.savefig(fname="{dir}/outputs/HRdiagram.png".format(dir=self.folder))
        plt.close()

    def definesequences(self, arr, point1, point2, point3, point4):
        self.i1 = self.getindex(arr, point1)
        self.i2 = self.getindex(arr, point2)
        self.i3 = self.getindex(arr, point3)
        self.i4 = self.getindex(arr, point4)

    def builddensity(self, show):
        """
        Build density and mass plot
        :return:
        """
        self.star_age = []
        self.logRho = []
        self.mass = []

        # print(self.profiles[0]["star_age"])
        # print(self.profiles[0]["logRho"])
        # print(self.profiles[0]["mass"])

        for i in self.profiles:
            self.star_age.append(np.mean(i["star_age"], dtype=float))
            self.logRho.append(np.mean(i["logRho"]))
            self.mass.append(np.mean(i["mass"]))

        self.i1 = self.getindex(self.star_age, 73895630)

        self.i2 = self.getindex(self.star_age, 11800000000)

        self.i3 = self.getindex(self.star_age, 12477950000)

        self.i4 = self.getindex(self.star_age, 12477980000)

        ZAMS = [self.mass[:self.i1], self.logRho[:self.i1]]
        TAMS = [self.mass[self.i1:self.i2], self.logRho[self.i1:self.i2]]
        RGB = [self.mass[self.i2-1:self.i3], self.logRho[self.i2-1:self.i3]]
        HeFlash = [self.mass[self.i3-1:self.i4], self.logRho[self.i3-1:self.i4]]

        fig, ax = plt.subplots()
        ax.plot(ZAMS[0], ZAMS[1])
        ax.plot(TAMS[0], TAMS[1])
        ax.plot(RGB[0], RGB[1])
        ax.plot(HeFlash[0], HeFlash[1])
        ax.set_title("Log(Density) vs Mass")
        ax.set_ylabel("log(Rho)")
        ax.set_xlabel("Mass")
        if show:
            plt.show()
        plt.close()

        fig, ax = plt.subplots()
        ax.plot(ZAMS[0], ZAMS[1])
        ax.set_title("Log(Density) vs Mass: ZAMS")
        ax.set_ylabel("log(Rho)")
        ax.set_xlabel("Mass")
        if not show:
            plt.show()
        plt.close()

        fig, ax = plt.subplots()
        ax.plot(TAMS[0], TAMS[1])
        ax.set_title("Log(Density) vs Mass: TAMS")
        ax.set_ylabel("log(Rho)")
        ax.set_xlabel("Mass")
        if not show:
            plt.show()
        plt.close()

        fig, ax = plt.subplots()
        ax.plot(RGB[0], RGB[1])
        ax.set_title("Log(Density) vs Mass: RGB")
        ax.set_ylabel("log(Rho)")
        ax.set_xlabel("Mass")
        if not show:
            plt.show()
        plt.close()

        fig, ax = plt.subplots()
        ax.plot(HeFlash[0], HeFlash[1])
        ax.set_title("Log(Density) vs Mass: He Flash")
        ax.set_ylabel("log(Rho)")
        ax.set_xlabel("Mass")
        if not show:
            plt.show()
        plt.close()


        # try:
        #     fig.savefig(fname="{dir}/outputs/logrho_mass.png".format(dir=self.folder))
        # except FileNotFoundError:
        #     os.mkdir("./{dir}/outputs".format(dir=self.folder))
        #     fig.savefig(fname="{dir}/outputs/logrho_mass.png".format(dir=self.folder))
        # plt.close()

    def buildtemp(self):
        """
        Build temp and mass plot
        :return:
        """

    def buildburn(self):
        """
        Build burn rate vs mass plots
        :return:
        """

    def buildabundance(self):
        """
        Build abundance vs mass plots
        :return:
        """

# def arrmean(arr):
if __name__ == "__main__":
    a = star("MESA1")
    a.cutdict()
    # a.buildHR(show=False)
    a.builddensity(show=True)
