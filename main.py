import numpy as np
import astropy as ap
import mesa_web as mw
import matplotlib.pyplot as plt

class star:
    """
    Class for a folder of data
    """

    def __init__(self, folder):
        self.f = folder
        self.profiles = []
        self.L = []
        self.Teff = []

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

    def buildHR(self, show):
        """
        Function for building an H-R diagram based of class properties

        :param show: Boolean responsible for showing the plot
        :return: None
        """
        for i in range(306):
            self.profiles.append(mw.read_profile(filename="./{folder}/profile{index}.data".format(folder = self.f,index = i+1), as_table=False))

        for i in self.profiles:
            self.L.append(self.meanarr(arr = i["luminosity"], log=True))
            self.Teff.append(self.meanarr(arr = i["Teff"], log=False))

        self.Teff = np.log(self.Teff)

        fig, ax = plt.subplots()
        ax.plot(self.Teff, self.L)
        ax.set_title("H-R Diagram")
        ax.set_xlabel("Teff [K]")
        ax.set_ylabel("Luminosity [Lsun]")
        ax.invert_xaxis()
        if show:
            plt.show()
        plt.close()

# def arrmean(arr):
if __name__ == "__main__":
    a = star("MESA1")
    a.buildHR(show=True)

