import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rcParams['figure.figsize'] = (8.0, 6.0)


class Plot(object):
    def __init__(self, z, function, label, function_err=np.array([]), slope=None, intercept=None, xlabel="z",
                 ylabel=""):
        self.z = z
        self.function = function
        self.function_err = function_err
        self.label = label
        self.slope = slope
        self.intercept = intercept
        self.ylabel = ylabel
        self.xlabel = xlabel

    def ylabel(self):
        return self.ylabel

    def xlabel(self):
        return self.xlabel

    def plot_f(self, xlo=5, xhi=15, ylo=None, yhi=None):
        if len(self.function_err) > 0:
            plt.errorbar(self.z, self.function, yerr=self.function_err, linestyle='-', marker='o',
                         label=r"$ " + str(self.label) + "$")
        else:
            plt.plot(self.z, self.function, linestyle='-', marker='o',
                     label=r"$ " + str(self.label) + "$")
        plt.xlabel(r"$ " + str(self.xlabel) + "$", size=20)
        plt.xlim(xlo, xhi)
        plt.ylabel(r"$ " + str(self.ylabel) + "$", size=20)
        if ylo or yhi:
            plt.ylim(ylo, yhi)
        plt.legend(prop={'size': 10})