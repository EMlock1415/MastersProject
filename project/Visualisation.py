import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import csv
import scipy
import pandas as pd

print("Importing complete")
frequencies = [30,44,70,100,143,217,353,545,857]




asymmetry = pd.read_csv("../../dump/asymetrytest.csv")



print(asymmetry['m1'].to_list())
print(asymmetry['m2'].to_list())




plt.plot(asymmetry['ksloc'],'o')
#plt.plot(asymmetry['k2'],'o')
plt.yscale("log")
plt.savefig("../../dump/Graphs/asymmetryksloclog")
plt.yscale("linear")
plt.savefig("../../dump/Graphs/asymmetryksloc")
plt.close()




'''kstest = pd.read_csv("../../dump/KSTEST.csv")




simdata = pd.read_csv("../../dump/masked_stats0100.csv")
regdata = pd.read_csv("../../dump/REALMAPS.csv")

regpoints = np.arange(0,90,10)

plt.plot(simdata["kurtosis"],"o")
plt.plot(regpoints,regdata["kurt"],"o", color="red")
plt.yscale("linear")
plt.savefig("../../dump/Graphs/kurt.png")'''




'''alpha = 0.0005
n = 50000000
crit_val = np.sqrt(-np.log(alpha/2) * (n/(n**2)))
for i in range (0,9):

    tempdata = pd.read_csv("../../dump/KSTEST{}.csv".format(frequencies[i]))

    plt.plot(tempdata["Test Result"], "o")
    plt.hlines(y = crit_val,xmin=0,xmax=50,colors="red",linestyles="--")
    plt.savefig("../../dump/Graphs/KSTEST{}".format(frequencies[i]))

    plt.close()

    plt.plot(tempdata["p-value"], "o")
    plt.yscale("log")
    plt.savefig("../../dump/Graphs/KSTESTpvalue{}".format(frequencies[i]))

    plt.close()
    plt.yscale("linear")
    plt.plot(tempdata["location"], "o")
    plt.hlines(y = 0,xmin = 0, xmax=50,colors="red", linestyles="--")
    plt.savefig("../../dump/Graphs/KSTESTlocation{}".format(frequencies[i]))

    plt.close()
'''



