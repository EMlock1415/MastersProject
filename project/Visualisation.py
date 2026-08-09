import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import csv
import scipy
import pandas as pd

print("Importing complete")
frequencies = [30,44,70,100,143,217,353,545,857]


font_size = 18

asymmetry = pd.read_csv("../../dump/asymetrytest.csv")


plt.locator_params(axis="y", nbins=3)
print(asymmetry['m1'].to_list())
print(asymmetry['m2'].to_list())




plt.plot(asymmetry['dagvalue1'],'o')
plt.plot(asymmetry['dagvalue2'],'o')
plt.yscale("log")
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.savefig("../../dump/Graphs/asymmetrydaglog")
plt.yscale("linear")
plt.savefig("../../dump/Graphs/asymmetrydag")
plt.close()


plt.locator_params(axis="y", nbins=3)
plt.plot(asymmetry['m1'],'o')
plt.plot(asymmetry['m1'],'o')
plt.yscale("log")
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.savefig("../../dump/Graphs/asymmetrymeanlog")
plt.yscale("linear")
plt.savefig("../../dump/Graphs/asymmetrymean")
plt.close()

plt.locator_params(axis="y", nbins=3)
plt.plot(asymmetry['s1'],'o')
plt.plot(asymmetry['s2'],'o')
plt.yscale("symlog")
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.savefig("../../dump/Graphs/asymmetryskewlog")
plt.yscale("linear")
plt.savefig("../../dump/Graphs/asymmetryskew")
plt.close()

plt.locator_params(axis="y", nbins=3)
plt.plot(asymmetry['k1'],'o')
plt.plot(asymmetry['k2'],'o')
plt.yscale("symlog")
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.savefig("../../dump/Graphs/asymmetrykurtlog")
plt.yscale("linear")
plt.savefig("../../dump/Graphs/asymmetrykurt")
plt.close()

plt.locator_params(axis="y", nbins=3)
plt.plot(asymmetry['ksvalue'],'o')
plt.yscale("log")
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.savefig("../../dump/Graphs/asymmetryksvaluelog")
plt.yscale("linear")
plt.savefig("../../dump/Graphs/asymmetryksvalue")
plt.close()

plt.locator_params(axis="y", nbins=3)
plt.plot(asymmetry['ksp'],'o')
plt.yscale("log")
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.savefig("../../dump/Graphs/asymmetryksplog")
plt.yscale("linear")
plt.savefig("../../dump/Graphs/asymmetryksp")
plt.close()

plt.locator_params(axis="y", nbins=3)
plt.plot(asymmetry['ksloc'],'o')
plt.yscale("log")
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.savefig("../../dump/Graphs/asymmetryksloclog")
plt.yscale("linear")
plt.savefig("../../dump/Graphs/asymmetryksloc")
plt.close()




kstest = pd.read_csv("../../dump/KSTEST.csv")




simdata = pd.read_csv("../../dump/masked_stats0100.csv")
regdata = pd.read_csv("../../dump/REALMAPS.csv")



regpoints = np.arange(0,90,10)

plt.locator_params(axis="y", nbins=3)
plt.plot(simdata["kurtosis"],"o")
plt.plot(regpoints,regdata["kurt"],"o", color="red")
plt.yscale("linear")
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.savefig("../../dump/Graphs/kurt.png")
plt.yscale("log")
plt.savefig("../../dump/Graphs/kurtlog.png")
plt.close()

plt.locator_params(axis="y", nbins=3)
plt.plot(simdata["skew"],"o")
plt.plot(regpoints,regdata["skew"],"o", color="red")
plt.yscale("linear")
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.savefig("../../dump/Graphs/skew.png")
plt.yscale("log")
plt.savefig("../../dump/Graphs/skewlog.png")
plt.close()

plt.locator_params(axis="y", nbins=3)
plt.plot(simdata["dagostino"],"o")
plt.plot(regpoints,regdata["dagostino"],"o", color="red")
plt.yscale("linear")
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.savefig("../../dump/Graphs/dagostino.png")
plt.yscale("log")
plt.savefig("../../dump/Graphs/dagostinolog.png")
plt.close()






alpha = 0.0005
n = 50000000
crit_val = np.sqrt(-np.log(alpha/2) * (n/(n**2)))
for i in range (0,9):

    tempdata = pd.read_csv("../../dump/KSTEST{}.csv".format(frequencies[i]))

    plt.locator_params(axis="y", nbins=3)
    plt.plot(tempdata["Test Result"], "o")
    plt.hlines(y = crit_val,xmin=0,xmax=50,colors="red",linestyles="--")
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.yscale("linear")
    plt.savefig("../../dump/Graphs/KSTEST{}".format(frequencies[i]))
    plt.yscale("log")
    plt.savefig("../../dump/Graphs/KSTESTLOG{}".format(frequencies[i]))


    plt.close()

    plt.locator_params(axis="y", nbins=3)
    plt.plot(tempdata["p-value"], "o")
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.yscale("log")
    plt.savefig("../../dump/Graphs/KSTESTpvalue{}".format(frequencies[i]))
    #plt.yscale("linear")
    plt.savefig("../../dump/Graphs/KSTESTLOGpvalue{}".format(frequencies[i]))


    plt.close()
    
    plt.locator_params(axis="y", nbins=3)
    plt.plot(tempdata["location"], "o")
    plt.hlines(y = 0,xmin = 0, xmax=50,colors="red", linestyles="--")
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.yscale("linear")
    plt.savefig("../../dump/Graphs/KSTESTlocation{}".format(frequencies[i]))
    plt.yscale("symlog")
    plt.savefig("../../dump/Graphs/KSTESTLOGlocation{}".format(frequencies[i]))

    plt.close()



'''plt.hist(simdata["dagostino"], bins=10,label="a")
plt.hist(regdata["dagostino"], bins=10, label="b")
plt.legend()
plt.savefig("../../dump/Graphs/daghist.png")'''
