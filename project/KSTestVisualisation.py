import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import csv
import scipy
import pandas as pd

dist1 = np.random.normal(0,1,20)
dist2 = np.random.normal(0,1,20)

a = scipy.stats.kstest(dist1,dist2)
print(a)



plt.ecdf(dist1)
plt.ecdf(dist2)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.vlines(a.statistic_location, ymin=0,ymax=1,colors="black",linestyles="--")
plt.savefig("../../dump/Graphs/aaaaaaaaaaaaaaaa")
