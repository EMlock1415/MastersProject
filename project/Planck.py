#%%
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import time
import scipy

import logging
log = logging.getLogger("healpy")

print("Importing worked!")

#%%

planckMaps = {
    30 : hp.read_map("../PlanckMaps/LFI_SkyMap_030_1024_R3.00_full.fits"),
    44 : hp.read_map("../PlanckMaps/LFI_SkyMap_044_1024_R3.00_full.fits"),
    70 : hp.read_map("../PlanckMaps/LFI_SkyMap_070_1024_R3.00_full.fits"),
    100 : hp.read_map("../PlanckMaps/HFI_SkyMap_100_2048_R3.01_full.fits"),
    143 : hp.read_map("../PlanckMaps/HFI_SkyMap_143_2048_R3.01_full.fits"),
    217 : hp.read_map("../PlanckMaps/HFI_SkyMap_217_2048_R3.01_full.fits"),
    353 : hp.read_map("../PlanckMaps/HFI_SkyMap_353_2048_R3.01_full.fits"),
    545 : hp.read_map("../PlanckMaps/HFI_SkyMap_545_2048_R3.01_full.fits"),
    857 : hp.read_map("../PlanckMaps/HFI_SkyMap_857_2048_R3.01_full.fits")
}

#%%

planckMaps[100][45581992] = np.mean(planckMaps[100])
planckMapFrequencies = [30,44,70,100,143,217,353,545,857]


#%%
'''
READING MAP (not needed)
'''
planckMap = planckMaps[100]
#planckMap = hp.ud_grade(hp.read_map("../PlanckMaps/HFI_SkyMap_100_2048_R3.01_full.fits"),2048)
plankMask = hp.ud_grade(hp.read_map("../PlanckMaps/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits"), 2048)
planckMapMasked = hp.ma(planckMap)
planckMapMasked.mask = np.logical_not(plankMask)
wmapScalePlanckMap = hp.ud_grade(planckMap, nside_out=512)
wmapMap = hp.read_map("../WMAPMaps/wmap_band_iqumap_r9_7yr_W_v4.fits")
wmapMask = hp.read_map("../WMAPMaps/wmap_temperature_analysis_mask_r9_7yr_v4.fits")
wmapMapMasked = hp.ma(wmapMap)
wmapMapMasked.mask = np.logical_not(wmapMask)
print(hp.get_nside(planckMap))
print(np.min(planckMapMasked))
#%%
'''mAKING ALL THE PICTUREs (not needed)'''

for i in range(len(planckMapFrequencies)):
    hp.mollview(planckMaps[planckMapFrequencies[i]],unit="K", norm="hist", title="Planck Map {}GHz".format(planckMapFrequencies[i]))
#%%
'''MAKING MASKED MAPS AND PICTURES (needed again)'''
maskedMaps = []
lowresplanckMask = hp.ud_grade(hp.read_map("../PlanckMaps/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits"), 1024)

for i in range(3):
    maskedMaps.append(hp.ma(planckMaps[planckMapFrequencies[i]]))
    maskedMaps[i].mask = np.logical_not(lowresplanckMask)

planckMask = hp.read_map("../PlanckMaps/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits")
for i in range(6):
    maskedMaps.append(hp.ma(planckMaps[planckMapFrequencies[i+3]]))
    maskedMaps[i+3].mask = np.logical_not(planckMask)

for i in range(len(planckMapFrequencies)):
    hp.mollview(maskedMaps[i],unit="K", norm="hist", title="Planck Map {}GHz".format(planckMapFrequencies[i]))

#%%
'''MAKING HISTOGRAMS'''
for i in range(len(planckMapFrequencies)):
    plt.hist(maskedMaps[i],bins=100)
    plt.title("Planck Histogram: {}GHz".format(planckMapFrequencies[i]))
    plt.show()
# %%
'''FINDING PROBLEMS'''
count = 0
problems = []
for i in range(len(planckMap)):
    if (i%1000 == 0):
        print("Pixel {}".format(len(planckMap) - i-1))
    if planckMap[len(planckMap) - i-1] < -10000:
        print(planckMap[len(planckMap) - i-1])
        print(len(planckMap) - i-1)
        count += 1
        problems.append(len(planckMap) - i - 1)
        
        
print(count)
print(problems)
# %%
'''
STATISTICS BELOW

'''
#%%
'''SKEW'''
for i in range(len(planckMapFrequencies)):
    a = scipy.stats.skew(maskedMaps[i])
    print("{}GHz frequency band skewness: {}".format(planckMapFrequencies[i], a))

# %%
'''KURTOSIS'''
for i in range(len(planckMapFrequencies)):
    a = scipy.stats.kurtosis(maskedMaps[i])
    print("{}GHz frequency band kurtosis: {}".format(planckMapFrequencies[i], a))

# %%
'''MEAN'''
for i in range(len(planckMapFrequencies)):
    a = np.mean(maskedMaps[i])
    print("{}GHz frequency band mean: {}".format(planckMapFrequencies[i], a))
#%%
'''VARIANCE'''
for i in range(len(planckMapFrequencies)):
    a = np.var(maskedMaps[i])
    print("{}GHz frequency band variance: {}".format(planckMapFrequencies[i], a))
#%%
'''D'AGOSTINO???'''
for i in range(len(planckMapFrequencies)):
    a = scipy.stats.normaltest(maskedMaps[i])
    print("{}GHz frequency band D'Agostino K^2 statistic: {}".format(planckMapFrequencies[i], a))

# %%

'''
DISC CREATION
'''
phi =  np.pi / 5 # north to south
theta =  np.pi/2
vector = hp.ang2vec(phi,theta)
NSIDE = 2048
pixels = hp.nside2npix(NSIDE)
disc = hp.query_disc(nside=NSIDE,vec=vector,radius=np.radians(2))
m = np.arange(pixels)
m[disc] = m.max()
hp.mollview(m)

# %%

a = scipy.stats.normaltest(maskedMaps[2][disc])
print("GHz frequency band D'Agostino K^2 statistic (on small piece of sky): {}".format(a))
a = scipy.stats.normaltest(maskedMaps[3][disc])
print("GHz frequency band D'Agostino K^2 statistic (on small piece of sky): {}".format(a))
#%%
a = scipy.stats.normaltest(maskedMaps[4][disc])
print("100GHz frequency band D'Agostino K^2 statistic (on small piece of sky): {}".format(a))
a = scipy.stats.normaltest(maskedMaps[5][disc])
print("143GHz frequency band D'Agostino K^2 statistic (on small piece of sky): {}".format(a))
a = scipy.stats.normaltest(maskedMaps[6][disc])
print("GHz frequency band D'Agostino K^2 statistic (on small piece of sky): {}".format(a))
# %%

mean = np.mean(maskedMaps[4][disc])
#print(mean)
n = len(maskedMaps[4][disc])
orderedMap = np.sort(maskedMaps[4][disc])
T = 0
S_squared = 0
for i in range(n):
    S_squared += (1/n) * ((orderedMap[i]) - mean)**2
for i in range(n):
    T += (i - 0.5 * (n + 1))*orderedMap[i]
print(T/((n^2) * np.sqrt(S_squared)))


# %%

frequency = input("Frequency")
usedMap = hp.read_map("../PlanckMaps/HFI_SkyMap_{}_2048_R3.01_full.fits".format(frequency))
hp.mollview(usedMap, unit="K", norm="hist", title="Map")



# %%
