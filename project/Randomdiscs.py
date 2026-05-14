# Disc Histogram hopefully

#%% 
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import time
import scipy

import logging
log = logging.getLogger("healpy")

print("Importing worked!")
# %%
NSIDE = 32
NPIX = hp.nside2npix(NSIDE)
print("Approximate resolution at NSIDE {} is {:.3} deg. {} total pixels".format(NSIDE, hp.nside2resol(NSIDE, arcmin=True)/60, NPIX))

#%%

wmap_map_I = hp.read_map("wmap_band_iqumap_r9_7yr_W_v4.fits")
O_NSIDE = hp.get_nside(wmap_map_I)
print(O_NSIDE)
O_NPIX = hp.nside2npix(O_NSIDE)
res_map = hp.ud_grade(wmap_map_I, nside_out=NSIDE)
hp.mollview(wmap_map_I, unit="mK", norm="hist")
hp.mollview(res_map, unit="mK", norm="hist")
hp.graticule()

#%%

mask = hp.read_map("wmap_temperature_analysis_mask_r9_7yr_v4.fits")

wmap_map_I_masked = hp.ma(wmap_map_I)
wmap_map_I_masked.mask = np.logical_not(mask)
hp.mollview(wmap_map_I_masked.filled())

res_mask = hp.ud_grade(mask, nside_out=NSIDE)

res_map_masked = hp.ma(res_map)
res_map_masked.mask = np.logical_not(res_mask)
hp.mollview(res_map_masked.filled())

#%%
plt.hist(wmap_map_I_masked, bins=100, label="Full Resolution Masked Historgram")
plt.title("Full Resolution Masked Historgram")
#%%
plt.hist(res_map_masked, bins=100, label="NSIDE = 32 Masked Histogram")
plt.title("NSIDE=32 Masked Historgram")
# %%

'''NSIDE = 32; making discs

Don't do this one'''

wmap_map_I = hp.read_map("wmap_band_iqumap_r9_7yr_W_v4.fits")
O_NSIDE = hp.get_nside(wmap_map_I)
O_NPIX = hp.nside2npix(O_NSIDE)
res_map = hp.ud_grade(wmap_map_I, nside_out=NSIDE)

res_map_used = res_map
angles = np.zeros(12)
placements = []
discs = []

for i in range(0,12):
    angles[i] = i * np.pi / 12
    placements.append(hp.ang2vec(angles[i],0))
    discs.append(hp.query_disc(nside=NSIDE, vec=placements[i],radius=np.radians(5)))
    res_map_used[discs[i]] = res_map_used.max()
print(angles/np.pi)



hp.mollview(res_map_used, title="Hopefully some discs will be yellow", norm="hist")

# %%

'''NSIDE = 512 (original); making discs

Or this one'''

wmap_map_I = hp.read_map("wmap_band_iqumap_r9_7yr_W_v4.fits")
O_NSIDE = hp.get_nside(wmap_map_I)
O_NPIX = hp.nside2npix(O_NSIDE)
res_map = hp.ud_grade(wmap_map_I, nside_out=NSIDE)

full_map_used = wmap_map_I
angles = np.zeros(12)
placements = []
discs = []


angles = [6 * np.pi / 12, 0 * np.pi / 6]
placements.append(hp.ang2vec(angles[0],angles[1]))
discs.append(hp.query_disc(nside=O_NSIDE, vec=placements[0],radius=np.radians(5)))
full_map_used[discs[0]] = full_map_used.max()

hp.mollview(full_map_used, title="Better Resolution?", norm="hist")

# %%

'''NSIDE = 512 (original res); making discs (and doing statistics???)


Do this one! Not the others! 
'''

wmap_map_I = hp.read_map("wmap_band_iqumap_r9_7yr_W_v4.fits")
O_NSIDE = hp.get_nside(wmap_map_I)
O_NPIX = hp.nside2npix(O_NSIDE)
res_map = hp.ud_grade(wmap_map_I, nside_out=NSIDE)

full_map_used = wmap_map_I_masked
v_angles = np.zeros(12)
h_angles = np.zeros(12)
placements = []
discs = []

for i in range(0,12):
    v_angles[i] = np.pi * i / 12
    h_angles[i] = np.pi * i / 6
for i in range(0,12):
    for j in range(0,12):
        placements.append(hp.ang2vec(v_angles[i], h_angles[j]))
        discs.append(hp.query_disc(nside=O_NSIDE, vec=placements[-1], radius=np.radians(5)))
        

# hp.mollview(full_map_used, title="Better Resolution?", norm="hist")


# %%

'''makes a plot of a disc'''
the_disc = discs[40]
#print(the_disc)
#print(full_map_used[the_disc])
plt.hist(full_map_used[the_disc],bins=100)

# %%

'''makes plots of all of the discs'''


time_start = time.time()
checkpoints=[]
plots = []
for i in range(len(discs)):
    plots.append(wmap_map_I_masked[discs[i]])
checkpoints.append(time.time())
print("List for plots made in {} seconds".format(checkpoints[0]-time_start))
'''disc_n = 60
plt.hist(plots[disc_n],bins=100)
plt.title("Histogram of disc no. {}".format(disc_n))'''
plt.rcParams.update({'font.size':1})
fig, axs =plt.subplots(12,12)
for i in range(12):
    for j in range(12):
        axs[j,i].hist(plots[j+12*i-1],bins=100)
        axs[j,i].set_title("Histogram of disc no. {}".format(j+12*i))
checkpoints.append(time.time())
print("Making plots took another {}".format(checkpoints[1]-checkpoints[0]))
print("Total time: {}".format(checkpoints[1]-time_start))




# %%

'''Skewness plot'''

disc_skew = []
new_plots = []
for i in range(11,len(plots)):
    new_plots.append(plots[i])

for i in range(len(new_plots)):
    disc_skew.append(scipy.stats.skew(new_plots[i], nan_policy='omit'))
    print(disc_skew[i])

plt.rcParams.update({'font.size':10})
plt.hist(disc_skew, bins=72)

#print(scipy.stats.skewtest(new_plots[50]))
print(disc_skew)

skew_mean = np.nanmean(disc_skew)
skew_var = np.nanvar(disc_skew)
skew_std = np.nanstd(disc_skew)
print("Mean of skew: {}".format(skew_mean))
print("Variation of the skew: {}".format(skew_var))
print("Standard deviation of the skew: {} = {}".format(np.sqrt(skew_var), skew_std))

# %%

'''Kurtosis plot?'''

disc_kurt = []
for i in range(11,len(plots)):
    new_plots.append(plots[i])

for i in range(len(new_plots)):
    disc_kurt.append(scipy.stats.kurtosis(new_plots[i], nan_policy='omit'))
    print(disc_kurt[i])

plt.hist(disc_kurt, bins=72)

kurt_mean = np.nanmean(disc_kurt)
kurt_var = np.nanvar(disc_kurt)
kurt_std = np.nanstd(disc_kurt)

print("Mean of kurtosis: {}".format(kurt_mean))
print("Variation of kurtotis: {}".format(kurt_var))
print("Standard deviation of kurtosis: {}".format(kurt_std))

# %%

'''Another statistic????'''



