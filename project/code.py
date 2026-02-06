#%%
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

import logging
log = logging.getLogger("healpy")

print("Importing worked!")
#%%



NSIDE = 32
NPIX = hp.nside2npix(NSIDE)
print("Approximate resolution at NSIDE {} is {:.3} deg. {} total pixels".format(NSIDE, hp.nside2resol(NSIDE, arcmin=True)/60, NPIX))
#%%
print("new cell")
# %%
m = np.arange(NPIX)
hp.mollview(m, title=("Mollview image RING"))
hp.graticule

# %%
'''Angle 1: theta runs from 0 to pi, north to south. 
Angle 2: phi runs from 0-2pi eastward from centre'''
vec = hp.ang2vec(np.pi /2, 3*np.pi /4)
print(vec)
# %%

'''making a little disk on the thingy'''

ipix_disc = hp.query_disc(nside=NSIDE, vec=vec, radius=np.radians(10)) # returns pixels defined by a vector (point) and radius around
m[ipix_disc] = m.max()
hp.mollview(m, title=("Mollview image RING"))

# %%
theta, phi = np.degrees(hp.pix2ang(nside=NSIDE, ipix=[0,1,2,3,4]))
print(theta)
print(phi)
# %%

hp.mollview(m, nest=True, title="Mollview image NESTED")

# %%

wmap_map_I = hp.read_map("wmap_band_iqumap_r9_7yr_W_v4.fits")
changed_map = hp.ud_grade(wmap_map_I,nside_out=NSIDE)
hp.get_nside(changed_map)



# %%

# Making some plots of WMAP 7yr data

hp.mollview(wmap_map_I,title=("temp map"), coord=["G","E"], unit="mK", norm="hist")
hp.mollview(wmap_map_I,title=("temp map"), unit="mK", norm="hist")
hp.graticule()
# %%
# gnomview

hp.gnomview(wmap_map_I, title=("gnom map"),rot=[0,0.3], unit="mK", format="%.2g")

# %%

mask = hp.read_map("wmap_temperature_analysis_mask_r9_7yr_v4.fits")
wmap_map_I_masked = hp.ma(wmap_map_I)
wmap_map_I_masked.mask = np.logical_not(mask)


# %%

hp.mollview(wmap_map_I_masked.filled())
# %%

# masked image histogram

plt.hist(wmap_map_I_masked,bins=1000)

# %%
plt.hist(wmap_map_I, bins=1000) # unmasked histogram

# %%
angles = np.zeros(12)
placements = []
discs = []

for i in range(0,12):
    angles[i] = i * np.pi / 12
    placements.append(hp.ang2vec(angles[i],0))
    discs.append(hp.query_disc(nside=NSIDE, vec=placements[i],radius=np.radians(5)))
    m[discs[i]] = m.max()
print(angles/np.pi)



hp.mollview(m, title="Hopefully some discs will be yellow")


# %%
