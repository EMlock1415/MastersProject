import numpy as np
import healpy as hp

np.random.seed(1)
NSIDE = 2048
pixels = hp.nside2npix(NSIDE)
fake = np.random.normal(loc=0,scale=1e-5,size=(pixels))
print(type(fake))