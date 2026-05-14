import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import time
import scipy
import random

import logging
log = logging.getLogger("healpy")

print("Importing worked!")


planckMaps = {
    30 : "../../PlanckMaps/LFI_SkyMap_030_1024_R3.00_full.fits",
    44 : "../../PlanckMaps/LFI_SkyMap_044_1024_R3.00_full.fits",
    70 : "../../PlanckMaps/LFI_SkyMap_070_1024_R3.00_full.fits",
    100 : "../../PlanckMaps/HFI_SkyMap_100_2048_R3.01_full.fits",
    143 : "../../PlanckMaps/HFI_SkyMap_143_2048_R3.01_full.fits",
    217 : "../../PlanckMaps/HFI_SkyMap_217_2048_R3.01_full.fits",
    353 : "../../PlanckMaps/HFI_SkyMap_353_2048_R3.01_full.fits",
    545 : "../../PlanckMaps/HFI_SkyMap_545_2048_R3.01_full.fits",
    857 : "../../PlanckMaps/HFI_SkyMap_857_2048_R3.01_full.fits"
}

inputting = True
map=[]
masked_map = []
disc = []
mask = hp.read_map("../../PlanckMaps/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits")




def function():
    global map
    frequency = input("Frequency ")
    try:
        print("Reading map with frequency {}".format(frequency))
        map = hp.read_map(planckMaps[int(frequency)])
    except:
        print("Invalid frequency probably")
    hp.mollview(map, norm="hist", title="current map")
    plt.savefig("../../dump/map1.png")


def load_map():
    global map
    frequency = input("Frequency ")
    try:
        print("Reading map with frequency {}".format(frequency))
        map = hp.ud_grade(hp.read_map(planckMaps[int(frequency)]),2048)
    except:
        print("Invalid frequency probably")

def mask_map():
    global masked_map
    masked_map = hp.ma(map)
    masked_map.mask = np.logical_not(mask)


def save_map():
    maps = {
        "map" : map,
        "masked_map" : masked_map
    }
    saved_map = maps[input("map/masked_map? ")]
    print(len(saved_map))
    location = input("Save file to where: ")
    hp.mollview(saved_map, norm="hist", title="current map")
    plt.savefig("../../dump/{}.png".format(location))
    plt.close()
    #plt.savefig("dump/mapp.png")

def histograms():
    maps = {
        "map" : map,
        "masked_map" : masked_map,
        "disc" : masked_map[disc]
    }
    working_map = maps[input("what are we working with here???? ")]
    x = plt.hist(working_map,bins=50)
    plt.title(input("Title? "))
    plt.savefig("../../dump/{}.png".format(input("File name ")))
    plt.close()

def skew():
    maps = {
        "map" : map,
        "masked_map" : masked_map,
        "disc" : masked_map[disc]
    }
    x = input("what are we working with here???? ")
    working_map = maps[x]
    a = scipy.stats.skew(working_map)
    print("Skew of {} : {}".format(x,a))


def kurt():
    maps = {
        "map" : map,
        "masked_map" : masked_map,
        "disc" : masked_map[disc]
    }
    x = input("what are we working with here???? ")
    working_map = maps[x]
    a = scipy.stats.kurtosis(working_map)
    print("Kurtosis of {} : {}".format(x,a))

def mean():
    maps = {
        "map" : map,
        "masked_map" : masked_map,
        "disc" : masked_map[disc]
    }
    x = input("what are we working with here???? ")
    working_map = maps[x]
    a = np.mean(working_map)
    print("Mean of {} : {}".format(x,a))

def var():
    maps = {
        "map" : map,
        "masked_map" : masked_map,
        "disc" : masked_map[disc]
    }
    x = input("what are we working with here???? ")
    working_map = maps[x]
    a = np.var(working_map)
    print("Variance of {} : {}".format(x,a))

def create_disc():
    phi =  float(input("Phi (north to south, [0,pi])" ) ) # north to south
    theta =  float(input("Theta (west to east, [0,2pi])" ) )
    vector = hp.ang2vec(phi,theta)
    NSIDE = 2048
    pixels = hp.nside2npix(NSIDE)
    d = hp.query_disc(nside=NSIDE,vec=vector,radius=np.radians(2))
    m = np.arange(pixels)
    m[d] = m.max()
    hp.mollview(m)
    plt.savefig("../../dump/disc.png")
    plt.close()

    global disc
    disc = d

def simulation():
    seed = input("Seed")
    np.random.seed(int(seed))
    NSIDE = 2048
    pixels = hp.nside2npix(NSIDE)
    fake = np.random.normal(loc=0,scale=1e-5,size=(pixels))

    lmax = 512

    cl = np.zeros(lmax + 1)
    ells = np.arange(lmax + 1)
    cl[1:] = 1/(ells[1:] * (ells[1:] + 1))

    alm = hp.synalm(cl,lmax=lmax)
    cmb_map = hp.alm2map(alm,NSIDE,lmax)


    hp.mollview(cmb_map)
    x = input("Simulation file name: ")
    plt.title("CMB sim: seed = {}, lmax = {}, nside = {}".format(seed,lmax,NSIDE))
    plt.savefig("../../dump/{}.png".format(x))




'''
Function dictionary
'''

function_dictionary = {
    "function" : function,
    "load_map" : load_map,
    "save_map" : save_map,
    "hist" : histograms,
    "mask_map" : mask_map,
    "skew" : skew,
    "kurt" : kurt,
    "mean" : mean,
    "variance" : var,
    "create_disc" : create_disc,
    "simulate" : simulation,
}

'''
Inputting loop
'''

while (inputting == True):

    inputt = input(" > ")
    if (inputt == "escape"):
        inputting = False
        break
    try:
        function_dictionary[inputt]()
    except:
        print("Something went wrong")
    






