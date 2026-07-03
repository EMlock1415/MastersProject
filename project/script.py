import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import time
import scipy
import random
import csv
import camb


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

NSIDE = 2048
pixels = hp.nside2npix(NSIDE)

inputting = True
map=[]
sim_map = []
masked_map = []
masked_sim_map = []
disc = []
mask = hp.read_map("../../PlanckMaps/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits")
maskk = np.zeros(pixels)




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
        map = hp.ud_grade(hp.read_map(planckMaps[int(frequency)]),2048)*(10**4) # Scaling to make variance = 1
        if int(frequency) == 100:
            map[45581992] = 0
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

    #print(saved_map[25000000])
    #plt.savefig("dump/mapp.png")

def histograms():
    maps = {
        "map" : map,
        "masked_map" : masked_map,
        "disc" : masked_map[disc]
    }
    working_map = maps[input("what are we working with here???? ")]
    x = plt.hist(working_map,bins=25)
    print(working_map.max())
    print(working_map.min())
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

def dagostino():
    maps = {
        "map" : map,
        "masked_map" : masked_map,
        "disc" : masked_map[disc]
    }
    x = input("what are we working with here???? ")
    working_map = maps[x]
    a = scipy.stats.normaltest(working_map)
    print("D'Agostino K^2 test of {}: {}".format(x,a))

def KStest():

    global masked_map
    global masked_sim_map


    working_masked_map1 = masked_map.compressed() - np.mean(masked_map)
    working_masked_sim_map1 = masked_sim_map.compressed() - np.mean(masked_sim_map)
   

    a = scipy.stats.kstest(working_masked_map1, working_masked_sim_map1)
    print(a)

    
    file_name = input("File name: ")
    plt.ecdf(working_masked_map1[::50], label="Real Map")
    plt.ecdf(working_masked_sim_map1[::50], label="Simulated Map")
    plt.legend()
    plt.savefig("../../dump/{}.png".format(file_name))
    plt.close()


def create_disc():
    phi =  float(input("Phi (north to south, [0,pi])" ) ) # north to south
    theta =  float(input("Theta (west to east, [0,2pi])" ) )
    vector = hp.ang2vec(phi,theta)
    d = hp.query_disc(nside=NSIDE,vec=vector,radius=np.radians(2))
    m = np.arange(pixels)
    m[d] = m.max()
    hp.mollview(m)
    plt.savefig("../../dump/disc.png")
    plt.close()

    global disc
    disc = d

def simulation():
    global sim_map
    global map
    seed = input("Seed")
    np.random.seed(int(seed))
    NSIDE = int(input("NSIDE: "))
    pixels = hp.nside2npix(NSIDE)
    fake = np.random.normal(loc=0,scale=1e-5,size=(pixels))

    lmax = int(input("Highest multipole? "))

    np.random.seed(int(seed))
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.5)
    pars.set_for_lmax(lmax)
    np.random.seed(int(seed))
    results = camb.get_results(pars)
    cl = results.get_cmb_power_spectra(
        params=pars,
        lmax=lmax,
        CMB_unit = 'K',
        raw_cl=True
    )['total'][:,0]


    '''cl = np.zeros(lmax + 1)
    ells = np.arange(lmax + 1)
    cl[1:] = 1/(ells[1:] * (ells[1:] + 1))'''
    np.random.seed(int(seed))

    alm = hp.synalm(cl,lmax=lmax)

    #floats = [float(i) for i in alm]

    
    cmb_map = hp.alm2map(alm,NSIDE,lmax)*(10**4)


    sim_map = cmb_map
    map = cmb_map

    global masked_sim_map
    masked_sim_map = hp.ma(sim_map)
    masked_sim_map.mask = maskk

    hp.mollview(masked_sim_map, unit="K", norm="hist")
    x = input("Simulation file name: ")
    plt.title("CMB sim: seed = {}, lmax = {}, nside = {}".format(seed,lmax,NSIDE))
    plt.savefig("../../dump/{}.png".format(x))
    plt.close()

def simulate_two():
    global map
    seed = input("Seed")
    np.random.seed(int(seed))
    NSIDE = 2048
    lmax = int(input("Highest multipole? "))
    lmin =  int(input("Lowest multipole? "))

    cl = np.zeros(lmax + 1)
    ells = np.arange(lmax + 1)
    cl[1:] = 1/(ells[1:] * (ells[1:] + 1))

    alm = hp.synalm(cl,lmax=lmax)
    
    cmb_map_1 = hp.alm2map(alm,NSIDE,lmax)

    np.random.seed(int(seed))

    cl = np.zeros(lmin + 1)
    ells = np.arange(lmin + 1)
    cl[1:] = 1/(ells[1:] * (ells[1:] + 1))

    alm = hp.synalm(cl,lmax=lmin)

    cmb_map_2 = hp.alm2map(alm,NSIDE,lmin)

    cmb_map = cmb_map_1 - cmb_map_2

    hp.mollview(cmb_map, unit="K", norm="hist")
    x = input("Simulation file name: ")
    plt.title("CMB sim: seed = {}, lmax = {}, nside = {}".format(seed,lmax,NSIDE))
    plt.savefig("../../dump/{}.png".format(x))
    plt.close()

    map = cmb_map


def bulk_simulate():
    start_seed = int(input("Starting Seed: "))
    amount = int(input("How many simulations? "))


    if (amount < 1000):

        with open("../../dump/statistics{}{}.csv".format(start_seed,amount), "a") as file:
            writer = csv.writer(file)
            writer.writerow(["number", "mean", "variance", "skew", "kurtosis", "dagostino"])
        
        for i in range(amount):
            np.random.seed(int(i+start_seed))
            NSIDE = 2048
            lmax = 1024

            pars = camb.CAMBparams()
            pars.set_cosmology(H0=67.5)
            pars.set_for_lmax(lmax)
            results = camb.get_results(pars)
            cl = results.get_cmb_power_spectra(
                params=pars,
                lmax=lmax,
                CMB_unit = 'K',
                raw_cl=True
            )['total'][:,0]
            
            '''cl = np.zeros(lmax + 1)
            ells = np.arange(lmax + 1)
            cl[1:] = 1/(ells[1:] * (ells[1:] + 1))'''

            alm = hp.synalm(cl,lmax=lmax)
            
            cmb_map = hp.alm2map(alm,NSIDE,lmax)

            mean = np.mean(cmb_map)
            variance = np.var(cmb_map)
            skew = scipy.stats.skew(cmb_map)
            kurt = scipy.stats.kurtosis(cmb_map)
            dagostino = scipy.stats.normaltest(cmb_map)

            statistic = [i,mean,variance,skew,kurt,dagostino]
            with open("../../dump/statistics{}{}.csv".format(start_seed,amount), "a") as file:
                writer = csv.writer(file)
                writer.writerow(statistic)
            print(statistic)
    else:
        print("too many simulations")

def bad_sim():
    global map
    fake = np.random.normal(loc=0,scale=1,size=(pixels))
    map = fake

def make_mask():
    global maskk
    res = 20
    size = 25
    vector = hp.ang2vec(np.pi/2, np.pi)
    maskkk = hp.query_disc(nside=NSIDE,vec=vector,radius=np.radians(size))
    for i in range(res):
        vector = hp.ang2vec(np.pi/2, np.pi+i*2*np.pi/res)
        maskkk = np.append(maskkk,hp.query_disc(nside=NSIDE,vec=vector,radius=np.radians(size)))
    maskk = np.zeros(pixels)
    maskk[maskkk] = 1
    maskk = np.logical_or(np.logical_not(mask),maskk)
    
    global masked_map
    masked_map = hp.ma(map)
    masked_map.mask = maskk

def alm_stats():
    maps = {
        "map" : map,
        "masked_map" : masked_map,
        "disc" : masked_map[disc]
    }
    x = input("what are we working with here???? ")
    working_map = maps[x]
    alms = hp.map2alm(working_map)
    print(alms[0])
    print(alms[1])
    

def masked_bulk_sim():
    start_seed = int(input("Starting Seed: "))
    amount = int(input("How many simulations? "))

    with open("../../dump/masked_stats{}{}.csv".format(start_seed,amount), "a") as file:
        writer = csv.writer(file)
        writer.writerow(["number", "mean", "variance", "skew", "kurtosis", "dagostino"])

    for i in range(amount):
        np.random.seed(int(i+start_seed))
        NSIDE = 2048
        lmax = 1024

        pars = camb.CAMBparams()
        pars.set_cosmology(H0=67.5)
        pars.set_for_lmax(lmax)
        results = camb.get_results(pars)
        cl = results.get_cmb_power_spectra(
            params=pars,
            lmax=lmax,
            CMB_unit = 'K',
            raw_cl=True
        )['total'][:,0]
        
        

        alm = hp.synalm(cl,lmax=lmax)
        
        cmb_map = hp.alm2map(alm,NSIDE,lmax)
        masked_cmb = hp.ma(cmb_map)
        masked_cmb.mask = maskk

        if (i < 3):
            hp.mollview(masked_cmb, norm="hist", title="current map")
            plt.savefig("../../dump/masked_cmb_{}.png".format(i))
            plt.close()
        
        mean = np.mean(masked_cmb)
        variance = np.var(masked_cmb)
        skew = scipy.stats.skew(masked_cmb)
        kurt = scipy.stats.kurtosis(masked_cmb)
        dagostino = scipy.stats.normaltest(masked_cmb)

        statistic = [i,mean,variance,skew,kurt,dagostino]
        with open("../../dump/masked_stats{}{}.csv".format(start_seed,amount), "a") as file:
            writer = csv.writer(file)
            writer.writerow(statistic)
        print(statistic)

        


    

'''
Function dictionary
'''

function_dictionary = {
    "function" : function,
    "load_map" : load_map,
    "save_map" : save_map,
    "hist" : histograms,
    "mask_map" : mask_map,
    "make_mask" : make_mask,
    "skew" : skew,
    "kurt" : kurt,
    "mean" : mean,
    "variance" : var,
    "dagostino" : dagostino,
    "kstest" : KStest,
    "create_disc" : create_disc,
    "simulate" : simulation,
    "simulate_two" : simulate_two,
    "bad_sim" : bad_sim,
    "bulk" : bulk_simulate,
    "masked_bulk" : masked_bulk_sim,
    "alm_stats" : alm_stats,
}

'''
Inputting loop
'''

while (inputting == True):

    inputt = input(" > ")
    if (inputt == "escape"):
        inputting = False
        break
    '''try:
        function_dictionary[inputt]()
    except:
        print("Something went wrong")'''
    
    function_dictionary[inputt]()
    






