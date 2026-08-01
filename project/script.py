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
    hp.mollview(saved_map, norm="hist", title="", margins=[0.1,0.1,0.1,0.1])
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


    working_masked_map1 = (masked_map.compressed() - np.mean(masked_map))/np.sqrt(np.var(masked_map))
    working_masked_sim_map1 = (masked_sim_map.compressed() - np.mean(masked_sim_map))/np.sqrt(np.var(masked_sim_map))

    print(np.mean(working_masked_map1))
    print(np.var(working_masked_map1))

    working_masked_map2 = (masked_map.compressed() - np.mean(masked_map))
    working_masked_sim_map2 = (masked_sim_map.compressed() - np.mean(masked_sim_map))
   

    a = scipy.stats.kstest(working_masked_map1, working_masked_sim_map1)
    b = scipy.stats.kstest(working_masked_map2, working_masked_sim_map2)
    print("var adjusted")
    print(a)
    print("not var adjusted")
    print(b)

    
    file_name = input("File name 1: ")
    plt.ecdf(working_masked_map1[::50], label="Real Map")
    plt.ecdf(working_masked_sim_map1[::50], label="Simulated Map")
    plt.legend()
    plt.savefig("../../dump/{}.png".format(file_name))
    plt.close()

    file_name = input("File name 2: ")
    plt.ecdf(working_masked_map2[::50], label="Real Map")
    plt.ecdf(working_masked_sim_map2[::50], label="Simulated Map")
    plt.legend()
    plt.savefig("../../dump/{}.png".format(file_name))
    plt.close()

def multivariate():
    global masked_map
    #global masked_sim_map


    working_masked_map1 = (masked_map.compressed() - np.mean(masked_map))/np.sqrt(np.var(masked_map))
    #working_masked_sim_map1 = (masked_sim_map.compressed() - np.mean(masked_sim_map))/np.sqrt(np.var(masked_sim_map))

    #data = np.array([working_masked_map1],[working_masked_map1])
    
    #S = np.cov(data)
    #S_inv = np.linalg.inv(S)

    

    D = np.linalg.matrix_transpose(working_masked_map1) @ 1/np.var(working_masked_map1) @ working_masked_map1

    skew = np.sum(D**3)/(np.len(D)**2)
    kurt = np.sum(np.diag(D)**2)/np.len(D)

    print("Multivariate skew: ")
    print(skew)
    print("Multivariate kurtosis: ")
    print(kurt)

    

    #scipy.spatial.distance.mahalanobis(masked_map,masked_map,)


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
        writer.writerow(["number", "mean", "variance", "skew", "kurtosis", "dagostino","dag_p"])

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

        

        statistic = [i,mean,variance,skew,kurt,dagostino.statistic,dagostino.pvalue]
        with open("../../dump/masked_stats{}{}.csv".format(start_seed,amount), "a") as file:
            writer = csv.writer(file)
            writer.writerow(statistic)
        print(statistic)

'''
Large scale production
'''

def bulk_ks():
    NSIDE = 2048
    pixels = hp.nside2npix(NSIDE)
    frequencies = [30,44,70,100,143,217,353,545,857]
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

    
    
    

    for j in range(0,9):

        with open("../../dump/KSTEST{}.csv".format(frequencies[j]), "a") as file:
            writer = csv.writer(file)
            writer.writerow(["Seed", "Test Result", "p-value", "location"])
        map = hp.ud_grade(hp.read_map(planckMaps[frequencies[j]]),2048)*(10**4) # Scaling to make variance = 1
        if frequencies[j] == 100:
            map[45581992] = 0

        masked_map = hp.ma(map)
        masked_map.mask = maskk

       

        

        x = 50

        

        for i in range(0,x-1):

            '''Simulation'''

            seed = i
            np.random.seed(int(seed))
            pixels = hp.nside2npix(NSIDE)
            

            lmax = 2048

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

            alm = hp.synalm(cl,lmax=lmax)
                
            sim_map = hp.alm2map(alm,NSIDE,lmax)

            masked_sim_map = hp.ma(sim_map)
            masked_sim_map.mask = maskk

            '''Tests'''

            working_masked_map1 = (masked_map.compressed() - np.mean(masked_map))/np.sqrt(np.var(masked_map))
            working_masked_sim_map1 = (masked_sim_map.compressed() - np.mean(masked_sim_map))/np.sqrt(np.var(masked_sim_map))

            a = scipy.stats.kstest(working_masked_map1, working_masked_sim_map1)

            b = [seed, a.statistic,a.pvalue,a.statistic_location]

            with open("../../dump/KSTEST{}.csv".format(frequencies[j]), "a") as file:
                writer = csv.writer(file)
                writer.writerow(b)
            print(a)


def bulk_ks_disks():

    NSIDE = 2048
    pixels = hp.nside2npix(NSIDE)
    frequencies = [30,44,70,100,143,217,353,545,857]
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

    with open("../../dump/KSTESTDISC.csv", "a") as file:
        writer = csv.writer(file)
        writer.writerow(["Disc", "Test Results"])
    
    for i in range (0,9):
        map = hp.ud_grade(hp.read_map(planckMaps[frequencies[i]]),2048)*(10**4) # Scaling to make variance = 1
        if frequencies[i] == 100:
            map[45581992] = 0

        masked_map = hp.ma(map)
        masked_map.mask = maskk

        
        

        x = 10
        y = 10
        for i in range(0,x):
            phi =  float(np.random.uniform(0,np.pi) ) # north to south
            theta =  float(np.random.uniform(0,np.pi*2))
            vector = hp.ang2vec(phi,theta)
            d = hp.query_disc(nside=NSIDE,vec=vector,radius=np.radians(2))

            with open("../../dump/KSTESTDISC.csv", "a") as file:
                    writer = csv.writer(file)
                    writer.writerow([vector,"-----------------"])
            

            for i in range(0,y):

                phi1 =  float(np.random.uniform(0,np.pi) ) # north to south
                theta1 =  float(np.random.uniform(0,np.pi*2))
                vector1 = hp.ang2vec(phi1,theta1)
                e = hp.query_disc(nside=NSIDE,vec = vector1, radius=np.radians(2) )
                a = scipy.stats.kstest(masked_map[d],masked_map[e])
                with open("../../dump/KSTESTDISC.csv", "a") as file:
                    writer = csv.writer(file)
                    writer.writerow([vector1,a])
                print(a)
        
def north_south_contrast():
    NSIDE = 2048
    pixels = hp.nside2npix(NSIDE)
    frequencies = [30,44,70,100,143,217,353,545,857]
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

    vector1 = hp.ang2vec(0,np.pi)
    vector2 = hp.ang2vec(np.pi,np.pi)

    with open("../../dump/asymetrytest.csv", "a") as file:
        writer = csv.writer(file)
        writer.writerow(["frequency","m1","m2", "s1", "s2", "k1", "k2","ksvalue","ksp,ksloc"])
    
    for i in range (0,9):
        map = hp.ud_grade(hp.read_map(planckMaps[frequencies[i]]),2048)*(10**4) # Scaling to make variance = 1
        if frequencies[i] == 100:
            map[45581992] = 0



        masked_map = hp.ma(map)
        masked_map.mask = maskk

        disk1 = hp.query_disc(nside=NSIDE,vec=vector1,radius=np.radians(30))
        disk2 = hp.query_disc(nside=NSIDE,vec=vector2,radius=np.radians(30))

        mean1 = np.mean(masked_map[disk1])
        mean2 = np.mean(masked_map[disk2])

        skew1 = scipy.stats.skew(masked_map[disk1])
        skew2 = scipy.stats.skew(masked_map[disk2])

        kurt1 = scipy.stats.kurtosis(masked_map[disk1])
        kurt2 = scipy.stats.kurtosis(masked_map[disk2])

        test1 = (masked_map[disk1]-np.mean(masked_map[disk1]))/np.var(masked_map[disk1])
        test2 = (masked_map[disk2]-np.mean(masked_map[disk2]))/np.var(masked_map[disk2])

        a = scipy.stats.kstest(test1,test2)
        print(a)

        with open("../../dump/asymetrytest.csv", "a") as file:
            writer = csv.writer(file)
            writer.writerow([frequencies[i],mean1,mean2,skew1,skew2,kurt1,kurt2,a.statistic,a.pvalue,a.statistic_location])

def masked_real_stats():
    NSIDE = 2048
    pixels = hp.nside2npix(NSIDE)
    frequencies = [30,44,70,100,143,217,353,545,857]
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

    with open("../../dump/REALMAPS.csv", "a") as file:
        writer = csv.writer(file)
        writer.writerow(["Frequency", "mean", "var", "skew", "kurt", "dagostino", "dag_p"])
    
    for i in range (0,9):
        map = hp.ud_grade(hp.read_map(planckMaps[frequencies[i]]),2048)*(10**4) # Scaling to make variance = 1
        if frequencies[i] == 100:
            map[45581992] = 0

        masked_map = hp.ma(map)
        masked_map.mask = maskk

        
        mean = np.mean(masked_map)
        var = np.var(masked_map)
        skew = scipy.stats.skew(masked_map)
        kurt = scipy.stats.kurtosis(masked_map)
        a = scipy.stats.normaltest(masked_map)

        with open("../../dump/REALMAPS.csv", "a") as file:
            writer = csv.writer(file)
            writer.writerow([frequencies[i], mean, var, skew, kurt, a.statistic,a.pvalue])





        
    
        


        

            
        
        

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
    "multivariate" : multivariate,
    "bulk_ks" : bulk_ks,
    "bulk_ks_disc" : bulk_ks_disks,
    "ns" : north_south_contrast,
    "real_stats" : masked_real_stats,
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
    
    #function_dictionary[inputt]()
    






