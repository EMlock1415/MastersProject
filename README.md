# Masters Project

The following is a brief readme file to introduce the code for this project (hosted at https://github.com/EMlock1415/MastersProject)


- Currently used maps are found at https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/scripts/get_planck_maps_primary.sh (as of 10/08/2026)

- Component separation common mask file found at https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/previews/COM_Mask_CMB-common-Mask-Int_2048_R3.00/index.html


Script.py contains a set of functions to be typed in the console that can load in maps, mask them, simulate CMB skies, and perform tests. Some functions only work when an unmasked *and* a  masked map have been prepared. Some functions that should be ran to initialise are

- load_map: Loads a map. Will prompt you to enter a number for the frequency in GHz of the map you want to load.
- simulate: Simulates a map. Prompts you to enter a number for the NSIDE and lmax parameter.
- mask_map: Masks the loaded map with the common mask
- make_mask: Masked the loaded map with the common mask + thick band across galactic plane

The rest of the functions are for running tests and are relatively self explanatory
