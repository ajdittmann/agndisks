# AGN disks
Calculate disk models for AGN. Currently only simple α-disk models are implemented (e.g. [Shakura and Sunyaev 1973](https://ui.adsabs.harvard.edu/abs/1973A%26A....24..337S/abstract)). 

If you would like to use tabluated opacities (from https://github.com/zhuzh1983/combined-opacity), then please clone this repository using 

`git clone --recurse-submodules https://github.com/ajdittmann/agndisks.git`. 

If you would like to use a simple Kramers' opacity, or use other opacity tables,
simply use

`git clone https://github.com/ajdittmann/agndisks.git`.

Resolution parameters are currently controlled using the `rmax`, `rmin`, and `Nr` variables. 
Free parameters in the accretion disk model are controlled by setting the `eta` (fraction of Eddington accretion rate), `M` (SMBH mass in M/M☉), `alpha` (SS73 α parameter), `X` (hydrogen mass fraction), and `Z` (metal mass fraction). `eps`, the efficiency at which rest mass is converted into energy during accretion onto the SMBH, is used to define the Eddington accretion rate.  

To compute a disk model, simply run `python diskcomp.py`. 
