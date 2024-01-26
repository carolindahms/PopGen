#Stats from SFS
Scripts to calculate pi, dxy, Tajima's D from site frequency spectra.

a) pi (genetic diversity) in the 'pi_tajima_n.R script (it also gives you Tajima's D, the number of segregating sites S, ThetaW (Whatterson's Theta)
b) dxy (absolute divergence) in the dxy.R script
 
This is written for calculating the stats from a list of SFS bootstraps, so please ignore the 'Boots' column in the output file. But it should of course work perfectly for just a single SFS, which you supply without any header.

Example on command line:
`Rscript dxy.R name_of_sfs_file n m`

where n is the number of individuals in pop1, m the number of individuals in pop2. SFS file in format pop1-pop2.sfs
