# Stats from SFS
Scripts to calculate pi, dxy, Tajima's D from site frequency spectra.

a) pi (genetic diversity) in the 'pi_tajima_n.R script (it also gives you Tajima's D, the number of segregating sites S, ThetaW (Whatterson's Theta)
b) dxy (absolute divergence) in the dxy.R script
 
This is written for calculating the stats from a list of SFS bootstraps, so please ignore the 'Boots' column in the output file. But it should of course work perfectly for just a single SFS, which you supply without any header.

Example on command line:
`Rscript dxy.R name_of_sfs_file n m`

where n is the number of individuals in pop1, m the number of individuals in pop2. SFS file in format pop1-pop2.sfs

## Compute all stats with one script + get marginalized 1D SFS
computes pi, theta, dxy, da for both populations, and their Fst, marginalizes into 1D SFS from moments formatted 2D SFS (unfolded or folded without mask).

Command:
`python allStats_moments_1Dsfs.py pop${pop1}_${pop2}.sfs ${pop1} ${pop2} pop1size pop2size`
