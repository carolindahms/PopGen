#!/usr/bin/env python

import matplotlib
matplotlib.use('PDF')
import moments
import sys
import os
import numpy as np
import pandas as pd
from moments import Spectrum

"""
Carolin Dahms, Sept 2025

Description:
 Computes pi, theta, TD, dxy, da, Fst, marginalised 1D site frequency spectra.
Input: bootstraps of moments formatted 2D SFS (unfolded or folded without mask). 
Expects sfs bootstraps consistign of 2 lines each: sfs + header line in format '2n+1 2m+1 unfolded' 
(note: always 'unfolded' irrespective of polarization to make this script work)

Execute: 
 python allStats_moments_1Dsfs.py pop1_pop2.sfs pop1 pop2 n m [fold]
n = number of individuals in pop1, m = number of individuals in pop2
The 'fold' arg is optional if you want the marginalised sfs to be folded.

"""
def process_sfs(infile, pop1, pop2, output_file, dxy, pi1, pi2, fold=False):
    data = Spectrum.from_file(infile, mask_corners=False)
    if fold:
        data = data.fold()
    ns = data.sample_sizes
    seg = np.around(data.S(), 2)
    sfs_sum = np.around(np.sum(data), 0)
    Fst = np.around(data.Fst(), 4)
    #sdata = moments.Spectrum.scramble_pop_ids(data, mask_corners=True)
    #sFst = np.around(sdata.Fst(), 4)

    sfs1 = moments.Spectrum.marginalize(data, [1])
    sfs2 = moments.Spectrum.marginalize(data, [0])
    sfs1_sum = np.sum(sfs1)
    sfs2_sum = np.sum(sfs2)

    pi1 = (np.around(sfs1.pi(), 2))/sfs_sum
    theta1 = np.around(sfs1.Watterson_theta(), 2)
    TD1 = np.around(sfs1.Tajima_D(), 4)

    pi2 = (np.around(sfs2.pi(), 2))/sfs_sum
    theta2 = np.around(sfs2.Watterson_theta(), 2)
    TD2 = np.around(sfs2.Tajima_D(), 4)
    
    da = dxy - (pi1 + pi2) / 2

    # Create a list of SFS1 values
    sfs1_values = [np.around(v, 2) for v in sfs1]
    #print(sfs1_values)
    # Create a list of SFS2 values
    sfs2_values = [np.around(v, 2) for v in sfs2]

    # Determine the number of values in each SFS
    sfs1_num = len(sfs1_values)
    sfs2_num = len(sfs2_values)

    # Build the output line with the additional columns
    output_line = f"{pop1}\t{pop2}\t{seg}\t{Fst}\t{pi1}\t{theta1/sfs_sum}\t{TD1}\t{pi2}\t{theta2/sfs_sum}\t{TD2}\t{sfs_sum}\t{dxy}\t{da}"
    for value in sfs1_values:
        output_line += f"\t{value}"
    for value in sfs2_values:
        output_line += f"\t{value}"
    
    # Write to the output file
    with open(output_file, 'a') as out_file:
        out_file.write(output_line + "\n")
    
    #print(sfs1)
    #print(sfs2)    
    #print(sfs_sum)
    #print (ns)
    return pi1, pi2

# Function to compute dxy
def compute_dxy(infile, pop1, pop2, n1, n2):

    # Read the 2D SFS data from the input file
    dat = pd.read_csv(infile, sep="\t", header=None)
    #data = Spectrum.from_file(infile, mask_corners=False)
    
    # Flatten the DataFrame to ensure proper shape
    dat = dat.values.flatten()

    # Adjust the sample sizes to match the format of the script you provided
    n = n1 * 2 + 1
    m = n2 * 2 + 1

    # Create the allele frequency vectors for the two populations
    P1 = np.linspace(0, 1, n)
    P2 = np.linspace(0, 1, m)

    # Create the weights
    weights = np.empty(n * m)
    c = 0
    for i in range(n):
        p1 = P1[i]
        for j in range(m):
            p2 = P2[j]
            weights[c] = p1 * (1 - p2) + p2 * (1 - p1)
            c += 1

    # Compute dxy
    for i in range(len(dat)):
        d = dat[i]  # Get the current count
        d_counts = np.fromstring(d, sep=' ')  # Convert space-separated string to numpy array
        nsites = int(np.sum(d_counts))  # Sum as integer
        dxy = np.sum(d_counts * weights) / nsites if nsites > 0 else 0  # Prevent division by zero    '''
    
    return dxy
    ''' 
    #Compute da
def compute_da(dxy, pi1, p2):
    
    da = dxy - (pi1+pi2)/2
    
    return da
    
'''
# Main execution
if __name__ == "__main__":
    # Get the input file and population names from command line arguments
    infile = sys.argv[1]
    pop1 = sys.argv[2]
    pop2 = sys.argv[3]
    n1 = int(sys.argv[4])
    n2 = int(sys.argv[5])
    fold = len(sys.argv) > 6 and sys.argv[6].lower() == 'fold'

    output_file = f"Sumstats_{pop1}_{pop2}_1D.txt"

    # Write the header to the output file
    with open(output_file, 'w') as out_file:
        out_file.write("pop1\tpop2\tsegregating_sites\tFst\tpi1\ttheta1\tTD1\tpi2\ttheta2\tTD2\tsfs_sum\tdxy\tda")
        for i in range(1, n1 * 2 + 2):
            out_file.write(f"\tsfs1_{i}")
        for i in range(1, n2 * 2 + 2):
            out_file.write(f"\tsfs2_{i}")
        out_file.write(f"\n")
        
    # Read all lines from the input file
    with open(infile, 'r') as f:
        lines = f.readlines()
        
    if len(lines) % 2 != 0:
        print(f"Error: The input file '{infile}' has an odd number of lines. Please check the file format.")
        sys.exit(1)
        
    # Loop through lines, processing every two lines as one SFS
    for i in range(0, len(lines), 2):
        #print(f"Processing lines {i} and {i+1}")
        
        # Create temporary SFS files for compute_dxy and process_sfs
        temp_dxy_file = "temp_dxy.txt"
        temp_sfs_file = "temp_sfs.txt"
        
        with open(temp_sfs_file, 'w') as temp_sfs:
            temp_sfs.write(lines[i])  # Write the header
            temp_sfs.write(lines[i + 1])  # Write the SFS data
        
        # Write the SFS data to the temporary files
        with open(temp_dxy_file, 'w') as temp_dxy:
            temp_dxy.write(lines[i + 1])
          

        # Compute dxy for the current SFS
        # dxy = compute_dxy(temp_dxy_file, pop1, pop2)
        dxy = compute_dxy(temp_dxy_file, pop1, pop2, n1, n2)
        
        # Process the temporary SFS file for other statistics
        pi1 = 0
        pi2 = 0

        pi1, pi2 = process_sfs(temp_sfs_file, pop1, pop2, output_file, dxy, pi1, pi2, fold)
        
        #da = compute_da(dxy, pi1, pi2)
        #print("absolute divergence =", da, pi1, pi2)
        

       

    # Optionally, remove the temporary files after processing
    os.remove(temp_dxy_file)
    os.remove(temp_sfs_file)

