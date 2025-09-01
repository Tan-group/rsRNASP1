

# rsRNASP1
rsRNASP1: a distance and torsion-angle dependent statistical potential with residue separation for RNA 3D structure evaluation.


# Install

## For Linux

Use Makefile (use gcc compiler):

```
make
```


Note: after compilation, You can manually set environment variable rsRNASP_RNA_HOME to the source directory of rsRNASP-RNA so the energy files can be detected, or run the following script 

bash install.sh 


# Usage

```
# Help is displayed if run rsRNASP1 without parameters

#######################################################
# Calculate rsRNASP1 score for a pdb or a list of pdbs
#######################################################
Usage: ./bin/rsRNASP1 pdb 
   or: ./bin/rsRNASP1 [ options ] 
Options:
   pdb [ pdb2 pdb3 ...], input RNA structures in pdb format
   -d directory,         OPTIONAL, override default directory of energyfiles
                         default: rsRNASP_rna/data/energyfiles
   -l pdblist,           A list of absolute paths to pdb files (plain text) UTF-8 encoding

```


# Example of usage

```
# Run in the example dir

../bin/rsRNASP1 *.pdb

#output:
#1a9nR.pdb -3146.575662
#1h4sT.pdb -7757.550000

```

If you have any questions about rsRNASP1, please contact us by the email: zjtan@whu.edu.cn .

Reference:                                      
[1] Lou E, Zheng C, Yu S, Tan YL, & Tan ZJ. rsRNASP1: a distance- and dihedral-dependent statistical potential for RNA 3D structure evaluation. Biophys J 2025;124:2740–53. doi:10.1016/j.bpj.2025.07.013


