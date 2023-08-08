# BHRefinementAnalysis
These are the python scripts I created to analyze my implementation of super-Lagrangian refinement around black holes in AREPO.
Note that the AREPO code with my BH refinement is not here; you can find that on a bh-refinement branch on the AREPO repository (as of writing).
Also note that the data is not available in this repository. Some data is archived on my HiPerGator account, though it is recommended that you run your own simulation instead.

Each script includes a description at the top, a path to the data, and comments throughout; they generally end with creating and saving a plot.
Be sure to adjust the directory paths according to your setup.

Extracting data from an AREPO run requires reading hdf5 (binary) files. During my time working on this project (2020-2023), `torreylabtools` were required to read in the data; this is written in python 2.7. 
To bring the analysis up to python3, Jonah Rose developed DataLoader. Both of these packages are imported in these scripts; feel free to use one, both, or your own hdf5 reader.
