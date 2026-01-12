# Step 2: Analysis

This directory contains analysis scripts used to generate the figures and tables
reported in the paper. The following is a list folders, the rationale behind each of the analyses, and some information about specific scripts and corresponding outputs. 

1. **figure4_density**: To determine how sequence variations modulate spatial organization within the micelle core, we quantified the z-axis density profiles of silica-binding segments relative to the ELP-linker center-of-mass (COM), revealing the degree of interfacial mixing versus segregation that correlates with distinct self-assembly and templating behaviors.

    - `get_densities_and_com.sh`: bash script to compute density profiles of the Car9/variant tags (`density.xvg`) along the z-axis and COM of the ELP domains (`com_ELP.xvg`).
    - `plot_density.py`: python script generate the figure (`Car9_density.png`) for the density profile ELP-linker-Car9 along the z-axis and the z-coordinate of the COM of the ELP-domains

2. **figureS10_energies**: To rationalize sequence-dependent self-assembly, we quantified interchain and chain–water LJ and Coulombic interactions, linking differences in energetic balance to hydration versus self-association tendencies observed across variants.  
    - `get_energies.sh`: bash script to 
        - make an index file containing protein and water indices together 
        - create a GROMACS mdp file with energy groups listed (`npt2_prot_water.mdp`)
        - rerun `grompp` and `mdrun` commands to compute energies, and
        - extract the protein-protein and protein-water Lennard Jones and Coulombic energies (`energy.xvg`) using gromacs. 

3. **figureS13_end_to_end_distances**: To assess peptide compactness and conformational variability, we computed KDEs of end-to-end distances for each tag across chains and trials, enabling direct comparison of sequence-dependent compactness within the polymer-brush environment.
   - `e2e_car9_plumed.dat`: plumed script to compute and print the distance between the COM of the first and last amino acid residues of Car9/variant tag domains in each of the 9 ELP-linker-tag chains placed in our box (`car9_COLVAR`).
   - `plot_kde.py`: Python script to read the plumed output, extract data for the last 100 ns of simulation trajectory, and plot the kernel density estimate (KDE) of the end-to-end distances of the 9 chains in our simulation box for each trial (`Car9_kde.png`). 

4. **figureS14_clustering**: To quantify how broadly each peptide explores conformational space, we performed RMSD-based clustering on concatenated trajectories and analyzed the cumulative distribution of cluster sizes across all chains and trials.
    -  `clustering.sh`: bash script to 
        - extract the last 100 ns of trajectory for each trial
        - separate the trajectory for each of the 9 chains in the system and center and fit it on the alpha Carbon of the peptide tag (Car9/variants)
        - create reference pdb files for clustering trajectory
        - concatenate the trajectories for the 9 chains for each trial
        - concatenate the trajectories from 3 different trials into one large 2.7 µs trajectory
        - perform the clustering analysis using the GROMOS algorithm. 
    - `plot_CDF.py`: plots the cumulative distribution function of cluster sizes for ELP-linker-Car9 and ELP-linker-P9AG10A at rmsd cutoff value of 0.15 nm (`cluster_size_CDF`). 

5. **tableS2_Shannon_entropy**: To obtain a single metric capturing conformational diversity, we computed the Shannon entropy of cluster populations, where higher entropy reflects greater structural heterogeneity and sampling breadth.
    - `calculate_shannon_entropy.py`: calculates the Shannon entropy of cluster sizes of each variant across the 2.7 µs trajectories at any given cutoff (in this case 0.15)


For more details, complete images, and inferences, go through the results, and supplementary section of the paper. 