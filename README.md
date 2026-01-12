# Molecular Dynamics Workflow and Analysis of ELP-Car9 variants controlling micelle-vesicle transitions.

This repository contains the molecular dynamics simulation setup and analysis pipelines used in our recent Soft Matter paper:
[Genetic control of morphological transitions in a coacervating protein template.](https://doi.org/10.1039/D5SM01047K)

---

## Repository Structure (High-Level)

Below is a concise description of the repository layout and the purpose of the main top-level folders and files.

- `Step1_setup_and_run/`: System setup and simulation protocols. This directory documents the complete molecular dynamics workflow used in this study.
The goal is to construct a polymer-brush-like system starting from a single ELP–peptide chain and simulate its behavior in explicit solvent.
    - `step1.1_initial_structure_prep/` 
        - Align and preprocess a single peptide chain
        - Define a consistent reference orientation
    - `step1.2_system_building/` 
        - Replicate the chain to form a multi-chain polymer-brush-like system
        - Pack, solvate, generate topology, and index files
    - `step1.3_simulation_protocol/` 
        - GROMACS `mdp` files and run scripts

- `Step2_analysis/`: Analysis scripts and pipelines used to generate the figures
    and tables reported in the paper and the [supplementary information.](https://pubs.rsc.org/en/content/articlelanding/2026/sm/d5sm01047k) 
    - `figure4_density/`
            - Computes z-axis density profiles of the tag domains relative to the
                    ELP-linker center-of-mass and plot density profiles and COM traces.
    - `figureS10_energies/`
            - Decomposes interaction energies (protein–protein and protein–water;
                    LJ and Coulombic) to rationalize sequence-dependent assembly.
    - `figureS13_end_to_end_distances/`
            - Computes end-to-end distances with PLUMED and plot kernel density
                    estimates for comparative compactness analyses.
    - `figureS14_clustering/`
            - Perform RMSD-based clustering across chains and trials and plot
                    cluster-size cumulative distribution functions (CDFs).
    - `tableS2_Shannon_entropy/`
            - Calculates Shannon entropy of cluster populations as a single metric
                    of conformational diversity.

---

## Key Tools & Methods

- **Molecular Dynamics:** GROMACS (2022.x)
- **System construction:** Packmol, tleap (AmberTools)
- **Trajectory analysis:** PLUMED, custom Python (NumPy, Matplotlib)
- **Structural alignment & visualization:** VMD (Tcl)
- **Automation:** Bash scripting

---

## Notes
- Each step contains a short README describing the purpose of the scripts.
- Paths, filenames, and some parameters are system-specific
- Scripts are provided as **reference implementations**
- Users are expected to adapt these scripts to their own systems
