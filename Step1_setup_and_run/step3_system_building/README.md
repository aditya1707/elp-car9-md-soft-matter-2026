## Step 3: Molecular Dynamics Simulation Protocol

This script documents the reference GROMACS protocol used in the paper.
It is provided for transparency and reproducibility, not as a prescriptive workflow.

### Key design choices
- Grid of 9 aligned ELP-Car9 chains
- Packmol-generated water followed by re-solvation
- Gradual relaxation with frozen protein chains
- Semi-isotropic NPT (xy coupled, z decoupled)

### Notes
- Users are expected to adapt this protocol to their own HPC environment
- mdp files reflect one reasonable choice, not a universal best practice
- GPU flags assume GROMACS ≥2022 and CUDA-enabled nodes