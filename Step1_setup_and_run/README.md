# System Setup and Simulation Workflow

This directory documents the complete molecular dynamics workflow used in this study.

The goal is to construct a polymer-brush-like system starting from a single ELP–peptide chain
and simulate its behavior in explicit solvent.

---

## Workflow Overview

1. **step1_initial_structure**
   - Align and preprocess a single peptide chain
   - Define a consistent reference orientation

2. **step2_system_building**
   - Replicate the chain to form a multi-chain brush-like system
   - Pack, solvate, generate topology, and index files

3. **step3_simulation_protocol**
   - Energy minimization
   - Equilibration (NVT/NPT)
   - Production MD

Each step contains a short README describing the purpose of the scripts.
