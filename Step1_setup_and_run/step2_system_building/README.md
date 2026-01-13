# Step 2 — System Building (GROMACS + PACKMOL)

This step builds a packed multi-chain system (9×9 grid) starting from a single aligned ELP chain PDB.

## What this step does

1. **Remove hydrogens** from the tleap-generated PDB so `pdb2gmx` can re-add them consistently.
2. Run **`gmx pdb2gmx` (Amber99SB)** on the single chain to generate:
   - a `.gro` (GROMACS structure)
   - a topology (`topol.top`)
3. Use **PACKMOL 20.15.1** to replicate/place the chain into a **9×9 grid** in a **5 nm × 5 nm × 12 nm** box and add water molecules.
4. Run **`gmx pdb2gmx` again** on the PACKMOL-packed PDB to produce the **final** GROMACS system:
   - final `.gro`
   - final topology files
5. Create a **GROMACS index file** `system.ndx` (Protein group by default).

> **Why remove H first?**  
> `tleap` can generate hydrogen naming/geometry that conflicts with what `pdb2gmx` expects.  
> Removing H ensures GROMACS rebuilds hydrogens according to **Amber99SB** rules.   
> We run `editconf` to give us the pdb with the correct hydrogens because `packmol` requires a pdb file as input. 

---

## Folder layout

- `input/`
  - `elp_car9_aligned.pdb` : aligned single-chain PDB (starting structure) for the elp-linker-car9 system generated from step 1. 
  - `water.pdb` : water geometry used by PACKMOL (TIP3P-compatible)
- `output/`
  - generated `.pdb`, `.gro`, topology files (`topol.top.1` and `topol.top`), and `system.ndx`
- `packmol/`
  - `create_packmol.sh` : generates `packmol_input.inp` (hardcodes box size, counts, distances, etc.)
  - `packmol_input.inp` : generated PACKMOL input file
- `scripts/`
  - `removeH.py` : removes hydrogens from PDB
  - `setup.sh` : runs the full Step 2 pipeline

---

## Requirements (HPC assumed)

- GROMACS available as `gmx`
- PACKMOL installed on the cluster (path may be cluster-specific)
- Python available for `removeH.py`

---

## How to run

From `step2_system_building/scripts/`:

```bash
bash setup.sh
```
Note: Install and load GROMACS, install and update the path to PACKMOL 20.15.1 in the script, update and check the paths to all the inputs and scripts that `setup.sh` runs (`elp_car9_aligned.pdb`, `water.pdb`, `removeH.py`,`create_packmol.sh`, `packmol_input.inp`) or run each command in `setup.sh` independently. 