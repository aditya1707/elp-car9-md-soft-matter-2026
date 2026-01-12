#!/usr/bin/env bash
# NOJUMP
set -euo pipefail

# === User-specific settings (edit if needed) ===
IDX=chains_ca.ndx                 # your index file with per-chain SBP C-alpha groups
TRIALS=(2 4 5)                    # the trials you want to use
CHAINS=(1 2 3 4 5 6 7 8 9)        # chain IDs to process
PEPTIDES=(car9 p9ag10a)           # Tag variants to process for clustering
XTC_NAME=npt2.xtc                 # trajectory filename inside each trial dir
TPR_NAME=npt2.tpr                 # tpr filename inside each trial dir
GROUP_OFFSET=25               # in $IDX, chain1..chain9 are groups 26..34

# Time window in ps (last 100 ns = 200 ns to 300 ns)
B_PS=200000
E_PS=300000

# === 1) Unwrap PBC and cut to 200–300 ns once per trial (full system) ===
for t in "${TRIALS[@]}"; do
    for pep in "${PEPTIDES[@]}";do

        xtc="${t}/${pep}/${XTC_NAME}"
        tpr="${t}/${pep}/${TPR_NAME}"

        echo "[Trial $t] PBC unwrap + slice ${B_PS}-${E_PS} ps -> ${t}/${pep}/nojump_200-300ns.xtc"
        # Select "0" (System) when prompted
        printf "0\n" | gmx trjconv -s "$tpr" -f "$xtc" \
            -o "${t}/${pep}/nojump_200-300ns.xtc" \
            -pbc nojump -b "$B_PS" -e "$E_PS"
    done
done
# === 2) For each chain: center+fit on its SBP Cα and output SBP Cα only (no subsampling) ===
for t in "${TRIALS[@]}"; do
    for pep in "${PEPTIDES[@]}";do
        tpr="${t}/${pep}/${TPR_NAME}"
        inxtc="${t}/${pep}/nojump_200-300ns.xtc"

        for c in "${CHAINS[@]}"; do
            grp_idx=$(( c + GROUP_OFFSET ))
            outxtc="${t}/${pep}/chain${c}_T${t}_200-300_fit.xtc"

            echo "[Trial $t | Peptide $pep | Chain $c] center+fit on group #$grp_idx -> $outxtc"
            # Prompts: (1) group for centering/fit (2) output group — use the same numeric group twice
            echo -e "${grp_idx}\n${grp_idx}\n${grp_idx}\n" | gmx trjconv -s "$tpr" -f "$inxtc" -n "$t/$pep/$IDX" \
                -center -ur compact -fit rot+trans -o "$outxtc"
        done
    done
done

# === 3) Create reference PDBs for clustering ===
# Global reference (from Chain 1, first frame at 200 ns) for pooling across all chains

REF_TRIAL="${TRIALS[0]}"

for pep in ${PEPTIDES[@]};do
    ref_grp_idx=$(( 1 + GROUP_OFFSET ))   # group for chain 1 (26)

    echo "[Ref] Writing global SBP_ref.pdb from chain1_T${REF_TRIAL}_200-300_fit.xtc at ${B_PS} ps"
   
    echo -e "$ref_grp_idx" | gmx trjconv -s "${REF_TRIAL}/${pep}/${TPR_NAME}" \
    -f "${REF_TRIAL}/${pep}/nojump_200-300ns.xtc" -n "${REF_TRIAL}/${pep}/$IDX" -dump "$B_PS" -o "${pep}_ref.pdb"

    # Per-chain references (for per-chain clustering)
    for c in "${CHAINS[@]}"; do
      grp_idx=$(( c + GROUP_OFFSET ))
      echo "[Ref] Writing chain${c}_ref.pdb from chain${c}_T${REF_TRIAL}_200-300_fit.xtc at ${B_PS} ps"
      echo -e "$grp_idx\n" | gmx trjconv -s "${REF_TRIAL}/${pep}/${TPR_NAME}" \
      -f "${REF_TRIAL}/${pep}/nojump_200-300ns.xtc" -n "${REF_TRIAL}/${pep}/$IDX" -dump "$B_PS" -o "${pep}_chain${c}_ref.pdb"

    done
done
# === 4) Concatenate all chains per trial (100 ns each) -> then across trials ===
# 4A) For each trial, pool chains 1..9 (each 100 ns) into a 0.9 µs file

for t in "${TRIALS[@]}"; do
  for pep in "${PEPTIDES[@]}";do
    > ${t}/${pep}/settime_allchains_T${t}.txt
    offset=0
    for _ in "${CHAINS[@]}"; do
      echo $offset >> ${t}/${pep}/settime_allchains_T${t}.txt
      offset=$(( offset + 100000 ))   # 100 ns in ps per chain block
    done

    echo "[Cat] trial ${t}: all chains (200–300 ns each) -> allchains_T${t}_0p9us.xtc"
    gmx trjcat -f $(printf "${t}/${pep}/chain%d_T${t}_200-300_fit.xtc " "${CHAINS[@]}") \
              -o ${t}/${pep}/allchains_T${t}_0p9us.xtc -cat -settime < ${t}/${pep}/settime_allchains_T${t}.txt
  done
# done

# 4B) Pool across trials -> one 2.7 µs trajectory

for pep in ${PEPTIDES[@]};do

  # Build per-peptide settime file: 3 blocks × 0.9 µs = 900000 ps
  > ${pep}_settime_alltrials.txt
  offset=0   
  for _ in "${TRIALS[@]}"; do
    echo ${offset} >> ${pep}_settime_alltrials.txt
    offset=$(( offset + 900000 ))
  done

  files=()
  for t in "${TRIALS[@]}";do
    files+=( "${t}/${pep}/allchains_T${t}_0p9us.xtc" )
  done


  echo "[Cat] ${pep}: all trials pooled (chains already pooled per trial) -> ${pep}_allchains_alltrials_2p7us.xtc"
  gmx trjcat -f "${files[@]}" \
             -o "${pep}_allchains_alltrials_2p7us.xtc" \
             -cat -settime < "${pep}_settime_alltrials.txt"

done

# === 5) Clustering of the final 2.7 µs trajectory ===
mkdir combined_clusters
for pep in "${PEPTIDES[@]}"; do
  mkdir combined_clusters/${pep}
  printf "0
  0
  "| gmx cluster -f "${pep}_allchains_alltrials_2p7us.xtc" \
              -s "${pep}_ref.pdb" \
              -method gromos -cutoff 0.15 -fit yes \
              -cl "combined_clusters/${pep}/${pep}_reps.pdb" \
              -sz "combined_clusters/${pep}/${pep}_cluster_sizes.xvg" \
              -clid "combined_clusters/${pep}/${pep}_clustid.xvg" \
              -dist "combined_clusters/${pep}/${pep}_rmsd.xvg"
done