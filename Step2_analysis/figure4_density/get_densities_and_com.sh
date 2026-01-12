#!/bin/bash
set -euo pipefail

num_atoms=580

# --- Car9 tag selection (script 1 logic) ---
tag_start=380
num_tag_atoms=194
# By default make_ndx has a total of 17 groups
group_num=17
# --- ELP+linker selection (script 2 logic) ---
elp_start=1
num_elp_atoms=378

build_union_expr () {
  local start="$1"
  local count="$2"
  local expr=""
  for i in {0..8}; do
    local lo=$(( start + i*num_atoms ))
    local hi=$(( start + count ))
    hi=$(( hi + i*num_atoms ))
    # append with OR between terms, no trailing '|'
    if [[ -z "$expr" ]]; then
      expr="a ${lo}-${hi}"
    else
      expr="${expr} | a ${lo}-${hi}"
    fi
  done
  echo "$expr"
}

tag_expr="$(build_union_expr "$tag_start" "$num_tag_atoms")"
elp_expr="$(build_union_expr "$elp_start" "$num_elp_atoms")"

# One make_ndx call creates both groups and names them reliably via "name 0"
input=""
input+="${tag_expr}\n"
input+="name  ${group_num} Car9_tag\n" # The tag should be group 18
input+="${elp_expr}\n"
group_num=$((group_num+1)) # ELP will group 19
input+="name ${group_num} ELP\n"
input+="q\n"

echo -e "$input" | gmx make_ndx -f npt2.gro -o tags_and_elp.ndx
echo "$group_num" | gmx traj -s npt2.tpr -f npt2.xtc -n tags_and_elp.ndx -com -ox com_ELP.xvg

echo "Car9_tag"| gmx density -f npt2.xtc -s npt2.tpr -o density.xvg -n tags_and_elp.ndx 
