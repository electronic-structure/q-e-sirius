# Copyright (C) 2001 Quantum ESPRESSO
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
#
# Output parser for the CRPA_METAL (crpa.x on metals) test-suite category.
# Reads the stdout capture ($1) and prints "label / value" line pairs.
# For crpa.x the values come from the CRPA-EXTRACT lines written by
# run-crpa-metal.sh.
#
fname=$1

# --- pw.x SCF / NSCF ---
e1=`grep ^! $fname | tail -1 | awk '{printf "%12.6f\n", $5}'`
n1=`grep 'convergence has' $fname | tail -1 | awk '{print $6}'`
# Fermi energy (eV): also the only checkable quantity of an NSCF run
ef=`grep 'the Fermi energy is' $fname | tail -1 | awk '{print $5}'`

# --- wannier90 localisation ---
omega=`grep "Omega Total" $fname | tail -1 | awk '{print $7}'`

# --- crpa.x digested quantities ---
for key in nnkpts eig1 vbare vscrd_cg vbare_q2 vscrd_cg_q2 vscrd_ln_q2 \
           eps_cg eps_cg_im eps_ln eps_ln_im wp2 ; do
  val=`grep "CRPA-EXTRACT ${key} " $fname | tail -1 | awk '{print $3}'`
  eval "$key=\"$val\""
done

for key in e1 n1 ef omega nnkpts eig1 vbare vscrd_cg vbare_q2 vscrd_cg_q2 \
           vscrd_ln_q2 eps_cg eps_cg_im eps_ln eps_ln_im wp2 ; do
  eval "val=\$$key"
  if test "$val" != ""; then echo $key; echo $val; fi
done
