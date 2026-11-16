# Copyright (C) 2001 Quantum ESPRESSO
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
#
# Output parser for the CRPA (crpa.x) test-suite category.
# Reads the stdout capture ($1) and prints "label / value" line pairs.
# For the crpa.x output the CRPA-EXTRACT helper in run-crpa.sh. 
#
fname=$1

# --- pw.x SCF ---
e1=`grep ^! $fname | tail -1 | awk '{printf "%12.6f\n", $5}'`
n1=`grep 'convergence has' $fname | tail -1 | awk '{print $6}'`

# --- wannier90 localisation ---
omega=`grep "Omega Total" $fname | tail -1 | awk '{print $7}'`

# --- crpa.x digested quantities ---
for key in nnkpts eig1 vbare vbare_ln vscrd_cg vscrd_ln \
           u_bare u_bare_ln u_scrd_cg u_scrd_ln \
           eps_cg eps_cg_im eps_ln eps_ln_im \
           wr_vbare wr_wscrd_cg wr_off_vbare wr_off_wscrd_cg ; do
  val=`grep "CRPA-EXTRACT ${key} " $fname | tail -1 | awk '{print $3}'`
  eval "$key=\"$val\""
done

if test "$e1" != ""; then echo e1; echo $e1; fi
if test "$n1" != ""; then echo n1; echo $n1; fi
if test "$omega" != ""; then echo omega; echo $omega; fi
if test "$nnkpts" != ""; then echo nnkpts; echo $nnkpts; fi
if test "$eig1" != ""; then echo eig1; echo $eig1; fi
if test "$vbare" != ""; then echo vbare; echo $vbare; fi
if test "$vbare_ln" != ""; then echo vbare_ln; echo $vbare_ln; fi
if test "$vscrd_cg" != ""; then echo vscrd_cg; echo $vscrd_cg; fi
if test "$vscrd_ln" != ""; then echo vscrd_ln; echo $vscrd_ln; fi
if test "$u_bare" != ""; then echo u_bare; echo $u_bare; fi
if test "$u_bare_ln" != ""; then echo u_bare_ln; echo $u_bare_ln; fi
if test "$u_scrd_cg" != ""; then echo u_scrd_cg; echo $u_scrd_cg; fi
if test "$u_scrd_ln" != ""; then echo u_scrd_ln; echo $u_scrd_ln; fi
if test "$eps_cg" != ""; then echo eps_cg; echo $eps_cg; fi
if test "$eps_cg_im" != ""; then echo eps_cg_im; echo $eps_cg_im; fi
if test "$eps_ln" != ""; then echo eps_ln; echo $eps_ln; fi
if test "$eps_ln_im" != ""; then echo eps_ln_im; echo $eps_ln_im; fi
if test "$wr_vbare" != ""; then echo wr_vbare; echo $wr_vbare; fi
if test "$wr_wscrd_cg" != ""; then echo wr_wscrd_cg; echo $wr_wscrd_cg; fi
if test "$wr_off_vbare" != ""; then echo wr_off_vbare; echo $wr_off_vbare; fi
if test "$wr_off_wscrd_cg" != ""; then echo wr_off_wscrd_cg; echo $wr_off_wscrd_cg; fi
