#!/bin/bash
#
# Copyright (C) 2001 Quantum ESPRESSO
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
#
# Driver for the CRPA_METAL (crpa.x on metals, SrVO3) test-suite category.
# The insulating (Si) suites use run-crpa.sh.
#
#   $1 = step token (see below)   $2 = input file
#   $3 = stdout capture           $4 = stderr capture
#
# Steps:
#   1  pw.x            SCF or NSCF
#   2  wannier90.x -pp (nnkp generation)
#   3  pw2wannier90.x  (amn/mmn/eig)
#   4  wannier90.x     (localisation)
#   5  crpa.x          Sternheimer/CG solver
#   6  crpa.x          non-Hermitian Lanczos solver
#
# For steps 5 and 6 a block of "CRPA-EXTRACT <key> <value>" lines is appended
# to the stdout capture ($3) so that extract-crpa-metal.sh (which only sees $3)
# can compare quantities that live in crpa.x's output files:
#   screened Coulomb: on-site U_11 at q#1 (vbare, vscrd_cg) and at q#2
#                     (vbare_q2, vscrd_cg_q2, vscrd_ln_q2);
#   dielectric:       Re/Im eps_11 from CG and Lanczos (eps_cg, eps_cg_im,
#                     eps_ln, eps_ln_im) and the Drude omega_p^2 (wp2).
# CG-vs-Lanczos differences are not emitted: a short metallic Lanczos chain is
# far from converged at small eta, so each solver is checked against its own
# reference.

# Pool count is hardcoded to 1 for reproducibility across launch layouts.
if [[ "$QE_USE_MPI" != "" ]]; then
  export PARA_PREFIX="mpirun -np $QE_USE_MPI"
  export PARA_SUFFIX="-npool 1"
else
  unset PARA_PREFIX
  unset PARA_SUFFIX
fi

# ----------------------------------------------------------------------------
# helper: on-site diagonal element  iw1=jw1=iw2=jw2=1, all R=0
#   $1 = file   $2 = column index of Re(V) (11 for *.bare, 13 for *.scrd)
crpa_diag () {
  local f=$1 col=$2
  [[ -f "$f" ]] || return 1
  awk -v c="$col" \
    '$1==1 && $2==1 && $3==0 && $4==0 && $5==0 && \
     $6==1 && $7==1 && $8==0 && $9==0 && $10==0 { print $c; exit }' "$f"
}

# helper: Re(e11) / Im(e11) from the first data row of a written epsilon file
eps_re11 () {
  [[ -f "$1" ]] || return 1
  grep -v '^#' "$1" | awk 'NF>=3 { print $3; exit }'
}
eps_im11 () {
  [[ -f "$1" ]] || return 1
  grep -v '^#' "$1" | awk 'NF>=4 { print $4; exit }'
}

emit () { echo "CRPA-EXTRACT $1 $2" >> "$3"; }
# ----------------------------------------------------------------------------

if [[ "$1" == "1" ]]
then
  echo "Running PW ..."
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]; then cat $3; fi
fi

if [[ "$1" == "2" ]]
then
  echo "Running Wannier90 PP ..." > $3
  seed=${2%.win}
  ${ESPRESSO_BUILD}/bin/wannier90.x -pp $2
  if [[ -e ${seed}.werr ]]; then cat ${seed}.werr; fi
  nb=$(grep -A1 'begin nnkpts' ${seed}.nnkp | tail -1 | awk '{print $1}')
  [[ -n "$nb" ]] && emit nnkpts "$nb" "$3"
fi

if [[ "$1" == "3" ]]
then
  echo "Running pw2wannier90 ..."
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw2wannier90.x < $2 > $3 2> $4
  if [[ -e CRASH ]]; then cat $3; fi
  # seedname from the pw2wannier90 namelist (quotes and blanks stripped)
  seed=$(grep -i 'seedname' $2 | head -1 | sed "s/.*=[[:space:]]*['\"]\{0,1\}\([^'\"[:space:]]*\).*/\1/")
  eig=$(head -1 ${seed}.eig | awk '{print $3}')
  [[ -n "$eig" ]] && emit eig1 "$eig" "$3"
fi

if [[ "$1" == "4" ]]
then
  echo "Running Wannier90 ..."
  ${ESPRESSO_BUILD}/bin/wannier90.x $2
  cp *.wout $3 2> /dev/null
  if [[ -e ${2%.win}.werr ]]; then cat ${2%.win}.werr; fi
fi

if [[ "$1" == "5" || "$1" == "6" ]]
then
  if [[ "$1" == "5" ]]; then
    echo "Running CRPA (Sternheimer/CG) ..."
  else
    echo "Running CRPA (Lanczos) ..."
  fi
  # Clean folder to avoid wrongfull restarts.
  rm -rf out/_crpa0 out/_epsw0 out/_ph0 diel_lanczos_chain.dat
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/crpa.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]; then cat $3; fi

  # ---- dielectric-function mode (with the Drude term) ---------------------
  if grep -qi 'crpa_epsil_only' "$2" && [[ -f epsilon_w.dat ]]; then
    wp=$(grep 'omega_p^2 (Ry^2)' "$3" | head -1 | awk '{print $4}')
    [[ -n "$wp" ]] && emit wp2 "$wp" "$3"
    vr=$(eps_re11 epsilon_w.dat)
    vi=$(eps_im11 epsilon_w.dat)
    if [[ "$1" == "5" ]]; then
      [[ -n "$vr" ]] && emit eps_cg    "$vr" "$3"
      [[ -n "$vi" ]] && emit eps_cg_im "$vi" "$3"
    else
      [[ -n "$vr" ]] && emit eps_ln    "$vr" "$3"
      [[ -n "$vi" ]] && emit eps_ln_im "$vi" "$3"
    fi
  fi

  # ---- screened-Coulomb mode: on-site U_11 at q#1 (if run) and q#2 --------
  if [[ "$1" == "5" ]]; then
    if [[ -f crpa_cg.bare1 ]]; then
      vb=$(crpa_diag crpa_cg.bare1 11); vs=$(crpa_diag crpa_cg.scrd1 13)
      [[ -n "$vb" ]] && emit vbare    "$vb" "$3"
      [[ -n "$vs" ]] && emit vscrd_cg "$vs" "$3"
    fi
    if [[ -f crpa_cg.bare2 ]]; then
      vb=$(crpa_diag crpa_cg.bare2 11); vs=$(crpa_diag crpa_cg.scrd2 13)
      [[ -n "$vb" ]] && emit vbare_q2    "$vb" "$3"
      [[ -n "$vs" ]] && emit vscrd_cg_q2 "$vs" "$3"
    fi
  else
    if [[ -f crpa_ln.scrd2 ]]; then
      vs=$(crpa_diag crpa_ln.scrd2 13)
      [[ -n "$vs" ]] && emit vscrd_ln_q2 "$vs" "$3"
    fi
  fi
fi
