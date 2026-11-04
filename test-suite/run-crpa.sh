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
# Driver for the CRPA (crpa.x) test-suite category.
#
#   $1 = step token (see below)   $2 = input file
#   $3 = stdout capture           $4 = stderr capture
#
# Steps:
#   1  pw.x            SCF
#   2  wannier90.x -pp (nnkp generation)
#   3  pw2wannier90.x  (amn/mmn/eig)
#   4  wannier90.x     (localisation)
#   5  crpa.x          Sternheimer/CG solver
#   6  crpa.x          non-Hermitian Lanczos solver
#
# For steps 5 and 6 a block of "CRPA-EXTRACT <key> <value>" lines is appended
# to the stdout capture ($3) so that extract-crpa.sh (which only sees $3) can
# compare the quantities that live in crpa.x's output files rather than stdout:
# the on-site Coulomb matrix element at q=Gamma (vbare/vscrd), the star-weighted
# Wannier-traced on-site U (u_bare/u_scrd, from the crpa_reduce_q irreducible q),
# and Re/Im of the macroscopic dielectric constant (eps_cg/eps_cg_im,
# eps_ln/eps_ln_im). Each solver is checked against its own benchmark only:
# CG-vs-Lanczos differences are not compared, so the Lanczos chains are short.

# Pool count is hardcoded to 1 for reproducibility across launch layouts.
# Change "-npool 1" here if you want crpa.x / pw.x to run with k-point pools.
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

# helper: star-weighted, Wannier-traced on-site Coulomb U
#         U = sum_iq w_iq * (1/N_w) sum_i U_ii(iq)
#   where w_iq are the star weights in <prefix>.q (col 4, after 2 header lines)
#   and U_ii(iq) is the diagonal on-site element (iw1=jw1=iw2=jw2=i, all R=0) of
#   <prefix>.<suffix><iq>. In crpa_mode='dHP' with crpa_reduce_q the trace is
#   star-invariant, so this equals the full-grid average while computing only
#   the irreducible q.
#   $1 = file prefix (crpa_cg | crpa_ln)  $2 = element suffix (bare|scrd)
#   $3 = Re(V) column (11 for *.bare, 13 for *.scrd)
crpa_wtrace () {
  local pref=$1 suf=$2 col=$3
  [[ -f "${pref}.q" ]] || return 1
  awk -v pref="$pref" -v suf="$suf" -v col="$col" '
    NR==1 { next }            # nq1 nq2 nq3
    NR==2 { next }            # number of (reduced) q
    {
      iq++; w=$4
      f = pref "." suf iq
      tr=0; n=0; w0first=""
      while ((getline line < f) > 0) {
        ns=split(line, a, " ")
        if (a[1] ~ /^#/ || ns < col) continue
        if (a[1]==a[2] && a[6]==a[7] && a[1]==a[6] &&
            a[3]==0 && a[4]==0 && a[5]==0 && a[8]==0 && a[9]==0 && a[10]==0) {
          if (suf == "scrd") {                      # one frequency only
            if (w0first == "") w0first = a[11]
            if (a[11] != w0first) continue
          }
          tr += a[col]; n++
        }
      }
      close(f)
      if (n>0) { acc += w * tr / n; wsum += w }
    }
    END { if (wsum>0) printf "%.9f\n", acc }
  ' "${pref}.q"
}

# helper: modulus of one element of the R3-resolved W (crpa_w_rspace), at
#         the lattice vector R3 = (1,0,0) modulo the q grid (the written label
#         is its Wigner-Seitz image, the value is the same). The modulus does
#         not depend on the phases of the Wannier functions, which wannier90's
#         conv_noise makes run-dependent and which rotate Re into Im for
#         off-diagonal pairs.
#   $1 = file prefix (crpa_cg | crpa_ln)  $2 = element file name (no .dat)
#   $3 = column of the real part (7 = V, 9 = W)
crpa_wr () {
  local f=$1.dir-WR/$2.dat
  [[ -f "$f" && -f "$1.q" ]] || return 1
  read -r n1 n2 n3 < "$1.q"
  awk -v c="$3" -v n1="$n1" -v n2="$n2" -v n3="$n3" '
    function m(x, n) { return ((x % n) + n) % n }
    $1 !~ /^#/ && m($1, n1) == 1 % n1 && m($2, n2) == 0 && m($3, n3) == 0 {
      printf "%.12f\n", sqrt($c * $c + $(c+1) * $(c+1)); exit }
  ' "$f"
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
  ${ESPRESSO_BUILD}/bin/wannier90.x -pp $2
  if [[ -e Si.werr ]]; then cat Si.werr; fi
  seed=${2%.win}
  nb=$(grep -A1 'begin nnkpts' ${seed}.nnkp | tail -1 | awk '{print $1}')
  [[ -n "$nb" ]] && emit nnkpts "$nb" "$3"
fi

if [[ "$1" == "3" ]]
then
  echo "Running pw2wannier90 ..."
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw2wannier90.x < $2 > $3 2> $4
  if [[ -e CRASH ]]; then cat $3; fi
  eig=$(head -1 Si.eig | awk '{print $3}')
  [[ -n "$eig" ]] && emit eig1 "$eig" "$3"
fi

if [[ "$1" == "4" ]]
then
  echo "Running Wannier90 ..."
  ${ESPRESSO_BUILD}/bin/wannier90.x $2
  cp *.wout $3 2> /dev/null
  if [[ -e Si.werr ]]; then cat Si.werr; fi
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

  # ---- dielectric-function mode -------------------------------------------
  if grep -qi 'crpa_epsil_only' "$2" && [[ -f epsilon_w.dat ]]; then
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

  # ---- screened-Coulomb mode --------------------------------------------
  if [[ "$1" == "5" && -f crpa_cg.bare1 ]]; then
    # on-site element at q#1 (= Gamma, always irreducible)
    vb=$(crpa_diag crpa_cg.bare1 11)
    vs=$(crpa_diag crpa_cg.scrd1 13)
    [[ -n "$vb" ]] && emit vbare    "$vb" "$3"
    [[ -n "$vs" ]] && emit vscrd_cg "$vs" "$3"
    # star-weighted Wannier-traced on-site U (the physical q-averaged value)
    ub=$(crpa_wtrace crpa_cg bare 11)
    uc=$(crpa_wtrace crpa_cg scrd 13)
    [[ -n "$ub" ]] && emit u_bare    "$ub" "$3"
    [[ -n "$uc" ]] && emit u_scrd_cg "$uc" "$3"
  fi
  if [[ "$1" == "6" && -f crpa_ln.bare1 ]]; then
    vb=$(crpa_diag crpa_ln.bare1 11)
    vs=$(crpa_diag crpa_ln.scrd1 13)
    [[ -n "$vb" ]] && emit vbare_ln "$vb" "$3"
    [[ -n "$vs" ]] && emit vscrd_ln "$vs" "$3"
    ub=$(crpa_wtrace crpa_ln bare 11)
    ul=$(crpa_wtrace crpa_ln scrd 13)
    [[ -n "$ub" ]] && emit u_bare_ln "$ub" "$3"
    [[ -n "$ul" ]] && emit u_scrd_ln "$ul" "$3"
  fi

  # ---- R3-resolved W (crpa_w_rspace) --------------------------------------
  # second pair translated by R3 = (1,0,0): a density-density element (bonds 1
  # and 2) and one with the off-diagonal matrix-element pair (1,2,0)
  if [[ "$1" == "5" && -d crpa_cg.dir-WR ]]; then
    wv=$(crpa_wr crpa_cg W_1_1_0_0_0_2_2_0_0_0 7)
    ww=$(crpa_wr crpa_cg W_1_1_0_0_0_2_2_0_0_0 9)
    [[ -n "$wv" ]] && emit wr_vbare    "$wv" "$3"
    [[ -n "$ww" ]] && emit wr_wscrd_cg "$ww" "$3"
    ov=$(crpa_wr crpa_cg W_1_2_0_0_0_1_1_0_0_0 7)
    ow=$(crpa_wr crpa_cg W_1_2_0_0_0_1_1_0_0_0 9)
    [[ -n "$ov" ]] && emit wr_off_vbare    "$ov" "$3"
    [[ -n "$ow" ]] && emit wr_off_wscrd_cg "$ow" "$3"
  fi
fi
