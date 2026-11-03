#!/bin/bash
#
# Copyright (C) 2020 Quantum ESPRESSO Foundation
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

if [[ "$QE_USE_MPI" != "" ]]; then
  export PARA_PREFIX="mpirun -np $QE_USE_MPI"
  export PARA_SUFFIX=" "
else
  unset PARA_PREFIX
  unset PARA_SUFFIX
fi

if [[ "$1" == "1" ]]
then
  # echo "Running PW ..."
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi  
elif [[ "$1" == "2" ]]
then
  # echo "Running PP ..."
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/ppacf.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "3" ]]
then
  # echo "Running pw2wannier90 ..."
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw2wannier90.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
  # In library mode Wannier90 runs inside pw2wannier90.x and reports to .wout,
  # so append it to the output that gets compared.
  cat *.wout >> $3 2> /dev/null
elif [[ "$1" == "4" ]]
then
  # echo "Running pw2wannier90 with k-point pools ..."
  # One pool per process. Without QE_USE_MPI this is the same as mode 3.
  if [[ "$QE_USE_MPI" != "" ]]; then
    export PARA_SUFFIX="-nk $QE_USE_MPI"
  fi
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw2wannier90.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
  cat *.wout >> $3 2> /dev/null
elif [[ "$1" == "5" ]]
then
  # echo "Running pw2wannier90 with two processes per k-point pool ..."
  # Odd process counts and runs without QE_USE_MPI use a single pool.
  if [[ "$QE_USE_MPI" != "" ]] && (( QE_USE_MPI % 2 == 0 )); then
    export PARA_SUFFIX="-nk $(( QE_USE_MPI / 2 ))"
  fi
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw2wannier90.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
  cat *.wout >> $3 2> /dev/null
elif [[ "$1" == "6" ]]
then
  # echo "Running pw2wannier90 decompose ..."
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw2wannier90.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
  # The decomposition is written to one .power file per Wannier function
  # and one _gc_ .coeff file per extra centre, so append the values, tagged,
  # to the output that gets compared.
  for f in *_[0-9][0-9][0-9][0-9][0-9].power
  do
    grep -v '^#' $f | awk '{print "POWER", $1}' >> $3
  done
  for f in *_gc_[0-9][0-9][0-9][0-9][0-9].coeff
  do
    grep -v '^#' $f | awk '{print "GCOEFF", $1}' >> $3
  done
elif [[ "$1" == "7" ]]
then
  # echo "Running Wannier90 -pp ..."
  ${ESPRESSO_BUILD}/bin/wannier90.x -pp $2
  cp *.nnkp $3
elif [[ "$1" == "8" ]]
then
  # echo "Running Wannier90 ..."
  ${ESPRESSO_BUILD}/bin/wannier90.x $2
  cp *.wout $3
elif [[ "$1" == "9" ]]
then
  # echo "Running pw2wannier90 standalone ..."
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw2wannier90.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
  # Standalone mode reports no physical quantity on stdout, so append the
  # eigenvalues it wrote, tagged, to the output that gets compared.
  awk '{print "EIG", $3}' *.eig >> $3
fi

