#!/bin/bash
#
# Copyright (C) 2001 Quantum ESPRESSO
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

if [[ "$QE_USE_MPI" != "" ]]; then
  export PARA_PREFIX="mpirun -np $QE_USE_MPI"
  export PARA_SUFFIX="-npool $QE_USE_MPI"
else
  unset PARA_PREFIX
  unset PARA_SUFFIX
fi

echo $0" "$@
# First SCF, no k pools
if [[ "$1" == "1" ]]
then
  echo "Running PW ..."
  #echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
# First SCF, with k pools
elif [[ "$1" == "2" ]]
then
  echo "Running PW ..."
  #echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x < $2 > $3 2> $4"  
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
# Calc of U, no k pools
elif [[ "$1" == "3" ]]
then
  echo "Running HP ..."
  #echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/hp.x < $2 > $3 2> $4"  
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/hp.x < $2 > $3 2> $4
  cp *.Hubbard_parameters.dat $3
  if [[ -e CRASH ]]
  then
    cat $3
  fi
# Calc of U, with k pools
elif [[ "$1" == "4" ]]
then
  echo "Running HP ..."
  #echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/hp.x < $2 > $3 2> $4"  
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/hp.x ${PARA_SUFFIX} < $2 > $3 2> $4
  cp *.Hubbard_parameters.dat $3
  if [[ -e CRASH ]]
  then
    cat $3
  fi
# cRPA (hp.x method='cRPA'), no k pools (band basis)
elif [[ "$1" == "5" ]]
then
  echo "Running HP (cRPA) ..."
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/hp.x < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
# cRPA (hp.x method='cRPA'), with k pools
elif [[ "$1" == "6" ]]
then
  echo "Running HP (cRPA) ..."
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/hp.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
# pw2wannier90 with wan_mode='library': Wannier90 runs inside pw2wannier90.x,
# writes the .chk read by the cRPA run and reports to .wout, which is appended
# to the output that gets compared
elif [[ "$1" == "7" ]]
then
  echo "Running pw2wannier90 (Wannier90 library) ..."
  # a failed run must not leave the previous .wout / .chk to be compared and read
  rm -f *.wout *.chk
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw2wannier90.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
  cat *.wout >> $3 2> /dev/null
fi

#rm -f input_tmp.in
