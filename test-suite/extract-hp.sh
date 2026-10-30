# Copyright (C) 2001 Quantum ESPRESSO
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
#
# SP: This can lead to issue if you reach the OS pipe buffer
#     You can increase the buffer in /proc/sys/fs/pipe-max-size
  
fname=$1
args=$(echo $fname | awk -F= '{print $NF}')

# cRPA Wannier setup (7 = pw2wannier90 in library mode, .wout appended): the
# final spreads of the one-shot (projection-only) Wannier functions
if test "$args" = "7"; then
	for o in 'Omega I ' 'Omega D ' 'Omega OD' 'Omega Total'; do
		val=`grep -- "$o" $fname | tail -1 | awk '{print $NF}'`
		lab=`echo "w90_$o" | tr -d ' ' | tr 'A-Z' 'a-z'`
		if test "$val" != ""; then echo $lab; echo $val; fi
	done
fi

##echo $fname > /home/sponce/program/espresso/test-suite/tmp.txt
###    echo $args >> /home/sponce/program/espresso/test-suite/tmp.txt

### if [[ "$args" == "1" ]]
##  then
# SCF
e1=`grep ! $fname | tail -1 | awk '{printf "%12.6f\n", $5}'`
n1=`grep 'convergence has' $fname | tail -1 | awk '{print $6}'`
f1=`grep "Total force" $fname | head -1 | awk '{printf "%8.4f\n", $4}'`
p1=`grep "P= " $fname | tail -1 | awk '{print $6}'`
### fi


# HP
u=`grep "Hubbard U (eV)" -A 1 $fname | tail -n 1 | awk '{print $8}'`
v1=`grep "Hubbard V (eV)" -A 2 $fname | tail -n 1 | awk '{print $6}'`
v2=`grep "Hubbard V (eV)" -A 3 $fname | tail -n 1 | awk '{print $6}'`


if test "$e1" != ""; then
	echo e1
	echo $e1
fi

if test "$n1" != ""; then
	echo n1
	echo $n1
fi

if test "$f1" != ""; then
	echo f1
	echo $f1
fi

if test "$p1" != ""; then
	echo p1
	echo $p1
fi

if test "$u" != ""; then
        echo u
        echo $u
fi

if test "$v1" != ""; then
        echo v1
        echo $v1
fi

if test "$v2" != ""; then
        echo v2
        echo $v2
fi

# cRPA (method='CRPA'): interaction summary printed by
# HP/src/crpa_postprocess_coulomb.f90. Each "label=" marker is a separate
# whitespace token, so its value is the next token (I4/F16.9 space-padded).
#
# Band basis: gauge-invariant pair-diagonal trace, one line per q point.
#   "     cRPA-U-TRACE q=   1  bare=     0.965796362  scrd=     0.668629601  bare_im=     0.000000000  scrd_im=     0.000000000"
awk '/cRPA-U-TRACE/ {
        q=""; b=""; s=""; bi=""; si="";
        for (i=1;i<=NF;i++) {
            if ($i=="q=")       q=$(i+1);
            if ($i=="bare=")    b=$(i+1);
            if ($i=="scrd=")    s=$(i+1);
            if ($i=="bare_im=") bi=$(i+1);
            if ($i=="scrd_im=") si=$(i+1);
        }
        print "utrace_bare_q" q; print b;
        print "utrace_scrd_q" q; print s;
        print "utrace_bare_im_q" q; print bi;
        print "utrace_scrd_im_q" q; print si;
     }' $fname

# Wannier pair lists (perturbation + matrix-element), printed once by
# crpa_prepare_perturbations_q. Emit the pair counts and each pair's WF-center
# distance (deterministic for projection-only Wannier); exercises the two
# distance cutoffs (wan_dist_thr / wan_dist_thr_mels).
#   "     cRPA-WFPAIRS pert npairs=     4"
#   "     cRPA-WFPAIR pert idx=    1  iw=   1  jw=   1  R=   0   0   0  dist=    0.00000000"
awk '/cRPA-WFPAIRS/ {
        n="";
        for (i=1;i<=NF;i++) if ($i=="npairs=") n=$(i+1);
        print "w90_npairs_" $2; print n;
     }' $fname
awk '/cRPA-WFPAIR / {
        tg=$2; idx=""; d="";
        for (i=1;i<=NF;i++) {
            if ($i=="idx=")  idx=$(i+1);
            if ($i=="dist=") d=$(i+1);
        }
        print "w90_wfp_" tg "_" idx "_dist"; print d;
     }' $fname

# Wannier basis: full on-site density-density matrix, one line per (q,i,j).
#   "     cRPA-V-IIJJ q=   1  i=   1  j=   1  bare=   2.520237691  scrd=   1.717595732  bare_im=   0.000000000  scrd_im=   0.000000000"
awk '/cRPA-V-IIJJ/ {
        q=""; ii=""; jj=""; b=""; s=""; bi=""; si="";
        for (i=1;i<=NF;i++) {
            if ($i=="q=")       q=$(i+1);
            if ($i=="i=")       ii=$(i+1);
            if ($i=="j=")       jj=$(i+1);
            if ($i=="bare=")    b=$(i+1);
            if ($i=="scrd=")    s=$(i+1);
            if ($i=="bare_im=") bi=$(i+1);
            if ($i=="scrd_im=") si=$(i+1);
        }
        print "viijj_bare_q" q "_" ii "_" jj; print b;
        print "viijj_scrd_q" q "_" ii "_" jj; print s;
        print "viijj_bare_im_q" q "_" ii "_" jj; print bi;
        print "viijj_scrd_im_q" q "_" ii "_" jj; print si;
     }' $fname
