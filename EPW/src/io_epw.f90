  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  !
  ! Copyright (C) 2002-2013 Quantum ESPRESSO group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  ! Code adapted from Modules/io_global.f90 - Quantum-ESPRESSO group
  !
  ! SP: I'm missing some that depend on QE like iuwfc, lrwfc .. Should ideally
  !     be included 
  !
  !----------------------------------------------------------------------------
  MODULE io_epw
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  PUBLIC :: lambda_phself, linewidth_phself, linewidth_elself, iospectral, &
            iua2ffil, iudosfil, iufillambda, iuqdos, iufe, iufilker, &
            iufilgap, iospectral_sup, iua2ftrfil, iufilgapFS, iufillambdaFS, &
            iospectral_cum, iuwanep, iuwane, iunukk, iudvscf, iuqpeig, iures
  PUBLIC :: epwdata, iundmedata, iunvmedata, iunksdata, iudyn, iukgmap, iuepb,&
            iufilfreq, iufilegnv, iufileph, iufilkqmap, &
            iufilikmap, iueig, iunepmatwp, iunepmatwe, iunkf, iunqf, &
            iufileig, iukmap, crystal, iunifc, iunimem, iunepmatwp2
  PUBLIC :: iuwinfil, iun_plot, iuukk, iuprojfil, iudecayH, iudecayP, &
            iudecaydyn, iudecayv, iummn, iubvec
  PUBLIC :: iufilsigma, iufilseebeck, iufilkappael, iufilkappa, iufilscatt_rate,&
            iufilFi_all, iufilsigma_all, iufiltau_all, iuindabs
  PUBLIC :: iunsparseq, iunsparsek, iunsparsei, iunsparsej, iunsparset, iunselecq, &
            iunsparseqcb, iunsparsekcb, iunsparseicb, iunsparsejcb, iunsparsetcb, &
            iunrestart, iufilibtev_sup, iunepmat, iunepmatcb

  !
  ! Output of physically relevant quantities (60-100)
  !    
  INTEGER :: iuqpeig         = 59  ! Reading quasi-particle eigenenergies from file  
  INTEGER :: lambda_phself   = 60  ! Lambda factor of the phonon self-energy
                                   ! [lambda.phself] 
  INTEGER :: linewidth_phself= 61  ! Imaginary part of the phonon self-energy
                                   ! [linewidth.phself]
  INTEGER :: linewidth_elself= 62  ! Imaginary part of the electron self-energy
                                   ! [linewidth.elself]
  INTEGER :: iospectral      = 63  ! Electronic spectral function [specfun.elself]
  INTEGER :: iospectral_sup  = 64  ! Support data for the spectral function 
                                   ! [specfun_sup.elself]
  INTEGER :: iua2ffil        = 65  ! Eliashberg a2f function [.a2f]
  INTEGER :: iudosfil        = 66  ! Phonon density of state [.phdos]
  INTEGER :: iufillambda     = 67  ! Electron-phonon coupling strength lambda
                                   ! [.lambda_X]
  INTEGER :: iuqdos          = 68  ! Quasiparticle density of states in the
                                   ! superconducting state [.qdos]
  INTEGER :: iufe            = 69  ! Free energy [.fe]
  INTEGER :: iufilker        = 70  ! Eliashberg kernel [.ker]
  INTEGER :: iufilgap        = 71  ! Eliashberg superconducting gap [.gapr]
  INTEGER :: iua2ftrfil      = 72  ! Eliashberg transport a2f function [.a2f_tr]
  INTEGER :: iufilgapFS      = 73  ! Eliashberg superconducting gap on FS with k-points  
  INTEGER :: iufillambdaFS   = 74  ! Electron-phonon coupling strength on FS with k-points
  INTEGER :: iospectral_cum  = 75  ! Electronic spectral function with the cumulant method
                                   ! [specfun_cum##.elself]
!DBSP : iukgmap was 96. Should be the same as set_kplusq.f90. 
  INTEGER :: iunukk          = 77  ! Unit with rotation matrix U(k) from wannier code
  INTEGER :: iures           = 78  ! Resistivity in metals using Ziman formula [.res]
  INTEGER :: iudvscf         = 80  ! Unit for the dvscf_q file
  INTEGER :: iudyn           = 81  ! Unit for the dynamical matrix file
  INTEGER :: iufilkqmap      = 82  ! Map of k+q
  INTEGER :: iukgmap         = 96  ! Map of folding G-vector indexes [.kgmap]
  INTEGER :: iuwanep         = 97  ! Spatial decay of e-p matrix elements in wannier basis 
                                   ! Electrons + phonons [epmat_wanep]
  INTEGER :: iuwane          = 98  ! Spatial decay of matrix elements in Wannier basis    
                                   ! [.epwane]  
  !
  ! Output of quantity for restarting purposes (101-200)
  ! Note that 100-102 are reserved Cray unit and cannot be used. 
  ! 
  INTEGER :: iunvmedata      = 103  ! Velocity matrix in wannier basis [vmedata.fmt]
  INTEGER :: iunksdata       = 104  ! Hamiltonian in wannier basis
  INTEGER :: iuepb           = 105  ! Electron-phonon matrix in Bloch 
                                    ! representation [.epb]
  INTEGER :: iufilfreq       = 108  ! Phonon frequency from a previous epw run
                                    ! [.freq]
  INTEGER :: iufilegnv       = 109  ! Eigenvalues from a previous epw run [.egnv]
  INTEGER :: iufileph        = 110  ! Electron-phonon matrix element in the
                                    ! Bloch representation on the fine mesh
                                    ! [.ephmat]
  INTEGER :: iufilikmap      = 112  ! Index of k+(sign)q on the irreducible k-mesh
                                    ! [.ikmap]
!  INTEGER :: iuetf           = 113  ! Interpolated hamiltonian eigenvalues
  INTEGER :: iueig           = 114  ! Temporary eig for interpolation    
  INTEGER :: iunepmatwp      = 115  ! The unit with the e-ph matrix in Wannier-Wannier representation
  INTEGER :: iunepmatwe      = 116  ! The unit with the e-ph matrix in Wannier-Bloch representation
  INTEGER :: iunkf           = 117  ! The unit with the fine k-point mesh in crystal coord.
  INTEGER :: iunqf           = 118  ! The unit with the fine q-point mesh in crystal coord. 
  INTEGER :: iufileig        = 119  ! The unit with eigenenergies [band.eig]
  INTEGER :: iukmap          = 120  ! Unit for the k-point map generation
  INTEGER :: crystal         = 121  ! Unit for crystal data
  INTEGER :: iunifc          = 122  ! Unit for the IFC produced by q2r.x
  INTEGER :: iunimem         = 123  ! Unit for reading memory information from the system status file
  INTEGER :: epwdata         = 124  ! EPW data [epwdata.fmt] 
  INTEGER :: iundmedata      = 125  ! Dipole matrix in wannier basis [dmedata.fmt]
  INTEGER :: iunepmatwp2     = 126  ! Opening the epmatwp file
  INTEGER :: iufilibtev_sup  = 127  ! Files containing velocities for IBTE
  INTEGER :: iunsparseq      = 128  ! Q-mapping for IBTE
  INTEGER :: iunsparsek      = 129  ! K-mapping for IBTE
  INTEGER :: iunsparsei      = 130  ! i band mapping for IBTE
  INTEGER :: iunsparsej      = 131  ! j band mapping for IBTE
  INTEGER :: iunsparset      = 132  ! temperature mapping for IBTE
  INTEGER :: iunsparseqcb    = 133  ! Q-mapping for IBTE of conduction band
  INTEGER :: iunsparsekcb    = 134  ! K-mapping for IBTE for conduction band
  INTEGER :: iunsparseicb    = 135  ! i band mapping for IBTE for conduction band
  INTEGER :: iunsparsejcb    = 136  ! j band mapping for IBTE for conduction band
  INTEGER :: iunsparsetcb    = 137  ! temperature mapping for IBTE for conduction band
  INTEGER :: iunselecq       = 138  ! file containing q-point inside the fsthick windows
  INTEGER :: iunrestart      = 139  ! restart file during writing of IBTE scattering elements
  INTEGER :: iunepmat        = 140  ! Opening the epmatkq files
  INTEGER :: iunepmatcb      = 141  ! Opening the epmatkqcb file

  !
  ! Output quantites related to Wannier (201-250)
  !  
  INTEGER :: iuwinfil        = 201  ! Wannier projectors and other quantities
! SP : Not used for now but could be in the future. Would require the amn as well.
  INTEGER :: iummn           = 202  ! Overlap of the cell periodic part of the Bloch 
                                    ! states <u_nmk|u_nk+b>
  INTEGER :: iun_plot        = 203  ! UNK file (needed by Wannier90 for plotting the 
                                    ! real space Wannier functions)
  INTEGER :: iuukk           = 204  ! Final ukk rotation matrix (the big U!)
  INTEGER :: iuprojfil       = 205  ! Unit for projector [.projw90]  
  INTEGER :: iudecayH        = 206  ! Hamiltonian decay in real space
  INTEGER :: iudecayP        = 207  ! Dipole decay in real space
  INTEGER :: iudecaydyn      = 208  ! Dynamical matrix decay in real space
  INTEGER :: iudecayv        = 209  ! Velocity matrix decay in real space
  INTEGER :: iubvec          = 206  ! b-vectors and their weight wb
  !
  ! Output quantites related to transport (251-300)
  INTEGER :: iufilsigma      = 251 ! Electrical conductivity
  INTEGER :: iufilseebeck    = 252 ! Seebeck coefficient
  INTEGER :: iufilkappael    = 253 ! Electronic contribution to thermal conductivity
  INTEGER :: iufilkappa      = 254 ! Electronic contribution to thermal conductivity
  INTEGER :: iufilscatt_rate = 255 ! scattering rate
  INTEGER :: iufilFi_all     = 256 ! Fi_all file to retart at X iteration
  INTEGER :: iufilsigma_all  = 257 ! Sigmar_all and Sigmai_all file to retart an interpolation
  INTEGER :: iufiltau_all    = 258 ! inv_tau_all file to retart an interpolation
  !
  ! Output quantities related to Indirect absorption (301-325)
  INTEGER :: iuindabs        = 301 ! Indirect absorption data
  ! 
END MODULE io_epw
