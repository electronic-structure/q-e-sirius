! dftd3 program for computing the dispersion energy and forces from cart
! and atomic numbers as described in
!
! S. Grimme, J. Antony, S. Ehrlich and H. Krieg
! J. Chem. Phys, 132 (2010), 154104
!
! S. Grimme, S. Ehrlich and L. Goerigk, J. Comput. Chem, 32 (2011), 1456
! (for BJ-damping)
!
! Copyright (C) 2009 - 2011 Stefan Grimme, University of Muenster, Germany
!
! Repackaging of the original code without any change in the functionality:
!
! Copyright (C) 2016, Bálint Aradi
!
! r2r4 data for Ba corrected by Vahid Askarpour, July 2022
! removal of unused routines in QE by Paolo Giannozzi, March 2022
! MPI parallelization  added by Paolo Giannozzi, June 2021
! OpenACC acceleration added by Ivan Carnimeo,   June 2021
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 1, or (at your option)
! any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! For the GNU General Public License, see <http://www.gnu.org/licenses/>
!

module dftd3_core
  use dftd3_sizes
  use dftd3_common
  use dftd3_pars
  implicit none

  !! number of processors within an image
  INTEGER :: nproc_dftd3 = 1
  !! index of the processor within an image
  INTEGER :: me_dftd3 = 0
  ! intra image communicator
  INTEGER :: comm_dftd3 = 0

  ! atomic <r^2>/<r^4> values
  real(wp) r2r4(max_elem)

  ! r2r4 =sqrt(0.5*r2r4(i)*dfloat(i)**0.5 ) with i=elementnumber
  ! the large number of digits is just to keep the results consistent
  ! with older versions. They should not imply any higher accuracy than
  ! the old values
  data r2r4 / &
    &2.00734898_wp, 1.56637132_wp, 5.01986934_wp, 3.85379032_wp, 3.64446594_wp,&
    &3.10492822_wp, 2.71175247_wp, 2.59361680_wp, 2.38825250_wp, 2.21522516_wp,&
    &6.58585536_wp, 5.46295967_wp, 5.65216669_wp, 4.88284902_wp, 4.29727576_wp,&
    &4.04108902_wp, 3.72932356_wp, 3.44677275_wp, 7.97762753_wp, 7.07623947_wp,&
    &6.60844053_wp, 6.28791364_wp, 6.07728703_wp, 5.54643096_wp, 5.80491167_wp,&
    &5.58415602_wp, 5.41374528_wp, 5.28497229_wp, 5.22592821_wp, 5.09817141_wp,&
    &6.12149689_wp, 5.54083734_wp, 5.06696878_wp, 4.87005108_wp, 4.59089647_wp,&
    &4.31176304_wp, 9.55461698_wp, 8.67396077_wp, 7.97210197_wp, 7.43439917_wp,&
    &6.58711862_wp, 6.19536215_wp, 6.01517290_wp, 5.81623410_wp, 5.65710424_wp,&
    &5.52640661_wp, 5.44263305_wp, 5.58285373_wp, 7.02081898_wp, 6.46815523_wp,&
    &5.98089120_wp, 5.81686657_wp, 5.53321815_wp, 5.25477007_wp,11.02204549_wp,&
    &10.15679528_wp,9.35167836_wp, 9.06926079_wp, 8.97241155_wp, 8.90092807_wp,&
    &8.85984840_wp, 8.81736827_wp, 8.79317710_wp, 7.89969626_wp, 8.80588454_wp,&
    &8.42439218_wp, 8.54289262_wp, 8.47583370_wp, 8.45090888_wp, 8.47339339_wp,&
    &7.83525634_wp, 8.20702843_wp, 7.70559063_wp, 7.32755997_wp, 7.03887381_wp,&
    &6.68978720_wp, 6.05450052_wp, 5.88752022_wp, 5.70661499_wp, 5.78450695_wp,&
    &7.79780729_wp, 7.26443867_wp, 6.78151984_wp, 6.67883169_wp, 6.39024318_wp,&
    &6.09527958_wp,11.79156076_wp,11.10997644_wp, 9.51377795_wp, 8.67197068_wp,&
    &8.77140725_wp, 8.65402716_wp, 8.53923501_wp, 8.85024712_wp /

  ! PBE0/def2-QZVP atomic values
  ! data r2r4 /
  ! . 8.0589, 3.4698, 29.0974, 14.8517, 11.8799, 7.8715, 5.5588,
  ! . 4.7566, 3.8025, 3.1036, 26.1552, 17.2304, 17.7210, 12.7442,
  ! . 9.5361, 8.1652, 6.7463, 5.6004, 29.2012, 22.3934, 19.0598,
  ! . 16.8590, 15.4023, 12.5589, 13.4788, 12.2309, 11.2809, 10.5569,
  ! . 10.1428, 9.4907, 13.4606, 10.8544, 8.9386, 8.1350, 7.1251,
  ! . 6.1971, 30.0162, 24.4103, 20.3537, 17.4780, 13.5528, 11.8451,
  ! . 11.0355, 10.1997, 9.5414, 9.0061, 8.6417, 8.9975, 14.0834,
  ! . 11.8333, 10.0179, 9.3844, 8.4110, 7.5152, 32.7622, 27.5708,
  ! . 23.1671, 21.6003, 20.9615, 20.4562, 20.1010, 19.7475, 19.4828,
  ! . 15.6013, 19.2362, 17.4717, 17.8321, 17.4237, 17.1954, 17.1631,
  ! . 14.5716, 15.8758, 13.8989, 12.4834, 11.4421, 10.2671, 8.3549,
  ! . 7.8496, 7.3278, 7.4820, 13.5124, 11.6554, 10.0959, 9.7340,
  ! . 8.8584, 8.0125, 29.8135, 26.3157, 19.1885, 15.8542, 16.1305,
  ! . 15.6161, 15.1226, 16.1576 /


  ! scale r4/r2 values of the atoms by sqrt(Z)
  ! sqrt is also globally close to optimum
  ! together with the factor 1/2 this yield reasonable
  ! c8 for he, ne and ar. for larger Z, C8 becomes too large
  ! which effectively mimics higher R^n terms neglected due
  ! to stability reasons

  ! r2r4 =sqrt(0.5*r2r4(i)*dfloat(i)**0.5 ) with i=elementnumber
  ! the large number of digits is just to keep the results consistent
  ! with older versions. They should not imply any higher accuracy than
  ! the old values
  !!data r2r4 / &
  !!    & 2.00734898, 1.56637132, 5.01986934, 3.85379032, 3.64446594, &
  !!    & 3.10492822, 2.71175247, 2.59361680, 2.38825250, 2.21522516, &
  !!    & 6.58585536, 5.46295967, 5.65216669, 4.88284902, 4.29727576, &
  !!    & 4.04108902, 3.72932356, 3.44677275, 7.97762753, 7.07623947, &
  !!    & 6.60844053, 6.28791364, 6.07728703, 5.54643096, 5.80491167, &
  !!    & 5.58415602, 5.41374528, 5.28497229, 5.22592821, 5.09817141, &
  !!    & 6.12149689, 5.54083734, 5.06696878, 4.87005108, 4.59089647, &
  !!    & 4.31176304, 9.55461698, 8.67396077, 7.97210197, 7.43439917, &
  !!    & 6.58711862, 6.19536215, 6.01517290, 5.81623410, 5.65710424, &
  !!    & 5.52640661, 5.44263305, 5.58285373, 7.02081898, 6.46815523, &
  !!    & 5.98089120, 5.81686657, 5.53321815, 5.25477007, 11.02204549, &
  !!    &10.15679528, 9.35167836, 9.06926079, 8.97241155, 8.90092807, &
  !!    & 8.85984840, 8.81736827, 8.79317710, 7.89969626, 8.80588454, &
  !!    & 8.42439218, 8.54289262, 8.47583370, 8.45090888, 8.47339339, &
  !!    & 7.83525634, 8.20702843, 7.70559063, 7.32755997, 7.03887381, &
  !!    & 6.68978720, 6.05450052, 5.88752022, 5.70661499, 5.78450695, &
  !!    & 7.79780729, 7.26443867, 6.78151984, 6.67883169, 6.39024318, &
  !!    & 6.09527958, 11.79156076, 11.10997644, 9.51377795, 8.67197068, &
  !!    & 8.77140725, 8.65402716, 8.53923501, 8.85024712 /


  real(wp) rcov(max_elem)

  ! covalent radii
  ! covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
  ! values for metals decreased by 10 %
  ! data rcov/
  ! . 0.32, 0.46, 1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67
  ! ., 1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96, 1.76, 1.54
  ! ., 1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09
  ! ., 1.12, 1.09, 1.15, 1.10, 1.14, 1.17, 1.89, 1.67, 1.47, 1.39
  ! ., 1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, 1.28, 1.26
  ! ., 1.26, 1.23, 1.32, 1.31, 2.09, 1.76, 1.62, 1.47, 1.58, 1.57
  ! ., 1.56, 1.55, 1.51, 1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53
  ! ., 1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32
  ! ., 1.30, 1.30, 1.36, 1.31, 1.38, 1.42, 2.01, 1.81, 1.67, 1.58
  ! ., 1.52, 1.53, 1.54, 1.55 /

  ! these new data are scaled with k2=4./3. and converted a_0 via
  ! autoang=0.52917726d0
  !!data rcov/ &
  !!    & 0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865, &
  !!    & 1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527, &
  !!    & 3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820, &
  !!    & 2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730, &
  !!    & 3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923, &
  !!    & 2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188, &
  !!    & 2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349, &
  !!    & 2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216, &
  !!    & 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717, &
  !!    & 2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967, &
  !!    & 3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625, &
  !!    & 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657, &
  !!    & 3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833, &
  !!    & 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098, &
  !!    & 3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878, &
  !!    & 2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790, &
  !!    & 3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584, &
  !!    & 3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289, &
  !!    & 3.82984466, 3.85504098, 3.88023730, 3.90543362 /

  ! these new data are scaled with k2=4./3. and converted a_0 via
  ! autoang=0.52917726d0
  data rcov/ &
  & 0.80628308_wp, 1.15903197_wp, 3.02356173_wp, 2.36845659_wp, 1.94011865_wp, &
  & 1.88972601_wp, 1.78894056_wp, 1.58736983_wp, 1.61256616_wp, 1.68815527_wp, &
  & 3.52748848_wp, 3.14954334_wp, 2.84718717_wp, 2.62041997_wp, 2.77159820_wp, &
  & 2.57002732_wp, 2.49443835_wp, 2.41884923_wp, 4.43455700_wp, 3.88023730_wp, &
  & 3.35111422_wp, 3.07395437_wp, 3.04875805_wp, 2.77159820_wp, 2.69600923_wp, &
  & 2.62041997_wp, 2.51963467_wp, 2.49443835_wp, 2.54483100_wp, 2.74640188_wp, &
  & 2.82199085_wp, 2.74640188_wp, 2.89757982_wp, 2.77159820_wp, 2.87238349_wp, &
  & 2.94797246_wp, 4.76210950_wp, 4.20778980_wp, 3.70386304_wp, 3.50229216_wp, &
  & 3.32591790_wp, 3.12434702_wp, 2.89757982_wp, 2.84718717_wp, 2.84718717_wp, &
  & 2.72120556_wp, 2.89757982_wp, 3.09915070_wp, 3.22513231_wp, 3.17473967_wp, &
  & 3.17473967_wp, 3.09915070_wp, 3.32591790_wp, 3.30072128_wp, 5.26603625_wp, &
  & 4.43455700_wp, 4.08180818_wp, 3.70386304_wp, 3.98102289_wp, 3.95582657_wp, &
  & 3.93062995_wp, 3.90543362_wp, 3.80464833_wp, 3.82984466_wp, 3.80464833_wp, &
  & 3.77945201_wp, 3.75425569_wp, 3.75425569_wp, 3.72905937_wp, 3.85504098_wp, &
  & 3.67866672_wp, 3.45189952_wp, 3.30072128_wp, 3.09915070_wp, 2.97316878_wp, &
  & 2.92277614_wp, 2.79679452_wp, 2.82199085_wp, 2.84718717_wp, 3.32591790_wp, &
  & 3.27552496_wp, 3.27552496_wp, 3.42670319_wp, 3.30072128_wp, 3.47709584_wp, &
  & 3.57788113_wp, 5.06446567_wp, 4.56053862_wp, 4.20778980_wp, 3.98102289_wp, &
  & 3.82984466_wp, 3.85504098_wp, 3.88023730_wp, 3.90543362_wp /


contains


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! set parameters
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine setfuncpar(func,version,TZ,s6,rs6,s18,rs18,alp)
    integer version
    real(wp) s6,rs6,s18,alp,rs18
    character*(*) func
    logical TZ
    ! double hybrid values revised according to procedure in the GMTKN30 pap


    if(version.eq.6)then
    s6  =1.0d0
    alp =14.0d0
! BJ damping with parameters from ...
    select case (func)
      case ("b2-plyp")
           rs6 =0.486434_wp
           s18 =0.672820_wp
           rs18=3.656466_wp
     case ("b3-lyp") 
           rs6 =0.278672_wp
           s18 =1.466677_wp
           rs18=4.606311_wp
      case ("b97-d")
           rs6 =0.240184_wp
           s18 =1.206988_wp
           rs18=3.864426_wp
      case ("b-lyp")
           rs6 =0.448486_wp
           s18 =1.875007_wp
           rs18=3.610679_wp
      case ("b-p")
           rs6 =0.821850_wp
           s18 =3.140281_wp
           rs18=2.728151_wp
      case ("pbe")
           rs6 =0.012092_wp
           s18 =0.358940_wp
           rs18=5.938951_wp
      case ("pbe0")
           rs6 =0.007912_wp
           s18 =0.528823_wp
           rs18=6.162326_wp
      case ("lc-wpbe")
           rs6 =0.563761_wp
           s18 =0.906564_wp
           rs18=3.593680_wp
      case DEFAULT
           call errore('dft-d3', 'functional name unknown', 1)
    end select
    endif

    if(version.eq.5)then
    s6  =1.0d0
    alp =14.0d0
! zero damping with parameters from ...
    select case (func)
      case ("b2-plyp")
           rs6 =1.313134_wp
           s18 =0.717543_wp
           rs18=0.016035_wp
           s6  =0.640000_wp
      case ("b3-lyp")
           rs6 =1.338153_wp
           s18 =1.532981_wp
           rs18=0.013988_wp
      case ("b97-d")
           rs6 =1.151808_wp
           s18 =1.020078_wp
           rs18=0.035964_wp
      case ("b-lyp")
           rs6 =1.279637_wp
           s18 =1.841686_wp
           rs18=0.014370_wp
      case ("b-p")
           rs6 =1.233460_wp
           s18 =1.945174_wp
           rs18=0.000000
      case ("pbe")
           rs6 =2.340218_wp
           s18 =0.000000
           rs18=0.129434_wp
      case ("pbe0")
           rs6 =2.077949_wp
           s18 =0.000081_wp
           rs18=0.116755_wp
      case ("lc-wpbe")
           rs6 =1.366361_wp
           s18 =1.280619_wp
           rs18=0.003160_wp
      case DEFAULT
           call errore('dft-d3', 'functional name unknown', 1)
    end select
    endif

    ! DFT-D3 with Becke-Johnson finite-damping, variant 2 with their radii
    ! SE: Alp is only used in 3-body calculations
    if (version.eq.4)then
      s6=1.0d0
      alp =14.0d0

      select case (func)
      case ("b-p")
        rs6 =0.3946_wp
        s18 =3.2822_wp
        rs18=4.8516_wp
      case ("b-lyp")
        rs6 =0.4298_wp
        s18 =2.6996_wp
        rs18=4.2359_wp
      case ("revpbe")
        rs6 =0.5238_wp
        s18 =2.3550_wp
        rs18=3.5016_wp
      case ("rpbe")
        rs6 =0.1820_wp
        s18 =0.8318_wp
        rs18=4.0094_wp
      case ("b97-d")
        rs6 =0.5545_wp
        s18 =2.2609_wp
        rs18=3.2297_wp
      case ("pbe")
        rs6 =0.4289_wp
        s18 =0.7875_wp
        rs18=4.4407_wp
      case ("rpw86-pbe")
        rs6 =0.4613_wp
        s18 =1.3845_wp
        rs18=4.5062_wp
      case ("b3-lyp")
        rs6 =0.3981_wp
        s18 =1.9889_wp
        rs18=4.4211_wp
      case ("tpss")
        rs6 =0.4535_wp
        s18 =1.9435_wp
        rs18=4.4752_wp
      case ("hf")
        rs6 =0.3385_wp
        s18 =0.9171_wp
        rs18=2.8830_wp
      case ("tpss0")
        rs6 =0.3768_wp
        s18 =1.2576_wp
        rs18=4.5865_wp
      case ("pbe0")
        rs6 =0.4145_wp
        s18 =1.2177_wp
        rs18=4.8593_wp
      case ("hse06")
        rs6 =0.383_wp
        s18 =2.310_wp
        rs18=5.685_wp
      case ("revpbe38")
        rs6 =0.4309_wp
        s18 =1.4760_wp
        rs18=3.9446_wp
      case ("pw6b95")
        rs6 =0.2076_wp
        s18 =0.7257_wp
        rs18=6.3750
      case ("b2-plyp")
        rs6 =0.3065_wp
        s18 =0.9147_wp
        rs18=5.0570_wp
        s6=0.64d0
      case ("dsd-blyp")
        rs6 =0.0000
        s18 =0.2130_wp
        rs18=6.0519_wp
        s6=0.50d0
      case ("dsd-blyp-fc")
        rs6 =0.0009_wp
        s18 =0.2112_wp
        rs18=5.9807_wp
        s6=0.50d0
      case ("bop")
        rs6 =0.4870_wp
        s18 =3.2950_wp
        rs18=3.5043_wp
      case ("mpwlyp")
        rs6 =0.4831_wp
        s18 =2.0077_wp
        rs18=4.5323_wp
      case ("o-lyp")
        rs6 =0.5299_wp
        s18 =2.6205_wp
        rs18=2.8065_wp
      case ("pbesol")
        rs6 =0.4466_wp
        s18 =2.9491_wp
        rs18=6.1742_wp
      case ("bpbe")
        rs6 =0.4567_wp
        s18 =4.0728_wp
        rs18=4.3908_wp
      case ("opbe")
        rs6 =0.5512_wp
        s18 =3.3816_wp
        rs18=2.9444_wp
      case ("ssb")
        rs6 =-0.0952
        s18 =-0.1744
        rs18=5.2170_wp
      case ("revssb")
        rs6 =0.4720_wp
        s18 =0.4389_wp
        rs18=4.0986_wp
      case ("otpss")
        rs6 =0.4634_wp
        s18 =2.7495_wp
        rs18=4.3153_wp
      case ("b3pw91")
        rs6 =0.4312_wp
        s18 =2.8524_wp
        rs18=4.4693_wp
      case ("bh-lyp")
        rs6 =0.2793_wp
        s18 =1.0354_wp
        rs18=4.9615_wp
      case ("revpbe0")
        rs6 =0.4679_wp
        s18 =1.7588_wp
        rs18=3.7619_wp
      case ("tpssh")
        rs6 =0.4529_wp
        s18 =2.2382_wp
        rs18=4.6550_wp
      case ("mpw1b95")
        rs6 =0.1955_wp
        s18 =1.0508_wp
        rs18=6.4177_wp
      case ("pwb6k")
        rs6 =0.1805_wp
        s18 =0.9383_wp
        rs18=7.7627_wp
      case ("b1b95")
        rs6 =0.2092_wp
        s18 =1.4507_wp
        rs18=5.5545_wp
      case ("bmk")
        rs6 =0.1940_wp
        s18 =2.0860_wp
        rs18=5.9197_wp
      case ("cam-b3lyp")
        rs6 =0.3708_wp
        s18 =2.0674_wp
        rs18=5.4743_wp
      case ("lc-wpbe")
        rs6 =0.3919_wp
        s18 =1.8541_wp
        rs18=5.0897_wp
      case ("b2gp-plyp")
        rs6 =0.0000
        s18 =0.2597_wp
        rs18=6.3332_wp
        s6=0.560_wp
      case ("ptpss")
        rs6 =0.0000
        s18 =0.2804_wp
        rs18=6.5745_wp
        s6=0.750
      case ("pwpb95")
        rs6 =0.0000
        s18 =0.2904_wp
        rs18=7.3141_wp
        s6=0.820_wp
        ! special HF/DFT with eBSSE correction
      case ("hf/mixed")
        rs6 =0.5607_wp
        s18 =3.9027_wp
        rs18=4.5622_wp
      case ("hf/sv")
        rs6 =0.4249_wp
        s18 =2.1849_wp
        rs18=4.2783_wp
      case ("hf/minis")
        rs6 =0.1702_wp
        s18 =0.9841_wp
        rs18=3.8506_wp
      case ("b3-lyp/6-31gd")
        rs6 =0.5014_wp
        s18 =4.0672_wp
        rs18=4.8409_wp
      case ("hcth120")
        rs6=0.3563_wp
        s18=1.0821_wp
        rs18=4.3359_wp
        ! DFTB3 old, deprecated parameters:
        ! case ("dftb3")
        ! rs6=0.7461
        ! s18=3.209
        ! rs18=4.1906
        ! special SCC-DFTB parametrization
        ! full third order DFTB, self consistent charges, hydrogen pair damping
        ! exponent 4.2
      case("dftb3")
        rs6=0.5719d0
        s18=0.5883d0
        rs18=3.6017d0
      case ("pw1pw")
        rs6 =0.3807d0
        s18 =2.3363d0
        rs18=5.8844d0
      case ("pwgga")
        rs6 =0.2211d0
        s18 =2.6910d0
        rs18=6.7278d0
      case ("hsesol")
        rs6 =0.4650d0
        s18 =2.9215d0
        rs18=6.2003d0
        ! special HF-D3-gCP-SRB/MINIX parametrization
      case ("hf3c")
        rs6=0.4171d0
        s18=0.8777d0
        rs18=2.9149d0
        ! special HF-D3-gCP-SRB2/ECP-2G parametrization
      case ("hf3cv")
        rs6=0.3063d0
        s18=0.5022d0
        rs18=3.9856d0
        ! special PBEh-D3-gCP/def2-mSVP parametrization
      case ("pbeh3c", "pbeh-3c")
        rs6=0.4860d0
        s18=0.0000d0
        rs18=4.5000d0
      case ("scan")
        ! Parameters from PRB 94, 115144 (2016); doi: 10.1103/PhysRevB.94.115144
        ! Table 1
        rs6 =0.538_wp
        s18 =0.0
        rs18=5.4200_wp
      case ("r2scan")
        ! Parameters from JCP, 154, 061101 (2021); doi: 10.1063/5.0041008
        ! Table 1
        rs6 =0.4948_wp
        s18 =0.7898_wp
        rs18=5.7308_wp

      case DEFAULT
        call errore('dft-d3', 'functional name unknown', 1)
      end select
    end if

    ! DFT-D3
    if (version.eq.3)then
      s6 =1.0d0
      alp =14.0d0
      rs18=1.0d0
      ! default def2-QZVP (almost basis set limit)
      if (.not.TZ) then
        select case (func)
        case ("slater-dirac-exchange")
          rs6 =0.999_wp
          s18 =-1.957
          rs18=0.697_wp
        case ("b-lyp")
          rs6=1.094_wp
          s18=1.682_wp
        case ("b-p")
          rs6=1.139_wp
          s18=1.683_wp
        case ("b97-d")
          rs6=0.892_wp
          s18=0.909_wp
        case ("revpbe")
          rs6=0.923_wp
          s18=1.010_wp
        case ("pbe")
          rs6=1.217_wp
          s18=0.722_wp
        case ("pbesol")
          rs6=1.345_wp
          s18=0.612_wp
        case ("rpw86-pbe")
          rs6=1.224_wp
          s18=0.901_wp
        case ("rpbe")
          rs6=0.872_wp
          s18=0.514_wp
        case ("tpss")
          rs6=1.166_wp
          s18=1.105_wp
        case ("b3-lyp")
          rs6=1.261_wp
          s18=1.703_wp
        case ("pbe0")
          rs6=1.287_wp
          s18=0.928_wp

        case ("hse06")
          rs6=1.129_wp
          s18=0.109_wp
        case ("revpbe38")
          rs6=1.021_wp
          s18=0.862_wp
        case ("pw6b95")
          rs6=1.532_wp
          s18=0.862_wp
        case ("tpss0")
          rs6=1.252_wp
          s18=1.242_wp
        case ("b2-plyp")
          rs6=1.427_wp
          s18=1.022_wp
          s6=0.64_wp
        case ("pwpb95")
          rs6=1.557_wp
          s18=0.705_wp
          s6=0.82_wp
        case ("b2gp-plyp")
          rs6=1.586_wp
          s18=0.760_wp
          s6=0.56_wp
        case ("ptpss")
          rs6=1.541_wp
          s18=0.879_wp
          s6=0.75
        case ("hf")
          rs6=1.158_wp
          s18=1.746_wp
        case ("mpwlyp")
          rs6=1.239_wp
          s18=1.098_wp
        case ("bpbe")
          rs6=1.087_wp
          s18=2.033_wp
        case ("bh-lyp")
          rs6=1.370_wp
          s18=1.442_wp
        case ("tpssh")
          rs6=1.223_wp
          s18=1.219_wp
        case ("pwb6k")
          rs6=1.660_wp
          s18=0.550_wp
        case ("b1b95")
          rs6=1.613_wp
          s18=1.868_wp
        case ("bop")
          rs6=0.929_wp
          s18=1.975_wp
        case ("o-lyp")
          rs6=0.806_wp
          s18=1.764_wp
        case ("o-pbe")
          rs6=0.837_wp
          s18=2.055_wp
        case ("ssb")
          rs6=1.215_wp
          s18=0.663_wp
        case ("revssb")
          rs6=1.221_wp
          s18=0.560_wp
        case ("otpss")
          rs6=1.128_wp
          s18=1.494_wp
        case ("b3pw91")
          rs6=1.176_wp
          s18=1.775_wp
        case ("revpbe0")
          rs6=0.949_wp
          s18=0.792_wp
        case ("pbe38")
          rs6=1.333_wp
          s18=0.998_wp
        case ("mpw1b95")
          rs6=1.605_wp
          s18=1.118_wp
        case ("mpwb1k")
          rs6=1.671_wp
          s18=1.061_wp
        case ("bmk")
          rs6=1.931_wp
          s18=2.168_wp
        case ("cam-b3lyp")
          rs6=1.378_wp
          s18=1.217_wp
        case ("lc-wpbe")
          rs6=1.355_wp
          s18=1.279_wp
        case ("m05")
          rs6=1.373_wp
          s18=0.595_wp
        case ("m052x")
          rs6=1.417_wp
          s18=0.000
        case ("m06l")
          rs6=1.581_wp
          s18=0.000
        case ("m06")
          rs6=1.325_wp
          s18=0.000
        case ("m062x")
          rs6=1.619_wp
          s18=0.000
        case ("m06hf")
          rs6=1.446_wp
          s18=0.000
          ! DFTB3 (zeta=4.0), old deprecated parameters
          ! case ("dftb3")
          ! rs6=1.235
          ! s18=0.673
        case ("hcth120")
          rs6=1.221_wp
          s18=1.206_wp
        case DEFAULT
          call errore('dft-d3', 'functional name unknown', 1)
        end select
      else
        ! special TZVPP parameter
        select case (func)
        case ("b-lyp")
          rs6=1.243_wp
          s18=2.022_wp
        case ("b-p")
          rs6=1.221_wp
          s18=1.838_wp
        case ("b97-d")
          rs6=0.921_wp
          s18=0.894_wp
        case ("revpbe")
          rs6=0.953_wp
          s18=0.989_wp
        case ("pbe")
          rs6=1.277_wp
          s18=0.777_wp
        case ("tpss")
          rs6=1.213_wp
          s18=1.176_wp
        case ("b3-lyp")
          rs6=1.314_wp
          s18=1.706_wp
        case ("pbe0")
          rs6=1.328_wp
          s18=0.926_wp
        case ("pw6b95")
          rs6=1.562_wp
          s18=0.821_wp
        case ("tpss0")
          rs6=1.282_wp
          s18=1.250
        case ("b2-plyp")
          rs6=1.551_wp
          s18=1.109_wp
          s6=0.5
        case DEFAULT
          call errore('dft-d3', 'functional name unknown (TZ case)', 1)
        end select
      end if
    end if
    ! DFT-D2
    if (version.eq.2)then
      rs6=1.1d0
      s18=0.0d0
      alp=20.0d0
      select case (func)
      case ("b-lyp")
        s6=1.2_wp
      case ("b-p")
        s6=1.05_wp
      case ("b97-d")
        s6=1.25
      case ("revpbe")
        s6=1.25
      case ("pbe")
        s6=0.75
      case ("tpss")
        s6=1.0
      case ("b3-lyp")
        s6=1.05_wp
      case ("pbe0")
        s6=0.6_wp
      case ("pw6b95")
        s6=0.5
      case ("tpss0")
        s6=0.85_wp
      case ("b2-plyp")
        s6=0.55_wp
      case ("b2gp-plyp")
        s6=0.4_wp
      case ("dsd-blyp")
        s6=0.41_wp
        alp=60.0d0
      case DEFAULT
        call errore('dft-d3', 'functional name unknown', 1)
      end select

    end if

  end subroutine setfuncpar


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! The N E W gradC6 routine C
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !
  subroutine get_dC6_dCNij(maxc,max_elem,c6ab,mxci,mxcj,cni,cnj, &
      & izi,izj,iat,jat,c6check,dc6i,dc6j)

    integer maxc,max_elem
    real(wp) c6ab(max_elem,max_elem,maxc,maxc,3)
    !mxc(iz(iat))
    integer mxci,mxcj
    real(wp) cni,cnj,term
    integer iat,jat,izi,izj
    real(wp) dc6i,dc6j,c6check


    integer i,j,a,b
    real(wp) zaehler,nenner,dzaehler_i,dnenner_i,dzaehler_j,dnenner_j
    real(wp) expterm,cn_refi,cn_refj,c6ref,r
    real(wp) c6mem,r_save



    c6mem=-1.d99
    r_save=9999.0
    zaehler=0.0d0
    nenner=0.0d0

    dzaehler_i=0.d0
    dnenner_i=0.d0
    dzaehler_j=0.d0
    dnenner_j=0.d0


    do a=1,mxci
      do b=1,mxcj
        c6ref=c6ab(izi,izj,a,b,1)
        if (c6ref.gt.0) then
          ! c6mem=c6ref
          cn_refi=c6ab(izi,izj,a,b,2)
          cn_refj=c6ab(izi,izj,a,b,3)
          r=(cn_refi-cni)*(cn_refi-cni)+(cn_refj-cnj)*(cn_refj-cnj)
          if (r.lt.r_save) then
            r_save=r
            c6mem=c6ref
          end if
          expterm=exp(k3*r)
          zaehler=zaehler+c6ref*expterm
          nenner=nenner+expterm
          expterm=expterm*2.d0*k3
          term=expterm*(cni-cn_refi)
          dzaehler_i=dzaehler_i+c6ref*term
          dnenner_i =dnenner_i + term

          term=expterm*(cnj-cn_refj)
          dzaehler_j=dzaehler_j+c6ref*term
          dnenner_j =dnenner_j + term
        end if
        !b
      end do
      !a
    end do

    if (nenner.gt.1.0d-99) then
      c6check=zaehler/nenner
      dc6i=((dzaehler_i*nenner)-(dnenner_i*zaehler)) &
          & /(nenner*nenner)
      dc6j=((dzaehler_j*nenner)-(dnenner_j*zaehler)) &
          & /(nenner*nenner)
    else
      c6check=c6mem
      dc6i=0.0d0
      dc6j=0.0d0
    end if
  end subroutine get_dC6_dCNij



  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! interpolate c6
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine getc6(maxc,max_elem,c6ab,mxc,iat,jat,nci,ncj,c6)
    integer maxc,max_elem
    integer iat,jat,i,j,mxc(max_elem)
    real(wp) nci,ncj,c6,c6mem
    real(wp) c6ab(max_elem,max_elem,maxc,maxc,3)
    ! the exponential is sensitive to numerics
    ! when nci or ncj is much larger than cn1/cn2
    real(wp) cn1,cn2,r,rsum,csum,tmp,tmp1
    real(wp) r_save

    c6mem=-1.d+99
    rsum=0.0d0
    csum=0.0d0
    c6 =0.0d0
    r_save=1.0d99
    do i=1,mxc(iat)
      do j=1,mxc(jat)
        c6=c6ab(iat,jat,i,j,1)
        if (c6.gt.0)then
          ! c6mem=c6
          cn1=c6ab(iat,jat,i,j,2)
          cn2=c6ab(iat,jat,i,j,3)
          ! distance
          r=(cn1-nci)**2+(cn2-ncj)**2
          if (r.lt.r_save) then
            r_save=r
            c6mem=c6
          end if
          tmp1=exp(k3*r)
          rsum=rsum+tmp1
          csum=csum+tmp1*c6
        end if
      end do
    end do

    if (rsum.gt.1.0d-99)then
      c6=csum/rsum
    else
      c6=c6mem
    end if

  end subroutine getc6

  integer function lin(i1,i2)
!$acc routine  seq
    integer i1,i2,idum1,idum2
    idum1=max(i1,i2)
    idum2=min(i1,i2)
    lin=idum2+idum1*(idum1-1)/2
    return
  end function lin


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! set cut-off radii
  ! in parts due to INTEL compiler bug
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine setr0ab(max_elem,autoang,r)
    integer max_elem,i,j,k
    real(wp) r(max_elem,max_elem),autoang
    real(wp) r0ab(4465)
    r0ab( 1: 70)=(/ &
        & 2.1823_wp, 1.8547_wp, 1.7347_wp, 2.9086_wp, 2.5732_wp, 3.4956_wp, 2.3550_wp &
        &, 2.5095_wp, 2.9802_wp, 3.0982_wp, 2.5141_wp, 2.3917_wp, 2.9977_wp, 2.9484_wp &
        &, 3.2160_wp, 2.4492_wp, 2.2527_wp, 3.1933_wp, 3.0214_wp, 2.9531_wp, 2.9103_wp &
        &, 2.3667_wp, 2.1328_wp, 2.8784_wp, 2.7660_wp, 2.7776_wp, 2.7063_wp, 2.6225_wp &
        &, 2.1768_wp, 2.0625_wp, 2.6395_wp, 2.6648_wp, 2.6482_wp, 2.5697_wp, 2.4846_wp &
        &, 2.4817_wp, 2.0646_wp, 1.9891_wp, 2.5086_wp, 2.6908_wp, 2.6233_wp, 2.4770_wp &
        &, 2.3885_wp, 2.3511_wp, 2.2996_wp, 1.9892_wp, 1.9251_wp, 2.4190_wp, 2.5473_wp &
        &, 2.4994_wp, 2.4091_wp, 2.3176_wp, 2.2571_wp, 2.1946_wp, 2.1374_wp, 2.9898_wp &
        &, 2.6397_wp, 3.6031_wp, 3.1219_wp, 3.7620_wp, 3.2485_wp, 2.9357_wp, 2.7093_wp &
        &, 2.5781_wp, 2.4839_wp, 3.7082_wp, 2.5129_wp, 2.7321_wp, 3.1052_wp, 3.2962_wp &
        &/)
    r0ab( 71: 140)=(/ &
        & 3.1331_wp, 3.2000_wp, 2.9586_wp, 3.0822_wp, 2.8582_wp, 2.7120_wp, 3.2570_wp &
        &, 3.4839_wp, 2.8766_wp, 2.7427_wp, 3.2776_wp, 3.2363_wp, 3.5929_wp, 3.2826_wp &
        &, 3.0911_wp, 2.9369_wp, 2.9030_wp, 2.7789_wp, 3.3921_wp, 3.3970_wp, 4.0106_wp &
        &, 2.8884_wp, 2.6605_wp, 3.7513_wp, 3.1613_wp, 3.3605_wp, 3.3325_wp, 3.0991_wp &
        &, 2.9297_wp, 2.8674_wp, 2.7571_wp, 3.8129_wp, 3.3266_wp, 3.7105_wp, 3.7917_wp &
        &, 2.8304_wp, 2.5538_wp, 3.3932_wp, 3.1193_wp, 3.1866_wp, 3.1245_wp, 3.0465_wp &
        &, 2.8727_wp, 2.7664_wp, 2.6926_wp, 3.4608_wp, 3.2984_wp, 3.5142_wp, 3.5418_wp &
        &, 3.5017_wp, 2.6190_wp, 2.4797_wp, 3.1331_wp, 3.0540_wp, 3.0651_wp, 2.9879_wp &
        &, 2.9054_wp, 2.8805_wp, 2.7330_wp, 2.6331_wp, 3.2096_wp, 3.5668_wp, 3.3684_wp &
        &, 3.3686_wp, 3.3180_wp, 3.3107_wp, 2.4757_wp, 2.4019_wp, 2.9789_wp, 3.1468_wp &
        &/)
    r0ab( 141: 210)=(/ &
        & 2.9768_wp, 2.8848_wp, 2.7952_wp, 2.7457_wp, 2.6881_wp, 2.5728_wp, 3.0574_wp &
        &, 3.3264_wp, 3.3562_wp, 3.2529_wp, 3.1916_wp, 3.1523_wp, 3.1046_wp, 2.3725_wp &
        &, 2.3289_wp, 2.8760_wp, 2.9804_wp, 2.9093_wp, 2.8040_wp, 2.7071_wp, 2.6386_wp &
        &, 2.5720_wp, 2.5139_wp, 2.9517_wp, 3.1606_wp, 3.2085_wp, 3.1692_wp, 3.0982_wp &
        &, 3.0352_wp, 2.9730_wp, 2.9148_wp, 3.2147_wp, 2.8315_wp, 3.8724_wp, 3.4621_wp &
        &, 3.8823_wp, 3.3760_wp, 3.0746_wp, 2.8817_wp, 2.7552_wp, 2.6605_wp, 3.9740_wp &
        &, 3.6192_wp, 3.6569_wp, 3.9586_wp, 3.6188_wp, 3.3917_wp, 3.2479_wp, 3.1434_wp &
        &, 4.2411_wp, 2.7597_wp, 3.0588_wp, 3.3474_wp, 3.6214_wp, 3.4353_wp, 3.4729_wp &
        &, 3.2487_wp, 3.3200_wp, 3.0914_wp, 2.9403_wp, 3.4972_wp, 3.7993_wp, 3.6773_wp &
        &, 3.8678_wp, 3.5808_wp, 3.8243_wp, 3.5826_wp, 3.4156_wp, 3.8765_wp, 4.1035_wp &
        &/)
    r0ab( 211: 280)=(/ &
        & 2.7361_wp, 2.9765_wp, 3.2475_wp, 3.5004_wp, 3.4185_wp, 3.4378_wp, 3.2084_wp &
        &, 3.2787_wp, 3.0604_wp, 2.9187_wp, 3.4037_wp, 3.6759_wp, 3.6586_wp, 3.8327_wp &
        &, 3.5372_wp, 3.7665_wp, 3.5310_wp, 3.3700_wp, 3.7788_wp, 3.9804_wp, 3.8903_wp &
        &, 2.6832_wp, 2.9060_wp, 3.2613_wp, 3.4359_wp, 3.3538_wp, 3.3860_wp, 3.1550_wp &
        &, 3.2300_wp, 3.0133_wp, 2.8736_wp, 3.4024_wp, 3.6142_wp, 3.5979_wp, 3.5295_wp &
        &, 3.4834_wp, 3.7140_wp, 3.4782_wp, 3.3170_wp, 3.7434_wp, 3.9623_wp, 3.8181_wp &
        &, 3.7642_wp, 2.6379_wp, 2.8494_wp, 3.1840_wp, 3.4225_wp, 3.2771_wp, 3.3401_wp &
        &, 3.1072_wp, 3.1885_wp, 2.9714_wp, 2.8319_wp, 3.3315_wp, 3.5979_wp, 3.5256_wp &
        &, 3.4980_wp, 3.4376_wp, 3.6714_wp, 3.4346_wp, 3.2723_wp, 3.6859_wp, 3.8985_wp &
        &, 3.7918_wp, 3.7372_wp, 3.7211_wp, 2.9230_wp, 2.6223_wp, 3.4161_wp, 2.8999_wp &
        &/)
    r0ab( 281: 350)=(/ &
        & 3.0557_wp, 3.3308_wp, 3.0555_wp, 2.8508_wp, 2.7385_wp, 2.6640_wp, 3.5263_wp &
        &, 3.0277_wp, 3.2990_wp, 3.7721_wp, 3.5017_wp, 3.2751_wp, 3.1368_wp, 3.0435_wp &
        &, 3.7873_wp, 3.2858_wp, 3.2140_wp, 3.1727_wp, 3.2178_wp, 3.4414_wp, 2.5490_wp &
        &, 2.7623_wp, 3.0991_wp, 3.3252_wp, 3.1836_wp, 3.2428_wp, 3.0259_wp, 3.1225_wp &
        &, 2.9032_wp, 2.7621_wp, 3.2490_wp, 3.5110_wp, 3.4429_wp, 3.3845_wp, 3.3574_wp &
        &, 3.6045_wp, 3.3658_wp, 3.2013_wp, 3.6110_wp, 3.8241_wp, 3.7090_wp, 3.6496_wp &
        &, 3.6333_wp, 3.0896_wp, 3.5462_wp, 2.4926_wp, 2.7136_wp, 3.0693_wp, 3.2699_wp &
        &, 3.1272_wp, 3.1893_wp, 2.9658_wp, 3.0972_wp, 2.8778_wp, 2.7358_wp, 3.2206_wp &
        &, 3.4566_wp, 3.3896_wp, 3.3257_wp, 3.2946_wp, 3.5693_wp, 3.3312_wp, 3.1670_wp &
        &, 3.5805_wp, 3.7711_wp, 3.6536_wp, 3.5927_wp, 3.5775_wp, 3.0411_wp, 3.4885_wp &
        &/)
    r0ab( 351: 420)=(/ &
        & 3.4421_wp, 2.4667_wp, 2.6709_wp, 3.0575_wp, 3.2357_wp, 3.0908_wp, 3.1537_wp &
        &, 2.9235_wp, 3.0669_wp, 2.8476_wp, 2.7054_wp, 3.2064_wp, 3.4519_wp, 3.3593_wp &
        &, 3.2921_wp, 3.2577_wp, 3.2161_wp, 3.2982_wp, 3.1339_wp, 3.5606_wp, 3.7582_wp &
        &, 3.6432_wp, 3.5833_wp, 3.5691_wp, 3.0161_wp, 3.4812_wp, 3.4339_wp, 3.4327_wp &
        &, 2.4515_wp, 2.6338_wp, 3.0511_wp, 3.2229_wp, 3.0630_wp, 3.1265_wp, 2.8909_wp &
        &, 3.0253_wp, 2.8184_wp, 2.6764_wp, 3.1968_wp, 3.4114_wp, 3.3492_wp, 3.2691_wp &
        &, 3.2320_wp, 3.1786_wp, 3.2680_wp, 3.1036_wp, 3.5453_wp, 3.7259_wp, 3.6090_wp &
        &, 3.5473_wp, 3.5327_wp, 3.0018_wp, 3.4413_wp, 3.3907_wp, 3.3593_wp, 3.3462_wp &
        &, 2.4413_wp, 2.6006_wp, 3.0540_wp, 3.1987_wp, 3.0490_wp, 3.1058_wp, 2.8643_wp &
        &, 2.9948_wp, 2.7908_wp, 2.6491_wp, 3.1950_wp, 3.3922_wp, 3.3316_wp, 3.2585_wp &
        &/)
    r0ab( 421: 490)=(/ &
        & 3.2136_wp, 3.1516_wp, 3.2364_wp, 3.0752_wp, 3.5368_wp, 3.7117_wp, 3.5941_wp &
        &, 3.5313_wp, 3.5164_wp, 2.9962_wp, 3.4225_wp, 3.3699_wp, 3.3370_wp, 3.3234_wp &
        &, 3.3008_wp, 2.4318_wp, 2.5729_wp, 3.0416_wp, 3.1639_wp, 3.0196_wp, 3.0843_wp &
        &, 2.8413_wp, 2.7436_wp, 2.7608_wp, 2.6271_wp, 3.1811_wp, 3.3591_wp, 3.3045_wp &
        &, 3.2349_wp, 3.1942_wp, 3.1291_wp, 3.2111_wp, 3.0534_wp, 3.5189_wp, 3.6809_wp &
        &, 3.5635_wp, 3.5001_wp, 3.4854_wp, 2.9857_wp, 3.3897_wp, 3.3363_wp, 3.3027_wp &
        &, 3.2890_wp, 3.2655_wp, 3.2309_wp, 2.8502_wp, 2.6934_wp, 3.2467_wp, 3.1921_wp &
        &, 3.5663_wp, 3.2541_wp, 3.0571_wp, 2.9048_wp, 2.8657_wp, 2.7438_wp, 3.3547_wp &
        &, 3.3510_wp, 3.9837_wp, 3.6871_wp, 3.4862_wp, 3.3389_wp, 3.2413_wp, 3.1708_wp &
        &, 3.6096_wp, 3.6280_wp, 3.6860_wp, 3.5568_wp, 3.4836_wp, 3.2868_wp, 3.3994_wp &
        &/)
    r0ab( 491: 560)=(/ &
        & 3.3476_wp, 3.3170_wp, 3.2950_wp, 3.2874_wp, 3.2606_wp, 3.9579_wp, 2.9226_wp &
        &, 2.6838_wp, 3.7867_wp, 3.1732_wp, 3.3872_wp, 3.3643_wp, 3.1267_wp, 2.9541_wp &
        &, 2.8505_wp, 2.7781_wp, 3.8475_wp, 3.3336_wp, 3.7359_wp, 3.8266_wp, 3.5733_wp &
        &, 3.3959_wp, 3.2775_wp, 3.1915_wp, 3.9878_wp, 3.8816_wp, 3.5810_wp, 3.5364_wp &
        &, 3.5060_wp, 3.8097_wp, 3.3925_wp, 3.3348_wp, 3.3019_wp, 3.2796_wp, 3.2662_wp &
        &, 3.2464_wp, 3.7136_wp, 3.8619_wp, 2.9140_wp, 2.6271_wp, 3.4771_wp, 3.1774_wp &
        &, 3.2560_wp, 3.1970_wp, 3.1207_wp, 2.9406_wp, 2.8322_wp, 2.7571_wp, 3.5455_wp &
        &, 3.3514_wp, 3.5837_wp, 3.6177_wp, 3.5816_wp, 3.3902_wp, 3.2604_wp, 3.1652_wp &
        &, 3.7037_wp, 3.6283_wp, 3.5858_wp, 3.5330_wp, 3.4884_wp, 3.5789_wp, 3.4094_wp &
        &, 3.3473_wp, 3.3118_wp, 3.2876_wp, 3.2707_wp, 3.2521_wp, 3.5570_wp, 3.6496_wp &
        &/)
    r0ab( 561: 630)=(/ &
        & 3.6625_wp, 2.7300_wp, 2.5870_wp, 3.2471_wp, 3.1487_wp, 3.1667_wp, 3.0914_wp &
        &, 3.0107_wp, 2.9812_wp, 2.8300_wp, 2.7284_wp, 3.3259_wp, 3.3182_wp, 3.4707_wp &
        &, 3.4748_wp, 3.4279_wp, 3.4182_wp, 3.2547_wp, 3.1353_wp, 3.5116_wp, 3.9432_wp &
        &, 3.8828_wp, 3.8303_wp, 3.7880_wp, 3.3760_wp, 3.7218_wp, 3.3408_wp, 3.3059_wp &
        &, 3.2698_wp, 3.2446_wp, 3.2229_wp, 3.4422_wp, 3.5023_wp, 3.5009_wp, 3.5268_wp &
        &, 2.6026_wp, 2.5355_wp, 3.1129_wp, 3.2863_wp, 3.1029_wp, 3.0108_wp, 2.9227_wp &
        &, 2.8694_wp, 2.8109_wp, 2.6929_wp, 3.1958_wp, 3.4670_wp, 3.4018_wp, 3.3805_wp &
        &, 3.3218_wp, 3.2815_wp, 3.2346_wp, 3.0994_wp, 3.3937_wp, 3.7266_wp, 3.6697_wp &
        &, 3.6164_wp, 3.5730_wp, 3.2522_wp, 3.5051_wp, 3.4686_wp, 3.4355_wp, 3.4084_wp &
        &, 3.3748_wp, 3.3496_wp, 3.3692_wp, 3.4052_wp, 3.3910_wp, 3.3849_wp, 3.3662_wp &
        &/)
    r0ab( 631: 700)=(/ &
        & 2.5087_wp, 2.4814_wp, 3.0239_wp, 3.1312_wp, 3.0535_wp, 2.9457_wp, 2.8496_wp &
        &, 2.7780_wp, 2.7828_wp, 2.6532_wp, 3.1063_wp, 3.3143_wp, 3.3549_wp, 3.3120_wp &
        &, 3.2421_wp, 3.1787_wp, 3.1176_wp, 3.0613_wp, 3.3082_wp, 3.5755_wp, 3.5222_wp &
        &, 3.4678_wp, 3.4231_wp, 3.1684_wp, 3.3528_wp, 3.3162_wp, 3.2827_wp, 3.2527_wp &
        &, 3.2308_wp, 3.2029_wp, 3.3173_wp, 3.3343_wp, 3.3092_wp, 3.2795_wp, 3.2452_wp &
        &, 3.2096_wp, 3.2893_wp, 2.8991_wp, 4.0388_wp, 3.6100_wp, 3.9388_wp, 3.4475_wp &
        &, 3.1590_wp, 2.9812_wp, 2.8586_wp, 2.7683_wp, 4.1428_wp, 3.7911_wp, 3.8225_wp &
        &, 4.0372_wp, 3.7059_wp, 3.4935_wp, 3.3529_wp, 3.2492_wp, 4.4352_wp, 4.0826_wp &
        &, 3.9733_wp, 3.9254_wp, 3.8646_wp, 3.9315_wp, 3.7837_wp, 3.7465_wp, 3.7211_wp &
        &, 3.7012_wp, 3.6893_wp, 3.6676_wp, 3.7736_wp, 4.0660_wp, 3.7926_wp, 3.6158_wp &
        &/)
    r0ab( 701: 770)=(/ &
        & 3.5017_wp, 3.4166_wp, 4.6176_wp, 2.8786_wp, 3.1658_wp, 3.5823_wp, 3.7689_wp &
        &, 3.5762_wp, 3.5789_wp, 3.3552_wp, 3.4004_wp, 3.1722_wp, 3.0212_wp, 3.7241_wp &
        &, 3.9604_wp, 3.8500_wp, 3.9844_wp, 3.7035_wp, 3.9161_wp, 3.6751_wp, 3.5075_wp &
        &, 4.1151_wp, 4.2877_wp, 4.1579_wp, 4.1247_wp, 4.0617_wp, 3.4874_wp, 3.9848_wp &
        &, 3.9280_wp, 3.9079_wp, 3.8751_wp, 3.8604_wp, 3.8277_wp, 3.8002_wp, 3.9981_wp &
        &, 3.7544_wp, 4.0371_wp, 3.8225_wp, 3.6718_wp, 4.3092_wp, 4.4764_wp, 2.8997_wp &
        &, 3.0953_wp, 3.4524_wp, 3.6107_wp, 3.6062_wp, 3.5783_wp, 3.3463_wp, 3.3855_wp &
        &, 3.1746_wp, 3.0381_wp, 3.6019_wp, 3.7938_wp, 3.8697_wp, 3.9781_wp, 3.6877_wp &
        &, 3.8736_wp, 3.6451_wp, 3.4890_wp, 3.9858_wp, 4.1179_wp, 4.0430_wp, 3.9563_wp &
        &, 3.9182_wp, 3.4002_wp, 3.8310_wp, 3.7716_wp, 3.7543_wp, 3.7203_wp, 3.7053_wp &
        &/)
    r0ab( 771: 840)=(/ &
        & 3.6742_wp, 3.8318_wp, 3.7631_wp, 3.7392_wp, 3.9892_wp, 3.7832_wp, 3.6406_wp &
        &, 4.1701_wp, 4.3016_wp, 4.2196_wp, 2.8535_wp, 3.0167_wp, 3.3978_wp, 3.5363_wp &
        &, 3.5393_wp, 3.5301_wp, 3.2960_wp, 3.3352_wp, 3.1287_wp, 2.9967_wp, 3.6659_wp &
        &, 3.7239_wp, 3.8070_wp, 3.7165_wp, 3.6368_wp, 3.8162_wp, 3.5885_wp, 3.4336_wp &
        &, 3.9829_wp, 4.0529_wp, 3.9584_wp, 3.9025_wp, 3.8607_wp, 3.3673_wp, 3.7658_wp &
        &, 3.7035_wp, 3.6866_wp, 3.6504_wp, 3.6339_wp, 3.6024_wp, 3.7708_wp, 3.7283_wp &
        &, 3.6896_wp, 3.9315_wp, 3.7250_wp, 3.5819_wp, 4.1457_wp, 4.2280_wp, 4.1130_wp &
        &, 4.0597_wp, 3.0905_wp, 2.7998_wp, 3.6448_wp, 3.0739_wp, 3.2996_wp, 3.5262_wp &
        &, 3.2559_wp, 3.0518_wp, 2.9394_wp, 2.8658_wp, 3.7514_wp, 3.2295_wp, 3.5643_wp &
        &, 3.7808_wp, 3.6931_wp, 3.4723_wp, 3.3357_wp, 3.2429_wp, 4.0280_wp, 3.5589_wp &
        &/)
    r0ab( 841: 910)=(/ &
        & 3.4636_wp, 3.4994_wp, 3.4309_wp, 3.6177_wp, 3.2946_wp, 3.2376_wp, 3.2050_wp &
        &, 3.1847_wp, 3.1715_wp, 3.1599_wp, 3.5555_wp, 3.8111_wp, 3.7693_wp, 3.5718_wp &
        &, 3.4498_wp, 3.3662_wp, 4.1608_wp, 3.7417_wp, 3.6536_wp, 3.6154_wp, 3.8596_wp &
        &, 3.0301_wp, 2.7312_wp, 3.5821_wp, 3.0473_wp, 3.2137_wp, 3.4679_wp, 3.1975_wp &
        &, 2.9969_wp, 2.8847_wp, 2.8110_wp, 3.6931_wp, 3.2076_wp, 3.4943_wp, 3.5956_wp &
        &, 3.6379_wp, 3.4190_wp, 3.2808_wp, 3.1860_wp, 3.9850_wp, 3.5105_wp, 3.4330_wp &
        &, 3.3797_wp, 3.4155_wp, 3.6033_wp, 3.2737_wp, 3.2145_wp, 3.1807_wp, 3.1596_wp &
        &, 3.1461_wp, 3.1337_wp, 3.4812_wp, 3.6251_wp, 3.7152_wp, 3.5201_wp, 3.3966_wp &
        &, 3.3107_wp, 4.1128_wp, 3.6899_wp, 3.6082_wp, 3.5604_wp, 3.7834_wp, 3.7543_wp &
        &, 2.9189_wp, 2.6777_wp, 3.4925_wp, 2.9648_wp, 3.1216_wp, 3.2940_wp, 3.0975_wp &
        &/)
    r0ab( 911: 980)=(/ &
        & 2.9757_wp, 2.8493_wp, 2.7638_wp, 3.6085_wp, 3.1214_wp, 3.4006_wp, 3.4793_wp &
        &, 3.5147_wp, 3.3806_wp, 3.2356_wp, 3.1335_wp, 3.9144_wp, 3.4183_wp, 3.3369_wp &
        &, 3.2803_wp, 3.2679_wp, 3.4871_wp, 3.1714_wp, 3.1521_wp, 3.1101_wp, 3.0843_wp &
        &, 3.0670_wp, 3.0539_wp, 3.3890_wp, 3.5086_wp, 3.5895_wp, 3.4783_wp, 3.3484_wp &
        &, 3.2559_wp, 4.0422_wp, 3.5967_wp, 3.5113_wp, 3.4576_wp, 3.6594_wp, 3.6313_wp &
        &, 3.5690_wp, 2.8578_wp, 2.6334_wp, 3.4673_wp, 2.9245_wp, 3.0732_wp, 3.2435_wp &
        &, 3.0338_wp, 2.9462_wp, 2.8143_wp, 2.7240_wp, 3.5832_wp, 3.0789_wp, 3.3617_wp &
        &, 3.4246_wp, 3.4505_wp, 3.3443_wp, 3.1964_wp, 3.0913_wp, 3.8921_wp, 3.3713_wp &
        &, 3.2873_wp, 3.2281_wp, 3.2165_wp, 3.4386_wp, 3.1164_wp, 3.1220_wp, 3.0761_wp &
        &, 3.0480_wp, 3.0295_wp, 3.0155_wp, 3.3495_wp, 3.4543_wp, 3.5260_wp, 3.4413_wp &
        &/)
    r0ab( 981:1050)=(/ &
        & 3.3085_wp, 3.2134_wp, 4.0170_wp, 3.5464_wp, 3.4587_wp, 3.4006_wp, 3.6027_wp &
        &, 3.5730_wp, 3.4945_wp, 3.4623_wp, 2.8240_wp, 2.5960_wp, 3.4635_wp, 2.9032_wp &
        &, 3.0431_wp, 3.2115_wp, 2.9892_wp, 2.9148_wp, 2.7801_wp, 2.6873_wp, 3.5776_wp &
        &, 3.0568_wp, 3.3433_wp, 3.3949_wp, 3.4132_wp, 3.3116_wp, 3.1616_wp, 3.0548_wp &
        &, 3.8859_wp, 3.3719_wp, 3.2917_wp, 3.2345_wp, 3.2274_wp, 3.4171_wp, 3.1293_wp &
        &, 3.0567_wp, 3.0565_wp, 3.0274_wp, 3.0087_wp, 2.9939_wp, 3.3293_wp, 3.4249_wp &
        &, 3.4902_wp, 3.4091_wp, 3.2744_wp, 3.1776_wp, 4.0078_wp, 3.5374_wp, 3.4537_wp &
        &, 3.3956_wp, 3.5747_wp, 3.5430_wp, 3.4522_wp, 3.4160_wp, 3.3975_wp, 2.8004_wp &
        &, 2.5621_wp, 3.4617_wp, 2.9154_wp, 3.0203_wp, 3.1875_wp, 2.9548_wp, 2.8038_wp &
        &, 2.7472_wp, 2.6530_wp, 3.5736_wp, 3.0584_wp, 3.3304_wp, 3.3748_wp, 3.3871_wp &
        &/)
    r0ab(1051:1120)=(/ &
        & 3.2028_wp, 3.1296_wp, 3.0214_wp, 3.8796_wp, 3.3337_wp, 3.2492_wp, 3.1883_wp &
        &, 3.1802_wp, 3.4050_wp, 3.0756_wp, 3.0478_wp, 3.0322_wp, 3.0323_wp, 3.0163_wp &
        &, 3.0019_wp, 3.3145_wp, 3.4050_wp, 3.4656_wp, 3.3021_wp, 3.2433_wp, 3.1453_wp &
        &, 3.9991_wp, 3.5017_wp, 3.4141_wp, 3.3520_wp, 3.5583_wp, 3.5251_wp, 3.4243_wp &
        &, 3.3851_wp, 3.3662_wp, 3.3525_wp, 2.7846_wp, 2.5324_wp, 3.4652_wp, 2.8759_wp &
        &, 3.0051_wp, 3.1692_wp, 2.9273_wp, 2.7615_wp, 2.7164_wp, 2.6212_wp, 3.5744_wp &
        &, 3.0275_wp, 3.3249_wp, 3.3627_wp, 3.3686_wp, 3.1669_wp, 3.0584_wp, 2.9915_wp &
        &, 3.8773_wp, 3.3099_wp, 3.2231_wp, 3.1600_wp, 3.1520_wp, 3.4023_wp, 3.0426_wp &
        &, 3.0099_wp, 2.9920_wp, 2.9809_wp, 2.9800_wp, 2.9646_wp, 3.3068_wp, 3.3930_wp &
        &, 3.4486_wp, 3.2682_wp, 3.1729_wp, 3.1168_wp, 3.9952_wp, 3.4796_wp, 3.3901_wp &
        &/)
    r0ab(1121:1190)=(/ &
        & 3.3255_wp, 3.5530_wp, 3.5183_wp, 3.4097_wp, 3.3683_wp, 3.3492_wp, 3.3360_wp &
        &, 3.3308_wp, 2.5424_wp, 2.6601_wp, 3.2555_wp, 3.2807_wp, 3.1384_wp, 3.1737_wp &
        &, 2.9397_wp, 2.8429_wp, 2.8492_wp, 2.7225_wp, 3.3875_wp, 3.4910_wp, 3.4520_wp &
        &, 3.3608_wp, 3.3036_wp, 3.2345_wp, 3.2999_wp, 3.1487_wp, 3.7409_wp, 3.8392_wp &
        &, 3.7148_wp, 3.6439_wp, 3.6182_wp, 3.1753_wp, 3.5210_wp, 3.4639_wp, 3.4265_wp &
        &, 3.4075_wp, 3.3828_wp, 3.3474_wp, 3.4071_wp, 3.3754_wp, 3.3646_wp, 3.3308_wp &
        &, 3.4393_wp, 3.2993_wp, 3.8768_wp, 3.9891_wp, 3.8310_wp, 3.7483_wp, 3.3417_wp &
        &, 3.3019_wp, 3.2250_wp, 3.1832_wp, 3.1578_wp, 3.1564_wp, 3.1224_wp, 3.4620_wp &
        &, 2.9743_wp, 2.8058_wp, 3.4830_wp, 3.3474_wp, 3.6863_wp, 3.3617_wp, 3.1608_wp &
        &, 3.0069_wp, 2.9640_wp, 2.8427_wp, 3.5885_wp, 3.5219_wp, 4.1314_wp, 3.8120_wp &
        &/)
    r0ab(1191:1260)=(/ &
        & 3.6015_wp, 3.4502_wp, 3.3498_wp, 3.2777_wp, 3.8635_wp, 3.8232_wp, 3.8486_wp &
        &, 3.7215_wp, 3.6487_wp, 3.4724_wp, 3.5627_wp, 3.5087_wp, 3.4757_wp, 3.4517_wp &
        &, 3.4423_wp, 3.4139_wp, 4.1028_wp, 3.8388_wp, 3.6745_wp, 3.5562_wp, 3.4806_wp &
        &, 3.4272_wp, 4.0182_wp, 3.9991_wp, 4.0007_wp, 3.9282_wp, 3.7238_wp, 3.6498_wp &
        &, 3.5605_wp, 3.5211_wp, 3.5009_wp, 3.4859_wp, 3.4785_wp, 3.5621_wp, 4.2623_wp &
        &, 3.0775_wp, 2.8275_wp, 4.0181_wp, 3.3385_wp, 3.5379_wp, 3.5036_wp, 3.2589_wp &
        &, 3.0804_wp, 3.0094_wp, 2.9003_wp, 4.0869_wp, 3.5088_wp, 3.9105_wp, 3.9833_wp &
        &, 3.7176_wp, 3.5323_wp, 3.4102_wp, 3.3227_wp, 4.2702_wp, 4.0888_wp, 3.7560_wp &
        &, 3.7687_wp, 3.6681_wp, 3.6405_wp, 3.5569_wp, 3.4990_wp, 3.4659_wp, 3.4433_wp &
        &, 3.4330_wp, 3.4092_wp, 3.8867_wp, 4.0190_wp, 3.7961_wp, 3.6412_wp, 3.5405_wp &
        &/)
    r0ab(1261:1330)=(/ &
        & 3.4681_wp, 4.3538_wp, 4.2136_wp, 3.9381_wp, 3.8912_wp, 3.9681_wp, 3.7909_wp &
        &, 3.6774_wp, 3.6262_wp, 3.5999_wp, 3.5823_wp, 3.5727_wp, 3.5419_wp, 4.0245_wp &
        &, 4.1874_wp, 3.0893_wp, 2.7917_wp, 3.7262_wp, 3.3518_wp, 3.4241_wp, 3.5433_wp &
        &, 3.2773_wp, 3.0890_wp, 2.9775_wp, 2.9010_wp, 3.8048_wp, 3.5362_wp, 3.7746_wp &
        &, 3.7911_wp, 3.7511_wp, 3.5495_wp, 3.4149_wp, 3.3177_wp, 4.0129_wp, 3.8370_wp &
        &, 3.7739_wp, 3.7125_wp, 3.7152_wp, 3.7701_wp, 3.5813_wp, 3.5187_wp, 3.4835_wp &
        &, 3.4595_wp, 3.4439_wp, 3.4242_wp, 3.7476_wp, 3.8239_wp, 3.8346_wp, 3.6627_wp &
        &, 3.5479_wp, 3.4639_wp, 4.1026_wp, 3.9733_wp, 3.9292_wp, 3.8667_wp, 3.9513_wp &
        &, 3.8959_wp, 3.7698_wp, 3.7089_wp, 3.6765_wp, 3.6548_wp, 3.6409_wp, 3.5398_wp &
        &, 3.8759_wp, 3.9804_wp, 4.0150_wp, 2.9091_wp, 2.7638_wp, 3.5066_wp, 3.3377_wp &
        &/)
    r0ab(1331:1400)=(/ &
        & 3.3481_wp, 3.2633_wp, 3.1810_wp, 3.1428_wp, 2.9872_wp, 2.8837_wp, 3.5929_wp &
        &, 3.5183_wp, 3.6729_wp, 3.6596_wp, 3.6082_wp, 3.5927_wp, 3.4224_wp, 3.2997_wp &
        &, 3.8190_wp, 4.1865_wp, 4.1114_wp, 4.0540_wp, 3.6325_wp, 3.5697_wp, 3.5561_wp &
        &, 3.5259_wp, 3.4901_wp, 3.4552_wp, 3.4315_wp, 3.4091_wp, 3.6438_wp, 3.6879_wp &
        &, 3.6832_wp, 3.7043_wp, 3.5557_wp, 3.4466_wp, 3.9203_wp, 4.2919_wp, 4.2196_wp &
        &, 4.1542_wp, 3.7573_wp, 3.7039_wp, 3.6546_wp, 3.6151_wp, 3.5293_wp, 3.4849_wp &
        &, 3.4552_wp, 3.5192_wp, 3.7673_wp, 3.8359_wp, 3.8525_wp, 3.8901_wp, 2.7806_wp &
        &, 2.7209_wp, 3.3812_wp, 3.4958_wp, 3.2913_wp, 3.1888_wp, 3.0990_wp, 3.0394_wp &
        &, 2.9789_wp, 2.8582_wp, 3.4716_wp, 3.6883_wp, 3.6105_wp, 3.5704_wp, 3.5059_wp &
        &, 3.4619_wp, 3.4138_wp, 3.2742_wp, 3.7080_wp, 3.9773_wp, 3.9010_wp, 3.8409_wp &
        &/)
    r0ab(1401:1470)=(/ &
        & 3.7944_wp, 3.4465_wp, 3.7235_wp, 3.6808_wp, 3.6453_wp, 3.6168_wp, 3.5844_wp &
        &, 3.5576_wp, 3.5772_wp, 3.5959_wp, 3.5768_wp, 3.5678_wp, 3.5486_wp, 3.4228_wp &
        &, 3.8107_wp, 4.0866_wp, 4.0169_wp, 3.9476_wp, 3.6358_wp, 3.5800_wp, 3.5260_wp &
        &, 3.4838_wp, 3.4501_wp, 3.4204_wp, 3.3553_wp, 3.6487_wp, 3.6973_wp, 3.7398_wp &
        &, 3.7405_wp, 3.7459_wp, 3.7380_wp, 2.6848_wp, 2.6740_wp, 3.2925_wp, 3.3386_wp &
        &, 3.2473_wp, 3.1284_wp, 3.0301_wp, 2.9531_wp, 2.9602_wp, 2.8272_wp, 3.3830_wp &
        &, 3.5358_wp, 3.5672_wp, 3.5049_wp, 3.4284_wp, 3.3621_wp, 3.3001_wp, 3.2451_wp &
        &, 3.6209_wp, 3.8299_wp, 3.7543_wp, 3.6920_wp, 3.6436_wp, 3.3598_wp, 3.5701_wp &
        &, 3.5266_wp, 3.4904_wp, 3.4590_wp, 3.4364_wp, 3.4077_wp, 3.5287_wp, 3.5280_wp &
        &, 3.4969_wp, 3.4650_wp, 3.4304_wp, 3.3963_wp, 3.7229_wp, 3.9402_wp, 3.8753_wp &
        &/)
    r0ab(1471:1540)=(/ &
        & 3.8035_wp, 3.5499_wp, 3.4913_wp, 3.4319_wp, 3.3873_wp, 3.3520_wp, 3.3209_wp &
        &, 3.2948_wp, 3.5052_wp, 3.6465_wp, 3.6696_wp, 3.6577_wp, 3.6388_wp, 3.6142_wp &
        &, 3.5889_wp, 3.3968_wp, 3.0122_wp, 4.2241_wp, 3.7887_wp, 4.0049_wp, 3.5384_wp &
        &, 3.2698_wp, 3.1083_wp, 2.9917_wp, 2.9057_wp, 4.3340_wp, 3.9900_wp, 4.6588_wp &
        &, 4.1278_wp, 3.8125_wp, 3.6189_wp, 3.4851_wp, 3.3859_wp, 4.6531_wp, 4.3134_wp &
        &, 4.2258_wp, 4.1309_wp, 4.0692_wp, 4.0944_wp, 3.9850_wp, 3.9416_wp, 3.9112_wp &
        &, 3.8873_wp, 3.8736_wp, 3.8473_wp, 4.6027_wp, 4.1538_wp, 3.8994_wp, 3.7419_wp &
        &, 3.6356_wp, 3.5548_wp, 4.8353_wp, 4.5413_wp, 4.3891_wp, 4.3416_wp, 4.3243_wp &
        &, 4.2753_wp, 4.2053_wp, 4.1790_wp, 4.1685_wp, 4.1585_wp, 4.1536_wp, 4.0579_wp &
        &, 4.1980_wp, 4.4564_wp, 4.2192_wp, 4.0528_wp, 3.9489_wp, 3.8642_wp, 5.0567_wp &
        &/)
    r0ab(1541:1610)=(/ &
        & 3.0630_wp, 3.3271_wp, 4.0432_wp, 4.0046_wp, 4.1555_wp, 3.7426_wp, 3.5130_wp &
        &, 3.5174_wp, 3.2884_wp, 3.1378_wp, 4.1894_wp, 4.2321_wp, 4.1725_wp, 4.1833_wp &
        &, 3.8929_wp, 4.0544_wp, 3.8118_wp, 3.6414_wp, 4.6373_wp, 4.6268_wp, 4.4750_wp &
        &, 4.4134_wp, 4.3458_wp, 3.8582_wp, 4.2583_wp, 4.1898_wp, 4.1562_wp, 4.1191_wp &
        &, 4.1069_wp, 4.0639_wp, 4.1257_wp, 4.1974_wp, 3.9532_wp, 4.1794_wp, 3.9660_wp &
        &, 3.8130_wp, 4.8160_wp, 4.8272_wp, 4.6294_wp, 4.5840_wp, 4.0770_wp, 4.0088_wp &
        &, 3.9103_wp, 3.8536_wp, 3.8324_wp, 3.7995_wp, 3.7826_wp, 4.2294_wp, 4.3380_wp &
        &, 4.4352_wp, 4.1933_wp, 4.4580_wp, 4.2554_wp, 4.1072_wp, 5.0454_wp, 5.1814_wp &
        &, 3.0632_wp, 3.2662_wp, 3.6432_wp, 3.8088_wp, 3.7910_wp, 3.7381_wp, 3.5093_wp &
        &, 3.5155_wp, 3.3047_wp, 3.1681_wp, 3.7871_wp, 3.9924_wp, 4.0637_wp, 4.1382_wp &
        &/)
    r0ab(1611:1680)=(/ &
        & 3.8591_wp, 4.0164_wp, 3.7878_wp, 3.6316_wp, 4.1741_wp, 4.3166_wp, 4.2395_wp &
        &, 4.1831_wp, 4.1107_wp, 3.5857_wp, 4.0270_wp, 3.9676_wp, 3.9463_wp, 3.9150_wp &
        &, 3.9021_wp, 3.8708_wp, 4.0240_wp, 4.1551_wp, 3.9108_wp, 4.1337_wp, 3.9289_wp &
        &, 3.7873_wp, 4.3666_wp, 4.5080_wp, 4.4232_wp, 4.3155_wp, 3.8461_wp, 3.8007_wp &
        &, 3.6991_wp, 3.6447_wp, 3.6308_wp, 3.5959_wp, 3.5749_wp, 4.0359_wp, 4.3124_wp &
        &, 4.3539_wp, 4.1122_wp, 4.3772_wp, 4.1785_wp, 4.0386_wp, 4.7004_wp, 4.8604_wp &
        &, 4.6261_wp, 2.9455_wp, 3.2470_wp, 3.6108_wp, 3.8522_wp, 3.6625_wp, 3.6598_wp &
        &, 3.4411_wp, 3.4660_wp, 3.2415_wp, 3.0944_wp, 3.7514_wp, 4.0397_wp, 3.9231_wp &
        &, 4.0561_wp, 3.7860_wp, 3.9845_wp, 3.7454_wp, 3.5802_wp, 4.1366_wp, 4.3581_wp &
        &, 4.2351_wp, 4.2011_wp, 4.1402_wp, 3.5381_wp, 4.0653_wp, 4.0093_wp, 3.9883_wp &
        &/)
    r0ab(1681:1750)=(/ &
        & 3.9570_wp, 3.9429_wp, 3.9112_wp, 3.8728_wp, 4.0682_wp, 3.8351_wp, 4.1054_wp &
        &, 3.8928_wp, 3.7445_wp, 4.3415_wp, 4.5497_wp, 4.3833_wp, 4.3122_wp, 3.8051_wp &
        &, 3.7583_wp, 3.6622_wp, 3.6108_wp, 3.5971_wp, 3.5628_wp, 3.5408_wp, 4.0780_wp &
        &, 4.0727_wp, 4.2836_wp, 4.0553_wp, 4.3647_wp, 4.1622_wp, 4.0178_wp, 4.5802_wp &
        &, 4.9125_wp, 4.5861_wp, 4.6201_wp, 2.9244_wp, 3.2241_wp, 3.5848_wp, 3.8293_wp &
        &, 3.6395_wp, 3.6400_wp, 3.4204_wp, 3.4499_wp, 3.2253_wp, 3.0779_wp, 3.7257_wp &
        &, 4.0170_wp, 3.9003_wp, 4.0372_wp, 3.7653_wp, 3.9672_wp, 3.7283_wp, 3.5630_wp &
        &, 4.1092_wp, 4.3347_wp, 4.2117_wp, 4.1793_wp, 4.1179_wp, 3.5139_wp, 4.0426_wp &
        &, 3.9867_wp, 3.9661_wp, 3.9345_wp, 3.9200_wp, 3.8883_wp, 3.8498_wp, 4.0496_wp &
        &, 3.8145_wp, 4.0881_wp, 3.8756_wp, 3.7271_wp, 4.3128_wp, 4.5242_wp, 4.3578_wp &
        &/)
    r0ab(1751:1820)=(/ &
        & 4.2870_wp, 3.7796_wp, 3.7318_wp, 3.6364_wp, 3.5854_wp, 3.5726_wp, 3.5378_wp &
        &, 3.5155_wp, 4.0527_wp, 4.0478_wp, 4.2630_wp, 4.0322_wp, 4.3449_wp, 4.1421_wp &
        &, 3.9975_wp, 4.5499_wp, 4.8825_wp, 4.5601_wp, 4.5950_wp, 4.5702_wp, 2.9046_wp &
        &, 3.2044_wp, 3.5621_wp, 3.8078_wp, 3.6185_wp, 3.6220_wp, 3.4019_wp, 3.4359_wp &
        &, 3.2110_wp, 3.0635_wp, 3.7037_wp, 3.9958_wp, 3.8792_wp, 4.0194_wp, 3.7460_wp &
        &, 3.9517_wp, 3.7128_wp, 3.5474_wp, 4.0872_wp, 4.3138_wp, 4.1906_wp, 4.1593_wp &
        &, 4.0973_wp, 3.4919_wp, 4.0216_wp, 3.9657_wp, 3.9454_wp, 3.9134_wp, 3.8986_wp &
        &, 3.8669_wp, 3.8289_wp, 4.0323_wp, 3.7954_wp, 4.0725_wp, 3.8598_wp, 3.7113_wp &
        &, 4.2896_wp, 4.5021_wp, 4.3325_wp, 4.2645_wp, 3.7571_wp, 3.7083_wp, 3.6136_wp &
        &, 3.5628_wp, 3.5507_wp, 3.5155_wp, 3.4929_wp, 4.0297_wp, 4.0234_wp, 4.2442_wp &
        &/)
    r0ab(1821:1890)=(/ &
        & 4.0112_wp, 4.3274_wp, 4.1240_wp, 3.9793_wp, 4.5257_wp, 4.8568_wp, 4.5353_wp &
        &, 4.5733_wp, 4.5485_wp, 4.5271_wp, 2.8878_wp, 3.1890_wp, 3.5412_wp, 3.7908_wp &
        &, 3.5974_wp, 3.6078_wp, 3.3871_wp, 3.4243_wp, 3.1992_wp, 3.0513_wp, 3.6831_wp &
        &, 3.9784_wp, 3.8579_wp, 4.0049_wp, 3.7304_wp, 3.9392_wp, 3.7002_wp, 3.5347_wp &
        &, 4.0657_wp, 4.2955_wp, 4.1705_wp, 4.1424_wp, 4.0800_wp, 3.4717_wp, 4.0043_wp &
        &, 3.9485_wp, 3.9286_wp, 3.8965_wp, 3.8815_wp, 3.8500_wp, 3.8073_wp, 4.0180_wp &
        &, 3.7796_wp, 4.0598_wp, 3.8470_wp, 3.6983_wp, 4.2678_wp, 4.4830_wp, 4.3132_wp &
        &, 4.2444_wp, 3.7370_wp, 3.6876_wp, 3.5935_wp, 3.5428_wp, 3.5314_wp, 3.4958_wp &
        &, 3.4730_wp, 4.0117_wp, 4.0043_wp, 4.2287_wp, 3.9939_wp, 4.3134_wp, 4.1096_wp &
        &, 3.9646_wp, 4.5032_wp, 4.8356_wp, 4.5156_wp, 4.5544_wp, 4.5297_wp, 4.5083_wp &
        &/)
    r0ab(1891:1960)=(/ &
        & 4.4896_wp, 2.8709_wp, 3.1737_wp, 3.5199_wp, 3.7734_wp, 3.5802_wp, 3.5934_wp &
        &, 3.3724_wp, 3.4128_wp, 3.1877_wp, 3.0396_wp, 3.6624_wp, 3.9608_wp, 3.8397_wp &
        &, 3.9893_wp, 3.7145_wp, 3.9266_wp, 3.6877_wp, 3.5222_wp, 4.0448_wp, 4.2771_wp &
        &, 4.1523_wp, 4.1247_wp, 4.0626_wp, 3.4530_wp, 3.9866_wp, 3.9310_wp, 3.9115_wp &
        &, 3.8792_wp, 3.8641_wp, 3.8326_wp, 3.7892_wp, 4.0025_wp, 3.7636_wp, 4.0471_wp &
        &, 3.8343_wp, 3.6854_wp, 4.2464_wp, 4.4635_wp, 4.2939_wp, 4.2252_wp, 3.7169_wp &
        &, 3.6675_wp, 3.5739_wp, 3.5235_wp, 3.5126_wp, 3.4768_wp, 3.4537_wp, 3.9932_wp &
        &, 3.9854_wp, 4.2123_wp, 3.9765_wp, 4.2992_wp, 4.0951_wp, 3.9500_wp, 4.4811_wp &
        &, 4.8135_wp, 4.4959_wp, 4.5351_wp, 4.5105_wp, 4.4891_wp, 4.4705_wp, 4.4515_wp &
        &, 2.8568_wp, 3.1608_wp, 3.5050_wp, 3.7598_wp, 3.5665_wp, 3.5803_wp, 3.3601_wp &
        &/)
    r0ab(1961:2030)=(/ &
        & 3.4031_wp, 3.1779_wp, 3.0296_wp, 3.6479_wp, 3.9471_wp, 3.8262_wp, 3.9773_wp &
        &, 3.7015_wp, 3.9162_wp, 3.6771_wp, 3.5115_wp, 4.0306_wp, 4.2634_wp, 4.1385_wp &
        &, 4.1116_wp, 4.0489_wp, 3.4366_wp, 3.9732_wp, 3.9176_wp, 3.8983_wp, 3.8659_wp &
        &, 3.8507_wp, 3.8191_wp, 3.7757_wp, 3.9907_wp, 3.7506_wp, 4.0365_wp, 3.8235_wp &
        &, 3.6745_wp, 4.2314_wp, 4.4490_wp, 4.2792_wp, 4.2105_wp, 3.7003_wp, 3.6510_wp &
        &, 3.5578_wp, 3.5075_wp, 3.4971_wp, 3.4609_wp, 3.4377_wp, 3.9788_wp, 3.9712_wp &
        &, 4.1997_wp, 3.9624_wp, 4.2877_wp, 4.0831_wp, 3.9378_wp, 4.4655_wp, 4.7974_wp &
        &, 4.4813_wp, 4.5209_wp, 4.4964_wp, 4.4750_wp, 4.4565_wp, 4.4375_wp, 4.4234_wp &
        &, 2.6798_wp, 3.0151_wp, 3.2586_wp, 3.5292_wp, 3.5391_wp, 3.4902_wp, 3.2887_wp &
        &, 3.3322_wp, 3.1228_wp, 2.9888_wp, 3.4012_wp, 3.7145_wp, 3.7830_wp, 3.6665_wp &
        &/)
    r0ab(2031:2100)=(/ &
        & 3.5898_wp, 3.8077_wp, 3.5810_wp, 3.4265_wp, 3.7726_wp, 4.0307_wp, 3.9763_wp &
        &, 3.8890_wp, 3.8489_wp, 3.2706_wp, 3.7595_wp, 3.6984_wp, 3.6772_wp, 3.6428_wp &
        &, 3.6243_wp, 3.5951_wp, 3.7497_wp, 3.6775_wp, 3.6364_wp, 3.9203_wp, 3.7157_wp &
        &, 3.5746_wp, 3.9494_wp, 4.2076_wp, 4.1563_wp, 4.0508_wp, 3.5329_wp, 3.4780_wp &
        &, 3.3731_wp, 3.3126_wp, 3.2846_wp, 3.2426_wp, 3.2135_wp, 3.7491_wp, 3.9006_wp &
        &, 3.8332_wp, 3.8029_wp, 4.1436_wp, 3.9407_wp, 3.7998_wp, 4.1663_wp, 4.5309_wp &
        &, 4.3481_wp, 4.2911_wp, 4.2671_wp, 4.2415_wp, 4.2230_wp, 4.2047_wp, 4.1908_wp &
        &, 4.1243_wp, 2.5189_wp, 2.9703_wp, 3.3063_wp, 3.6235_wp, 3.4517_wp, 3.3989_wp &
        &, 3.2107_wp, 3.2434_wp, 3.0094_wp, 2.8580_wp, 3.4253_wp, 3.8157_wp, 3.7258_wp &
        &, 3.6132_wp, 3.5297_wp, 3.7566_wp, 3.5095_wp, 3.3368_wp, 3.7890_wp, 4.1298_wp &
        &/)
    r0ab(2101:2170)=(/ &
        & 4.0190_wp, 3.9573_wp, 3.9237_wp, 3.2677_wp, 3.8480_wp, 3.8157_wp, 3.7656_wp &
        &, 3.7317_wp, 3.7126_wp, 3.6814_wp, 3.6793_wp, 3.6218_wp, 3.5788_wp, 3.8763_wp &
        &, 3.6572_wp, 3.5022_wp, 3.9737_wp, 4.3255_wp, 4.1828_wp, 4.1158_wp, 3.5078_wp &
        &, 3.4595_wp, 3.3600_wp, 3.3088_wp, 3.2575_wp, 3.2164_wp, 3.1856_wp, 3.8522_wp &
        &, 3.8665_wp, 3.8075_wp, 3.7772_wp, 4.1391_wp, 3.9296_wp, 3.7772_wp, 4.2134_wp &
        &, 4.7308_wp, 4.3787_wp, 4.3894_wp, 4.3649_wp, 4.3441_wp, 4.3257_wp, 4.3073_wp &
        &, 4.2941_wp, 4.1252_wp, 4.2427_wp, 3.0481_wp, 2.9584_wp, 3.6919_wp, 3.5990_wp &
        &, 3.8881_wp, 3.4209_wp, 3.1606_wp, 3.1938_wp, 2.9975_wp, 2.8646_wp, 3.8138_wp &
        &, 3.7935_wp, 3.7081_wp, 3.9155_wp, 3.5910_wp, 3.4808_wp, 3.4886_wp, 3.3397_wp &
        &, 4.1336_wp, 4.1122_wp, 3.9888_wp, 3.9543_wp, 3.8917_wp, 3.5894_wp, 3.8131_wp &
        &/)
    r0ab(2171:2240)=(/ &
        & 3.7635_wp, 3.7419_wp, 3.7071_wp, 3.6880_wp, 3.6574_wp, 3.6546_wp, 3.9375_wp &
        &, 3.6579_wp, 3.5870_wp, 3.6361_wp, 3.5039_wp, 4.3149_wp, 4.2978_wp, 4.1321_wp &
        &, 4.1298_wp, 3.8164_wp, 3.7680_wp, 3.7154_wp, 3.6858_wp, 3.6709_wp, 3.6666_wp &
        &, 3.6517_wp, 3.8174_wp, 3.8608_wp, 4.1805_wp, 3.9102_wp, 3.8394_wp, 3.8968_wp &
        &, 3.7673_wp, 4.5274_wp, 4.6682_wp, 4.3344_wp, 4.3639_wp, 4.3384_wp, 4.3162_wp &
        &, 4.2972_wp, 4.2779_wp, 4.2636_wp, 4.0253_wp, 4.1168_wp, 4.1541_wp, 2.8136_wp &
        &, 3.0951_wp, 3.4635_wp, 3.6875_wp, 3.4987_wp, 3.5183_wp, 3.2937_wp, 3.3580_wp &
        &, 3.1325_wp, 2.9832_wp, 3.6078_wp, 3.8757_wp, 3.7616_wp, 3.9222_wp, 3.6370_wp &
        &, 3.8647_wp, 3.6256_wp, 3.4595_wp, 3.9874_wp, 4.1938_wp, 4.0679_wp, 4.0430_wp &
        &, 3.9781_wp, 3.3886_wp, 3.9008_wp, 3.8463_wp, 3.8288_wp, 3.7950_wp, 3.7790_wp &
        &/)
    r0ab(2241:2310)=(/ &
        & 3.7472_wp, 3.7117_wp, 3.9371_wp, 3.6873_wp, 3.9846_wp, 3.7709_wp, 3.6210_wp &
        &, 4.1812_wp, 4.3750_wp, 4.2044_wp, 4.1340_wp, 3.6459_wp, 3.5929_wp, 3.5036_wp &
        &, 3.4577_wp, 3.4528_wp, 3.4146_wp, 3.3904_wp, 3.9014_wp, 3.9031_wp, 4.1443_wp &
        &, 3.8961_wp, 4.2295_wp, 4.0227_wp, 3.8763_wp, 4.4086_wp, 4.7097_wp, 4.4064_wp &
        &, 4.4488_wp, 4.4243_wp, 4.4029_wp, 4.3842_wp, 4.3655_wp, 4.3514_wp, 4.1162_wp &
        &, 4.2205_wp, 4.1953_wp, 4.2794_wp, 2.8032_wp, 3.0805_wp, 3.4519_wp, 3.6700_wp &
        &, 3.4827_wp, 3.5050_wp, 3.2799_wp, 3.3482_wp, 3.1233_wp, 2.9747_wp, 3.5971_wp &
        &, 3.8586_wp, 3.7461_wp, 3.9100_wp, 3.6228_wp, 3.8535_wp, 3.6147_wp, 3.4490_wp &
        &, 3.9764_wp, 4.1773_wp, 4.0511_wp, 4.0270_wp, 3.9614_wp, 3.3754_wp, 3.8836_wp &
        &, 3.8291_wp, 3.8121_wp, 3.7780_wp, 3.7619_wp, 3.7300_wp, 3.6965_wp, 3.9253_wp &
        &/)
    r0ab(2311:2380)=(/ &
        & 3.6734_wp, 3.9733_wp, 3.7597_wp, 3.6099_wp, 4.1683_wp, 4.3572_wp, 4.1862_wp &
        &, 4.1153_wp, 3.6312_wp, 3.5772_wp, 3.4881_wp, 3.4429_wp, 3.4395_wp, 3.4009_wp &
        &, 3.3766_wp, 3.8827_wp, 3.8868_wp, 4.1316_wp, 3.8807_wp, 4.2164_wp, 4.0092_wp &
        &, 3.8627_wp, 4.3936_wp, 4.6871_wp, 4.3882_wp, 4.4316_wp, 4.4073_wp, 4.3858_wp &
        &, 4.3672_wp, 4.3485_wp, 4.3344_wp, 4.0984_wp, 4.2036_wp, 4.1791_wp, 4.2622_wp &
        &, 4.2450_wp, 2.7967_wp, 3.0689_wp, 3.4445_wp, 3.6581_wp, 3.4717_wp, 3.4951_wp &
        &, 3.2694_wp, 3.3397_wp, 3.1147_wp, 2.9661_wp, 3.5898_wp, 3.8468_wp, 3.7358_wp &
        &, 3.9014_wp, 3.6129_wp, 3.8443_wp, 3.6054_wp, 3.4396_wp, 3.9683_wp, 4.1656_wp &
        &, 4.0394_wp, 4.0158_wp, 3.9498_wp, 3.3677_wp, 3.8718_wp, 3.8164_wp, 3.8005_wp &
        &, 3.7662_wp, 3.7500_wp, 3.7181_wp, 3.6863_wp, 3.9170_wp, 3.6637_wp, 3.9641_wp &
        &/)
    r0ab(2381:2450)=(/ &
        & 3.7503_wp, 3.6004_wp, 4.1590_wp, 4.3448_wp, 4.1739_wp, 4.1029_wp, 3.6224_wp &
        &, 3.5677_wp, 3.4785_wp, 3.4314_wp, 3.4313_wp, 3.3923_wp, 3.3680_wp, 3.8698_wp &
        &, 3.8758_wp, 4.1229_wp, 3.8704_wp, 4.2063_wp, 3.9987_wp, 3.8519_wp, 4.3832_wp &
        &, 4.6728_wp, 4.3759_wp, 4.4195_wp, 4.3952_wp, 4.3737_wp, 4.3551_wp, 4.3364_wp &
        &, 4.3223_wp, 4.0861_wp, 4.1911_wp, 4.1676_wp, 4.2501_wp, 4.2329_wp, 4.2208_wp &
        &, 2.7897_wp, 3.0636_wp, 3.4344_wp, 3.6480_wp, 3.4626_wp, 3.4892_wp, 3.2626_wp &
        &, 3.3344_wp, 3.1088_wp, 2.9597_wp, 3.5804_wp, 3.8359_wp, 3.7251_wp, 3.8940_wp &
        &, 3.6047_wp, 3.8375_wp, 3.5990_wp, 3.4329_wp, 3.9597_wp, 4.1542_wp, 4.0278_wp &
        &, 4.0048_wp, 3.9390_wp, 3.3571_wp, 3.8608_wp, 3.8056_wp, 3.7899_wp, 3.7560_wp &
        &, 3.7400_wp, 3.7081_wp, 3.6758_wp, 3.9095_wp, 3.6552_wp, 3.9572_wp, 3.7436_wp &
        &/)
    r0ab(2451:2520)=(/ &
        & 3.5933_wp, 4.1508_wp, 4.3337_wp, 4.1624_wp, 4.0916_wp, 3.6126_wp, 3.5582_wp &
        &, 3.4684_wp, 3.4212_wp, 3.4207_wp, 3.3829_wp, 3.3586_wp, 3.8604_wp, 3.8658_wp &
        &, 4.1156_wp, 3.8620_wp, 4.1994_wp, 3.9917_wp, 3.8446_wp, 4.3750_wp, 4.6617_wp &
        &, 4.3644_wp, 4.4083_wp, 4.3840_wp, 4.3625_wp, 4.3439_wp, 4.3253_wp, 4.3112_wp &
        &, 4.0745_wp, 4.1807_wp, 4.1578_wp, 4.2390_wp, 4.2218_wp, 4.2097_wp, 4.1986_wp &
        &, 2.8395_wp, 3.0081_wp, 3.3171_wp, 3.4878_wp, 3.5360_wp, 3.5145_wp, 3.2809_wp &
        &, 3.3307_wp, 3.1260_wp, 2.9940_wp, 3.4741_wp, 3.6675_wp, 3.7832_wp, 3.6787_wp &
        &, 3.6156_wp, 3.8041_wp, 3.5813_wp, 3.4301_wp, 3.8480_wp, 3.9849_wp, 3.9314_wp &
        &, 3.8405_wp, 3.8029_wp, 3.2962_wp, 3.7104_wp, 3.6515_wp, 3.6378_wp, 3.6020_wp &
        &, 3.5849_wp, 3.5550_wp, 3.7494_wp, 3.6893_wp, 3.6666_wp, 3.9170_wp, 3.7150_wp &
        &/)
    r0ab(2521:2590)=(/ &
        & 3.5760_wp, 4.0268_wp, 4.1596_wp, 4.1107_wp, 3.9995_wp, 3.5574_wp, 3.5103_wp &
        &, 3.4163_wp, 3.3655_wp, 3.3677_wp, 3.3243_wp, 3.2975_wp, 3.7071_wp, 3.9047_wp &
        &, 3.8514_wp, 3.8422_wp, 3.8022_wp, 3.9323_wp, 3.7932_wp, 4.2343_wp, 4.4583_wp &
        &, 4.3115_wp, 4.2457_wp, 4.2213_wp, 4.1945_wp, 4.1756_wp, 4.1569_wp, 4.1424_wp &
        &, 4.0620_wp, 4.0494_wp, 3.9953_wp, 4.0694_wp, 4.0516_wp, 4.0396_wp, 4.0280_wp &
        &, 4.0130_wp, 2.9007_wp, 2.9674_wp, 3.8174_wp, 3.5856_wp, 3.6486_wp, 3.5339_wp &
        &, 3.2832_wp, 3.3154_wp, 3.1144_wp, 2.9866_wp, 3.9618_wp, 3.8430_wp, 3.9980_wp &
        &, 3.8134_wp, 3.6652_wp, 3.7985_wp, 3.5756_wp, 3.4207_wp, 4.4061_wp, 4.2817_wp &
        &, 4.1477_wp, 4.0616_wp, 3.9979_wp, 3.6492_wp, 3.8833_wp, 3.8027_wp, 3.7660_wp &
        &, 3.7183_wp, 3.6954_wp, 3.6525_wp, 3.9669_wp, 3.8371_wp, 3.7325_wp, 3.9160_wp &
        &/)
    r0ab(2591:2660)=(/ &
        & 3.7156_wp, 3.5714_wp, 4.6036_wp, 4.4620_wp, 4.3092_wp, 4.2122_wp, 3.8478_wp &
        &, 3.7572_wp, 3.6597_wp, 3.5969_wp, 3.5575_wp, 3.5386_wp, 3.5153_wp, 3.7818_wp &
        &, 4.1335_wp, 4.0153_wp, 3.9177_wp, 3.8603_wp, 3.9365_wp, 3.7906_wp, 4.7936_wp &
        &, 4.7410_wp, 4.5461_wp, 4.5662_wp, 4.5340_wp, 4.5059_wp, 4.4832_wp, 4.4604_wp &
        &, 4.4429_wp, 4.2346_wp, 4.4204_wp, 4.3119_wp, 4.3450_wp, 4.3193_wp, 4.3035_wp &
        &, 4.2933_wp, 4.1582_wp, 4.2450_wp, 2.8559_wp, 2.9050_wp, 3.8325_wp, 3.5442_wp &
        &, 3.5077_wp, 3.4905_wp, 3.2396_wp, 3.2720_wp, 3.0726_wp, 2.9467_wp, 3.9644_wp &
        &, 3.8050_wp, 3.8981_wp, 3.7762_wp, 3.6216_wp, 3.7531_wp, 3.5297_wp, 3.3742_wp &
        &, 4.3814_wp, 4.2818_wp, 4.1026_wp, 4.0294_wp, 3.9640_wp, 3.6208_wp, 3.8464_wp &
        &, 3.7648_wp, 3.7281_wp, 3.6790_wp, 3.6542_wp, 3.6117_wp, 3.8650_wp, 3.8010_wp &
        &/)
    r0ab(2661:2730)=(/ &
        & 3.6894_wp, 3.8713_wp, 3.6699_wp, 3.5244_wp, 4.5151_wp, 4.4517_wp, 4.2538_wp &
        &, 4.1483_wp, 3.8641_wp, 3.7244_wp, 3.6243_wp, 3.5589_wp, 3.5172_wp, 3.4973_wp &
        &, 3.4715_wp, 3.7340_wp, 4.0316_wp, 3.9958_wp, 3.8687_wp, 3.8115_wp, 3.8862_wp &
        &, 3.7379_wp, 4.7091_wp, 4.7156_wp, 4.5199_wp, 4.5542_wp, 4.5230_wp, 4.4959_wp &
        &, 4.4750_wp, 4.4529_wp, 4.4361_wp, 4.1774_wp, 4.3774_wp, 4.2963_wp, 4.3406_wp &
        &, 4.3159_wp, 4.3006_wp, 4.2910_wp, 4.1008_wp, 4.1568_wp, 4.0980_wp, 2.8110_wp &
        &, 2.8520_wp, 3.7480_wp, 3.5105_wp, 3.4346_wp, 3.3461_wp, 3.1971_wp, 3.2326_wp &
        &, 3.0329_wp, 2.9070_wp, 3.8823_wp, 3.7928_wp, 3.8264_wp, 3.7006_wp, 3.5797_wp &
        &, 3.7141_wp, 3.4894_wp, 3.3326_wp, 4.3048_wp, 4.2217_wp, 4.0786_wp, 3.9900_wp &
        &, 3.9357_wp, 3.6331_wp, 3.8333_wp, 3.7317_wp, 3.6957_wp, 3.6460_wp, 3.6197_wp &
        &/)
    r0ab(2731:2800)=(/ &
        & 3.5779_wp, 3.7909_wp, 3.7257_wp, 3.6476_wp, 3.5729_wp, 3.6304_wp, 3.4834_wp &
        &, 4.4368_wp, 4.3921_wp, 4.2207_wp, 4.1133_wp, 3.8067_wp, 3.7421_wp, 3.6140_wp &
        &, 3.5491_wp, 3.5077_wp, 3.4887_wp, 3.4623_wp, 3.6956_wp, 3.9568_wp, 3.8976_wp &
        &, 3.8240_wp, 3.7684_wp, 3.8451_wp, 3.6949_wp, 4.6318_wp, 4.6559_wp, 4.4533_wp &
        &, 4.4956_wp, 4.4641_wp, 4.4366_wp, 4.4155_wp, 4.3936_wp, 4.3764_wp, 4.1302_wp &
        &, 4.3398_wp, 4.2283_wp, 4.2796_wp, 4.2547_wp, 4.2391_wp, 4.2296_wp, 4.0699_wp &
        &, 4.1083_wp, 4.0319_wp, 3.9855_wp, 2.7676_wp, 2.8078_wp, 3.6725_wp, 3.4804_wp &
        &, 3.3775_wp, 3.2411_wp, 3.1581_wp, 3.1983_wp, 2.9973_wp, 2.8705_wp, 3.8070_wp &
        &, 3.7392_wp, 3.7668_wp, 3.6263_wp, 3.5402_wp, 3.6807_wp, 3.4545_wp, 3.2962_wp &
        &, 4.2283_wp, 4.1698_wp, 4.0240_wp, 3.9341_wp, 3.8711_wp, 3.5489_wp, 3.7798_wp &
        &/)
    r0ab(2801:2870)=(/ &
        & 3.7000_wp, 3.6654_wp, 3.6154_wp, 3.5882_wp, 3.5472_wp, 3.7289_wp, 3.6510_wp &
        &, 3.6078_wp, 3.5355_wp, 3.5963_wp, 3.4480_wp, 4.3587_wp, 4.3390_wp, 4.1635_wp &
        &, 4.0536_wp, 3.7193_wp, 3.6529_wp, 3.5512_wp, 3.4837_wp, 3.4400_wp, 3.4191_wp &
        &, 3.3891_wp, 3.6622_wp, 3.8934_wp, 3.8235_wp, 3.7823_wp, 3.7292_wp, 3.8106_wp &
        &, 3.6589_wp, 4.5535_wp, 4.6013_wp, 4.3961_wp, 4.4423_wp, 4.4109_wp, 4.3835_wp &
        &, 4.3625_wp, 4.3407_wp, 4.3237_wp, 4.0863_wp, 4.2835_wp, 4.1675_wp, 4.2272_wp &
        &, 4.2025_wp, 4.1869_wp, 4.1774_wp, 4.0126_wp, 4.0460_wp, 3.9815_wp, 3.9340_wp &
        &, 3.8955_wp, 2.6912_wp, 2.7604_wp, 3.6037_wp, 3.4194_wp, 3.3094_wp, 3.1710_wp &
        &, 3.0862_wp, 3.1789_wp, 2.9738_wp, 2.8427_wp, 3.7378_wp, 3.6742_wp, 3.6928_wp &
        &, 3.5512_wp, 3.4614_wp, 3.4087_wp, 3.4201_wp, 3.2607_wp, 4.1527_wp, 4.0977_wp &
        &/)
    r0ab(2871:2940)=(/ &
        & 3.9523_wp, 3.8628_wp, 3.8002_wp, 3.4759_wp, 3.7102_wp, 3.6466_wp, 3.6106_wp &
        &, 3.5580_wp, 3.5282_wp, 3.4878_wp, 3.6547_wp, 3.5763_wp, 3.5289_wp, 3.5086_wp &
        &, 3.5593_wp, 3.4099_wp, 4.2788_wp, 4.2624_wp, 4.0873_wp, 3.9770_wp, 3.6407_wp &
        &, 3.5743_wp, 3.5178_wp, 3.4753_wp, 3.3931_wp, 3.3694_wp, 3.3339_wp, 3.6002_wp &
        &, 3.8164_wp, 3.7478_wp, 3.7028_wp, 3.6952_wp, 3.7669_wp, 3.6137_wp, 4.4698_wp &
        &, 4.5488_wp, 4.3168_wp, 4.3646_wp, 4.3338_wp, 4.3067_wp, 4.2860_wp, 4.2645_wp &
        &, 4.2478_wp, 4.0067_wp, 4.2349_wp, 4.0958_wp, 4.1543_wp, 4.1302_wp, 4.1141_wp &
        &, 4.1048_wp, 3.9410_wp, 3.9595_wp, 3.8941_wp, 3.8465_wp, 3.8089_wp, 3.7490_wp &
        &, 2.7895_wp, 2.5849_wp, 3.6484_wp, 3.0162_wp, 3.1267_wp, 3.2125_wp, 3.0043_wp &
        &, 2.9572_wp, 2.8197_wp, 2.7261_wp, 3.7701_wp, 3.2446_wp, 3.5239_wp, 3.4696_wp &
        &/)
    r0ab(2941:3010)=(/ &
        & 3.4261_wp, 3.3508_wp, 3.1968_wp, 3.0848_wp, 4.1496_wp, 3.6598_wp, 3.5111_wp &
        &, 3.4199_wp, 3.3809_wp, 3.5382_wp, 3.2572_wp, 3.2100_wp, 3.1917_wp, 3.1519_wp &
        &, 3.1198_wp, 3.1005_wp, 3.5071_wp, 3.5086_wp, 3.5073_wp, 3.4509_wp, 3.3120_wp &
        &, 3.2082_wp, 4.2611_wp, 3.8117_wp, 3.6988_wp, 3.5646_wp, 3.6925_wp, 3.6295_wp &
        &, 3.5383_wp, 3.4910_wp, 3.4625_wp, 3.4233_wp, 3.4007_wp, 3.2329_wp, 3.6723_wp &
        &, 3.6845_wp, 3.6876_wp, 3.6197_wp, 3.4799_wp, 3.3737_wp, 4.4341_wp, 4.0525_wp &
        &, 3.9011_wp, 3.8945_wp, 3.8635_wp, 3.8368_wp, 3.8153_wp, 3.7936_wp, 3.7758_wp &
        &, 3.4944_wp, 3.4873_wp, 3.9040_wp, 3.7110_wp, 3.6922_wp, 3.6799_wp, 3.6724_wp &
        &, 3.5622_wp, 3.6081_wp, 3.5426_wp, 3.4922_wp, 3.4498_wp, 3.3984_wp, 3.4456_wp &
        &, 2.7522_wp, 2.5524_wp, 3.5742_wp, 2.9508_wp, 3.0751_wp, 3.0158_wp, 2.9644_wp &
        &/)
    r0ab(3011:3080)=(/ &
        & 2.8338_wp, 2.7891_wp, 2.6933_wp, 3.6926_wp, 3.1814_wp, 3.4528_wp, 3.4186_wp &
        &, 3.3836_wp, 3.2213_wp, 3.1626_wp, 3.0507_wp, 4.0548_wp, 3.5312_wp, 3.4244_wp &
        &, 3.3409_wp, 3.2810_wp, 3.4782_wp, 3.1905_wp, 3.1494_wp, 3.1221_wp, 3.1128_wp &
        &, 3.0853_wp, 3.0384_wp, 3.4366_wp, 3.4562_wp, 3.4638_wp, 3.3211_wp, 3.2762_wp &
        &, 3.1730_wp, 4.1632_wp, 3.6825_wp, 3.5822_wp, 3.4870_wp, 3.6325_wp, 3.5740_wp &
        &, 3.4733_wp, 3.4247_wp, 3.3969_wp, 3.3764_wp, 3.3525_wp, 3.1984_wp, 3.5989_wp &
        &, 3.6299_wp, 3.6433_wp, 3.4937_wp, 3.4417_wp, 3.3365_wp, 4.3304_wp, 3.9242_wp &
        &, 3.7793_wp, 3.7623_wp, 3.7327_wp, 3.7071_wp, 3.6860_wp, 3.6650_wp, 3.6476_wp &
        &, 3.3849_wp, 3.3534_wp, 3.8216_wp, 3.5870_wp, 3.5695_wp, 3.5584_wp, 3.5508_wp &
        &, 3.4856_wp, 3.5523_wp, 3.4934_wp, 3.4464_wp, 3.4055_wp, 3.3551_wp, 3.3888_wp &
        &/)
    r0ab(3081:3150)=(/ &
        & 3.3525_wp, 2.7202_wp, 2.5183_wp, 3.4947_wp, 2.8731_wp, 3.0198_wp, 3.1457_wp &
        &, 2.9276_wp, 2.7826_wp, 2.7574_wp, 2.6606_wp, 3.6090_wp, 3.0581_wp, 3.3747_wp &
        &, 3.3677_wp, 3.3450_wp, 3.1651_wp, 3.1259_wp, 3.0147_wp, 3.9498_wp, 3.3857_wp &
        &, 3.2917_wp, 3.2154_wp, 3.1604_wp, 3.4174_wp, 3.0735_wp, 3.0342_wp, 3.0096_wp &
        &, 3.0136_wp, 2.9855_wp, 2.9680_wp, 3.3604_wp, 3.4037_wp, 3.4243_wp, 3.2633_wp &
        &, 3.1810_wp, 3.1351_wp, 4.0557_wp, 3.5368_wp, 3.4526_wp, 3.3699_wp, 3.5707_wp &
        &, 3.5184_wp, 3.4085_wp, 3.3595_wp, 3.3333_wp, 3.3143_wp, 3.3041_wp, 3.1094_wp &
        &, 3.5193_wp, 3.5745_wp, 3.6025_wp, 3.4338_wp, 3.3448_wp, 3.2952_wp, 4.2158_wp &
        &, 3.7802_wp, 3.6431_wp, 3.6129_wp, 3.5853_wp, 3.5610_wp, 3.5406_wp, 3.5204_wp &
        &, 3.5036_wp, 3.2679_wp, 3.2162_wp, 3.7068_wp, 3.4483_wp, 3.4323_wp, 3.4221_wp &
        &/)
    r0ab(3151:3220)=(/ &
        & 3.4138_wp, 3.3652_wp, 3.4576_wp, 3.4053_wp, 3.3618_wp, 3.3224_wp, 3.2711_wp &
        &, 3.3326_wp, 3.2950_wp, 3.2564_wp, 2.5315_wp, 2.6104_wp, 3.2734_wp, 3.2299_wp &
        &, 3.1090_wp, 2.9942_wp, 2.9159_wp, 2.8324_wp, 2.8350_wp, 2.7216_wp, 3.3994_wp &
        &, 3.4475_wp, 3.4354_wp, 3.3438_wp, 3.2807_wp, 3.2169_wp, 3.2677_wp, 3.1296_wp &
        &, 3.7493_wp, 3.8075_wp, 3.6846_wp, 3.6104_wp, 3.5577_wp, 3.2052_wp, 3.4803_wp &
        &, 3.4236_wp, 3.3845_wp, 3.3640_wp, 3.3365_wp, 3.3010_wp, 3.3938_wp, 3.3624_wp &
        &, 3.3440_wp, 3.3132_wp, 3.4035_wp, 3.2754_wp, 3.8701_wp, 3.9523_wp, 3.8018_wp &
        &, 3.7149_wp, 3.3673_wp, 3.3199_wp, 3.2483_wp, 3.2069_wp, 3.1793_wp, 3.1558_wp &
        &, 3.1395_wp, 3.4097_wp, 3.5410_wp, 3.5228_wp, 3.5116_wp, 3.4921_wp, 3.4781_wp &
        &, 3.4690_wp, 4.0420_wp, 4.1759_wp, 4.0078_wp, 4.0450_wp, 4.0189_wp, 3.9952_wp &
        &/)
    r0ab(3221:3290)=(/ &
        & 3.9770_wp, 3.9583_wp, 3.9434_wp, 3.7217_wp, 3.8228_wp, 3.7826_wp, 3.8640_wp &
        &, 3.8446_wp, 3.8314_wp, 3.8225_wp, 3.6817_wp, 3.7068_wp, 3.6555_wp, 3.6159_wp &
        &, 3.5831_wp, 3.5257_wp, 3.2133_wp, 3.1689_wp, 3.1196_wp, 3.3599_wp, 2.9852_wp &
        &, 2.7881_wp, 3.5284_wp, 3.3493_wp, 3.6958_wp, 3.3642_wp, 3.1568_wp, 3.0055_wp &
        &, 2.9558_wp, 2.8393_wp, 3.6287_wp, 3.5283_wp, 4.1511_wp, 3.8259_wp, 3.6066_wp &
        &, 3.4527_wp, 3.3480_wp, 3.2713_wp, 3.9037_wp, 3.8361_wp, 3.8579_wp, 3.7311_wp &
        &, 3.6575_wp, 3.5176_wp, 3.5693_wp, 3.5157_wp, 3.4814_wp, 3.4559_wp, 3.4445_wp &
        &, 3.4160_wp, 4.1231_wp, 3.8543_wp, 3.6816_wp, 3.5602_wp, 3.4798_wp, 3.4208_wp &
        &, 4.0542_wp, 4.0139_wp, 4.0165_wp, 3.9412_wp, 3.7698_wp, 3.6915_wp, 3.6043_wp &
        &, 3.5639_wp, 3.5416_wp, 3.5247_wp, 3.5153_wp, 3.5654_wp, 4.2862_wp, 4.0437_wp &
        &/)
    r0ab(3291:3360)=(/ &
        & 3.8871_wp, 3.7741_wp, 3.6985_wp, 3.6413_wp, 4.2345_wp, 4.3663_wp, 4.3257_wp &
        &, 4.0869_wp, 4.0612_wp, 4.0364_wp, 4.0170_wp, 3.9978_wp, 3.9834_wp, 3.9137_wp &
        &, 3.8825_wp, 3.8758_wp, 3.9143_wp, 3.8976_wp, 3.8864_wp, 3.8768_wp, 3.9190_wp &
        &, 4.1613_wp, 4.0566_wp, 3.9784_wp, 3.9116_wp, 3.8326_wp, 3.7122_wp, 3.6378_wp &
        &, 3.5576_wp, 3.5457_wp, 4.3127_wp, 3.1160_wp, 2.8482_wp, 4.0739_wp, 3.3599_wp &
        &, 3.5698_wp, 3.5366_wp, 3.2854_wp, 3.1039_wp, 2.9953_wp, 2.9192_wp, 4.1432_wp &
        &, 3.5320_wp, 3.9478_wp, 4.0231_wp, 3.7509_wp, 3.5604_wp, 3.4340_wp, 3.3426_wp &
        &, 4.3328_wp, 3.8288_wp, 3.7822_wp, 3.7909_wp, 3.6907_wp, 3.6864_wp, 3.5793_wp &
        &, 3.5221_wp, 3.4883_wp, 3.4649_wp, 3.4514_wp, 3.4301_wp, 3.9256_wp, 4.0596_wp &
        &, 3.8307_wp, 3.6702_wp, 3.5651_wp, 3.4884_wp, 4.4182_wp, 4.2516_wp, 3.9687_wp &
        &/)
    r0ab(3361:3430)=(/ &
        & 3.9186_wp, 3.9485_wp, 3.8370_wp, 3.7255_wp, 3.6744_wp, 3.6476_wp, 3.6295_wp &
        &, 3.6193_wp, 3.5659_wp, 4.0663_wp, 4.2309_wp, 4.0183_wp, 3.8680_wp, 3.7672_wp &
        &, 3.6923_wp, 4.5240_wp, 4.4834_wp, 4.1570_wp, 4.3204_wp, 4.2993_wp, 4.2804_wp &
        &, 4.2647_wp, 4.2481_wp, 4.2354_wp, 3.8626_wp, 3.8448_wp, 4.2267_wp, 4.1799_wp &
        &, 4.1670_wp, 3.8738_wp, 3.8643_wp, 3.8796_wp, 4.0575_wp, 4.0354_wp, 3.9365_wp &
        &, 3.8611_wp, 3.7847_wp, 3.7388_wp, 3.6826_wp, 3.6251_wp, 3.5492_wp, 4.0889_wp &
        &, 4.2764_wp, 3.1416_wp, 2.8325_wp, 3.7735_wp, 3.3787_wp, 3.4632_wp, 3.5923_wp &
        &, 3.3214_wp, 3.1285_wp, 3.0147_wp, 2.9366_wp, 3.8527_wp, 3.5602_wp, 3.8131_wp &
        &, 3.8349_wp, 3.7995_wp, 3.5919_wp, 3.4539_wp, 3.3540_wp, 4.0654_wp, 3.8603_wp &
        &, 3.7972_wp, 3.7358_wp, 3.7392_wp, 3.8157_wp, 3.6055_wp, 3.5438_wp, 3.5089_wp &
        &/)
    r0ab(3431:3500)=(/ &
        & 3.4853_wp, 3.4698_wp, 3.4508_wp, 3.7882_wp, 3.8682_wp, 3.8837_wp, 3.7055_wp &
        &, 3.5870_wp, 3.5000_wp, 4.1573_wp, 4.0005_wp, 3.9568_wp, 3.8936_wp, 3.9990_wp &
        &, 3.9433_wp, 3.8172_wp, 3.7566_wp, 3.7246_wp, 3.7033_wp, 3.6900_wp, 3.5697_wp &
        &, 3.9183_wp, 4.0262_wp, 4.0659_wp, 3.8969_wp, 3.7809_wp, 3.6949_wp, 4.2765_wp &
        &, 4.2312_wp, 4.1401_wp, 4.0815_wp, 4.0580_wp, 4.0369_wp, 4.0194_wp, 4.0017_wp &
        &, 3.9874_wp, 3.8312_wp, 3.8120_wp, 3.9454_wp, 3.9210_wp, 3.9055_wp, 3.8951_wp &
        &, 3.8866_wp, 3.8689_wp, 3.9603_wp, 3.9109_wp, 3.9122_wp, 3.8233_wp, 3.7438_wp &
        &, 3.7436_wp, 3.6981_wp, 3.6555_wp, 3.5452_wp, 3.9327_wp, 4.0658_wp, 4.1175_wp &
        &, 2.9664_wp, 2.8209_wp, 3.5547_wp, 3.3796_wp, 3.3985_wp, 3.3164_wp, 3.2364_wp &
        &, 3.1956_wp, 3.0370_wp, 2.9313_wp, 3.6425_wp, 3.5565_wp, 3.7209_wp, 3.7108_wp &
        &/)
    r0ab(3501:3570)=(/ &
        & 3.6639_wp, 3.6484_wp, 3.4745_wp, 3.3492_wp, 3.8755_wp, 4.2457_wp, 3.7758_wp &
        &, 3.7161_wp, 3.6693_wp, 3.6155_wp, 3.5941_wp, 3.5643_wp, 3.5292_wp, 3.4950_wp &
        &, 3.4720_wp, 3.4503_wp, 3.6936_wp, 3.7392_wp, 3.7388_wp, 3.7602_wp, 3.6078_wp &
        &, 3.4960_wp, 3.9800_wp, 4.3518_wp, 4.2802_wp, 3.8580_wp, 3.8056_wp, 3.7527_wp &
        &, 3.7019_wp, 3.6615_wp, 3.5768_wp, 3.5330_wp, 3.5038_wp, 3.5639_wp, 3.8192_wp &
        &, 3.8883_wp, 3.9092_wp, 3.9478_wp, 3.7995_wp, 3.6896_wp, 4.1165_wp, 4.5232_wp &
        &, 4.4357_wp, 4.4226_wp, 4.4031_wp, 4.3860_wp, 4.3721_wp, 4.3580_wp, 4.3466_wp &
        &, 4.2036_wp, 4.2037_wp, 3.8867_wp, 4.2895_wp, 4.2766_wp, 4.2662_wp, 4.2598_wp &
        &, 3.8408_wp, 3.9169_wp, 3.8681_wp, 3.8250_wp, 3.7855_wp, 3.7501_wp, 3.6753_wp &
        &, 3.5499_wp, 3.4872_wp, 3.5401_wp, 3.8288_wp, 3.9217_wp, 3.9538_wp, 4.0054_wp &
        &/)
    r0ab(3571:3640)=(/ &
        & 2.8388_wp, 2.7890_wp, 3.4329_wp, 3.5593_wp, 3.3488_wp, 3.2486_wp, 3.1615_wp &
        &, 3.1000_wp, 3.0394_wp, 2.9165_wp, 3.5267_wp, 3.7479_wp, 3.6650_wp, 3.6263_wp &
        &, 3.5658_wp, 3.5224_wp, 3.4762_wp, 3.3342_wp, 3.7738_wp, 4.0333_wp, 3.9568_wp &
        &, 3.8975_wp, 3.8521_wp, 3.4929_wp, 3.7830_wp, 3.7409_wp, 3.7062_wp, 3.6786_wp &
        &, 3.6471_wp, 3.6208_wp, 3.6337_wp, 3.6519_wp, 3.6363_wp, 3.6278_wp, 3.6110_wp &
        &, 3.4825_wp, 3.8795_wp, 4.1448_wp, 4.0736_wp, 4.0045_wp, 3.6843_wp, 3.6291_wp &
        &, 3.5741_wp, 3.5312_wp, 3.4974_wp, 3.4472_wp, 3.4034_wp, 3.7131_wp, 3.7557_wp &
        &, 3.7966_wp, 3.8005_wp, 3.8068_wp, 3.8015_wp, 3.6747_wp, 4.0222_wp, 4.3207_wp &
        &, 4.2347_wp, 4.2191_wp, 4.1990_wp, 4.1811_wp, 4.1666_wp, 4.1521_wp, 4.1401_wp &
        &, 3.9970_wp, 3.9943_wp, 3.9592_wp, 4.0800_wp, 4.0664_wp, 4.0559_wp, 4.0488_wp &
        &/)
    r0ab(3641:3710)=(/ &
        & 3.9882_wp, 4.0035_wp, 3.9539_wp, 3.9138_wp, 3.8798_wp, 3.8355_wp, 3.5359_wp &
        &, 3.4954_wp, 3.3962_wp, 3.5339_wp, 3.7595_wp, 3.8250_wp, 3.8408_wp, 3.8600_wp &
        &, 3.8644_wp, 2.7412_wp, 2.7489_wp, 3.3374_wp, 3.3950_wp, 3.3076_wp, 3.1910_wp &
        &, 3.0961_wp, 3.0175_wp, 3.0280_wp, 2.8929_wp, 3.4328_wp, 3.5883_wp, 3.6227_wp &
        &, 3.5616_wp, 3.4894_wp, 3.4241_wp, 3.3641_wp, 3.3120_wp, 3.6815_wp, 3.8789_wp &
        &, 3.8031_wp, 3.7413_wp, 3.6939_wp, 3.4010_wp, 3.6225_wp, 3.5797_wp, 3.5443_wp &
        &, 3.5139_wp, 3.4923_wp, 3.4642_wp, 3.5860_wp, 3.5849_wp, 3.5570_wp, 3.5257_wp &
        &, 3.4936_wp, 3.4628_wp, 3.7874_wp, 3.9916_wp, 3.9249_wp, 3.8530_wp, 3.5932_wp &
        &, 3.5355_wp, 3.4757_wp, 3.4306_wp, 3.3953_wp, 3.3646_wp, 3.3390_wp, 3.5637_wp &
        &, 3.7053_wp, 3.7266_wp, 3.7177_wp, 3.6996_wp, 3.6775_wp, 3.6558_wp, 3.9331_wp &
        &/)
    r0ab(3711:3780)=(/ &
        & 4.1655_wp, 4.0879_wp, 4.0681_wp, 4.0479_wp, 4.0299_wp, 4.0152_wp, 4.0006_wp &
        &, 3.9883_wp, 3.8500_wp, 3.8359_wp, 3.8249_wp, 3.9269_wp, 3.9133_wp, 3.9025_wp &
        &, 3.8948_wp, 3.8422_wp, 3.8509_wp, 3.7990_wp, 3.7570_wp, 3.7219_wp, 3.6762_wp &
        &, 3.4260_wp, 3.3866_wp, 3.3425_wp, 3.5294_wp, 3.7022_wp, 3.7497_wp, 3.7542_wp &
        &, 3.7494_wp, 3.7370_wp, 3.7216_wp, 3.4155_wp, 3.0522_wp, 4.2541_wp, 3.8218_wp &
        &, 4.0438_wp, 3.5875_wp, 3.3286_wp, 3.1682_wp, 3.0566_wp, 2.9746_wp, 4.3627_wp &
        &, 4.0249_wp, 4.6947_wp, 4.1718_wp, 3.8639_wp, 3.6735_wp, 3.5435_wp, 3.4479_wp &
        &, 4.6806_wp, 4.3485_wp, 4.2668_wp, 4.1690_wp, 4.1061_wp, 4.1245_wp, 4.0206_wp &
        &, 3.9765_wp, 3.9458_wp, 3.9217_wp, 3.9075_wp, 3.8813_wp, 3.9947_wp, 4.1989_wp &
        &, 3.9507_wp, 3.7960_wp, 3.6925_wp, 3.6150_wp, 4.8535_wp, 4.5642_wp, 4.4134_wp &
        &/)
    r0ab(3781:3850)=(/ &
        & 4.3688_wp, 4.3396_wp, 4.2879_wp, 4.2166_wp, 4.1888_wp, 4.1768_wp, 4.1660_wp &
        &, 4.1608_wp, 4.0745_wp, 4.2289_wp, 4.4863_wp, 4.2513_wp, 4.0897_wp, 3.9876_wp &
        &, 3.9061_wp, 5.0690_wp, 5.0446_wp, 4.6186_wp, 4.6078_wp, 4.5780_wp, 4.5538_wp &
        &, 4.5319_wp, 4.5101_wp, 4.4945_wp, 4.1912_wp, 4.2315_wp, 4.5534_wp, 4.4373_wp &
        &, 4.4224_wp, 4.4120_wp, 4.4040_wp, 4.2634_wp, 4.7770_wp, 4.6890_wp, 4.6107_wp &
        &, 4.5331_wp, 4.4496_wp, 4.4082_wp, 4.3095_wp, 4.2023_wp, 4.0501_wp, 4.2595_wp &
        &, 4.5497_wp, 4.3056_wp, 4.1506_wp, 4.0574_wp, 3.9725_wp, 5.0796_wp, 3.0548_wp &
        &, 3.3206_wp, 3.8132_wp, 3.9720_wp, 3.7675_wp, 3.7351_wp, 3.5167_wp, 3.5274_wp &
        &, 3.3085_wp, 3.1653_wp, 3.9500_wp, 4.1730_wp, 4.0613_wp, 4.1493_wp, 3.8823_wp &
        &, 4.0537_wp, 3.8200_wp, 3.6582_wp, 4.3422_wp, 4.5111_wp, 4.3795_wp, 4.3362_wp &
        &/)
    r0ab(3851:3920)=(/ &
        & 4.2751_wp, 3.7103_wp, 4.1973_wp, 4.1385_wp, 4.1129_wp, 4.0800_wp, 4.0647_wp &
        &, 4.0308_wp, 4.0096_wp, 4.1619_wp, 3.9360_wp, 4.1766_wp, 3.9705_wp, 3.8262_wp &
        &, 4.5348_wp, 4.7025_wp, 4.5268_wp, 4.5076_wp, 3.9562_wp, 3.9065_wp, 3.8119_wp &
        &, 3.7605_wp, 3.7447_wp, 3.7119_wp, 3.6916_wp, 4.1950_wp, 4.2110_wp, 4.3843_wp &
        &, 4.1631_wp, 4.4427_wp, 4.2463_wp, 4.1054_wp, 4.7693_wp, 5.0649_wp, 4.7365_wp &
        &, 4.7761_wp, 4.7498_wp, 4.7272_wp, 4.7076_wp, 4.6877_wp, 4.6730_wp, 4.4274_wp &
        &, 4.5473_wp, 4.5169_wp, 4.5975_wp, 4.5793_wp, 4.5667_wp, 4.5559_wp, 4.3804_wp &
        &, 4.6920_wp, 4.6731_wp, 4.6142_wp, 4.5600_wp, 4.4801_wp, 4.0149_wp, 3.8856_wp &
        &, 3.7407_wp, 4.1545_wp, 4.2253_wp, 4.4229_wp, 4.1923_wp, 4.5022_wp, 4.3059_wp &
        &, 4.1591_wp, 4.7883_wp, 4.9294_wp, 3.3850_wp, 3.4208_wp, 3.7004_wp, 3.8800_wp &
        &/)
    r0ab(3921:3990)=(/ &
        & 3.9886_wp, 3.9040_wp, 3.6719_wp, 3.6547_wp, 3.4625_wp, 3.3370_wp, 3.8394_wp &
        &, 4.0335_wp, 4.2373_wp, 4.3023_wp, 4.0306_wp, 4.1408_wp, 3.9297_wp, 3.7857_wp &
        &, 4.1907_wp, 4.3230_wp, 4.2664_wp, 4.2173_wp, 4.1482_wp, 3.6823_wp, 4.0711_wp &
        &, 4.0180_wp, 4.0017_wp, 3.9747_wp, 3.9634_wp, 3.9383_wp, 4.1993_wp, 4.3205_wp &
        &, 4.0821_wp, 4.2547_wp, 4.0659_wp, 3.9359_wp, 4.3952_wp, 4.5176_wp, 4.3888_wp &
        &, 4.3607_wp, 3.9583_wp, 3.9280_wp, 3.8390_wp, 3.7971_wp, 3.7955_wp, 3.7674_wp &
        &, 3.7521_wp, 4.1062_wp, 4.3633_wp, 4.2991_wp, 4.2767_wp, 4.4857_wp, 4.3039_wp &
        &, 4.1762_wp, 4.6197_wp, 4.8654_wp, 4.6633_wp, 4.5878_wp, 4.5640_wp, 4.5422_wp &
        &, 4.5231_wp, 4.5042_wp, 4.4901_wp, 4.3282_wp, 4.3978_wp, 4.3483_wp, 4.4202_wp &
        &, 4.4039_wp, 4.3926_wp, 4.3807_wp, 4.2649_wp, 4.6135_wp, 4.5605_wp, 4.5232_wp &
        &/)
    r0ab(3991:4060)=(/ &
        & 4.4676_wp, 4.3948_wp, 4.0989_wp, 3.9864_wp, 3.8596_wp, 4.0942_wp, 4.2720_wp &
        &, 4.3270_wp, 4.3022_wp, 4.5410_wp, 4.3576_wp, 4.2235_wp, 4.6545_wp, 4.7447_wp &
        &, 4.7043_wp, 3.0942_wp, 3.2075_wp, 3.5152_wp, 3.6659_wp, 3.8289_wp, 3.7459_wp &
        &, 3.5156_wp, 3.5197_wp, 3.3290_wp, 3.2069_wp, 3.6702_wp, 3.8448_wp, 4.0340_wp &
        &, 3.9509_wp, 3.8585_wp, 3.9894_wp, 3.7787_wp, 3.6365_wp, 4.1425_wp, 4.1618_wp &
        &, 4.0940_wp, 4.0466_wp, 3.9941_wp, 3.5426_wp, 3.8952_wp, 3.8327_wp, 3.8126_wp &
        &, 3.7796_wp, 3.7635_wp, 3.7356_wp, 4.0047_wp, 3.9655_wp, 3.9116_wp, 4.1010_wp &
        &, 3.9102_wp, 3.7800_wp, 4.2964_wp, 4.3330_wp, 4.2622_wp, 4.2254_wp, 3.8195_wp &
        &, 3.7560_wp, 3.6513_wp, 3.5941_wp, 3.5810_wp, 3.5420_wp, 3.5178_wp, 3.8861_wp &
        &, 4.1459_wp, 4.1147_wp, 4.0772_wp, 4.3120_wp, 4.1207_wp, 3.9900_wp, 4.4733_wp &
        &/)
    r0ab(4061:4130)=(/ &
        & 4.6157_wp, 4.4580_wp, 4.4194_wp, 4.3954_wp, 4.3739_wp, 4.3531_wp, 4.3343_wp &
        &, 4.3196_wp, 4.2140_wp, 4.2339_wp, 4.1738_wp, 4.2458_wp, 4.2278_wp, 4.2158_wp &
        &, 4.2039_wp, 4.1658_wp, 4.3595_wp, 4.2857_wp, 4.2444_wp, 4.1855_wp, 4.1122_wp &
        &, 3.7839_wp, 3.6879_wp, 3.5816_wp, 3.8633_wp, 4.1585_wp, 4.1402_wp, 4.1036_wp &
        &, 4.3694_wp, 4.1735_wp, 4.0368_wp, 4.5095_wp, 4.5538_wp, 4.5240_wp, 4.4252_wp &
        &, 3.0187_wp, 3.1918_wp, 3.5127_wp, 3.6875_wp, 3.7404_wp, 3.6943_wp, 3.4702_wp &
        &, 3.4888_wp, 3.2914_wp, 3.1643_wp, 3.6669_wp, 3.8724_wp, 3.9940_wp, 4.0816_wp &
        &, 3.8054_wp, 3.9661_wp, 3.7492_wp, 3.6024_wp, 4.0428_wp, 4.1951_wp, 4.1466_wp &
        &, 4.0515_wp, 4.0075_wp, 3.5020_wp, 3.9158_wp, 3.8546_wp, 3.8342_wp, 3.8008_wp &
        &, 3.7845_wp, 3.7549_wp, 3.9602_wp, 3.8872_wp, 3.8564_wp, 4.0793_wp, 3.8835_wp &
        &/)
    r0ab(4131:4200)=(/ &
        & 3.7495_wp, 4.2213_wp, 4.3704_wp, 4.3300_wp, 4.2121_wp, 3.7643_wp, 3.7130_wp &
        &, 3.6144_wp, 3.5599_wp, 3.5474_wp, 3.5093_wp, 3.4853_wp, 3.9075_wp, 4.1115_wp &
        &, 4.0473_wp, 4.0318_wp, 4.2999_wp, 4.1050_wp, 3.9710_wp, 4.4320_wp, 4.6706_wp &
        &, 4.5273_wp, 4.4581_wp, 4.4332_wp, 4.4064_wp, 4.3873_wp, 4.3684_wp, 4.3537_wp &
        &, 4.2728_wp, 4.2549_wp, 4.2032_wp, 4.2794_wp, 4.2613_wp, 4.2491_wp, 4.2375_wp &
        &, 4.2322_wp, 4.3665_wp, 4.3061_wp, 4.2714_wp, 4.2155_wp, 4.1416_wp, 3.7660_wp &
        &, 3.6628_wp, 3.5476_wp, 3.8790_wp, 4.1233_wp, 4.0738_wp, 4.0575_wp, 4.3575_wp &
        &, 4.1586_wp, 4.0183_wp, 4.4593_wp, 4.5927_wp, 4.4865_wp, 4.3813_wp, 4.4594_wp &
        &, 2.9875_wp, 3.1674_wp, 3.4971_wp, 3.6715_wp, 3.7114_wp, 3.6692_wp, 3.4446_wp &
        &, 3.4676_wp, 3.2685_wp, 3.1405_wp, 3.6546_wp, 3.8579_wp, 3.9637_wp, 4.0581_wp &
        &/)
    r0ab(4201:4270)=(/ &
        & 3.7796_wp, 3.9463_wp, 3.7275_wp, 3.5792_wp, 4.0295_wp, 4.1824_wp, 4.1247_wp &
        &, 4.0357_wp, 3.9926_wp, 3.4827_wp, 3.9007_wp, 3.8392_wp, 3.8191_wp, 3.7851_wp &
        &, 3.7687_wp, 3.7387_wp, 3.9290_wp, 3.8606_wp, 3.8306_wp, 4.0601_wp, 3.8625_wp &
        &, 3.7269_wp, 4.2062_wp, 4.3566_wp, 4.3022_wp, 4.1929_wp, 3.7401_wp, 3.6888_wp &
        &, 3.5900_wp, 3.5350_wp, 3.5226_wp, 3.4838_wp, 3.4594_wp, 3.8888_wp, 4.0813_wp &
        &, 4.0209_wp, 4.0059_wp, 4.2810_wp, 4.0843_wp, 3.9486_wp, 4.4162_wp, 4.6542_wp &
        &, 4.5005_wp, 4.4444_wp, 4.4196_wp, 4.3933_wp, 4.3741_wp, 4.3552_wp, 4.3406_wp &
        &, 4.2484_wp, 4.2413_wp, 4.1907_wp, 4.2656_wp, 4.2474_wp, 4.2352_wp, 4.2236_wp &
        &, 4.2068_wp, 4.3410_wp, 4.2817_wp, 4.2479_wp, 4.1921_wp, 4.1182_wp, 3.7346_wp &
        &, 3.6314_wp, 3.5168_wp, 3.8582_wp, 4.0927_wp, 4.0469_wp, 4.0313_wp, 4.3391_wp &
        &/)
    r0ab(4271:4340)=(/ &
        & 4.1381_wp, 3.9962_wp, 4.4429_wp, 4.5787_wp, 4.4731_wp, 4.3588_wp, 4.4270_wp &
        &, 4.3957_wp, 2.9659_wp, 3.1442_wp, 3.4795_wp, 3.6503_wp, 3.6814_wp, 3.6476_wp &
        &, 3.4222_wp, 3.4491_wp, 3.2494_wp, 3.1209_wp, 3.6324_wp, 3.8375_wp, 3.9397_wp &
        &, 3.8311_wp, 3.7581_wp, 3.9274_wp, 3.7085_wp, 3.5598_wp, 4.0080_wp, 4.1641_wp &
        &, 4.1057_wp, 4.0158_wp, 3.9726_wp, 3.4667_wp, 3.8802_wp, 3.8188_wp, 3.7989_wp &
        &, 3.7644_wp, 3.7474_wp, 3.7173_wp, 3.9049_wp, 3.8424_wp, 3.8095_wp, 4.0412_wp &
        &, 3.8436_wp, 3.7077_wp, 4.1837_wp, 4.3366_wp, 4.2816_wp, 4.1686_wp, 3.7293_wp &
        &, 3.6709_wp, 3.5700_wp, 3.5153_wp, 3.5039_wp, 3.4684_wp, 3.4437_wp, 3.8663_wp &
        &, 4.0575_wp, 4.0020_wp, 3.9842_wp, 4.2612_wp, 4.0643_wp, 3.9285_wp, 4.3928_wp &
        &, 4.6308_wp, 4.4799_wp, 4.4244_wp, 4.3996_wp, 4.3737_wp, 4.3547_wp, 4.3358_wp &
        &/)
    r0ab(4341:4410)=(/ &
        & 4.3212_wp, 4.2275_wp, 4.2216_wp, 4.1676_wp, 4.2465_wp, 4.2283_wp, 4.2161_wp &
        &, 4.2045_wp, 4.1841_wp, 4.3135_wp, 4.2562_wp, 4.2226_wp, 4.1667_wp, 4.0932_wp &
        &, 3.7134_wp, 3.6109_wp, 3.4962_wp, 3.8352_wp, 4.0688_wp, 4.0281_wp, 4.0099_wp &
        &, 4.3199_wp, 4.1188_wp, 3.9768_wp, 4.4192_wp, 4.5577_wp, 4.4516_wp, 4.3365_wp &
        &, 4.4058_wp, 4.3745_wp, 4.3539_wp, 2.8763_wp, 3.1294_wp, 3.5598_wp, 3.7465_wp &
        &, 3.5659_wp, 3.5816_wp, 3.3599_wp, 3.4024_wp, 3.1877_wp, 3.0484_wp, 3.7009_wp &
        &, 3.9451_wp, 3.8465_wp, 3.9873_wp, 3.7079_wp, 3.9083_wp, 3.6756_wp, 3.5150_wp &
        &, 4.0829_wp, 4.2780_wp, 4.1511_wp, 4.1260_wp, 4.0571_wp, 3.4865_wp, 3.9744_wp &
        &, 3.9150_wp, 3.8930_wp, 3.8578_wp, 3.8402_wp, 3.8073_wp, 3.7977_wp, 4.0036_wp &
        &, 3.7604_wp, 4.0288_wp, 3.8210_wp, 3.6757_wp, 4.2646_wp, 4.4558_wp, 4.2862_wp &
        &/)
    r0ab(4411:4465)=(/ &
        & 4.2122_wp, 3.7088_wp, 3.6729_wp, 3.5800_wp, 3.5276_wp, 3.5165_wp, 3.4783_wp &
        &, 3.4539_wp, 3.9553_wp, 3.9818_wp, 4.2040_wp, 3.9604_wp, 4.2718_wp, 4.0689_wp &
        &, 3.9253_wp, 4.4869_wp, 4.7792_wp, 4.4918_wp, 4.5342_wp, 4.5090_wp, 4.4868_wp &
        &, 4.4680_wp, 4.4486_wp, 4.4341_wp, 4.2023_wp, 4.3122_wp, 4.2710_wp, 4.3587_wp &
        &, 4.3407_wp, 4.3281_wp, 4.3174_wp, 4.1499_wp, 4.3940_wp, 4.3895_wp, 4.3260_wp &
        &, 4.2725_wp, 4.1961_wp, 3.7361_wp, 3.6193_wp, 3.4916_wp, 3.9115_wp, 3.9914_wp &
        &, 3.9809_wp, 3.9866_wp, 4.3329_wp, 4.1276_wp, 3.9782_wp, 4.5097_wp, 4.6769_wp &
        &, 4.5158_wp, 4.3291_wp, 4.3609_wp, 4.3462_wp, 4.3265_wp, 4.4341_wp &
        &/)

    k=0
    do i=1,max_elem
      do j=1,i
        k=k+1
        r(i,j)=r0ab(k)/autoang
        r(j,i)=r0ab(k)/autoang
      end do
    end do

  end subroutine setr0ab


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! Returns the number of a given element string (h-pu, 1-94)
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  elemental subroutine ELEM(KEY1, NAT)
    CHARACTER(*), intent(in) :: KEY1
    INTEGER, intent(out) :: NAT

    character(2), parameter :: ELEMNT(94) = [&
        & 'h ','he', &
        & 'li','be','b ','c ','n ','o ','f ','ne', &
        & 'na','mg','al','si','p ','s ','cl','ar', &
        & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu', &
        & 'zn','ga','ge','as','se','br','kr', &
        & 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag', &
        & 'cd','in','sn','sb','te','i ','xe', &
        & 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy', &
        & 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt', &
        & 'au','hg','tl','pb','bi','po','at','rn', &
        & 'fr','ra','ac','th','pa','u ','np','pu' ]

    CHARACTER(2) :: E
    integer :: k, j, n, i

    nat=0
    e=' '
    k=1
    do J=1,len(key1)
      if (k.gt.2)exit
      N=ICHAR(key1(J:J))
      if (n.ge.ichar('A') .and. n.le.ichar('Z') )then
        e(k:k)=char(n+ICHAR('a')-ICHAR('A'))
        k=k+1
      end if
      if (n.ge.ichar('a') .and. n.le.ichar('z') )then
        e(k:k)=key1(j:j)
        k=k+1
      end if
    end do

    do I=1,94
      if (e.eq.elemnt(i))then
        NAT=I
        RETURN
      end if
    end do


  end subroutine ELEM



  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! B E G I N O F P B C P A R T
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! compute coordination numbers by adding an inverse damping function
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine pbcncoord(natoms,rcov,iz,xyz,cn,lat,rep_cn,crit_cn)
    integer,intent(in) :: natoms,iz(*)
    real(wp),intent(in) :: rcov(94)
    real(wp), intent(in):: xyz(3,*),lat(3,3)
    real(wp), intent(in) :: crit_cn
    real(wp), intent(out):: cn(*)

    integer i,max_elem,rep_cn(3)

    integer iat,taux,tauy,tauz
    real(wp) dx,dy,dz,r,damp,xn,rr,rco,tau(3)

    do i=1,natoms
      xn=0.0d0
      do iat=1,natoms
        do taux=-rep_cn(1),rep_cn(1)
          do tauy=-rep_cn(2),rep_cn(2)
            do tauz=-rep_cn(3),rep_cn(3)
              if (iat.eq.i .and. taux.eq.0 .and. tauy.eq.0 .and.&
                  & tauz.eq.0) cycle
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
              dx=xyz(1,iat)-xyz(1,i)+tau(1)
              dy=xyz(2,iat)-xyz(2,i)+tau(2)
              dz=xyz(3,iat)-xyz(3,i)+tau(3)
              r=(dx*dx+dy*dy+dz*dz)
              if (r.gt.crit_cn) cycle
              r=sqrt(r)
              ! covalent distance in Bohr
              rco=rcov(iz(i))+rcov(iz(iat))
              rr=rco/r
              ! counting function exponential has a better long-range behavior than MH
              damp=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
              xn=xn+damp
              ! print '("cn(",I2,I2,"): ",E14.8)',i,iat,damp

            end do
          end do
        end do
      end do
      cn(i)=xn
    end do

  end subroutine pbcncoord

  subroutine pbcncoord_new(natoms,rcov,iz,xyz,cn,lat,rep_cn,crit_cn,xyz_hstep)
    integer,intent(in) :: natoms,iz(*)
    real(wp),intent(in) :: rcov(94)
    real(wp), intent(in):: xyz(3,*),lat(3,3)
    real(wp), intent(in) :: crit_cn
    real(wp), intent(out):: cn(*)
    real(wp), intent(in) :: xyz_hstep(3,*) ! displaced geometry
  
    integer i,max_elem,rep_cn(3)

    integer iat,taux,tauy,tauz
    real(wp) dx,dy,dz,r,damp,xn,rr,rco,tau(3)

    do i=1,natoms
      xn=0.0d0
      do iat=1,natoms
        do taux=-rep_cn(1),rep_cn(1)
          do tauy=-rep_cn(2),rep_cn(2)
            do tauz=-rep_cn(3),rep_cn(3)
              if (iat.eq.i .and. taux.eq.0 .and. tauy.eq.0 .and.&
                  & tauz.eq.0) cycle
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

              if(taux.eq.0.and.tauy.eq.0.and.tauz.eq.0.) then 
                dx=xyz_hstep(1,iat)-xyz_hstep(1,i) ! both in the unit cell, use the displaced geometry
                dy=xyz_hstep(2,iat)-xyz_hstep(2,i) 
                dz=xyz_hstep(3,iat)-xyz_hstep(3,i) 
              else
                dx=xyz(1,iat)-xyz_hstep(1,i)+tau(1) ! only iat in the unit cell, use the undisplaced geometry for i
                dy=xyz(2,iat)-xyz_hstep(2,i)+tau(2) 
                dz=xyz(3,iat)-xyz_hstep(3,i)+tau(3) 
              end if 

              r=(dx*dx+dy*dy+dz*dz)
              if (r.gt.crit_cn) cycle
              r=sqrt(r)
              ! covalent distance in Bohr
              rco=rcov(iz(i))+rcov(iz(iat))
              rr=rco/r
              ! counting function exponential has a better long-range behavior than MH
              damp=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
              xn=xn+damp
              ! print '("cn(",I2,I2,"): ",E14.8)',i,iat,damp

            end do
          end do
        end do
      end do
      cn(i)=xn
    end do

  end subroutine pbcncoord_new



  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! compute energy
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,&
      & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,&
      & e6,e8,e10,e12,e63,lat,rthr,rep_vdw,cn_thr,rep_cn)

    USE mp,           ONLY : mp_sum
    integer :: mykey, na_s, na_e

    integer, intent(in) :: max_elem,maxc
    real(wp), intent(in):: r2r4(max_elem),rcov(max_elem)
    real(wp), intent(in)::  rs6,rs8,rs10,alp6,alp8,alp10
    real(wp), intent(in):: rthr,cn_thr
    integer, intent(in):: rep_vdw(3),rep_cn(3)
    integer, intent(in):: n,iz(*),version,mxc(max_elem)
    ! integer rep_v(3)=rep_vdw!,rep_cn(3)
    real(wp), intent(in):: xyz(3,*),r0ab(max_elem,max_elem),lat(3,3)
    ! real(wp) rs6,rs8,rs10,alp6,alp8,alp10,rcov(max_elem)
    real(wp), intent(in):: c6ab(max_elem,max_elem,maxc,maxc,3)
    real(wp), intent(out) :: e6, e8, e10, e12, e63
    logical, intent(in):: noabc

    integer iat,jat,kat
    real(wp) :: crit_cn
    real(wp) r,r2,r6,r8,tmp,dx,dy,dz,c6,c8,c10,ang,rav,R0
    real(wp) damp6,damp8,damp10,rr,thr,c9,r42,c12,r10,c14
    real(wp) cn(n),rxyz(3),dxyz(3)
    real(wp) cc6ab(n*n),dmp(n*n),d2(3),t1,t2,t3,tau(3)
    integer ij,ik,jk
    integer taux,tauy,tauz,counter
    real(wp) a1,a2
    real(wp) bj_dmp6,bj_dmp8
    real(wp) tmp1,tmp2

    e6 =0
    e8 =0
    e10=0
    e12=0
    e63=0
    tau=(/0.0,0.0,0.0/)
    counter=0
    crit_cn=cn_thr
    cc6ab(:) = 0.0_wp
    ! Becke-Johnson parameters
    a1=rs6
    a2=rs8


    CALL block_distribute( n, me_dftd3, nproc_dftd3, na_s, na_e, mykey )
    IF ( mykey == 0 ) THEN


    ! DFT-D2
    if (version.eq.2)then


      do iat=na_s,min(na_e,n-1)
        do jat=iat+1,n
          c6=c6ab(iz(jat),iz(iat),1,1,1)
          do taux=-rep_vdw(1),rep_vdw(1)
            do tauy=-rep_vdw(2),rep_vdw(2)
              do tauz=-rep_vdw(3),rep_vdw(3)
                tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
                dx=xyz(1,iat)-xyz(1,jat)+tau(1)
                dy=xyz(2,iat)-xyz(2,jat)+tau(2)
                dz=xyz(3,iat)-xyz(3,jat)+tau(3)
                r2=dx*dx+dy*dy+dz*dz
                if (r2.gt.rthr) cycle
                r=sqrt(r2)
                damp6=1./(1.+exp(-alp6*(r/(rs6*r0ab(iz(jat),iz(iat)))-1.)))
                r6=r2**3
                e6 =e6+c6*damp6/r6
              end do
            end do
          end do
        end do
      end do

      do iat=na_s,na_e
        jat=iat
        c6=c6ab(iz(jat),iz(iat),1,1,1)
        do taux=-rep_vdw(1),rep_vdw(1)
          do tauy=-rep_vdw(2),rep_vdw(2)
            do tauz=-rep_vdw(3),rep_vdw(3)
              if (taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
              dx=tau(1)
              dy=tau(2)
              dz=tau(3)
              r2=dx*dx+dy*dy+dz*dz
              if (r2.gt.rthr) cycle
              r=sqrt(r2)
              damp6=1./(1.+exp(-alp6*(r/(rs6*r0ab(iz(jat),iz(iat)))-1.)))
              r6=r2**3
              e6 =e6+c6*damp6/r6*0.50d0
            end do
          end do
        end do
      end do



    else if ((version.eq.3).or.(version.eq.5)) then
      ! DFT-D3(zero-damping)

      call pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,crit_cn)

      do iat=na_s,min(na_e,n-1)
        do jat=iat+1,n
          ! get C6
          call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),&
              & cn(iat),cn(jat),c6)

          if (.not.noabc)then
            ij=lin(jat,iat)
            ! store C6 for C9, calc as sqrt
            cc6ab(ij)=sqrt(c6)
          end if
          do taux=-rep_vdw(1),rep_vdw(1)
            do tauy=-rep_vdw(2),rep_vdw(2)
              do tauz=-rep_vdw(3),rep_vdw(3)
                tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

                dx=xyz(1,iat)-xyz(1,jat)+tau(1)
                dy=xyz(2,iat)-xyz(2,jat)+tau(2)
                dz=xyz(3,iat)-xyz(3,jat)+tau(3)
                r2=dx*dx+dy*dy+dz*dz
                ! cutoff

                if (r2.gt.rthr) cycle
                r =sqrt(r2)
                R0=r0ab(iz(jat),iz(iat))
                rr=R0/r
                ! damping
                if(version.eq.3)then
                  ! DFT-D3 zero-damp
                  tmp=rs6*rr
                  damp6 =1.d0/( 1.d0+6.d0*tmp**alp6 )
                  tmp=rs8*rr
                  damp8 =1.d0/( 1.d0+6.d0*tmp**alp8 )
                else
                  ! DFT-D3M zero-damp
                  tmp=(r/(rs6*R0))+rs8*R0
                  damp6 =1.d0/( 1.d0+6.d0*tmp**(-alp6) )
                  tmp=(r/R0)+rs8*R0
                  damp8 =1.d0/( 1.d0+6.d0*tmp**(-alp8) )
                endif

                r6=r2**3
                e6 =e6+damp6/r6* c6
                ! write(*,*)'e6: ',c6*damp6/r6*autokcal

                ! stored in main as sqrt
                c8 =3.0d0*r2r4(iz(iat))*r2r4(iz(jat))*c6
                r8 =r6*r2

                e8 =e8+c8*damp8/r8

              end do
            end do
          end do
        end do
      end do

      do iat=na_s,na_e
        jat=iat
        ! get C6
        call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),&
            & cn(iat),cn(jat),c6)

        if (.not.noabc)then
          ij=lin(jat,iat)
          ! store C6 for C9, calc as sqrt
          cc6ab(ij)=sqrt(c6)
        end if
        do taux=-rep_vdw(1),rep_vdw(1)
          do tauy=-rep_vdw(2),rep_vdw(2)
            do tauz=-rep_vdw(3),rep_vdw(3)
              if (taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

              dx=tau(1)
              dy=tau(2)
              dz=tau(3)
              r2=dx*dx+dy*dy+dz*dz
              ! cutoff
              if (r2.gt.rthr) cycle
              r =sqrt(r2)
              R0=r0ab(iz(jat),iz(iat))
              rr=R0/r

              ! damping 
              if(version.eq.3)then
                ! DFT-D3 zero-damp
                tmp=rs6*rr
                damp6 =1.d0/( 1.d0+6.d0*tmp**alp6 )
                tmp=rs8*rr
                damp8 =1.d0/( 1.d0+6.d0*tmp**alp8 )
              else
                ! DFT-D3M zero-damp
                tmp=(r/(rs6*R0))+rs8*R0
                damp6 =1.d0/( 1.d0+6.d0*tmp**(-alp6) )
                tmp=(r/R0)+rs8*R0
                damp8 =1.d0/( 1.d0+6.d0*tmp**(-alp8) )
              endif


              r6=r2**3

              e6 =e6+damp6/r6*0.50d0 *C6

              ! stored in main as sqrt
              c8 =3.0d0*r2r4(iz(iat))*r2r4(iz(jat)) *C6
              r8 =r6*r2

              e8 =e8+c8*damp8/r8*0.50d0
              counter=counter+1

            end do
          end do
        end do
      end do
      ! write(*,*)'counter(edisp): ',counter
    else if((version.eq.4).or.(version.eq.6)) then


      ! DFT-D3(BJ-damping)
      call pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,crit_cn)

      do iat=na_s,na_e
        do jat=iat+1,n
          ! get C6
          call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),&
              & cn(iat),cn(jat),c6)

          rxyz=xyz(:,iat)-xyz(:,jat)
          r42=r2r4(iz(iat))*r2r4(iz(jat))
          bj_dmp6=(a1*dsqrt(3.0d0*r42)+a2)**6
          bj_dmp8=(a1*dsqrt(3.0d0*r42)+a2)**8

          if (.not.noabc)then
            ij=lin(jat,iat)
            ! store C6 for C9, calc as sqrt
            cc6ab(ij)=sqrt(c6)
          end if
          do taux=-rep_vdw(1),rep_vdw(1)
            do tauy=-rep_vdw(2),rep_vdw(2)
              do tauz=-rep_vdw(3),rep_vdw(3)
                tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

                dxyz=rxyz+tau

                r2=sum(dxyz*dxyz)
                ! cutoff
                if (r2.gt.rthr) cycle
                r =sqrt(r2)
                rr=r0ab(iz(jat),iz(iat))/r


                r6=r2**3

                e6 =e6+c6/(r6+bj_dmp6)

                ! stored in main as sqrt
                c8 =3.0d0*c6*r42
                r8 =r6*r2

                e8 =e8+c8/(r8+bj_dmp8)

                counter=counter+1

              end do
            end do
          end do
        end do

        ! Now the self interaction
        jat=iat
        ! get C6
        call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),&
            & cn(iat),cn(jat),c6)
        r42=r2r4(iz(iat))*r2r4(iz(iat))
        bj_dmp6=(a1*dsqrt(3.0d0*r42)+a2)**6
        bj_dmp8=(a1*dsqrt(3.0d0*r42)+a2)**8

        if (.not.noabc)then
          ij=lin(jat,iat)
          ! store C6 for C9, calc as sqrt
          cc6ab(ij)=dsqrt(c6)
        end if

        do taux=-rep_vdw(1),rep_vdw(1)
          do tauy=-rep_vdw(2),rep_vdw(2)
            do tauz=-rep_vdw(3),rep_vdw(3)
              if (taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

              r2=sum(tau*tau)
              ! cutoff
              if (r2.gt.rthr) cycle
              r =sqrt(r2)
              rr=r0ab(iz(jat),iz(iat))/r


              r6=r2**3

              e6 =e6+c6/(r6+bj_dmp6)*0.50d0

              ! stored in main as sqrt
              c8 =3.0d0*c6*r42
              r8 =r6*r2

              e8 =e8+c8/(r8+bj_dmp8)*0.50d0
              counter=counter+1

            end do
          end do
        end do
      end do


    end if


    ENDIF
    CALL mp_sum ( e6 , comm_dftd3 )
    CALL mp_sum ( e8 , comm_dftd3 )

    if (.not.noabc) then 

    ! compute non-additive third-order energy using averaged C6
    CALL mp_sum ( cc6ab , comm_dftd3 )
    call pbcthreebody(max_elem,xyz,lat,n,iz,rep_cn,crit_cn,&
        & cc6ab,r0ab,e63)

    end if 

  end subroutine pbcedisp


  subroutine pbcthreebody(max_elem,xyz,lat,n,iz,repv,cnthr,cc6ab,&
      & r0ab,eabc)

    USE mp,           ONLY : mp_sum
    integer :: mykey, na_s, na_smax, na_e

    integer max_elem
    INTEGER :: n,i,j,k,jtaux,jtauy,jtauz,iat,jat,kat
    INTEGER :: ktaux,ktauy,ktauz,counter,ij,ik,jk,idum
    REAL(WP) :: dx,dy,dz,rij2,rik2,rjk2,c9,rr0ij,rr0ik
    REAL(WP) :: rr0jk,geomean,fdamp,rik,rjk,rij
    REAL(WP) :: r0ij,r0ik,r0jk
    REAL(WP),INTENT(OUT)::eabc
    REAL(WP) :: tmp,tmp1,tmp2,tmp3,tmp4,ang

    REAL(WP) ,DIMENSION(3,3),INTENT(IN)::lat
    REAL(WP) ,DIMENSION(3,*),INTENT(IN) :: xyz
    INTEGER,DIMENSION(*),INTENT(IN)::iz
    REAL(WP),DIMENSION(3):: jtau,ktau,jxyz,kxyz,ijvec,ikvec,jkvec,dumvec
    INTEGER,DIMENSION(3):: repv
    REAL(WP),INTENT(IN) ::cnthr
    REAL(WP),DIMENSION(n*n),INTENT(IN)::cc6ab
    REAL(WP),DIMENSION(max_elem,max_elem),INTENT(IN):: r0ab
    REAL(WP),PARAMETER::sr9=0.75d0
    REAL(WP),PARAMETER::alp9=-16.0d0
    REAL(WP) :: abcthr
    INTEGER :: repmin1, repmin2, repmin3, repmax1, repmax2, repmax3    
    INTEGER :: repv1, repv2, repv3 
    REAL(WP)  :: ijvec1, ijvec2, ijvec3
    REAL(WP)  :: ikvec1, ikvec2, ikvec3
    REAL(WP)  :: jkvec1, jkvec2, jkvec3
    REAL(WP)  :: jtau1, jtau2, jtau3
    REAL(WP)  :: ktau1, ktau2, ktau3
    REAL(WP)  :: dumvec11, dumvec12, dumvec13  
    REAL(WP)  :: dumvec21, dumvec22, dumvec23  
    ! REAL(WP) :: time1,time2

    counter=0
    eabc=0.0d0
    abcthr=cnthr
    ! abcthr=1.0d99
    ! write(*,*)'thr:',(abcthr)

    ! call cpu_time(time1)

    repv1 = repv(1)
    repv2 = repv(2)
    repv3 = repv(3)

    CALL block_distribute( n, me_dftd3, nproc_dftd3, na_s, na_e, mykey )
    IF ( mykey == 0 ) THEN
    na_smax = max(3,na_s)

!$acc data copyin(xyz(1:3,1:n),iz(1:n),cc6ab(1:n*n),lat(1:3,1:3),r0ab(1:max_elem,1:max_elem)) 
!$acc kernels  vector_length(32)
!$acc loop collapse(2) gang private(ijvec1,ijvec2,ijvec3, ikvec1,ikvec2,ikvec3, jkvec1,jkvec2,jkvec3, c9, r0ij,r0ik,r0jk, &
!$acc&                              repmin1,repmin2,repmin3, repmax1,repmax2,repmax3, jtau1,jtau2,jtau3, &
!$acc&                              dumvec11,dumvec12,dumvec13,rij2,rr0ij)  reduction(+:eabc)
    do iat=na_smax,na_e
     do jat = 2, n
        if(jat.ge.iat) cycle
        ijvec1=xyz(1,jat)-xyz(1,iat)
        ijvec2=xyz(2,jat)-xyz(2,iat)
        ijvec3=xyz(3,jat)-xyz(3,iat)
        ij=lin(iat,jat)
        r0ij=r0ab(iz(iat),iz(jat))
!$acc loop seq
        do kat=1,jat-1
          ik=lin(iat,kat)
          jk=lin(jat,kat)
          ikvec1=xyz(1,kat)-xyz(1,iat)
          ikvec2=xyz(2,kat)-xyz(2,iat)
          ikvec3=xyz(3,kat)-xyz(3,iat)
          jkvec1=xyz(1,kat)-xyz(1,jat)
          jkvec2=xyz(2,kat)-xyz(2,jat)
          jkvec3=xyz(3,kat)-xyz(3,jat)
          c9=-1.0d0*(cc6ab(ij)*cc6ab(ik)*cc6ab(jk))

          r0ik=r0ab(iz(iat),iz(kat))
          r0jk=r0ab(iz(jat),iz(kat))

          do jtaux=-repv1,repv1
            do jtauy=-repv2,repv2
              do jtauz=-repv3,repv3
                repmin1=max(-repv1,jtaux-repv1)
                repmax1=min(repv1,jtaux+repv1)
                repmin2=max(-repv2,jtauy-repv2)
                repmax2=min(repv2,jtauy+repv2)
                repmin3=max(-repv3,jtauz-repv3)
                repmax3=min(repv3,jtauz+repv3)
                jtau1=jtaux*lat(1,1)+jtauy*lat(1,2)+jtauz*lat(1,3)
                jtau2=jtaux*lat(2,1)+jtauy*lat(2,2)+jtauz*lat(2,3)
                jtau3=jtaux*lat(3,1)+jtauy*lat(3,2)+jtauz*lat(3,3)
                dumvec11=(ijvec1+jtau1)*(ijvec1+jtau1)
                dumvec12=(ijvec2+jtau2)*(ijvec2+jtau2)
                dumvec13=(ijvec3+jtau3)*(ijvec3+jtau3)
                rij2=dumvec11+dumvec12+dumvec13
                if (rij2.gt.abcthr)cycle

                rr0ij=DSQRT(rij2)/r0ij

!$acc loop vector collapse(3) private(ktau1,ktau2,ktau3,dumvec21,dumvec22,dumvec23,rik2,rr0ik,rjk2,rr0jk,geomean,fdamp,tmp1,tmp2,tmp3,tmp4,ang) reduction(+:eabc)
                do ktaux=repmin1,repmax1
                  do ktauy=repmin2,repmax2
                    do ktauz=repmin3,repmax3
                      ktau1=ktaux*lat(1,1)+ktauy*lat(1,2)+ktauz*lat(1,3)
                      ktau2=ktaux*lat(2,1)+ktauy*lat(2,2)+ktauz*lat(2,3)
                      ktau3=ktaux*lat(3,1)+ktauy*lat(3,2)+ktauz*lat(3,3)
                      dumvec21=(ikvec1+ktau1)*(ikvec1+ktau1)
                      dumvec22=(ikvec2+ktau2)*(ikvec2+ktau2)
                      dumvec23=(ikvec3+ktau3)*(ikvec3+ktau3)
                      rik2=dumvec21+dumvec22+dumvec23
                      if (rik2.gt.abcthr)cycle
                      rr0ik=DSQRT(rik2)/r0ik
                      dumvec21=jkvec1+ktau1-jtau1
                      dumvec22=jkvec2+ktau2-jtau2
                      dumvec23=jkvec3+ktau3-jtau3
                      rjk2=dumvec21*dumvec21+dumvec22*dumvec22+dumvec23*dumvec23
                      if (rjk2.gt.abcthr)cycle
                      rr0jk=DSQRT(rjk2)/r0jk


                      geomean=(rr0ij*rr0ik*rr0jk)**(1.0d0/3.0d0)
                      ! write(*,*)'geomean:',geomean
                      fdamp=1./(1.+6.*(sr9*geomean)**alp9)
                      tmp1 = (rij2+rjk2-rik2)
                      tmp2 = (rij2+rik2-rjk2)
                      tmp3 = (rik2+rjk2-rij2)
                      tmp4=rij2*rjk2*rik2
                      ang=(0.375d0*tmp1*tmp2*tmp3/tmp4+1.0d0)/tmp4**1.5d0

                      eabc=eabc+ang*c9*fdamp

                    end do
                  end do
                end do

              end do
            end do
          end do

        end do
      end do
    end do
!$acc end kernels  
!$acc end data 

    do iat=max(2,na_s),na_e
      jat=iat
      ij=lin(iat,jat)
      ijvec=0.0d0
      r0ij=r0ab(iz(iat),iz(jat))
      do kat=1,iat-1
        jk=lin(jat,kat)
        ik=jk
        ikvec=xyz(:,kat)-xyz(:,iat)
        jkvec=ikvec
        c9=-(cc6ab(ij)*cc6ab(ik)*cc6ab(jk))

        r0ik=r0ab(iz(iat),iz(kat))
        r0jk=r0ab(iz(jat),iz(kat))
        do jtaux=-repv1,repv1
          do jtauy=-repv2,repv2
            do jtauz=-repv3,repv3
              repmin1=max(-repv1,jtaux-repv1)
              repmax1=min(repv1,jtaux+repv1)
              repmin2=max(-repv2,jtauy-repv2)
              repmax2=min(repv2,jtauy+repv2)
              repmin3=max(-repv3,jtauz-repv3)
              repmax3=min(repv3,jtauz+repv3)
              if (jtaux.eq.0 .and. jtauy.eq.0 .and. jtauz.eq.0) cycle
              jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
              dumvec=ijvec+jtau
              dumvec=dumvec*dumvec
              rij2=SUM(dumvec)
              if (rij2.gt.abcthr)cycle

              rr0ij=DSQRT(rij2)/r0ij

              do ktaux=repmin1,repmax1
                do ktauy=repmin2,repmax2
                  do ktauz=repmin3,repmax3
                    ! every result * 0.5
                    ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
                    dumvec=ikvec+ktau
                    dumvec=dumvec*dumvec
                    rik2=SUM(dumvec)
                    if (rik2.gt.abcthr)cycle
                    rr0ik=DSQRT(rik2)/r0ik

                    dumvec=jkvec+ktau-jtau
                    dumvec=dumvec*dumvec
                    rjk2=SUM(dumvec)
                    if (rjk2.gt.abcthr)cycle
                    rr0jk=DSQRT(rjk2)/r0jk


                    geomean=(rr0ij*rr0ik*rr0jk)**(1./3.)
                    fdamp=1./(1.+6.*(sr9*geomean)**alp9)
                    tmp1 = (rij2+rjk2-rik2)
                    tmp2 = (rij2+rik2-rjk2)
                    tmp3 = (rik2+rjk2-rij2)
                    tmp4=rij2*rjk2*rik2
                    ang=(0.375d0*tmp1*tmp2*tmp3/tmp4+1.0d0)/tmp4**1.5d0

                    eabc=eabc+c9*fdamp*ang/2.0
                  end do
                end do
              end do

            end do
          end do
        end do
      end do
    end do

    do iat=max(2,na_s),na_e
      do jat=1,iat-1
        kat=jat
        ij=lin(iat,jat)
        jk=lin(jat,kat)
        ik=ij
        ikvec=xyz(:,kat)-xyz(:,iat)
        ijvec=ikvec
        jkvec=0.0d0
        c9=-(cc6ab(ij)*cc6ab(ik)*cc6ab(jk))

        r0ij=r0ab(iz(iat),iz(jat))
        r0ik=r0ij
        r0jk=r0ab(iz(jat),iz(kat))

        do jtaux=-repv1,repv1
          do jtauy=-repv2,repv2
            do jtauz=-repv3,repv3
              repmin1=max(-repv1,jtaux-repv1)
              repmax1=min(repv1,jtaux+repv1)
              repmin2=max(-repv2,jtauy-repv2)
              repmax2=min(repv2,jtauy+repv2)
              repmin3=max(-repv3,jtauz-repv3)
              repmax3=min(repv3,jtauz+repv3)
              jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
              dumvec=ijvec+jtau
              dumvec=dumvec*dumvec
              rij2=SUM(dumvec)
              if (rij2.gt.abcthr)cycle

              rr0ij=DSQRT(rij2)/r0ij

              do ktaux=repmin1,repmax1
                do ktauy=repmin2,repmax2
                  do ktauz=repmin3,repmax3
                    ! every result * 0.5
                    if (jtaux.eq.ktaux .and. jtauy.eq.ktauy&
                        & .and. jtauz.eq.ktauz) cycle
                    ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
                    dumvec=ikvec+ktau
                    dumvec=dumvec*dumvec
                    rik2=SUM(dumvec)
                    if (rik2.gt.abcthr)cycle
                    rr0ik=DSQRT(rik2)/r0ik

                    dumvec=jkvec+ktau-jtau
                    dumvec=dumvec*dumvec
                    rjk2=SUM(dumvec)
                    if (rjk2.gt.abcthr)cycle
                    rr0jk=DSQRT(rjk2)/r0jk


                    geomean=(rr0ij*rr0ik*rr0jk)**(1./3.)
                    fdamp=1./(1.+6.*(sr9*geomean)**alp9)
                    tmp1 = (rij2+rjk2-rik2)
                    tmp2 = (rij2+rik2-rjk2)
                    tmp3 = (rik2+rjk2-rij2)
                    tmp4=rij2*rjk2*rik2
                    ang=(0.375d0*tmp1*tmp2*tmp3/tmp4+1.0d0)/tmp4**1.5d0

                    eabc=eabc+c9*fdamp*ang/2.0
                  end do
                end do
              end do

            end do
          end do
        end do
      end do
    end do


    ! And finally the self interaction iat=jat=kat all

    idum=0
    do iat=na_s,na_e
      jat=iat
      kat=iat
      ijvec=0.0d0
      ij=lin(iat,iat)
      ik=ij
      jk=ij
      ikvec=ijvec
      jkvec=ikvec
      c9=-(cc6ab(ij)*cc6ab(ik)*cc6ab(jk))

      r0ij=r0ab(iz(iat),iz(iat))
      r0ik=r0ij
      r0jk=r0ij
      do jtaux=-repv1,repv1
        do jtauy=-repv2,repv2
          do jtauz=-repv3,repv3
            repmin1=max(-repv1,jtaux-repv1)
            repmax1=min(repv1,jtaux+repv1)
            repmin2=max(-repv2,jtauy-repv2)
            repmax2=min(repv2,jtauy+repv2)
            repmin3=max(-repv3,jtauz-repv3)
            repmax3=min(repv3,jtauz+repv3)
            if (jtaux.eq.0 .and. jtauy.eq.0 .and. jtauz.eq.0) cycle
            jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
            dumvec=jtau
            dumvec=dumvec*dumvec
            rij2=SUM(dumvec)
            if (rij2.gt.abcthr)cycle
            rr0ij=DSQRT(rij2)/r0ij

            do ktaux=repmin1,repmax1
              do ktauy=repmin2,repmax2
                do ktauz=repmin3,repmax3
                  if ((ktaux.eq.0) .and.( ktauy.eq.0) .and.( ktauz.eq.0))cycle
                  if ((ktaux.eq.jtaux) .and. (ktauy.eq.jtauy)&
                      & .and. (ktauz.eq.jtauz)) cycle

                  ! every result * 1/6 becaues every triple is counted twice due to the tw
                  !
                  !plus 1/3 becaues every triple is three times in each unitcell
                  ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
                  dumvec=ktau
                  dumvec=dumvec*dumvec
                  rik2=SUM(dumvec)
                  if (rik2.gt.abcthr)cycle
                  rr0ik=DSQRT(rik2)/r0ik

                  dumvec=jkvec+ktau-jtau
                  dumvec=dumvec*dumvec
                  rjk2=SUM(dumvec)
                  if (rjk2.gt.abcthr)cycle
                  rr0jk=DSQRT(rjk2)/r0jk


                  geomean=(rr0ij*rr0ik*rr0jk)**(1./3.)
                  fdamp=1./(1.+6.*(sr9*geomean)**alp9)
                  tmp1 = (rij2+rjk2-rik2)
                  tmp2 = (rij2+rik2-rjk2)
                  tmp3 = (rik2+rjk2-rij2)
                  tmp4=rij2*rjk2*rik2
                  ang=(0.375d0*tmp1*tmp2*tmp3/tmp4+1.0d0)/tmp4**1.5d0

                  eabc=eabc+c9*fdamp*ang/6.0d0

                end do
              end do
            end do
          end do
        end do
      end do

    end do

   ENDIF
   CALL mp_sum ( eabc , comm_dftd3 )

  end subroutine pbcthreebody


  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! compute gradient
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  subroutine pbcgdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,&
      & rcov,s6,s18,rs6,rs8,rs10,alp6,alp8,alp10,noabc,num,&
      & version,g,disp,gnorm,stress,lat,rep_v,rep_cn,&
      & crit_vdw,echo,crit_cn)


    USE mp,           ONLY : mp_sum
    integer :: mykey, na_s, na_smax, na_e

    integer n,iz(*),max_elem,maxc,version,mxc(max_elem)
    real(wp) xyz(3,*),r0ab(max_elem,max_elem),r2r4(*)
    real(wp) c6ab(max_elem,max_elem,maxc,maxc,3)
    real(wp) g(3,*),s6,s18,rcov(max_elem)
    real(wp) rs6,rs8,rs10,alp10,alp8,alp6
    real(wp) a1,a2
    real(wp) bj_dmp6,bj_dmp8
    logical noabc,num,echo
    ! conversion factors

    integer iat,jat,i,j,kat,my,ny,a,b,idum,tau2
    real(wp) R0,C6,alp,R42,disp,x1,y1,z1,x2,y2,z2,rr,e6abc,fdum
    real(wp) dx,dy,dz,r2,r,r4,r6,r8,r10,r12,t6,t8,t10,damp1
    real(wp) damp6,damp8,damp9,e6,e8,e10,e12,gnorm,tmp1
    real(wp) s10,s8,gC6(3),term,step,dispr,displ,r235,tmp2
    real(wp) cn(n),gx1,gy1,gz1,gx2,gy2,gz2,rthr,testsum
    real(wp), DIMENSION(3,3) :: lat,stress,sigma,virialstress,lat_1
    real(wp), DIMENSION(3,3) :: gC6_stress
    integer, DIMENSION(3) :: rep_v,rep_cn
    real(wp) crit_vdw,crit_cn
    integer taux,tauy,tauz
    real(wp), DIMENSION(3) :: tau,vec12,dxyz,dxyz0
    real(wp) ::outpr(3,3)
    real(wp), DIMENSION(3,3):: outerprod

    real(wp) rij(3),rik(3),rjk(3),r7,r9
    real(wp) rik_dist,rjk_dist
    real(wp) drik,drjk
    real(wp) rcovij
    real(wp) dc6,c6chk
    real(wp) expterm,dcni
    real(wp), allocatable,dimension(:,:,:,:) :: drij
    real(wp), allocatable,dimension(:,:,:,:) :: dcn
    real(wp) dcnn
    real(wp) :: dc6_rest
    real(wp) vec(3),vec2(3),dummy
    real(wp) dc6i(n)
    real(wp) dc6ij(n,n)
    real(wp) dc6_rest_sum(n*(n+1)/2)
    integer linij,linik,linjk
    real(wp) abc(3,n)

    real(wp) eabc
    real(wp) gabc(3,n),glatabc(3,3)
    real(wp) sigma_abc(3,3)
    real(wp) labc,rabc
    real(wp) ,dimension(3) ::ijvec,ikvec,jkvec,jtau,ktau,dumvec
    integer jtaux,jtauy,jtauz,ktaux,ktauy,ktauz,mtaux,mtauy,mtauz
    integer,dimension(3) :: taumin,taumax
    integer mat,linim,linjm,linkm
    real(wp) rij2,rik2,rjk2,c9,c6ij,c6ik,c6jk,geomean,geomean3
    real(wp) rr0ij,rr0jk,rr0ik,dc6iji,dc6ijj
    real(wp) :: sr9=0.75d0
    real(wp), parameter :: alp9=-16.0d0
    real(wp),DIMENSION(n*(n+1)) ::c6save
    real(wp) abcthr,time1,time2,geomean2,r0av,dc9,dfdmp,dang,ang
    integer,dimension(3) ::repv,repmin,repmax
    integer :: rep_v1, rep_v2, rep_v3  
    integer :: rep_cn1, rep_cn2, rep_cn3  
    integer :: repmin1, repmin2, repmin3   
    integer :: repmax1, repmax2, repmax3   
    real(wp) :: dumvec1, dumvec2, dumvec3  
    real(wp) :: ijvec1, ijvec2, ijvec3   
    real(wp) :: ikvec1, ikvec2, ikvec3   
    real(wp) :: jkvec1, jkvec2, jkvec3   
    real(wp) :: jtau1, jtau2, jtau3  
    real(wp) :: ktau1, ktau2, ktau3  

    ! R^2 cut-off
    rthr=crit_vdw
    abcthr=crit_cn
    ! write(*,*)'abcthr:', abcthr**(1./1.)
    sigma=0.0d0
    virialstress=0.0d0
    stress=0.0d0
    gabc=0.0d0
    glatabc=0.0d0

    ! testsum=0.0d0

    if (echo)write(*,*)

    if (num) then
      if (echo) &
          & write(*,*) 'doing numerical gradient O(N^3) ...'

      call pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,&
          & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,&
          & e6,e8,e10,e12,e6abc,lat,rthr,rep_v,crit_cn,rep_cn)


      disp=-s6*e6-s18*e8-e6abc

      step=2.d-5

      do i=1,n
        do j=1,3
          xyz(j,i)=xyz(j,i)+step
          call pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,&
              & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,&
              & e6,e8,e10,e12,e6abc,lat,rthr,rep_v,crit_cn,rep_cn)

          dispr=-s6*e6-s18*e8-e6abc
          rabc=e6abc
          xyz(j,i)=xyz(j,i)-2*step
          call pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,&
              & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,&
              & e6,e8,e10,e12,e6abc,lat,rthr,rep_v,crit_cn,rep_cn)

          displ=-s6*e6-s18*e8-e6abc
          labc=e6abc
          gabc(j,i)=0.5*(rabc-labc)/step
          g(j,i)=0.5*(dispr-displ)/step
          xyz(j,i)=xyz(j,i)+step
        end do
      end do
      if (echo) write(*,*)'Doing numerical stresstensor...'

      call xyz_to_abc(xyz,abc,lat,n)
      step=2.d-5
      if (echo) write(*,*)'step: ',step
      do i=1,3
        do j=1,3
          lat(j,i)=lat(j,i)+step
          call abc_to_xyz(abc,xyz,lat,n)
          call pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,&
              & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,&
              & e6,e8,e10,e12,e6abc,lat,rthr,rep_v,crit_cn,rep_cn)

          dispr=-s6*e6-s18*e8-e6abc
          labc=e6abc


          lat(j,i)=lat(j,i)-2*step
          call abc_to_xyz(abc,xyz,lat,n)
          call pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,&
              & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,&
              & e6,e8,e10,e12,e6abc,lat,rthr,rep_v,crit_cn,rep_cn)

          displ=-s6*e6-s18*e8-e6abc
          rabc=e6abc
          stress(j,i)=(dispr-displ)/(step*2.0)
          glatabc(j,i)=(rabc-labc)/(step*2.0d0)

          lat(j,i)=lat(j,i)+step
          call abc_to_xyz(abc,xyz,lat,n)

        end do
      end do

      sigma=0.0d0
      call inv_cell(lat,lat_1)
      do a=1,3
        do b=1,3
          do my=1,3
            sigma(a,b)=sigma(a,b)-stress(a,my)*lat(b,my)
          end do
        end do
      end do

      goto 999

    end if


    if (version.eq.2)then
      if (echo)write(*,*) 'doing analytical gradient D-old O(N^2) ...'
      disp=0
      stress=0.0d0
      do iat=1,n-1
        do jat=iat+1,n
          R0=r0ab(iz(jat),iz(iat))*rs6
          c6=c6ab(iz(jat),iz(iat),1,1,1)*s6
          do taux=-rep_v(1),rep_v(1)
            do tauy=-rep_v(2),rep_v(2)
              do tauz=-rep_v(3),rep_v(3)
                tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
                dxyz=xyz(:,iat)-xyz(:,jat)+tau
                r2 =sum(dxyz*dxyz)
                if (r2.gt.rthr) cycle
                r235=r2**3.5
                r =dsqrt(r2)
                damp6=exp(-alp6*(r/R0-1.0d0))
                damp1=1.+damp6
                tmp1=damp6/(damp1*damp1*r235*R0)
                tmp2=6./(damp1*r*r235)

                term=alp6*tmp1-tmp2
                g(:,iat)=g(:,iat)-term*dxyz*c6
                g(:,jat)=g(:,jat)+term*dxyz*c6
                disp=disp+c6*(1./damp1)/r2**3

                do ny=1,3
                  do my=1,3
                    sigma(my,ny)=sigma(my,ny)+term*dxyz(ny)*dxyz(my)*c6
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      ! and now the self interaction, only for convenient energy in dispersion
      do iat=1,n
        jat=iat
        R0=r0ab(iz(jat),iz(iat))*rs6
        c6=c6ab(iz(jat),iz(iat),1,1,1)*s6
        do taux=-rep_v(1),rep_v(1)
          do tauy=-rep_v(2),rep_v(2)
            do tauz=-rep_v(3),rep_v(3)
              if (taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0) cycle
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

              dxyz=tau
              ! vec12=(/ dx,dy,dz /)
              r2 =sum(dxyz*dxyz)
              if (r2.gt.rthr) cycle
              r235=r2**3.5
              r =dsqrt(r2)
              damp6=exp(-alp6*(r/R0-1.0d0))
              damp1=1.+damp6
              tmp1=damp6/(damp1*damp1*r235*R0)
              tmp2=6./(damp1*r*r235)
              disp=disp+(c6*(1./damp1)/r2**3)*0.50d0
              term=alp6*tmp1-tmp2
              do ny=1,3
                do my=1,3
                  sigma(my,ny)=sigma(my,ny)+term*dxyz(ny)*dxyz(my)*c6*0.5d0
                end do
              end do


            end do
          end do
        end do
      end do

      call inv_cell(lat,lat_1)
      do a=1,3
        do b=1,3
          do my=1,3
            stress(a,b)=stress(a,b)-sigma(a,my)*lat_1(b,my)
          end do
        end do
      end do

      disp=-disp
      ! sigma=virialstress
      goto 999
    end if

    CALL block_distribute( n, me_dftd3, nproc_dftd3, na_s, na_e, mykey )

    if ((version.eq.3).or.(version.eq.5)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! begin ZERO DAMPING GRADIENT
      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (echo)&
          & write(*,*) 'doing analytical gradient O(N^2) ...'
      ! precompute for analytical part
      call pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,crit_cn)


      s8 =s18
      s10=s18
      allocate(drij(-rep_v(3):rep_v(3),-rep_v(2):rep_v(2),&
          & -rep_v(1):rep_v(1),n*(n+1)/2))

      disp=0

      drij=0.0d0
      dc6_rest=0.0d0
      dc6_rest_sum=0.0d0
      c6save=0.0d0
      kat=0
      dc6i=0.0d0
      dc6ij=0.0d0

      IF ( mykey == 0 ) THEN

      do iat=na_s, na_e
        call get_dC6_dCNij(maxc,max_elem,c6ab,mxc(iz(iat)),&
            & mxc(iz(iat)),cn(iat),cn(iat),iz(iat),iz(iat),iat,iat,&
            & c6,dc6iji,dc6ijj)

        c6save(lin(iat,iat))=c6
        dc6ij(iat,iat)=dc6iji
        r0=r0ab(iz(iat),iz(iat))
        r42=r2r4(iz(iat))*r2r4(iz(iat))
        rcovij=rcov(iz(iat))+rcov(iz(iat))


        do taux=-rep_v(1),rep_v(1)
          do tauy=-rep_v(2),rep_v(2)
            do tauz=-rep_v(3),rep_v(3)
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)


              !first dE/d(tau) saved in drij(i,i,counter)
              rij=tau
              r2=sum(rij*rij)
              ! if (r2.gt.rthr) cycle

              if (r2.gt.0.1.and.r2.lt.rthr) then


                r=dsqrt(r2)
                r6=r2*r2*r2
                r7=r6*r
                r8=r6*r2
                r9=r8*r

                !
                ! Calculates damping functions:
                if (version.eq.3) then
                  t6 = (r/(rs6*R0))**(-alp6)
                  damp6 =1.d0/( 1.d0+6.d0*t6 )
                  t8 = (r/(rs8*R0))**(-alp8)
                  damp8 =1.d0/( 1.d0+6.d0*t8 )

                  drij(tauz,tauy,taux,lin(iat,iat))=drij(tauz,tauy,taux,lin(iat,&
                      & iat))&
                      & +(-s6*(6.0/(r7)*C6*damp6)&
                      & -s8*(24.0/(r9)*C6*r42*damp8))*0.5d0


                  drij(tauz,tauy,taux,lin(iat,iat))=drij(tauz,tauy,taux,lin(iat,&
                      & iat))&
                      & +(s6*C6/r7*6.d0*alp6*t6*damp6*damp6&
                      & +s8*C6*r42/r9*18.d0*alp8*t8*damp8*damp8)*0.5d0
                else !version.eq.5
                  t6 = (r/(rs6*R0)+R0*rs8)**(-alp6)
                  damp6 =1.d0/( 1.d0+6.d0*t6 )
                  t8 = (r/(R0)+R0*rs8)**(-alp8)
                  damp8 =1.d0/( 1.d0+6.d0*t8 )
  
                  tmp1=s6*6.d0*damp6*C6/r7
                  tmp2=s8*6.d0*C6*r42*damp8/r9
                  drij(tauz,tauy,taux,lin(iat,iat))=drij(tauz,tauy,taux,lin(iat, &
                      & iat)) - (tmp1 +4.d0*tmp2)*0.5d0               ! d(r^(-6))/d(r_ij)
  
  
                  drij(tauz,tauy,taux,lin(iat,iat))=drij(tauz,tauy,taux,lin(iat, &
                      & iat)) +(tmp1*alp6*t6*damp6*r/(r+rs6*R0*R0*rs8) &
                      & +3.d0*tmp2*alp8*t8*damp8*r/(r+R0*R0*rs8))*0.5d0  !d(f_dmp)/d(r_ij)
                endif
                !
                ! in dC6_rest all terms BUT C6-term is saved for the kat-loop
                !
                dc6_rest=&
                    & (s6/r6*damp6+3.d0*s8*r42/r8*damp8)*0.50d0


                disp=disp-dc6_rest*c6

                dc6i(iat)=dc6i(iat)+dc6_rest*(dc6iji+dc6ijj)
                ! if (r2.lt.crit_cn)
                dc6_rest_sum(lin(iat,iat))=dc6_rest_sum(lin(iat,iat))+dc6_rest


              else
                drij(tauz,tauy,taux,lin(iat,iat))=0.0d0
              end if


            end do
          end do
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! B E G I N jat L O O P
!!!!!!!!!!!!!!!!!!!!!!!!!!
        do jat=1,iat-1
          !
          ! get_dC6_dCNij calculates the derivative dC6(iat,jat)/dCN(iat) and
          ! dC6(iat,jat)/dCN(jat). these are saved in dC6ij for the kat loop
          !
          call get_dC6_dCNij(maxc,max_elem,c6ab,mxc(iz(iat)),&
              & mxc(iz(jat)),cn(iat),cn(jat),iz(iat),iz(jat),iat,jat,&
              & c6,dc6iji,dc6ijj)

          r0=r0ab(iz(jat),iz(iat))
          r42=r2r4(iz(iat))*r2r4(iz(jat))
          rcovij=rcov(iz(iat))+rcov(iz(jat))
          linij=lin(iat,jat)

          dc6ij(iat,jat)=dc6iji
          dc6ij(jat,iat)=dc6ijj
          c6save(linij)=c6
          do taux=-rep_v(1),rep_v(1)
            do tauy=-rep_v(2),rep_v(2)
              do tauz=-rep_v(3),rep_v(3)
                tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)


                rij=xyz(:,jat)-xyz(:,iat)+tau
                r2=sum(rij*rij)
                if (r2.gt.rthr) cycle


                r=dsqrt(r2)
                r6=r2*r2*r2
                r7=r6*r
                r8=r6*r2
                r9=r8*r

                !
                ! Calculates damping functions:
                if (version.eq.3) then
                  t6 = (r/(rs6*R0))**(-alp6)
                  damp6 =1.d0/( 1.d0+6.d0*t6 )
                  t8 = (r/(rs8*R0))**(-alp8)
                  damp8 =1.d0/( 1.d0+6.d0*t8 )

                  drij(tauz,tauy,taux,linij)=drij(tauz,tauy,taux,&
                      & linij)&
                      & -s6*(6.0/(r7)*C6*damp6)&
                      & -s8*(24.0/(r9)*C6*r42*damp8)

                  drij(tauz,tauy,taux,linij)=drij(tauz,tauy,taux,&
                      & linij)&
                      & +s6*C6/r7*6.d0*alp6*t6*damp6*damp6&
                      & +s8*C6*r42/r9*18.d0*alp8*t8*damp8*damp8
                else !version.eq.5
                  t6 = (r/(rs6*R0)+R0*rs8)**(-alp6)
                  damp6 =1.d0/( 1.d0+6.d0*t6 )
                  t8 = (r/(R0)+R0*rs8)**(-alp8)
                  damp8 =1.d0/( 1.d0+6.d0*t8 )
  
                  tmp1=s6*6.d0*damp6*C6/r7
                  tmp2=s8*6.d0*C6*r42*damp8/r9
                  drij(tauz,tauy,taux,linij)=drij(tauz,tauy,taux, &
                      & linij) - (tmp1 +4.d0*tmp2)  ! d(r^(-6))/d(r_ij)
  
  
                  drij(tauz,tauy,taux,linij)=drij(tauz,tauy,taux,linij) &
                      & +(tmp1*alp6*t6*damp6*r/(r+rs6*R0*R0*rs8) & 
                      & +3.d0*tmp2*alp8*t8*damp8*r/(r+R0*R0*rs8)) !d(f_dmp)/d(r_ij)
                endif
                !
                ! in dC6_rest all terms BUT C6-term is saved for the kat-loop
                !
                dc6_rest=&
                    & (s6/r6*damp6+3.d0*s8*r42/r8*damp8)


                disp=disp-dc6_rest*c6

                dc6i(iat)=dc6i(iat)+dc6_rest*dc6iji
                dc6i(jat)=dc6i(jat)+dc6_rest*dc6ijj
                ! if (r2.lt.crit_cn)
                dc6_rest_sum(linij)=dc6_rest_sum(linij)&
                    & +dc6_rest


              end do
            end do
          end do

        end do

      end do
      END IF

    elseif ((version.eq.4).or.(version.eq.6)) then



!!!!!!!!!!!!!!!!!!!!!!!
      ! NOW THE BJ Gradient !
!!!!!!!!!!!!!!!!!!!!!!!


      if (echo) write(*,*) 'doing analytical gradient O(N^2) ...'
      call pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,crit_cn)

      a1 =rs6
      a2 =rs8
      s8 =s18

      allocate(drij(-rep_v(3):rep_v(3),-rep_v(2):rep_v(2),&
          & -rep_v(1):rep_v(1),n*(n+1)/2))
      disp=0
      drij=0.0d0
      dc6_rest=0.0d0
      dc6_rest_sum=0.0d0
      c6save=0.0d0
      kat=0
      dc6i=0.0d0
      dc6ij=0.0d0

      IF ( mykey == 0 ) THEN

      do iat=na_s, na_e
        call get_dC6_dCNij(maxc,max_elem,c6ab,mxc(iz(iat)),&
            & mxc(iz(iat)),cn(iat),cn(iat),iz(iat),iz(iat),iat,iat,&
            & c6,dc6iji,dc6ijj)

        dc6ij(iat,iat)=dc6iji
        c6save(lin(iat,iat))=c6
        r42=r2r4(iz(iat))*r2r4(iz(iat))
        rcovij=rcov(iz(iat))+rcov(iz(iat))

        R0=a1*sqrt(3.0d0*r42)+a2

        do taux=-rep_v(1),rep_v(1)
          do tauy=-rep_v(2),rep_v(2)
            do tauz=-rep_v(3),rep_v(3)
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

              !first dE/d(tau) saved in drij(i,i,counter)
              rij=tau
              r2=sum(rij*rij)
              ! if (r2.gt.rthr) cycle

              ! if (r2.gt.0.1) then
              if (r2.gt.0.1.and.r2.lt.rthr) then
                !
                ! get_dC6_dCNij calculates the derivative dC6(iat,jat)/dCN(iat) and
                ! dC6(iat,jat)/dCN(jat). these are saved in dC6ij for the kat loop
                !
                r=dsqrt(r2)
                r4=r2*r2
                r6=r4*r2
                r7=r6*r
                r8=r6*r2
                r9=r8*r

                !
                ! Calculates damping functions:

                t6=(r6+R0**6)
                t8=(r8+R0**8)

                drij(tauz,tauy,taux,lin(iat,iat))=drij(tauz,tauy,taux,lin(iat, &
                    & iat))&
                    & -s6*C6*6.0d0*r4*r/(t6*t6)*0.5d0&
                    & -s8*C6*24.0d0*r42*r7/(t8*t8)*0.5d0


                !
                ! in dC6_rest all terms BUT C6-term is saved for the kat-loop
                !
                dc6_rest=&
                    & (s6/t6+3.d0*s8*r42/t8)*0.50d0


                disp=disp-dc6_rest*c6

                dc6i(iat)=dc6i(iat)+dc6_rest*(dc6iji+dc6ijj)
                ! if (r2.lt.crit_cn)
                dc6_rest_sum(lin(iat,iat))=dc6_rest_sum(lin(iat,iat))+&
                    & dc6_rest


              else
                drij(tauz,tauy,taux,lin(iat,iat))=0.0d0
              end if


            end do
          end do
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! B E G I N jat L O O P
!!!!!!!!!!!!!!!!!!!!!!!!!!
        do jat=1,iat-1
          !
          ! get_dC6_dCNij calculates the derivative dC6(iat,jat)/dCN(iat) and
          ! dC6(iat,jat)/dCN(jat). these are saved in dC6ij for the kat loop
          !
          call get_dC6_dCNij(maxc,max_elem,c6ab,mxc(iz(iat)),&
              & mxc(iz(jat)),cn(iat),cn(jat),iz(iat),iz(jat),iat,jat,&
              & c6,dc6iji,dc6ijj)

          r42=r2r4(iz(iat))*r2r4(iz(jat))
          rcovij=rcov(iz(iat))+rcov(iz(jat))

          R0=a1*dsqrt(3.0d0*r42)+a2

          linij=lin(iat,jat)
          dc6ij(iat,jat)=dc6iji
          dc6ij(jat,iat)=dc6ijj
          c6save(linij)=c6
          do taux=-rep_v(1),rep_v(1)
            do tauy=-rep_v(2),rep_v(2)
              do tauz=-rep_v(3),rep_v(3)
                tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)


                rij=xyz(:,jat)-xyz(:,iat)+tau
                r2=sum(rij*rij)
                if (r2.gt.rthr) cycle


                r=dsqrt(r2)
                r4=r2*r2
                r6=r4*r2
                r7=r6*r
                r8=r6*r2
                r9=r8*r

                !
                ! Calculates damping functions:
                t6=(r6+R0**6)
                t8=(r8+R0**8)


                drij(tauz,tauy,taux,linij)=drij(tauz,tauy,taux,&
                    & linij)&
                    & -s6*C6*6.0d0*r4*r/(t6*t6)&
                    & -s8*C6*24.0d0*r42*r7/(t8*t8)

                !
                ! in dC6_rest all terms BUT C6-term is saved for the kat-loop
                !
                dc6_rest=&
                    & (s6/t6+3.d0*s8*r42/t8)


                disp=disp-dc6_rest*c6

                dc6i(iat)=dc6i(iat)+dc6_rest*dc6iji
                dc6i(jat)=dc6i(jat)+dc6_rest*dc6ijj
                ! if (r2.lt.crit_cn)
                dc6_rest_sum(lin(iat,jat))=dc6_rest_sum(linij)&
                    & +dc6_rest


              end do
            end do
          end do

        end do

      end do
      END IF

    end if

!!!!!!!!!!!!!!!!!!!!!!!
    !! BEGIN Threebody gradient
!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.noabc) then

      ! write(*,*)'!!!!!!!!!! THREEBODY GRADIENT !!!!!!!!!!'
      sr9=0.75d0
      eabc=0.0d0
      abcthr=crit_cn
      repv=rep_cn
      ! write(*,*)'thr:',sqrt(abcthr)

      call cpu_time(time1)

      rep_cn1 = rep_cn(1) 
      rep_cn2 = rep_cn(2) 
      rep_cn3 = rep_cn(3) 
      rep_v1 = rep_v(1)
      rep_v2 = rep_v(2)
      rep_v3 = rep_v(3)

      CALL mp_sum ( c6save , comm_dftd3 )
      CALL mp_sum ( dc6ij  , comm_dftd3 )
      IF ( mykey == 0 ) THEN
      na_smax = max(3,na_s)

!$acc data copyin(xyz(1:3,1:n),iz(1:n),lat(1:3,1:3),r0ab(1:max_elem,1:max_elem),c6save(1:n*(n+1)),dc6ij(1:n,1:n)) &
!$acc&            copy(dc6i(1:n),drij(-rep_v3:rep_v3,-rep_v2:rep_v2,-rep_v1:rep_v1,1:n*(n+1)/2)) 
!$acc parallel vector_length(32) 
!$acc loop collapse(3) gang  private(ijvec1,ijvec2,ijvec3, ikvec1,ikvec2,ikvec3, jkvec1,jkvec2,jkvec3, c6ij,c6ik,c6jk,c9, linij,linik,linjk, &
!$acc&                              jtau1,jtau2,jtau3, rij2,rr0ij, repmin1,repmin2,repmin3,repmax1,repmax2,repmax3 ) &
!$acc&                       reduction(+:eabc) 
      do iat=na_smax,na_e
        do jat=2, n
          do kat=1, n 
            if((jat.ge.iat).or.(kat.ge.jat)) cycle    
            linij=lin(iat,jat)
            ijvec1=xyz(1,jat)-xyz(1,iat)
            ijvec2=xyz(2,jat)-xyz(2,iat)
            ijvec3=xyz(3,jat)-xyz(3,iat)

            c6ij=c6save(linij)

            linik=lin(iat,kat)
            linjk=lin(jat,kat)
            ikvec1=xyz(1,kat)-xyz(1,iat)
            ikvec2=xyz(2,kat)-xyz(2,iat)
            ikvec3=xyz(3,kat)-xyz(3,iat)
            jkvec1=xyz(1,kat)-xyz(1,jat)
            jkvec2=xyz(2,kat)-xyz(2,jat)
            jkvec3=xyz(3,kat)-xyz(3,jat)

            c6ik=c6save(linik)
            c6jk=c6save(linjk)
            c9=-1.0d0*dsqrt(c6ij*c6ik*c6jk)
!$acc loop seq independent
            do jtaux=-rep_cn1,rep_cn1
              repmin1=max(-rep_cn1,jtaux-rep_cn1)
              repmax1=min(rep_cn1,jtaux+rep_cn1)
!$acc loop seq independent   
              do jtauy=-rep_cn2,rep_cn2
                repmin2=max(-rep_cn2,jtauy-rep_cn2)
                repmax2=min(rep_cn2,jtauy+rep_cn2)
!$acc loop seq independent 
                do jtauz=-rep_cn3,rep_cn3
                  repmin3=max(-rep_cn3,jtauz-rep_cn3)
                  repmax3=min(rep_cn3,jtauz+rep_cn3)
                  jtau1=jtaux*lat(1,1)+jtauy*lat(1,2)+jtauz*lat(1,3)
                  jtau2=jtaux*lat(2,1)+jtauy*lat(2,2)+jtauz*lat(2,3)
                  jtau3=jtaux*lat(3,1)+jtauy*lat(3,2)+jtauz*lat(3,3)
                  rij2= (ijvec1+jtau1)*(ijvec1+jtau1) + (ijvec2+jtau2)*(ijvec2+jtau2) + (ijvec3+jtau3)*(ijvec3+jtau3) 
                  if (rij2.gt.abcthr)cycle

                  rr0ij=DSQRT(rij2)/r0ab(iz(iat),iz(jat))

!$acc loop vector collapse(3) private(ktau1,ktau2,ktau3, dumvec1,dumvec2,dumvec3, rik2,rjk2,rr0ik,rr0jk, &
!$acc&                                geomean,geomean2,geomean3,r0av,damp9,ang,dc6_rest,dfdmp,r,dang,tmp1,dc9) &
!$acc&                        reduction(+:eabc) 
                  do ktaux=repmin1,repmax1
                    do ktauy=repmin2,repmax2
                      do ktauz=repmin3,repmax3
                        ktau1=ktaux*lat(1,1)+ktauy*lat(1,2)+ktauz*lat(1,3)
                        ktau2=ktaux*lat(2,1)+ktauy*lat(2,2)+ktauz*lat(2,3)
                        ktau3=ktaux*lat(3,1)+ktauy*lat(3,2)+ktauz*lat(3,3)
                        rik2=(ikvec1+ktau1)*(ikvec1+ktau1)+(ikvec2+ktau2)*(ikvec2+ktau2)+(ikvec3+ktau3)*(ikvec3+ktau3)
                        if (rik2.gt.abcthr)cycle

                        dumvec1=jkvec1+ktau1-jtau1
                        dumvec2=jkvec2+ktau2-jtau2
                        dumvec3=jkvec3+ktau3-jtau3
                        rjk2=dumvec1*dumvec1+dumvec2*dumvec2+dumvec3*dumvec3
                        if (rjk2.gt.abcthr)cycle
                        rr0ik=dsqrt(rik2)/r0ab(iz(iat),iz(kat))
                        rr0jk=dsqrt(rjk2)/r0ab(iz(jat),iz(kat))
                        geomean2=(rij2*rjk2*rik2)
                        ! first calculate the three components for the energy calculation fdmp
                        ! and ang
                        r0av=(rr0ij*rr0ik*rr0jk)**(1.0d0/3.0d0)
                        damp9=1./(1.+6.*(sr9*r0av)**alp9)

                        geomean=dsqrt(geomean2)
                        geomean3=geomean*geomean2
                        ang=0.375d0*(rij2+rjk2-rik2)*(rij2-rjk2+rik2)&
                            & *(-rij2+rjk2+rik2)/(geomean3*geomean2)&
                            & +1.0d0/(geomean3)

                        dc6_rest=ang*damp9
                        eabc=eabc+dc6_rest*c9
                        !
                        !start calculating the gradient components dfdmp, dang and dc9

                        !dfdmp is the same for all three distances
                        dfdmp=2.d0*alp9*(0.75d0*r0av)**(alp9)*damp9*damp9

                        !start calculating the derivatives of each part w.r.t. r_ij
                        r=dsqrt(rij2)


                        dang=-0.375d0*(rij2**3+rij2**2*(rjk2+rik2)&
                            & +rij2*(3.0d0*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2)&
                            & -5.0*(rjk2-rik2)**2*(rjk2+rik2))&
                            & /(r*geomean3*geomean2)

                        tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
!$acc atomic update
                        drij(jtauz,jtauy,jtaux,linij)= drij(jtauz,jtauy,jtaux,linij)-tmp1
!$acc end atomic 
                        !start calculating the derivatives of each part w.r.t. r_ik

                        r=dsqrt(rik2)


                        dang=-0.375d0*(rik2**3+rik2**2*(rjk2+rij2)&
                            & +rik2*(3.0d0*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2)&
                            & -5.0*(rjk2-rij2)**2*(rjk2+rij2))&
                            & /(r*geomean3*geomean2)

                        tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                        ! tmp1=-dc9
!$acc atomic update
                        drij(ktauz,ktauy,ktaux,linik)= drij(ktauz,ktauy,ktaux,linik)-tmp1
!$acc end atomic

                        !
                        !start calculating the derivatives of each part w.r.t. r_jk

                        r=dsqrt(rjk2)

                        dang=-0.375d0*(rjk2**3+rjk2**2*(rik2+rij2)&
                            & +rjk2*(3.0d0*rik2**2+2.0*rik2*rij2+3.0*rij2**2)&
                            & -5.0*(rik2-rij2)**2*(rik2+rij2))&
                            & /(r*geomean3*geomean2)

                        tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
!$acc atomic update
                        drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)= drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)-tmp1
!$acc end atomic

                        !calculating the CN derivative dE_disp(ijk)/dCN(i)

                        dc9=dc6ij(iat,jat)/c6ij+dc6ij(iat,kat)/c6ik
                        dc9=0.5d0*c9*dc9
!$acc atomic update 
                        dc6i(iat) = dc6i(iat) + dc6_rest*dc9   
!$acc end atomic

                        dc9=dc6ij(jat,iat)/c6ij+dc6ij(jat,kat)/c6jk
                        dc9=0.5d0*c9*dc9
!$acc atomic update 
                        dc6i(jat) = dc6i(jat) + dc6_rest*dc9  
!$acc end atomic

                        dc9=dc6ij(kat,iat)/c6ik+dc6ij(kat,jat)/c6jk
                        dc9=0.5d0*c9*dc9
!$acc atomic update 
                        dc6i(kat) = dc6i(kat) + dc6_rest*dc9   
!$acc end atomic
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
!$acc end parallel 
!$acc end data 

      ! Now the interaction with jat=iat of the triples iat,iat,kat
      do iat=max(2,na_s),na_e
        jat=iat
        linij=lin(iat,jat)
        ijvec=0.0d0

        c6ij=c6save(linij)
        do kat=1,iat-1
          linjk=lin(jat,kat)
          linik=linjk

          c6ik=c6save(linik)
          c6jk=c6ik
          ikvec=xyz(:,kat)-xyz(:,iat)
          jkvec=ikvec
          c9=-dsqrt(c6ij*c6ik*c6jk)
          do jtaux=-repv(1),repv(1)
            repmin(1)=max(-repv(1),jtaux-repv(1))
            repmax(1)=min(repv(1),jtaux+repv(1))
            do jtauy=-repv(2),repv(2)
              repmin(2)=max(-repv(2),jtauy-repv(2))
              repmax(2)=min(repv(2),jtauy+repv(2))
              do jtauz=-repv(3),repv(3)
                repmin(3)=max(-repv(3),jtauz-repv(3))
                repmax(3)=min(repv(3),jtauz+repv(3))
                if (jtaux.eq.0 .and. jtauy.eq.0 .and. jtauz.eq.0) cycle
                jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
                dumvec=jtau
                rij2=SUM(dumvec*dumvec)
                if (rij2.gt.abcthr)cycle

                rr0ij=DSQRT(rij2)/r0ab(iz(iat),iz(jat))

                do ktaux=repmin(1),repmax(1)
                  do ktauy=repmin(2),repmax(2)
                    do ktauz=repmin(3),repmax(3)
                      ! every result * 0.5

                      ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
                      dumvec=ikvec+ktau
                      dumvec=dumvec*dumvec
                      rik2=SUM(dumvec)
                      if (rik2.gt.abcthr)cycle

                      dumvec=jkvec+ktau-jtau
                      dumvec=dumvec*dumvec
                      rjk2=SUM(dumvec)
                      if (rjk2.gt.abcthr)cycle
                      rr0ik=DSQRT(rik2)/r0ab(iz(iat),iz(kat))
                      rr0jk=DSQRT(rjk2)/r0ab(iz(jat),iz(kat))


                      geomean2=(rij2*rjk2*rik2)
                      r0av=(rr0ij*rr0ik*rr0jk)**(1.0d0/3.0d0)
                      damp9=1./(1.+6.*(sr9*r0av)**alp9)

                      geomean=dsqrt(geomean2)
                      geomean3=geomean*geomean2
                      ang=0.375d0*(rij2+rjk2-rik2)*(rij2-rjk2+rik2)&
                          & *(-rij2+rjk2+rik2)/(geomean3*geomean2)&
                          & +1.0d0/(geomean3)


                      dc6_rest=ang*damp9/2.0d0
                      eabc=eabc+dc6_rest*c9

                      ! iat=jat
                      dfdmp=2.d0*alp9*(0.75d0*r0av)**(alp9)*damp9*damp9

                      !start calculating the derivatives of each part w.r.t. r_ij
                      r=dsqrt(rij2)

                      dang=-0.375d0*(rij2**3+rij2**2*(rjk2+rik2) &
                          & +rij2*(3.0d0*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2)&
                          & -5.0*(rjk2-rik2)**2*(rjk2+rik2))&
                          & /(r*geomean3*geomean2)

                      tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                      drij(jtauz,jtauy,jtaux,linij)=&
                          & drij(jtauz,jtauy,jtaux,linij)-tmp1/2.0

                      !start calculating the derivatives of each part w.r.t. r_ik
                      r=dsqrt(rik2)


                      dang=-0.375d0*(rik2**3+rik2**2*(rjk2+rij2)&
                          & +rik2*(3.0d0*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2)&
                          & -5.0*(rjk2-rij2)**2*(rjk2+rij2))&
                          & /(r*geomean3*geomean2)

                      tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                      drij(ktauz,ktauy,ktaux,linik)=&
                          & drij(ktauz,ktauy,ktaux,linik)-tmp1/2.0
                      !
                      !start calculating the derivatives of each part w.r.t. r_ik
                      r=dsqrt(rjk2)

                      dang=-0.375d0*(rjk2**3+rjk2**2*(rik2+rij2)&
                          & +rjk2*(3.0d0*rik2**2+2.0*rik2*rij2+3.0*rij2**2)&
                          & -5.0*(rik2-rij2)**2*(rik2+rij2))&
                          & /(r*geomean3*geomean2)

                      tmp1=-dang*c9*damp9+dfdmp/r*c9*ang

                      drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)=&
                          & drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)-tmp1/2.0

                      dc9=dc6ij(iat,jat)/c6ij+dc6ij(iat,kat)/c6ik
                      dc9=0.5d0*c9*dc9
                      dc6i(iat)=dc6i(iat)+dc6_rest*dc9

                      dc9=dc6ij(jat,iat)/c6ij+dc6ij(jat,kat)/c6jk
                      dc9=0.5d0*c9*dc9
                      dc6i(jat)=dc6i(jat)+dc6_rest*dc9

                      dc9=dc6ij(kat,iat)/c6ik+dc6ij(kat,jat)/c6jk
                      dc9=0.5d0*c9*dc9
                      dc6i(kat)=dc6i(kat)+dc6_rest*dc9




                    end do
                  end do
                end do

              end do
            end do
          end do
        end do
      end do

      do iat=max(2,na_s),na_e
        do jat=1,iat-1
          kat=jat
          linij=lin(iat,jat)
          linjk=lin(jat,kat)
          linik=linij

          c6ij=c6save(linij)
          c6ik=c6ij

          c6jk=c6save(linjk)
          ikvec=xyz(:,kat)-xyz(:,iat)
          ijvec=ikvec
          jkvec=0.0d0

          c9=-1.0d0*dsqrt(c6ij*c6ik*c6jk)
          do jtaux=-repv(1),repv(1)
            repmin(1)=max(-repv(1),jtaux-repv(1))
            repmax(1)=min(repv(1),jtaux+repv(1))
            do jtauy=-repv(2),repv(2)
              repmin(2)=max(-repv(2),jtauy-repv(2))
              repmax(2)=min(repv(2),jtauy+repv(2))
              do jtauz=-repv(3),repv(3)
                repmin(3)=max(-repv(3),jtauz-repv(3))
                repmax(3)=min(repv(3),jtauz+repv(3))

                jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
                dumvec=ijvec+jtau
                dumvec=dumvec*dumvec
                rij2=SUM(dumvec)
                if (rij2.gt.abcthr)cycle

                rr0ij=SQRT(rij2)/r0ab(iz(iat),iz(jat))

                do ktaux=repmin(1),repmax(1)
                  do ktauy=repmin(2),repmax(2)
                    do ktauz=repmin(3),repmax(3)
                      ! every result * 0.5
                      if (jtaux.eq.ktaux .and. jtauy.eq.ktauy&
                          & .and. jtauz.eq.ktauz) cycle
                      ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
                      dumvec=ikvec+ktau
                      dumvec=dumvec*dumvec
                      rik2=SUM(dumvec)
                      if (rik2.gt.abcthr)cycle
                      rr0ik=SQRT(rik2)/r0ab(iz(iat),iz(kat))

                      dumvec=jkvec+ktau-jtau
                      dumvec=dumvec*dumvec
                      rjk2=SUM(dumvec)
                      if (rjk2.gt.abcthr)cycle
                      rr0jk=SQRT(rjk2)/r0ab(iz(jat),iz(kat))

                      ! if (rij*rjk*rik.gt.abcthr)cycle

                      geomean2=(rij2*rjk2*rik2)
                      r0av=(rr0ij*rr0ik*rr0jk)**(1.0d0/3.0d0)
                      damp9=1./(1.+6.d0*(sr9*r0av)**alp9)

                      geomean=dsqrt(geomean2)
                      geomean3=geomean*geomean2
                      ang=0.375d0*(rij2+rjk2-rik2)*(rij2-rjk2+rik2)&
                          & *(-rij2+rjk2+rik2)/(geomean2*geomean3)&
                          & +1.0d0/(geomean3)
                      dc6_rest=ang*damp9/2.0d0
                      eabc=eabc+dc6_rest*c9


                      ! jat=kat
                      dfdmp=2.d0*alp9*(0.75d0*r0av)**(alp9)*damp9*damp9
                      !start calculating the derivatives of each part w.r.t. r_ij
                      r=dsqrt(rij2)

                      dang=-0.375d0*(rij2**3+rij2**2*(rjk2+rik2)&
                          & +rij2*(3.0d0*rjk2**2+2.0d0*rjk2*rik2+3.0d0*rik2**2)&
                          & -5.0d0*(rjk2-rik2)**2*(rjk2+rik2))&
                          & /(r*geomean3*geomean2)

                      tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                      drij(jtauz,jtauy,jtaux,linij)=&
                          & drij(jtauz,jtauy,jtaux,linij)-tmp1/2.0d0

                      !start calculating the derivatives of each part w.r.t. r_ik
                      r=dsqrt(rik2)


                      dang=-0.375d0*(rik2**3+rik2**2*(rjk2+rij2)&
                          & +rik2*(3.0d0*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2)&
                          & -5.0*(rjk2-rij2)**2*(rjk2+rij2))&
                          & /(r*geomean3*geomean2)

                      tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                      ! tmp1=-dc9
                      drij(ktauz,ktauy,ktaux,linik)=&
                          & drij(ktauz,ktauy,ktaux,linik)-tmp1/2.0d0
                      !
                      !start calculating the derivatives of each part w.r.t. r_jk
                      r=dsqrt(rjk2)

                      dang=-0.375d0*(rjk2**3+rjk2**2*(rik2+rij2)&
                          & +rjk2*(3.0d0*rik2**2+2.0*rik2*rij2+3.0*rij2**2)&
                          & -5.0d0*(rik2-rij2)**2*(rik2+rij2))&
                          & /(r*geomean3*geomean2)

                      tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                      drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)=&
                          & drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)-tmp1/2.0d0

                      !calculating the CN derivative dE_disp(ijk)/dCN(i)

                      dc9=dc6ij(iat,jat)/c6ij+dc6ij(iat,kat)/c6ik
                      dc9=0.5d0*c9*dc9
                      dc6i(iat)=dc6i(iat)+dc6_rest*dc9

                      dc9=dc6ij(jat,iat)/c6ij+dc6ij(jat,kat)/c6jk
                      dc9=0.5d0*c9*dc9
                      dc6i(jat)=dc6i(jat)+dc6_rest*dc9

                      dc9=dc6ij(kat,iat)/c6ik+dc6ij(kat,jat)/c6jk
                      dc9=0.5d0*c9*dc9
                      dc6i(kat)=dc6i(kat)+dc6_rest*dc9




                    end do
                  end do
                end do

              end do
            end do
          end do
        end do
      end do


      ! And finally the self interaction iat=jat=kat all

      idum=0
      do iat=na_s,na_e
        jat=iat
        kat=iat
        ijvec=0.0d0
        linij=lin(iat,jat)
        linik=lin(iat,kat)
        linjk=lin(jat,kat)
        ikvec=ijvec
        jkvec=ikvec
        c6ij=c6save(linij)
        c6ik=c6ij
        c6jk=c6ij
        c9=-(DSQRT(c6ij*c6ij*c6ij))

        do jtaux=-repv(1),repv(1)
          repmin(1)=max(-repv(1),jtaux-repv(1))
          repmax(1)=min(repv(1),jtaux+repv(1))
          do jtauy=-repv(2),repv(2)
            repmin(2)=max(-repv(2),jtauy-repv(2))
            repmax(2)=min(repv(2),jtauy+repv(2))
            do jtauz=-repv(3),repv(3)
              repmin(3)=max(-repv(3),jtauz-repv(3))
              repmax(3)=min(repv(3),jtauz+repv(3))
              if ((jtaux.eq.0) .and.(jtauy.eq.0) .and.(jtauz.eq.0))cycle
              jtau=jtaux*lat(:,1)+jtauy*lat(:,2)+jtauz*lat(:,3)
              dumvec=jtau
              dumvec=dumvec*dumvec
              rij2=SUM(dumvec)
              if (rij2.gt.abcthr)cycle
              rr0ij=SQRT(rij2)/r0ab(iz(iat),iz(jat))

              do ktaux=repmin(1),repmax(1)
                do ktauy=repmin(2),repmax(2)
                  do ktauz=repmin(3),repmax(3)
                    if ((ktaux.eq.0) .and.( ktauy.eq.0) .and.( ktauz.eq.0))cycle
                    if ((ktaux.eq.jtaux) .and. (ktauy.eq.jtauy)&
                        & .and. (ktauz.eq.jtauz)) cycle

                    ! every result * 1/6 becaues every triple is counted twice due to the tw
                    !
                    !plus 1/3 becaues every triple is three times in each unitcell
                    ktau=ktaux*lat(:,1)+ktauy*lat(:,2)+ktauz*lat(:,3)
                    dumvec=ktau
                    dumvec=dumvec*dumvec
                    rik2=SUM(dumvec)
                    if (rik2.gt.abcthr)cycle
                    rr0ik=SQRT(rik2)/r0ab(iz(iat),iz(kat))

                    dumvec=jkvec+ktau-jtau
                    dumvec=dumvec*dumvec
                    rjk2=SUM(dumvec)
                    if (rjk2.gt.abcthr)cycle
                    rr0jk=SQRT(rjk2)/r0ab(iz(jat),iz(kat))

                    geomean2=(rij2*rjk2*rik2)
                    r0av=(rr0ij*rr0ik*rr0jk)**(1.0d0/3.0d0)
                    damp9=1./(1.+6.*(sr9*r0av)**alp9)

                    geomean=dsqrt(geomean2)
                    geomean3=geomean*geomean2
                    ang=0.375d0*(rij2+rjk2-rik2)*(rij2-rjk2+rik2)&
                        & *(-rij2+rjk2+rik2)/(geomean2*geomean3)&
                        & +1.0d0/(geomean3)
                    dc6_rest=ang*damp9/6.0d0
                    eabc=eabc+c9*dc6_rest

                    ! iat=jat=kat
                    dfdmp=2.d0*alp9*(0.75d0*r0av)**(alp9)*damp9*damp9
                    !start calculating the derivatives of each part w.r.t. r_ij

                    r=dsqrt(rij2)
                    dang=-0.375d0*(rij2**3+rij2**2*(rjk2+rik2)&
                        & +rij2*(3.0d0*rjk2**2+2.0*rjk2*rik2+3.0*rik2**2)&
                        & -5.0*(rjk2-rik2)**2*(rjk2+rik2))&
                        & /(r*geomean3*geomean2)


                    tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                    drij(jtauz,jtauy,jtaux,linij)=&
                        & drij(jtauz,jtauy,jtaux,linij)-tmp1/6.0d0

                    !start calculating the derivatives of each part w.r.t. r_ik

                    r=dsqrt(rik2)

                    dang=-0.375d0*(rik2**3+rik2**2*(rjk2+rij2)&
                        & +rik2*(3.0d0*rjk2**2+2.0*rjk2*rij2+3.0*rij2**2)&
                        & -5.0*(rjk2-rij2)**2*(rjk2+rij2))&
                        & /(r*geomean3*geomean2)

                    tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                    drij(ktauz,ktauy,ktaux,linik)=&
                        & drij(ktauz,ktauy,ktaux,linik)-tmp1/6.0d0
                    !
                    !start calculating the derivatives of each part w.r.t. r_jk

                    r=dsqrt(rjk2)
                    dang=-0.375d0*(rjk2**3+rjk2**2*(rik2+rij2)&
                        & +rjk2*(3.0d0*rik2**2+2.0*rik2*rij2+3.0*rij2**2)&
                        & -5.0*(rik2-rij2)**2*(rik2+rij2))&
                        & /(r*geomean3*geomean2)

                    tmp1=-dang*c9*damp9+dfdmp/r*c9*ang
                    drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)=&
                        & drij(ktauz-jtauz,ktauy-jtauy,ktaux-jtaux,linjk)-tmp1/6.0d0


                    !calculating the CN derivative dE_disp(ijk)/dCN(i)

                    dc9=dc6ij(iat,jat)/c6ij+dc6ij(iat,kat)/c6ik
                    dc9=0.5d0*c9*dc9
                    dc6i(iat)=dc6i(iat)+dc6_rest*dc9

                    dc9=dc6ij(jat,iat)/c6ij+dc6ij(jat,kat)/c6jk
                    dc9=0.5d0*c9*dc9
                    dc6i(jat)=dc6i(jat)+dc6_rest*dc9

                    dc9=dc6ij(kat,iat)/c6ik+dc6ij(kat,jat)/c6jk
                    dc9=0.5d0*c9*dc9
                    dc6i(kat)=dc6i(kat)+dc6_rest*dc9





                  end do
                end do
              end do
            end do
          end do
          !jtaux
        end do

      end do

      END IF
      CALL mp_sum ( eabc , comm_dftd3 )

      call cpu_time(time2)

      ! write(*,*)' eabc(gdisp): ',eabc
      ! write(*,'('' time(abc) '',f6.1)')time2-time1
      disp=disp-eabc
      ! write(*,*)'gdisp:',disp
    end if

    CALL mp_sum ( drij , comm_dftd3 )
    CALL mp_sum ( dc6i , comm_dftd3 )

    sigma_abc=0.0d0
    sigma=0.0d0

    ! After calculating all derivatives dE/dr_ij w.r.t. distances,
    ! the grad w.r.t. the coordinates is calculated dE/dr_ij * dr_ij/dxyz_i
    do iat=2,n
      do jat=1,iat-1
        linij=lin(iat,jat)
        rcovij=rcov(iz(iat))+rcov(iz(jat))
        do taux=-rep_v(1),rep_v(1)
          do tauy=-rep_v(2),rep_v(2)
            do tauz=-rep_v(3),rep_v(3)
              tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)

              rij=xyz(:,jat)-xyz(:,iat)+tau
              r2=sum(rij*rij)
              if (r2.gt.rthr.or.r2.lt.0.5) cycle
              r=dsqrt(r2)

              if (r2.lt.crit_cn) then
                expterm=exp(-k1*(rcovij/r-1.d0))
                dcnn=-k1*rcovij*expterm/&
                    & (r2*(expterm+1.d0)*(expterm+1.d0))
              else
                dcnn=0.0d0
              end if
              x1=drij(tauz,tauy,taux,linij)+dcnn*(dc6i(iat)+dc6i(jat))
              vec=x1*rij/r
              g(:,iat)=g(:,iat)+vec
              g(:,jat)=g(:,jat)-vec
              do i=1,3
                do j=1,3
                  sigma(j,i)=sigma(j,i)+vec(j)*rij(i)
                end do
              end do



            end do
          end do
        end do
      end do
    end do

    do iat=1,n
      rcovij=rcov(iz(iat))+rcov(iz(iat))
      do taux=-rep_v(1),rep_v(1)
        do tauy=-rep_v(2),rep_v(2)
          do tauz=-rep_v(3),rep_v(3)
            if (taux.eq.0.and.tauy.eq.0.and.tauz.eq.0) cycle

            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
            r2=(sum(tau*tau))
            r=dsqrt(r2)
            if (r2.lt.crit_cn) then
              expterm=exp(-k1*(rcovij/r-1.d0))
              dcnn=-k1*rcovij*expterm/&
                  & (r2*(expterm+1.d0)*(expterm+1.d0))
            else
              dcnn=0.0d0
            end if
            x1=drij(tauz,tauy,taux,lin(iat,iat))+dcnn*dc6i(iat)
            vec=x1*tau/r
            vec2(1)=taux
            vec2(2)=tauy
            vec2(3)=tauz
            do i=1,3
              do j=1,3
                sigma(j,i)=sigma(j,i)+vec(j)*tau(i)
              end do
            end do


          end do
        end do
      end do



    end do



    stress=0.0d0
    glatabc=0.0d0
    call inv_cell(lat,lat_1)
    do a=1,3
      do b=1,3
        do my=1,3
          stress(a,b)=stress(a,b)-sigma(a,my)*lat_1(b,my)
        end do
      end do
    end do



    ! write(*,*)'drij:',drij(lin(iat,jat),:)
    ! write(*,*)'g:',g(1,1:3)
    ! write(*,*)'dcn:',sum(dcn(lin(2,1),:))



    deallocate(drij)




999 continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !This is where the D2 gradient and the numerical gradient jump.
    !
!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do i=1,n
    ! write(*,'(83F17.12)') g(1:3,i)
    ! end do
    gnorm=sum(abs(g(1:3,1:n)))
    if (echo)then
      ! write(*,*)'testsum:',testsum*autoev/autoang
      write(*,*)'|G(force)| =',gnorm
      gnorm=sum(abs(stress(1:3,1:3)))
      write(*,*)'|G(stress)|=',gnorm
    end if

  end subroutine pbcgdisp

  subroutine pbcgdisp_new(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,&
      & rcov,s6,s18,rs6,rs8,rs10,alp6,alp8,alp10,noabc,num,&
      & version,g,disp,gnorm,lat,rep_v,rep_cn,&
      & crit_vdw,echo,crit_cn, hstep, ia, ix, is, g_supercell_)
    !
    USE mp,           ONLY : mp_sum
    !
    ! input/output variables
    !   
    real(wp), intent(in) :: c6ab(max_elem,max_elem,maxc,maxc,3)
    real(wp), intent(in) :: s6,s18,rcov(max_elem)
    real(wp), intent(in) :: rs6,rs8,rs10,alp10,alp8,alp6
    integer,  intent(in) :: n,iz(*),max_elem,maxc,version,mxc(max_elem)
    real(wp), intent(in) :: xyz(3,*),r0ab(max_elem,max_elem),r2r4(*)
    logical,  intent(in) :: noabc,num,echo
    real(wp), dimension(3,3), intent(in) :: lat
    integer,  dimension(3),   intent(in) :: rep_v,rep_cn
    real(wp), intent(in) :: crit_vdw,crit_cn
    integer,  intent(in) :: ia, ix, is
    real(wp), intent(in) :: hstep ! step (in Bohr) for atom displacement for numerical Hessian calculation
    real(wp), intent(inout) :: g(3,*), disp
    real(wp), intent(out) :: gnorm
    real(wp), target, contiguous, intent(out) :: g_supercell_(:,:,:,:,:) 
    ! 
    ! local variables
    !
    real(wp), pointer :: g_supercell(:,:,:,:,:) 
    integer, allocatable :: ns(:)
    real(wp), allocatable :: xyz_hstep(:,:) ! displaced geometry
    real(wp) :: iat_jat_fact ! weight for diagonal/off-diagonal gradient terms
    logical  :: unit_cell     ! a quick flag to check whether a certain taux,tauy,tauz points to the unit cell or to some image
    integer  :: mykey, na_s, na_e
    real(wp) :: hdisp,a1,a2,R0,C6,R42,r2,s10,s8,cn(n),rthr,rcovij,dcnn
    integer  :: iat,jat, taux,tauy,tauz, linij
    real(wp), DIMENSION(3) :: tau, rij, vec
    real(wp), allocatable,dimension(:,:,:,:) :: drij, drij_hstep(:,:,:,:)
    real(wp), allocatable,dimension(:,:,:,:) :: dcn
    real(wp) :: dc6i(n), dc6i_hstep(n), dc6ij(n,n), dc6iji,dc6ijj
    real(wp) :: dc6_rest_sum(n*(n+1)/2)
    real(wp),DIMENSION(n*(n+1)) ::c6save
    real(wp) :: res1, res2, gnorm_supercell

    if(num) Call errore('pbcgdisp', 'Atom displacement not implemented with numerical forces', 1)
    if(.not.noabc) Call errore('pbcgdisp', 'Atom displacement not implemented with the threebody term ' // &
                                 ' (set dftd3_threebody=.false. for phonon calculations)', 1)

    ns = shape(g_supercell_)
    g_supercell( -ns(1)/2:ns(1)/2, -ns(2)/2:ns(2)/2, -ns(3)/2:ns(3)/2, 1:ns(4), 1:ns(5) ) => g_supercell_
    g_supercell(:,:,:,:,:) = 0.0_wp 
    hdisp = dble(is) * hstep
    allocate(xyz_hstep(3,n))
    xyz_hstep(1:3,1:n) = xyz(1:3,1:n)
    xyz_hstep(ix, ia) = xyz_hstep(ix, ia) + hdisp

    ! R^2 cut-off
    rthr=crit_vdw

    if (echo)write(*,*)
    if (echo) write(*,*) 'doing analytical gradient for version...', version

    CALL block_distribute( n, me_dftd3, nproc_dftd3, na_s, na_e, mykey )

    if ((version.eq.3).or.(version.eq.5).or.(version.eq.4).or.(version.eq.6)) then

      ! precompute for analytical part
      call pbcncoord_new(n,rcov,iz,xyz,cn,lat,rep_cn,crit_cn,xyz_hstep)


      a1 =rs6
      a2 =rs8
      s8 =s18
      s10=s18
      allocate(drij(-rep_v(3):rep_v(3),-rep_v(2):rep_v(2),-rep_v(1):rep_v(1),n*(n+1)/2), &
               drij_hstep(-rep_v(3):rep_v(3),-rep_v(2):rep_v(2),-rep_v(1):rep_v(1),n*(n+1)/2) )

      disp=0
      drij=0.0d0
      drij_hstep=0.0d0
      dc6_rest_sum=0.0d0
      c6save=0.0d0
      dc6i=0.0d0
      dc6i_hstep=0.0d0
      dc6ij=0.0d0

      IF ( mykey == 0 ) THEN

      do taux=-rep_v(1),rep_v(1)
        do tauy=-rep_v(2),rep_v(2)
          do tauz=-rep_v(3),rep_v(3)

            unit_cell = taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0  
            tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
            r2=sqrt(sum(tau*tau))
            if(unit_cell .and. (r2.gt.0.000010d0)) Call errore('pbcgdisp', &
                                                      'non zero traslation vector found for the unit cell', 1)

            do iat=na_s, na_e
              call get_dC6_dCNij(maxc,max_elem,c6ab,mxc(iz(iat)),&
                  & mxc(iz(iat)),cn(iat),cn(iat),iz(iat),iz(iat),iat,iat,&
                  & c6,dc6iji,dc6ijj)

              c6save(lin(iat,iat))=c6
              dc6ij(iat,iat)=dc6iji
              r42=r2r4(iz(iat))*r2r4(iz(iat))
              rcovij=rcov(iz(iat))+rcov(iz(iat))
              if((version.eq.3).or.(version.eq.5)) then 
                R0=r0ab(iz(iat),iz(iat))
              elseif((version.eq.4).or.(version.eq.6)) then
                R0=a1*sqrt(3.0d0*r42)+a2
              end if 

              iat_jat_fact = 0.50d0 ! each diagonal contribution is weighted 1/2

              rij=tau
              r2=sum(rij*rij)
              if (r2.gt.0.1.and.r2.lt.rthr) then
                !
                Call gkernel1 (version, r2, R0, s6, rs6, alp6, s8, rs8, alp8, C6, r42, iat_jat_fact, res1, res2)
                drij(tauz,tauy,taux,lin(iat,iat))=drij(tauz,tauy,taux,lin(iat,iat)) + res1
                disp=disp-res2*c6
                dc6i(iat)=dc6i(iat)+res2*(dc6iji+dc6ijj)
                dc6_rest_sum(lin(iat,iat))=dc6_rest_sum(lin(iat,iat))+res2
                !
                if(.not. unit_cell .and. (iat.eq.ia) ) then 
                  !
                  disp=disp+res2*c6 ! remove the previous contribution 
                  dc6_rest_sum(lin(iat,iat))=dc6_rest_sum(lin(iat,iat))-res2 ! remove the previous contribution
                  !
                  rij = xyz(:,iat) - xyz_hstep(:,iat) + tau
                  r2=sum(rij*rij)
                  Call gkernel1 (version, r2, R0, s6, rs6, alp6, s8, rs8, alp8, C6, r42, iat_jat_fact, res1, res2)
                  drij_hstep(tauz,tauy,taux,lin(iat,iat))=drij_hstep(tauz,tauy,taux,lin(iat,iat)) + res1
                  disp=disp-res2*c6
                  dc6i_hstep(iat)=dc6i_hstep(iat)+res2*(dc6iji+dc6ijj)
                  dc6_rest_sum(lin(iat,iat))=dc6_rest_sum(lin(iat,iat))+res2
                end if 
              else
                drij(tauz,tauy,taux,lin(iat,iat))=0.0d0
                if(.not.unit_cell) drij_hstep(tauz,tauy,taux,lin(iat,iat))=0.0d0
              end if

              do jat=1,iat-1
                !
                ! get_dC6_dCNij calculates the derivative dC6(iat,jat)/dCN(iat) and
                ! dC6(iat,jat)/dCN(jat)
                !
                call get_dC6_dCNij(maxc,max_elem,c6ab,mxc(iz(iat)),&
                    & mxc(iz(jat)),cn(iat),cn(jat),iz(iat),iz(jat),iat,jat,&
                    & c6,dc6iji,dc6ijj)
      
                r42=r2r4(iz(iat))*r2r4(iz(jat))
                rcovij=rcov(iz(iat))+rcov(iz(jat))
                linij=lin(iat,jat)
                dc6ij(iat,jat)=dc6iji
                dc6ij(jat,iat)=dc6ijj
                c6save(linij)=c6
                if((version.eq.3).or.(version.eq.5)) then 
                  R0=r0ab(iz(iat),iz(iat))
                elseif((version.eq.4).or.(version.eq.6)) then
                  R0=a1*dsqrt(3.0d0*r42)+a2 
                end if 

                iat_jat_fact = 1.0d0 ! off diagonal contributions for iat and jat are equivalent: they are weighted twice the diagonal ones 

                rij=xyz(:,jat)-xyz(:,iat)+tau

                if (r2.gt.rthr) cycle

                r2=sum(rij*rij)
                Call gkernel1 (version, r2, R0, s6, rs6, alp6, s8, rs8, alp8, C6, r42, iat_jat_fact, res1, res2)
                drij(tauz,tauy,taux,linij)=drij(tauz,tauy,taux,linij) + res1
                disp = disp - res2 * C6
                dc6i(iat)=dc6i(iat)+res2 * dc6iji
                dc6i(jat)=dc6i(jat)+res2 * dc6ijj
                dc6_rest_sum(linij) = dc6_rest_sum(linij) + res2

                if( iat.eq.ia ) then 
                  !
                  ! remove previous contributions from disp and dc6_rest_sum
                  disp = disp + res2 * C6 
                  dc6_rest_sum(linij) = dc6_rest_sum(linij) - res2
                  !
                  rij = ( xyz(:,jat) + tau ) - xyz_hstep(:,iat) 
                  r2=sum(rij*rij)
                  Call gkernel1 (version, r2, R0, s6, rs6, alp6, s8, rs8, alp8, C6, r42, iat_jat_fact, res1, res2)
                  drij_hstep(tauz,tauy,taux,linij)=drij_hstep(tauz,tauy,taux,linij) + res1 
                  disp = disp - res2 * C6
                  dc6i_hstep(iat)=dc6i_hstep(iat)+res2*dc6iji
                  dc6i_hstep(jat)=dc6i_hstep(jat)+res2*dc6ijj
                  dc6_rest_sum(linij) = dc6_rest_sum(linij) + res2
                  !
                elseif( jat.eq.ia ) then 
                  !
                  ! remove previous contributions from disp and dc6_rest_sum
                  disp = disp + res2 * C6 
                  dc6_rest_sum(linij) = dc6_rest_sum(linij) - res2
                  !
                  rij = xyz_hstep(:,jat) - ( xyz(:,iat) - tau )  
                  r2=sum(rij*rij)
                  Call gkernel1 (version, r2, R0, s6, rs6, alp6, s8, rs8, alp8, C6, r42, iat_jat_fact, res1, res2)
                  drij_hstep(tauz,tauy,taux,linij)=drij_hstep(tauz,tauy,taux,linij) + res1 
                  disp = disp - res2 * C6
                  dc6i_hstep(iat)=dc6i_hstep(iat)+res2*dc6iji
                  dc6i_hstep(jat)=dc6i_hstep(jat)+res2*dc6ijj
                  dc6_rest_sum(linij) = dc6_rest_sum(linij) + res2
                end if 
                
              end do ! jat 
            end do ! iat

          end do ! tauz 
        end do ! tauy
      end do ! taux

      END IF ! mykey == 0

      CALL mp_sum ( drij , comm_dftd3 )
      CALL mp_sum ( dc6i , comm_dftd3 )
      CALL mp_sum ( drij_hstep , comm_dftd3 )
      CALL mp_sum ( dc6i_hstep , comm_dftd3 )

    end if ! version

    ! After calculating all derivatives dE/dr_ij w.r.t. distances,
    ! the grad w.r.t. the coordinates is calculated dE/dr_ij * dr_ij/dxyz_i

    do taux=-rep_v(1),rep_v(1)
      do tauy=-rep_v(2),rep_v(2)
        do tauz=-rep_v(3),rep_v(3)
          !
          unit_cell = taux.eq.0 .and. tauy.eq.0 .and. tauz.eq.0  
          tau=taux*lat(:,1)+tauy*lat(:,2)+tauz*lat(:,3)
          r2=sqrt(sum(tau*tau))
          if(unit_cell .and. (r2.gt.0.000010d0)) Call errore('pbcgdisp', &
                                                      'non zero traslation vector found for the unit cell', 1)
          !
          do iat = 1, n 
            !
            if(.not. unit_cell .and. (iat.eq.ia) ) then 
              !
              linij=lin(iat,iat)
              rcovij=rcov(iz(iat))+rcov(iz(iat))
              !
              rij = xyz(:,iat) - xyz_hstep(:,iat) + tau
              r2=sum(rij*rij)
              if(version.eq.2) then 
                iat_jat_fact = 0.5_wp
                R0=r0ab(iz(iat),iz(iat))*rs6
                c6=c6ab(iz(iat),iz(iat),1,1,1)*s6
                Call gkernel3 ( iat_jat_fact, rij, r2, alp6, R0, c6, vec  )
              else
                Call gkernel2(rij, r2, crit_cn, rcovij, drij_hstep(tauz,tauy,taux,linij), dc6i(iat), dc6i_hstep(iat), vec)
              endif 
              g(:,iat)=g(:,iat)+vec
              g_supercell(0,0,0,1:3,iat)          = g_supercell(0,0,0,1:3,iat)          + vec
              g_supercell(tauz,tauy,taux,1:3,iat) = g_supercell(tauz,tauy,taux,1:3,iat) - vec
              !
            !else 
            !  whenever iat is not displaced (i.e. iat.ne.ia), forces from -tau and +tau cancel out each other 
            !  whenever iat is in the unit_cell, there is no force contribution from iat with itself
            end if 
            !
            do jat=1,iat-1
              !
              linij=lin(iat,jat)
              rcovij=rcov(iz(iat))+rcov(iz(jat))
              if(version.eq.2) then 
                iat_jat_fact = 1.0_wp
                R0=r0ab(iz(jat),iz(iat))*rs6
                c6=c6ab(iz(jat),iz(iat),1,1,1)*s6
              end if 
              !
              if(unit_cell) then 
                !
                if(iat.eq.ia) then
                  rij=xyz(:,jat)-xyz_hstep(:,iat) 
                  r2=sum(rij*rij)
                  if (r2.gt.rthr.or.r2.lt.0.5) cycle
                  if(version.eq.2) then 
                    Call gkernel3 ( iat_jat_fact, rij, r2, alp6, R0, c6, vec  )
                  else
                    Call gkernel2(rij, r2, crit_cn, rcovij, drij_hstep(tauz,tauy,taux,linij), dc6i_hstep(iat), dc6i(jat), vec)
                  endif 
                elseif(jat.eq.ia) then
                  rij=xyz_hstep(:,jat)-xyz(:,iat) 
                  r2=sum(rij*rij)
                  if (r2.gt.rthr.or.r2.lt.0.5) cycle
                  if(version.eq.2) then 
                    Call gkernel3 ( iat_jat_fact, rij, r2, alp6, R0, c6, vec  )
                  else
                    Call gkernel2(rij, r2, crit_cn, rcovij, drij_hstep(tauz,tauy,taux,linij), dc6i(iat), dc6i_hstep(jat), vec)
                  end if
                else
                  rij=xyz(:,jat)-xyz(:,iat) 
                  r2=sum(rij*rij)
                  if (r2.gt.rthr.or.r2.lt.0.5) cycle
                  if(version.eq.2) then 
                    Call gkernel3 ( iat_jat_fact, rij, r2, alp6, R0, c6, vec )
                  else
                    Call gkernel2(rij, r2, crit_cn, rcovij, drij(tauz,tauy,taux,linij), dc6i(iat), dc6i(jat), vec)
                  end if 
                end if 
                g(:,iat)=g(:,iat)+vec
                g(:,jat)=g(:,jat)-vec
                g_supercell(tauz,tauy,taux,1:3,iat) = g_supercell(tauz,tauy,taux,1:3,iat) + vec
                g_supercell(tauz,tauy,taux,1:3,jat) = g_supercell(tauz,tauy,taux,1:3,jat) - vec
                !
              elseif(.not. unit_cell) then 
                !
                if(iat.eq.ia) then 
                  !
                  rij = ( xyz(:,jat) + tau ) - xyz_hstep(:,iat) ! displaced 
                  r2=sum(rij*rij)                
                  if (r2.gt.rthr.or.r2.lt.0.5) cycle
                  if(version.eq.2) then 
                    Call gkernel3 ( iat_jat_fact, rij, r2, alp6, R0, c6, vec  )
                  else
                    Call gkernel2(rij, r2, crit_cn, rcovij, drij_hstep(tauz,tauy,taux,linij), dc6i_hstep(iat), dc6i(jat), vec)
                  end if 
                  g_supercell(tauz,tauy,taux,1:3,jat) = g_supercell(tauz,tauy,taux,1:3,jat) - vec       
                  g_supercell(0   ,0   ,0   ,1:3,iat) = g_supercell(0   ,0   ,0   ,1:3,iat) + vec
                  g(:,iat)=g(:,iat)+vec
                  !
                  rij = ( xyz(:,jat) + tau ) - xyz(:,iat)        ! undisplaced
                  r2=sum(rij*rij)                
                  if(version.eq.2) then 
                    Call gkernel3 ( iat_jat_fact, rij, r2, alp6, R0, c6, vec  )
                  else
                    Call gkernel2(rij, r2, crit_cn, rcovij, drij(tauz,tauy,taux,linij), dc6i(iat), dc6i(jat), vec)
                  end if  
                  g_supercell(0   ,0   ,0   ,1:3,jat) = g_supercell(0   ,0   ,0   ,1:3,jat) - vec       
                  g_supercell(-tauz,-tauy,-taux,1:3,iat) = g_supercell(-tauz,-tauy,-taux,1:3,iat) + vec
                  g(:,jat)=g(:,jat)-vec
                  !
                elseif(jat.eq.ia) then 
                  !
                  rij = xyz_hstep(:,jat) - ( xyz(:,iat) - tau )  ! displaced
                  r2=sum(rij*rij)                
                  if (r2.gt.rthr.or.r2.lt.0.5) cycle
                  if(version.eq.2) then 
                    Call gkernel3 ( iat_jat_fact, rij, r2, alp6, R0, c6, vec )
                  else
                    Call gkernel2(rij, r2, crit_cn, rcovij, drij_hstep(tauz,tauy,taux,linij), dc6i(iat), dc6i_hstep(jat), vec)
                  endif 
                  !vec = -vec
                  g_supercell(0   ,0   ,0   ,1:3,jat) = g_supercell(0   ,0   ,0   ,1:3,jat) - vec       
                  g_supercell(-tauz,-tauy,-taux,1:3,iat) = g_supercell(-tauz,-tauy,-taux,1:3,iat) + vec
                  g(:,jat)=g(:,jat)-vec
                  !
                  rij = xyz(:,jat)  - ( xyz(:,iat) - tau )        ! undisplaced
                  r2=sum(rij*rij)                
                  if(version.eq.2) then 
                    Call gkernel3 ( iat_jat_fact, rij, r2, alp6, R0, c6, vec )
                  else
                    Call gkernel2(rij, r2, crit_cn, rcovij, drij(tauz,tauy,taux,linij), dc6i(iat), dc6i(jat), vec)
                  endif 
                  g_supercell(tauz,tauy,taux,1:3,jat) = g_supercell(tauz,tauy,taux,1:3,jat) - vec       
                  g_supercell(0   ,0   ,0   ,1:3,iat) = g_supercell(0   ,0   ,0   ,1:3,iat) + vec
                  g(:,iat)=g(:,iat)+vec
                  !
                else
                  !
                  rij = ( xyz(:,jat) + tau ) - xyz(:,iat)        ! undisplaced
                  r2=sum(rij*rij)                
                  if (r2.gt.rthr.or.r2.lt.0.5) cycle
                  if(version.eq.2) then 
                    Call gkernel3 ( iat_jat_fact, rij, r2, alp6, R0, c6, vec )
                  else
                    Call gkernel2(rij, r2, crit_cn, rcovij, drij(tauz,tauy,taux,linij), dc6i(iat), dc6i(jat), vec)
                  endif 
                  g_supercell(tauz,tauy,taux,1:3,jat) = g_supercell(tauz,tauy,taux,1:3,jat) - vec       
                  g_supercell(-tauz,-tauy,-taux,1:3,iat) = g_supercell(-tauz,-tauy,-taux,1:3,iat) + vec
                  g_supercell(0   ,0   ,0   ,1:3,jat) = g_supercell(0   ,0   ,0   ,1:3,jat) - vec
                  g_supercell(0   ,0   ,0   ,1:3,iat) = g_supercell(0   ,0   ,0   ,1:3,iat) + vec
                  g(:,jat)=g(:,jat)-vec
                  g(:,iat)=g(:,iat)+vec
                  !
                end if 
                !
              end if ! unit_cell 
              !
            end do
          end do
        end do
      end do
    end do
    !
    if(allocated(drij)) deallocate(drij)
    if(allocated(drij_hstep)) deallocate(drij_hstep)
    !
    gnorm_supercell=sum(abs(g_supercell(:,:,:,1:3,:)))
    write(*,*)'|G(force)| =',gnorm_supercell, ' supercell gradient'
    gnorm_supercell=sum(abs(g_supercell(0,0,0,1:3,:)))
    write(*,*)'|G(force)| =',gnorm_supercell, ' supercell gradient for atoms in the unit cell '
    gnorm=sum(abs(g(1:3,1:n)))
    write(*,*)'|G(force)| =',gnorm, ' unit cell gradient'
    deallocate(xyz_hstep, ns)
    !
    return
    !
  end subroutine pbcgdisp_new


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine copyc6(fname,maxc,max_elem,c6ab,maxci, minc6,minc6list,maxc6,&
      & maxc6list)
    integer maxc,max_elem,maxci(max_elem),mima
    real(wp)  c6ab(max_elem,max_elem,maxc,maxc,3)
    character*(*) fname
    logical minc6,maxc6,minc6list(max_elem),maxc6list(max_elem)

    character*1  atmp
    character*80 btmp
    real(wp)  x,y,f,cn1,cn2,cmax,xx(10)
    integer iat,jat,i,n,l,j,k,il,iadr,jadr,nn,kk
    logical special

    call init_pars()
    c6ab=-1
    maxci=0
    ! read file
    kk=1
    !only use values for cn=minimum
    if (minc6.or.maxc6) then
      do i=1,94
        if (minc6list(i))then
          c6ab(i,:,1,:,2)=10000000.0
          c6ab(:,i,:,1,3)=10000000.0
        end if
      end do


      do nn=1,nlines
        special=.false.
        iat=int(pars(kk+1))
        jat=int(pars(kk+2))
        call limit(iat,jat,iadr,jadr)

        !only CN=minimum for iat
        if (minc6list(iat)) then
          special=.true.
          maxci(iat)=1
          maxci(jat)=max(maxci(jat),jadr)

          if (pars(kk+3).le.c6ab(iat,jat,1,jadr,2)) then

            c6ab(iat,jat,1,jadr,1)=pars(kk)
            c6ab(iat,jat,1,jadr,2)=pars(kk+3)
            c6ab(iat,jat,1,jadr,3)=pars(kk+4)

            c6ab(jat,iat,jadr,1,1)=pars(kk)
            c6ab(jat,iat,jadr,1,2)=pars(kk+4)
            c6ab(jat,iat,jadr,1,3)=pars(kk+3)
          end if
        end if

        !only CN=minimum for jat
        if (minc6list(jat)) then
          special=.true.
          maxci(iat)=max(maxci(iat),iadr)
          maxci(jat)=1

          if (pars(kk+4).le.c6ab(iat,jat,iadr,1,3)) then

            c6ab(iat,jat,iadr,1,1)=pars(kk)
            c6ab(iat,jat,iadr,1,2)=pars(kk+3)
            c6ab(iat,jat,iadr,1,3)=pars(kk+4)

            c6ab(jat,iat,1,iadr,1)=pars(kk)
            c6ab(jat,iat,1,iadr,2)=pars(kk+4)
            c6ab(jat,iat,1,iadr,3)=pars(kk+3)
          end if
        end if


        !only CN=minimum for
        if (minc6list(iat).and.minc6list(jat)) then
          special=.true.
          maxci(jat)=1
          maxci(iat)=1

          if (pars(kk+4).le.c6ab(iat,jat,1,1,3).and.&
              &      pars(kk+3).le.c6ab(iat,jat,1,1,2)) then

            c6ab(iat,jat,1,1,1)=pars(kk)
            c6ab(iat,jat,1,1,2)=pars(kk+3)
            c6ab(iat,jat,1,1,3)=pars(kk+4)

            c6ab(jat,iat,1,1,1)=pars(kk)
            c6ab(jat,iat,1,1,2)=pars(kk+4)
            c6ab(jat,iat,1,1,3)=pars(kk+3)
          end if
        end if



        !only CN=maximum for iat
        if (maxc6list(iat)) then
          special=.true.

          maxci(iat)=1
          maxci(jat)=max(maxci(jat),jadr)

          if (pars(kk+3).ge.c6ab(iat,jat,1,jadr,2)) then

            c6ab(iat,jat,1,jadr,1)=pars(kk)
            c6ab(iat,jat,1,jadr,2)=pars(kk+3)
            c6ab(iat,jat,1,jadr,3)=pars(kk+4)

            c6ab(jat,iat,jadr,1,1)=pars(kk)
            c6ab(jat,iat,jadr,1,2)=pars(kk+4)
            c6ab(jat,iat,jadr,1,3)=pars(kk+3)
          end if
        end if
        !only CN=maximum for jat
        if (maxc6list(jat)) then
          special=.true.

          maxci(jat)=1
          maxci(iat)=max(maxci(iat),iadr)

          if (pars(kk+4).ge.c6ab(iat,jat,iadr,1,3)) then

            c6ab(iat,jat,iadr,1,1)=pars(kk)
            c6ab(iat,jat,iadr,1,2)=pars(kk+3)
            c6ab(iat,jat,iadr,1,3)=pars(kk+4)

            c6ab(jat,iat,1,iadr,1)=pars(kk)
            c6ab(jat,iat,1,iadr,2)=pars(kk+4)
            c6ab(jat,iat,1,iadr,3)=pars(kk+3)
          end if
        end if

        !only CN=maximum for
        if (maxc6list(iat).and.maxc6list(jat)) then
          special=.true.
          maxci(jat)=1
          maxci(iat)=1

          if (pars(kk+4).ge.c6ab(iat,jat,1,1,3).and.&
              &      pars(kk+3).ge.c6ab(iat,jat,1,1,2)) then

            c6ab(iat,jat,1,1,1)=pars(kk)
            c6ab(iat,jat,1,1,2)=pars(kk+3)
            c6ab(iat,jat,1,1,3)=pars(kk+4)

            c6ab(jat,iat,1,1,1)=pars(kk)
            c6ab(jat,iat,1,1,2)=pars(kk+4)
            c6ab(jat,iat,1,1,3)=pars(kk+3)
          end if
        end if

        !only CN=minimum for
        if (minc6list(iat).and.maxc6list(jat)) then
          !and CN=maximum jat
          special=.true.
          maxci(jat)=1
          maxci(iat)=1

          if (pars(kk+4).ge.c6ab(iat,jat,1,1,3).and.&
              &      pars(kk+3).le.c6ab(iat,jat,1,1,2)) then

            c6ab(iat,jat,1,1,1)=pars(kk)
            c6ab(iat,jat,1,1,2)=pars(kk+3)
            c6ab(iat,jat,1,1,3)=pars(kk+4)

            c6ab(jat,iat,1,1,1)=pars(kk)
            c6ab(jat,iat,1,1,2)=pars(kk+4)
            c6ab(jat,iat,1,1,3)=pars(kk+3)
          end if
        end if

        !only CN=maximum for
        if (maxc6list(iat).and.minc6list(jat)) then
          !  and CN=minimum fo
          special=.true.
          maxci(jat)=1
          maxci(iat)=1

          if (pars(kk+4).le.c6ab(iat,jat,1,1,3).and.&
              &      pars(kk+3).ge.c6ab(iat,jat,1,1,2)) then

            c6ab(iat,jat,1,1,1)=pars(kk)
            c6ab(iat,jat,1,1,2)=pars(kk+3)
            c6ab(iat,jat,1,1,3)=pars(kk+4)

            c6ab(jat,iat,1,1,1)=pars(kk)
            c6ab(jat,iat,1,1,2)=pars(kk+4)
            c6ab(jat,iat,1,1,3)=pars(kk+3)
          end if
        end if

        if (.not.special) then

          maxci(iat)=max(maxci(iat),iadr)
          maxci(jat)=max(maxci(jat),jadr)

          c6ab(iat,jat,iadr,jadr,1)=pars(kk)
          c6ab(iat,jat,iadr,jadr,2)=pars(kk+3)
          c6ab(iat,jat,iadr,jadr,3)=pars(kk+4)

          c6ab(jat,iat,jadr,iadr,1)=pars(kk)
          c6ab(jat,iat,jadr,iadr,2)=pars(kk+4)
          c6ab(jat,iat,jadr,iadr,3)=pars(kk+3)
        end if
        kk=(nn*5)+1
      end do



      !no min/max at all
    else
      do nn=1,nlines
        iat=int(pars(kk+1))
        jat=int(pars(kk+2))
        call limit(iat,jat,iadr,jadr)
        maxci(iat)=max(maxci(iat),iadr)
        maxci(jat)=max(maxci(jat),jadr)

        c6ab(iat,jat,iadr,jadr,1)=pars(kk)
        c6ab(iat,jat,iadr,jadr,2)=pars(kk+3)
        c6ab(iat,jat,iadr,jadr,3)=pars(kk+4)

        c6ab(jat,iat,jadr,iadr,1)=pars(kk)
        c6ab(jat,iat,jadr,iadr,2)=pars(kk+4)
        c6ab(jat,iat,jadr,iadr,3)=pars(kk+3)
        kk=(nn*5)+1
      end do
    end if

  end subroutine copyc6



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SET_CRITERIA(rthr,lat,tau_max)

    REAL(WP) :: r_cutoff,rthr
    REAL(WP) :: lat(3,3)
    REAL(WP) :: tau_max(3)
    REAL(WP) :: norm1(3),norm2(3),norm3(3)
    REAL(WP) :: cos10,cos21,cos32

    r_cutoff=sqrt(rthr)
    ! write(*,*) 'lat',lat
    call kreuzprodukt(lat(:,2),lat(:,3),norm1)
    call kreuzprodukt(lat(:,3),lat(:,1),norm2)
    call kreuzprodukt(lat(:,1),lat(:,2),norm3)
    ! write(*,*) 'norm2',norm2
    norm1=norm1/VECTORSIZE(norm1)
    norm2=norm2/VECTORSIZE(norm2)
    norm3=norm3/VECTORSIZE(norm3)
    ! write(*,*) 'norm2_',norm2
    cos10=SUM(norm1*lat(:,1))
    cos21=SUM(norm2*lat(:,2))
    cos32=SUM(norm3*lat(:,3))
    tau_max(1)=abs(r_cutoff/cos10)
    tau_max(2)=abs(r_cutoff/cos21)
    tau_max(3)=abs(r_cutoff/cos32)
    ! write(*,'(3f8.4)')tau_max(1),tau_max(2),tau_max(3)
  end subroutine SET_CRITERIA


  subroutine kreuzprodukt(A,B,C)

    REAL(WP) :: A(3),B(3)
    REAL(WP) :: X,Y,Z
    REAL(WP) :: C(3)

    X=A(2)*B(3)-B(2)*A(3)
    Y=A(3)*B(1)-B(3)*A(1)
    Z=A(1)*B(2)-B(1)*A(2)
    C=(/X,Y,Z/)
  end subroutine kreuzprodukt

  function VECTORSIZE(VECT)

    REAL(WP) :: VECT(3)
    REAL(WP) :: SVECT(3)
    REAL(WP) :: VECTORSIZE

    SVECT=VECT*VECT
    VECTORSIZE=SUM(SVECT)
    VECTORSIZE=VECTORSIZE**(0.5)
  end function VECTORSIZE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function determinant(aa) result(det)
    real(wp), intent(in) :: aa(:,:)
    real(wp) :: det

    det = aa(1,1) * (aa(2,2) * aa(3,3) - aa(3,2) * aa(2,3))&
        & - aa(1,2) * (aa(2,1) * aa(3,3) - aa(3,1) * aa(2,3))&
        & + aa(1,3) * (aa(2,1) * aa(3,2) - aa(3,1) * aa(2,2))
    
  end function determinant


  subroutine inv_cell(x,a)
    real(wp), intent(in) :: x(3,3)
    real(wp), intent(out) :: a(3,3)
    integer i
    real(wp) det

    a = 0.0
    det = determinant(x)
    ! write(*,*)'Det:',det
    a(1,1)=x(2,2)*x(3,3)-x(2,3)*x(3,2)
    a(2,1)=x(2,3)*x(3,1)-x(2,1)*x(3,3)
    a(3,1)=x(2,1)*x(3,2)-x(2,2)*x(3,1)
    a(1,2)=x(1,3)*x(3,2)-x(1,2)*x(3,3)
    a(2,2)=x(1,1)*x(3,3)-x(1,3)*x(3,1)
    a(3,2)=x(1,2)*x(3,1)-x(1,1)*x(3,2)
    a(1,3)=x(1,2)*x(2,3)-x(1,3)*x(2,2)
    a(2,3)=x(1,3)*x(2,1)-x(1,1)*x(2,3)
    a(3,3)=x(1,1)*x(2,2)-x(1,2)*x(2,1)
    a=a/det
  end subroutine inv_cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xyz_to_abc(xyz,abc,lat,n)
    integer,intent(in) :: n
    real(wp), INTENT(in) :: xyz(3,n)
    real(wp), intent(in) :: lat(3,3)
    real(wp), intent(out) :: abc(3,n)

    real(wp) lat_1(3,3)
    integer i,j,k

    call inv_cell(lat,lat_1)

    abc(:,:n)=0.0d0
    do i=1,n
      do j=1,3
        do k=1,3
          abc(j,i)=abc(j,i)+lat_1(j,k)*xyz(k,i)
        end do
        abc(j,i)=dmod(abc(j,i),1.0d0)
      end do
    end do

  end subroutine xyz_to_abc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine abc_to_xyz(abc,xyz,lat,n)
    real(wp), INTENT(in) :: abc(3,*)
    real(wp), intent(in) :: lat(3,3)
    real(wp), intent(out) :: xyz(3,*)
    integer,intent(in) :: n

    integer i,j,k

    xyz(:,:n)=0.0d0
    do i=1,n
      do j=1,3
        do k=1,3
          xyz(j,i)=xyz(j,i)+lat(j,k)*abc(k,i)
        end do
      end do
    end do

  end subroutine abc_to_xyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gkernel1 ( version, r2, R0, s6, rs6, alp6, s8, rs8, alp8, C6, r42, &
                          fact, res1, res2 ) 
  implicit none
  integer, intent(in)   :: version
  real(wp), intent(in)  :: fact ! overall scaling factor for res1 and res2
  real(wp), intent(in)  :: r2, R0, s6, rs6, alp6, s8, rs8, alp8, C6, r42
  real(wp), intent(out) :: res1, res2

  real(wp) :: r, r4, r6, r7, r8, r9, tmp1, tmp2
  real(wp) :: t6, damp6, t8, damp8 , dc6_rest

  r=dsqrt(r2)
  r4=r2*r2
  r6=r4*r2
  r7=r6*r
  r8=r6*r2
  r9=r8*r

  if (version.eq.3) then
    t6 = (r/(rs6*R0))**(-alp6)
    damp6 =1.d0/( 1.d0+6.d0*t6 )
    t8 = (r/(rs8*R0))**(-alp8)
    damp8 =1.d0/( 1.d0+6.d0*t8 )
    res1 = -s6*(6.0/(r7)*C6*damp6) -s8*(24.0/(r9)*C6*r42*damp8) &
            +s6*C6/r7*6.d0*alp6*t6*damp6*damp6 +s8*C6*r42/r9*18.d0*alp8*t8*damp8*damp8 
    res2 = s6/r6*damp6+3.d0*s8*r42/r8*damp8 

  elseif(version.eq.5) then  
    t6 = (r/(rs6*R0)+R0*rs8)**(-alp6)
    damp6 =1.d0/( 1.d0+6.d0*t6 )
    t8 = (r/(R0)+R0*rs8)**(-alp8)
    damp8 =1.d0/( 1.d0+6.d0*t8 )
    tmp1=s6*6.d0*damp6*C6/r7
    tmp2=s8*6.d0*C6*r42*damp8/r9
    res1 = - (tmp1 +4.d0*tmp2) +(tmp1*alp6*t6*damp6*r/(r+rs6*R0*R0*rs8) +3.d0*tmp2*alp8*t8*damp8*r/(r+R0*R0*rs8)) 
    res2 = s6/r6*damp6+3.d0*s8*r42/r8*damp8 

  elseif((version.eq.4).or.(version.eq.6)) then  
    t6=(r6+R0**6)
    t8=(r8+R0**8)
    res1 = -s6*C6*6.0d0*r4*r/(t6*t6) -s8*C6*24.0d0*r42*r7/(t8*t8) 
    res2 = s6/t6+3.d0*s8*r42/t8

  endif

  res1 = fact * res1
  res2 = fact * res2

  return

  end subroutine gkernel1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gkernel2 ( rij, r2, crit_cn, rcovij, drij_linij, dc6i_iat, dc6i_jat, vec ) 
  implicit none
  real(wp), intent(in)  :: rij(3) 
  real(wp), intent(in)  :: r2, crit_cn, rcovij, drij_linij, dc6i_iat, dc6i_jat
  real(wp), intent(out) :: vec(3) 

  real(wp) :: r, expterm, dcnn, x1
  
  r=dsqrt(r2)

  if (r2.lt.crit_cn) then
    expterm=exp(-k1*(rcovij/r-1.d0))
    dcnn=-k1*rcovij*expterm/&
        & (r2*(expterm+1.d0)*(expterm+1.d0))
  else
    dcnn=0.0d0
  end if

  x1=drij_linij + dcnn*(dc6i_iat + dc6i_jat )
  vec=x1*rij/r

  return 

  end subroutine gkernel2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine gkernel3 ( fact, dxyz, r2, alp6, R0, c6, vec )
  implicit none
  real(wp), intent(in) :: fact, dxyz(3), r2, alp6, R0, c6
  real(wp), intent(out) :: vec(3)

  real(wp) :: r235, r, damp6, damp1, tmp1, tmp2, term

  r235=r2**3.5
  r =dsqrt(r2)
  damp6=exp(-alp6*(r/R0-1.0d0))
  damp1=1.+damp6
  tmp1=damp6/(damp1*damp1*r235*R0)
  tmp2=6./(damp1*r*r235)

  term=alp6*tmp1-tmp2

  vec(1:3) = fact * term*dxyz(1:3)*c6

  return
  
  end subroutine gkernel3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module dftd3_core
