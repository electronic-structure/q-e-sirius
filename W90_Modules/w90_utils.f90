!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE w90_utils
  !----------------------------------------------------------------------------
  !! Module for writing real-space Wannier functions.
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  SAVE
  !
  INTERFACE w90_resize_in_place
    !! Resize a (n1, ..., nN) array in place by adding n_add indices to the last index.
    !! Input: (n1, ..., nN) array
    !! Output: (n1, ..., nN + n_add) array
    !! The data of the original array in the indices 1 to nN are kept.
    !! The new indices are uninitialized.
    MODULE PROCEDURE w90_resize_in_place_i2, w90_resize_in_place_r3
  END INTERFACE
  !
  CONTAINS
  !
  !------------------------------------------------------------------------------------
  SUBROUTINE w90_resize_in_place_i2(arr, n_add)
    !----------------------------------------------------------------------------
    !! Resize the array in place by adding n_add columns.
    !! The data of the original array are kept.
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: arr(:, :)
    !! Array to be resized.
    !! Input: size (n1, n2)
    !! Output: size (n1, n2, n + n_add)
    INTEGER, INTENT(IN) :: n_add
    !! Number of columns to add
    !
    INTEGER :: n1, n2
    !! Size of the original array
    INTEGER, ALLOCATABLE :: arr_copy(:, :)
    !! Temporary array to hold the original data
    !
    n1 = SIZE(arr, 1)
    n2 = SIZE(arr, 2)
    ALLOCATE(arr_copy(n1, n2 + n_add))
    arr_copy(:, 1:n2) = arr(:, 1:n2)
    !
    CALL move_alloc(arr_copy, arr)  ! Move allocation from arr_copy to arr
    !
    ! More verbose form with an explicit copy that could replace the above line :
    ! DEALLOCATE(arr)
    ! ALLOCATE(arr(n1, n2 + n_add))
    ! arr(:, :) = arr_copy(:, :)
    ! DEALLOCATE(arr_copy)
    !
  !------------------------------------------------------------------------------------
  END SUBROUTINE w90_resize_in_place_i2
  !------------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------------
  SUBROUTINE w90_resize_in_place_r3(arr, n_add)
    USE kinds, ONLY : DP
    IMPLICIT NONE
    !
    REAL(DP), ALLOCATABLE, INTENT(INOUT) :: arr(:, :, :)
    INTEGER, INTENT(IN) :: n_add
    !
    INTEGER :: n1, n2, n3
    REAL(DP), ALLOCATABLE :: arr_copy(:, :, :)
    !
    n1 = SIZE(arr, 1)
    n2 = SIZE(arr, 2)
    n3 = SIZE(arr, 3)
    ALLOCATE(arr_copy(n1, n2, n3 + n_add))
    arr_copy(:, :, 1:n3) = arr(:, :, 1:n3)
    CALL move_alloc(arr_copy, arr)  ! Move allocation from arr_copy to arr
  !------------------------------------------------------------------------------------
  END SUBROUTINE w90_resize_in_place_r3
  !------------------------------------------------------------------------------------
END MODULE w90_utils
