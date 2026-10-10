# To use wannier to interpolate wavefunction
* The following code should be run for setup. (only once)
    ```fortran
    USE w90_interface,     ONLY : w90_read
    USE w90_real_space,    ONLY : w90_write_real_space
    USE w90_wan_Rr_buffer, ONLY : w90_wan_Rr_save_buffer
    !
    ! Read Wannier90 information from checkpoint file
    !
    CALL w90_read(seedname)
    !
    ! Compute and write real-space Wannier functions to file in collected form
    !
    IF (write_wan_Rr) CALL w90_write_real_space(tmp_dir_wfc, folder_wan_Rr)
    !
    ! Read the collected real-space Wannier functions from file to buffer
    !
    CALL w90_wan_Rr_save_buffer(folder_wan_Rr, tmp_dir_buffer)
    !
    ```
    * `tmp_dir_wfc`: `outdir` for the NSCF calculations (containing `prefix.save/wfc###.dat` files for Wannierization)
    * `folder_wan_Rr`: folder where collected real-space Wannier function data are stored. (This data is transferrable.)
    * `tmp_dir_buffer`: temporary directory where the buffers will be stored.
* Then, use `CALL w90_interpolate_wfc(ik, evc_kp)`. It uses k point data from module `klist`.
* A buffer unit `iuwan_r = 100` (in module `w90_wan_Rr_buffer`) is opened in `w90_wan_Rr_save_buffer`. It must be kept open.
    * Unit collision may happen...
* To close the buffer, run `CALL w90_wan_Rr_close_buffer()`
* `io_level` in module `control_flags` is used to control I/O vs memory.


# TODO
* Use buffer for k point wavefunctions in the first block of `w90_write_real_space`. Make it a separate subroutine.
* Wrap the above setup blocks to a single subroutine.
* Allow reuse of buffer. (For pool parallelization, different pools can share same `prefix.wan_r#` buffer.)
    * This is doable, something like below. However, this conflicts too much with how buffer library of QE works. So I did not implement this.
    ```fortran
    !
    ! The wan_Rr data is distributed over plane waves, but are identical over k points and q
    ! points. Thus, we can reuse the same buffer for all pools and images. This reduces
    ! the disk usage. We set write_buffer to .FALSE. if the buffer is not written, only read.
    !
    write_buffer = .TRUE.
    !
    ! If io_level == 0, the buffer is not created and everything is stored in memory, so
    ! we must create the buffer. So we set write_buffer to .FALSE. only if io_level > 0.
    !
    IF (io_level > 0 .AND. ((my_pool_id /= 0) .OR. (my_image_id /= 0))) THEN
      write_buffer = .FALSE.
    ENDIF
    !
    iuwan_r = 100
    lrwan_r = dffts%nnr * npol * num_wann
    !
    IF (write_buffer) THEN
      !
      ! Open a file to write the Wannier functions in real space
      !
      CALL open_buffer(iuwan_r, 'wan_r', lrwan_r, io_level, exst_mem, exst, tmp_dir_buffer)
      !
      ! Read wan_Rr wavefunction from file and write to the buffer
      !
      ALLOCATE(wan_Rr(dffts%nnr, npol, num_wann))
      !
      DO iR = 1, nR_ws
        CALL w90_readwrite_wan_Rr(iRlist_ws(:, iR), wan_Rr, -1, tmp_dir_collected)
        CALL save_buffer(wan_Rr, lrwan_r, iuwan_r, iR)
      ENDDO ! iR
      !
      DEALLOCATE(wan_Rr)
      !
    ENDIF
    !
    ! Wait for all buffers to be written before proceeding
    !
    CALL mp_barrier(world_comm)
    !
    IF (.NOT. write_buffer) THEN
      !
      ! Open the buffer written by the first image/pool.
      !
      ! Taken from diropn in Modules/io_files.f90, changed the postfix nd_nmbr to be
      ! the index of the processor in the intra_bgrp_comm.
      !
      filename = TRIM(tmp_dir_buffer) // TRIM(prefix) // ".wan_r" // TRIM(int_to_char(me_bgrp+1))
      !
      INQUIRE(FILE=filename, EXIST=exst)
      !
      IF (.NOT. exst) CALL errore('w90_wan_Rr_get_buffer', &
          'buffer file '// TRIM(filename) // ' does not exist', 1)
      !
      ! the  record length in direct-access I/O is given by the number of
      ! real*8 words times direct_io_factor (may depend on the compiler)
      !
      INQUIRE(IOLENGTH=direct_io_factor) dummy
      unf_recl = direct_io_factor * INT(2*lrwan_r, KIND=KIND(unf_recl))
      !
      OPEN(UNIT=iuwan_r, FILE=TRIM(ADJUSTL(filename)), IOSTAT=ios, FORM='unformatted', &
          STATUS='old', ACTION='read', ACCESS='direct', RECL=unf_recl)
      !
      IF (ios /= 0) call errore('w90_wan_Rr_save_buffer', &
          'error opening ' // TRIM(filename), ABS(ios))
      !
    ENDIF
    ```
