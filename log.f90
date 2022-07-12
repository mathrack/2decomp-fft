!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
! Copyright (C) 2021               the University of Edinburgh (UoE)
!
!=======================================================================

submodule (decomp_2d) smod_log

  implicit none

  contains

  !
  ! Generate a IO unit for the decomp_2d listing
  !
  module function d2d_get_iounit() result(iounit)

     use iso_fortran_env, only : output_unit

     implicit none

     integer :: iounit
     integer :: ierror
     character(10) :: fname

#ifdef DEBUG
    if (nrank_loc >= 0) then
       write(fname, "(I3.3,'_',I3.3)") nrank, nrank_loc
    else
       write(fname, "(I3.3)") nrank
    endif
    open(newunit=iounit, file='decomp_2d_setup_'//trim(fname)//'.log', iostat=ierror)
#else
    if (nrank == 0) then
       open(newunit=iounit, file="decomp_2d_setup.log", iostat=ierror)
    else
       iounit = output_unit
       ierror = 0
    endif
#endif
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "Could not open log file")

  end function d2d_get_iounit
  !
  ! Print some information about decomp_2d
  !
  module subroutine d2d_listing(given_io_unit)

    use iso_fortran_env, only : output_unit, compiler_version, compiler_options

    implicit none

    ! Argument
    integer, intent(in), optional :: given_io_unit

    ! Local variable
    integer :: io_unit, ierror

    !
    ! Default : only rank 0 (or local master on node 0) will print a listing
    !
    ! In DEBUG mode, all ranks will print a listing
    !
#ifndef DEBUG
    if (nrank /= 0 .or. nrank_loc > 0) return
#endif

    ! If no IO unit provided, use stdout
    if (present(given_io_unit)) then
            io_unit = given_io_unit
    else
            io_unit = output_unit
    endif

    ! Header
    write (io_unit, *) '==========================================================='
    write (io_unit, *) '=================== Decomp2D - log ========================'
    write (io_unit, *) '==========================================================='

    ! Git hash if available
#if defined(VERSION)
    write (io_unit, *) 'Git version        : ', VERSION
#else
    write (io_unit, *) 'Git version        : unknown'
#endif

    ! Basic info
#ifdef DEBUG
    write (io_unit, *) 'I am mpi rank ', nrank
    write (io_unit, *) 'Local rank on the node ', nrank_loc
#endif
    write (io_unit, *) 'Total ranks ', nproc
    write (io_unit, *) 'Total ranks inside the node ', nproc_loc
    write (io_unit, *) 'Global data size : ', nx_global, ny_global, nz_global
    write (io_unit, *) 'p_row, p_col : ', dims(1), dims(2)
    write (io_unit, *) 'Periodicity : ', periodic_x, periodic_y, periodic_z
    write (io_unit, *) 'Intranode : ', d2d_intranode
    write (io_unit, *) 'Number of bytes / float number : ', mytype_bytes
    write (io_unit, *) '==========================================================='

    ! Show detected flags, compiler options, version of the MPI library
    write (io_unit, *) 'Compile flags detected :'
#ifdef DOUBLE_PREC
#ifdef SAVE_SINGLE
    write (io_unit, *) 'Numerical precision: Double, saving in single'
#else
    write (io_unit, *) 'Numerical precision: Double'
#endif
#else
    write (io_unit, *) 'Numerical precision: Single'
#endif
    write (io_unit, *) 'Compiled with ', compiler_version()
    write (io_unit, *) 'Compiler options : ', compiler_options()
    write (io_unit, '(" Version of the MPI library : ",I0,".",I0)') MPI_VERSION, MPI_SUBVERSION
#ifdef DEBUG
    write (io_unit, *) 'Compile flag DEBUG detected'
#endif
#ifdef EVEN
    write (io_unit, *) 'Compile flag EVEN detected'
#endif
#ifdef OCC
    write (io_unit, *) 'Compile flag OCC detected'
#endif
#ifdef OVERWRITE
    write (io_unit, *) 'Compile flag OVERWRITE detected'
#endif
#ifdef HALO_DEBUG
    write (io_unit, *) 'Compile flag HALO_DEBUG detected'
#endif
#ifdef _GPU
    write (io_unit, *) 'Compile flag _GPU detected'
#endif
#ifdef _NCCL
    write (io_unit, *) 'Compile flag _NCCL detected'
#endif
    write (io_unit, *) '==========================================================='
    ! Info about each decomp_info object
    call decomp_info_print(decomp_main, io_unit, "decomp_main")
    call decomp_info_print(phG, io_unit, "phG")
    call decomp_info_print(ph1, io_unit, "ph1")
    call decomp_info_print(ph2, io_unit, "ph2")
    call decomp_info_print(ph3, io_unit, "ph3")
    call decomp_info_print(ph4, io_unit, "ph4")
    write (io_unit, *) '==========================================================='
    write (io_unit, *) '==========================================================='

    ! Close the unit if needed
    if (io_unit /= output_unit) then
       close(io_unit, iostat=ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "Could not close log file")
    endif

  end subroutine d2d_listing

  !
  ! Print some information about given decomp_info object
  !
  module subroutine decomp_info_print(d2d, io_unit, d2dname)

    implicit none

    ! Arguments
    type(decomp_info), intent(in) :: d2d
    integer, intent(in) :: io_unit
    character(len=*), intent(in) :: d2dname

    ! Nothing to print if not initialized
    if (.not.allocated(d2d%x1dist)) then
      write (io_unit, *) 'Uninitialized decomp_info ', d2dname
      return
    endif

    !
    ! If DEBUG mode, print everything
    ! Otherwise, print only global size
    !
    write (io_unit, *) 'Decomp_info : ', d2dname
    write (io_unit, *) '   Global size : ', d2d%xsz(1), d2d%ysz(1), d2d%zsz(1)
#ifdef DEBUG
    write (io_unit, *) '   xsz, xst, xen : ', d2d%xsz, d2d%xst, d2d%xen
    write (io_unit, *) '   ysz, yst, yen : ', d2d%ysz, d2d%yst, d2d%yen
    write (io_unit, *) '   zsz, zst, zen : ', d2d%zsz, d2d%zst, d2d%zen
    write (io_unit, *) '   xsz_loc, xst_loc, xen_loc : ', d2d%xsz_loc, d2d%xst_loc, d2d%xen_loc
    write (io_unit, *) '   ysz_loc, yst_loc, yen_loc : ', d2d%ysz_loc, d2d%yst_loc, d2d%yen_loc
    write (io_unit, *) '   zsz_loc, zst_loc, zen_loc : ', d2d%zsz_loc, d2d%zst_loc, d2d%zen_loc
    write (io_unit, *) '   x1dist : ', d2d%x1dist
    write (io_unit, *) '   y1dist : ', d2d%y1dist
    write (io_unit, *) '   y2dist : ', d2d%y2dist
    write (io_unit, *) '   z2dist : ', d2d%z2dist
    write (io_unit, *) '   x1cnts : ', d2d%x1cnts
    write (io_unit, *) '   y1cnts : ', d2d%y1cnts
    write (io_unit, *) '   y2cnts : ', d2d%y2cnts
    write (io_unit, *) '   z2cnts : ', d2d%z2cnts
    write (io_unit, *) '   x1disp : ', d2d%x1disp
    write (io_unit, *) '   y1disp : ', d2d%y1disp
    write (io_unit, *) '   y2disp : ', d2d%y2disp
    write (io_unit, *) '   z2disp : ', d2d%z2disp
    write (io_unit, *) '   x1count : ', d2d%x1count
    write (io_unit, *) '   y1count : ', d2d%y1count
    write (io_unit, *) '   y2count : ', d2d%y2count
    write (io_unit, *) '   z2count : ', d2d%z2count
    write (io_unit, *) '   even : ', d2d%even
#endif

  end subroutine decomp_info_print

end submodule smod_log
