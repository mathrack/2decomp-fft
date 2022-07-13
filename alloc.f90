!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2021 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This is the submodule dedicated to routines allocating memory

submodule (decomp_2d) smod_alloc

  use, intrinsic :: iso_c_binding, only : c_loc, c_ptr, c_f_pointer

  implicit none

  contains

  !
  ! Generic subroutine to allocate a global 3D real array
  !
  subroutine alloc_glob_real(var, s1, e1, s2, e2, s3, e3)

    implicit none

    ! Arguments
    integer, intent(in) :: s1, e1, s2, e2, s3, e3
    real(mytype), dimension(:,:,:), pointer :: var

    ! Local variables
    integer :: alloc_stat, errorcode

    if (associated(var)) nullify(var)
    allocate(var(s1:e1,s2:e2,s3:e3), stat=alloc_stat)

    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

  end subroutine alloc_glob_real

  !
  ! Generic subroutine to allocate a global 3D complex array                                            
  !
  subroutine alloc_glob_cplx(var, s1, e1, s2, e2, s3, e3)

    implicit none

    ! Arguments
    integer, intent(in) :: s1, e1, s2, e2, s3, e3
    complex(mytype), dimension(:,:,:), pointer :: var

    ! Local variables
    integer :: alloc_stat, errorcode

    if (associated(var)) nullify(var)
    allocate(var(s1:e1,s2:e2,s3:e3), stat=alloc_stat)

    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

  end subroutine alloc_glob_cplx

  !
  ! Generic subroutine to allocate a local 3D real array
  !
  subroutine alloc_loc_real(var, n1, n2, n3, var2d)

    implicit none

    ! Arguments
    integer, intent(in) :: n1, n2 ,n3
    real(mytype), dimension(:,:,:), pointer :: var
    real(mytype), dimension(:,:), pointer, optional :: var2d

    ! Local variables
    integer :: alloc_stat, errorcode

    ! Allocate memory
    if (associated(var)) nullify(var)
    allocate(var(n1, n2, n3), stat=alloc_stat)
    if (alloc_stat /= 0) then 
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

    ! Also provide a 2D map if needed
    if (present(var2d)) then
       if (associated(var2d)) nullify(var2d)
       call c_f_pointer(c_loc(var), var2d, (/n1, n2 * n3/))
    endif

  end subroutine alloc_loc_real

  !
  ! Generic subroutine to allocate a local 3D complex array
  !
  subroutine alloc_loc_cplx(var, n1, n2, n3, var2d)

    implicit none

    ! Arguments
    integer, intent(in) :: n1, n2 ,n3
    complex(mytype), dimension(:,:,:), pointer :: var 
    complex(mytype), dimension(:,:), pointer, optional :: var2d

    ! Local variables
    integer :: alloc_stat, errorcode

    ! Allocate memory
    if (associated(var)) nullify(var)
    allocate(var(n1, n2, n3), stat=alloc_stat)
    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

    ! Also provide a 2D map if needed
    if (present(var2d)) then
       if (associated(var2d)) nullify(var2d)
       call c_f_pointer(c_loc(var), var2d, (/n1, n2 * n3/))
    endif

  end subroutine alloc_loc_cplx

  !
  ! Generic subroutine to allocate a MPI3 unified memory real array
  !
  subroutine alloc_mpi3_real(var, n1, n2, n3, var2d, nloc1, nloc2, nloc3, win)

    implicit none

    ! Arguments
    integer, intent(out) :: win
    integer, intent(in) :: n1, n2, n3, nloc1, nloc2, nloc3
    real(mytype), dimension(:,:,:), pointer :: var
    real(mytype), dimension(:,:), pointer :: var2d

    ! Local variables
    type(c_ptr) :: baseptr
    integer(kind=MPI_ADDRESS_KIND) :: winsize, tmpsize, tmpdispunit
    integer :: info, ierror

    ! Size of the memory to allocate on each CPU
    winsize = nloc1 * nloc2 * nloc3

    ! Same disp_unit everywhere
    call MPI_INFO_CREATE(info, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_INFO_CREATE")
    call MPI_INFO_SET(info, "same_disp_unit", "true", ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_INFO_SET")

    ! Create a window and allow direct memory access inside DECOMP_2D_LOCALCOMM
    win = MPI_WIN_NULL
    call MPI_WIN_ALLOCATE_SHARED(winsize * mytype_bytes, &
                                 mytype_bytes, &
                                 info, &
                                 DECOMP_2D_LOCALCOMM, &
                                 baseptr, &
                                 win, &
                                 ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_ALLOCATE_SHARED")

    ! Get the CPU-level 2D memory
    if (associated(var2d)) nullify(var2d)
    call C_F_POINTER(baseptr, var2d, (/nloc1, nloc2 * nloc3/))

    ! Get the node-level 3D memory
    call MPI_WIN_SHARED_QUERY(win, 0, tmpsize, tmpdispunit, baseptr, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_SHARED_QUERY")
    if (associated(var)) nullify(var)
    call C_F_POINTER(baseptr, var, (/n1, n2, n3/))

    call MPI_INFO_FREE(info, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_INFO_FREE")

  end subroutine alloc_mpi3_real

  !
  ! Generic subroutine to allocate a MPI3 unified memory complex array
  !
  subroutine alloc_mpi3_cplx(var, n1, n2, n3, var2d, nloc1, nloc2, nloc3, win)

    implicit none

    ! Arguments
    integer, intent(out) :: win
    integer, intent(in) :: n1, n2, n3, nloc1, nloc2, nloc3
    complex(mytype), dimension(:,:,:), pointer :: var
    complex(mytype), dimension(:,:), pointer :: var2d

    ! Local variables
    type(c_ptr) :: baseptr
    integer(kind=MPI_ADDRESS_KIND) :: winsize, tmpsize, tmpdispunit
    integer :: info, ierror

    ! Size of the memory to allocate on each CPU
    winsize = nloc1 * nloc2 * nloc3

    ! Same disp_unit everywhere
    call MPI_INFO_CREATE(info, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_INFO_CREATE")
    call MPI_INFO_SET(info, "same_disp_unit", "true", ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_INFO_SET")

    ! Create a window and allow direct memory access inside DECOMP_2D_LOCALCOMM
    win = MPI_WIN_NULL
    call MPI_WIN_ALLOCATE_SHARED(winsize * mytype_bytes * 2, &
                                 mytype_bytes * 2, &
                                 info, &
                                 DECOMP_2D_LOCALCOMM, &
                                 baseptr, &
                                 win, &
                                 ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_ALLOCATE_SHARED")
       
    ! Get the CPU-level 2D memory
    if (associated(var2d)) nullify(var2d)
    call C_F_POINTER(baseptr, var2d, (/nloc1, nloc2 * nloc3/))
                                    
    ! Get the node-level 3D memory  
    call MPI_WIN_SHARED_QUERY(win, 0, tmpsize, tmpdispunit, baseptr, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_SHARED_QUERY")
    if (associated(var)) nullify(var)
    call C_F_POINTER(baseptr, var, (/n1, n2, n3/))

    call MPI_INFO_FREE(info, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_INFO_FREE")

  end subroutine alloc_mpi3_cplx

  ! X-pencil real arrays
  module subroutine alloc_x_real(var, var2d, opt_decomp, opt_global, win)

    implicit none

    real(mytype), dimension(:,:,:), pointer :: var
    real(mytype), dimension(:,:), pointer, optional :: var2d
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global
    integer, intent(out), optional :: win

    TYPE(DECOMP_INFO) :: decomp
    logical :: global

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    ! Safety check
#ifdef DEBUG
    if (d2d_intranode.and.(global.or..not.present(win))) &
       call decomp_2d_abort(__FILE__, __LINE__, -1, "Incorrect alloc arguments")
#endif

    if (d2d_intranode) then
       ! MPI3 shared memory, the caller should free the MPI window
       call alloc_mpi3_real(var, decomp%xsz(1), decomp%xsz(2), decomp%xsz(3), &
                            var2d, decomp%xsz_loc(1), decomp%xsz_loc(2), decomp%xsz_loc(3), win)
    else
       ! No MPI3 shared memory
       if (global) then
          call alloc_glob_real(var, decomp%xst(1), decomp%xen(1), &
                                    decomp%xst(2), decomp%xen(2), &
                                    decomp%xst(3), decomp%xen(3))
       else
          call alloc_loc_real(var, decomp%xsz(1), decomp%xsz(2), decomp%xsz(3), var2d)
       end if
    endif
    
  end subroutine alloc_x_real

  ! X-pencil complex arrays
  module subroutine alloc_x_complex(var, var2d, opt_decomp, opt_global, win)

    implicit none

    complex(mytype), dimension(:,:,:), pointer :: var
    complex(mytype), dimension(:,:), pointer, optional :: var2d
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global
    integer, intent(out), optional :: win

    TYPE(DECOMP_INFO) :: decomp
    logical :: global

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    ! Safety check
#ifdef DEBUG
    if (d2d_intranode.and.(global.or..not.present(win))) &
       call decomp_2d_abort(__FILE__, __LINE__, -1, "Incorrect alloc arguments")
#endif

    if (d2d_intranode) then
       ! MPI3 shared memory, the caller should free the MPI window
       call alloc_mpi3_cplx(var, decomp%xsz(1), decomp%xsz(2), decomp%xsz(3), &
                            var2d, decomp%xsz_loc(1), decomp%xsz_loc(2), decomp%xsz_loc(3), win)
    else
       ! No MPI3 shared memory
       if (global) then
          call alloc_glob_cplx(var, decomp%xst(1), decomp%xen(1), &
                                    decomp%xst(2), decomp%xen(2), &
                                    decomp%xst(3), decomp%xen(3))
       else
          call alloc_loc_cplx(var, decomp%xsz(1), decomp%xsz(2), decomp%xsz(3), var2d)
       end if
    endif
    
  end subroutine alloc_x_complex

  ! Y-pencil real arrays
  module subroutine alloc_y_real(var, var2d, opt_decomp, opt_global, win)

    implicit none

    real(mytype), dimension(:,:,:), pointer :: var
    real(mytype), dimension(:,:), pointer, optional :: var2d
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global
    integer, intent(out), optional :: win

    TYPE(DECOMP_INFO) :: decomp
    logical :: global

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    ! Safety check
#ifdef DEBUG
    if (d2d_intranode.and.(global.or..not.present(win))) &
       call decomp_2d_abort(__FILE__, __LINE__, -1, "Incorrect alloc arguments")
#endif

    if (d2d_intranode) then
       ! MPI3 shared memory, the caller should free the MPI window
       call alloc_mpi3_real(var, decomp%ysz(1), decomp%ysz(2), decomp%ysz(3), &
                            var2d, decomp%ysz_loc(1), decomp%ysz_loc(2), decomp%ysz_loc(3), win)
    else
       ! No MPI3 shared memory
       if (global) then
          call alloc_glob_real(var, decomp%yst(1), decomp%yen(1), &
                                    decomp%yst(2), decomp%yen(2), &
                                    decomp%yst(3), decomp%yen(3))
       else
          call alloc_loc_real(var, decomp%ysz(1), decomp%ysz(2), decomp%ysz(3), var2d)
       end if
    endif
    
  end subroutine alloc_y_real

  ! Y-pencil complex arrays
  module subroutine alloc_y_complex(var, var2d, opt_decomp, opt_global, win)

    implicit none

    complex(mytype), dimension(:,:,:), pointer :: var
    complex(mytype), dimension(:,:), pointer, optional :: var2d
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global
    integer, intent(out), optional :: win

    TYPE(DECOMP_INFO) :: decomp
    logical :: global

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    ! Safety check
#ifdef DEBUG
    if (d2d_intranode.and.(global.or..not.present(win))) &
       call decomp_2d_abort(__FILE__, __LINE__, -1, "Incorrect alloc arguments")
#endif

    if (d2d_intranode) then
       ! MPI3 shared memory, the caller should free the MPI window
       call alloc_mpi3_cplx(var, decomp%ysz(1), decomp%ysz(2), decomp%ysz(3), &
                            var2d, decomp%ysz_loc(1), decomp%ysz_loc(2), decomp%ysz_loc(3), win)
    else
       ! No MPI3 shared memory
       if (global) then
          call alloc_glob_cplx(var, decomp%yst(1), decomp%yen(1), &
                                    decomp%yst(2), decomp%yen(2), &
                                    decomp%yst(3), decomp%yen(3))
       else
          call alloc_loc_cplx(var, decomp%ysz(1), decomp%ysz(2), decomp%ysz(3), var2d)
       end if
    endif
    
  end subroutine alloc_y_complex

  ! Z-pencil real arrays
  module subroutine alloc_z_real(var, var2d, opt_decomp, opt_global, win)

    implicit none

    real(mytype), dimension(:,:,:), pointer :: var
    real(mytype), dimension(:,:), pointer, optional :: var2d
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global
    integer, intent(out), optional :: win

    TYPE(DECOMP_INFO) :: decomp
    logical :: global

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    ! Safety check
#ifdef DEBUG
    if (d2d_intranode.and.(global.or..not.present(win))) &
       call decomp_2d_abort(__FILE__, __LINE__, -1, "Incorrect alloc arguments")
#endif

    if (d2d_intranode) then
       ! MPI3 shared memory, the caller should free the MPI window
       call alloc_mpi3_real(var, decomp%zsz(1), decomp%zsz(2), decomp%zsz(3), &
                            var2d, decomp%zsz_loc(1), decomp%zsz_loc(2), decomp%zsz_loc(3), win)
    else
       ! No MPI3 shared memory
       if (global) then
          call alloc_glob_real(var, decomp%zst(1), decomp%zen(1), &
                                    decomp%zst(2), decomp%zen(2), &
                                    decomp%zst(3), decomp%zen(3))
       else
          call alloc_loc_real(var, decomp%zsz(1), decomp%zsz(2), decomp%zsz(3), var2d)
       end if
    endif
    
    return
  end subroutine alloc_z_real

  ! Z-pencil complex arrays
  module subroutine alloc_z_complex(var, var2d, opt_decomp, opt_global, win)

    implicit none

    complex(mytype), dimension(:,:,:), pointer :: var
    complex(mytype), dimension(:,:), pointer, optional :: var2d
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global
    integer, intent(out), optional :: win

    TYPE(DECOMP_INFO) :: decomp
    logical :: global

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    ! Safety check
#ifdef DEBUG
    if (d2d_intranode.and.(global.or..not.present(win))) &
       call decomp_2d_abort(__FILE__, __LINE__, -1, "Incorrect alloc arguments")
#endif

    if (d2d_intranode) then
       ! MPI3 shared memory, the caller should free the MPI window
       call alloc_mpi3_cplx(var, decomp%zsz(1), decomp%zsz(2), decomp%zsz(3), &
                            var2d, decomp%zsz_loc(1), decomp%zsz_loc(2), decomp%zsz_loc(3), win)
    else
       ! No MPI3 shared memory
       if (global) then
          call alloc_glob_cplx(var, decomp%zst(1), decomp%zen(1), &
                                    decomp%zst(2), decomp%zen(2), &
                                    decomp%zst(3), decomp%zen(3))
       else
          call alloc_loc_cplx(var, decomp%zsz(1), decomp%zsz(2), decomp%zsz(3), var2d)
       end if
    endif

  end subroutine alloc_z_complex

end submodule smod_alloc
