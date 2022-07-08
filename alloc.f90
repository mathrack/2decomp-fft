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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Utility routine to help allocate 3D arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! X-pencil real arrays
  subroutine alloc_x_real(var, var2d, opt_decomp, opt_global, win)

    implicit none

    real(mytype), dimension(:,:,:), pointer :: var
    real(mytype), dimension(:,:), pointer, optional :: var2d
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global
    integer, intent(out), optional :: win

    TYPE(DECOMP_INFO) :: decomp
    type(c_ptr) :: baseptr
    integer(kind=MPI_ADDRESS_KIND) :: winsize, tmpsize, tmpdispunit
    logical :: global
    integer :: alloc_stat, errorcode, ierror

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

    if (DECOMP_2D_LOCALCOMM == MPI_COMM_NULL) then
       ! No MPI3 shared memory
       if (global) then
          allocate(var(decomp%xst(1):decomp%xen(1), &
               decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)), &
               stat=alloc_stat)
       else
          allocate(var(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)), &
               stat=alloc_stat)
       end if
       if (present(var2d)) nullify(var2d)
    else
       ! MPI3 shared memory, the caller should free the MPI window
       !
       ! Safety check
       if (.not.present(win) .or. .not.present(var2d) .or. global) then
          call decomp_2d_abort(__FILE__, __LINE__, 0, "Allocation is not possible")
       endif
       alloc_stat = 0
       !
       ! Size of the memory to allocate on each CPU
       winsize = decomp%xsz_loc(1) * decomp%xsz_loc(2) * decomp%xsz_loc(3)
       !
       ! Create a window and allow direct memory access inside DECOMP_2D_LOCALCOMM
       call MPI_WIN_ALLOCATE_SHARED(winsize * mytype_bytes, &
                                    mytype_bytes, &
                                    MPI_INFO_NULL, &
                                    DECOMP_2D_LOCALCOMM, &
                                    baseptr, &
                                    win, &
                                    ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_ALLOCATE_SHARED")
       !
       ! Get the CPU-level 2D memory
       if (present(var2d)) call C_F_POINTER(baseptr, var2d, (/decomp%xsz_loc(1),decomp%xsz_loc(2)*decomp%xsz_loc(3)/))
       !
       ! Get the node-level 3D memory
       call MPI_WIN_SHARED_QUERY(win, 0, tmpsize, tmpdispunit, baseptr, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_SHARED_QUERY")
       call C_F_POINTER(baseptr, var, (/decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)/))
    endif
    
    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

    return
  end subroutine alloc_x_real

  ! X-pencil complex arrays
  subroutine alloc_x_complex(var, var2d, opt_decomp, opt_global, win)

    implicit none

    complex(mytype), dimension(:,:,:), pointer :: var
    complex(mytype), dimension(:,:), pointer, optional :: var2d
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global
    integer, intent(out), optional :: win

    TYPE(DECOMP_INFO) :: decomp
    TYPE(C_PTR) :: baseptr ! No need to keep this
    integer(kind=MPI_ADDRESS_KIND) :: winsize, tmpsize, tmpdispunit ! No need to keep this
    logical :: global
    integer :: alloc_stat, errorcode, ierror

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

    if (DECOMP_2D_LOCALCOMM == MPI_COMM_NULL) then
       ! No MPI3 shared memory
       if (global) then
          allocate(var(decomp%xst(1):decomp%xen(1), &
               decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)), &
               stat=alloc_stat)
       else
          allocate(var(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)), &
               stat=alloc_stat)
       end if
       if (present(var2d)) nullify(var2d)
    else
       ! MPI3 shared memory, the caller should free the MPI window
       !
       ! Safety check
       if (.not.present(win) .or. .not.present(var2d) .or. global) then
          call decomp_2d_abort(__FILE__, __LINE__, 0, "Allocation is not possible")
       endif
       alloc_stat = 0
       !
       ! Size of the memory available on each CPU
       winsize = decomp%xsz_loc(1) * decomp%xsz_loc(2) * decomp%xsz_loc(3)
       !
       ! Create a window and allow direct memory access inside DECOMP_2D_LOCALCOMM
       call MPI_WIN_ALLOCATE_SHARED(winsize * mytype_bytes * 2, &
                                    mytype_bytes * 2, &
                                    MPI_INFO_NULL, &
                                    DECOMP_2D_LOCALCOMM, &
                                    baseptr, &
                                    win, &
                                    ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_ALLOCATE_SHARED")
       !
       ! Get the CPU-level 2D memory
       if (present(var2d)) call C_F_POINTER(baseptr, var2d, (/decomp%xsz_loc(1),decomp%xsz_loc(2)*decomp%xsz_loc(3)/))
       !
       ! Get the node-level 3D memory
       call MPI_WIN_SHARED_QUERY(win, 0, tmpsize, tmpdispunit, baseptr, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_SHARED_QUERY")
       call C_F_POINTER(baseptr, var, (/decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)/))
    endif
    
    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

    return
  end subroutine alloc_x_complex

  ! Y-pencil real arrays
  subroutine alloc_y_real(var, var2d, opt_decomp, opt_global, win)

    implicit none

    real(mytype), dimension(:,:,:), pointer :: var
    real(mytype), dimension(:,:), pointer, optional :: var2d
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global
    integer, intent(out), optional :: win

    TYPE(DECOMP_INFO) :: decomp
    type(c_ptr) :: baseptr
    integer(kind=MPI_ADDRESS_KIND) :: winsize, tmpsize, tmpdispunit
    logical :: global
    integer :: alloc_stat, errorcode, ierror

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

    if (DECOMP_2D_LOCALCOMM == MPI_COMM_NULL) then
       ! No MPI3 shared memory
       if (global) then
          allocate(var(decomp%yst(1):decomp%yen(1), &
               decomp%yst(2):decomp%yen(2), decomp%yst(3):decomp%yen(3)), &
               stat=alloc_stat)
       else
          allocate(var(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)), &
               stat=alloc_stat)
       end if
    else
       ! MPI3 shared memory, the caller should free the MPI window
       !
       ! Safety check
       if (.not.present(win) .or. .not.present(var2d) .or. global) then
          call decomp_2d_abort(__FILE__, __LINE__, 0, "Allocation is not possible")
       endif
       alloc_stat = 0
       !
       ! Size of the memory to allocate on each CPU
       winsize = decomp%ysz_loc(1) * decomp%ysz_loc(2) * decomp%ysz_loc(3)
       !
       ! Create a window and allow direct memory access inside DECOMP_2D_LOCALCOMM
       call MPI_WIN_ALLOCATE_SHARED(winsize * mytype_bytes, &
                                    mytype_bytes, &
                                    MPI_INFO_NULL, &
                                    DECOMP_2D_LOCALCOMM, &
                                    baseptr, &
                                    win, &
                                    ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_ALLOCATE_SHARED")
       !
       ! Get the CPU-level 2D memory
       if (present(var2d)) call C_F_POINTER(baseptr, var2d, (/decomp%ysz_loc(1),decomp%ysz_loc(2)*decomp%ysz_loc(3)/))
       !
       ! Get the node-level 3D memory
       call MPI_WIN_SHARED_QUERY(win, 0, tmpsize, tmpdispunit, baseptr, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_SHARED_QUERY")
       call C_F_POINTER(baseptr, var, (/decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)/))
    endif
    
    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

    return
  end subroutine alloc_y_real

  ! Y-pencil complex arrays
  subroutine alloc_y_complex(var, var2d, opt_decomp, opt_global, win)

    implicit none

    complex(mytype), dimension(:,:,:), pointer :: var
    complex(mytype), dimension(:,:), pointer, optional :: var2d
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global
    integer, intent(out), optional :: win

    TYPE(DECOMP_INFO) :: decomp
    type(c_ptr) :: baseptr
    integer(kind=MPI_ADDRESS_KIND) :: winsize, tmpsize, tmpdispunit
    logical :: global
    integer :: alloc_stat, errorcode, ierror

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

    if (DECOMP_2D_LOCALCOMM == MPI_COMM_NULL) then
       ! No MPI3 shared memory
       if (global) then
          allocate(var(decomp%yst(1):decomp%yen(1), &
               decomp%yst(2):decomp%yen(2), decomp%yst(3):decomp%yen(3)), &
               stat=alloc_stat)
       else
          allocate(var(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)), &
               stat=alloc_stat)
       end if
    else
       ! MPI3 shared memory, the caller should free the MPI window
       !
       ! Safety check
       if (.not.present(win) .or. .not.present(var2d) .or. global) then
          call decomp_2d_abort(__FILE__, __LINE__, 0, "Allocation is not possible")
       endif
       alloc_stat = 0
       !
       ! Size of the memory to allocate on each CPU
       winsize = decomp%ysz_loc(1) * decomp%ysz_loc(2) * decomp%ysz_loc(3)
       !
       ! Create a window and allow direct memory access inside DECOMP_2D_LOCALCOMM
       call MPI_WIN_ALLOCATE_SHARED(winsize * mytype_bytes * 2, &
                                    mytype_bytes * 2, &
                                    MPI_INFO_NULL, &
                                    DECOMP_2D_LOCALCOMM, &
                                    baseptr, &
                                    win, &
                                    ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_ALLOCATE_SHARED")
       !
       ! Get the CPU-level 2D memory
       if (present(var2d)) call C_F_POINTER(baseptr, var2d, (/decomp%ysz_loc(1),decomp%ysz_loc(2)*decomp%ysz_loc(3)/))
       !
       ! Get the node-level 3D memory
       call MPI_WIN_SHARED_QUERY(win, 0, tmpsize, tmpdispunit, baseptr, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_SHARED_QUERY")
       call C_F_POINTER(baseptr, var, (/decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)/))
    endif
    
    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

    return
  end subroutine alloc_y_complex

  ! Z-pencil real arrays
  subroutine alloc_z_real(var, var2d, opt_decomp, opt_global, win)

    implicit none

    real(mytype), dimension(:,:,:), pointer :: var
    real(mytype), dimension(:,:), pointer, optional :: var2d
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global
    integer, intent(out), optional :: win

    TYPE(DECOMP_INFO) :: decomp
    type(c_ptr) :: baseptr
    integer(kind=MPI_ADDRESS_KIND) :: winsize, tmpsize, tmpdispunit
    logical :: global
    integer :: alloc_stat, errorcode, ierror

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

    if (DECOMP_2D_LOCALCOMM == MPI_COMM_NULL) then
       ! No MPI3 shared memory
       if (global) then
          allocate(var(decomp%zst(1):decomp%zen(1), &
               decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)), &
               stat=alloc_stat)
       else
          allocate(var(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)), &
               stat=alloc_stat)
       end if
    else
       ! MPI3 shared memory, the caller should free the MPI window
       !
       ! Safety check
       if (.not.present(win) .or. .not.present(var2d) .or. global) then
          call decomp_2d_abort(__FILE__, __LINE__, 0, "Allocation is not possible")
       endif
       alloc_stat = 0
       !
       ! Size of the memory to allocate on each CPU
       winsize = decomp%zsz_loc(1) * decomp%zsz_loc(2) * decomp%zsz_loc(3)
       !
       ! Create a window and allow direct memory access inside DECOMP_2D_LOCALCOMM
       call MPI_WIN_ALLOCATE_SHARED(winsize * mytype_bytes, &
                                    mytype_bytes, &
                                    MPI_INFO_NULL, &
                                    DECOMP_2D_LOCALCOMM, &
                                    baseptr, &
                                    win, &
                                    ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_ALLOCATE_SHARED")
       !
       ! Get the CPU-level 2D memory
       if (present(var2d)) call C_F_POINTER(baseptr, var2d, (/decomp%zsz_loc(1),decomp%zsz_loc(2)*decomp%zsz_loc(3)/))
       !
       ! Get the node-level 3D memory
       call MPI_WIN_SHARED_QUERY(win, 0, tmpsize, tmpdispunit, baseptr, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_SHARED_QUERY")
       call C_F_POINTER(baseptr, var, (/decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)/))
    endif
    
    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

    return
  end subroutine alloc_z_real

  ! Z-pencil complex arrays
  subroutine alloc_z_complex(var, var2d, opt_decomp, opt_global, win)

    implicit none

    complex(mytype), dimension(:,:,:), pointer :: var
    complex(mytype), dimension(:,:), pointer, optional :: var2d
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global
    integer, intent(out), optional :: win

    TYPE(DECOMP_INFO) :: decomp
    type(c_ptr) :: baseptr
    integer(kind=MPI_ADDRESS_KIND) :: winsize, tmpsize, tmpdispunit
    logical :: global
    integer :: alloc_stat, errorcode, ierror

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

    if (DECOMP_2D_LOCALCOMM == MPI_COMM_NULL) then
       ! No MPI3 shared memory
       if (global) then
          allocate(var(decomp%zst(1):decomp%zen(1), &
               decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)), &
               stat=alloc_stat)
       else
          allocate(var(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)), &
               stat=alloc_stat)
       end if
    else
       ! MPI3 shared memory, the caller should free the MPI window
       !
       ! Safety check
       if (.not.present(win) .or. .not.present(var2d) .or. global) then
          call decomp_2d_abort(__FILE__, __LINE__, 0, "Allocation is not possible")
       endif
       alloc_stat = 0
       !
       ! Size of the memory to allocate on each CPU
       winsize = decomp%zsz_loc(1) * decomp%zsz_loc(2) * decomp%zsz_loc(3)
       !
       ! Create a window and allow direct memory access inside DECOMP_2D_LOCALCOMM
       call MPI_WIN_ALLOCATE_SHARED(winsize * mytype_bytes * 2, &
                                    mytype_bytes * 2, &
                                    MPI_INFO_NULL, &
                                    DECOMP_2D_LOCALCOMM, &
                                    baseptr, &
                                    win, &
                                    ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_ALLOCATE_SHARED")
       !
       ! Get the CPU-level 2D memory
       if (present(var2d)) call C_F_POINTER(baseptr, var2d, (/decomp%zsz_loc(1),decomp%zsz_loc(2)*decomp%zsz_loc(3)/))
       !
       ! Get the node-level 3D memory
       call MPI_WIN_SHARED_QUERY(win, 0, tmpsize, tmpdispunit, baseptr, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_SHARED_QUERY")
       call C_F_POINTER(baseptr, var, (/decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)/))
    endif
    
    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

    return
  end subroutine alloc_z_complex
