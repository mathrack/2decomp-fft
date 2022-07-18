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

! This is the main 2D pencil decomposition module

module decomp_2d

  use MPI
  use, intrinsic :: iso_fortran_env, only : real32, real64
#if defined(_GPU)
  use cudafor
#if defined(_NCCL)
  use nccl
#endif
#endif

  implicit none

  private        ! Make everything private unless declared public

#ifdef DOUBLE_PREC
  integer, parameter, public :: mytype = KIND(0._real64)
  integer, parameter, public :: real_type = MPI_DOUBLE_PRECISION
  integer, parameter, public :: real2_type = MPI_2DOUBLE_PRECISION
  integer, parameter, public :: complex_type = MPI_DOUBLE_COMPLEX
#ifdef SAVE_SINGLE
  integer, parameter, public :: mytype_single = KIND(0._real32)
  integer, parameter, public :: real_type_single = MPI_REAL
#else
  integer, parameter, public :: mytype_single = KIND(0._real64)
  integer, parameter, public :: real_type_single = MPI_DOUBLE_PRECISION
#endif
#else
  integer, parameter, public :: mytype = KIND(0._real32)
  integer, parameter, public :: real_type = MPI_REAL
  integer, parameter, public :: real2_type = MPI_2REAL
  integer, parameter, public :: complex_type = MPI_COMPLEX
  integer, parameter, public :: mytype_single = KIND(0._real32)
  integer, parameter, public :: real_type_single = MPI_REAL
#endif

  integer, save, public :: mytype_bytes

  ! some key global variables
  integer, save, public :: nx_global, ny_global, nz_global  ! global size

  integer, save, public :: nrank  ! local MPI rank 
  integer, save, public :: nproc  ! total number of processors
  integer, save, public :: nrank_loc ! intranode MPI rank
  integer, save, public :: nproc_loc ! intranode number of processors

  ! parameters for 2D Cartesian topology 
  integer, save, dimension(2) :: dims, coord, dims_loc, coord_loc
  logical, save, dimension(2) :: periodic
  integer, save, public :: DECOMP_2D_COMM, DECOMP_2D_LOCALCOMM
  integer, save, public :: DECOMP_2D_COMM_CART_X, &
       DECOMP_2D_COMM_CART_Y, DECOMP_2D_COMM_CART_Z 
  integer, save :: DECOMP_2D_COMM_ROW, DECOMP_2D_COMM_COL

  ! define neighboring blocks (to be used in halo-cell support)
  !  first dimension 1=X-pencil, 2=Y-pencil, 3=Z-pencil
  ! second dimension 1=east, 2=west, 3=north, 4=south, 5=top, 6=bottom 
  integer, save, dimension(3,6) :: neighbour 

  ! flags for periodic condition in three dimensions
  logical, save :: periodic_x, periodic_y, periodic_z

  ! flag for intranode communicator
  logical, save, public :: d2d_intranode

#if defined(_GPU)
#if defined(_NCCL)
  integer, save :: row_rank, col_rank
#endif
#endif

  ! derived type to store generic data
  type, public :: decomp_data
     ! flag to determine if real or complex data
     logical :: is_cplx = .true.
     ! associated decomp_info object and direction
     type(decomp_info), pointer :: decomp => null()
     integer :: idir = 0 ! 1-2-3 for x-y-z decomp_info data
     ! 3D array exposing the memory
     real(mytype), dimension(:,:,:), pointer :: var => null()
     complex(mytype), dimension(:,:,:), pointer :: cvar => null()
     ! 2D array exposing the memory
     real(mytype), dimension(:,:), pointer :: var2d => null()
     complex(mytype), dimension(:,:), pointer :: cvar2d => null()
     ! Associated window if MPI3 shared memory
     logical :: shm = .false.
     integer :: win = MPI_WIN_NULL
     contains
        procedure :: init => decomp_data_init
        procedure :: copy => decomp_data_init_copy
        procedure :: fin => decomp_data_fin
  end type decomp_data

  ! derived type to store decomposition info for a given global data size
  TYPE, public :: DECOMP_INFO
     ! staring/ending index and size of data held by current processor
     integer, dimension(3) :: xst, xen, xsz  ! x-pencil
     integer, dimension(3) :: yst, yen, ysz  ! y-pencil
     integer, dimension(3) :: zst, zen, zsz  ! z-pencil

     ! local size and index in case of MPI3 shared memory
     integer, dimension(3) :: xst_loc, xen_loc, xsz_loc
     integer, dimension(3) :: yst_loc, yen_loc, ysz_loc
     integer, dimension(3) :: zst_loc, zen_loc, zsz_loc

     ! in addition to local information, processors also need to know 
     ! some global information for global communications to work 

     ! how each dimension is distributed along pencils
     integer, allocatable, dimension(:) :: &
          x1dist, y1dist, y2dist, z2dist

     ! send/receive buffer counts and displacements for MPI_ALLTOALLV
     integer, allocatable, dimension(:) :: &
          x1cnts, y1cnts, y2cnts, z2cnts
     integer, allocatable, dimension(:) :: &
          x1disp, y1disp, y2disp, z2disp

     ! buffer counts for MPI_ALLTOALL: either for evenly distributed data
     ! or for padded-alltoall
     integer :: x1count, y1count, y2count, z2count

     ! evenly distributed data
     logical :: even

     ! node-level map for intra-node split and merge
     integer, allocatable, dimension(:,:,:) :: intramap_split, intramap_merge

  END TYPE DECOMP_INFO

  ! main (default) decomposition information for global size nx*ny*nz
  TYPE(DECOMP_INFO), save, target, public :: decomp_main
  TYPE(DECOMP_INFO), save, target, public :: phG,ph1,ph2,ph3,ph4

  ! staring/ending index and size of data held by current processor
  ! duplicate 'decomp_main', needed by apps to define data structure 
  integer, save, dimension(3), public :: xstart, xend, xsize  ! x-pencil
  integer, save, dimension(3), public :: ystart, yend, ysize  ! y-pencil
  integer, save, dimension(3), public :: zstart, zend, zsize  ! z-pencil

  ! These are the buffers used by MPI_ALLTOALL(V) calls
  integer, save :: decomp_buf_size = 0
  integer, save :: work1_r_win, work2_r_win, work1_c_win, work2_c_win
  real(mytype), pointer, dimension(:) :: work1_r, work2_r
  complex(mytype), pointer, dimension(:) :: work1_c, work2_c

#if defined(_GPU)
  integer, save :: work1_r_d_win, work2_r_d_win, work1_c_d_win, work2_c_d_win
  real(mytype), pointer, dimension(:), device :: work1_r_d, work2_r_d
  complex(mytype), pointer, dimension(:), device :: work1_c_d, work2_c_d

#if defined(_NCCL)
  integer col_comm_size, row_comm_size
  type(ncclResult) :: nccl_stat
  integer, allocatable, dimension(:) :: local_to_global_col, local_to_global_row
  type(ncclUniqueId) :: nccl_uid_2decomp
  type(ncclComm) :: nccl_comm_2decomp
  integer cuda_stat
  integer(kind=cuda_stream_kind) :: cuda_stream_2decomp
#endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! To define smaller arrays using every several mesh points
  integer, save, dimension(3), public :: xszS,yszS,zszS,xstS,ystS,zstS,xenS,yenS,zenS
  integer, save, dimension(3), public :: xszV,yszV,zszV,xstV,ystV,zstV,xenV,yenV,zenV
  integer, save, dimension(3), public :: xszP,yszP,zszP,xstP,ystP,zstP,xenP,yenP,zenP
  logical, save :: coarse_mesh_starts_from_1
  integer, save :: iskipS, jskipS, kskipS
  integer, save :: iskipV, jskipV, kskipV
  integer, save :: iskipP, jskipP, kskipP

  ! public user routines
  public :: decomp_2d_init, decomp_2d_finalize, &
       decomp_2d_win_fence, decomp_2d_win_free, &
       transpose_x_to_y, transpose_y_to_z, &
       transpose_z_to_y, transpose_y_to_x, &
       decomp_info_init, decomp_info_finalize, &
       decomp_info_print, &
       init_coarser_mesh_statS,fine_to_coarseS,&
       init_coarser_mesh_statV,fine_to_coarseV,&
       init_coarser_mesh_statP,fine_to_coarseP,&
       alloc_x, alloc_y, alloc_z, &
       update_halo, decomp_2d_abort, &
       decomp_2d_warning, get_decomp_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These are routines to perform global data transpositions
  ! 
  !   Four combinations are available, enough to cover all situations
  !    - transpose_x_to_y (X-pencil --> Y-pencil)
  !    - transpose_y_to_z (Y-pencil --> Z-pencil)
  !    - transpose_z_to_y (Z-pencil --> Y-pencil)
  !    - transpose_y_to_x (Y-pencil --> X-pencil)
  !
  !   Generic interface provided here to support multiple data types
  !    - real and complex types supported through generic interface
  !    - single/double precision supported through pre-processing
  !       * see 'mytype' variable at the beginning
  !    - an optional argument can be supplied to transpose data whose 
  !      global size is not the default nx*ny*nz 
  !       * as the case in fft r2c/c2r interface 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface transpose_x_to_y
     module procedure transpose_x_to_y_data
  end interface transpose_x_to_y

  interface transpose_y_to_z
     module procedure transpose_y_to_z_data
  end interface transpose_y_to_z

  interface transpose_z_to_y
     module procedure transpose_z_to_y_data
  end interface transpose_z_to_y

  interface transpose_y_to_x
     module procedure transpose_y_to_x_data
  end interface transpose_y_to_x

  interface update_halo
     module procedure update_halo_real
     module procedure update_halo_complex
  end interface update_halo

  ! Submodule alloc

  interface
     module subroutine alloc_x_real(var, var2d, opt_decomp, opt_global, win)
        real(mytype), dimension(:,:,:), pointer :: var
        real(mytype), dimension(:,:), pointer, optional :: var2d
        TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
        logical, intent(IN), optional :: opt_global
        integer, intent(out), optional :: win
     end subroutine alloc_x_real
     module subroutine alloc_x_complex(var, var2d, opt_decomp, opt_global, win)
        complex(mytype), dimension(:,:,:), pointer :: var
        complex(mytype), dimension(:,:), pointer, optional :: var2d
        TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
        logical, intent(IN), optional :: opt_global
        integer, intent(out), optional :: win
     end subroutine alloc_x_complex
  end interface
  interface alloc_x
     module procedure alloc_x_real, alloc_x_complex
  end interface alloc_x

  interface
     module subroutine alloc_y_real(var, var2d, opt_decomp, opt_global, win)
        real(mytype), dimension(:,:,:), pointer :: var
        real(mytype), dimension(:,:), pointer, optional :: var2d
        TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
        logical, intent(IN), optional :: opt_global
        integer, intent(out), optional :: win
     end subroutine alloc_y_real
     module subroutine alloc_y_complex(var, var2d, opt_decomp, opt_global, win)
        complex(mytype), dimension(:,:,:), pointer :: var
        complex(mytype), dimension(:,:), pointer, optional :: var2d
        TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
        logical, intent(IN), optional :: opt_global
        integer, intent(out), optional :: win
     end subroutine alloc_y_complex
  end interface
  interface alloc_y
     module procedure alloc_y_real, alloc_y_complex
  end interface alloc_y

  interface
     module subroutine alloc_z_real(var, var2d, opt_decomp, opt_global, win)
        real(mytype), dimension(:,:,:), pointer :: var
        real(mytype), dimension(:,:), pointer, optional :: var2d
        TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
        logical, intent(IN), optional :: opt_global
        integer, intent(out), optional :: win
     end subroutine alloc_z_real
     module subroutine alloc_z_complex(var, var2d, opt_decomp, opt_global, win)
        complex(mytype), dimension(:,:,:), pointer :: var 
        complex(mytype), dimension(:,:), pointer, optional :: var2d
        TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
        logical, intent(IN), optional :: opt_global
        integer, intent(out), optional :: win
     end subroutine alloc_z_complex
  end interface
  interface alloc_z
     module procedure alloc_z_real, alloc_z_complex
  end interface alloc_z

  interface decomp_2d_abort
     module procedure decomp_2d_abort_basic
     module procedure decomp_2d_abort_file_line
  end interface decomp_2d_abort

  interface decomp_2d_warning
     module procedure decomp_2d_warning_basic
     module procedure decomp_2d_warning_file_line
  end interface decomp_2d_warning

  interface

     ! Submodule buffer

     module subroutine decomp_buffer_alloc(buf_size)
        integer, intent(inout) :: buf_size
     end subroutine decomp_buffer_alloc

     module subroutine decomp_buffer_init
     end subroutine decomp_buffer_init

     module subroutine decomp_buffer_free
     end subroutine decomp_buffer_free

     ! Submodule coarse

     module subroutine init_coarser_mesh_statS(i_skip,j_skip,k_skip,from1)
        integer, intent(IN) :: i_skip,j_skip,k_skip
        logical, intent(IN) :: from1
     end subroutine init_coarser_mesh_statS

     module subroutine init_coarser_mesh_statV(i_skip,j_skip,k_skip,from1)
        integer, intent(IN) :: i_skip,j_skip,k_skip
        logical, intent(IN) :: from1
     end subroutine init_coarser_mesh_statV

     module subroutine init_coarser_mesh_statP(i_skip,j_skip,k_skip,from1)
        integer, intent(IN) :: i_skip,j_skip,k_skip
        logical, intent(IN) :: from1
     end subroutine init_coarser_mesh_statP

     module subroutine fine_to_coarseS(ipencil,var_fine,var_coarse)
        integer, intent(IN) :: ipencil
        real(mytype), dimension(:,:,:) :: var_fine
        real(mytype), dimension(:,:,:) :: var_coarse
     end subroutine fine_to_coarseS

     module subroutine fine_to_coarseV(ipencil,var_fine,var_coarse)
        integer, intent(IN) :: ipencil
        real(mytype), dimension(:,:,:) :: var_fine
        real(mytype), dimension(:,:,:) :: var_coarse
     end subroutine fine_to_coarseV

     module subroutine fine_to_coarseP(ipencil,var_fine,var_coarse)
        integer, intent(IN) :: ipencil
        real(mytype), dimension(:,:,:) :: var_fine
        real(mytype), dimension(:,:,:) :: var_coarse
     end subroutine fine_to_coarseP

     ! Submodule decomp_data

     module subroutine decomp_data_init(self, is_cplx, idir, decomp, rwk, cwk)
        class(decomp_data), intent(out) :: self
        logical, intent(in) :: is_cplx
        integer, intent(in) :: idir
        type(decomp_info), intent(in), target, optional :: decomp
        real(mytype), dimension(:,:,:), target, optional :: rwk
        complex(mytype), dimension(:,:,:), target, optional :: cwk
     end subroutine decomp_data_init

     module subroutine decomp_data_init_copy(self, dat)
        class(decomp_data), intent(out) :: self
        type(decomp_data), intent(in) :: dat
     end subroutine decomp_data_init_copy

     module subroutine decomp_data_fin(self)
        class(decomp_data), intent(inout) :: self
     end subroutine decomp_data_fin

     ! Submodule decomp_info

     module subroutine decomp_info_init_reshapeyz(decomp)
        type(decomp_info), intent(inout) :: decomp
     end subroutine decomp_info_init_reshapeyz

     module subroutine get_decomp_info(decomp)
        type(decomp_info), intent(out) :: decomp
     end subroutine get_decomp_info

     module subroutine decomp_info_init(nx,ny,nz,decomp)
        integer, intent(IN) :: nx,ny,nz
        type(decomp_info), intent(inout) :: decomp
     end subroutine decomp_info_init

     module subroutine decomp_info_init_local(decomp)
        type(decomp_info), intent(inout) :: decomp
     end subroutine decomp_info_init_local

     module subroutine decomp_info_finalize(decomp)
        type(decomp_info), intent(inout) :: decomp
     end subroutine decomp_info_finalize

     ! Submodule log

     module function d2d_get_iounit() result(iounit)
        integer :: iounit
     end function d2d_get_iounit

     module subroutine d2d_listing(given_io_unit)
        integer, intent(in), optional :: given_io_unit
     end subroutine d2d_listing

     module subroutine decomp_info_print(d2d, io_unit, d2dname)
        type(decomp_info), intent(in) :: d2d
        integer, intent(in) :: io_unit
        character(len=*), intent(in) :: d2dname
     end subroutine decomp_info_print

     ! Submodule wrapper

     module subroutine decomp_2d_win_transpose_start_reading(src_win)
        integer, intent(in) :: src_win
     end subroutine decomp_2d_win_transpose_start_reading

     module subroutine decomp_2d_win_transpose_stop_reading(src_win)                                 
        integer, intent(in) :: src_win                                                               
     end subroutine decomp_2d_win_transpose_stop_reading

     module subroutine decomp_2d_win_transpose_start_writing(dst_win)
        integer, intent(in) :: dst_win
     end subroutine decomp_2d_win_transpose_start_writing

     module subroutine decomp_2d_win_transpose_stop_writing(dst_win)
        integer, intent(in) :: dst_win
     end subroutine decomp_2d_win_transpose_stop_writing

     module subroutine decomp_2d_win_fence(win, assert)
        integer, intent(in) :: win
        integer, intent(in), optional :: assert
     end subroutine decomp_2d_win_fence

     module subroutine decomp_2d_win_free(win)    
        integer, intent(inout) :: win
     end subroutine decomp_2d_win_free

  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine to be called by applications to initialise this library
  !   INPUT:
  !     nx, ny, nz   - global data dimension
  !     p_row, p_col - 2D processor grid
  !   OUTPUT:
  !     all internal data structures initialised properly
  !     library ready to use
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_init(nx,ny,nz,p_row,p_col,periodic_bc,glob_comm,local_comm)

    use iso_fortran_env, only : output_unit

    implicit none

    integer, intent(IN) :: nx,ny,nz
    integer, intent(INOUT) :: p_row,p_col
    logical, dimension(3), intent(IN), optional :: periodic_bc
    integer, intent(in), optional :: glob_comm, local_comm

    integer :: errorcode, ierror, row, col

    if (present(glob_comm)) then
       DECOMP_2D_COMM = glob_comm
    else
       DECOMP_2D_COMM = MPI_COMM_WORLD
    endif

    call decomp_buffer_init()

    if (DECOMP_2D_COMM /= MPI_COMM_WORLD .and. present(local_comm)) then
       ! MPI3 shared memory
       d2d_intranode = .true.
       DECOMP_2D_LOCALCOMM = local_comm
       ! Only local masters will perform MPI operations
       if (DECOMP_2D_COMM == MPI_COMM_NULL) then
          nrank = -1
          nproc = -1
          ! Intra-node CPU map
          call decomp_2d_map_local()
          ! Initialize the main decomp_info object
          call decomp_info_init(nx, ny, nz, decomp_main)
          ! Receive some decomp_2d stuff from local master
          call decomp_2d_init_local()
          ! Print decomp_2d setup
          call d2d_listing(d2d_get_iounit())
          return
       endif
    else
       ! No MPI3 shared memory
       d2d_intranode = .false.
       DECOMP_2D_LOCALCOMM = MPI_COMM_NULL
       nrank_loc = 0
       nproc_loc = -1
    endif

    !
    ! Get global rank and comm size and mytype_bytes
    !
    call MPI_COMM_RANK(DECOMP_2D_COMM, nrank, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_RANK")
    call MPI_COMM_SIZE(DECOMP_2D_COMM, nproc, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SIZE")
    ! determine the number of bytes per float number
    ! do not use 'mytype' which is compiler dependent
    ! also possible to use inquire(iolength=...) 
    call MPI_TYPE_SIZE(real_type,mytype_bytes,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_SIZE")

    nx_global = nx
    ny_global = ny
    nz_global = nz

    if (present(periodic_bc)) then
       periodic_x = periodic_bc(1)
       periodic_y = periodic_bc(2)
       periodic_z = periodic_bc(3)
    else
       periodic_x = .false.
       periodic_y = .false.
       periodic_z = .false.
    end if

    if (p_row==0 .and. p_col==0) then
       ! determine the best 2D processor grid
       call best_2d_grid(nproc, row, col)
       p_row = row
       p_col = col
    else
       if (nproc /= p_row*p_col) then
          errorcode = 1
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Invalid 2D processor grid - nproc /= p_row*p_col')
       else
          row = p_row
          col = p_col
       end if
    end if

    ! Create 2D Catersian topology
    ! Note that in order to support periodic B.C. in the halo-cell code,
    ! need to create multiple topology objects: DECOMP_2D_COMM_CART_?,
    ! corresponding to three pencil orientations. They contain almost
    ! identical topological information but allow different combinations
    ! of periodic conditions.
    dims(1) = row
    dims(2) = col
    periodic(1) = periodic_y
    periodic(2) = periodic_z
    call MPI_CART_CREATE(DECOMP_2D_COMM,2,dims,periodic, &
         .false., &  ! do not reorder rank
         DECOMP_2D_COMM_CART_X, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_CREATE")
    periodic(1) = periodic_x
    periodic(2) = periodic_z
    call MPI_CART_CREATE(DECOMP_2D_COMM,2,dims,periodic, &
         .false., DECOMP_2D_COMM_CART_Y, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_CREATE")
    periodic(1) = periodic_x
    periodic(2) = periodic_y
    call MPI_CART_CREATE(DECOMP_2D_COMM,2,dims,periodic, &
         .false., DECOMP_2D_COMM_CART_Z, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_CREATE")

    call MPI_CART_COORDS(DECOMP_2D_COMM_CART_X,nrank,2,coord,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_COORDS")

    ! derive communicators defining sub-groups for ALLTOALL(V)
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.true.,.false./), &
         DECOMP_2D_COMM_COL,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SUB")
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X,(/.false.,.true./), &
         DECOMP_2D_COMM_ROW,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SUB")

    ! MPI3 shared memory : intra-node CPU map
    if (d2d_intranode) call decomp_2d_map_local()

    ! gather information for halo-cell support code
    call init_neighbour

    ! actually generate all 2D decomposition information
    call decomp_info_init(nx,ny,nz,decomp_main)

    ! make a copy of the decomposition information associated with the
    ! default global size in these global variables so applications can
    ! use them to create data structures 
    xstart = decomp_main%xst
    ystart = decomp_main%yst
    zstart = decomp_main%zst
    xend   = decomp_main%xen
    yend   = decomp_main%yen
    zend   = decomp_main%zen
    xsize  = decomp_main%xsz
    ysize  = decomp_main%ysz
    zsize  = decomp_main%zsz

#ifdef EVEN
    if (nrank==0.and.nrank_loc<=0) write(*,*) 'Padded ALLTOALL optimisation on'
#endif 

#if defined(_GPU)
#if defined(_NCCL)
    call MPI_COMM_RANK(DECOMP_2D_COMM_COL,col_rank,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_RANK")
    call MPI_COMM_RANK(DECOMP_2D_COMM_ROW,row_rank,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_RANK")
    call MPI_COMM_SIZE(DECOMP_2D_COMM_COL,col_comm_size,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SIZE")
    call MPI_COMM_SIZE(DECOMP_2D_COMM_ROW,row_comm_size,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SIZE")

    allocate(local_to_global_col(col_comm_size), local_to_global_row(row_comm_size))
    
    local_to_global_col(:) = 0
    local_to_global_row(:) = 0
    local_to_global_col(col_rank+1) = nrank
    local_to_global_row(row_rank+1) = nrank
    
    call mpi_allreduce(MPI_IN_PLACE,local_to_global_col,col_comm_size,MPI_INTEGER,MPI_SUM,DECOMP_2D_COMM_COL,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLREDUCE")
    call mpi_allreduce(MPI_IN_PLACE,local_to_global_row,row_comm_size,MPI_INTEGER,MPI_SUM,DECOMP_2D_COMM_ROW,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLREDUCE")

    if (nrank .eq. 0) then
       nccl_stat = ncclGetUniqueId(nccl_uid_2decomp)
    end if
    call MPI_Bcast(nccl_uid_2decomp, int(sizeof(ncclUniqueId)), MPI_BYTE, 0, DECOMP_2D_COMM, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")

    nccl_stat = ncclCommInitRank(nccl_comm_2decomp, nproc, nccl_uid_2decomp, nrank)
    cuda_stat = cudaStreamCreate(cuda_stream_2decomp)
#endif
#endif

    ! If MPI3 shared memory, local master should broadcast some decomp_2d stuff
    if (d2d_intranode) call decomp_2d_init_local()

    !
    ! Print decomp_2d setup
    !
    call d2d_listing(d2d_get_iounit())

    return
  end subroutine decomp_2d_init

  !
  ! If MPI3 shared memory, local master should broadcast decomp_2d stuff
  !
  subroutine decomp_2d_init_local()

     implicit none

     integer :: ierror

     !
     ! Broadcast
     !

     call MPI_BCAST(xstart, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(xend, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(xsize, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")

     call MPI_BCAST(ystart, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(yend, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(ysize, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")

     call MPI_BCAST(zstart, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(zend, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(zsize, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")

  end subroutine decomp_2d_init_local

  !
  ! MPI3 shared memory : intranode CPU map (1D slab decomposition)
  !
  subroutine decomp_2d_map_local()

     implicit none

     integer :: ierror, TMP_COMM_CART

     call MPI_BCAST(nx_global, 1, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(ny_global, 1, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(nz_global, 1, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")

     call MPI_BCAST(nrank, 1, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(nproc, 1, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(dims, 2, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(coord, 2, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")

     !
     ! Get local rank and comm size and mytype_bytes
     call MPI_COMM_RANK(DECOMP_2D_LOCALCOMM, nrank_loc, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_RANK")
     call MPI_COMM_SIZE(DECOMP_2D_LOCALCOMM, nproc_loc, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SIZE")
     call MPI_BCAST(mytype_bytes, 1, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     !
     ! Build a 1D CPU grid (arrays dims_loc and coord_loc)
     dims_loc(1) = 1
     dims_loc(2) = nproc_loc
     call MPI_CART_CREATE(DECOMP_2D_LOCALCOMM, 2, dims_loc, &
             (/.false.,.false./), & ! no periodicity
             .false., &  ! do not reorder rank
             TMP_COMM_CART, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_CREATE")
     call MPI_CART_COORDS(TMP_COMM_CART, nrank_loc, 2, coord_loc, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_COORDS")
     !
     ! Free the cartesian MPI comm
     call MPI_COMM_FREE(TMP_COMM_CART, ierror)
     if (ierror /= 0) call decomp_2d_warning(__FILE__, __LINE__, ierror, "MPI_COMM_FREE")

  end subroutine decomp_2d_map_local

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine to be called by applications to clean things up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_finalize()
    
    implicit none
 
    integer :: ierror

    call decomp_buffer_free()

    ! Release localcomm if possible
    if (DECOMP_2D_LOCALCOMM /= MPI_COMM_NULL .and. DECOMP_2D_LOCALCOMM /= MPI_COMM_WORLD) then       
       call MPI_COMM_FREE(DECOMP_2D_LOCALCOMM, ierror)                                               
       if (ierror /= 0) call decomp_2d_warning(__FILE__, __LINE__, ierror, "MPI_COMM_FREE")          
    endif 

    ! Nothing more to do if MPI3 shared memory and the rank is not local master
    if (DECOMP_2D_COMM == MPI_COMM_NULL) return

    ierror = 0
    if (DECOMP_2D_COMM /= MPI_COMM_WORLD) call MPI_COMM_FREE(DECOMP_2D_COMM, ierror)
    if (ierror /= 0) call decomp_2d_warning(__FILE__, __LINE__, ierror, "MPI_COMM_FREE")
    call MPI_COMM_FREE(DECOMP_2D_COMM_ROW, ierror)
    if (ierror /= 0) call decomp_2d_warning(__FILE__, __LINE__, ierror, "MPI_COMM_FREE")
    call MPI_COMM_FREE(DECOMP_2D_COMM_COL, ierror)
    if (ierror /= 0) call decomp_2d_warning(__FILE__, __LINE__, ierror, "MPI_COMM_FREE")
    call MPI_COMM_FREE(DECOMP_2D_COMM_CART_X, ierror)
    if (ierror /= 0) call decomp_2d_warning(__FILE__, __LINE__, ierror, "MPI_COMM_FREE")
    call MPI_COMM_FREE(DECOMP_2D_COMM_CART_Y, ierror)
    if (ierror /= 0) call decomp_2d_warning(__FILE__, __LINE__, ierror, "MPI_COMM_FREE")
    call MPI_COMM_FREE(DECOMP_2D_COMM_CART_Z, ierror)
    if (ierror /= 0) call decomp_2d_warning(__FILE__, __LINE__, ierror, "MPI_COMM_FREE")
    call decomp_info_finalize(decomp_main)

  end subroutine decomp_2d_finalize


#ifdef OCC
  ! For non-blocking communication code, progress the comminication stack
  subroutine transpose_test(handle)

    implicit none

    integer :: handle, ierror

    call NBC_TEST(handle,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "NBC_TEST")

    return
  end subroutine transpose_test
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Transposition routines 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "transpose_x_to_y.f90"
#include "transpose_y_to_z.f90"
#include "transpose_z_to_y.f90"
#include "transpose_y_to_x.f90"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Auto-tuning algorithm to select the best 2D processor grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine best_2d_grid(iproc, best_p_row, best_p_col)

    implicit none

    integer, intent(IN) :: iproc
    integer, intent(OUT) :: best_p_row, best_p_col

    integer, allocatable, dimension(:) :: factors
    integer :: nfact, i, col, i_best

    if (nrank==0.and.nrank_loc<=0) write(*,*) 'In auto-tuning mode......'

    i = int(sqrt(real(iproc))) + 10  ! enough space to save all factors 
    allocate(factors(i))
    call findfactor(iproc, factors, nfact)
    if (nrank==0.and.nrank_loc<=0) write(*,*) 'factors: ', (factors(i), i=1,nfact)

    i_best=nfact/2+1
    col=factors(i_best)

    best_p_col = col
    best_p_row=iproc/col
    if (nrank==0.and.nrank_loc<=0) print *,'p_row x p_col', best_p_row, best_p_col
    if (best_p_col==1.and.nrank==0.and.nrank_loc<=0) then
       print *,'WARNING: current 2D DECOMP set-up might not work'
    endif
    
    deallocate(factors)

    return
  end subroutine best_2d_grid

#include "factor.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Halo cell support
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "halo.f90"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Error handling
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_abort_basic(errorcode, msg)

    use iso_fortran_env, only : error_unit

    implicit none

    integer, intent(IN) :: errorcode
    character(len=*), intent(IN) :: msg

    integer :: ierror

    if (nrank==0.and.nrank_loc<=0) then
       write(*,*) '2DECOMP&FFT ERROR - errorcode: ', errorcode
       write(*,*) 'ERROR MESSAGE: ' // msg
       write(error_unit,*) '2DECOMP&FFT ERROR - errorcode: ', errorcode
       write(error_unit,*) 'ERROR MESSAGE: ' // msg
    end if
    call MPI_ABORT(DECOMP_2D_COMM,errorcode,ierror)

  end subroutine decomp_2d_abort_basic

  subroutine decomp_2d_abort_file_line(file, line, errorcode, msg)

    use iso_fortran_env, only : error_unit

    implicit none

    integer, intent(IN) :: errorcode, line
    character(len=*), intent(IN) :: msg, file

    integer :: ierror

    if (nrank==0.and.nrank_loc<=0) then
       write(*,*) '2DECOMP&FFT / X3D ERROR'
       write(*,*) '  errorcode:     ', errorcode
       write(*,*) '  error in file  ' // file
       write(*,*) '           line  ', line
       write(*,*) '  error message: ' // msg
       write(error_unit,*) '2DECOMP&FFT / X3D ERROR'
       write(error_unit,*) '  errorcode:     ', errorcode
       write(error_unit,*) '  error in file  ' // file
       write(error_unit,*) '           line  ', line
       write(error_unit,*) '  error message: ' // msg
    end if
    call MPI_ABORT(DECOMP_2D_COMM,errorcode,ierror)

  end subroutine decomp_2d_abort_file_line

  subroutine decomp_2d_warning_basic(errorcode, msg)

    use iso_fortran_env, only : error_unit

    implicit none

    integer, intent(IN) :: errorcode
    character(len=*), intent(IN) :: msg

    if (nrank==0.and.nrank_loc<=0) then
       write(*,*) '2DECOMP&FFT WARNING - errorcode: ', errorcode
       write(*,*) 'ERROR MESSAGE: ' // msg
       write(error_unit,*) '2DECOMP&FFT WARNING - errorcode: ', errorcode
       write(error_unit,*) 'ERROR MESSAGE: ' // msg
    end if

  end subroutine decomp_2d_warning_basic

  subroutine decomp_2d_warning_file_line(file, line, errorcode, msg)

    use iso_fortran_env, only : error_unit

    implicit none

    integer, intent(IN) :: errorcode, line
    character(len=*), intent(IN) :: msg, file

    if (nrank==0.and.nrank_loc<=0) then
       write(*,*) '2DECOMP&FFT / X3D WARNING'
       write(*,*) '  errorcode:     ', errorcode
       write(*,*) '  error in file  ' // file
       write(*,*) '           line  ', line
       write(*,*) '  error message: ' // msg
       write(error_unit,*) '2DECOMP&FFT / X3D WARNING'
       write(error_unit,*) '  errorcode:     ', errorcode
       write(error_unit,*) '  error in file  ' // file
       write(error_unit,*) '           line  ', line
       write(error_unit,*) '  error message: ' // msg
    end if

  end subroutine decomp_2d_warning_file_line

end module decomp_2d

