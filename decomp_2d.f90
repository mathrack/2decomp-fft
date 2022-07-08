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
  use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
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

#if defined(_GPU)
#if defined(_NCCL)
  integer, save :: row_rank, col_rank
#endif
#endif

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

  END TYPE DECOMP_INFO

  ! main (default) decomposition information for global size nx*ny*nz
  TYPE(DECOMP_INFO), save, public :: decomp_main
  TYPE(DECOMP_INFO), save, public :: phG,ph1,ph2,ph3,ph4

  ! staring/ending index and size of data held by current processor
  ! duplicate 'decomp_main', needed by apps to define data structure 
  integer, save, dimension(3), public :: xstart, xend, xsize  ! x-pencil
  integer, save, dimension(3), public :: ystart, yend, ysize  ! y-pencil
  integer, save, dimension(3), public :: zstart, zend, zsize  ! z-pencil

  ! These are the buffers used by MPI_ALLTOALL(V) calls
  integer, save :: decomp_buf_size = 0
  real(mytype),    allocatable, dimension(:) :: work1_r, work2_r
  complex(mytype), allocatable, dimension(:) :: work1_c, work2_c

#if defined(_GPU)
  real(mytype), allocatable, dimension(:), device :: work1_r_d, work2_r_d
  complex(mytype), allocatable, dimension(:), device :: work1_c_d, work2_c_d

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
       transpose_x_to_y, transpose_y_to_z, &
       transpose_z_to_y, transpose_y_to_x, &
       decomp_info_init, decomp_info_finalize, partition, &
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
     module procedure transpose_x_to_y_real
     module procedure transpose_x_to_y_complex
  end interface transpose_x_to_y

  interface transpose_y_to_z
     module procedure transpose_y_to_z_real
     module procedure transpose_y_to_z_complex
  end interface transpose_y_to_z

  interface transpose_z_to_y
     module procedure transpose_z_to_y_real
     module procedure transpose_z_to_y_complex
  end interface transpose_z_to_y

  interface transpose_y_to_x
     module procedure transpose_y_to_x_real
     module procedure transpose_y_to_x_complex
  end interface transpose_y_to_x

  interface update_halo
     module procedure update_halo_real
     module procedure update_halo_complex
  end interface update_halo

  interface alloc_x
     module procedure alloc_x_real
     module procedure alloc_x_complex
  end interface alloc_x

  interface alloc_y
     module procedure alloc_y_real
     module procedure alloc_y_complex
  end interface alloc_y

  interface alloc_z
     module procedure alloc_z_real
     module procedure alloc_z_complex
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

    if (DECOMP_2D_COMM /= MPI_COMM_WORLD .and. present(local_comm)) then
       ! MPI3 shared memory
       DECOMP_2D_LOCALCOMM = local_comm
       ! Only local masters will perform MPI operations
       if (DECOMP_2D_COMM == MPI_COMM_NULL) then
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
       DECOMP_2D_LOCALCOMM = MPI_COMM_NULL
       nrank_loc = -1
       nproc_loc = -1
    endif

    !
    ! Get global rank and comm size
    !
    call MPI_COMM_RANK(DECOMP_2D_COMM, nrank, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_RANK")
    call MPI_COMM_SIZE(DECOMP_2D_COMM, nproc, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SIZE")

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
    if (DECOMP_2D_LOCALCOMM /= MPI_COMM_NULL) call decomp_2d_map_local()

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

    ! determine the number of bytes per float number
    ! do not use 'mytype' which is compiler dependent
    ! also possible to use inquire(iolength=...) 
    call MPI_TYPE_SIZE(real_type,mytype_bytes,ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_SIZE")

#ifdef EVEN
    if (nrank==0) write(*,*) 'Padded ALLTOALL optimisation on'
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
    if (DECOMP_2D_COMM /= MPI_COMM_WORLD) call decomp_2d_init_local()

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

     call MPI_BCAST(mytype_bytes, 1, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")

     call MPI_BCAST(nx_global, 1, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(ny_global, 1, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(nz_global, 1, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")

     call MPI_BCAST(dims, 2, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(coord, 2, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")

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

     !
     ! Get local rank and comm size
     call MPI_COMM_RANK(DECOMP_2D_LOCALCOMM, nrank_loc, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_RANK")
     call MPI_COMM_SIZE(DECOMP_2D_LOCALCOMM, nproc_loc, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SIZE")
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

  !
  ! Reshape Y and Z pencils so that the leading dimension is the working one
  !
  subroutine decomp_info_init_reshapeyz(decomp)

     implicit none

     type(decomp_info), intent(inout) :: decomp

     call d2d_swap(decomp%yst, 1, 2)
     call d2d_swap(decomp%yen, 1, 2)
     call d2d_swap(decomp%ysz, 1, 2)

     ! reshape z-pencil like x-pencil
     call d2d_swap(decomp%zst, 2, 3)
     call d2d_swap(decomp%zst, 1, 2)
     call d2d_swap(decomp%zen, 2, 3)
     call d2d_swap(decomp%zen, 1, 2)
     call d2d_swap(decomp%zsz, 2, 3)
     call d2d_swap(decomp%zsz, 1, 2)

  end subroutine decomp_info_init_reshapeyz
  !
  ! Swap integers
  !
  subroutine d2d_swap(array, i, j)

     implicit none

     integer, intent(inout) :: array(3)
     integer, intent(in) :: i, j
     integer :: tmp

     tmp = array(i)
     array(i) = array(j)
     array(j) = tmp

  end subroutine d2d_swap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine to be called by applications to clean things up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_finalize
    
    implicit none
 
    integer :: ierror

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

    decomp_buf_size = 0
    deallocate(work1_r, work2_r, work1_c, work2_c)
#if defined(_GPU)
    deallocate(work1_r_d, work2_r_d, work1_c_d, work2_c_d)

#if defined(_NCCL)
    nccl_stat = ncclCommDestroy(nccl_comm_2decomp)
#endif
#endif

    return
  end subroutine decomp_2d_finalize


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Return the default decomposition object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_decomp_info(decomp)

    implicit none

    TYPE(DECOMP_INFO), intent(OUT) :: decomp

    decomp = decomp_main

    return
  end subroutine get_decomp_info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Advanced Interface allowing applications to define globle domain of
  ! any size, distribute it, and then transpose data among pencils.
  !  - generate 2D decomposition details as defined in DECOMP_INFO
  !  - the default global data size is nx*ny*nz
  !  - a different global size nx/2+1,ny,nz is used in FFT r2c/c2r
  !  - multiple global sizes can co-exist in one application, each
  !    using its own DECOMP_INFO object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_info_init(nx,ny,nz,decomp)

    implicit none

    integer, intent(IN) :: nx,ny,nz
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    integer :: buf_size, status, errorcode

    ! In case of MPI3 shared memory and proc is not local master
    if (DECOMP_2D_COMM == MPI_COMM_NULL) then
       call decomp_info_init_local(decomp)
       call decomp_info_init_reshapeyz(decomp)
       call partition(decomp%xsz(1),1,decomp%xsz(2)*decomp%xsz(3), &
            (/1,2,3/), dims_loc, coord_loc, &
            decomp%xst_loc,decomp%xen_loc,decomp%xsz_loc)
       call partition(decomp%ysz(1),1,decomp%ysz(2)*decomp%ysz(3), &
            (/1,2,3/), dims_loc, coord_loc, &
            decomp%yst_loc,decomp%yen_loc,decomp%ysz_loc)
       call partition(decomp%zsz(1),1,decomp%zsz(2)*decomp%zsz(3), &
            (/1,2,3/), dims_loc, coord_loc, &
            decomp%zst_loc,decomp%zen_loc,decomp%zsz_loc)
       return
    endif

    ! verify the global size can actually be distributed as pencils
    if (nx_global<dims(1) .or. ny_global<dims(1) .or. ny_global<dims(2) .or. nz_global<dims(2)) then
       errorcode = 6
       call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
            'Invalid 2D processor grid. ' // &
            'Make sure that min(nx,ny) >= p_row and ' // &
            'min(ny,nz) >= p_col')
    end if

    if (mod(nx,dims(1))==0 .and. mod(ny,dims(1))==0 .and. &
         mod(ny,dims(2))==0 .and. mod(nz,dims(2))==0) then
       decomp%even = .true.
    else
       decomp%even = .false.
    end if

    ! distribute mesh points
    allocate(decomp%x1dist(0:dims(1)-1),decomp%y1dist(0:dims(1)-1), &
         decomp%y2dist(0:dims(2)-1),decomp%z2dist(0:dims(2)-1))
    call get_dist(nx,ny,nz,decomp)

    ! generate partition information - starting/ending index etc.
    call partition(nx, ny, nz, (/ 1,2,3 /), dims, coord, &
         decomp%xst, decomp%xen, decomp%xsz)
    call partition(nx, ny, nz, (/ 2,1,3 /), dims, coord, &
         decomp%yst, decomp%yen, decomp%ysz)
    call partition(nx, ny, nz, (/ 2,3,1 /), dims, coord, &
         decomp%zst, decomp%zen, decomp%zsz)

    ! prepare send/receive buffer displacement and count for ALLTOALL(V)
    allocate(decomp%x1cnts(0:dims(1)-1),decomp%y1cnts(0:dims(1)-1), &
         decomp%y2cnts(0:dims(2)-1),decomp%z2cnts(0:dims(2)-1))
    allocate(decomp%x1disp(0:dims(1)-1),decomp%y1disp(0:dims(1)-1), &
         decomp%y2disp(0:dims(2)-1),decomp%z2disp(0:dims(2)-1))
    call prepare_buffer(decomp)

    ! allocate memory for the MPI_ALLTOALL(V) buffers
    ! define the buffers globally for performance reason

    buf_size = max(decomp%xsz(1)*decomp%xsz(2)*decomp%xsz(3), &
         max(decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3), &
         decomp%zsz(1)*decomp%zsz(2)*decomp%zsz(3)) )
#ifdef EVEN
    ! padded alltoall optimisation may need larger buffer space
    buf_size = max(buf_size, &
         max(decomp%x1count*dims(1),decomp%y2count*dims(2)) ) 
#endif

    ! check if additional memory is required
    ! *** TODO: consider how to share the real/complex buffers 
    if (buf_size > decomp_buf_size) then
       decomp_buf_size = buf_size
#if defined(_GPU)
       if (allocated(work1_r_d)) deallocate(work1_r_d)
       if (allocated(work2_r_d)) deallocate(work2_r_d)
       if (allocated(work1_c_d)) deallocate(work1_c_d)
       if (allocated(work2_c_d)) deallocate(work2_c_d)
       allocate(work1_r_d(buf_size), STAT=status)
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
       allocate(work1_c_d(buf_size), STAT=status)
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
       allocate(work2_r_d(buf_size), STAT=status)
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
       allocate(work2_c_d(buf_size), STAT=status)
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
#endif
       if (allocated(work1_r)) deallocate(work1_r)
       if (allocated(work2_r)) deallocate(work2_r)
       if (allocated(work1_c)) deallocate(work1_c)
       if (allocated(work2_c)) deallocate(work2_c)
       allocate(work1_r(buf_size), STAT=status)
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
       allocate(work2_r(buf_size), STAT=status)
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
       allocate(work1_c(buf_size), STAT=status)
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
       allocate(work2_c(buf_size), STAT=status)
       if (status /= 0) then
          errorcode = 2
          call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
               'Out of memory when allocating 2DECOMP workspace')
       end if
    end if

    ! In case of MPI3 shared memory and proc is local master
    if (DECOMP_2D_LOCALCOMM /= MPI_COMM_NULL) call decomp_info_init_local(decomp)

    ! Reshape Y and Z pencils
    call decomp_info_init_reshapeyz(decomp)

    ! Compute the local index
    if (DECOMP_2D_LOCALCOMM == MPI_COMM_NULL) then
       ! No MPI3 shared memory
       decomp%xst_loc = decomp%xst
       decomp%xen_loc = decomp%xen
       decomp%xsz_loc = decomp%xsz
       decomp%yst_loc = decomp%yst                                                                   
       decomp%yen_loc = decomp%yen                                                                   
       decomp%ysz_loc = decomp%ysz
       decomp%zst_loc = decomp%zst                                                                   
       decomp%zen_loc = decomp%zen                                                                   
       decomp%zsz_loc = decomp%zsz
    else
       ! MPI3 shared memory
       call partition(decomp%xsz(1),1,decomp%xsz(2)*decomp%xsz(3), & 
            (/1,2,3/), dims_loc, coord_loc, &
            decomp%xst_loc,decomp%xen_loc,decomp%xsz_loc)
       call partition(decomp%ysz(1),1,decomp%ysz(2)*decomp%ysz(3), &
            (/1,2,3/), dims_loc, coord_loc, &
            decomp%yst_loc,decomp%yen_loc,decomp%ysz_loc)
       call partition(decomp%zsz(1),1,decomp%zsz(2)*decomp%zsz(3), &
            (/1,2,3/), dims_loc, coord_loc, &
            decomp%zst_loc,decomp%zen_loc,decomp%zsz_loc)
    endif

    return
  end subroutine decomp_info_init

  !
  ! If MPI3 shared memory, local master should broadcast decomp_info stuff
  !
  subroutine decomp_info_init_local(decomp)

     implicit none

     type(decomp_info), intent(inout) :: decomp

     integer :: ierror

     call MPI_BCAST(decomp%xst, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(decomp%xen, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(decomp%xsz, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")

     call MPI_BCAST(decomp%yst, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(decomp%yen, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(decomp%ysz, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")

     call MPI_BCAST(decomp%zst, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(decomp%zen, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(decomp%zsz, 3, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")

     call MPI_BCAST(decomp%x1count, 1, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(decomp%y1count, 1, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(decomp%y2count, 1, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
     call MPI_BCAST(decomp%z2count, 1, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)

     call MPI_BCAST(decomp%even, 1, MPI_LOGICAL, 0, DECOMP_2D_LOCALCOMM, ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")

  end subroutine decomp_info_init_local

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Release memory associated with a DECOMP_INFO object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_info_finalize(decomp)

    implicit none

    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    ! If MPI3 shared memory and rank is not local master
    if (DECOMP_2D_COMM == MPI_COMM_NULL) return

    if (allocated(decomp%x1dist)) deallocate(decomp%x1dist)
    if (allocated(decomp%y1dist)) deallocate(decomp%y1dist)
    if (allocated(decomp%y2dist)) deallocate(decomp%y2dist)
    if (allocated(decomp%z2dist)) deallocate(decomp%z2dist)
    if (allocated(decomp%x1cnts)) deallocate(decomp%x1cnts)
    if (allocated(decomp%y1cnts)) deallocate(decomp%y1cnts)
    if (allocated(decomp%y2cnts)) deallocate(decomp%y2cnts)
    if (allocated(decomp%z2cnts)) deallocate(decomp%z2cnts)
    if (allocated(decomp%x1disp)) deallocate(decomp%x1disp)
    if (allocated(decomp%y1disp)) deallocate(decomp%y1disp)
    if (allocated(decomp%y2disp)) deallocate(decomp%y2disp)
    if (allocated(decomp%z2disp)) deallocate(decomp%z2disp)

    return
  end subroutine decomp_info_finalize


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Coarser mesh support for statistic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_coarser_mesh_statS(i_skip,j_skip,k_skip,from1)

    implicit none

    integer, intent(IN) :: i_skip,j_skip,k_skip
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
    ! .false. - save n,2n,3n...

    integer, dimension(3) :: skip
    integer :: i

    coarse_mesh_starts_from_1 = from1
    iskipS = i_skip
    jskipS = j_skip
    kskipS = k_skip

    skip(1)=iskipS
    skip(2)=jskipS
    skip(3)=kskipS

    do i=1,3
       if (from1) then
          xstS(i) = (xstart(i)+skip(i)-1)/skip(i)
          if (mod(xstart(i)+skip(i)-1,skip(i))/=0) xstS(i)=xstS(i)+1
          xenS(i) = (xend(i)+skip(i)-1)/skip(i)
       else
          xstS(i) = xstart(i)/skip(i)
          if (mod(xstart(i),skip(i))/=0) xstS(i)=xstS(i)+1
          xenS(i) = xend(i)/skip(i)
       end if
       xszS(i) = xenS(i)-xstS(i)+1
    end do

    do i=1,3
       if (from1) then
          ystS(i) = (ystart(i)+skip(i)-1)/skip(i)
          if (mod(ystart(i)+skip(i)-1,skip(i))/=0) ystS(i)=ystS(i)+1
          yenS(i) = (yend(i)+skip(i)-1)/skip(i)
       else
          ystS(i) = ystart(i)/skip(i)
          if (mod(ystart(i),skip(i))/=0) ystS(i)=ystS(i)+1
          yenS(i) = yend(i)/skip(i)
       end if
       yszS(i) = yenS(i)-ystS(i)+1
    end do

    do i=1,3
       if (from1) then
          zstS(i) = (zstart(i)+skip(i)-1)/skip(i)
          if (mod(zstart(i)+skip(i)-1,skip(i))/=0) zstS(i)=zstS(i)+1
          zenS(i) = (zend(i)+skip(i)-1)/skip(i)
       else
          zstS(i) = zstart(i)/skip(i)
          if (mod(zstart(i),skip(i))/=0) zstS(i)=zstS(i)+1
          zenS(i) = zend(i)/skip(i)
       end if
       zszS(i) = zenS(i)-zstS(i)+1
    end do

    return
  end subroutine init_coarser_mesh_statS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Coarser mesh support for visualization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_coarser_mesh_statV(i_skip,j_skip,k_skip,from1)

    implicit none

    integer, intent(IN) :: i_skip,j_skip,k_skip
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
    ! .false. - save n,2n,3n...

    integer, dimension(3) :: skip
    integer :: i

    coarse_mesh_starts_from_1 = from1
    iskipV = i_skip
    jskipV = j_skip
    kskipV = k_skip

    skip(1)=iskipV
    skip(2)=jskipV
    skip(3)=kskipV

    do i=1,3
       if (from1) then
          xstV(i) = (xstart(i)+skip(i)-1)/skip(i)
          if (mod(xstart(i)+skip(i)-1,skip(i))/=0) xstV(i)=xstV(i)+1
          xenV(i) = (xend(i)+skip(i)-1)/skip(i)
       else
          xstV(i) = xstart(i)/skip(i)
          if (mod(xstart(i),skip(i))/=0) xstV(i)=xstV(i)+1
          xenV(i) = xend(i)/skip(i)
       end if
       xszV(i) = xenV(i)-xstV(i)+1
    end do

    do i=1,3
       if (from1) then
          ystV(i) = (ystart(i)+skip(i)-1)/skip(i)
          if (mod(ystart(i)+skip(i)-1,skip(i))/=0) ystV(i)=ystV(i)+1
          yenV(i) = (yend(i)+skip(i)-1)/skip(i)
       else
          ystV(i) = ystart(i)/skip(i)
          if (mod(ystart(i),skip(i))/=0) ystV(i)=ystV(i)+1
          yenV(i) = yend(i)/skip(i)
       end if
       yszV(i) = yenV(i)-ystV(i)+1
    end do

    do i=1,3
       if (from1) then
          zstV(i) = (zstart(i)+skip(i)-1)/skip(i)
          if (mod(zstart(i)+skip(i)-1,skip(i))/=0) zstV(i)=zstV(i)+1
          zenV(i) = (zend(i)+skip(i)-1)/skip(i)
       else
          zstV(i) = zstart(i)/skip(i)
          if (mod(zstart(i),skip(i))/=0) zstV(i)=zstV(i)+1
          zenV(i) = zend(i)/skip(i)
       end if
       zszV(i) = zenV(i)-zstV(i)+1
    end do

    return
  end subroutine init_coarser_mesh_statV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Coarser mesh support for probe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_coarser_mesh_statP(i_skip,j_skip,k_skip,from1)

    implicit none

    integer, intent(IN) :: i_skip,j_skip,k_skip
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
    ! .false. - save n,2n,3n...

    integer, dimension(3) :: skip
    integer :: i

    coarse_mesh_starts_from_1 = from1
    iskipP = i_skip
    jskipP = j_skip
    kskipP = k_skip

    skip(1)=iskipP
    skip(2)=jskipP
    skip(3)=kskipP

    do i=1,3
       if (from1) then
          xstP(i) = (xstart(i)+skip(i)-1)/skip(i)
          if (mod(xstart(i)+skip(i)-1,skip(i))/=0) xstP(i)=xstP(i)+1
          xenP(i) = (xend(i)+skip(i)-1)/skip(i)
       else
          xstP(i) = xstart(i)/skip(i)
          if (mod(xstart(i),skip(i))/=0) xstP(i)=xstP(i)+1
          xenP(i) = xend(i)/skip(i)
       end if
       xszP(i) = xenP(i)-xstP(i)+1
    end do

    do i=1,3
       if (from1) then
          ystP(i) = (ystart(i)+skip(i)-1)/skip(i)
          if (mod(ystart(i)+skip(i)-1,skip(i))/=0) ystP(i)=ystP(i)+1
          yenP(i) = (yend(i)+skip(i)-1)/skip(i)
       else
          ystP(i) = ystart(i)/skip(i)
          if (mod(ystart(i),skip(i))/=0) ystP(i)=ystP(i)+1
          yenP(i) = yend(i)/skip(i)
       end if
       yszP(i) = yenP(i)-ystP(i)+1
    end do

    do i=1,3
       if (from1) then
          zstP(i) = (zstart(i)+skip(i)-1)/skip(i)
          if (mod(zstart(i)+skip(i)-1,skip(i))/=0) zstP(i)=zstP(i)+1
          zenP(i) = (zend(i)+skip(i)-1)/skip(i)
       else
          zstP(i) = zstart(i)/skip(i)
          if (mod(zstart(i),skip(i))/=0) zstP(i)=zstP(i)+1
          zenP(i) = zend(i)/skip(i)
       end if
       zszP(i) = zenP(i)-zstP(i)+1
    end do

    return
  end subroutine init_coarser_mesh_statP

  ! Copy data from a fine-resolution array to a coarse one for statistic
  subroutine fine_to_coarseS(ipencil,var_fine,var_coarse)

    implicit none

    real(mytype), dimension(:,:,:) :: var_fine
    real(mytype), dimension(:,:,:) :: var_coarse
    integer, intent(IN) :: ipencil

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer :: i,j,k

    if (ipencil==1) then
       allocate(wk(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
       allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=xstS(3),xenS(3)
             do j=xstS(2),xenS(2)
                do i=xstS(1),xenS(1)
                   wk(i,j,k) = wk2((i-1)*iskipS+1,(j-1)*jskipS+1,(k-1)*kskipS+1)
                end do
             end do
          end do
       else
          do k=xstS(3),xenS(3)
             do j=xstS(2),xenS(2)
                do i=xstS(1),xenS(1)
                   wk(i,j,k) = wk2(i*iskipS,j*jskipS,k*kskipS)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==2) then
       allocate(wk(ystS(1):yenS(1),ystS(2):yenS(2),ystS(3):yenS(3)))
       allocate(wk2(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=ystS(3),yenS(3)
             do j=ystS(2),yenS(2)
                do i=ystS(1),yenS(1)
                   wk(i,j,k) = wk2((i-1)*iskipS+1,(j-1)*jskipS+1,(k-1)*kskipS+1)
                end do
             end do
          end do
       else
          do k=ystS(3),yenS(3)
             do j=ystS(2),yenS(2)
                do i=ystS(1),yenS(1)
                   wk(i,j,k) = wk2(i*iskipS,j*jskipS,k*kskipS)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==3) then
       allocate(wk(zstS(1):zenS(1),zstS(2):zenS(2),zstS(3):zenS(3)))
       allocate(wk2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=zstS(3),zenS(3)
             do j=zstS(2),zenS(2)
                do i=zstS(1),zenS(1)
                   wk(i,j,k) = wk2((i-1)*iskipS+1,(j-1)*jskipS+1,(k-1)*kskipS+1)
                end do
             end do
          end do
       else
          do k=zstS(3),zenS(3)
             do j=zstS(2),zenS(2)
                do i=zstS(1),zenS(1)
                   wk(i,j,k) = wk2(i*iskipS,j*jskipS,k*kskipS)
                end do
             end do
          end do
       end if
       var_coarse=wk
    end if

    deallocate(wk,wk2)

    return
  end subroutine fine_to_coarseS

  ! Copy data from a fine-resolution array to a coarse one for visualization
  subroutine fine_to_coarseV(ipencil,var_fine,var_coarse)

    implicit none

    real(mytype), dimension(:,:,:) :: var_fine
    real(mytype), dimension(:,:,:) :: var_coarse
    integer, intent(IN) :: ipencil

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer :: i,j,k

    if (ipencil==1) then
       allocate(wk(xstV(1):xenV(1),xstV(2):xenV(2),xstV(3):xenV(3)))
       allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=xstV(3),xenV(3)
             do j=xstV(2),xenV(2)
                do i=xstV(1),xenV(1)
                   wk(i,j,k) = wk2((i-1)*iskipV+1,(j-1)*jskipV+1,(k-1)*kskipV+1)
                end do
             end do
          end do
       else
          do k=xstV(3),xenV(3)
             do j=xstV(2),xenV(2)
                do i=xstV(1),xenV(1)
                   wk(i,j,k) = wk2(i*iskipV,j*jskipV,k*kskipV)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==2) then
       allocate(wk(ystV(1):yenV(1),ystV(2):yenV(2),ystV(3):yenV(3)))
       allocate(wk2(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=ystV(3),yenV(3)
             do j=ystV(2),yenV(2)
                do i=ystV(1),yenV(1)
                   wk(i,j,k) = wk2((i-1)*iskipV+1,(j-1)*jskipV+1,(k-1)*kskipV+1)
                end do
             end do
          end do
       else
          do k=ystV(3),yenV(3)
             do j=ystV(2),yenV(2)
                do i=ystV(1),yenV(1)
                   wk(i,j,k) = wk2(i*iskipV,j*jskipV,k*kskipV)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==3) then
       allocate(wk(zstV(1):zenV(1),zstV(2):zenV(2),zstV(3):zenV(3)))
       allocate(wk2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=zstV(3),zenV(3)
             do j=zstV(2),zenV(2)
                do i=zstV(1),zenV(1)
                   wk(i,j,k) = wk2((i-1)*iskipV+1,(j-1)*jskipV+1,(k-1)*kskipV+1)
                end do
             end do
          end do
       else
          do k=zstV(3),zenV(3)
             do j=zstV(2),zenV(2)
                do i=zstV(1),zenV(1)
                   wk(i,j,k) = wk2(i*iskipV,j*jskipV,k*kskipV)
                end do
             end do
          end do
       end if
       var_coarse=wk
    end if

    deallocate(wk,wk2)

    return
  end subroutine fine_to_coarseV

  ! Copy data from a fine-resolution array to a coarse one for probe
  subroutine fine_to_coarseP(ipencil,var_fine,var_coarse)

    implicit none

    real(mytype), dimension(:,:,:) :: var_fine
    real(mytype), dimension(:,:,:) :: var_coarse
    integer, intent(IN) :: ipencil

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer :: i,j,k

    if (ipencil==1) then
       allocate(wk(xstP(1):xenP(1),xstP(2):xenP(2),xstP(3):xenP(3)))
       allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=xstP(3),xenP(3)
             do j=xstP(2),xenP(2)
                do i=xstP(1),xenP(1)
                   wk(i,j,k) = wk2((i-1)*iskipP+1,(j-1)*jskipP+1,(k-1)*kskipP+1)
                end do
             end do
          end do
       else
          do k=xstP(3),xenP(3)
             do j=xstP(2),xenP(2)
                do i=xstP(1),xenP(1)
                   wk(i,j,k) = wk2(i*iskipP,j*jskipP,k*kskipP)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==2) then
       allocate(wk(ystP(1):yenP(1),ystP(2):yenP(2),ystP(3):yenP(3)))
       allocate(wk2(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=ystP(3),yenP(3)
             do j=ystP(2),yenP(2)
                do i=ystP(1),yenP(1)
                   wk(i,j,k) = wk2((i-1)*iskipP+1,(j-1)*jskipP+1,(k-1)*kskipP+1)
                end do
             end do
          end do
       else
          do k=ystP(3),yenP(3)
             do j=ystP(2),yenP(2)
                do i=ystP(1),yenP(1)
                   wk(i,j,k) = wk2(i*iskipP,j*jskipP,k*kskipP)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==3) then
       allocate(wk(zstP(1):zenP(1),zstP(2):zenP(2),zstP(3):zenP(3)))
       allocate(wk2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=zstP(3),zenP(3)
             do j=zstP(2),zenP(2)
                do i=zstP(1),zenP(1)
                   wk(i,j,k) = wk2((i-1)*iskipP+1,(j-1)*jskipP+1,(k-1)*kskipP+1)
                end do
             end do
          end do
       else
          do k=zstP(3),zenP(3)
             do j=zstP(2),zenP(2)
                do i=zstP(1),zenP(1)
                   wk(i,j,k) = wk2(i*iskipP,j*jskipP,k*kskipP)
                end do
             end do
          end do
       end if
       var_coarse=wk
    end if

    deallocate(wk,wk2)

    return
  end subroutine fine_to_coarseP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Find sub-domain information held by current processor
  !   INPUT: 
  !     nx, ny, nz - global data dimension
  !     pdim(3)    - number of processor grid in each dimension, 
  !                  valid values: 1 - distibute locally; 
  !                                2 - distribute across p_row; 
  !                                3 - distribute across p_col
  !     dims(2)    - (p_row, p_col)
  !     coord(2)   - coordinates of the CPU on the grid (p_row, p_col)
  !   OUTPUT:
  !     lstart(3)  - starting index
  !     lend(3)    - ending index
  !     lsize(3)   - size of the sub-block (redundant) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine partition(nx, ny, nz, pdim, dims, coord, lstart, lend, lsize)

    implicit none

    integer, intent(IN) :: nx, ny, nz
    integer, dimension(3), intent(IN) :: pdim
    integer, dimension(2), intent(IN) :: dims, coord
    integer, dimension(3), intent(OUT) :: lstart, lend, lsize

    integer, allocatable, dimension(:) :: st,en,sz
    integer :: i, gsize

    do i = 1, 3

       if (i==1) then
          gsize = nx
       else if (i==2) then
          gsize = ny
       else if (i==3) then
          gsize = nz
       end if

       if (pdim(i) == 1) then        ! all local
          lstart(i) = 1
          lend(i)   = gsize
          lsize(i)  = gsize
       elseif (pdim(i) == 2) then    ! distribute across dims(1)
          allocate(st(0:dims(1)-1))
          allocate(en(0:dims(1)-1))
          allocate(sz(0:dims(1)-1))
          call distribute(gsize,dims(1),st,en,sz)
          lstart(i) = st(coord(1))
          lend(i)   = en(coord(1))
          lsize(i)  = sz(coord(1))
          deallocate(st,en,sz)
       elseif (pdim(i) == 3) then    ! distribute across dims(2)
          allocate(st(0:dims(2)-1))
          allocate(en(0:dims(2)-1))
          allocate(sz(0:dims(2)-1))
          call distribute(gsize,dims(2),st,en,sz)
          lstart(i) = st(coord(2))
          lend(i)   = en(coord(2))
          lsize(i)  = sz(coord(2))
          deallocate(st,en,sz)
       end if

    end do
    return   

  end subroutine partition

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   - distibutes grid points in one dimension
  !   - handles uneven distribution properly 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine distribute(data1,proc,st,en,sz)

    implicit none
    ! data1 -- data size in any dimension to be partitioned
    ! proc  -- number of processors in that dimension
    ! st    -- array of starting index
    ! en    -- array of ending index
    ! sz    -- array of local size  (redundent)
    integer data1,proc,st(0:proc-1),en(0:proc-1),sz(0:proc-1)
    integer i,size1,nl,nu

    size1=data1/proc
    nu = data1 - size1 * proc
    nl = proc - nu
    st(0) = 1
    sz(0) = size1
    en(0) = size1
    do i=1,nl-1
       st(i) = st(i-1) + size1
       sz(i) = size1
       en(i) = en(i-1) + size1
    end do
    size1 = size1 + 1
    do i=nl,proc-1
       st(i) = en(i-1) + 1
       sz(i) = size1
       en(i) = en(i-1) + size1
    end do
    en(proc-1)= data1 
    sz(proc-1)= data1-st(proc-1)+1

    return
  end subroutine distribute

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Define how each dimension is distributed across processors
  !    e.g. 17 meshes across 4 processor would be distibuted as (4,4,4,5)
  !    such global information is required locally at MPI_ALLTOALLV time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_dist(nx,ny,nz,decomp)

    implicit none

    integer, intent(IN) :: nx, ny, nz
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp
    integer, allocatable, dimension(:) :: st,en

    allocate(st(0:dims(1)-1))
    allocate(en(0:dims(1)-1))
    call distribute(nx,dims(1),st,en,decomp%x1dist)
    call distribute(ny,dims(1),st,en,decomp%y1dist)
    deallocate(st,en)

    allocate(st(0:dims(2)-1))
    allocate(en(0:dims(2)-1))
    call distribute(ny,dims(2),st,en,decomp%y2dist)
    call distribute(nz,dims(2),st,en,decomp%z2dist)
    deallocate(st,en)

    return
  end subroutine get_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Prepare the send / receive buffers for MPI_ALLTOALLV communications
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine prepare_buffer(decomp)

    implicit none

    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    integer :: i

    !LG : AJOUTS "bidons" pour eviter un plantage en -O3 avec gcc9.3
    !       * la fonction sortait des valeurs 'aleatoires'
    !         et le calcul plantait dans MPI_ALLTOALLV
    !       * pas de plantage en O2
    
    character(len=100) :: tmp_char
    if (nrank==0) then
       open(newunit=i, file='temp.dat', form='unformatted')
       write(i) decomp%x1dist,decomp%y1dist,decomp%y2dist,decomp%z2dist, &
                decomp%xsz,decomp%ysz,decomp%zsz
       close(i, status='delete')
    endif

    ! MPI_ALLTOALLV buffer information

    do i=0, dims(1)-1
       decomp%x1cnts(i) = decomp%x1dist(i)*decomp%xsz(2)*decomp%xsz(3)
       decomp%y1cnts(i) = decomp%ysz(1)*decomp%y1dist(i)*decomp%ysz(3)
       if (i==0) then
          decomp%x1disp(i) = 0  ! displacement is 0-based index
          decomp%y1disp(i) = 0
       else
          decomp%x1disp(i) = decomp%x1disp(i-1) + decomp%x1cnts(i-1)
          decomp%y1disp(i) = decomp%y1disp(i-1) + decomp%y1cnts(i-1)
       end if
    end do

    do i=0, dims(2)-1
       decomp%y2cnts(i) = decomp%ysz(1)*decomp%y2dist(i)*decomp%ysz(3)
       decomp%z2cnts(i) = decomp%zsz(1)*decomp%zsz(2)*decomp%z2dist(i)
       if (i==0) then
          decomp%y2disp(i) = 0  ! displacement is 0-based index
          decomp%z2disp(i) = 0
       else
          decomp%y2disp(i) = decomp%y2disp(i-1) + decomp%y2cnts(i-1)
          decomp%z2disp(i) = decomp%z2disp(i-1) + decomp%z2cnts(i-1)
       end if
    end do

    ! MPI_ALLTOALL buffer information

    ! For evenly distributed data, following is an easier implementation.
    ! But it should be covered by the more general formulation below.
    !decomp%x1count = decomp%xsz(1)*decomp%xsz(2)*decomp%xsz(3)/dims(1)
    !decomp%y1count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(1) 
    !decomp%y2count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(2)
    !decomp%z2count = decomp%zsz(1)*decomp%zsz(2)*decomp%zsz(3)/dims(2)

    ! For unevenly distributed data, pad smaller messages. Note the 
    ! last blocks along pencils always get assigned more mesh points
    ! for X <=> Y transposes
    decomp%x1count = decomp%x1dist(dims(1)-1) * &
         decomp%y1dist(dims(1)-1) * decomp%xsz(3)
    decomp%y1count = decomp%x1count
    ! for Y <=> Z transposes
    decomp%y2count = decomp%y2dist(dims(2)-1) * &
         decomp%z2dist(dims(2)-1) * decomp%zsz(1)
    decomp%z2count = decomp%y2count

    return
  end subroutine prepare_buffer


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

    if (nrank==0) write(*,*) 'In auto-tuning mode......'

    i = int(sqrt(real(iproc))) + 10  ! enough space to save all factors 
    allocate(factors(i))
    call findfactor(iproc, factors, nfact)
    if (nrank==0) write(*,*) 'factors: ', (factors(i), i=1,nfact)

    i_best=nfact/2+1
    col=factors(i_best)

    best_p_col = col
    best_p_row=iproc/col
    if (nrank==0) print *,'p_row x p_col', best_p_row, best_p_col
    if ((best_p_col==1).and.(nrank==0)) then
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

    if (nrank==0) then
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

    if (nrank==0) then
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

    if (nrank==0) then
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

    if (nrank==0) then
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Utility routines to help allocate 3D arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "alloc.f90"
    
  
end module decomp_2d

