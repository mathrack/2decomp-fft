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

submodule (decomp_2d) smod_decomp_info

  implicit none

contains

  !
  ! Reshape Y and Z pencils so that the leading dimension is the working one
  !
  module subroutine decomp_info_init_reshapeyz(decomp)

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Return the default decomposition object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  module subroutine get_decomp_info(decomp)

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
  module subroutine decomp_info_init(nx,ny,nz,decomp)

    implicit none

    integer, intent(IN) :: nx,ny,nz
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    integer :: buf_size, status, errorcode

    ! In case of MPI3 shared memory and proc is not local master
    if (d2d_intranode .and. nrank_loc/=0) then
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
    if (d2d_intranode) call decomp_info_init_local(decomp)

    ! Reshape Y and Z pencils
    call decomp_info_init_reshapeyz(decomp)

    ! Compute the local index
    if (.not.d2d_intranode) then
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
  module subroutine decomp_info_init_local(decomp)

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
  module subroutine decomp_info_finalize(decomp)

    implicit none

    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    ! If MPI3 shared memory and rank is not local master
    if (d2d_intranode .and. nrank_loc > 0) return

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

end submodule smod_decomp_info
