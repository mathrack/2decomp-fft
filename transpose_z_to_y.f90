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

! This file contains the routines that transpose data from Z to Y pencil

  subroutine transpose_z_to_y_data(src, dst)

    implicit none

    ! Arguments
    type(decomp_data), intent(in) :: src
    type(decomp_data), intent(inout) :: dst

    ! Local variables
    integer :: s1, s2, s3, d1, d2, d3, ierror
#if defined(_GPU)
#if defined(_NCCL)
    integer :: row_rank_id
#endif
#endif

    ! Safety check
#ifdef DEBUG
    if (src%idir /= 3 .or. dst%idir /= 2 .or. &
        (src%is_cplx.neqv.dst%is_cplx) .or. &
        (src%shm.neqv.dst%shm) .or. &
        .not.associated(src%decomp, dst%decomp)) &
       call decomp_2d_abort(__FILE__, __LINE__, -1, "Impossible transpose operation")
#endif

    s1 = src%decomp%zsz(1)
    s2 = src%decomp%zsz(2)
    s3 = src%decomp%zsz(3)
    d1 = src%decomp%ysz(1)
    d2 = src%decomp%ysz(2)
    d3 = src%decomp%ysz(3)

    ! In case of MPI3 shared memory, local master starts reading
    if (d2d_intranode) then
       call decomp_2d_win_transpose_start_reading(src%win)
       if (src%is_cplx) then
          call decomp_2d_win_transpose_start_writing(work1_c_win)
       else
          call decomp_2d_win_transpose_start_writing(work1_r_win)
       endif
    endif

    ! rearrange source array as send buffer
    if (src%is_cplx) then
#if defined(GPU)
       istat = cudaMemcpy( work1_c_d, src%cvar, s1*s2*s3 )
#else
       call mem_split_zy_complex(src, s1, s2, s3, work1_c, dims(2), &
            src%decomp%z2dist, src%decomp)
#endif
    else
#if defined(GPU)
       istat = cudaMemcpy( work1_r_d, src%var, s1*s2*s3 )
#else
       call mem_split_zy_real(src, s1, s2, s3, work1_r, dims(2), &
            src%decomp%z2dist, src%decomp)
#endif
    endif

    ! In case of MPI3 shared memory, local master is done reading
    if (d2d_intranode) then
       call decomp_2d_win_transpose_stop_reading(src%win)
       if (src%is_cplx) then
          call decomp_2d_win_transpose_stop_writing(work1_c_win)
       else
          call decomp_2d_win_transpose_stop_writing(work1_r_win)
       endif
    endif

    ! transpose using MPI_ALLTOALL(V)
    if (nrank_loc <= 0) then
    if (src%is_cplx) then
#ifdef EVEN
       call MPI_ALLTOALL(work1_c, src%decomp%z2count, &
            complex_type, work2_c, src%decomp%y2count, &
            complex_type, DECOMP_2D_COMM_ROW, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
#else

#if defined(_GPU)
       call MPI_ALLTOALLV(work1_c_d, src%decomp%z2cnts, src%decomp%z2disp, &
            complex_type, work2_c_d, src%decomp%y2cnts, src%decomp%y2disp, &
            complex_type, DECOMP_2D_COMM_ROW, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#else
       call MPI_ALLTOALLV(work1_c, src%decomp%z2cnts, src%decomp%z2disp, &
            complex_type, work2_c, src%decomp%y2cnts, src%decomp%y2disp, &
            complex_type, DECOMP_2D_COMM_ROW, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif
#endif
    else
#ifdef EVEN
       call MPI_ALLTOALL(work1_r, src%decomp%z2count, &
            real_type, work2_r, src%decomp%y2count, &
            real_type, DECOMP_2D_COMM_ROW, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
#else
#if defined(_GPU)
#if defined(_NCCL)
       nccl_stat = ncclGroupStart()
       do row_rank_id = 0, (row_comm_size - 1)
           nccl_stat = ncclSend(work1_r_d( src%decomp%z2disp(row_rank_id)+1 ), src%decomp%z2cnts(row_rank_id), ncclDouble, local_to_global_row(row_rank_id+1), nccl_comm_2decomp, cuda_stream_2decomp)
           nccl_stat = ncclRecv(work2_r_d( src%decomp%y2disp(row_rank_id)+1 ), src%decomp%y2cnts(row_rank_id), ncclDouble, local_to_global_row(row_rank_id+1), nccl_comm_2decomp, cuda_stream_2decomp)
       end do
       nccl_stat = ncclGroupEnd()
       cuda_stat = cudaStreamSynchronize(cuda_stream_2decomp)
#else
       call MPI_ALLTOALLV(work1_r_d, src%decomp%z2cnts, src%decomp%z2disp, &
            real_type, work2_r_d, src%decomp%y2cnts, src%decomp%y2disp, &
            real_type, DECOMP_2D_COMM_ROW, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif
#else
       call MPI_ALLTOALLV(work1_r, src%decomp%z2cnts, src%decomp%z2disp, &
            real_type, work2_r, src%decomp%y2cnts, src%decomp%y2disp, &
            real_type, DECOMP_2D_COMM_ROW, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif
#endif
    endif
    endif

    ! In case of MPI3 shared memory, local master starts writing
    if (d2d_intranode) then
       call decomp_2d_win_transpose_start_writing(dst%win)
       if (src%is_cplx) then
          call decomp_2d_win_transpose_start_reading(work2_c_win)
       else
          call decomp_2d_win_transpose_start_reading(work2_r_win)
       endif
    endif

    ! rearrange receive buffer
    if (src%is_cplx) then
#if defined (_GPU)
       call mem_merge_zy_complex(work2_c_d, d1, d2, d3, dst%cvar, dims(2), &
            src%decomp%y2dist, src%decomp)
#else
       call mem_merge_zy_complex(work2_c, d1, d2, d3, dst, dims(2), &
            src%decomp%y2dist, src%decomp)
#endif
    else
#if defined (_GPU)
       call mem_merge_zy_real(work2_r_d, d1, d2, d3, dst%var, dims(2), &
            src%decomp%y2dist, src%decomp)
#else
       call mem_merge_zy_real(work2_r, d1, d2, d3, dst, dims(2), &
            src%decomp%y2dist, src%decomp)
#endif
    endif

    ! In case of MPI3 shared memory, local master is done writing
    if (d2d_intranode) then
       call decomp_2d_win_transpose_stop_writing(dst%win)
       if (src%is_cplx) then
          call decomp_2d_win_transpose_stop_reading(work2_c_win)
       else
          call decomp_2d_win_transpose_stop_reading(work2_r_win)
       endif
    endif

  end subroutine transpose_z_to_y_data

  subroutine transpose_z_to_y_real(src, dst, opt_decomp, src_win, dst_win)

    implicit none
 
    ! Arguments
    real(mytype), dimension(:,:,:), pointer, intent(IN) :: src
    real(mytype), dimension(:,:,:), pointer, intent(INOUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    integer, intent(in), optional :: src_win, dst_win

    ! Local variables
    type(decomp_data) :: data_src, data_dst

    ! Init using given array
    call data_src%init(is_cplx = .false., idir = 3, decomp = opt_decomp, rwk = src)
    call data_dst%init(is_cplx = .false., idir = 2, decomp = opt_decomp, rwk = dst)
    if (present(src_win)) data_src%win = src_win
    if (present(dst_win)) data_dst%win = dst_win
    ! Transpose
    call transpose_z_to_y(data_src, data_dst)
    ! Clean
    nullify(data_src%decomp)
    nullify(data_src%var)
    nullify(data_dst%decomp)
    nullify(data_dst%var)
    data_src%win = MPI_WIN_NULL
    data_dst%win = MPI_WIN_NULL

  end subroutine transpose_z_to_y_real

  subroutine transpose_z_to_y_complex(src, dst, opt_decomp, src_win, dst_win)

    implicit none

    ! Arguments
    complex(mytype), dimension(:,:,:), pointer, intent(IN) :: src
    complex(mytype), dimension(:,:,:), pointer, intent(INOUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    integer, intent(in), optional :: src_win, dst_win

    ! Local variables
    type(decomp_data) :: data_src, data_dst
    
    ! Init using given array
    call data_src%init(is_cplx = .true., idir = 3, decomp = opt_decomp, cwk = src)
    call data_dst%init(is_cplx = .true., idir = 2, decomp = opt_decomp, cwk = dst)
    if (present(src_win)) data_src%win = src_win
    if (present(dst_win)) data_dst%win = dst_win
    ! Transpose
    call transpose_z_to_y(data_src, data_dst)
    ! Clean
    nullify(data_src%decomp)
    nullify(data_src%cvar)
    nullify(data_dst%decomp)
    nullify(data_dst%cvar)
    data_src%win = MPI_WIN_NULL
    data_dst%win = MPI_WIN_NULL
    
  end subroutine transpose_z_to_y_complex

  ! pack/unpack ALLTOALL(V) buffers

  subroutine mem_split_zy_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    type(decomp_data), intent(in) :: in
    real(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2,pos

    do m=0,iproc-1
       if (m==0) then 
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef EVEN
       pos = m * decomp%z2count + 1
#else
       pos = decomp%z2disp(m) + 1
#endif
       if (d2d_intranode .and. associated(in%var2d)) then
          do j = 1, decomp%zsz_loc(3)
             do k = i1, i2
                out(decomp%intramap_split(4,k,j)) = in%var2d(k,j)
             enddo
          enddo
       else
          do j=1,n3 !k=i1,i2
             do i=1,n2 !j=1,n3
                do k=i1,i2 !i=1,n2
                   out(pos + i-1 + n2*(j-1) + n2*n3*(k-i1)) = in%var(k,i,j)
                end do
             end do
          end do
       endif
    end do

    return
  end subroutine mem_split_zy_real


  subroutine mem_split_zy_complex(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    type(decomp_data), intent(in) :: in
    complex(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2, pos

    do m=0,iproc-1
       if (m==0) then 
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef EVEN
       pos = m * decomp%z2count + 1
#else
       pos = decomp%z2disp(m) + 1
#endif
       if (d2d_intranode .and. associated(in%cvar2d)) then
          do j = 1, decomp%zsz_loc(3)
             do k = i1, i2
                out(decomp%intramap_split(4,k,j)) = in%cvar2d(k,j)
             enddo
          enddo
       else
          do j=1,n3
             do i=1,n2
                do k=i1,i2
                   out(pos + i-1 + n2*(j-1) + n2*n3*(k-i1)) = in%cvar(k,i,j)
                end do
             end do
          end do
       endif
    end do

    return
  end subroutine mem_split_zy_complex


  subroutine mem_merge_zy_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(*), intent(IN) :: in
    type(decomp_data), intent(inout) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
    attributes(device) :: in
    integer :: istat
#endif

    integer :: i,j,k, m,i1,i2,i3, pos

    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if
       i3 = i2 - i1 + 1

#ifdef EVEN
       pos = m * decomp%y2count + 1
#else
       pos = decomp%y2disp(m) + 1
#endif

#if defined(_GPU)
       istat = cudaMemcpy2D( out(1,i1,1), n1*n2, in(pos), n1*(i2-i1+1), n1*(i2-i1+1), n3, cudaMemcpyDeviceToDevice )
#else
       if (d2d_intranode .and. associated(out%var2d)) then
          do k = 1, decomp%ysz_loc(3)
             do j = i1, i2
                out%var2d(j,k) = in(decomp%intramap_merge(4,j,k))
             enddo
          enddo
       else
          do k=1,n3
             do i=1,n2
                do j=i1,i2
                   out%var(j,i,k) = in(pos + i-1 + n2*(j-i1) + n2*i3*(k-1))
                end do
             end do
          end do
       endif
#endif
    end do

    return
  end subroutine mem_merge_zy_real


  subroutine mem_merge_zy_complex(in,n1,n2,n3,out,iproc,dist,decomp)
    
    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(*), intent(IN) :: in
    type(decomp_data), intent(inout) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
#if defined(_GPU)
    attributes(device) :: in
    integer :: istat
#endif

    integer :: i,j,k, m,i1,i2,i3, pos

    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if
       i3 = i2 - i1 + 1

#ifdef EVEN
       pos = m * decomp%y2count + 1
#else
       pos = decomp%y2disp(m) + 1
#endif

#if defined(_GPU)
       istat = cudaMemcpy2D( out(1,i1,1), n1*n2, in(pos), n1*(i2-i1+1), n1*(i2-i1+1), n3, cudaMemcpyDeviceToDevice )
#else
       if (d2d_intranode .and. associated(out%cvar2d)) then
          do k = 1, decomp%ysz_loc(3)
             do j = i1, i2
                out%cvar2d(j,k) = in(decomp%intramap_merge(4,j,k))
             enddo
          enddo
       else
          do k=1,n3
             do i=1,n2
                do j=i1,i2
                   out%cvar(j,i,k) = in(pos + i-1 + n2*(j-i1) + n2*i3*(k-1))
                end do
             end do
          end do
       endif
#endif
    end do

    return
  end subroutine mem_merge_zy_complex
