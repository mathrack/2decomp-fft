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

! This file contains the routines that transpose data from X to Y pencil

  subroutine transpose_x_to_y_data(src, dst)

    implicit none

    ! Arguments
    type(decomp_data), intent(in) :: src
    type(decomp_data), intent(inout) :: dst

    ! Local variables
    integer :: s1, s2, s3, d1, d2, d3, ierror
#if defined(_GPU)
#if defined(_NCCL)
    integer :: col_rank_id
#endif
#endif

    ! Safety check
#ifdef DEBUG
    if (src%idir /= 1 .or. dst%idir /= 2 .or. &
        (src%is_cplx.neqv.dst%is_cplx) .or. &
        (src%shm.neqv.dst%shm) .or. &
        .not.associated(src%decomp, dst%decomp)) &
       call decomp_2d_abort(__FILE__, __LINE__, -1, "Impossible transpose operation")
#endif

    s1 = src%decomp%xsz(1)
    s2 = src%decomp%xsz(2)
    s3 = src%decomp%xsz(3)
    d1 = src%decomp%ysz(1)
    d2 = src%decomp%ysz(2)
    d3 = src%decomp%ysz(3)

    ! In case of MPI3 shared memory
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
#if defined(_GPU)
       call mem_split_xy_complex(src%cvar, s1, s2, s3, work1_c_d, dims(1), &
            src%decomp%x1dist, src%decomp)
#else
       call mem_split_xy_complex(src, s1, s2, s3, work1_c, dims(1), &
            src%decomp%x1dist, src%decomp)
#endif
    else
#if defined(_GPU)
       call mem_split_xy_real(src%var, s1, s2, s3, work1_r_d, dims(1), &
            src%decomp%x1dist, src%decomp)
#else
       call mem_split_xy_real(src, s1, s2, s3, work1_r, dims(1), &
            src%decomp%x1dist, src%decomp)
#endif
    endif

    ! In case of MPI3 shared memory
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
       call MPI_ALLTOALL(work1_c, src%decomp%x1count, &
            complex_type, work2_c, src%decomp%y1count, &
            complex_type, DECOMP_2D_COMM_COL, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
#else

#if defined(_GPU)
       call MPI_ALLTOALLV(work1_c_d, src%decomp%x1cnts, src%decomp%x1disp, &
            complex_type, work2_c_d, src%decomp%y1cnts, src%decomp%y1disp, &
            complex_type, DECOMP_2D_COMM_COL, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#else
       call MPI_ALLTOALLV(work1_c, src%decomp%x1cnts, src%decomp%x1disp, &
            complex_type, work2_c, src%decomp%y1cnts, src%decomp%y1disp, &
            complex_type, DECOMP_2D_COMM_COL, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif

#endif
    else
#ifdef EVEN
       call MPI_ALLTOALL(work1_r, src%decomp%x1count, &
            real_type, work2_r, src%decomp%y1count, &
            real_type, DECOMP_2D_COMM_COL, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
#else

#if defined(_GPU)
#if defined(_NCCL)
       nccl_stat = ncclGroupStart()
       do col_rank_id = 0, (col_comm_size - 1)
           nccl_stat = ncclSend(work1_r_d( src%decomp%x1disp(col_rank_id)+1 ), src%decomp%x1cnts(col_rank_id), ncclDouble, local_to_global_col(col_rank_id+1), nccl_comm_2decomp, cuda_stream_2decomp)
           nccl_stat = ncclRecv(work2_r_d( src%decomp%y1disp(col_rank_id)+1 ), src%decomp%y1cnts(col_rank_id), ncclDouble, local_to_global_col(col_rank_id+1), nccl_comm_2decomp, cuda_stream_2decomp)
       end do
       nccl_stat = ncclGroupEnd()
       cuda_stat = cudaStreamSynchronize(cuda_stream_2decomp)
#else
       call MPI_ALLTOALLV(work1_r_d, src%decomp%x1cnts, src%decomp%x1disp, &
            real_type, work2_r_d, src%decomp%y1cnts, src%decomp%y1disp, &
            real_type, DECOMP_2D_COMM_COL, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif
#else
       call MPI_ALLTOALLV(work1_r, src%decomp%x1cnts, src%decomp%x1disp, &
            real_type, work2_r, src%decomp%y1cnts, src%decomp%y1disp, &
            real_type, DECOMP_2D_COMM_COL, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif

#endif
    endif
    endif

    ! In case of MPI3 shared memory
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
#if defined(_GPU)          
       call mem_merge_xy_complex(work2_c_d, d1, d2, d3, dst%cvar, dims(1), & 
            src%decomp%y1dist, src%decomp)
#else                 
       call mem_merge_xy_complex(work2_c, d1, d2, d3, dst, dims(1), & 
            src%decomp%y1dist, src%decomp)
#endif
    else
#if defined(_GPU)
       call mem_merge_xy_real(work2_r_d, d1, d2, d3, dst%var, dims(1), &
            src%decomp%y1dist, src%decomp)
#else
       call mem_merge_xy_real(work2_r, d1, d2, d3, dst, dims(1), &
            src%decomp%y1dist, src%decomp)
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

  end subroutine transpose_x_to_y_data

  subroutine mem_split_xy_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    ! Arguments
    integer, intent(IN) :: n1,n2,n3
    type(decomp_data), intent(in) :: in
    real(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    ! Local variables
#if defined(_GPU)
    attributes(device) :: out
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
       pos = m * decomp%x1count + 1
#else
       pos = decomp%x1disp(m) + 1
#endif

#if defined(_GPU)
       istat = cudaMemcpy2D( out(pos), i2-i1+1, in(i1,1,1), n1, i2-i1+1, n2*n3, cudaMemcpyDeviceToDevice )
#else
       if (d2d_intranode) then
          do k = 1, decomp%xsz_loc(3)
             do i = i1, i2
                out(decomp%intramap_split(1,i,k)) = in%var2d(i,k)
             enddo
          enddo
       else
          do k=1,n3
             do j=1,n2
                do i=i1,i2
                   out(pos + i-i1 + i3*(j-1) + i3*n2*(k-1)) = in%var(i,j,k)
                end do
             end do
          end do
       endif
#endif
    end do

    return
  end subroutine mem_split_xy_real


  subroutine mem_split_xy_complex(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    ! Arguments
    integer, intent(IN) :: n1,n2,n3
    type(decomp_data), intent(in) :: in
    complex(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    ! Local variables
#if defined(_GPU)
    attributes(device) :: out
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
       pos = m * decomp%x1count + 1
#else
       pos = decomp%x1disp(m) + 1
#endif

#if defined(_GPU)
       istat = cudaMemcpy2D( out(pos), i2-i1+1, in(i1,1,1), n1, i2-i1+1, n2*n3, cudaMemcpyDeviceToDevice )
#else
       if (d2d_intranode) then
          do k = 1, decomp%xsz_loc(3)
             do i = i1, i2
                out(decomp%intramap_split(1,i,k)) = in%cvar2d(i,k)
             enddo
          enddo
       else
          do k=1,n3
             do j=1,n2
                do i=i1,i2
                   out(pos + i-i1 + i3*(j-1) + i3*n2*(k-1)) = in%cvar(i,j,k)
                end do
             end do
          end do
       endif
#endif
    end do

    return
  end subroutine mem_split_xy_complex


  subroutine mem_merge_xy_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    ! Arguments
    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(*), intent(IN) :: in
    type(decomp_data), intent(inout) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    ! Local variables
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
       pos = m * decomp%y1count + 1
#else
       pos = decomp%y1disp(m) + 1
#endif

#if defined(_GPU)
       istat = cudaMemcpy2D( out(1,i1,1), n1*n2, in(pos), n1*(i2-i1+1), n1*(i2-i1+1), n3, cudaMemcpyDeviceToDevice )
#else
       if (d2d_intranode) then
          do k = 1, decomp%ysz_loc(3)
             do j = i1, i2
                out%var2d(j,k) = in(decomp%intramap_merge(1,j,k))
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
  end subroutine mem_merge_xy_real


  subroutine mem_merge_xy_complex(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    ! Arguments
    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(*), intent(IN) :: in
    type(decomp_data), intent(inout) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    ! Local variables
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
       pos = m * decomp%y1count + 1
#else
       pos = decomp%y1disp(m) + 1
#endif

#if defined(_GPU)
       istat = cudaMemcpy2D( out(1,i1,1), n1*n2, in(pos), n1*(i2-i1+1), n1*(i2-i1+1), n3, cudaMemcpyDeviceToDevice )
#else
       if (d2d_intranode) then
          do k = 1, decomp%ysz_loc(3)
             do j = i1, i2
                out%cvar2d(j,k) = in(decomp%intramap_merge(1,j,k))
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
  end subroutine mem_merge_xy_complex
