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

! This file contains the routines that transpose data from Y to X pencil

  subroutine transpose_y_to_x_real(src, dst, opt_decomp, src_win, dst_win)

    implicit none
    
    real(mytype), dimension(:,:,:), pointer, intent(IN) :: src
    real(mytype), dimension(:,:,:), pointer, intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    integer, intent(in), optional :: src_win, dst_win

    TYPE(DECOMP_INFO) :: decomp

#if defined(_GPU)
#if defined(_NCCL)
    integer :: col_rank_id
#endif
#endif

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    ! In case of MPI3 shared memory
    if (d2d_intranode.and.nrank_loc>0) then
       ! Local master will read(write) from everybody at start(end)
       call decomp_2d_win_transpose_start_reading(src_win)
       call decomp_2d_win_transpose_stop_reading(src_win)
       call decomp_2d_win_transpose_start_writing(dst_win)
       call decomp_2d_win_transpose_stop_writing(dst_win)
       return
    endif

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = decomp%ysz(1)
    s2 = decomp%ysz(2)
    s3 = decomp%ysz(3)
    d1 = decomp%xsz(1)
    d2 = decomp%xsz(2)
    d3 = decomp%xsz(3)

    ! In case of MPI3 shared memory, local master starts reading
    if (d2d_intranode) call decomp_2d_win_transpose_start_reading(src_win)

    ! rearrange source array as send buffer
#if defined(_GPU)
    call mem_split_yx_real(src, s1, s2, s3, work1_r_d, dims(1), &
         decomp%y1dist, decomp)
#else
    call mem_split_yx_real(src, s1, s2, s3, work1_r, dims(1), &
         decomp%y1dist, decomp)
#endif

    ! In case of MPI3 shared memory, local master is done reading
    if (d2d_intranode) call decomp_2d_win_transpose_stop_reading(src_win)

    ! transpose using MPI_ALLTOALL(V)
#ifdef EVEN
    call MPI_ALLTOALL(work1_r, decomp%y1count, &
         real_type, work2_r, decomp%x1count, &
         real_type, DECOMP_2D_COMM_COL, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
#else

#if defined(_GPU)
#if defined(_NCCL)
    nccl_stat = ncclGroupStart()
    do col_rank_id = 0, (col_comm_size - 1)
     nccl_stat = ncclSend(work1_r_d( decomp%y1disp(col_rank_id)+1 ), decomp%y1cnts(col_rank_id), ncclDouble, local_to_global_col(col_rank_id+1), nccl_comm_2decomp, cuda_stream_2decomp)
     nccl_stat = ncclRecv(work2_r_d( decomp%x1disp(col_rank_id)+1 ), decomp%x1cnts(col_rank_id), ncclDouble, local_to_global_col(col_rank_id+1), nccl_comm_2decomp, cuda_stream_2decomp)
    end do
    nccl_stat = ncclGroupEnd()
    cuda_stat = cudaStreamSynchronize(cuda_stream_2decomp)
#else
    call MPI_ALLTOALLV(work1_r_d, decomp%y1cnts, decomp%y1disp, &
         real_type, work2_r_d, decomp%x1cnts, decomp%x1disp, &
         real_type, DECOMP_2D_COMM_COL, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif
#else
    call MPI_ALLTOALLV(work1_r, decomp%y1cnts, decomp%y1disp, &
         real_type, work2_r, decomp%x1cnts, decomp%x1disp, &
         real_type, DECOMP_2D_COMM_COL, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif

#endif

    ! In case of MPI3 shared memory, local master starts writing
    if (d2d_intranode) call decomp_2d_win_transpose_start_writing(dst_win)

    ! rearrange receive buffer
#if defined(_GPU)
    call mem_merge_yx_real(work2_r_d, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#else
    call mem_merge_yx_real(work2_r, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#endif

    ! In case of MPI3 shared memory, local master is done writing
    if (d2d_intranode) call decomp_2d_win_transpose_stop_writing(dst_win)

    return
  end subroutine transpose_y_to_x_real


  subroutine transpose_y_to_x_complex(src, dst, opt_decomp, src_win, dst_win)

    implicit none
    
    complex(mytype), dimension(:,:,:), pointer, intent(IN) :: src
    complex(mytype), dimension(:,:,:), pointer, intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    integer, intent(in), optional :: src_win, dst_win

    TYPE(DECOMP_INFO) :: decomp

#if defined(_GPU)
#if defined(_NCCL)
    integer :: col_rank_id
#endif
#endif

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    ! In case of MPI3 shared memory
    if (d2d_intranode.and.nrank_loc>0) then
       ! Local master will read(write) from everybody at start(end)
       call decomp_2d_win_transpose_start_reading(src_win)
       call decomp_2d_win_transpose_stop_reading(src_win)
       call decomp_2d_win_transpose_start_writing(dst_win)
       call decomp_2d_win_transpose_stop_writing(dst_win)
       return
    endif

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = decomp%ysz(1)
    s2 = decomp%ysz(2)
    s3 = decomp%ysz(3)
    d1 = decomp%xsz(1)
    d2 = decomp%xsz(2)
    d3 = decomp%xsz(3)
 
    ! In case of MPI3 shared memory, local master starts reading
    if (d2d_intranode) call decomp_2d_win_transpose_start_reading(src_win)

    ! rearrange source array as send buffer
#if defined(_GPU)
    call mem_split_yx_complex(src, s1, s2, s3, work1_c_d, dims(1), &
         decomp%y1dist, decomp)
#else
    call mem_split_yx_complex(src, s1, s2, s3, work1_c, dims(1), &
         decomp%y1dist, decomp)
#endif

    ! In case of MPI3 shared memory, local master is done reading
    if (d2d_intranode) call decomp_2d_win_transpose_stop_reading(src_win)

    ! transpose using MPI_ALLTOALL(V)
#ifdef EVEN
    call MPI_ALLTOALL(work1_c, decomp%y1count, &
         complex_type, work2_c, decomp%x1count, &
         complex_type, DECOMP_2D_COMM_COL, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
#else

#if defined(_GPU)
    call MPI_ALLTOALLV(work1_c_d, decomp%y1cnts, decomp%y1disp, &
         complex_type, work2_c_d, decomp%x1cnts, decomp%x1disp, &
         complex_type, DECOMP_2D_COMM_COL, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#else
    call MPI_ALLTOALLV(work1_c, decomp%y1cnts, decomp%y1disp, &
         complex_type, work2_c, decomp%x1cnts, decomp%x1disp, &
         complex_type, DECOMP_2D_COMM_COL, ierror)
    if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif

#endif

    ! In case of MPI3 shared memory, local master starts writing
    if (d2d_intranode) call decomp_2d_win_transpose_start_writing(dst_win)

    ! rearrange receive buffer
#if defined(_GPU)
    call mem_merge_yx_complex(work2_c_d, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#else
    call mem_merge_yx_complex(work2_c, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#endif

    ! In case of MPI3 shared memory, local master is done writing
    if (d2d_intranode) call decomp_2d_win_transpose_stop_writing(dst_win)

    return
  end subroutine transpose_y_to_x_complex


  subroutine mem_split_yx_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(n1,n2,n3), intent(IN) :: in
    real(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
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
       pos = m * decomp%y1count + 1
#else
       pos = decomp%y1disp(m) + 1
#endif

#if defined(_GPU)
       istat = cudaMemcpy2D( out(pos), n1*(i2-i1+1), in(1,i1,1), n1*n2, n1*(i2-i1+1), n3, cudaMemcpyDeviceToDevice )
#else
       do k=1,n3
          do i=1,n2
             do j=i1,i2
                out(pos + i-1 + n2*(j-i1) + n2*i3*(k-1)) = in(j,i,k)
             end do
          end do
       end do
#endif
    end do

    return
  end subroutine mem_split_yx_real


  subroutine mem_split_yx_complex(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(n1,n2,n3), intent(IN) :: in
    complex(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
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
       pos = m * decomp%y1count + 1
#else
       pos = decomp%y1disp(m) + 1
#endif

#if defined(_GPU)
       istat = cudaMemcpy2D( out(pos), n1*(i2-i1+1), in(1,i1,1), n1*n2, n1*(i2-i1+1), n3, cudaMemcpyDeviceToDevice )
#else
       do k=1,n3
          do i=1,n2
             do j=i1,i2
                out(pos + i-1 + n2*(j-i1) + n2*i3*(k-1)) = in(j,i,k)
             end do
          end do
       end do
#endif
    end do

    return
  end subroutine mem_split_yx_complex


  subroutine mem_merge_yx_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(*), intent(IN) :: in
    real(mytype), dimension(n1,n2,n3), intent(OUT) :: out
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
       pos = m * decomp%x1count + 1
#else
       pos = decomp%x1disp(m) + 1
#endif

#if defined(_GPU)
       istat = cudaMemcpy2D( out(i1,1,1), n1, in(pos), i2-i1+1, i2-i1+1, n2*n3, cudaMemcpyDeviceToDevice )
#else
       do k=1,n3
          do j=1,n2
             do i=i1,i2
                out(i,j,k) = in(pos + i-i1 + i3*(j-1) + i3*n2*(k-1))
             end do
          end do
       end do
#endif
    end do

    return
  end subroutine mem_merge_yx_real


  subroutine mem_merge_yx_complex(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(*), intent(IN) :: in
    complex(mytype), dimension(n1,n2,n3), intent(OUT) :: out
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
       pos = m * decomp%x1count + 1
#else
       pos = decomp%x1disp(m) + 1
#endif

#if defined(_GPU)
       istat = cudaMemcpy2D( out(i1,1,1), n1, in(pos), i2-i1+1, i2-i1+1, n2*n3, cudaMemcpyDeviceToDevice )
#else
       do k=1,n3
          do j=1,n2
             do i=i1,i2
                out(i,j,k) = in(pos + i-i1 + i3*(j-1) + i3*n2*(k-1))
             end do
          end do
       end do
#endif
    end do

    return
  end subroutine mem_merge_yx_complex
