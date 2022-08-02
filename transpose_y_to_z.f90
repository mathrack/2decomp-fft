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

! This file contains the routines that transpose data from Y to Z pencil

  subroutine transpose_y_to_z_data(src, dst, no_prepro)

    implicit none

    ! Arguments
    type(decomp_data), intent(in) :: src
    type(decomp_data), intent(inout) :: dst
    logical, optional :: no_prepro

    ! Local variables
    logical :: prepro
    integer :: s1, s2, s3, d1, d2, d3, ierror
#if defined(_GPU)
#if defined(_NCCL)
    integer :: row_rank_id
#endif 
#endif

    ! Safety check
#ifdef DEBUG
    if (src%idir /= 2 .or. dst%idir /= 3 .or. &
        (src%is_cplx.neqv.dst%is_cplx) .or. &
        (src%shm.neqv.dst%shm) .or. &
        .not.associated(src%decomp, dst%decomp)) &
       call decomp_2d_abort(__FILE__, __LINE__, -1, "Impossible transpose operation")
#endif

    ! Default value : prepro is needed
    if (.not.present(no_prepro)) then
       prepro = .true.
    else
       prepro = .not.no_prepro
    endif

    s1 = src%decomp%ysz(1)
    s2 = src%decomp%ysz(2)
    s3 = src%decomp%ysz(3)
    d1 = src%decomp%zsz(1)
    d2 = src%decomp%zsz(2)
    d3 = src%decomp%zsz(3)

    ! In case of MPI3 shared memory
    if (d2d_intranode .and. prepro) then
       call decomp_2d_win_transpose_start_reading(src%win)
       if (src%is_cplx) then
          call decomp_2d_win_transpose_start_writing(work1_c_win)
       else
          call decomp_2d_win_transpose_start_writing(work1_r_win)
       endif
    endif

    ! rearrange source array as send buffer
    if (prepro .or. (.not.d2d_intranode)) then
    if (src%is_cplx) then
#if defined(_GPU)
       call mem_split_yz_complex(src%cvar, s1, s2, s3, work1_c_d, dims(2), &
            src%decomp%y2dist, src%decomp)
#else
       call mem_split_yz_complex(src, s1, s2, s3, work1_c, dims(2), &
            src%decomp%y2dist, src%decomp)
#endif
    else
#if defined(_GPU)
       call mem_split_yz_real(src%var, s1, s2, s3, work1_r_d, dims(2), &
            src%decomp%y2dist, src%decomp)
#else
       call mem_split_yz_real(src, s1, s2, s3, work1_r, dims(2), &
            src%decomp%y2dist, src%decomp)
#endif
    endif
    endif

    ! In case of MPI3 shared memory
    if (d2d_intranode) then
       if (prepro) call decomp_2d_win_transpose_stop_reading(src%win)
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
       call MPI_ALLTOALL(work1_c, src%decomp%y2count, &
            complex_type, work2_c, src%decomp%z2count, &
            complex_type, DECOMP_2D_COMM_ROW, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
#else

#if defined(_GPU)
       call MPI_ALLTOALLV(work1_c_d, src%decomp%y2cnts, src%decomp%y2disp, &
            complex_type, work2_c_d, src%decomp%z2cnts, src%decomp%z2disp, &
            complex_type, DECOMP_2D_COMM_ROW, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#else
       call MPI_ALLTOALLV(work1_c, src%decomp%y2cnts, src%decomp%y2disp, &
            complex_type, work2_c, src%decomp%z2cnts, src%decomp%z2disp, &
            complex_type, DECOMP_2D_COMM_ROW, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif

#endif
    else
#ifdef EVEN
       call MPI_ALLTOALL(work1_r, src%decomp%y2count, &
            real_type, work2_r, src%decomp%z2count, &
            real_type, DECOMP_2D_COMM_ROW, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALL")
#else

#if defined(_GPU)
#if defined(_NCCL)
       nccl_stat = ncclGroupStart()
       do row_rank_id = 0, (row_comm_size - 1)
           nccl_stat = ncclSend(work1_r_d( src%decomp%y2disp(row_rank_id)+1 ), src%decomp%y2cnts(row_rank_id), ncclDouble, local_to_global_row(row_rank_id+1), nccl_comm_2decomp, cuda_stream_2decomp)
           nccl_stat = ncclRecv(work2_r_d( src%decomp%z2disp(row_rank_id)+1 ), src%decomp%z2cnts(row_rank_id), ncclDouble, local_to_global_row(row_rank_id+1), nccl_comm_2decomp, cuda_stream_2decomp)
       end do
       nccl_stat = ncclGroupEnd()
       cuda_stat = cudaStreamSynchronize(cuda_stream_2decomp)
#else
       call MPI_ALLTOALLV(work1_r_d, src%decomp%y2cnts, src%decomp%y2disp, &
            real_type, work2_r_d, src%decomp%z2cnts, src%decomp%z2disp, &
            real_type, DECOMP_2D_COMM_ROW, ierror)
       if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_ALLTOALLV")
#endif
#else
       call MPI_ALLTOALLV(work1_r, src%decomp%y2cnts, src%decomp%y2disp, &
            real_type, work2_r, src%decomp%z2cnts, src%decomp%z2disp, &
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
#if defined(_GPU)
       istat = cudaMemcpy( dst%cvar, work2_c_d, d1*d2*d3, cudaMemcpyDeviceToDevice )
#else
       call mem_merge_yz_complex(work2_c, d1, d2, d3, dst, dims(2), &
         src%decomp%z2dist, src%decomp)
#endif
    else
#if defined(_GPU)
       istat = cudaMemcpy( dst%var, work2_r_d, d1*d2*d3, cudaMemcpyDeviceToDevice )
#else
       call mem_merge_yz_real(work2_r, d1, d2, d3, dst, dims(2), &
         src%decomp%z2dist, src%decomp)
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

  end subroutine transpose_y_to_z_data

  subroutine mem_split_yz_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    type(decomp_data), intent(in) :: in
    real(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
    attributes(device) :: out
    integer :: istat
#endif

    integer :: i,j,k, m,i1,i2,i3, pos

    if (d2d_intranode) then
       do k = 1, decomp%ysz_loc(3)
          do j = 1, decomp%ysz_loc(1)
             out(decomp%intramap_split(j,k,3)) = in%var2d(j,k)
          enddo
       enddo
       return
    endif

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
       istat = cudaMemcpy2D( out(pos), n1*(i2-i1+1), in(1,i1,1), n1*n2, n1*(i2-i1+1), n3, cudaMemcpyDeviceToDevice )
#else
       do k=1,n3
          do i=1,n2
             do j=i1,i2
                out(pos + i-1 + n2*(j-i1) + n2*i3*(k-1)) = in%var(j,i,k)
             end do
          end do
       end do
#endif
    end do

    return
  end subroutine mem_split_yz_real


  subroutine mem_split_yz_complex(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    type(decomp_data), intent(in) :: in
    complex(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
#if defined(_GPU)
    attributes(device) :: out
    integer :: istat
#endif

    integer :: i,j,k, m,i1,i2,i3, pos

    if (d2d_intranode) then
       do k = 1, decomp%ysz_loc(3)
          do j = 1, decomp%ysz_loc(1)
             out(decomp%intramap_split(j,k,3)) = in%cvar2d(j,k)
          enddo
       enddo
       return
    endif

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
       istat = cudaMemcpy2D( out(pos), n1*(i2-i1+1), in(1,i1,1), n1*n2, n1*(i2-i1+1), n3, cudaMemcpyDeviceToDevice )
#else
       do k=1,n3
          do i=1,n2
             do j=i1,i2
                out(pos + i-1 + n2*(j-i1) + n2*i3*(k-1)) = in%cvar(j,i,k)
             end do
          end do
       end do
#endif
    end do

    return
  end subroutine mem_split_yz_complex


  subroutine mem_merge_yz_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(*), intent(IN) :: in
    type(decomp_data), intent(inout) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2, pos

    if (d2d_intranode) then
       do j = 1, decomp%zsz_loc(3)
          do k = 1, decomp%zsz_loc(1)
             out%var2d(k,j) = in(decomp%intramap_merge(k,j,3))
          enddo
       enddo
       return
    endif

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
       do j=1,n3
          do i=1,n2
             do k=i1,i2
                out%var(k,i,j) = in(pos + i-1 + n2*(j-1) + n2*n3*(k-i1))
             end do
          end do
       end do
    end do

    return
  end subroutine mem_merge_yz_real


  subroutine mem_merge_yz_complex(in,n1,n2,n3,out,iproc,dist,decomp)
    
    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(*), intent(IN) :: in
    type(decomp_data), intent(inout) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2, pos

    if (d2d_intranode) then
       do j = 1, decomp%zsz_loc(3)
          do k = 1, decomp%zsz_loc(1)
             out%cvar2d(k,j) = in(decomp%intramap_merge(k,j,3))
          enddo
       enddo
       return
    endif

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
       do j=1,n3
          do i=1,n2
             do k=i1,i2
                out%cvar(k,i,j) = in(pos + i-1 + n2*(j-1) + n2*n3*(k-i1))
             end do
          end do
       end do
    end do

    return
  end subroutine mem_merge_yz_complex
