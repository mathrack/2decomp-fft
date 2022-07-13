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

! This is the submodule dedicated to log subroutines

submodule(decomp_2d) smod_buffer

   implicit none

contains

   !
   ! Free and allocate the buffers if needed
   !
   module subroutine decomp_buffer_alloc(buf_size)

      implicit none

      ! Argument
      integer, intent(inout) :: buf_size

      ! Local variables
      integer :: ierror, status, errorcode

      ! Local master broadcast buf_size if needed
      if (d2d_intranode) then
         call MPI_BCAST(buf_size, 1, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
      end if

      if (buf_size <= decomp_buf_size) return

      decomp_buf_size = buf_size
      call _buffer_free()

#if defined(_GPU)
      allocate (work1_r_d(buf_size), STAT=status)
      if (status /= 0) then
         errorcode = 2
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                              'Out of memory when allocating 2DECOMP workspace')
      end if
      allocate (work1_c_d(buf_size), STAT=status)
      if (status /= 0) then
         errorcode = 2
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                              'Out of memory when allocating 2DECOMP workspace')
      end if
      allocate (work2_r_d(buf_size), STAT=status)
      if (status /= 0) then
         errorcode = 2
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                              'Out of memory when allocating 2DECOMP workspace')
      end if
      allocate (work2_c_d(buf_size), STAT=status)
      if (status /= 0) then
         errorcode = 2
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                              'Out of memory when allocating 2DECOMP workspace')
      end if
#endif
      allocate (work1_r(buf_size), STAT=status)
      if (status /= 0) then
         errorcode = 2
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                              'Out of memory when allocating 2DECOMP workspace')
      end if
      allocate (work2_r(buf_size), STAT=status)
      if (status /= 0) then
         errorcode = 2
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                              'Out of memory when allocating 2DECOMP workspace')
      end if
      allocate (work1_c(buf_size), STAT=status)
      if (status /= 0) then
         errorcode = 2
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                              'Out of memory when allocating 2DECOMP workspace')
      end if
      allocate (work2_c(buf_size), STAT=status)
      if (status /= 0) then
         errorcode = 2
         call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                              'Out of memory when allocating 2DECOMP workspace')
      end if

   end subroutine decomp_buffer_alloc

   !
   ! Public routine to free the buffers
   !
   module subroutine decomp_buffer_free

      implicit none

      ! Local variable
#if defined(_GPU)
#if defined(_NCCL)
      integer :: nccl_stat
#endif
#endif

      decomp_buf_size = 0
      call _buffer_free()
#if defined(_GPU)
#if defined(_NCCL)
      nccl_stat = ncclCommDestroy(nccl_comm_2decomp)
#endif
#endif

   end subroutine decomp_buffer_free

   !
   ! Internal routine to free the buffers
   !
   module subroutine _buffer_free

      implicit none

      if (allocated(work1_r)) deallocate (work1_r)
      if (allocated(work2_r)) deallocate (work2_r)
      if (allocated(work1_c)) deallocate (work1_c)
      if (allocated(work2_c)) deallocate (work2_c)
#if defined(_GPU)
      if (allocated(work1_r_d)) deallocate (work1_r_d)
      if (allocated(work2_r_d)) deallocate (work2_r_d)
      if (allocated(work1_c_d)) deallocate (work1_c_d)
      if (allocated(work2_c_d)) deallocate (work2_c_d)
#endif

   end subroutine _buffer_free

end submodule smod_buffer
