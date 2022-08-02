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

   use, intrinsic :: iso_c_binding, only: c_loc, c_ptr, c_f_pointer

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
      integer :: ierror

      ! Local master broadcast buf_size if needed
      if (d2d_intranode) then
         call MPI_BCAST(buf_size, 1, MPI_INT, 0, DECOMP_2D_LOCALCOMM, ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_BCAST")
      end if

      if (buf_size <= decomp_buf_size) return

      decomp_buf_size = buf_size
      call smod_buffer_free()

#if defined(_GPU)
      call smod_real_array_alloc(work1_r_d, buf_size, work1_r_d_win)
      call smod_cplx_array_alloc(work1_c_d, buf_size, work1_c_d_win)
      call smod_real_array_alloc(work2_r_d, buf_size, work2_r_d_win)
      call smod_cplx_array_alloc(work2_c_d, buf_size, work2_c_d_win)
#endif
      call smod_real_array_alloc(work1_r, buf_size, work1_r_win)
      call smod_cplx_array_alloc(work1_c, buf_size, work1_c_win)
      call smod_real_array_alloc(work2_r, buf_size, work2_r_win)
      call smod_cplx_array_alloc(work2_c, buf_size, work2_c_win)

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
      call smod_buffer_free()
#if defined(_GPU)
#if defined(_NCCL)
      nccl_stat = ncclCommDestroy(nccl_comm_2decomp)
#endif
#endif

   end subroutine decomp_buffer_free

   !
   ! Initialize the window and buffers to null
   !
   module subroutine decomp_buffer_init()

      implicit none

      work1_r_win = MPI_WIN_NULL
      work2_r_win = MPI_WIN_NULL
      work1_c_win = MPI_WIN_NULL
      work2_c_win = MPI_WIN_NULL
      work1_r => null()
      work2_r => null()
      work1_c => null()
      work2_c => null()
#if defined(_GPU)
      work1_r_d_win = MPI_WIN_NULL
      work2_r_d_win = MPI_WIN_NULL
      work1_c_d_win = MPI_WIN_NULL
      work2_c_d_win = MPI_WIN_NULL
      work1_r_d => null()
      work2_r_d => null()
      work1_c_d => null()
      work2_c_d => null()
#endif

   end subroutine decomp_buffer_init

   !
   ! Internal routine to allocate memory
   !
   subroutine smod_real_array_alloc(array, buf_size, win)

      implicit none

      ! Arguments
      real(mytype), dimension(:), pointer, intent(out) :: array
      integer, intent(in) :: buf_size
      integer, intent(out) :: win

      ! Local variables
      integer :: status, errorcode, ierror, tmpdispunit, info
      type(c_ptr) :: baseptr
      integer(kind=MPI_ADDRESS_KIND) :: winsize, tmpsize, tmpptr

      if (d2d_intranode) then
         ! Same disp_unit everywhere
         call MPI_INFO_CREATE(info, ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_INFO_CREATE")
         call MPI_INFO_SET(info, "same_disp_unit", "true", ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_INFO_SET")
         call MPI_INFO_SET(info, "alloc_shared_noncontig", "true", ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_INFO_SET")

         ! Allocate shared memory on local master
         if (nrank_loc == 0) then
            winsize = buf_size
         else
            winsize = 0
         end if
         call MPI_WIN_ALLOCATE_SHARED(winsize*mytype_bytes, &
                                      mytype_bytes, &
                                      info, &
                                      DECOMP_2D_LOCALCOMM, &
                                      baseptr, &
                                      win, &
                                      ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_ALLOCATE_SHARED")

         ! Get the memory on local master
         call MPI_WIN_SHARED_QUERY(win, 0, tmpsize, tmpdispunit, tmpptr, ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_SHARED_QUERY")
         baseptr = transfer(tmpptr, baseptr)
         call C_F_POINTER(baseptr, array, (/buf_size/))

         ! info is no longer needed
         call MPI_INFO_FREE(info, ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_INFO_FREE")

      else
         ! Allocate memory
         allocate (array(buf_size), STAT=status)
         if (status /= 0) then
            errorcode = 2
            call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                                 'Out of memory when allocating 2DECOMP workspace')
         end if
      end if

   end subroutine smod_real_array_alloc
   subroutine smod_cplx_array_alloc(array, buf_size, win)

      implicit none

      ! Arguments
      complex(mytype), dimension(:), pointer, intent(out) :: array
      integer, intent(in) :: buf_size
      integer, intent(out) :: win

      ! Local variables
      integer :: status, errorcode, ierror, tmpdispunit, info
      type(c_ptr) :: baseptr
      integer(kind=MPI_ADDRESS_KIND) :: winsize, tmpsize, tmpptr

      if (d2d_intranode) then
         ! Same disp_unit everywhere
         call MPI_INFO_CREATE(info, ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_INFO_CREATE")
         call MPI_INFO_SET(info, "same_disp_unit", "true", ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_INFO_SET")
         call MPI_INFO_SET(info, "alloc_shared_noncontig", "true", ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_INFO_SET")

         ! Allocate shared memory on local master
         if (nrank_loc == 0) then
            winsize = buf_size
         else
            winsize = 0
         end if
         call MPI_WIN_ALLOCATE_SHARED(winsize*mytype_bytes*2, &
                                      mytype_bytes*2, &
                                      info, &
                                      DECOMP_2D_LOCALCOMM, &
                                      baseptr, &
                                      win, &
                                      ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_ALLOCATE_SHARED")

         ! Get the memory on local master
         call MPI_WIN_SHARED_QUERY(win, 0, tmpsize, tmpdispunit, tmpptr, ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_SHARED_QUERY")
         baseptr = transfer(tmpptr, baseptr)
         call C_F_POINTER(baseptr, array, (/buf_size/))

         ! info is no longer needed
         call MPI_INFO_FREE(info, ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_INFO_FREE")

      else
         ! Allocate memory
         allocate (array(buf_size), STAT=status)
         if (status /= 0) then
            errorcode = 2
            call decomp_2d_abort(__FILE__, __LINE__, errorcode, &
                                 'Out of memory when allocating 2DECOMP workspace')
         end if
      end if

   end subroutine smod_cplx_array_alloc

   !
   ! Internal routine to free the buffers
   !
   subroutine smod_buffer_free

      implicit none

      call decomp_2d_win_free(work1_r_win)
      call decomp_2d_win_free(work2_r_win)
      call decomp_2d_win_free(work1_c_win)
      call decomp_2d_win_free(work2_c_win)
      if (d2d_intranode) then
         if (associated(work1_r)) nullify (work1_r)
         if (associated(work2_r)) nullify (work2_r)
         if (associated(work1_c)) nullify (work1_c)
         if (associated(work2_c)) nullify (work2_c)
      else
         if (associated(work1_r)) deallocate (work1_r)
         if (associated(work2_r)) deallocate (work2_r)
         if (associated(work1_c)) deallocate (work1_c)
         if (associated(work2_c)) deallocate (work2_c)
      end if
#if defined(_GPU)
      call decomp_2d_win_free(work1_r_d_win)
      call decomp_2d_win_free(work2_r_d_win)
      call decomp_2d_win_free(work1_c_d_win)
      call decomp_2d_win_free(work2_c_d_win)
      if (d2d_intranode) then
         if (associated(work1_r_d)) nullify (work1_r_d)
         if (associated(work2_r_d)) nullify (work2_r_d)
         if (associated(work1_c_d)) nullify (work1_c_d)
         if (associated(work2_c_d)) nullify (work2_c_d)
      else
         if (associated(work1_r_d)) deallocate (work1_r_d)
         if (associated(work2_r_d)) deallocate (work2_r_d)
         if (associated(work1_c_d)) deallocate (work1_c_d)
         if (associated(work2_c_d)) deallocate (work2_c_d)
      end if
#endif

   end subroutine smod_buffer_free

end submodule smod_buffer
