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

! This is the submodule dedicated to wrappers

submodule(decomp_2d) smod_wrapper

   implicit none

   logical, parameter :: transpose_win_sync_with_fence = .true.

contains

   !
   ! Windows sync. at the beginning of the transpose operation, reading part
   !
   module subroutine decomp_2d_win_transpose_start_reading(src_win)

      implicit none

      integer, intent(in) :: src_win

      if (transpose_win_sync_with_fence) then
         ! Basic sync. using mpi_win_fence
         call decomp_2d_win_fence(src_win, MPI_MODE_NOPUT+MPI_MODE_NOSTORE)
      else
         ! Advanced sync. using post-start-complete-wait
      endif

   end subroutine decomp_2d_win_transpose_start_reading

   !
   ! Windows sync. during the transpose operation, reading part
   !
   module subroutine decomp_2d_win_transpose_stop_reading(src_win)

      implicit none

      integer, intent(in) :: src_win

      if (transpose_win_sync_with_fence) then
         ! Basic sync. using mpi_win_fence
         call decomp_2d_win_fence(src_win, MPI_MODE_NOPUT+MPI_MODE_NOSTORE)
      else
         ! Advanced sync. using post-start-complete-wait
      endif

   end subroutine decomp_2d_win_transpose_stop_reading

   !
   ! Windows sync. during the transpose operation, writing part
   !
   module subroutine decomp_2d_win_transpose_start_writing(dst_win)

      implicit none

      integer, intent(in) :: dst_win

      if (transpose_win_sync_with_fence) then
         ! Basic sync. using mpi_win_fence
         call decomp_2d_win_fence(dst_win)
      else
         ! Advanced sync. using post-start-complete-wait
      endif

   end subroutine decomp_2d_win_transpose_start_writing

   !
   ! Windows sync. at the end of the transpose operation, writing part
   !
   module subroutine decomp_2d_win_transpose_stop_writing(dst_win)

      implicit none

      integer, intent(in) :: dst_win

      if (transpose_win_sync_with_fence) then
         ! Basic sync. using mpi_win_fence
         call decomp_2d_win_fence(dst_win)
      else
         ! Advanced sync. using post-start-complete-wait
      endif

   end subroutine decomp_2d_win_transpose_stop_writing

   !
   ! Fence the given MPI window
   !
   module subroutine decomp_2d_win_fence(win, assert)

      implicit none

      ! Argument
      integer, intent(in) :: win
      integer, intent(in), optional :: assert

      ! Local variable
      integer :: ierror

      if (.not.d2d_intranode) return

      ! Safety check only in debug mode
#ifdef DEBUG
      if (win == MPI_WIN_NULL) call decomp_2d_abort(-1,"Error")
#endif

      if (present(assert)) then
         call MPI_WIN_FENCE(assert, win, ierror)
      else
         call MPI_WIN_FENCE(0, win, ierror)
      endif
      if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_WIN_FENCE")

   end subroutine decomp_2d_win_fence

   !
   ! Free the given MPI window
   !
   module subroutine decomp_2d_win_free(win)

      implicit none

      ! Argument
      integer, intent(inout) :: win

      ! Local variable
      integer :: ierror

      ! Safety check only in debug mode
#ifdef DEBUG
      if (.not.d2d_intranode) call decomp_2d_abort(-1,"Error")
      if (win == MPI_WIN_NULL) call decomp_2d_abort(-1,"Error")
#endif

      call MPI_WIN_FREE(win, ierror)
      if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_WIN_FREE")

   end subroutine decomp_2d_win_free

end submodule smod_wrapper
