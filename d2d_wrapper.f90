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

submodule(decomp_2d) d2d_wrapper

   implicit none

contains

   !
   ! Fence the given MPI window
   !
   module subroutine decomp_2d_win_fence(win)

      implicit none

      ! Argument
      integer, intent(in) :: win

      ! Local variable
      integer :: ierror

      if (win == MPI_WIN_NULL) return

      call MPI_WIN_FENCE(0, win, ierror)
      if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_WIN_FENCE")

   end subroutine decomp_2d_win_fence

   !
   ! Free the given MPI window
   !
   module subroutine decomp_2d_win_free(win)

      implicit none

      ! Argument
      integer, intent(in) :: win

      ! Local variable
      integer :: ierror

      if (win == MPI_WIN_NULL) return

      call MPI_WIN_FREE(win, ierror)
      if (ierror /= 0) call decomp_2d_abort(ierror, "MPI_WIN_FREE")

   end subroutine decomp_2d_win_free

end submodule d2d_wrapper
