subroutine boot_xcompact3d(hostcomm, localcomm)

   use MPI
   use decomp_2d, only: nproc, decomp_2d_abort

   implicit none

   integer, intent(inout) :: hostcomm, localcomm

   integer :: localrank, key, code

   ! Initialise MPI
   call MPI_INIT(code)
   if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_INIT")

   ! Split MPI_COMM_WORLD ?
   if (.true.) then
      write(*,*) "SPLIT"
      ! Split using shared memory
      call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, localcomm, code)
      if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_COMM_SPLIT_TYPE")
      ! Each CPU has a rank in the local communicator
      call MPI_COMM_RANK(localcomm, localrank, code)                                                         
      if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_COMM_RANK") 
      ! Use the global rank for the key
      call MPI_COMM_RANK(MPI_COMM_WORLD, key, code)
      if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_COMM_RANK")
      ! Create a new communicator using all masters
      if (localrank == 0) then
         call MPI_COMM_SPLIT(MPI_COMM_WORLD, 0, key, hostcomm, code)
      else
         call MPI_COMM_SPLIT(MPI_COMM_WORLD, MPI_UNDEFINED, key, hostcomm, code)
      endif
      if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_COMM_SPLIT")
   else
      write(*,*) "NOSPLIT"
      hostcomm = MPI_COMM_WORLD
      localcomm = MPI_COMM_NULL
   endif

end subroutine boot_xcompact3d

subroutine init_xcompact3d(hostcomm, localcomm)

   use decomp_2d, only: decomp_2d_init
   use decomp_2d_fft, only : PHYSICAL_IN_X, decomp_2d_fft_init

   implicit none

   integer, intent(in) :: hostcomm, localcomm

   integer :: nx, ny, nz, p_row, p_col

   nx = 2
   ny = 3
   nz = 4
   p_row = 0
   p_col = 0

   call decomp_2d_init(nx, ny, nz, p_row, p_col, &
           glob_comm = hostcomm, &
           local_comm = localcomm)
   call decomp_2d_fft_init(PHYSICAL_IN_X)

end subroutine init_xcompact3d

subroutine finalise_xcompact3d()

   use MPI
   use decomp_2d, only: decomp_2d_finalize
   use decomp_2d_fft, only : decomp_2d_fft_finalize

   implicit none

   integer :: ierr

   call decomp_2d_fft_finalize()
   call decomp_2d_finalize()
   CALL MPI_FINALIZE(ierr)

end subroutine finalise_xcompact3d

program xcompact3d

   use MPI
   use decomp_2d
   use decomp_2d_fft

   implicit none

   ! inter-node and intra-node communicators
   integer :: hostcomm, localcomm

   ! local variables
   integer :: i, j, k, pos, code

   type(decomp_data) :: ux1, ref1, t1, ux3, ref3, t3

   type(decomp_data) :: cux1, cref1, ct1, cux3, cref3, ct3, wk

   call boot_xcompact3d(hostcomm, localcomm)

   call init_xcompact3d(hostcomm, localcomm)

   call ux1%init(is_cplx=.false., idir=1)
   call ref1%init(is_cplx=.false., idir=1, contig=.true.)
   call t1%copy(ux1)

   call ux3%init(is_cplx=.false., idir=3)
   call ref3%init(is_cplx=.false., idir=3, contig=.true.)
   call t3%copy(ux3)

   call wk%init(is_cplx=.true., idir=3, decomp=sp)

   call cux1%init(is_cplx=.true., idir=1, decomp=ph)
   call cref1%init(is_cplx=.true., idir=1, decomp=ph, contig=.true.)
   call ct1%copy(cux1)

   call cux3%init(is_cplx=.true., idir=3, decomp=ph)
   call cref3%init(is_cplx=.true., idir=3, decomp=ph, contig=.true.)
   call ct3%copy(cux3)

   call decomp_2d_win_fence(ref1%win)
   call decomp_2d_win_fence(cref1%win)
   if (nrank_loc <= 0) then
      do k = 1, decomp_main%xsz(3)
      do j = 1, decomp_main%xsz(2)
      do i = 1, decomp_main%xsz(1)
         pos = i-1 + decomp_main%xst(1)-1 &
             + nx_global * (j-1 + decomp_main%xst(2)-1) &
             + nx_global * ny_global * (k-1 + decomp_main%xst(3)-1)
         ref1%var(i,j,k) = pos
         cref1%cvar(i,j,k) = cmplx(pos, -pos)
      enddo
      enddo
      enddo
   endif
   call decomp_2d_win_fence(ref1%win)
   call decomp_2d_win_fence(cref1%win)

   call decomp_2d_win_fence(ref3%win)
   call decomp_2d_win_fence(cref3%win)
   if (nrank_loc <= 0) then
      do k = 1, decomp_main%zsz(3)
      do j = 1, decomp_main%zsz(2)
      do i = 1, decomp_main%zsz(1)
         pos = nx_global * ny_global * (i-1 + decomp_main%zst(1)-1) &
             + j-1 + decomp_main%zst(2)-1 &
             + nx_global * (k-1 + decomp_main%zst(3)-1)
         ref3%var(i,j,k) = pos
         cref3%cvar(i,j,k) = cmplx(pos, -pos)
      enddo
      enddo
      enddo 
   endif
   call decomp_2d_win_fence(ref3%win)
   call decomp_2d_win_fence(cref3%win)

   ! Test r2c => c2r
   call decomp_2d_fft_3d(ref1, wk)
   wk%cvar2d = wk%cvar2d / real(xsize(1), kind=mytype) / real(ysize(1), kind=mytype) / real(zsize(1), kind=mytype)
   call decomp_2d_fft_3d(wk, ux1)
   t1%var2d = ref1%var2d - ux1%var2d
   print *, "z => x => z ?", maxval(abs(t1%var2d))

   ! Test c2c => c2c
   cux1%cvar2d = cref1%cvar2d
   call decomp_2d_fft_3d(cux1, cux3, DECOMP_2D_FFT_FORWARD)
   cux3%cvar2d = cux3%cvar2d / real(xsize(1), kind=mytype) / real(ysize(1), kind=mytype) / real(zsize(1), kind=mytype)
   call decomp_2d_fft_3d(cux3, cux1, DECOMP_2D_FFT_BACKWARD)
   ct1%cvar2d = cref1%cvar2d - cux1%cvar2d
   print *, "z =c=> x =c=> z ?", maxval(abs(real(ct1%cvar2d))), maxval(abs(aimag(ct1%cvar2d)))

   call ux1%fin
   call ref1%fin
   call t1%fin

   call ux3%fin
   call ref3%fin
   call t3%fin

   call wk%fin

   call cux1%fin
   call cref1%fin
   call ct1%fin

   call cux3%fin
   call cref3%fin
   call ct3%fin

   call finalise_xcompact3d()

end program xcompact3d
