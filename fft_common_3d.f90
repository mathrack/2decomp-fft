!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2011 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This file contains 3D c2c/r2c/c2r transform subroutines which are
! identical for several FFT engines 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D FFT - complex to complex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_3d_c2c(in, out, isign, in_2d, out_2d, in_win, out_win)

implicit none

complex(mytype), dimension(:,:,:), pointer, intent(INOUT) :: in
complex(mytype), dimension(:,:,:), pointer, intent(INOUT) :: out
integer, intent(IN) :: isign
complex(mytype), dimension(:,:), pointer, intent(INOUT), optional :: in_2d
complex(mytype), dimension(:,:), pointer, intent(INOUT), optional :: out_2d
integer, intent(in), optional :: in_win, out_win

integer :: i, j, k, ierror

#ifndef OVERWRITE
integer :: wk1_win
complex(mytype), pointer, dimension(:,:,:) :: wk1
complex(mytype), pointer, dimension(:,:) :: wk1_2d
#endif

! Safety check
if (d2d_intranode) then
   if (.not.present(in_2d) .or. .not.present(out_2d) .or. &
       .not.present(in_win) .or. .not.present(out_win)) &
      call decomp_2d_abort(__FILE__, __LINE__, 0, "Incorrect arguments")
endif

if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_FORWARD .OR.  &
format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_BACKWARD) then

! ===== 1D FFTs in X =====
#ifdef OVERWRITE
if (d2d_intranode) then
call c2c_1m_x(in_2d,isign,ph)
else
call c2c_1m_x(in,isign,ph)
endif
#else
if (d2d_intranode) then
call alloc_x(wk1, var2d=wk1_2d, win=wk1_win)
do concurrent (k=1:ph%xsz_loc(3), i=1:ph%xsz_loc(1))
wk1_2d(i,k) = in_2d(i,k)
end do
call c2c_1m_x(wk1_2d,isign,ph)
else
call alloc_x(wk1)
do concurrent (k=1:ph%xsz(3), j=1:ph%xsz(2), i=1:ph%xsz(1))
wk1(i,j,k) = in(i,j,k)
end do
call c2c_1m_x(wk1,isign,ph)
endif
#endif

! ===== Swap X --> Y; 1D FFTs in Y =====

#ifdef OVERWRITE
if (d2d_intranode) then
call MPI_WIN_FENCE(0, in_win, ierror)
if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_FENCE")
endif
call transpose_x_to_y(in,wk2_c2c,ph)
#else
if (d2d_intranode) then
call MPI_WIN_FENCE(0, wk1_win, ierror)
if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_FENCE")                   
endif
call transpose_x_to_y(wk1,wk2_c2c,ph)
#endif
if (d2d_intranode) then
call MPI_WIN_FENCE(0, wk2_c2c_win, ierror)
if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_FENCE")
call c2c_1m_y(wk2_c2c_2d,isign,ph)
else
call c2c_1m_y(wk2_c2c,isign,ph)
endif

! ===== Swap Y --> Z; 1D FFTs in Z =====
if (d2d_intranode) then
call MPI_WIN_FENCE(0, wk2_c2c_win, ierror)
if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_FENCE")
endif
call transpose_y_to_z(wk2_c2c,out,ph)
if (d2d_intranode) then
call MPI_WIN_FENCE(0, out_win, ierror)
if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_FENCE")
call c2c_1m_z(out_2d,isign,ph)
else
call c2c_1m_z(out,isign,ph)
endif

else if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_BACKWARD &
.OR. & 
format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_FORWARD) then

! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
if (d2d_intranode) then
call c2c_1m_z(in_2d,isign,ph)
else
call c2c_1m_z(in,isign,ph)
endif
#else
if (d2d_intranode) then
call alloc_z(wk1, var2d=wk1_2d, win=wk1_win)
do concurrent (k=1:ph%zsz_loc(3), i=1:ph%zsz_loc(1))
wk1_2d(i,k) = in_2d(i,k)
end do
call c2c_1m_z(wk1_2d,isign,ph)
else
call alloc_z(wk1)
do concurrent (k=1:ph%zsz(3), j=1:ph%zsz(2), i=1:ph%zsz(1))
wk1(i,j,k) = in(i,j,k)
end do
call c2c_1m_z(wk1,isign,ph)
endif
#endif

! ===== Swap Z --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
if (d2d_intranode) then
call MPI_WIN_FENCE(0, in_win, ierror)
if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_FENCE")
endif
call transpose_z_to_y(in,wk2_c2c,ph)
#else
if (d2d_intranode) then
call MPI_WIN_FENCE(0, wk1_win, ierror)
if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_FENCE")
endif
call transpose_z_to_y(wk1,wk2_c2c,ph)
#endif
if (d2d_intranode) then
call MPI_WIN_FENCE(0, wk2_c2c_win, ierror)
if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_FENCE")
call c2c_1m_y(wk2_c2c_2d,isign,ph)
else
call c2c_1m_y(wk2_c2c,isign,ph)
endif

! ===== Swap Y --> X; 1D FFTs in X =====
if (d2d_intranode) then
call MPI_WIN_FENCE(0, wk2_c2c_win, ierror)
if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_FENCE")
endif
call transpose_y_to_x(wk2_c2c,out,ph)
if (d2d_intranode) then
call MPI_WIN_FENCE(0, out_win, ierror)
if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_FENCE")
call c2c_1m_x(out_2d,isign,ph)
else
call c2c_1m_x(out,isign,ph)
endif

end if

#ifndef OVERWRITE
! Free memory
if (associated(wk1)) nullify(wk1)
if (associated(wk1_2d)) nullify(wk1_2d)
if (d2d_intranode) then
   call MPI_WIN_FREE(wk1_win, ierror)
   if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_WIN_FREE")
endif
#endif

return
end subroutine fft_3d_c2c


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D forward FFT - real to complex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_3d_r2c(in_r, out_c)

implicit none

real(mytype), dimension(:,:,:), pointer, intent(IN) :: in_r
complex(mytype), dimension(:,:,:), pointer, intent(OUT) :: out_c

if (format==PHYSICAL_IN_X) then

! ===== 1D FFTs in X =====
call r2c_1m_x(in_r,wk13)

! ===== Swap X --> Y; 1D FFTs in Y =====
call transpose_x_to_y(wk13,wk2_r2c,sp)
call c2c_1m_y(wk2_r2c,-1,sp)

! ===== Swap Y --> Z; 1D FFTs in Z =====
call transpose_y_to_z(wk2_r2c,out_c,sp)
call c2c_1m_z(out_c,-1,sp)

else if (format==PHYSICAL_IN_Z) then

! ===== 1D FFTs in Z =====
call r2c_1m_z(in_r,wk13)

! ===== Swap Z --> Y; 1D FFTs in Y =====
call transpose_z_to_y(wk13,wk2_r2c,sp)
call c2c_1m_y(wk2_r2c,-1,sp)

! ===== Swap Y --> X; 1D FFTs in X =====
call transpose_y_to_x(wk2_r2c,out_c,sp)
call c2c_1m_x(out_c,-1,sp)

end if

return
end subroutine fft_3d_r2c


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D inverse FFT - complex to real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_3d_c2r(in_c, out_r)

implicit none

complex(mytype), dimension(:,:,:), pointer, intent(INOUT) :: in_c
real(mytype), dimension(:,:,:), pointer, intent(INOUT) :: out_r

integer :: i, j, k

#ifndef OVERWRITE
complex(mytype), pointer, dimension(:,:,:) :: wk1
#endif

if (format==PHYSICAL_IN_X) then

! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
call c2c_1m_z(in_c,1,sp)       
#else
allocate(wk1(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
do concurrent (k=1:sp%zsz(3), j=1:sp%zsz(2), i=1:sp%zsz(1))
wk1(i,j,k) = in_c(i,j,k)
end do
call c2c_1m_z(wk1,1,sp)
#endif

! ===== Swap Z --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
call transpose_z_to_y(in_c,wk2_r2c,sp)
#else
call transpose_z_to_y(wk1,wk2_r2c,sp)
#endif
call c2c_1m_y(wk2_r2c,1,sp)

! ===== Swap Y --> X; 1D FFTs in X =====
call transpose_y_to_x(wk2_r2c,wk13,sp)
call c2r_1m_x(wk13,out_r)

else if (format==PHYSICAL_IN_Z) then

! ===== 1D FFTs in X =====
#ifdef OVERWRITE
call c2c_1m_x(in_c,1,sp)
#else
allocate(wk1(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
do concurrent (k=1:sp%xsz(3), j=1:sp%xsz(2), i=1:sp%xsz(1))
wk1(i,j,k) = in_c(i,j,k)
end do
call c2c_1m_x(wk1,1,sp)
#endif

! ===== Swap X --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
call transpose_x_to_y(in_c,wk2_r2c,sp)
#else
call transpose_x_to_y(wk1,wk2_r2c,sp)
#endif
call c2c_1m_y(wk2_r2c,1,sp)

! ===== Swap Y --> Z; 1D FFTs in Z =====
call transpose_y_to_z(wk2_r2c,wk13,sp)
call c2r_1m_z(wk13,out_r)

end if

#ifndef OVERWRITE
! Free memory
if (associated(wk1)) nullify(wk1)
#endif

return
end subroutine fft_3d_c2r
