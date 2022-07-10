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

integer :: i, j, k

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
call decomp_2d_win_fence(in_win)
call c2c_1m_x(in_2d,isign,ph,in_win)
else
call c2c_1m_x(in,isign,ph)
endif
#else
if (d2d_intranode) then
call decomp_2d_win_fence(in_win)
call alloc_x(wk1, var2d=wk1_2d, win=wk1_win)
do concurrent (k=1:ph%xsz_loc(3), i=1:ph%xsz_loc(1))
wk1_2d(i,k) = in_2d(i,k)
end do
call c2c_1m_x(wk1_2d,isign,ph,wk1_win)
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
call transpose_x_to_y(in,wk2_c2c,ph)
#else
call transpose_x_to_y(wk1,wk2_c2c,ph)
#endif
if (d2d_intranode) then
call decomp_2d_win_fence(wk2_c2c_win)
call c2c_1m_y(wk2_c2c_2d,isign,ph,wk2_c2c_win)
else
call c2c_1m_y(wk2_c2c,isign,ph)
endif

! ===== Swap Y --> Z; 1D FFTs in Z =====
call transpose_y_to_z(wk2_c2c,out,ph)
if (d2d_intranode) then
call decomp_2d_win_fence(out_win)
call c2c_1m_z(out_2d,isign,ph,out_win)
else
call c2c_1m_z(out,isign,ph)
endif

else if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_BACKWARD &
.OR. & 
format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_FORWARD) then

! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
if (d2d_intranode) then
call decomp_2d_win_fence(in_win)
call c2c_1m_z(in_2d,isign,ph,in_win)
else
call c2c_1m_z(in,isign,ph)
endif
#else
if (d2d_intranode) then
call decomp_2d_win_fence(in_win)
call alloc_z(wk1, var2d=wk1_2d, win=wk1_win)
do concurrent (k=1:ph%zsz_loc(3), i=1:ph%zsz_loc(1))
wk1_2d(i,k) = in_2d(i,k)
end do
call c2c_1m_z(wk1_2d,isign,ph,wk1_win)
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
call transpose_z_to_y(in,wk2_c2c,ph)
#else
call transpose_z_to_y(wk1,wk2_c2c,ph)
#endif
if (d2d_intranode) then
call decomp_2d_win_fence(wk2_c2c_win)
call c2c_1m_y(wk2_c2c_2d,isign,ph,wk2_c2c_win)
else
call c2c_1m_y(wk2_c2c,isign,ph)
endif

! ===== Swap Y --> X; 1D FFTs in X =====
call transpose_y_to_x(wk2_c2c,out,ph)
if (d2d_intranode) then
call decomp_2d_win_fence(out_win)
call c2c_1m_x(out_2d,isign,ph,out_win)
else
call c2c_1m_x(out,isign,ph)
endif

end if

#ifndef OVERWRITE
! Free memory
if (associated(wk1)) nullify(wk1)
if (associated(wk1_2d)) nullify(wk1_2d)
if (d2d_intranode) call decomp_2d_win_free(wk1_win)
#endif

return
end subroutine fft_3d_c2c


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D forward FFT - real to complex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_3d_r2c(in_r, out_c, in_r_2d, out_c_2d, in_r_win, out_c_win)

implicit none

real(mytype), dimension(:,:,:), pointer, intent(INOUT) :: in_r
complex(mytype), dimension(:,:,:), pointer, intent(INOUT) :: out_c
real(mytype), dimension(:,:), pointer, intent(INOUT), optional :: in_r_2d
complex(mytype), dimension(:,:), pointer, intent(INOUT), optional :: out_c_2d
integer, intent(in), optional :: in_r_win, out_c_win

! Safety check
if (d2d_intranode) then
   if (.not.present(in_r_2d) .or. .not.present(out_c_2d) .or. &
       .not.present(in_r_win) .or. .not.present(out_c_win)) &
      call decomp_2d_abort(__FILE__, __LINE__, 0, "Incorrect arguments")
endif

if (format==PHYSICAL_IN_X) then

! ===== 1D FFTs in X =====
if (d2d_intranode) then
call decomp_2d_win_fence(in_r_win)
call r2c_1m_x(in_r_2d,wk13_2d,wk13_win)
else
call r2c_1m_x(in_r,wk13)
endif

! ===== Swap X --> Y; 1D FFTs in Y =====
call transpose_x_to_y(wk13,wk2_r2c,sp)
if (d2d_intranode) then
call decomp_2d_win_fence(wk2_r2c_win)
call c2c_1m_y(wk2_r2c_2d,-1,sp,wk2_r2c_win)
else
call c2c_1m_y(wk2_r2c,-1,sp)
endif

! ===== Swap Y --> Z; 1D FFTs in Z =====
call transpose_y_to_z(wk2_r2c,out_c,sp)
if (d2d_intranode) then
call decomp_2d_win_fence(out_c_win)
call c2c_1m_z(out_c_2d,-1,sp,out_c_win)
else
call c2c_1m_z(out_c,-1,sp)
endif

else if (format==PHYSICAL_IN_Z) then

! ===== 1D FFTs in Z =====
if (d2d_intranode) then
call decomp_2d_win_fence(in_r_win)
call r2c_1m_z(in_r_2d,wk13_2d,wk13_win)
else
call r2c_1m_z(in_r,wk13)
endif

! ===== Swap Z --> Y; 1D FFTs in Y =====
call transpose_z_to_y(wk13,wk2_r2c,sp)
if (d2d_intranode) then
call decomp_2d_win_fence(wk2_r2c_win)
call c2c_1m_y(wk2_r2c_2d,-1,sp,wk2_r2c_win)
else
call c2c_1m_y(wk2_r2c,-1,sp)
endif

! ===== Swap Y --> X; 1D FFTs in X =====
call transpose_y_to_x(wk2_r2c,out_c,sp)
if (d2d_intranode) then
call decomp_2d_win_fence(out_c_win)
call c2c_1m_x(out_c_2d,-1,sp,out_c_win)
else
call c2c_1m_x(out_c,-1,sp)
endif

endif

end subroutine fft_3d_r2c


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D inverse FFT - complex to real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_3d_c2r(in_c, out_r, in_c_2d, out_r_2d, in_c_win, out_r_win)

implicit none

complex(mytype), dimension(:,:,:), pointer, intent(INOUT) :: in_c
real(mytype), dimension(:,:,:), pointer, intent(INOUT) :: out_r
complex(mytype), dimension(:,:), pointer, intent(INOUT), optional :: in_c_2d
real(mytype), dimension(:,:), pointer, intent(INOUT), optional :: out_r_2d
integer, intent(in), optional :: in_c_win, out_r_win

integer :: i, j, k

#ifndef OVERWRITE
complex(mytype), pointer, dimension(:,:,:) :: wk1
complex(mytype), pointer, dimension(:,:) :: wk1_2d
integer :: wk1_win
#endif

if (format==PHYSICAL_IN_X) then

! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
if (d2d_intranode) then
call decomp_2d_win_fence(in_c_win)
call c2c_1m_z(in_c_2d,1,sp,in_c_win)
else
call c2c_1m_z(in_c,1,sp)
endif
#else
if (d2d_intranode) then
call decomp_2d_win_fence(in_c_win)
call alloc_z(wk1, var2d=wk1_2d, opt_decomp=sp, win=wk1_win)
do concurrent (k=1:sp%zsz_loc(3), i=1:sp%zsz_loc(1))
wk1_2d(i,k) = in_c_2d(i,k)
end do
call c2c_1m_z(wk1_2d,1,sp,wk1_win)
else
call alloc_z(wk1, opt_decomp=sp)
do concurrent (k=1:sp%zsz(3), j=1:sp%zsz(2), i=1:sp%zsz(1))
wk1(i,j,k) = in_c(i,j,k)
end do
call c2c_1m_z(wk1,1,sp)
endif
#endif

! ===== Swap Z --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
call transpose_z_to_y(in_c,wk2_r2c,sp)
#else
call transpose_z_to_y(wk1,wk2_r2c,sp)
#endif
if (d2d_intranode) then
call decomp_2d_win_fence(wk2_r2c_win)
call c2c_1m_y(wk2_r2c_2d,1,sp,wk2_r2c_win)
else
call c2c_1m_y(wk2_r2c,1,sp)
endif

! ===== Swap Y --> X; 1D FFTs in X =====
call transpose_y_to_x(wk2_r2c,wk13,sp)
if (d2d_intranode) then
call decomp_2d_win_fence(wk13_win)
call c2r_1m_x(wk13_2d,out_r_2d,out_r_win)
else
call c2r_1m_x(wk13,out_r)
endif

else if (format==PHYSICAL_IN_Z) then

! ===== 1D FFTs in X =====
#ifdef OVERWRITE
if (d2d_intranode) then
call decomp_2d_win_fence(in_c_win)
call c2c_1m_x(in_c_2d,1,sp,in_c_win)
else
call c2c_1m_x(in_c,1,sp)
endif
#else
if (d2d_intranode) then
call decomp_2d_win_fence(in_c_win)
call alloc_x(wk1, var2d=wk1_2d, opt_decomp=sp, win=wk1_win)
do concurrent (k=1:sp%xsz_loc(3), i=1:sp%xsz_loc(1))
wk1_2d(i,k) = in_c_2d(i,k)
end do
call c2c_1m_x(wk1_2d,1,sp,wk1_win)
else
call alloc_x(wk1, opt_decomp=sp)
do concurrent (k=1:sp%xsz(3), j=1:sp%xsz(2), i=1:sp%xsz(1))
wk1(i,j,k) = in_c(i,j,k)
end do
call c2c_1m_x(wk1,1,sp)
endif
#endif

! ===== Swap X --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
call transpose_x_to_y(in_c,wk2_r2c,sp)
#else
call transpose_x_to_y(wk1,wk2_r2c,sp)
#endif
if (d2d_intranode) then
call decomp_2d_win_fence(wk2_r2c_win)
call c2c_1m_y(wk2_r2c_2d,1,sp,wk2_r2c_win)
else
call c2c_1m_y(wk2_r2c,1,sp)
endif

! ===== Swap Y --> Z; 1D FFTs in Z =====
call transpose_y_to_z(wk2_r2c,wk13,sp)
if (d2d_intranode) then
call decomp_2d_win_fence(wk13_win)
call c2r_1m_z(wk13_2d,out_r_2d,out_r_win)
else
call c2r_1m_z(wk13,out_r)
endif

endif

#ifndef OVERWRITE
! Free memory
if (associated(wk1)) nullify(wk1)
if (associated(wk1_2d)) nullify(wk1_2d)
if (d2d_intranode) call decomp_2d_win_free(wk1_win)
#endif

return
end subroutine fft_3d_c2r
