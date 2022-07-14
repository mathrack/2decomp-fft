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

subroutine fft_3d_c2c(in, out, isign)

implicit none

! Arguments
type(decomp_data), intent(INOUT) :: in
type(decomp_data), intent(INOUT) :: out
integer, intent(IN) :: isign

! Local variables
integer :: i, j, k

#ifndef OVERWRITE
type(decomp_data) :: wk1
#endif

if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_FORWARD .OR.  &
format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_BACKWARD) then

! ===== 1D FFTs in X =====
#ifdef OVERWRITE
if (d2d_intranode) then
call c2c_1m_x(in%cvar2d,isign,ph)
else
call c2c_1m_x(in%cvar,isign,ph)
endif
#else
call wk1%init(is_cplx=.true., idir=1, decomp=ph)
if (d2d_intranode) then
wk1%cvar2d = in%cvar2d
call c2c_1m_x(wk1%cvar2d,isign,ph)
else
wk1%cvar = in%cvar
call c2c_1m_x(wk1%cvar,isign,ph)
endif
#endif

! ===== Swap X --> Y; 1D FFTs in Y =====

#ifdef OVERWRITE
call transpose_x_to_y(in,wk2_c2c)
#else
call transpose_x_to_y(wk1,wk2_c2c)
#endif
if (d2d_intranode) then
call c2c_1m_y(wk2_c2c%cvar2d,isign,ph)
else
call c2c_1m_y(wk2_c2c%cvar,isign,ph)
endif

! ===== Swap Y --> Z; 1D FFTs in Z =====
call transpose_y_to_z(wk2_c2c,out)
if (d2d_intranode) then
call c2c_1m_z(out%cvar2d,isign,ph)
else
call c2c_1m_z(out%cvar,isign,ph)
endif

else if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_BACKWARD &
.OR. & 
format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_FORWARD) then

! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
if (d2d_intranode) then
call c2c_1m_z(in%cvar2d,isign,ph)
else
call c2c_1m_z(in%cvar,isign,ph)
endif
#else
call wk1%init(is_cplx=.true., idir=3, decomp=ph)
if (d2d_intranode) then
wk1%cvar2d = in%cvar2d
call c2c_1m_z(wk1%cvar2d,isign,ph)
else
wk1%cvar = in%cvar
call c2c_1m_z(wk1%cvar,isign,ph)
endif
#endif

! ===== Swap Z --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
call transpose_z_to_y(in,wk2_c2c)
#else
call transpose_z_to_y(wk1,wk2_c2c)
#endif
if (d2d_intranode) then
call c2c_1m_y(wk2_c2c%cvar2d,isign,ph)
else
call c2c_1m_y(wk2_c2c%cvar,isign,ph)
endif

! ===== Swap Y --> X; 1D FFTs in X =====
call transpose_y_to_x(wk2_c2c,out)
if (d2d_intranode) then
call c2c_1m_x(out%cvar2d,isign,ph)
else
call c2c_1m_x(out%cvar,isign,ph)
endif

end if

#ifndef OVERWRITE
! Free memory
call wk1%fin
#endif

return
end subroutine fft_3d_c2c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D forward FFT - real to complex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_3d_rc2rc(in, out)

implicit none

type(decomp_data), intent(INOUT) :: in
type(decomp_data), intent(INOUT) :: out

if (in%is_cplx .and. (.not.out%is_cplx)) call fft_3d_c2r(in, out)

if ((.not.in%is_cplx) .and. out%is_cplx) call fft_3d_r2c(in, out)

end subroutine fft_3d_rc2rc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D forward FFT - real to complex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_3d_r2c(in_r, out_c)

implicit none

type(decomp_data), intent(INOUT) :: in_r
type(decomp_data), intent(INOUT) :: out_c

if (format==PHYSICAL_IN_X) then

! ===== 1D FFTs in X =====
if (d2d_intranode) then
call r2c_1m_x(in_r%var2d,wk13%cvar2d)
else
call r2c_1m_x(in_r%var,wk13%cvar)
endif

! ===== Swap X --> Y; 1D FFTs in Y =====
call transpose_x_to_y(wk13,wk2_r2c)
if (d2d_intranode) then
call c2c_1m_y(wk2_r2c%cvar2d,-1,sp)
else
call c2c_1m_y(wk2_r2c%cvar,-1,sp)
endif

! ===== Swap Y --> Z; 1D FFTs in Z =====
call transpose_y_to_z(wk2_r2c,out_c)
if (d2d_intranode) then
call c2c_1m_z(out_c%cvar2d,-1,sp)
else
call c2c_1m_z(out_c%cvar,-1,sp)
endif

else if (format==PHYSICAL_IN_Z) then

! ===== 1D FFTs in Z =====
if (d2d_intranode) then
call r2c_1m_z(in_r%var2d,wk13%cvar2d)
else
call r2c_1m_z(in_r%var,wk13%cvar)
endif

! ===== Swap Z --> Y; 1D FFTs in Y =====
call transpose_z_to_y(wk13,wk2_r2c)
if (d2d_intranode) then
call c2c_1m_y(wk2_r2c%cvar2d,-1,sp)
else
call c2c_1m_y(wk2_r2c%cvar,-1,sp)
endif

! ===== Swap Y --> X; 1D FFTs in X =====
call transpose_y_to_x(wk2_r2c,out_c)
if (d2d_intranode) then
call c2c_1m_x(out_c%cvar2d,-1,sp)
else
call c2c_1m_x(out_c%cvar,-1,sp)
endif

endif

end subroutine fft_3d_r2c


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D inverse FFT - complex to real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_3d_c2r(in_c, out_r)

implicit none

! Arguments
type(decomp_data), intent(INOUT) :: in_c
type(decomp_data), intent(INOUT) :: out_r

! local variables
#ifndef OVERWRITE
type(decomp_data) :: wk1
#endif

if (format==PHYSICAL_IN_X) then

! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
if (d2d_intranode) then
call c2c_1m_z(in_c%cvar2d,1,sp)
else
call c2c_1m_z(in_c%cvar,1,sp)
endif
#else
call wk1%init(is_cplx=.true., idir=1, decomp=sp)
if (d2d_intranode) then
wk1%cvar2d = in_c%cvar2d
call c2c_1m_z(wk1%cvar2d,1,sp)
else
wk1%cvar = in_c%cvar
call c2c_1m_z(wk1%cvar,1,sp)
endif
#endif

! ===== Swap Z --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
call transpose_z_to_y(in_c,wk2_r2c)
#else
call transpose_z_to_y(wk1,wk2_r2c)
#endif
if (d2d_intranode) then
call c2c_1m_y(wk2_r2c%cvar2d,1,sp)
else
call c2c_1m_y(wk2_r2c%cvar,1,sp)
endif

! ===== Swap Y --> X; 1D FFTs in X =====
call transpose_y_to_x(wk2_r2c,wk13)
if (d2d_intranode) then
call c2r_1m_x(wk13%cvar2d,out_r%var2d)
else
call c2r_1m_x(wk13%cvar,out_r%var)
endif

else if (format==PHYSICAL_IN_Z) then

! ===== 1D FFTs in X =====
#ifdef OVERWRITE
if (d2d_intranode) then
call c2c_1m_x(in_c%cvar2d,1,sp)
else
call c2c_1m_x(in_c%cvar,1,sp)
endif
#else
call wk1%init(is_cplx=.true., idir=1, decomp=sp)
if (d2d_intranode) then
wk1%cvar2d = in_c%cvar2d
call c2c_1m_x(wk1%cvar2d,1,sp)
else
wk1%cvar = in_c%cvar
call c2c_1m_x(wk1%cvar,1,sp)
endif
#endif

! ===== Swap X --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
call transpose_x_to_y(in_c,wk2_r2c)
#else
call transpose_x_to_y(wk1,wk2_r2c)
#endif
if (d2d_intranode) then
call c2c_1m_y(wk2_r2c%cvar2d,1,sp)
else
call c2c_1m_y(wk2_r2c%cvar,1,sp)
endif

! ===== Swap Y --> Z; 1D FFTs in Z =====
call transpose_y_to_z(wk2_r2c,wk13)
if (d2d_intranode) then
call c2r_1m_z(wk13%cvar2d,out_r%var2d)
else
call c2r_1m_z(wk13%cvar,out_r%var)
endif

endif

#ifndef OVERWRITE
! Free memory
call wk1%fin
#endif

return
end subroutine fft_3d_c2r
