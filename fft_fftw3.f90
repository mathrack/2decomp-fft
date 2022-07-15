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

! This is the FFTW (version 3.x) implementation of the FFT library

module decomp_2d_fft

  use decomp_2d  ! 2D decomposition module

  implicit none

  include "fftw3.f"

  private        ! Make everything private unless declared public

  ! engine-specific global variables
  integer, save :: plan_type = FFTW_MEASURE

  ! FFTW plans
  ! j=1,2,3 corresponds to the 1D FFTs in X,Y,Z direction, respectively
  ! For c2c transforms: 
  !     use plan(-1,j) for  forward transform; 
  !     use plan( 1,j) for backward transform;
  ! For r2c/c2r transforms:
  !     use plan(0,j) for r2c transforms;
  !     use plan(2,j) for c2r transforms;
  integer*8, save :: plan(-1:2,3)

  ! common code used for all engines, including global variables, 
  ! generic interface definitions and several subroutines
#include "fft_common.f90"

  ! Return a FFTW3 plan for multiple 1D c2c FFTs
  subroutine c2c_1m_plan(plan1, isign, nn)

    implicit none

    ! Arguments
    integer*8, intent(OUT) :: plan1
    integer, intent(in) :: isign, nn(3)

    ! Local variables
    complex(mytype), allocatable, dimension(:,:) :: a2d

    allocate(a2d(nn(1),nn(2)*nn(3)))

#ifdef DOUBLE_PREC
    call dfftw_plan_many_dft(plan1, 1, nn(1), &
         nn(2)*nn(3), a2d, nn(1), 1, &
         nn(1), a2d, nn(1), 1, nn(1), &
         isign, plan_type)
#else
    call sfftw_plan_many_dft(plan1, 1, nn(1), &
         nn(2)*nn(3), a2d, nn(1), 1, &
         nn(1), a2d, nn(1), 1, nn(1), &
         isign, plan_type)
#endif

    deallocate(a2d)

  end subroutine c2c_1m_plan

  ! Return a FFTW3 plan for multiple 1D c2c FFTs in X direction
  subroutine c2c_1m_x_plan(plan1, decomp, isign)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: isign

    call c2c_1m_plan(plan1, isign, decomp%xsz_loc)

  end subroutine c2c_1m_x_plan

  ! Return a FFTW3 plan for multiple 1D c2c FFTs in Y direction
  subroutine c2c_1m_y_plan(plan1, decomp, isign)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: isign

    call c2c_1m_plan(plan1, isign, decomp%ysz_loc)

  end subroutine c2c_1m_y_plan

  ! Return a FFTW3 plan for multiple 1D c2c FFTs in Z direction
  subroutine c2c_1m_z_plan(plan1, decomp, isign)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: isign

    call c2c_1m_plan(plan1, isign, decomp%zsz_loc)

  end subroutine c2c_1m_z_plan

  ! Return a FFTW3 plan for multiple 1D r2c FFT
  subroutine r2c_1m_plan(plan1, nr, nc)

    implicit none

    ! Arguments
    integer*8, intent(OUT) :: plan1
    integer, intent(in) :: nr(3), nc(3)

    ! Local variables
    real(mytype), allocatable, dimension(:,:) :: b1
    complex(mytype), allocatable, dimension(:,:) :: b2

    allocate(b1(nr(1), nr(2) * nr(3))) 
    allocate(b2(nc(1), nc(2) * nc(3)))

#ifdef DOUBLE_PREC
    call dfftw_plan_many_dft_r2c(plan1, 1, nr(1), &
         nr(2)*nr(3), b1, nr(1), 1, &
         nr(1), b2, nc(1), 1, nc(1), &
         plan_type)
#else      
    call sfftw_plan_many_dft_r2c(plan1, 1, nr(1), &           
         nr(2)*nr(3), b1, nr(1), 1, &                                 
         nr(1), b2, nc(1), 1, nc(1), &                                 
         plan_type)
#endif

    deallocate(b1, b2)

  end subroutine r2c_1m_plan

  ! Return a FFTW3 plan for multiple 1D c2r FFT
  subroutine c2r_1m_plan(plan1, nc, nr)

    implicit none

    ! Arguments
    integer*8, intent(OUT) :: plan1
    integer, intent(in) :: nc(3), nr(3)

    ! Local variables
    complex(mytype), allocatable, dimension(:,:) :: b1
    real(mytype), allocatable, dimension(:,:) :: b2

    allocate(b1(nc(1), nc(2) * nc(3)))
    allocate(b2(nr(1), nr(2) * nr(3)))

#ifdef DOUBLE_PREC
    call dfftw_plan_many_dft_c2r(plan1, 1, nr(1), &
         nr(2)*nr(3), b1, nc(1), 1, &
         nc(1), b2, nr(1), 1, nr(1), &
         plan_type)
#else
    call sfftw_plan_many_dft_c2r(plan1, 1, nr(1), &
         nr(2)*nr(3), b1, nc(1), 1, &
         nc(1), b2, nr(1), 1, nr(1), &
         plan_type)
#endif

    deallocate(b1, b2)

  end subroutine c2r_1m_plan

  ! Return a FFTW3 plan for multiple 1D r2c FFTs in X direction
  subroutine r2c_1m_x_plan(plan1, decomp_ph, decomp_sp)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp

    call r2c_1m_plan(plan1, decomp_ph%xsz_loc, decomp_sp%xsz_loc)

  end subroutine r2c_1m_x_plan


  ! Return a FFTW3 plan for multiple 1D c2r FFTs in X direction
  subroutine c2r_1m_x_plan(plan1, decomp_sp, decomp_ph)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph

    call c2r_1m_plan(plan1, decomp_sp%xsz_loc, decomp_ph%xsz_loc)

  end subroutine c2r_1m_x_plan


  ! Return a FFTW3 plan for multiple 1D r2c FFTs in Z direction
  subroutine r2c_1m_z_plan(plan1, decomp_ph, decomp_sp)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp

    call r2c_1m_plan(plan1, decomp_ph%zsz_loc, decomp_sp%zsz_loc)

  end subroutine r2c_1m_z_plan


  ! Return a FFTW3 plan for multiple 1D c2r FFTs in Z direction
  subroutine c2r_1m_z_plan(plan1, decomp_sp, decomp_ph)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph

    call c2r_1m_plan(plan1, decomp_sp%zsz_loc, decomp_ph%zsz_loc)

  end subroutine c2r_1m_z_plan


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time initialisations for the FFT engine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_fft_engine

    implicit none

    if (nrank==0 .and. nrank_loc<=0) then
       write(*,*) ' '
       write(*,*) '***** Using the FFTW (version 3.x) engine *****'
       write(*,*) ' '
    end if

    if (format == PHYSICAL_IN_X) then

       ! For C2C transforms
       call c2c_1m_x_plan(plan(-1,1), ph, FFTW_FORWARD )
       call c2c_1m_y_plan(plan(-1,2), ph, FFTW_FORWARD )
       call c2c_1m_z_plan(plan(-1,3), ph, FFTW_FORWARD )
       call c2c_1m_z_plan(plan( 1,3), ph, FFTW_BACKWARD)
       call c2c_1m_y_plan(plan( 1,2), ph, FFTW_BACKWARD)
       call c2c_1m_x_plan(plan( 1,1), ph, FFTW_BACKWARD)

       ! For R2C/C2R tranforms
       call r2c_1m_x_plan(plan(0,1), ph, sp)
       call c2c_1m_y_plan(plan(0,2), sp, FFTW_FORWARD )
       call c2c_1m_z_plan(plan(0,3), sp, FFTW_FORWARD )
       call c2c_1m_z_plan(plan(2,3), sp, FFTW_BACKWARD)
       call c2c_1m_y_plan(plan(2,2), sp, FFTW_BACKWARD)
       call c2r_1m_x_plan(plan(2,1), sp, ph)

    else if (format == PHYSICAL_IN_Z) then

       ! For C2C transforms
       call c2c_1m_z_plan(plan(-1,3), ph, FFTW_FORWARD )
       call c2c_1m_y_plan(plan(-1,2), ph, FFTW_FORWARD ) 
       call c2c_1m_x_plan(plan(-1,1), ph, FFTW_FORWARD )
       call c2c_1m_x_plan(plan( 1,1), ph, FFTW_BACKWARD)
       call c2c_1m_y_plan(plan( 1,2), ph, FFTW_BACKWARD)
       call c2c_1m_z_plan(plan( 1,3), ph, FFTW_BACKWARD)

       ! For R2C/C2R tranforms
       call r2c_1m_z_plan(plan(0,3), ph, sp)
       call c2c_1m_y_plan(plan(0,2), sp, FFTW_FORWARD )
       call c2c_1m_x_plan(plan(0,1), sp, FFTW_FORWARD )
       call c2c_1m_x_plan(plan(2,1), sp, FFTW_BACKWARD)
       call c2c_1m_y_plan(plan(2,2), sp, FFTW_BACKWARD)
       call c2r_1m_z_plan(plan(2,3), sp, ph)

    end if

  end subroutine init_fft_engine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time finalisations for the FFT engine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine finalize_fft_engine

    implicit none

    integer :: i,j

    do j=1,3
       do i=-1,2
#ifdef DOUBLE_PREC
          call dfftw_destroy_plan(plan(i,j))
#else
          call sfftw_destroy_plan(plan(i,j))
#endif
       end do
    end do

  end subroutine finalize_fft_engine


  ! Following routines calculate multiple one-dimensional FFTs to form 
  ! the basis of three-dimensional FFTs.

  ! c2c transform, multiple 1D FFTs in x direction
  subroutine c2c_1m_x_2d(inout, isign, plan1)

    implicit none

    complex(mytype), dimension(:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    integer*8, intent(IN) :: plan1

#ifdef DOUBLE_PREC
    call dfftw_execute_dft(plan1, inout, inout)
#else
    call sfftw_execute_dft(plan1, inout, inout)
#endif

  end subroutine c2c_1m_x_2d


  ! c2c transform, multiple 1D FFTs in y direction
  subroutine c2c_1m_y_2d(inout, isign, plan1)

    implicit none

    complex(mytype), dimension(:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    integer*8, intent(IN) :: plan1

#ifdef DOUBLE_PREC
    call dfftw_execute_dft(plan1, inout, inout)
#else
    call sfftw_execute_dft(plan1, inout, inout)
#endif

  end subroutine c2c_1m_y_2d

  ! c2c transform, multiple 1D FFTs in z direction
  subroutine c2c_1m_z_2d(inout, isign, plan1)

    implicit none

    complex(mytype), dimension(:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    integer*8, intent(IN) :: plan1

#ifdef DOUBLE_PREC
    call dfftw_execute_dft(plan1, inout, inout)
#else
    call sfftw_execute_dft(plan1, inout, inout)
#endif

  end subroutine c2c_1m_z_2d

  ! r2c transform, multiple 1D FFTs in x direction
  subroutine r2c_1m_x_2d(input, output)

    implicit none

    real(mytype), dimension(:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:), intent(OUT) :: output

#ifdef DOUBLE_PREC
    call dfftw_execute_dft_r2c(plan(0,1), input, output)
#else
    call sfftw_execute_dft_r2c(plan(0,1), input, output)
#endif    

  end subroutine r2c_1m_x_2d

  ! r2c transform, multiple 1D FFTs in z direction
  subroutine r2c_1m_z_2d(input, output)

    implicit none

    real(mytype), dimension(:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:), intent(OUT) :: output

#ifdef DOUBLE_PREC
    call dfftw_execute_dft_r2c(plan(0,3), input, output)
#else
    call sfftw_execute_dft_r2c(plan(0,3), input, output)
#endif

  end subroutine r2c_1m_z_2d

  ! c2r transform, multiple 1D FFTs in x direction
  subroutine c2r_1m_x_2d(input, output)

    implicit none

    complex(mytype), dimension(:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:), intent(OUT) :: output

#ifdef DOUBLE_PREC
    call dfftw_execute_dft_c2r(plan(2,1), input, output)
#else
    call sfftw_execute_dft_c2r(plan(2,1), input, output)
#endif

  end subroutine c2r_1m_x_2d

  ! c2r transform, multiple 1D FFTs in z direction
  subroutine c2r_1m_z_2d(input, output)

    implicit none

    complex(mytype), dimension(:,:), intent(IN) :: input
    real(mytype), dimension(:,:), intent(OUT) :: output

#ifdef DOUBLE_PREC
    call dfftw_execute_dft_c2r(plan(2,3), input, output)
#else
    call sfftw_execute_dft_c2r(plan(2,3), input, output)
#endif    

  end subroutine c2r_1m_z_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D FFT - complex to complex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_c2c(in, out, isign)

    implicit none

    type(decomp_data), intent(inout) :: in
    type(decomp_data), intent(inout) :: out
    integer, intent(IN) :: isign

#ifndef OVERWRITE
    type(decomp_data) :: wk1
#endif

    if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_FORWARD .OR.  &
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_BACKWARD) then

       ! ===== 1D FFTs in X =====
#ifdef OVERWRITE
       call c2c_1m_x(in%cvar2d,isign,plan(isign,1))
#else
       call wk1%copy(in)
       wk1%cvar2d = in%cvar2d
       call c2c_1m_x(wk1%cvar2d,isign,plan(isign,1))
#endif

       ! ===== Swap X --> Y; 1D FFTs in Y =====

#ifdef OVERWRITE
       call transpose_x_to_y(in,wk2_c2c)
#else
       call transpose_x_to_y(wk1,wk2_c2c)
#endif
       call c2c_1m_y(wk2_c2c%cvar2d,isign,plan(isign,2))

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       call transpose_y_to_z(wk2_c2c,out)
       call c2c_1m_z(out%cvar2d,isign,plan(isign,3))

    else if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_BACKWARD &
         .OR. & 
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_FORWARD) then

       ! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
       call c2c_1m_z(in%cvar2d,isign,plan(isign,3))
#else
       call wk1%copy(in)
       wk1%cvar2d = in%cvar2d
       call c2c_1m_z(wk1%cvar2d,isign,plan(isign,3))
#endif

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
       call transpose_z_to_y(in,wk2_c2c)
#else
       call transpose_z_to_y(wk1,wk2_c2c)
#endif
       call c2c_1m_y(wk2_c2c%cvar2d,isign,plan(isign,2))

       ! ===== Swap Y --> X; 1D FFTs in X =====
       call transpose_y_to_x(wk2_c2c,out)
       call c2c_1m_x(out%cvar2d,isign,plan(isign,1))

    end if

#ifndef OVERWRITE
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

    type(decomp_data), intent(in) :: in_r
    type(decomp_data), intent(inout) :: out_c

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in X =====
       call r2c_1m_x(in_r%var2d,wk13%cvar2d)

       ! ===== Swap X --> Y; 1D FFTs in Y =====
       call transpose_x_to_y(wk13,wk2_r2c)
       call c2c_1m_y(wk2_r2c%cvar2d,-1,plan(0,2))

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       call transpose_y_to_z(wk2_r2c,out_c)
       call c2c_1m_z(out_c%cvar2d,-1,plan(0,3))

    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in Z =====
       call r2c_1m_z(in_r%var2d,wk13%cvar2d)

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
       call transpose_z_to_y(wk13,wk2_r2c)
       call c2c_1m_y(wk2_r2c%cvar2d,-1,plan(0,2))

       ! ===== Swap Y --> X; 1D FFTs in X =====
       call transpose_y_to_x(wk2_r2c,out_c)
       call c2c_1m_x(out_c%cvar2d,-1,plan(0,1))

    end if

    return
  end subroutine fft_3d_r2c


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D inverse FFT - complex to real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_c2r(in_c, out_r)

    implicit none

    type(decomp_data), intent(inout) :: in_c
    type(decomp_data), intent(inout) :: out_r

#ifndef OVERWRITE
    type(decomp_data) :: wk1
#endif

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
       call c2c_1m_z(in_c%cvar2d,1,plan(2,3))
#else
       call wk1%copy(in_c)
       wk1%cvar2d = in_c%cvar2d
       call c2c_1m_z(wk1%cvar2d,1,plan(2,3))
#endif

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
       call transpose_z_to_y(in_c,wk2_r2c)
#else
       call transpose_z_to_y(wk1,wk2_r2c)
#endif
       call c2c_1m_y(wk2_r2c%cvar2d,1,plan(2,2))

       ! ===== Swap Y --> X; 1D FFTs in X =====
       call transpose_y_to_x(wk2_r2c,wk13)
       call c2r_1m_x(wk13%cvar2d,out_r%var2d)

    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in X =====
#ifdef OVERWRITE
       call c2c_1m_x(in_c%cvar2d,1,plan(2,1))
#else
       call wk1%copy(in_c)
       wk1%cvar2d = in_c%cvar2d
       call c2c_1m_x(wk1%cvar2d,1,plan(2,1))
#endif

       ! ===== Swap X --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
       call transpose_x_to_y(in_c,wk2_r2c)
#else
       call transpose_x_to_y(wk1,wk2_r2c)
#endif
       call c2c_1m_y(wk2_r2c%cvar2d,1,plan(2,2))

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       call transpose_y_to_z(wk2_r2c,wk13)
       call c2r_1m_z(wk13%cvar2d,out_r%var2d)

    end if

#ifndef OVERWRITE
    call wk1%fin
#endif

    return
  end subroutine fft_3d_c2r


end module decomp_2d_fft
