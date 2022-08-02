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

! This is the 'generic' implementation of the FFT library

module decomp_2d_fft

  use decomp_2d  ! 2D decomposition module
  use glassman

  implicit none

  private        ! Make everything private unless declared public

  ! engine-specific global variables
  complex(mytype), allocatable, dimension(:) :: buf, scratch

  ! common code used for all engines, including global variables, 
  ! generic interface definitions and several subroutines
#include "fft_common.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time initialisations for the FFT engine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_fft_engine

    implicit none

    integer :: cbuf_size

    if (nrank==0 .and. nrank_loc<=0) then
       write(*,*) ' '
       write(*,*) '***** Using the generic FFT engine *****'
       write(*,*) ' '
    end if

    cbuf_size = max(ph%xsz(1), ph%ysz(1))
    cbuf_size = max(cbuf_size, ph%zsz(1))
    allocate(buf(cbuf_size))
    allocate(scratch(cbuf_size))

    return
  end subroutine init_fft_engine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time finalisations for the FFT engine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine finalize_fft_engine

    implicit none

    if (allocated(buf)) deallocate(buf)
    if (allocated(scratch)) deallocate(scratch)

    return
  end subroutine finalize_fft_engine


  ! Following routines calculate multiple one-dimensional FFTs to form 
  ! the basis of three-dimensional FFTs.

  ! c2c transform, multiple 1D FFTs
  subroutine c2c_1m_2d(inout, isign, n1, n2, idir)

    !$acc routine(spcfft) seq

    implicit none

    complex(mytype), dimension(:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign, n1, n2
    integer, intent(in), optional :: idir

    integer :: i,j

    ! This is needed before writing to the MPI3 shared memory
    if (present(idir).and.d2d_intranode) &
       call decomp_2d_win_transpose_start_writing(work1_c_win)

    !$acc parallel loop gang vector private(buf, scratch)
    do j=1,n2
       !
       ! Compute 1D FFT
       do i=1,n1
          buf(i) = inout(i,j)
       end do
       call spcfft(buf,n1,isign,scratch)
       do i=1,n1
          inout(i,j) = buf(i)
       end do
       !
       ! Write to the MPI3 shared memory
       ! Construction of the sending buffer in progress
       if (present(idir).and.d2d_intranode) then
          select case (idir)
          case (12)
             ! x => y, ph
             do i = 1, n1
                work1_c(ph%intramap_split(i,j,1)) = buf(i)
             enddo
          case (120)
             ! x => y, sp
             do i = 1, n1
                work1_c(sp%intramap_split(i,j,1)) = buf(i)
             enddo
          case (23)
             ! y => z, ph
             do i = 1, n1
                work1_c(ph%intramap_split(i,j,3)) = buf(i)
             enddo
          case (230)
             ! y => z, sp
             do i = 1, n1
                work1_c(sp%intramap_split(i,j,3)) = buf(i)
             enddo
          case (32)
             ! z => y, ph
             do i = 1, n1
                work1_c(ph%intramap_split(i,j,4)) = buf(i)
             enddo
          case(320)
             ! z => y, sp
             do i = 1, n1
                work1_c(sp%intramap_split(i,j,4)) = buf(i)
             enddo
          case (21)
             ! y => x, ph
             do i = 1, n1
                work1_c(ph%intramap_split(i,j,2)) = buf(i)
             enddo
          case (210)
             ! y => x, sp
             do i = 1, n1
                work1_c(sp%intramap_split(i,j,2)) = buf(i)
             enddo
          end select
       endif
    end do
    !$acc end parallel loop

  end subroutine c2c_1m_2d

  ! c2c transform, multiple 1D FFTs in x direction
  subroutine c2c_1m_x_2d(inout, isign, decomp, idir)

    implicit none

    complex(mytype), dimension(:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(in), optional :: idir

    call c2c_1m_2d(inout, isign, decomp%xsz_loc(1), decomp%xsz_loc(2)*decomp%xsz_loc(3), idir)

  end subroutine c2c_1m_x_2d

  ! c2c transform, multiple 1D FFTs in y direction
  subroutine c2c_1m_y_2d(inout, isign, decomp, idir)

    implicit none

    complex(mytype), dimension(:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(in), optional :: idir

    call c2c_1m_2d(inout, isign, decomp%ysz_loc(1), decomp%ysz_loc(2)*decomp%ysz_loc(3), idir)

  end subroutine c2c_1m_y_2d

  ! c2c transform, multiple 1D FFTs in z direction
  subroutine c2c_1m_z_2d(inout, isign, decomp, idir)

    implicit none

    complex(mytype), dimension(:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(in), optional :: idir

    call c2c_1m_2d(inout, isign, decomp%zsz_loc(1), decomp%zsz_loc(2)*decomp%zsz_loc(3), idir)

  end subroutine c2c_1m_z_2d

  ! r2c transform, multiple 1D FFTs in x direction
  subroutine r2c_1m_x_2d(input, output, idir)

    !$acc routine(spcfft) seq

    implicit none

    real(mytype), dimension(:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:), intent(OUT) :: output
    integer, intent(in), optional :: idir

    integer :: i,j, s1,s2, d1

    s1 = size(input,1)
    s2 = size(input,2)
    d1 = size(output,1)

    ! This is needed before writing to the MPI3 shared memory
    if (present(idir).and.d2d_intranode) &
       call decomp_2d_win_transpose_start_writing(work1_c_win)

    !$acc parallel loop gang vector private(buf, scratch)
    do j=1,s2
       ! Glassman's FFT is c2c only, 
       ! needing some pre- and post-processing for r2c
       ! pack real input in complex storage
       do i=1,s1
          buf(i) = cmplx(input(i,j),0._mytype, kind=mytype)
       end do
       call spcfft(buf,s1,-1,scratch)
       ! note d1 ~ s1/2+1
       ! simply drop the redundant part of the complex output
       do i=1,d1
          output(i,j) = buf(i)
       end do
       !
       ! Write to the MPI3 shared memory
       ! Construction of the sending buffer in progress
       if (present(idir).and.d2d_intranode) then
          if (idir==120) then
             ! x => y, sp
             do i = 1, d1
                work1_c(sp%intramap_split(i,j,1)) = buf(i)
             enddo
          else
             ! z => y, sp
             do i = 1, d1
                work1_c(sp%intramap_split(i,j,4)) = buf(i)
             enddo
          endif
       endif
    end do
    !$acc end parallel loop

  end subroutine r2c_1m_x_2d

  ! r2c transform, multiple 1D FFTs in z direction
  subroutine r2c_1m_z_2d(input, output, idir)

    implicit none

    real(mytype), dimension(:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:), intent(OUT) :: output
    integer, intent(in), optional :: idir

    call r2c_1m_x(input, output, idir)

  end subroutine r2c_1m_z_2d

  ! c2r transform, multiple 1D FFTs in x direction
  subroutine c2r_1m_x_2d(input, output)

    !$acc routine(spcfft) seq

    implicit none

    complex(mytype), dimension(:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:), intent(OUT) :: output

    integer :: i,j, d1,d2

    d1 = size(output,1)
    d2 = size(output,2)

    !$acc parallel loop gang vector collapse(2) private(buf, scratch)
    do j=1,d2
       ! Glassman's FFT is c2c only, 
       ! needing some pre- and post-processing for c2r
       do i=1,d1/2+1
          buf(i) = input(i,j)
       end do
       ! expanding to a full-size complex array
       ! For odd N, the storage is:
       !  1, 2, ...... N/2+1   integer division rounded down
       !     N, ...... N/2+2   => a(i) is conjugate of a(N+2-i)
       ! For even N, the storage is:
       !  1, 2, ...... N/2  , N/2+1
       !     N, ...... N/2+2  again a(i) conjugate of a(N+2-i)
       do i=d1/2+2,d1
          buf(i) =  conjg(buf(d1+2-i))
       end do
       call spcfft(buf,d1,1,scratch)
       do i=1,d1
          ! simply drop imaginary part
          output(i,j) = real(buf(i), kind=mytype)
       end do
    end do
    !$acc end parallel loop

  end subroutine c2r_1m_x_2d

  ! c2r transform, multiple 1D FFTs in z direction
  subroutine c2r_1m_z_2d(input, output)

    implicit none

    complex(mytype), dimension(:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:), intent(OUT) :: output

    call c2r_1m_x(input, output)

  end subroutine c2r_1m_z_2d


#include "fft_common_3d.f90"


end module decomp_2d_fft
