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

! This is the submodule dedicated to decomp_data objects

submodule(decomp_2d) smod_data

   implicit none

contains

   !
   ! Initialize the object and allocate memory
   !
   module subroutine decomp_data_init(self, is_cplx, idir, decomp, rwk, cwk, contig)

      implicit none

      ! Arguments
      class(decomp_data), intent(out) :: self
      logical, intent(in) :: is_cplx
      integer, intent(in) :: idir
      type(decomp_info), intent(in), target, optional :: decomp
      real(mytype), dimension(:, :, :), target, optional :: rwk
      complex(mytype), dimension(:, :, :), target, optional :: cwk
      logical, optional :: contig

      ! Local variables
      logical, parameter :: glob = .false.

      ! Safety check
#ifdef DEBUG
      if (idir <= 0 .or. idir > 3 .or. &
          associated(self%var) .or. associated(self%cvar) .or. &
          (is_cplx .and. present(rwk)) .or. &
          ((.not. is_cplx) .and. present(cwk))) &
         call decomp_2d_abort(__FILE__, &
                              __LINE__, &
                              -1, &
                              "Incorrect arguments for decomp_data object init")
#endif

      ! Store key parameters
      self%is_cplx = is_cplx
      if (present(decomp)) then
         self%decomp => decomp
      else
         self%decomp => decomp_main
      end if
      self%idir = idir
      self%shm = d2d_intranode
      if (present(contig)) then
         self%contig = contig
      else
         self%contig = .false.
      endif

      ! Specific case : a 3D array was provided, no alloc needed
      if (present(rwk) .or. present(cwk)) then
         if (present(rwk)) self%var => rwk
         if (present(cwk)) self%cvar => cwk
         self%shm = .false.
         return
      end if

      ! x-y-z pencil decomp_info data
      if (self%idir == 1) then
         ! x-pencil
         if (self%shm) then
            ! MPI3 unified memory
            if (self%is_cplx) then
               call alloc_x(self%cvar, self%cvar2d, self%decomp, glob, self%win, self%contig)
            else
               call alloc_x(self%var, self%var2d, self%decomp, glob, self%win, self%contig)
            end if
         else
            ! Regular local object
            if (self%is_cplx) then
               call alloc_x(self%cvar, self%cvar2d, self%decomp, glob)
            else
               call alloc_x(self%var, self%var2d, self%decomp, glob)
            end if
         end if
      elseif (self%idir == 2) then
         ! y-pencil
         if (self%shm) then
            ! MPI3 unified memory
            if (self%is_cplx) then
               call alloc_y(self%cvar, self%cvar2d, self%decomp, glob, self%win, self%contig)
            else
               call alloc_y(self%var, self%var2d, self%decomp, glob, self%win, self%contig)
            end if
         else
            ! Regular local object
            if (self%is_cplx) then
               call alloc_y(self%cvar, self%cvar2d, self%decomp, glob)
            else
               call alloc_y(self%var, self%var2d, self%decomp, glob)
            end if
         end if
      else
         ! z-pencil
         if (self%shm) then
            ! MPI3 unified memory
            if (self%is_cplx) then
               call alloc_z(self%cvar, self%cvar2d, self%decomp, glob, self%win, self%contig)
            else
               call alloc_z(self%var, self%var2d, self%decomp, glob, self%win, self%contig)
            end if
         else
            ! Regular local object
            if (self%is_cplx) then
               call alloc_z(self%cvar, self%cvar2d, self%decomp, glob)
            else
               call alloc_z(self%var, self%var2d, self%decomp, glob)
            end if
         end if
      end if

   end subroutine decomp_data_init

   !
   ! Initialize the object using an existing one
   !
   ! THIS DOES NOT COPY THE CONTENT OF THE EXISTING OBJECT
   !
   module subroutine decomp_data_init_copy(self, dat)

      implicit none

      ! Arguments
      class(decomp_data), intent(out) :: self
      type(decomp_data), intent(in) :: dat

      call self%init(is_cplx = dat%is_cplx, &
                     idir = dat%idir, &
                     decomp = dat%decomp, &
                     contig = dat%contig)

   end subroutine decomp_data_init_copy

   !
   ! Free memory
   !
   module subroutine decomp_data_fin(self)

      implicit none

      ! Argument
      class(decomp_data), intent(inout) :: self

      self%is_cplx = .true.

      if (associated(self%decomp)) nullify(self%decomp)
      self%idir = 0

      if (self%shm) then
         call decomp_2d_win_free(self%win)
         if (associated(self%cvar)) nullify (self%cvar)
         if (associated(self%var)) nullify (self%var)
      else
         if (associated(self%cvar)) deallocate (self%cvar)
         if (associated(self%var)) deallocate (self%var)
      end if

      if (associated(self%cvar2d)) nullify (self%cvar2d)
      if (associated(self%var2d)) nullify (self%var2d)

      self%shm = .false.
      self%contig = .false.

   end subroutine decomp_data_fin

end submodule smod_data
