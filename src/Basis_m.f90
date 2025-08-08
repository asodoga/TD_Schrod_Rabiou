MODULE Basis_m
   USE QDUtil_m
   USE NDindex_m
   USE polyortho_m

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: Basis_t
   PUBLIC :: Read_Basis, Write_Basis
   PUBLIC :: Basis_IS_allocated, Basis_IS_allocatedtot
   PUBLIC :: Deallocate_Basis
   
   PUBLIC :: Calc_dngg_grid, Calc_dngg_grid_0
   PUBLIC :: Calc_tranpose_d0gb
   PUBLIC :: test_basitogridgridtobasis
   
   PUBLIC :: GridTOBasis_nD_cplx, BasisTOGrid_nD_cplx
   
   PUBLIC :: Calc_Q_grid, Calc_index
   PUBLIC :: Rdensity_alloc, Rdensity_Writing
   
   PUBLIC :: Scale_Basis, init_Basis1_TO_Basis2
   PUBLIC :: construct_primitive_basis, construct_primitive_basis_temp
   PUBLIC :: Construct_Basis_el,Construct_Basis_Fourier,Construct_Basis_Sin
   PUBLIC :: Construct_Basis_Ho,Construct_Basis_Hag
   PUBLIC :: CheckOrtho_Basis
   
   PUBLIC :: pdv2psi_nD
   
   PUBLIC :: Get_Basis_Parameters, Change_Basis_Parameters
   PUBLIC :: Calc_Basis_dPtQtAt
   
   PUBLIC :: Calc_d0d1d2W
   
   PUBLIC :: Calc_reduced_Density_surf, Calc_reduced_density
   PUBLIC :: REDUCED_DENSIRY_t
   PUBLIC :: calc_tab_Iq0 ,Complete_construct_Basis

   TYPE :: Basis_t
      integer                             :: nb_basis = ZERO
      integer                             :: nb ,nq
      real(kind=Rkind)                    :: A , B ,imp_k , scaleQ,Q0
      complex(kind=Rkind)                 :: Alpha = CZERO
      character(len=:),    allocatable    :: Basis_name
      real(kind=Rkind),    allocatable    :: x(:)
      real(kind=Rkind),    allocatable    :: w(:)
      complex(kind=Rkind), allocatable    :: d0gb(:, :) 
      complex(kind=Rkind), allocatable    :: d0bgw(:, :) 
      complex(kind=Rkind), allocatable    :: d1gb(:, :, :)
      complex(kind=Rkind), allocatable    :: d1gg(:, :, :)  
      complex(kind=Rkind), allocatable    :: d2gb(:, :, :, :)
      complex(kind=Rkind), allocatable    :: d2gg(:, :, :, :)
      complex(kind=Rkind), allocatable    :: dagb(:, :)
      complex(kind=Rkind), allocatable    :: dagg(:, :)
      complex(kind=Rkind), allocatable    :: dq0gb(:, :)
      complex(kind=Rkind), allocatable    :: dq0gg(:, :)
      complex(kind=Rkind), allocatable    :: dp0gb(:, :)
      complex(kind=Rkind), allocatable    :: dp0gg(:, :)
      complex(kind=Rkind), allocatable    :: S(:, :)
      TYPE(NDindex_t)                     :: NDindexq
      TYPE(NDindex_t)                     :: NDindexb
      TYPE(Basis_t), allocatable          :: tab_basis(:)     

   END TYPE Basis_t


   TYPE REDUCED_DENSIRY_t
      
      real(kind=Rkind),    allocatable    :: Norm(:)
      real(kind=Rkind),    allocatable    :: prob(:)
      TYPE(REDUCED_DENSIRY_t),allocatable :: tab_prob(:)
   END TYPE REDUCED_DENSIRY_t

CONTAINS

   RECURSIVE FUNCTION Basis_IS_allocated(Basis) RESULT(alloc)

      TYPE(Basis_t), intent(in)  :: Basis
      logical                      :: alloc
      integer                      :: i

      IF (Basis%Basis_name == 'el' .AND. Basis%nb > 0) Then
         alloc = .TRUE.
         RETURN
      END IF

      alloc = allocated(Basis%tab_basis)
      IF (allocated(Basis%tab_basis)) THEN
         Do i = 1, size(Basis%tab_basis)
            alloc = alloc .and. Basis_IS_allocated(Basis%tab_basis(i))
         END DO
      ELSE
         alloc = allocated(Basis%x)
         alloc = alloc .AND. allocated(Basis%w)
         alloc = alloc .AND. allocated(Basis%d0gb)
         alloc = alloc .AND. allocated(Basis%d1gb)
         alloc = alloc .AND. allocated(Basis%d2gb)
         !alloc = alloc .AND. allocated(Basis%d0bgw)
      END IF
   END FUNCTION Basis_IS_allocated

   RECURSIVE FUNCTION Basis_IS_allocatedtot(Basis) RESULT(alloc)

      TYPE(Basis_t), intent(in)    :: Basis
      logical                      :: alloc
      integer                      :: i

      alloc = allocated(Basis%tab_basis)
      IF (allocated(Basis%tab_basis)) THEN
         Do i = 1, size(Basis%tab_basis)
            alloc = alloc .and. Basis_IS_allocated(Basis%tab_basis(i))
         END DO
      ELSE
         alloc = allocated(Basis%x)
         alloc = alloc .AND. allocated(Basis%w)
         alloc = alloc .AND. allocated(Basis%d0gb)
         alloc = alloc .AND. allocated(Basis%d1gb)
         alloc = alloc .AND. allocated(Basis%d2gb)
         alloc = alloc .AND. allocated(Basis%d1gg)
         alloc = alloc .AND. allocated(Basis%d2gg)

      END IF

   END FUNCTION Basis_IS_allocatedtot

   RECURSIVE SUBROUTINE Write_Basis(Basis)
      USE QDUtil_m

      TYPE(Basis_t), intent(in)  :: Basis
      integer                          :: i

      !write(out_unit,*) '---------------------------------------------------------------------'
      !write(out_unit,*) 'Write_Basis'
      write(out_unit,*) "Basis_name=", Basis%Basis_name
      !write(out_unit,*) "n_basis=", Basis%nb_basis
      !write(out_unit,*) 'nb,nq',Basis%nb,Basis%nq
      !write(out_unit,*) "Q0,scaleQ ", Basis%Q0,Basis%scaleQ
      !write(out_unit,*)  'A,B',Basis%A,Basis%B
      !write(out_unit,*) 'nb_basis',Basis%nb_basis

      IF (.NOT. allocated(Basis%x)) THEN
         write(out_unit,*)' Basis table x is not allocated.'
      ELSE
         call  Write_VecMat(Basis%x, out_unit, 5,  info='x')
      END IF
      write (out_unit, *)
      IF (.NOT. allocated(Basis%W)) THEN
           write(out_unit,*)' Basis table w is not allocated.'
      ELSE
         call  Write_VecMat(Basis%w, out_unit, 5,  info='w')
      END IF
      write (out_unit, *)
       IF (.NOT.allocated(Basis%d0gb)) THEN
        write(out_unit,*)' Basis table d0gb is not allocated.'
      ELSE
      call   Write_VecMat(Basis%d0gb, out_unit, 5,  info='d0gb')
       END IF
          

       write (out_unit, *)
       IF (.NOT.allocated(Basis%dp0gb)) THEN
        write(out_unit,*)' Basis table dp0gb is not allocated.'
      ELSE
      call   Write_VecMat(Basis%dp0gb, out_unit, 5,  info='dp0gb')
       END IF
      
       
       write (out_unit, *)
       IF (.NOT.allocated(Basis%dq0gb)) THEN
        write(out_unit,*)' Basis table dq0gb is not allocated.'
      ELSE
      call   Write_VecMat(Basis%dq0gb, out_unit, 5,  info='dq0gb')
       END IF
       
       write (out_unit, *)
       IF (.NOT.allocated(Basis%dagb)) THEN
        write(out_unit,*)' Basis table dagb is not allocated.'
      ELSE
      call   Write_VecMat(Basis%dagb, out_unit, 5,  info='dagb')
       END IF

       write (out_unit, *)
       IF (.NOT.allocated(Basis%dagg)) THEN
        write(out_unit,*)' Basis table dagg is not allocated.'
      ELSE
      call   Write_VecMat(Basis%dagg, out_unit, 5,  info='dagg')
       END IF

       write (out_unit, *)
       IF (.NOT.allocated(Basis%dp0gg)) THEN
        write(out_unit,*)' Basis table dp0gg is not allocated.'
      ELSE
      call   Write_VecMat(Basis%dp0gg, out_unit, 5,  info='dp0gg')
       END IF
      
       write (out_unit, *)
       IF (.NOT.allocated(Basis%dq0gg)) THEN
        write(out_unit,*)' Basis table dq0gg is not allocated.'
      ELSE
      call   Write_VecMat(Basis%dq0gg, out_unit, 5,  info='dq0gg')
       END IF

      write(out_unit,*)
        IF (.NOT.allocated(Basis%d0bgw)) THEN
           write(out_unit,*)' Basis table d0bgw is not allocated.'
         ELSE
      call   Write_VecMat(Basis%d0bgw,out_unit,5, info='d0gbw')
        END IF

         write(out_unit,*)
         IF (.NOT.allocated(Basis%d1gb)) THEN
           write(out_unit,*)' Basis table d1gb is not allocated.'
         ELSE
      call   Write_VecMat(Basis%d1gb(:,:,1),out_unit,5, info='d1gb')
         END IF
         write(out_unit,*)
         IF (.NOT.allocated(Basis%d1gg)) THEN
           write(out_unit,*)' Basis table d1gb is not allocated.'
         ELSE
      call   Write_VecMat(Basis%d1gg(:,:,1),out_unit,5, info='d1gg')
         END IF

         write(out_unit,*)
         IF (.NOT.allocated(Basis%d2gb)) THEN
           write(out_unit,*)' Basis table d1gb is not allocated.'
         ELSE
      call   Write_VecMat(Basis%d2gb(:,:,1,1),out_unit,5, info='d2gb')
         END IF

         write(out_unit,*)
         IF (.NOT.allocated(Basis%d2gg)) THEN
           write(out_unit,*)' Basis table d2gg is not allocated.'
         ELSE
      call   Write_VecMat(Basis%d2gg(:,:,1,1),out_unit,5, info='d2gg')
         END IF
          

         write(out_unit,*)
         IF (.NOT.allocated(Basis%S)) THEN
            write(out_unit,*)' Basis table S is not allocated.'
         ELSE
      !call  Write_VecMat(Basis%S,out_unit,5, info='S')
         END IF  
         IF (allocated(Basis%tab_basis)) THEN
            DO i = 1, size(Basis%tab_basis)
              CALL Write_Basis(Basis%tab_basis(i))
            END DO
         END IF
      write (out_unit, *) '--------------------------------------------------------------------------'

   END SUBROUTINE Write_Basis


   !> Recursive subroutine to initialize Basis2 from Basis1
   !! It supports both composite (multi-dimensional) and 1D primitive basis structures.
   RECURSIVE SUBROUTINE init_Basis1_TO_Basis2(Basis2, Basis1)
   USE QDUtil_m
   IMPLICIT NONE
     
     TYPE(Basis_t), intent(in)           :: Basis1   ! Input basis structure
     TYPE(Basis_t), intent(inout)        :: Basis2   ! Output (initialized) basis structure
     integer                             :: ib        ! Loop index
     
     IF (allocated(Basis1%tab_basis)) THEN
        ! Composite basis: recursively initialize each sub-basis
        call Deallocate_Basis(Basis2)
        Basis2%Basis_name  = Basis1%Basis_name
        Basis2%nb_basis    = Basis1%nb_basis
        allocate(Basis2%tab_basis(Basis2%nb_basis))
     
        DO ib = 1, Basis1%nb_basis
           call init_Basis1_TO_Basis2(Basis2%tab_basis(ib), Basis1%tab_basis(ib))
        END DO
     
        ! Recalculate total number of basis functions and quadrature points
        Basis2%nb = 1
        Basis2%nq = 1
        DO ib = 1, Basis1%nb_basis
           IF (Basis2%tab_basis(ib)%Basis_name == 'el') CYCLE
           Basis2%nb = Basis2%nb * Basis2%tab_basis(ib)%nb
           Basis2%nq = Basis2%nq * Basis2%tab_basis(ib)%nq
        END DO
     
     ELSE
        ! Elementary basis: directly copy parameters
        Basis2%Basis_name  = Basis1%Basis_name
        Basis2%nb_basis    = Basis1%nb_basis
        Basis2%nb          = Basis1%nb
        Basis2%nq          = Basis1%nq
        Basis2%Q0          = Basis1%Q0
        Basis2%SCALEQ      = Basis1%SCALEQ
        Basis2%Imp_k       = Basis1%Imp_k
        Basis2%Alpha       = Basis1%Alpha
        Basis2%A           = Basis1%A
        Basis2%B           = Basis1%B
        Basis2%S           = Basis1%S
     END IF
     
     END SUBROUTINE init_Basis1_TO_Basis2

   !> Recursive subroutine to safely deallocate all allocated components
   !! of a Basis_t object, including nested basis structures.
   RECURSIVE SUBROUTINE Deallocate_Basis(Basis)
   USE QDUtil_m
   IMPLICIT NONE
   
   TYPE(Basis_t), intent(inout) :: Basis
   integer                      :: i
   
   ! Deallocate ND indexing arrays if present
   IF (allocated(Basis%tab_basis)) THEN
      IF (allocated(Basis%NDindexq%Tab0)) deallocate(Basis%NDindexq%Tab0)
      IF (allocated(Basis%NDindexb%Tab0)) deallocate(Basis%NDindexb%Tab0)
      deallocate(Basis%tab_basis)
   END IF
   
   ! Deallocate 1D basis arrays if present
   IF (allocated(Basis%x)) THEN
      deallocate(Basis%x)
      write(out_unit,*) 'Basis table x is now deallocated.'
   END IF
   
   IF (allocated(Basis%d0gb)) THEN
      deallocate(Basis%d0gb)
      write(out_unit,*) 'Basis table d0gb is now deallocated.'
   END IF
   
   IF (allocated(Basis%d0bgw)) THEN
      deallocate(Basis%d0bgw)
      write(out_unit,*) 'Basis table d0bgw is now deallocated.'
   END IF
   
   IF (allocated(Basis%dagb)) THEN
      deallocate(Basis%dagb)
      write(out_unit,*) 'Basis table dagb is now deallocated.'
   END IF
   
   IF (allocated(Basis%dagg)) THEN
      deallocate(Basis%dagg)
      write(out_unit,*) 'Basis table dagg is now deallocated.'
   END IF
   
   IF (allocated(Basis%dq0gb)) THEN
      deallocate(Basis%dq0gb)
      write(out_unit,*) 'Basis table dq0gb is now deallocated.'
   END IF
   
   IF (allocated(Basis%dq0gg)) THEN
      deallocate(Basis%dq0gg)
      write(out_unit,*) 'Basis table dq0gg is now deallocated.'
   END IF
   
   IF (allocated(Basis%dp0gg)) THEN
      deallocate(Basis%dp0gg)
      write(out_unit,*) 'Basis table dp0gg is now deallocated.'
   END IF
   
   IF (allocated(Basis%d1gb)) THEN
      deallocate(Basis%d1gb)
      write(out_unit,*) 'Basis table d1gb is now deallocated.'
   END IF
   
   IF (allocated(Basis%d1gg)) THEN
      deallocate(Basis%d1gg)
      write(out_unit,*) 'Basis table d1gg is now deallocated.'
   END IF
   
   IF (allocated(Basis%d2gb)) THEN
      deallocate(Basis%d2gb)
      write(out_unit,*) 'Basis table d2gb is now deallocated.'
   END IF
   
   IF (allocated(Basis%d2gg)) THEN
      deallocate(Basis%d2gg)
      write(out_unit,*) 'Basis table d2gg is now deallocated.'
   END IF
   
   ! Recursively deallocate sub-basis if composite
   IF (allocated(Basis%tab_basis)) THEN
      DO i = 1, size(Basis%tab_basis)
         CALL Deallocate_Basis(Basis%tab_basis(i))
      END DO
   END IF
   
   END SUBROUTINE Deallocate_Basis

   RECURSIVE SUBROUTINE Read_Basis(Basis, nio)
      USE QDUtil_m
      logical, parameter                       :: debug = .true.
      !logical,             parameter          ::debug = .false.
      TYPE(Basis_t), intent(inout)             :: Basis
      integer, intent(in)                      :: nio
      integer                                  :: err_io, nb, nq, i, j, nb_basis, ib
      character(len=Name_len)                  :: name
      real(kind=Rkind)                         :: A, B, scaleQ, Q0, d0, d2, X1, W1,Imp_k
      complex(kind=Rkind)                      :: Alpha

      NAMELIST /basis_nD/ name, nb_basis, nb, nq, A, B, scaleQ, Q0,Imp_k,Alpha
      nb_basis = 0
      nb = 0
      nq = 0
      A = ZERO
      B = ZERO
      Q0 = ZERO
      scaleQ = ONE
      name = '0'
      Imp_k= ZERO
      Alpha = CONE

      read (nio, nml=basis_nD, IOSTAT=err_io)
      write (out_unit, nml=basis_nD)
      IF (err_io < 0) THEN
         write (out_unit, basis_nD)
         write (out_unit, *) ' ERROR in Read_Basis'
         write (out_unit, *) ' while reading the namelist "basis_nD"'
         write (out_unit, *) ' end of file or end of record'
         write (out_unit, *) ' Probably, you forget a basis set ...'
         write (out_unit, *) ' Check your data !!'
         STOP ' ERROR in Read_Basis: problems with the namelist.'
      END IF
      IF (err_io > 0) THEN
         write (out_unit, basis_nD)
         write (out_unit, *) ' ERROR in Read_Basis'
         write (out_unit, *) ' while reading the namelist "basis_nD"'
         write (out_unit, *) ' Probably, some arguments of namelist are wrong.'
         write (out_unit, *) ' Check your data !!'
         STOP ' ERROR in Read_Basis: problems with the namelist.'
      END IF

      IF (nb_basis > 1) THEN
         Basis%Basis_name = 'Dp'
         Basis%nb_basis = nb_basis
         call string_uppercase_TO_lowercase(Basis%Basis_name)
         allocate (Basis%tab_basis(nb_basis))
         DO i = 1, nb_basis
            CALL Read_Basis(Basis%tab_basis(i), nio)
         END DO
         Basis%nb = 1
         Basis%nq = 1
         DO i = 1, nb_basis
            if (Basis%tab_basis(i)%Basis_name == 'el') cycle
            Basis%nb = Basis%nb*Basis%tab_basis(i)%nb
            Basis%nq = Basis%nq*Basis%tab_basis(i)%nq
         END DO
      ELSE
         Basis%nb_basis = nb_basis
         Basis%nb = nb
         Basis%nq = nq
         Basis%Q0 = Q0
         Basis%SCALEQ = SCALEQ
         Basis%A = A
         Basis%B = B
         Basis%Imp_k = Imp_k
         Basis%alpha = alpha
         Basis%Basis_name = trim(adjustl(name))
         allocate (Basis%S(nb, nb))
         Basis%S(:, :) = ZERO
         do ib = 1, Basis%nb
            Basis%S(ib, ib) =CONE
         end do

         CALL string_uppercase_TO_lowercase(Basis%Basis_name)
      END IF
   END SUBROUTINE Read_Basis

   
   !> Calculates the index table Tab_iq for multi-dimensional quadrature indexing
   !! Tab_iq(:, iq) contains the index vector for the iq-th point in the total basis grid.
   SUBROUTINE Calc_tab_Iq0(Tab_Iq, Basis)
      USE QDUtil_m
      ! Purpose: Calculate the table of quantum indices (Tab_Iq) for a given basis
      
      ! Input/Output variables:
      TYPE(Basis_t), intent(in), target :: Basis       ! Basis set information
      integer, allocatable, intent(inout) :: Tab_iq(:,:)  ! Output table of quantum indices
      
      ! Local variables:
      integer, allocatable :: Tab_iq0(:)    ! Temporary array for current quantum indices
      integer :: ndim                      ! Number of dimensions
      integer :: iq                        ! Current quantum state index
      integer :: nq                        ! Total number of quantum states
      logical :: Endloop                   ! Flag to indicate end of loop
      
      ! Get dimensions from the basis set:
      ndim = size(Basis%tab_basis) - 1     ! Number of dimensions (subtracting 1 for 0-based or other reason)
      nq = Basis%nq                        ! Total number of quantum states
      
      ! Allocate memory for the output table and temporary array:
      allocate(Tab_iq(ndim, nq), Tab_iq0(ndim))
      
      ! Initialize the temporary index array:
      Call Init_tab_ind(Tab_iq0, Basis%NDindexq)
      
      ! Loop through all quantum states:
      Iq = 0
      DO
          Iq = Iq + 1  ! Increment quantum state counter
          
          ! Get the next set of quantum numbers:
          CALL increase_NDindex(Tab_Iq0, Basis%NDindexq, Endloop)
          
          ! Exit if we've processed all quantum states:
          IF (Endloop) exit
          
          ! Store the current quantum numbers in the output table:
          Tab_iq(:, Iq) = Tab_Iq0
      END DO
      
      ! Clean up temporary array:
      deallocate(Tab_Iq0)
  END SUBROUTINE Calc_tab_Iq0

   RECURSIVE SUBROUTINE construct_primitive_basis_temp(Basis)
   USE QDUtil_m
   IMPLICIT NONE

   !==================== Parameters =====================
   LOGICAL, PARAMETER :: debug = .true.
   !LOGICAL, PARAMETER :: debug = .false.  ! Uncomment to disable debug mode

   !==================== Interface ======================
   TYPE(Basis_t), INTENT(INOUT) :: Basis   ! Basis structure to construct

   !==================== Local variables ================
   INTEGER :: ndim     ! Number of sub-bases
   INTEGER :: ib       ! Index for loop over sub-bases

   !==================== Logic ==========================
   ndim = Basis%nb_basis

   ! If the basis is composite (i.e., tab_basis is allocated),
   ! recursively construct the primitive basis of each subcomponent
   IF (ALLOCATED(Basis%tab_basis)) THEN
      DO ib = 1, ndim
         CALL construct_primitive_basis_temp(Basis%tab_basis(ib))
      END DO

   ELSE
      ! Otherwise, construct the basis based on its name/type
      SELECT CASE (Basis%Basis_name)

      CASE ('el')
         ! Electronic basis: no spatial degree of freedom
         WRITE(6, *) 'Electronic basis. Electronic state number:', Basis%nb
         Basis%nq = 0

      CASE ('boxab')
         ! Particle-in-a-box basis between [A, B]
         CALL Construct_Basis_Sin(Basis)
         Basis%Q0     = Basis%A
         Basis%scaleQ = pi / (Basis%B - Basis%A)
         CALL Complete_construct_Basis(Basis)

      CASE ('fourier')
         ! Plane wave / Fourier basis (unbounded domain)
         Basis%Q0     = ZERO
         Basis%scaleQ = ONE
         CALL Construct_Basis_Fourier(Basis)
         CALL Complete_construct_Basis(Basis)

      CASE ('hag', 'ho')
         ! Harmonic or Hagedorn basis (semiclassical Gaussian functions)
         !CALL Construct_Basis_Ho(Basis)  ! Uncomment if needed
         CALL Construct_Basis_Hag(Basis)
         CALL Complete_construct_Basis(Basis)

      CASE DEFAULT
         ! Unknown or unsupported basis type
         STOP 'ERROR  Nothing to construct'

      END SELECT
   END IF

   ! Optional debug output
   !WRITE(out_unit,*) 'End construct primitive Basis'

END SUBROUTINE construct_primitive_basis_temp

   RECURSIVE SUBROUTINE construct_primitive_basis(Basis)
      USE QDUtil_m
      logical, parameter                     :: debug = .true.
      !logical,             parameter        ::debug = .false.
      TYPE(Basis_t), intent(inout)           :: Basis
      integer, allocatable                   :: NDend_q(:)
      integer, allocatable                   :: NDend_b(:)
      integer                                :: ndim, ib

       ndim = Basis%nb_basis - 1
       allocate (NDend_q(ndim))
       allocate (NDend_b(ndim)) 
       DO ib = 1, ndim
          NDend_q(ib) = Basis%tab_basis(ib)%nq
          NDend_b(ib) = Basis%tab_basis(ib)%nb
       END DO 
       call Init_NDindex(Basis%NDindexq, NDend_q, ndim)
       call Init_NDindex(Basis%NDindexb, NDend_b, ndim)

       call construct_primitive_basis_temp(Basis)

        
   END SUBROUTINE 

    SUBROUTINE Construct_Basis_Sin(Basis) ! sin : boxAB with A=0 and B=pi
      USE QDUtil_m

      TYPE(Basis_t), intent(inout)        :: Basis
      real(kind=Rkind)                    :: dx
      integer                             :: ib, iq, nb, nq


      if (allocated(Basis%x)) deallocate (Basis%x)
      if (allocated(Basis%w)) deallocate (Basis%w)
      if (allocated(Basis%d0gb)) deallocate (Basis%d0gb)
      if (allocated(Basis%d1gb)) deallocate (Basis%d1gb)
      if (allocated(Basis%d2gb)) deallocate (Basis%d2gb)
      if (allocated(Basis%dagb)) deallocate (Basis%dagb)
      if (allocated(Basis%dp0gb)) deallocate (Basis%dp0gb)
      if (allocated(Basis%dq0gb)) deallocate (Basis%dq0gb)

      nb = Basis%nb
      nq = Basis%nq
      dx = pi/nq

      ! grid and weight
      Basis%x = [(dx*(iq - HALF), iq=1, nq)]
      Basis%w = [(dx, iq=1, nq)]

      allocate (Basis%d0gb(nq, nb))
      allocate (Basis%d1gb(nq, nb, 1))
      allocate (Basis%d2gb(nq, nb, 1, 1))


      DO ib = 1, nb
         DO iq = 1, nq
            call d0d1d2box(Basis%x(iq),Basis%d0gb(iq, ib),Basis%d1gb(iq, ib, 1)&
                           &,Basis%d2gb(iq, ib, 1, 1),ib)
         END DO
      END DO

      IF (nb == nq) THEN
         Basis%d0gb(:, nb) = Basis%d0gb(:, nb)/sqrt(TWO)
         Basis%d1gb(:, nb, :) = Basis%d1gb(:, nb, :)/sqrt(TWO)
         Basis%d2gb(:, nb, :, :) = Basis%d2gb(:, nb, :, :)/sqrt(TWO)
      END IF

   END SUBROUTINE Construct_Basis_Sin




   !> Construct the Fourier basis in the interval [-π, π]
   !! This routine allocates the quadrature grid and computes the Fourier
   !! basis functions and their derivatives.
   SUBROUTINE Construct_Basis_Fourier(Basis)
      USE QDUtil_m
      IMPLICIT NONE
   
      TYPE(Basis_t), intent(inout) :: Basis
      real(kind=Rkind)             :: dx
      integer                      :: ib, iq, nb, nq
   
      !---------------------------------------------------------------------------
      ! Deallocate existing arrays if already allocated
      IF (allocated(Basis%x))       deallocate(Basis%x)
      IF (allocated(Basis%w))       deallocate(Basis%w)
      IF (allocated(Basis%d0gb))    deallocate(Basis%d0gb)
      IF (allocated(Basis%d1gb))    deallocate(Basis%d1gb)
      IF (allocated(Basis%d2gb))    deallocate(Basis%d2gb)
      IF (allocated(Basis%dagb))    deallocate(Basis%dagb)
      IF (allocated(Basis%dp0gb))   deallocate(Basis%dp0gb)
      IF (allocated(Basis%dq0gb))   deallocate(Basis%dq0gb)
      !---------------------------------------------------------------------------
   
      nb = Basis%nb
      nq = Basis%nq
   
      ! Step size for the uniform grid over [-π, π]
      dx = TWO * pi / nq
   
      !---------------------------------------------------------------------------
      ! Allocate and initialize grid points and weights
      allocate(Basis%x(nq))
      allocate(Basis%w(nq))
   
      Basis%x = [(iq * dx - dx / 2 - pi, iq = 1, nq)]
      Basis%w = [(dx, iq = 1, nq)]
      !---------------------------------------------------------------------------
   
      !---------------------------------------------------------------------------
      ! Allocate derivative matrices
      allocate(Basis%d0gb(nq, nb))
      allocate(Basis%d1gb(nq, nb, 1))
      allocate(Basis%d2gb(nq, nb, 1, 1))
      !---------------------------------------------------------------------------
   
      !---------------------------------------------------------------------------
      ! Compute the Fourier basis functions and their derivatives
      DO ib = 1, nb
         DO iq = 1, nq
            CALL d0d1d2fourier(Basis%x(iq), &
                               Basis%d0gb(iq, ib), &
                               Basis%d1gb(iq, ib, 1), &
                               Basis%d2gb(iq, ib, 1, 1), ib)
         END DO
      END DO
      !---------------------------------------------------------------------------
   
      !---------------------------------------------------------------------------
      ! Normalize the last basis function if conditions are met
      IF (nb == nq .AND. MOD(nb, 2) == 0) THEN
         Basis%d0gb(:, nb)        = Basis%d0gb(:, nb) / SQRT(TWO)
         Basis%d1gb(:, nb, :)     = Basis%d1gb(:, nb, :) / SQRT(TWO)
         Basis%d2gb(:, nb, :, :)  = Basis%d2gb(:, nb, :, :) / SQRT(TWO)
      END IF
      !---------------------------------------------------------------------------
   
   END SUBROUTINE Construct_Basis_Fourier





    !> @brief Constructs the 'el' (electronic) type basis
   SUBROUTINE Construct_Basis_el(Basis)
      USE QDUtil_m
      IMPLICIT NONE
    
      TYPE(Basis_t), INTENT(inout) :: Basis
    
      ! Set number of grid points to zero for 'el' basis
      Basis%nq = 0
    
      RETURN
    END SUBROUTINE Construct_Basis_el

   !> @brief Constructs the 'ho' (harmonic oscillator) type basis
   !! This subroutine initializes the basis functions and corresponding
   !! grid values for the harmonic oscillator using Hermite polynomials.
   !! It allocates and computes the grid points, weights, and the 0th, 1st, and 2nd
   !! derivatives of the Hermite polynomials at each grid point.
   SUBROUTINE Construct_Basis_Ho(Basis)
      USE QDUtil_m
      IMPLICIT NONE
    
      TYPE(Basis_t), INTENT(inout) :: Basis
      INTEGER                      :: iq, ib
      INTEGER                      :: nb, nq
    
      ! Retrieve the number of basis functions and grid points
      nb = Basis%nb
      nq = Basis%nq
    
      ! Deallocate previously allocated arrays if necessary
      IF (ALLOCATED(Basis%x))       DEALLOCATE(Basis%x)
      IF (ALLOCATED(Basis%w))       DEALLOCATE(Basis%w)
      IF (ALLOCATED(Basis%d0gb))    DEALLOCATE(Basis%d0gb)
      IF (ALLOCATED(Basis%d1gb))    DEALLOCATE(Basis%d1gb)
      IF (ALLOCATED(Basis%d2gb))    DEALLOCATE(Basis%d2gb)
      IF (ALLOCATED(Basis%dagb))    DEALLOCATE(Basis%dagb)
      IF (ALLOCATED(Basis%dp0gb))   DEALLOCATE(Basis%dp0gb)
      IF (ALLOCATED(Basis%dq0gb))   DEALLOCATE(Basis%dq0gb)
    
      ! Allocate grid points and weights
      ALLOCATE(Basis%x(nq))
      ALLOCATE(Basis%w(nq))
      CALL hercom(nq, Basis%x(:), Basis%w(:))  ! Hermite quadrature points and weights
    
      ! Allocate arrays for derivatives
      ALLOCATE(Basis%d0gb(nq, nb))
      ALLOCATE(Basis%d1gb(nq, nb, 1))
      ALLOCATE(Basis%d2gb(nq, nb, 1, 1))
    
      ! Compute Hermite polynomials and their derivatives at each grid point
      DO iq = 1, nq
        DO ib = 1, nb
          CALL d0d1d2poly_Hermite_exp(Basis%x(iq), ib - 1,                        &
               Basis%d0gb(iq, ib), Basis%d1gb(iq, ib, 1), Basis%d2gb(iq, ib, 1, 1), .TRUE.)
        END DO
      END DO
    
    END SUBROUTINE Construct_Basis_Ho


   SUBROUTINE Construct_Basis_Hag(Basis)
      USE QDUtil_m
      TYPE(Basis_t), intent(inout) :: Basis
  
      integer :: iq, ib
      integer :: nb, nq
      real(kind=Rkind) :: q0, sq, p0, a, x
      complex(kind=Rkind) :: alpha
  
      !----------------------------------------------------------
      ! Retrieve basis parameters from the Basis structure
      !----------------------------------------------------------
      q0    = Basis%Q0
      sq    = Basis%scaleQ
      alpha = Basis%alpha
      p0    = Basis%imp_k
  
      nb = Basis%nb
      nq = Basis%nq
  
      a = REAL(alpha, kind=Rkind)  ! Re(alpha)
  
      !----------------------------------------------------------
      ! Deallocate previous basis arrays if allocated
      !----------------------------------------------------------
      if (allocated(Basis%x)) deallocate (Basis%x)
      if (allocated(Basis%w)) deallocate (Basis%w)
      if (allocated(Basis%d0gb)) deallocate (Basis%d0gb)
      if (allocated(Basis%d1gb)) deallocate (Basis%d1gb)
      if (allocated(Basis%d2gb)) deallocate (Basis%d2gb)
  
      !----------------------------------------------------------
      ! Allocate grids and quadrature weights (Hermite quadrature)
      !----------------------------------------------------------
      allocate (Basis%x(nq))
      allocate (Basis%w(nq))
      call hercom(nq, Basis%x(:), Basis%w(:))
  
      !----------------------------------------------------------
      ! Allocate basis function arrays
      !----------------------------------------------------------
      allocate (Basis%d0gb(nq, nb))
      allocate (Basis%d1gb(nq, nb, 1))
      allocate (Basis%d2gb(nq, nb, 1, 1))
  
      !----------------------------------------------------------
      ! Construct Hagedorn basis functions directly
      !----------------------------------------------------------
      DO iq = 1, nq
          DO ib = 1, nb
              call d0d1d2poly_Hag_exp(Basis%x(iq), ib - 1, &
                   Basis%d0gb(iq, ib), Basis%d1gb(iq, ib, 1), Basis%d2gb(iq, ib, 1, 1), &
                   .true., alpha, p0)
          END DO
      END DO
  
  END SUBROUTINE Construct_Basis_Hag


   !> Check orthonormality of the basis set
   !! This subroutine checks if the basis functions are orthonormal (⟨b_i | b_j⟩ ≈ δ_ij)
   !! It also checks the overlaps with derivatives if requested (up to second order).
   SUBROUTINE CheckOrtho_Basis(Basis, nderiv)
      USE QDUtil_m
      IMPLICIT NONE
   
      TYPE(Basis_t), intent(in)          :: Basis
      integer, intent(in)                :: nderiv
      integer                            :: ib, jb
      real(kind=Rkind), ALLOCATABLE      :: S(:, :)
      real(kind=Rkind)                   :: Sii, Sij
   
      !-------------------------------------------------------------------------
      ! Check for incompatible basis type
      IF (Basis%Basis_name == 'el') THEN
         PRINT *, 'This routine is not applicable to Basis ''el'''
         RETURN
      END IF
      !-------------------------------------------------------------------------
   
      !-------------------------------------------------------------------------
      ! Proceed only if basis data is allocated
      IF (Basis_IS_allocated(Basis)) THEN
   
         ! Compute overlap matrix: S = <d0b|d0b>
         S = matmul(conjg(Basis%d0bgw), Basis%d0gb)
   
         ! Optional debug loop (currently disabled)
         ! DO ib = 1, Basis%nb
         !    DO jb = 1, Basis%nb
         !       WRITE(115,*) ib, jb, S(ib, ib), S(ib, jb)
         !    END DO
         ! END DO
   
         ! Check deviation from orthonormality
         Sii = ZERO
         DO ib = 1, Basis%nb
            IF (ABS(S(ib, ib) - ONE) > Sii) Sii = ABS(S(ib, ib) - ONE)
            S(ib, ib) = ZERO
         END DO
         Sij = MAXVAL(ABS(S))
   
         WRITE(out_unit, *) 'Deviation from orthonormality: max|Sii - 1| =', Sii, ', max|Sij| =', Sij
   
         !-------------------------------------------------------------------------
         ! Check overlap with first derivatives if required
         IF (nderiv > 0) THEN
            WRITE(out_unit, *)
            S = matmul(conjg(Basis%d0bgw), Basis%d1gb(:, :, 1))
            ! CALL Write_VecMat(S, out_unit, 5, info='<d0b|d1b>')
         END IF
   
         !-------------------------------------------------------------------------
         ! Check overlap with second derivatives if required
         IF (nderiv > 1) THEN
            WRITE(out_unit, *)
            S = matmul(conjg(Basis%d0bgw), Basis%d2gb(:, :, 1, 1))
            ! CALL Write_VecMat(S, out_unit, 5, info='<d0b|d2b>')
         END IF
   
      ELSE
         WRITE(out_unit, *) 'WARNING in CheckOrtho_Basis: The basis is not allocated.'
      END IF
   
   END SUBROUTINE CheckOrtho_Basis
   SUBROUTINE Scale_Basis(Basis, q0, sq)
      USE QDUtil_m
      TYPE(Basis_t), intent(inout)  :: Basis
      real(kind=Rkind), intent(in)  :: q0, sq
   
      !----------------------------------------------------------------------
      ! SUBROUTINE Scale_Basis(Basis, q0, sq)
      !  --> Perform the change of variable: x = sqrt(a)*(q - q0)
      !      with a = Re(alpha), sq = sqrt(a) provided directly as input.
      !
      ! Inputs:
      !   Basis : type(Basis_t) containing grid points x, weights w,
      !            basis functions d0gb, first derivatives d1gb, second derivatives d2gb
      !   q0    : translation point q0
      !   sq    : sqrt(a), scaling factor given as input (do not recompute a)
      !
      ! The subroutine applies:
      !   - Translation and rescaling of points x -> q
      !   - Adjusts integration weights w
      !   - Applies normalization factors using sq directly:
      !         a^{1/4}  = sqrt(sq)
      !         a^{3/4}  = sq * sqrt(sq)
      !         a^{5/4}  = sq^2 * sqrt(sq)
      !----------------------------------------------------------------------
   
      REAL(KIND=Rkind) :: sq_half, sq_3half, sq_5half
   
      IF (Basis%nq == 0) RETURN
   
      IF (abs(sq) > ONETENTH**6 .AND. Basis_IS_allocated(Basis)) THEN
         sq_half  = sqrt(sq)
         sq_3half = sq * sq_half
         sq_5half = sq * sq * sq_half
   
         ! Transform grid points from x to q
         Basis%x(:) = q0 + Basis%x(:)/sq
         ! Adjust integration weights
         Basis%w(:) = Basis%w(:)/sq
   
         ! Apply normalization factor a^{1/4} = sqrt(sq)
         Basis%d0gb(:, :) = Basis%d0gb(:, :) * sq_half
   
         ! Scale derivatives accordingly
         Basis%d1gb(:, :, :) = Basis%d1gb(:, :, :) * sq_3half
         Basis%d2gb(:, :, :, :) = Basis%d2gb(:, :, :, :) * sq_5half
      ELSE
         write (out_unit, *) ' ERROR in Scale_Basis'
         write (out_unit, *) ' sq is too small or the basis is not allocated.'
         STOP 'ERROR in Scale_Basis'
      END IF
   END SUBROUTINE Scale_Basis

   SUBROUTINE Complete_construct_Basis(Basis)
      USE QDUtil_m
      TYPE(Basis_t), intent(inout) :: Basis
      real(kind=Rkind) :: Q0, SQ0, a
      complex(kind=Rkind) :: alpha
  
      !--------------------------------------------------------------
      ! Retrieve alpha and compute scaling factor SQ0 = sqrt(Re(alpha))
      !--------------------------------------------------------------
      alpha = Basis%alpha
      a     = REAL(alpha, kind=Rkind)
      SQ0   = SQRT(a)
  
      !--------------------------------------------------------------
      ! Retrieve the translation center Q0
      !--------------------------------------------------------------
      Q0 = Basis%Q0
  
      !--------------------------------------------------------------
      ! Apply scaling transformation: x = SQ0 * (q - Q0)
      !--------------------------------------------------------------
      call Scale_Basis(Basis, Q0, SQ0)
      !--------------------------------------------------------------
      ! Compute the transpose of d0gb (for efficiency in projections)
      !--------------------------------------------------------------
      call Calc_tranpose_d0gb(Basis)
      !--------------------------------------------------------------
      call Calc_dngg_grid_0(Basis)
      !--------------------------------------------------------------
      ! Check orthonormality of the constructed basis up to 2nd derivative
      !--------------------------------------------------------------
      call CheckOrtho_Basis(Basis, nderiv=2)
  
  END SUBROUTINE 
 

    !> Compute the transpose of d0gb and apply the quadrature weights
    !! This subroutine computes the weighted transpose of the d0gb 
    !! matrix and stores it in d0bgw.
    SUBROUTINE Calc_tranpose_d0gb(Basis)
       USE QDUtil_m
       IMPLICIT NONE
    
       TYPE(Basis_t), intent(inout)      :: Basis
       INTEGER                           :: ib, nb
    
       !-------------------------------------------------------------------------
       ! Get the number of basis functions
       nb = Basis%nb
    
       !-------------------------------------------------------------------------
       ! If the output array already exists, deallocate it
       IF (allocated(Basis%d0bgw)) DEALLOCATE(Basis%d0bgw)
    
       !-------------------------------------------------------------------------
       ! Special case: skip computation if Basis is of type 'el'
       IF (Basis%Basis_name == 'el') RETURN
    
       !-------------------------------------------------------------------------
       ! Compute the transpose and apply weights
       Basis%d0bgw = TRANSPOSE(Basis%d0gb)
    
       DO ib = 1, nb
          Basis%d0bgw(ib, :) = Basis%d0bgw(ib, :) * Basis%w(:)
       END DO
    
    END SUBROUTINE Calc_tranpose_d0gb


    SUBROUTINE Calc_Basis_dPtQtAt(Basis)

    USE QDUtil_m
    TYPE(Basis_t), intent(inout)    :: Basis
    integer                         :: ib,iq,nb,nq
    !logical, parameter             :: debug = .true.
    logical, parameter              :: debug = .false.
    real(kind=Rkind)                :: a,Q0 


    Q0=Basis%Q0
    a = Basis%SCALEQ
    nb = Basis%nb
    nq = Basis%nq

    IF (debug) THEN
       write (out_unit, *) 'BEGINNING Calc_dngg_grid'
       CALL Write_Basis(Basis)
       flush (out_unit)
    END IF

    if (allocated(Basis%dagb)) deallocate (Basis%dagb)
    if (allocated(Basis%dq0gb)) deallocate (Basis%dq0gb)
    if (allocated(Basis%dp0gb)) deallocate (Basis%dp0gb)
     
    allocate (Basis%dagb(nq, nb))
    allocate (Basis%dq0gb(nq, nb))
    allocate (Basis%dp0gb(nq, nb))

   if (Basis%Basis_name == 'herm' .or. Basis%Basis_name == 'ho') then

    do ib = 1,nb
       do iq = 1,nq
          Basis%dagb(iq, ib)  = (ONE/a)*(HALF*Basis%d0gb(iq, ib)+(Basis%x(iq)-Q0)*Basis%d1gb(iq, ib, 1))
          Basis%dq0gb(iq, ib) = -Basis%d1gb(iq, ib, 1)
          Basis%dp0gb(iq, ib) = EYE*(Basis%x(iq)-Q0)*Basis%d0gb(iq, ib)
       end do
    end do

    Elseif(Basis%Basis_name == 'boxab') then
   do ib = 1,nb
     do iq = 1,nq
       Basis%dagb(iq, ib)  = (ONE/a)*(HALF*Basis%d0gb(iq, ib)+(Basis%x(iq)-Q0)*Basis%d1gb(iq, ib, 1))
       Basis%dq0gb(iq, ib) = -Basis%d1gb(iq, ib, 1)
       Basis%dp0gb(iq, ib) =CZERO
     end do
   end do
     Elseif(Basis%Basis_name == 'fourier') then

     Basis%dagb(:, :)   = CZERO
     Basis%dq0gb(:,:) = CZERO
     Basis%dp0gb(:,:) = CZERO

     End If

 END SUBROUTINE

    !> Computes several grid-based derivative matrices in the gg (grid-grid) representation
   !! using the precomputed gb (grid-basis) data and the weighted conjugate transpose (d0bgw).
   SUBROUTINE Calc_dngg_grid(Basis)
      USE QDUtil_m
      IMPLICIT NONE
   
      TYPE(Basis_t), intent(inout) :: Basis
      INTEGER                      :: ib, iq, nb, nq
      LOGICAL, PARAMETER           :: debug = .false.  ! Set to true for debugging output
   
      !-------------------------------------------------------------------------
      ! Retrieve number of basis functions and grid points
      nb = Basis%nb
      nq = Basis%nq
   
      !-------------------------------------------------------------------------
      ! Debug output (if enabled)
      IF (debug) THEN
         WRITE(out_unit, *) 'BEGINNING Calc_dngg_grid'
         CALL Write_Basis(Basis)
         FLUSH(out_unit)
      END IF
   
      !-------------------------------------------------------------------------
      ! Deallocate previous gg matrices if already allocated
      IF (allocated(Basis%d1gg)) DEALLOCATE(Basis%d1gg)
      IF (allocated(Basis%d2gg)) DEALLOCATE(Basis%d2gg)
      IF (allocated(Basis%dagg)) DEALLOCATE(Basis%dagg)
      IF (allocated(Basis%dq0gg)) DEALLOCATE(Basis%dq0gg)
      IF (allocated(Basis%dp0gg)) DEALLOCATE(Basis%dp0gg)
   
      !-------------------------------------------------------------------------
      ! Allocate matrices for gg representations
      ALLOCATE(Basis%d1gg(nq, nq, 1))
      ALLOCATE(Basis%d2gg(nq, nq, 1, 1))
      ALLOCATE(Basis%dagg(nq, nq))
      ALLOCATE(Basis%dq0gg(nq, nq))
      ALLOCATE(Basis%dp0gg(nq, nq))
   
      !-------------------------------------------------------------------------
      ! Optional debug output: show d0bgw matrix
      IF (debug) THEN
         CALL Write_VecMat(Basis%d0bgw, out_unit, 5, info='d0bgw')
         WRITE(out_unit, *)
      END IF
   
      !-------------------------------------------------------------------------
      ! Skip computation for 'el' basis type
      IF (Basis%Basis_name == 'el') RETURN
   
      !-------------------------------------------------------------------------
      ! Compute grid-grid matrices via matrix multiplications
      Basis%d1gg(:, :, 1)     = MATMUL(Basis%d1gb(:, :, 1), CONJG(Basis%d0bgw))
      Basis%d2gg(:, :, 1, 1)  = MATMUL(Basis%d2gb(:, :, 1, 1), CONJG(Basis%d0bgw))
      Basis%dagg              = MATMUL(Basis%dagb, CONJG(Basis%d0bgw))
      Basis%dq0gg             = MATMUL(Basis%dq0gb, CONJG(Basis%d0bgw))
      Basis%dp0gg             = MATMUL(Basis%dp0gb, CONJG(Basis%d0bgw))
   
      !-------------------------------------------------------------------------
      ! Debug output: print calculated matrices
      IF (debug) THEN
         CALL Write_VecMat(Basis%d1gg(:, :, 1), out_unit, 5, info='d1gg')
         WRITE(out_unit, *)
         CALL Write_VecMat(Basis%d2gg(:, :, 1, 1), out_unit, 5, info='d2gg')
         WRITE(out_unit, *)
         CALL Write_VecMat(Basis%dp0gg, out_unit, 5, info='dp0gg')
         WRITE(out_unit, *)
         CALL Write_VecMat(Basis%dq0gg, out_unit, 5, info='dq0gg')
         WRITE(out_unit, *)
         CALL Write_VecMat(Basis%dagg, out_unit, 5, info='dagg')
         WRITE(out_unit, *) 'END Calc_dngg_grid'
         FLUSH(out_unit)
      END IF
   
   END SUBROUTINE Calc_dngg_grid
   
   
   !> Computes only the first and second derivative operators (d1gg and d2gg)
   !! in the grid-grid representation, using the basis-grid derivative data.
   SUBROUTINE Calc_dngg_grid_0(Basis)
      USE QDUtil_m
      IMPLICIT NONE
   
      TYPE(Basis_t), intent(inout) :: Basis
      INTEGER                      :: ib, iq, nb, nq
      LOGICAL, PARAMETER           :: debug = .false.  ! Set to true for debug printing
   
      !-------------------------------------------------------------------------
      ! Retrieve the number of basis functions and grid points
      nb = Basis%nb
      nq = Basis%nq
   
      !-------------------------------------------------------------------------
      ! Optional debug print of basis setup
      IF (debug) THEN
         WRITE(out_unit, *) 'BEGINNING Calc_dngg_grid_0'
         CALL Write_Basis(Basis)
         FLUSH(out_unit)
      END IF
   
      !-------------------------------------------------------------------------
      ! Deallocate previous versions of d1gg and d2gg if already allocated
      IF (ALLOCATED(Basis%d1gg)) DEALLOCATE(Basis%d1gg)
      IF (ALLOCATED(Basis%d2gg)) DEALLOCATE(Basis%d2gg)
   
      !-------------------------------------------------------------------------
      ! Allocate arrays for the derivative operators in the grid representation
      ALLOCATE(Basis%d1gg(nq, nq, 1))
      ALLOCATE(Basis%d2gg(nq, nq, 1, 1))
   
      !-------------------------------------------------------------------------
      ! Optional debug print of basis-weighted transpose
      IF (debug) THEN
         CALL Write_VecMat(Basis%d0bgw, out_unit, 5, info='d0bgw')
         WRITE(out_unit, *)
      END IF
   
      !-------------------------------------------------------------------------
      ! Skip computation if basis is of type 'el'
      IF (Basis%Basis_name == 'el') RETURN
   
      !-------------------------------------------------------------------------
      ! Compute derivative matrices in grid-grid representation
      Basis%d1gg(:, :, 1) = MATMUL(Basis%d1gb(:, :, 1), CONJG(Basis%d0bgw))
      Basis%d2gg(:, :, 1, 1) = MATMUL(Basis%d2gb(:, :, 1, 1), CONJG(Basis%d0bgw))
   
      !-------------------------------------------------------------------------
      ! Debug output of resulting matrices
      IF (debug) THEN
         CALL Write_VecMat(Basis%d1gg(:, :, 1), out_unit, 5, info='d1gg')
         WRITE(out_unit, *)
         CALL Write_VecMat(Basis%d2gg(:, :, 1, 1), out_unit, 5, info='d2gg')
         WRITE(out_unit, *)
         WRITE(out_unit, *) 'END Calc_dngg_grid_0'
         FLUSH(out_unit)
      END IF
   
   END SUBROUTINE Calc_dngg_grid_0

   SUBROUTINE test_basitogridgridtobasis(Basis)
      USE QDUtil_m
      TYPE(Basis_t), intent(in)       :: Basis
      logical, parameter                 :: debug = .true.
      complex(kind=Rkind), allocatable   :: G1(:), B1(:)!,G(:), Hpsi(:)
      complex(kind=Rkind), allocatable   :: G2(:), B2(:)
      complex(kind=Rkind), allocatable   :: B(:)
      real(kind=Rkind), allocatable       :: diff_g(:), diff_b(:)
      !REAL(KIND=Rkind)                   :: Norm0,Norm1,min_diff,max_diff
      integer                             :: iq, ndim
      ndim = size(Basis%tab_basis)

      IF (debug) THEN
         write (out_unit, *) 'BEGINNING Test'
         flush (out_unit)
      END IF

      allocate (B(Basis%nb*Basis%tab_basis(ndim)%nb))
      allocate (G1(Basis%nq*Basis%tab_basis(ndim)%nb))
      allocate (diff_g(Basis%nq*Basis%tab_basis(ndim)%nb))
      allocate (diff_b(Basis%nb*Basis%tab_basis(ndim)%nb))
      allocate (G2(Basis%nq*Basis%tab_basis(ndim)%nb))
      allocate (B1(Basis%nb*Basis%tab_basis(ndim)%nb))
      allocate (B2(Basis%nb*Basis%tab_basis(ndim)%nb))

      B(:) = CZERO
      B1(:) = CONE
      G1(:) = CONE
      G2(:) = CZERO
      !print *, 'G in', G1
      Call GridTOBasis_nD_cplx(B, G1, Basis)
      Call BasisTOGrid_nD_cplx(G2, B, Basis)
      !print *, 'G out', G2

      print *, '---------------------------------------------------------------------'

      diff_g(:) = ABS(G1(:) - G2(:))

      print *, '---------------------------------------------------------------------'
      Write (out_unit, *) 'maxval(diff_g(:))=', maxval(diff_g(:))
      Write (out_unit, *) 'MINVAL(diff_g(:))=', MINVAL(diff_g(:))
      print *, '---------------------------------------------------------------------'
      G2(:) = CZERO
      B2(:) = CZERO
      B1(:) = CONE
      !print *, 'B in', B1
      Call BasisTOGrid_nD_cplx(G2, B1, Basis)
      Call GridTOBasis_nD_cplx(B2, G2, Basis)
      !print *, 'B out', B2
      diff_b(:) = ABS(B2(:) - B1(:))
      print *, '---------------------------------------------------------------------'
      Write (out_unit, *) 'maxval(diff_b(:))=', maxval(diff_b(:))
      Write (out_unit, *) 'MINVAL(diff_b(:))=', MINVAL(diff_b(:))
      print *, '---------------------------------------------------------------------'

      IF (debug) THEN
         write (out_unit, *) 'END Test'
         flush (out_unit)
      END IF

   END SUBROUTINE test_basitogridgridtobasis

   !> Transform a function from grid representation (GG) to basis representation (BB)
   !! in one dimension, using complex arithmetic.
   !! The transformation is: BB = <basis|grid> * GG
   SUBROUTINE GridTOBasis_1D_cplx(BB, GG, Basis)
      USE QDUtil_m
      IMPLICIT NONE
   
      TYPE(Basis_t), intent(in), target        :: Basis
      COMPLEX(kind=Rkind), intent(inout)       :: BB(:, :, :)  ! Output: Basis representation
      COMPLEX(kind=Rkind), intent(in)          :: GG(:, :, :)  ! Input: Grid representation
      LOGICAL, PARAMETER                       :: debug = .true.
      INTEGER                                  :: i1, i3
   
      !-------------------------------------------------------------------------
      ! Optional debug flushing
      IF (debug) THEN
         FLUSH(out_unit)
      END IF
   
      !-------------------------------------------------------------------------
      ! Initialize the BB array to zero
      BB = CZERO
   
      !-------------------------------------------------------------------------
      ! Loop over all dimensions (assumed to be 1D data with extra indices)
      DO i3 = 1, UBOUND(GG, dim=3)
         DO i1 = 1, UBOUND(GG, dim=1)
            ! Matrix multiplication: BB(i1,:,i3) = <basis|grid> * GG(i1,:,i3)
            BB(i1, :, i3) = MATMUL(CONJG(Basis%d0bgw), GG(i1, :, i3))
         END DO
      END DO
   
      !-------------------------------------------------------------------------
      ! Optional debug flushing
      IF (debug) THEN
         FLUSH(out_unit)
      END IF
   
   END SUBROUTINE GridTOBasis_1D_cplx

   !> Transform a function from basis representation (BB) to grid representation (GB)
   !! in one dimension, using complex arithmetic.
   !! The transformation is: GB = <x|basis> * BB
   SUBROUTINE BasisTOGrid_1D_cplx(GB, BB, Basis)
      USE QDUtil_m
      IMPLICIT NONE
   
      TYPE(Basis_t), intent(in), target         :: Basis
      COMPLEX(kind=Rkind), intent(inout)        :: GB(:, :, :)  ! Output: Grid representation
      COMPLEX(kind=Rkind), intent(in)           :: BB(:, :, :)  ! Input: Basis representation
      LOGICAL, PARAMETER                        :: debug = .true.
      INTEGER                                   :: i1, i3
   
      !-------------------------------------------------------------------------
      ! Optional debug flushing before processing
      IF (debug) THEN
         FLUSH(out_unit)
      END IF
   
      !-------------------------------------------------------------------------
      ! Initialize the GB array to zero
      GB = CZERO
   
      !-------------------------------------------------------------------------
      ! Loop over all dimensions (i1 and i3 are loop indices for external dimensions)
      DO i3 = 1, UBOUND(BB, dim=3)
         DO i1 = 1, UBOUND(BB, dim=1)
            ! Matrix multiplication:
            ! GB(i1,:,i3) = <x|basis> * BB(i1,:,i3)
            ! d0gb contains basis functions evaluated on the grid: φ_b(x_q)
            GB(i1, :, i3) = MATMUL(Basis%d0gb(:, :), BB(i1, :, i3))
         END DO
      END DO
   
      !-------------------------------------------------------------------------
      ! Optional debug flushing after processing
      IF (debug) THEN
         FLUSH(out_unit)
      END IF
   
   END SUBROUTINE BasisTOGrid_1D_cplx



   !> Calculates indexing parameters (Ib1, Ib2, Ib3, Iq1, Iq2, Iq3) 
   !! used for reshaping and accessing multidimensional basis/grid data.
   !! This subroutine is typically used in multidimensional simulations (e.g., 3D quantum grids).
   SUBROUTINE Calc_index(Ib1, Ib2, Ib3, Iq1, Iq2, Iq3, Basis)
   TYPE(Basis_t), intent(in), target              :: Basis
   INTEGER, intent(inout), allocatable, optional  :: Iq1(:), Iq2(:), Iq3(:)
   INTEGER, intent(inout), allocatable, optional  :: Ib1(:), Ib2(:), Ib3(:)
   INTEGER                                        :: Ndim      ! Number of spatial dimensions
   INTEGER                                        :: inb       ! Loop index

   ! Determine the number of spatial dimensions (tab_basis has Ndim+1 entries)
   Ndim = SIZE(Basis%tab_basis) - 1

   ! Allocate index arrays if they are present in the arguments
   IF (PRESENT(Ib3)) ALLOCATE(Ib3(Ndim))
   IF (PRESENT(Ib2)) ALLOCATE(Ib2(Ndim))
   IF (PRESENT(Ib1)) ALLOCATE(Ib1(Ndim))
   IF (PRESENT(Iq3)) ALLOCATE(Iq3(Ndim))
   IF (PRESENT(Iq2)) ALLOCATE(Iq2(Ndim))
   IF (PRESENT(Iq1)) ALLOCATE(Iq1(Ndim))

   ! Loop over each dimension to calculate indexing terms
   DO inb = 1, Ndim
      IF (inb == 1) THEN
         ! First dimension
         IF (PRESENT(Iq1)) Iq1(1) = 1
         IF (PRESENT(Ib1)) Ib1(1) = 1
         IF (PRESENT(Iq2)) Iq2(1) = Basis%tab_basis(1)%nq
         IF (PRESENT(Ib2)) Ib2(1) = Basis%tab_basis(1)%nb
         IF (PRESENT(Iq3)) Iq3(1) = PRODUCT(Basis%tab_basis(2:Ndim)%nq) * Basis%tab_basis(Ndim + 1)%nb
         IF (PRESENT(Ib3)) Ib3(1) = PRODUCT(Basis%tab_basis(2:Ndim + 1)%nb)

      ELSE IF (inb == Ndim) THEN
         ! Last dimension
         IF (PRESENT(Iq1)) Iq1(inb) = PRODUCT(Basis%tab_basis(1:Ndim - 1)%nq)
         IF (PRESENT(Ib1)) Ib1(inb) = PRODUCT(Basis%tab_basis(1:Ndim - 1)%nb)
         IF (PRESENT(Iq2)) Iq2(inb) = Basis%tab_basis(Ndim)%nq
         IF (PRESENT(Ib2)) Ib2(inb) = Basis%tab_basis(Ndim)%nb
         IF (PRESENT(Iq3)) Iq3(inb) = Basis%tab_basis(Ndim + 1)%nb
         IF (PRESENT(Ib3)) Ib3(inb) = Basis%tab_basis(Ndim + 1)%nb

      ELSE
         ! Intermediate dimensions
         IF (PRESENT(Iq1)) Iq1(inb) = PRODUCT(Basis%tab_basis(1:inb - 1)%nq)
         IF (PRESENT(Ib1)) Ib1(inb) = PRODUCT(Basis%tab_basis(1:inb - 1)%nb)
         IF (PRESENT(Iq2)) Iq2(inb) = Basis%tab_basis(inb)%nq
         IF (PRESENT(Ib2)) Ib2(inb) = Basis%tab_basis(inb)%nb
         IF (PRESENT(Iq3)) Iq3(inb) = PRODUCT(Basis%tab_basis(inb + 1:Ndim)%nq) * Basis%tab_basis(Ndim + 1)%nb
         IF (PRESENT(Ib3)) Ib3(inb) = PRODUCT(Basis%tab_basis(inb + 1:Ndim + 1)%nb)
      END IF
   END DO

END SUBROUTINE Calc_index

SUBROUTINE BasisTOGrid_nD_cplx(G, B, Basis)
   USE QDUtil_m
   USE NDindex_m

   ! Toggle debugging output
   LOGICAL, PARAMETER :: debug = .false.

   ! Input/Output arguments
   TYPE(Basis_t), INTENT(IN)                  :: Basis          ! Basis definition and metadata
   COMPLEX(KIND=Rkind), INTENT(IN), TARGET    :: B(:)           ! Input in Basis representation
   COMPLEX(KIND=Rkind), INTENT(INOUT), TARGET :: G(:)           ! Output in Grid representation

   ! Pointers for reshaped views of arrays
   COMPLEX(KIND=Rkind), POINTER :: BBG(:, :, :), BBB(:, :, :), GGB(:, :, :)
   COMPLEX(KIND=Rkind), POINTER :: GBB(:, :, :)

   ! Indexing variables (flattening multidimensional data)
   INTEGER, ALLOCATABLE :: Ib3, Iq1, Iq2, Ib1, Ib2, Iq3

   ! Temporary working arrays
   COMPLEX(KIND=Rkind), ALLOCATABLE, TARGET :: GBB1(:), GGB2(:)

   ! Integer variables
   INTEGER :: Ib, Iq, nq, nb, inb, ndim
   INTEGER :: jb, jb1, jb2

   ! Debugging: print inputs and basis
   IF (debug) THEN
       WRITE(out_unit, *) 'BEGINNING BasisTOGrid_nD_cplx'
       WRITE(out_unit, *) 'intent(in) :: B(:)', B
       CALL Write_Basis(Basis)
       FLUSH(out_unit)
   END IF

   ! Check that basis is allocated
   IF (.NOT. Basis_IS_Allocated(Basis)) THEN
       WRITE(out_unit, *) ' ERROR: the basis is not allocated.'
       STOP 'ERROR BasisTOGrid_Basis: the basis is not allocated.'
   END IF

   ! Check that the sizes of B and G are compatible with the basis structure
   IF (SIZE(B) /= Basis%nb * Basis%tab_basis(SIZE(Basis%tab_basis))%nb) THEN
       WRITE(out_unit, *) ' ERROR: the size of B does not match expected nb.'
       STOP 'ERROR in BasisTOGrid_Basis: wrong B size.'
   END IF
   IF (SIZE(G) /= Basis%nq * Basis%tab_basis(SIZE(Basis%tab_basis))%nb) THEN
       WRITE(out_unit, *) ' ERROR: the size of G does not match expected nq.'
       STOP 'ERROR in BasisTOGrid_Basis: wrong G size.'
   END IF

   ! Determine number of dimensions (excluding electronic dimension)
   ndim = SIZE(Basis%tab_basis) - 1

   ! --- Special Case: 1D Transform ---
   IF (ndim == 1) THEN
       ! Set index factors for single-dimension case
       Iq1 = 1
       Iq2 = Basis%tab_basis(1)%nq
       Iq3 = Basis%tab_basis(ndim + 1)%nb

       Ib1 = 1
       Ib2 = Basis%tab_basis(1)%nb
       Ib3 = PRODUCT(Basis%tab_basis(2:ndim + 1)%nb)

       ! Reshape input B and output G into 3D slices
       BBB(1:Ib1, 1:Ib2, 1:Ib3) => B
       GGB(1:Iq1, 1:Iq2, 1:Iq3) => G

       ! Initialize output
       G = CZERO

       ! Perform transformation for 1D case
       CALL BasisTOGrid_1D_cplx(GGB, BBB, Basis%tab_basis(1))

   ELSE
   ! --- General Case: Multi-Dimensional Transform ---
       ! Initialize first transformation (dimension 1)
       Iq1 = 1
       Iq2 = Basis%tab_basis(1)%nq
       Iq3 = Basis%tab_basis(ndim + 1)%nb

       Ib1 = 1
       Ib2 = Basis%tab_basis(1)%nb
       Ib3 = PRODUCT(Basis%tab_basis(2:ndim + 1)%nb)

       ALLOCATE(GBB1(Iq1 * Iq2 * Ib3))
       BBB(1:Iq1, 1:Ib2, 1:Ib3) => B
       GBB(1:Iq1, 1:Iq2, 1:Ib3) => GBB1
       GBB1 = CZERO

       ! First 1D transformation
       CALL BasisTOGrid_1D_cplx(GBB, BBB, Basis%tab_basis(1))

       ! Loop over intermediate dimensions (2 to ndim-1)
       DO inb = 2, ndim - 1
           Iq1 = PRODUCT(Basis%tab_basis(1:inb - 1)%nq)
           Ib1 = PRODUCT(Basis%tab_basis(1:inb - 1)%nb)
           Iq2 = Basis%tab_basis(inb)%nq
           Ib2 = Basis%tab_basis(inb)%nb
           Ib3 = PRODUCT(Basis%tab_basis(inb + 1:ndim + 1)%nb)
           Iq3 = PRODUCT(Basis%tab_basis(inb + 1:ndim)%nq) * Basis%tab_basis(ndim + 1)%nb

           ! Allocate temporary output
           ALLOCATE(GGB2(Iq1 * Iq2 * Ib3))

           ! Set reshaped views for input/output
           BBB(1:Iq1, 1:Ib2, 1:Ib3) => GBB1
           GBB(1:Iq1, 1:Iq2, 1:Ib3) => GGB2

           ! Perform 1D transformation for dimension `inb`
           CALL BasisTOGrid_1D_cplx(GBB, BBB, Basis%tab_basis(inb))

           ! Move result to input buffer for next iteration
           GBB1 = GGB2
           DEALLOCATE(GGB2)
       END DO

       ! Final transformation for last dimension
       Iq1 = PRODUCT(Basis%tab_basis(1:ndim - 1)%nq)
       Ib1 = PRODUCT(Basis%tab_basis(1:ndim - 1)%nb)
       Ib2 = Basis%tab_basis(ndim)%nb
       Iq2 = Basis%tab_basis(ndim)%nq
       Ib3 = Basis%tab_basis(ndim + 1)%nb
       Iq3 = Basis%tab_basis(ndim + 1)%nb

       BBB(1:Iq1, 1:Ib2, 1:Ib3) => GBB1
       GBB(1:Iq1, 1:Iq2, 1:Ib3) => G

       CALL BasisTOGrid_1D_cplx(GBB, BBB, Basis%tab_basis(ndim))
   END IF

   ! Debugging: output result
   IF (debug) THEN
       WRITE(out_unit, *) 'intent(INOUT) :: G(:)', G
       WRITE(out_unit, *) 'END BasisTOGrid_nD_cplx'
       FLUSH(out_unit)
   END IF

END SUBROUTINE

SUBROUTINE GridTOBasis_nD_cplx(B, G, Basis)
   USE QDUtil_m
   Logical, parameter                               :: debug = .false.
   TYPE(Basis_t), intent(in), target                :: Basis
   complex(kind=Rkind), intent(in), target          :: G(:)         ! Input: Grid representation
   complex(kind=Rkind), intent(inout), target       :: B(:)         ! Output: Basis representation
   complex(kind=Rkind), pointer                     :: BBB(:, :, :), GGB(:, :, :)
   complex(kind=Rkind), allocatable, target         :: BGG1(:), BGG2(:)
   complex(kind=Rkind), pointer                     :: GGG(:, :, :)
   Integer                                          :: ib, i1, i3, inb, Ndim, iq
   Integer                                          :: Ib1, Ib2, Iq3, Iq1, Iq2, Ib3

   ! Optional debug output
   IF (debug) THEN
      write (out_unit, *) 'BEGINNING GridTOBasis_nD_cplx'
      write (out_unit, *) 'intent(in) :: G(:)', G
      flush (out_unit)
   END IF

   ! Check that basis is properly allocated
   IF (.NOT. Basis_IS_Allocated(Basis)) THEN
      write (out_unit, *) ' ERROR in BasisTOGrid_nD_cplx'
      STOP "ERROR BasisTOGrid_Basis: the basis is not Allocated."
   END IF

   ! Check input and output array sizes
   IF (size(B) /= Basis%nb * Basis%tab_basis(size(Basis%tab_basis))%nb) THEN
      write (out_unit, *) ' ERROR: Wrong size for B'
      STOP 'ERROR in GridTOBasis_Basis: wrong B size.'
   END IF

   IF (size(G) /= Basis%nq * Basis%tab_basis(size(Basis%tab_basis))%nb) THEN
      write (out_unit, *) ' ERROR: Wrong size for G'
      STOP 'ERROR in GridTOBasis_Basis: wrong G size'
   END IF

   ndim = size(Basis%tab_basis) - 1   ! Number of spatial dimensions

   ! Special case: 1D
   IF (ndim == 1) THEN
      Iq1 = 1
      Iq2 = Basis%tab_basis(1)%nq
      Iq3 = Basis%tab_basis(ndim + 1)%nb

      Ib1 = 1
      Ib2 = Basis%tab_basis(1)%nb
      Ib3 = Product(Basis%tab_basis(2:ndim+1)%nb)

      B = CZERO
      GGB(1:Iq1, 1:Iq2, 1:Iq3) => G
      BBB(1:Ib1, 1:Ib2, 1:Ib3) => B
      Call GridTOBasis_1D_cplx(BBB, GGB, Basis%tab_basis(1))

   ELSE
      ! Initialization for first dimension
      Iq1 = 1
      Iq2 = Basis%tab_basis(1)%nq
      Iq3 = Product(Basis%tab_basis(2:ndim)%nq) * Basis%tab_basis(ndim + 1)%nb

      Ib1 = 1
      Ib2 = Basis%tab_basis(1)%nb
      Ib3 = Product(Basis%tab_basis(2:ndim+1)%nb)

      allocate (BGG1(Ib1*Ib2*Iq3))
      BGG1 = CZERO
      GGG(1:Ib1, 1:Iq2, 1:Iq3) => G
      GGB(1:Ib1, 1:Ib2, 1:Iq3) => BGG1

      Call GridTOBasis_1D_cplx(GGB, GGG, Basis%tab_basis(1))

      ! Iterative transform for inner dimensions
      DO inb = 2, ndim - 1
         Iq1 = Product(Basis%tab_basis(1:inb - 1)%nq)
         Ib1 = Product(Basis%tab_basis(1:inb - 1)%nb)
         Iq2 = Basis%tab_basis(inb)%nq
         Ib2 = Basis%tab_basis(inb)%nb
         Ib3 = Product(Basis%tab_basis(inb + 1:ndim + 1)%nb)
         Iq3 = Product(Basis%tab_basis(inb + 1:ndim)%nq) * Basis%tab_basis(ndim + 1)%nb 

         allocate (BGG2(Ib1*Ib2*Iq3))
         BGG2 = CZERO
         GGG(1:Ib1, 1:Iq2, 1:Iq3) => BGG1
         GGB(1:Ib1, 1:Ib2, 1:Iq3) => BGG2

         Call GridTOBasis_1D_cplx(GGB, GGG, Basis%tab_basis(inb))
         BGG1 = BGG2
         deallocate (BGG2)
      END DO

      ! Final step for last dimension
      B(:) = CZERO

      Iq1 = Product(Basis%tab_basis(1:ndim - 1)%nq)
      Ib1 = Product(Basis%tab_basis(1:ndim - 1)%nb)
      Ib2 = Basis%tab_basis(ndim)%nb
      Iq2 = Basis%tab_basis(ndim)%nq
      Ib3 = Basis%tab_basis(ndim + 1)%nb
      Iq3 = Basis%tab_basis(ndim + 1)%nb

      GGG(1:Ib1, 1:Iq2, 1:Iq3) => BGG1
      GGB(1:Ib1, 1:Ib2, 1:Iq3) => B
      Call GridTOBasis_1D_cplx(GGB, GGG, Basis%tab_basis(ndim))
   END IF

   IF (debug) THEN
      write (out_unit, *) 'END GridTOBasis_nD_cplx'
      flush (out_unit)
   END IF
END SUBROUTINE
SUBROUTINE Calc_Q_grid(Q, Basis)
   IMPLICIT NONE

   ! Input
   TYPE(Basis_t), INTENT(IN)                 :: Basis

   ! Local variables
   INTEGER, ALLOCATABLE                      :: Tab_iq(:)
   INTEGER                                   :: inb, ndim, iq
   REAL(KIND=Rkind), INTENT(INOUT), ALLOCATABLE :: Q(:, :)
   LOGICAL                                   :: Endloop

   ! Determine the number of dimensions excluding the electronic basis
   ndim = SIZE(Basis%tab_basis) - 1

   ! Allocate output array Q and index array Tab_iq
   ALLOCATE(Q(Basis%nq, ndim), Tab_iq(ndim))

   ! Initialize the index array for multidimensional looping
   CALL Init_tab_ind(Tab_iq, Basis%NDindexq)

   iq = 0
   DO
       iq = iq + 1

       ! Update the multi-dimensional index
       CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop)
       IF (Endloop) EXIT

       ! Fill the Q array with the appropriate basis grid values
       DO inb = 1, ndim
           Q(iq, inb) = Basis%tab_basis(inb)%x(Tab_iq(inb))
       END DO

       ! Optional debug print
       ! PRINT*, iq, Q(iq, :)

   END DO

   ! Deallocate local array to avoid memory leak
   DEALLOCATE(Tab_iq)


END SUBROUTINE Calc_Q_grid


 SUBROUTINE Hermite_double_product_func(Hf, x, B1,B2,j, i)

      USE QDUtil_m
      real(kind=Rkind), intent(in)                  :: x,B1(:),B2(:)
      complex(kind=Rkind), intent(inout)            :: Hf !f(q)
      integer, intent(in)                           :: i, j
      real(kind=Rkind)                              :: qi, qj
      complex(kind=Rkind)                           ::  Hfi, Hfj

      ! Hfi = He(I,qi) = HermiteH[I,qi/Sqrt[2]] / Sqrt [ 2^I ]
      qi = (x - B1(1))*B1(3)
     ! Hfi =  sqrt(B1(3))*poly_Hermite_exp_cplx(qi,B1(2)/B1(3), i- 1)

      ! Hfj = He(J,qj) = HermiteH[J,qj/Sqrt[2]] / Sqrt [ 2^J ]
       qj = (x - B2(1))*B2(3)
       !Hfj= sqrt(B2(3))*poly_Hermite_exp_cplx(qj,B2(2)/B2(3), j- 1)
       

      !  Hf = Hfi *Hfj = He(J,x)*He(J,x)

      Hf = conjg(Hfj)*Hfi

   END SUBROUTINE Hermite_double_product_func


   SUBROUTINE Calc_Hermitec_int(f, Q, B1,B2,j, i)
      USE QDUtil_m
      real(kind=Rkind), intent(in)                  :: Q(:),B1(:),B2(:)
      complex(kind=Rkind), intent(inout)            :: f(:) 
      integer, intent(in)                           :: i,j
      integer                                       :: iq


      Do iq = 1,size(Q)
        ! call  Hermite_double_product_func(f(iq), Q(iq), B1,B2, j, i)
      End Do
   END SUBROUTINE Calc_Hermitec_int

   SUBROUTINE Hagedorn_ovelp_mat(Mat, nb, nq,B1,B2)
      USE QDUtil_m
      real(kind=Rkind), intent(in)        :: B1(:),B2(:)
      complex(kind=Rkind),allocatable, intent(inout)  :: Mat(:,:) 
      integer, intent(in)                 :: nq, nb
      real(kind=Rkind), allocatable       :: Q(:), w(:)
      complex(kind=Rkind),allocatable     :: Hf(:)
      real(kind=Rkind)                    :: SQeq,Qeq 
      integer                             :: iq,ib,jb

      if (allocated(Mat)) deallocate(Mat)
      allocate (Mat(nb,nb))
      allocate (Q(nq),w(nq),Hf(nq))
      call hercom(nq,Q,w)

    SQeq = sqrt(B1(3)*B1(3) + B2(3)*B2(3))/sqrt(TWO)
    Qeq = (B1(3)*B1(3)*B1(1) + B2(3)*B2(3)*B2(1))/(B1(3)*B1(3) + B2(3)*B2(3))
    w(:) = w(:)/SQeq
     Q(:) = Qeq+Q(:)/SQeq

    !print*,"Q",Q
      !print*,"w",w
      Do jb = 1,nb
         Do ib = 1,nb
            !call Calc_Hermitec_int(Hf, Q, B1,B2, jb, ib)
            Mat(jb,ib) = dot_product(Hf, w)
         End Do
      End Do

      deallocate (Hf,q,w)
      print*,"first basis parameters",B1
      print*,"second basis parameters",B2
      CALL  Write_VecMat(Mat, out_unit, 5,  info='Mat')

   END SUBROUTINE 
   
    SUBROUTINE d2psi1D_cplx(d2GGB,GGB,Basis)
    USE QDUtil_m
    complex(kind=Rkind), intent(in)                            :: GGB(:,:,:)
    complex(kind=Rkind), intent(inout)                         :: d2GGB(:,:,:)
    TYPE(Basis_t),intent(in)                                   :: Basis
    
    !locals variables---------------------------------------------------------
    
    
    logical, parameter                                         :: debug = .true.
  
    integer                                                    :: i1,i3
    
    
    !debuging----------------------------------------------------------------
    
    
    IF (debug) THEN
       flush (out_unit)
    END IF
    
      
      DO i3 = 1, ubound(GGB, dim=3)
         DO i1 = 1, ubound(GGB, dim=1)
           d2GGB(i1, :, i3) = d2GGB(i1, :, i3)+ matmul( Basis%d2gg(:, :, 1, 1) , GGB(i1,:,i3) )
         END DO
      END DO
    
    
    IF (debug) THEN
       flush (out_unit)
    END IF
    
   
 END SUBROUTINE

 SUBROUTINE pdv2psi_nD(d2psi_g, psi_g, Basis,ib)

    USE  QDUtil_m
    TYPE(Basis_t), intent(in), target               :: Basis
    complex(kind=Rkind), intent(in), target         :: psi_g(:)
    complex(kind=Rkind), intent(inout), target      :: d2psi_g(:)
    integer ,intent(in)                             :: ib
    complex(kind=Rkind), pointer                    :: psi_ggb(:, :, :)
    complex(kind=Rkind), pointer                    :: d2gg(:, :)
    complex(kind=Rkind), pointer                    :: d2psi_ggb(:, :, :)
    logical, parameter                              :: debug = .true.
    integer                                         :: iq, i1, i3, Ndim
    integer, allocatable                            :: Iq1(:), Iq2(:), Iq3(:), Ib1(:), Ib2(:), Ib3(:)

    IF (debug) THEN
   
       flush (out_unit)

    END IF

    Ndim = size(Basis%tab_basis)

    call Calc_index(Ib1, Ib2, Ib3, Iq1, Iq2, Iq3, Basis)
    d2psi_g(:) = CZERO


       d2psi_ggb(1:Iq1(ib), 1:Iq2(ib), 1:Iq3(ib)) => d2psi_g
       psi_ggb(1:Iq1(ib), 1:Iq2(ib), 1:Iq3(ib)) => psi_g
       d2gg(1:Iq2(ib), 1:Iq2(ib)) => Basis%tab_basis(ib)%d2gg

       DO i3 = 1, ubound(psi_ggb, dim=3)

          DO i1 = 1, ubound(psi_ggb, dim=1)

             d2psi_ggb(i1, :, i3) = d2psi_ggb(i1, :, i3) + matmul(d2gg, psi_ggb(i1, :, i3))

          END DO

       END DO

    Deallocate (Ib1, Ib2, Ib3, Iq1, Iq2, Iq3)

    IF (debug) THEN

       flush (out_unit)

    END IF

 END SUBROUTINE pdv2psi_nD

 SUBROUTINE psi_per_surf(Gsurf,G,Basis,ibe)

     IMPLICIT NONE
     TYPE(Basis_t), intent(in)                              :: Basis
     complex(kind=Rkind), intent(in) ,target                :: G(:)
     complex(kind=Rkind), intent(inout) ,allocatable        :: Gsurf(:)
     integer  , intent(in)                                  :: ibe
     
     complex(kind=Rkind), pointer                           :: psi_gb(:, :)
     integer                                                :: nq,nsurf

     IF(Allocated(Gsurf)) deallocate(Gsurf)

     nq = Basis%nq
     nsurf = Basis%tab_basis(size(Basis%tab_basis))%nb
     allocate(Gsurf(nq))
     Gsurf = CZERO

      psi_gb(1:nq, 1:nsurf) => G
    
      Gsurf(1:nq) =  psi_gb(1:nq, ibe) 


  END SUBROUTINE


  SUBROUTINE Calc_reduced_Density_surf(Rdensity,G,Basis)

     USE QDUtil_m

    TYPE(Basis_t), intent(in)                                    :: Basis
    complex(Kind=Rkind),intent(in)                               :: G(:)                                         
   TYPE(REDUCED_DENSIRY_t),intent(inout)                         :: Rdensity           

    !Locals variables ----------------------------------------------------------------------------

    integer, allocatable                                         :: Tab_iq(:)
    integer                                                      :: ndim, Iq,nq,inb,ib
    real(Kind=Rkind)                                             :: W
    logical                                                      :: Endloop
      ndim  = SIZE(Basis%tab_basis) - 1
      allocate (Tab_iq(ndim))
      call Init_tab_ind(Tab_iq, Basis%NDindexq) 
         Iq = 0
         DO
            Iq = Iq + 1

            CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop)
            IF (Endloop) exit
             W = ONE
            DO inb = 1, ndim  
               W = W*Basis%tab_basis(inb)%w(Tab_iq(inb))
            END DO
            If(ndim ==1) then
               Rdensity%prob(Iq) = Rdensity%prob(Iq) + conjg(G(Iq))*G(Iq)*W
            Else
               DO ib = 1, ndim
                   Rdensity%tab_prob(ib)%prob(Tab_iq(ib)) = Rdensity%tab_prob(ib)%prob(Tab_iq(ib)) + conjg(G(Iq))*G(Iq)*W
               END DO
            End If
         END DO            
      Deallocate(Tab_iq)
 END SUBROUTINE 

 SUBROUTINE Rdensity_alloc(Rdensity,Basis)
     
   USE  QDUtil_m  
   TYPE(REDUCED_DENSIRY_t),intent(inout)         :: Rdensity
   TYPE(Basis_t), intent(in), target             :: Basis
   integer                                       :: ib,ndim
   ndim  = SIZE(Basis%tab_basis) - 1
   If(allocated(Rdensity%tab_prob)) deallocate(Rdensity%tab_prob)
   allocate(Rdensity%tab_prob(ndim))
   IF(ndim==1) THEN
     allocate(Rdensity%prob(Basis%tab_basis(1)%nq))
     Rdensity%prob(:) = ZERO
    ! print*,'Rdensity%prob(:)',Rdensity%prob(:)
   ELSE
      DO ib=1,ndim       
        allocate(Rdensity%tab_prob(ib)%prob(Basis%tab_basis(ib)%nq))
        Rdensity%tab_prob(ib)%prob(:) =ZERO
        !print*,' ib,Rdensity%prob(:)',ib,Rdensity%tab_prob(ib)%prob(:) 
      END DO
   END IF        
 END SUBROUTINE



 SUBROUTINE Rdensity_Writing(Rdensity,Basis,nio,ib,t)
    
  USE  QDUtil_m  
  TYPE(REDUCED_DENSIRY_t),intent(inout)         :: Rdensity
  TYPE(Basis_t), intent(in), target             :: Basis
  integer  ,intent(in)                          :: nio
    integer  ,intent(in),optional               :: ib
  real(kind=Rkind),intent(in),optional          :: t 
  integer                                       :: ndim,nq,iq

  ndim  = SIZE(Basis%tab_basis) - 1
   If(present(ib)) then 
     nq= Basis%tab_basis(ib)%nq 
   Else
    nq= Basis%tab_basis(1)%nq 
   End If  

   DO Iq=1,nq
      If(ndim ==1) then
        If(present(t)) then
           write(nio,*) t, Basis%tab_basis(1)%x(Iq), Rdensity%prob(Iq)
        Else
           write(nio,*)  Basis%tab_basis(1)%x(Iq), Rdensity%prob(Iq)
        End If 
      Else 
        If(present(ib)) then        
          If(present(t)) then
            write(nio,*) t, Basis%tab_basis(ib)%X(Iq), Rdensity%tab_prob(ib)%prob(Iq)
          Else
            write(nio,*)  Basis%tab_basis(ib)%X(Iq), Rdensity%tab_prob(ib)%prob(Iq)
          End If 
          Else
          
           If(present(t)) then
             write(nio,*) t, Basis%tab_basis(1)%X(Iq),Basis%tab_basis(2)%X(Iq), &
             & Rdensity%tab_prob(1)%prob(Iq), Rdensity%tab_prob(2)%prob(Iq)
           Else
              write(nio,*) Basis%tab_basis(1)%X(Iq),Basis%tab_basis(2)%X(Iq), &
              & Rdensity%tab_prob(1)%prob(Iq), Rdensity%tab_prob(2)%prob(Iq)
           End If 
         End If
      End If    
   END DO
       
END SUBROUTINE

 SUBROUTINE Calc_Rdensity_Tot(Rdensity,Basis)
    
  USE  QDUtil_m  
  TYPE(REDUCED_DENSIRY_t),intent(inout)         :: Rdensity
  TYPE(Basis_t), intent(in), target             :: Basis

  integer                                       :: ndim,Ib

  ndim  = SIZE(Basis%tab_basis) - 1

  allocate(Rdensity%Norm(ndim))
      

       Rdensity%Norm(:) = ZERO
  
     DO Ib=1,ndim
       
           Rdensity%Norm(Ib) = sum(Rdensity%tab_prob(Ib)%prob(:))
  
     END DO

     Write(*,*) 'RedTot',Rdensity%Norm(:)


END SUBROUTINE


 SUBROUTINE Rdensity_Dealloc(Rdensity,Basis)
    
  USE  QDUtil_m  
  TYPE(REDUCED_DENSIRY_t),intent(inout)         :: Rdensity
  TYPE(Basis_t), intent(in), target             :: Basis
  integer                                       :: ib,ndim

  ndim  = SIZE(Basis%tab_basis) - 1

   IF (allocated(Rdensity%prob)) Deallocate(Rdensity%prob)

   IF (allocated(Rdensity%tab_prob)) then

   Do ib = 1,ndim
      IF (allocated(Rdensity%tab_prob(ib)%prob)) Deallocate(Rdensity%tab_prob(ib)%prob)
   End  Do

   End If

  IF (allocated(Rdensity%tab_prob)) Deallocate(Rdensity%tab_prob)
       
END SUBROUTINE

 SUBROUTINE Calc_reduced_density(Rdensity,B,Basis)
      USE  QDUtil_m
      TYPE(Basis_t), intent(in), target               :: Basis
      complex(kind=Rkind), intent(in), target         :: B(:)
      TYPE(REDUCED_DENSIRY_t),intent(inout)           :: Rdensity

        TYPE(REDUCED_DENSIRY_t)                       :: Rdensitytemp
        real(kind=Rkind)                              :: Norm 
        complex(kind=Rkind),allocatable               :: G(:),Gsurf(:)
        integer                                       :: nsurf,nq,ibe,iq,ndim,ib

         nsurf = Basis%tab_basis(size(Basis%tab_basis))%nb
         nq = nsurf*Basis%nq
         ndim = size(Basis%tab_basis)-1

         call Rdensity_Dealloc(Rdensity,Basis)
         call Rdensity_alloc(Rdensity,Basis)          
         allocate(G(nq))
         CALL BasisTOGrid_nD_cplx(G,B,Basis)
        
       DO ibe = 1,nsurf

          CALL psi_per_surf(Gsurf,G,Basis,ibe)
          call  Rdensity_alloc(Rdensitytemp,Basis)
          CALL Calc_reduced_Density_surf(Rdensitytemp,Gsurf,Basis)
           If(ndim==1) then
             Rdensity%prob(:)  = Rdensity%prob(:)+Rdensitytemp%prob(:)
          Else
           Do ib =1, ndim
              Rdensity%tab_prob(ib)%prob(:)  = Rdensity%tab_prob(ib)%prob(:)+Rdensitytemp%tab_prob(ib)%prob(:)
           END DO
          End If
          call   Rdensity_Dealloc(Rdensitytemp,Basis)
          

       END DO   
      deallocate(G)

 END SUBROUTINE


 SUBROUTINE Calc_sum_psi(Int,G1,G2,Basis)
   USE QDUtil_m
   logical, parameter                            :: debug = .false.
   TYPE(Basis_t),intent(in)                      :: Basis
   complex(kind=Rkind), intent(in)   ,target     :: G1(:),G2(:)
   complex(kind=Rkind) ,intent(inout)            :: Int
   !Locals variabls ----------------------------------------------------------
   complex(kind=Rkind), pointer                  :: psi_gb1(:, :),psi_gb2(:, :)
   logical                                       :: Endloop_q
   complex(kind=Rkind), allocatable              :: Intel(:)
   real(kind=Rkind), allocatable                 :: Nel(:)
   real(kind=Rkind)                              :: W
   integer, allocatable                          :: Tab_iq(:)
   integer                                       :: iq, inbe,inb,nq,nsurf,ndim
   IF (debug) THEN
      write (out_unit, *) 'Beging Evaluating integral'
      flush (out_unit)
   END IF
   nq = Basis%nq
   nsurf = Basis%tab_basis(size(Basis%tab_basis))%nb
   ndim = size(Basis%tab_basis) - 1
   allocate (Nel(nsurf))
   allocate (Intel(nsurf)) 
   allocate (Tab_iq(ndim))
   psi_gb1(1:nq,1:nsurf) => G1
   psi_gb2(1:nq,1:nsurf) => G2
 Nel(:) = ZERO
 Intel(:)=ZERO
 Int = ZERO
 DO inbe = 1, Basis%tab_basis(size(Basis%tab_basis))%nb !electronic state
    Intel(inbe) = CZERO
    Call Init_tab_ind(Tab_iq, Basis%NDindexq)
    iq = 0
    DO
       iq = iq + 1
       CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop_q)
       IF (Endloop_q) exit
      
        W = ONE
        DO inb = 1, ndim
           W = W*Basis%tab_basis(inb)%w(tab_iq(inb))
        END DO
       Nel(inbe) = Nel(inbe) + conjg(psi_gb1(iq, inbe))*psi_gb1(iq, inbe)*W
       Intel(inbe) = Intel(inbe) + psi_gb1(iq, inbe)*psi_gb2(iq, inbe)*W
    END DO
 END DO
 Int = sum(Intel)/(Sum(Nel)**2)
 Deallocate (Tab_iq)
 Deallocate(Intel,Nel)
 IF (debug) THEN
    write (out_unit, *) 'END Evaluating integral'
    flush (out_unit)
 END IF
END SUBROUTINE

SUBROUTINE Calc_d0d1d2W_temp(Q, Qt, SQt, Bt, Pt, d0W, d1W, d2W, nderiv)
   IMPLICIT NONE
   REAL(kind=Rkind),    INTENT(IN)     :: Q       ! Position
   REAL(kind=Rkind),    INTENT(IN)     :: Qt      ! Center of Gaussian
   REAL(kind=Rkind),    INTENT(IN)     :: SQt     ! Width (std dev)
   REAL(kind=Rkind),    INTENT(IN)     :: Bt      ! Related to complex width
   REAL(kind=Rkind),    INTENT(IN)     :: Pt      ! Momentum
   LOGICAL,             INTENT(IN)     :: nderiv  ! Flag to compute derivatives
 
   COMPLEX(kind=Rkind), INTENT(INOUT)  :: d0W     ! complex exponential part of the gaussian value
   COMPLEX(kind=Rkind), INTENT(INOUT)  :: d1W     ! First derivative of the complex exponential part of the gaussian value
   COMPLEX(kind=Rkind), INTENT(INOUT)  :: d2W     ! Second derivative of the complex exponential part of the gaussian value
 
   ! === Local variables ===
   COMPLEX(kind=Rkind) :: cst, cst1, cst2
   REAL(kind=Rkind)    :: DQ
 
   ! === Compute Q - Qt
   DQ = Q - Qt
 
   ! === Compute exponential part of the wavefunction
   cst = -HALF * EYE * (Bt / (SQt * SQt)) * Q * Q + EYE * (Pt / SQt) * Q
 
   ! === Derivatives of the exponential (symbolic form)
   cst1 = -EYE * (Bt / (SQt * SQt)) * Q + EYE * (Pt / SQt)
   cst2 = -EYE * (Bt / (SQt * SQt)) + cst1 * cst1
 
   ! === Set output values
   IF (nderiv) THEN
     d0W = EXP(cst)
     d1W = cst1 * d0W
     d2W = cst2 * d0W
   ELSE
     d0W = EXP(cst)
     d1W = CZERO
     d2W = CZERO
   END IF
 
   ! Uncomment for debugging:
   ! PRINT *, 'd0W = ', d0W, ' d1W = ', d1W, ' d2W = ', d2W
 
 END SUBROUTINE Calc_d0d1d2W_temp


 SUBROUTINE Calc_d0d1d2W(Q, Qt, SQt, Bt, Pt, d0W, d1W, d2W, nderiv)
   IMPLICIT NONE
 
   ! === Input ===
   REAL(kind=Rkind),    INTENT(IN)     :: Q        ! Position
   REAL(kind=Rkind),    INTENT(IN)     :: Qt       ! Center of the wave packet
   REAL(kind=Rkind),    INTENT(IN)     :: SQt      ! Width of the Gaussian (not used here)
   REAL(kind=Rkind),    INTENT(IN)     :: Bt       ! Imaginary width parameter
   REAL(kind=Rkind),    INTENT(IN)     :: Pt       ! Momentum
   LOGICAL,             INTENT(IN)     :: nderiv   ! Flag to compute derivatives
 
   ! === Output ===
   COMPLEX(kind=Rkind), INTENT(INOUT)  :: d0W      ! complex exponential part of the gaussian value
   COMPLEX(kind=Rkind), INTENT(INOUT)  :: d1W      ! First derivative of the complex exponential part of the gaussian value
   COMPLEX(kind=Rkind), INTENT(INOUT)  :: d2W      ! Second derivative of the complex exponential part of the gaussian value
 
   ! === Local variables ===
   COMPLEX(kind=Rkind) :: cst, cst1, cst2
   REAL(kind=Rkind)    :: DQ  ! Displacement Q - Qt
 
   ! === Compute displacement from the center
   DQ = Q - Qt
 
   ! === Compute wavefunction and its derivatives
   IF (nderiv) THEN
      ! Compute complex phase for wavefunction
      cst  = -EYE * HALF * Bt * DQ * DQ + EYE * Pt * DQ
      ! First derivative of the exponential argument
      cst1 = -EYE * Bt * DQ + EYE * Pt
      ! Second derivative of the exponential argument
      cst2 = -EYE * Bt + cst1 * cst1
 
      ! Compute wavefunction and derivatives
      d0W = EXP(cst)
      d1W = cst1 * d0W
      d2W = cst2 * d0W
 
   ELSE
      ! Only compute the exponential part, derivatives are set to zero
      cst  = -EYE * HALF * Bt * DQ * DQ + EYE * Pt * DQ
      d0W  = EXP(cst)
      d1W  = CZERO
      d2W  = CZERO
   END IF
 
   ! Uncomment this line for debugging:
   ! PRINT *, 'd0W = ', d0W, ' d1W = ', d1W, ' d2W = ', d2W
 
 END SUBROUTINE Calc_d0d1d2W

 SUBROUTINE Change_Basis_Parameters(Basis, Qt, SQt, At, Pt)
   IMPLICIT NONE
 
   ! === Input/Output ===
   TYPE(Basis_t), INTENT(INOUT)          :: Basis   ! Basis object to modify
   REAL(kind=Rkind), INTENT(IN)          :: Qt(:)    ! New centers (position)
   REAL(kind=Rkind), INTENT(IN)          :: SQt(:)   ! New scaling factors
   REAL(kind=Rkind), INTENT(IN)          :: Pt(:)    ! New impulses (momentum)
   COMPLEX(kind=Rkind), INTENT(IN)       :: At(:)    ! New complex width parameters
 
   ! === Local variables ===
   INTEGER                               :: ndim     ! Number of degrees of freedom
   INTEGER                               :: ib       ! Basis index
 
   ! === Compute number of spatial dimensions (assume last index is for electronic state)
   ndim = SIZE(Basis%tab_basis) - 1
 
   ! === Update basis parameters for each spatial dimension
   DO ib = 1, ndim
     Basis%tab_basis(ib)%Q0      = Qt(ib)
     Basis%tab_basis(ib)%SCALEQ  = SQt(ib)
     Basis%tab_basis(ib)%Imp_k   = Pt(ib)
     Basis%tab_basis(ib)%alpha   = At(ib)
   END DO
 
 END SUBROUTINE Change_Basis_Parameters

 SUBROUTINE Get_Basis_Parameters(Basis, Qt, SQt, At, Pt)
   IMPLICIT NONE
 
   ! === Input ===
   TYPE(Basis_t), INTENT(IN)           :: Basis      ! basis structure
 
   ! === Output ===
   REAL(kind=Rkind), INTENT(OUT)       :: Qt(:)      ! position centers
   REAL(kind=Rkind), INTENT(OUT)       :: SQt(:)     ! scaling (widths)
   REAL(kind=Rkind), INTENT(OUT)       :: Pt(:)      ! momentum components
   COMPLEX(kind=Rkind), INTENT(OUT)    :: At(:)      ! complex parameters alpha
 
   ! === Local variables ===
   INTEGER                             :: ndim       ! Number of spatial dimensions
   INTEGER                             :: ib         ! Loop index
 
   ! Determine the number of spatial dimensions (excluding electronic state)
   ndim = SIZE(Basis%tab_basis) - 1
 
   ! Extract the parameters from the basis
   DO ib = 1, ndim
     Qt(ib)  = Basis%tab_basis(ib)%Q0
     SQt(ib) = Basis%tab_basis(ib)%SCALEQ
     Pt(ib)  = Basis%tab_basis(ib)%Imp_k
     At(ib)  = Basis%tab_basis(ib)%alpha
   END DO
 
 END SUBROUTINE Get_Basis_Parameters

END MODULE Basis_m
