MODULE Basis_m
   USE QDUtil_m
   USE NDindex_m
   USE polyortho_m

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: Basis_t, Read_Basis, Write_Basis, Basis_IS_allocated, Deallocate_Basis, Basis_IS_allocatedtot
   PUBLIC :: Calc_dngg_grid, Calc_tranpose_d0gb, test_basitogridgridtobasis
   PUBLIC :: GridTOBasis_nD_cplx, BasisTOGrid_nD_cplx
   PUBLIC :: Calc_Q_grid, Calc_index
   PUBLIC :: Scale_Basis, init_Basis1_TO_Basis2
   PUBLIC :: construct_primitive_basis, construct_primitive_basis0, construct_primitive_basis1
   PUBLIC :: Hermite_double_product_func, Hagedorn_ovelp_mat
   PUBLIC :: Calc_dp0G,Calc_dq0G,Calc_daG

   TYPE :: Basis_t
      integer                             :: nb_basis = 0
      integer                             :: nb = 0
      integer                             :: nq = 0
      real(kind=Rkind)                    :: Q0 = 0_Rkind
      real(kind=Rkind)                    :: scaleQ = 0_Rkind
      real(kind=Rkind)                    :: A = 0_Rkind
      real(kind=Rkind)                    :: B = 0_Rkind
      real(kind=Rkind)                    :: imp_k = 0_Rkind 
      character(len=:),    allocatable    :: Basis_name
      real(kind=Rkind),    allocatable    :: x(:)
      real(kind=Rkind),    allocatable    :: w(:)
      complex(kind=Rkind),    allocatable    :: d0gb(:, :)         ! basis functions d0gb(nq,nb)
      complex(kind=Rkind),    allocatable    :: d0bgw(:, :)        ! transpose of basis functions d0gb(nb,nq)
      complex(kind=Rkind),    allocatable    :: d1gb(:, :, :)      ! basis functions d2gb(nq,nb,1)
      complex(kind=Rkind),    allocatable    :: d1gg(:, :, :)      ! basis functions d2gg(nq,nq,1)
      complex(kind=Rkind),    allocatable    :: d2gb(:, :, :, :)   ! basis functions d2gb(nq,nb,1,1)
      complex(kind=Rkind),    allocatable    :: d2gg(:, :, :, :)   ! basis functions d2gg(nq,nq,1,1)
      complex(kind=Rkind),    allocatable :: d0gb_cplx(:, :)         ! basis functions d0gb(nq,nb)
      complex(kind=Rkind), allocatable    :: dagb(:, :)         !  basis for dérivative respect to alpha   d0gb(nq,nb)
      complex(kind=Rkind), allocatable    :: dagg(:, :)         ! dérivative respect to alpha on the grid  dagg(nq,nq)     
      complex(kind=Rkind), allocatable    :: dq0gb(:, :)        !  basis for dérivative respect to q0  dq0gb(nq,nb)
      complex(kind=Rkind), allocatable    :: dq0gg(:, :)        ! dérivative respect to q0 on the grid dq0gq(nq,nq)
      complex(kind=Rkind), allocatable    :: dp0gb(:, :)        !  basis for dérivative respect to p0  dq0gb(nq,nb)
      complex(kind=Rkind), allocatable    :: dp0gg(:, :)        ! dérivative respect to p0 on the grid dq0gq(nq,nq)    
      complex(kind=Rkind), allocatable    :: S(:, :)            ! for Hagedorn transformation
      TYPE(NDindex_t)                     :: NDindexq
      TYPE(NDindex_t)                     :: NDindexb
      TYPE(Basis_t), allocatable          :: tab_basis(:)       !  for more than one Basis.

   END TYPE Basis_t

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
      !write(out_unit,*) "Basis_name=", Basis%Basis_name
      !write(out_unit,*) "n_basis=", Basis%nb_basis
      !write(out_unit,*) 'nb,nq',Basis%nb,Basis%nq
      !write(out_unit,*) "Q0,scaleQ ", Basis%Q0,Basis%scaleQ
      !write(out_unit,*)  'A,B',Basis%A,Basis%B
      !write(out_unit,*) 'nb_basis',Basis%nb_basis

      IF (.NOT. allocated(Basis%x)) THEN
         write(out_unit,*)' Basis table x is not allocated.'
      ELSE
         CALL  Write_VecMat(Basis%x, out_unit, 5,  info='x')
      END IF
      write (out_unit, *)
      IF (.NOT. allocated(Basis%W)) THEN
           write(out_unit,*)' Basis table w is not allocated.'
      ELSE
         CALL  Write_VecMat(Basis%w, out_unit, 5,  info='w')
      END IF
      write (out_unit, *)
       IF (.NOT.allocated(Basis%d0gb)) THEN
        write(out_unit,*)' Basis table d0gb is not allocated.'
      ELSE
      CALL   Write_VecMat(Basis%d0gb, out_unit, 5,  info='d0gb')
       END IF
          

       write (out_unit, *)
       IF (.NOT.allocated(Basis%dp0gb)) THEN
        write(out_unit,*)' Basis table dp0gb is not allocated.'
      ELSE
      CALL   Write_VecMat(Basis%dp0gb, out_unit, 5,  info='dp0gb')
       END IF
      
       
       write (out_unit, *)
       IF (.NOT.allocated(Basis%dq0gb)) THEN
        write(out_unit,*)' Basis table dq0gb is not allocated.'
      ELSE
      CALL   Write_VecMat(Basis%dq0gb, out_unit, 5,  info='dq0gb')
       END IF
       
       write (out_unit, *)
       IF (.NOT.allocated(Basis%dagb)) THEN
        write(out_unit,*)' Basis table dagb is not allocated.'
      ELSE
      CALL   Write_VecMat(Basis%dagb, out_unit, 5,  info='dagb')
       END IF

       write (out_unit, *)
       IF (.NOT.allocated(Basis%dagg)) THEN
        write(out_unit,*)' Basis table dagg is not allocated.'
      ELSE
      CALL   Write_VecMat(Basis%dagg, out_unit, 5,  info='dagg')
       END IF

       write (out_unit, *)
       IF (.NOT.allocated(Basis%dp0gg)) THEN
        write(out_unit,*)' Basis table dp0gg is not allocated.'
      ELSE
      CALL   Write_VecMat(Basis%dp0gg, out_unit, 5,  info='dp0gg')
       END IF
      
       write (out_unit, *)
       IF (.NOT.allocated(Basis%dq0gg)) THEN
        write(out_unit,*)' Basis table dq0gg is not allocated.'
      ELSE
      CALL   Write_VecMat(Basis%dq0gg, out_unit, 5,  info='dq0gg')
       END IF

      write(out_unit,*)
        IF (.NOT.allocated(Basis%d0bgw)) THEN
           write(out_unit,*)' Basis table d0bgw is not allocated.'
         ELSE
      CALL   Write_VecMat(Basis%d0bgw,out_unit,5, info='d0gbw')
        END IF

         write(out_unit,*)
         IF (.NOT.allocated(Basis%d1gb)) THEN
           write(out_unit,*)' Basis table d1gb is not allocated.'
         ELSE
      CALL   Write_VecMat(Basis%d1gb(:,:,1),out_unit,5, info='d1gb')
         END IF
         write(out_unit,*)
         IF (.NOT.allocated(Basis%d1gg)) THEN
           write(out_unit,*)' Basis table d1gb is not allocated.'
         ELSE
      CALL   Write_VecMat(Basis%d1gg(:,:,1),out_unit,5, info='d1gg')
         END IF

         write(out_unit,*)
         IF (.NOT.allocated(Basis%d2gb)) THEN
           write(out_unit,*)' Basis table d1gb is not allocated.'
         ELSE
      CALL   Write_VecMat(Basis%d2gb(:,:,1,1),out_unit,5, info='d2gb')
         END IF

         write(out_unit,*)
         IF (.NOT.allocated(Basis%d2gg)) THEN
           write(out_unit,*)' Basis table d2gg is not allocated.'
         ELSE
      CALL   Write_VecMat(Basis%d2gg(:,:,1,1),out_unit,5, info='d2gg')
         END IF
          

         write(out_unit,*)
         IF (.NOT.allocated(Basis%S)) THEN
            write(out_unit,*)' Basis table S is not allocated.'
         ELSE
      CALL   Write_VecMat(Basis%S,out_unit,5, info='S')
         END IF  
         IF (allocated(Basis%tab_basis)) THEN
            DO i = 1, size(Basis%tab_basis)
            IF (Basis%tab_basis(i)%Basis_name /= 'el') CALL Write_Basis(Basis%tab_basis(i))
           END DO
        END IF
      write (out_unit, *) '--------------------------------------------------------------------------'

   END SUBROUTINE Write_Basis

   RECURSIVE SUBROUTINE init_Basis1_TO_Basis2(Basis2, Basis1)
      USE QDUtil_m
      TYPE(Basis_t), intent(in)           :: Basis1
      TYPE(Basis_t), intent(inout)        :: Basis2
      integer                             :: ib

      IF (allocated(Basis1%tab_basis)) THEN
         call Deallocate_Basis(Basis2)
         Basis2%Basis_name = Basis1%Basis_name
         Basis2%nb_basis = Basis1%nb_basis
         allocate (Basis2%tab_basis(Basis2%nb_basis))
         DO ib = 1, Basis1%nb_basis
            CALL init_Basis1_TO_Basis2(Basis2%tab_basis(ib), Basis1%tab_basis(ib))
         END DO
         Basis2%nb = 1
         Basis2%nq = 1
         DO ib = 1, Basis1%nb_basis
            if (Basis2%tab_basis(ib)%Basis_name == 'el') cycle
            Basis2%nb = Basis2%nb*Basis2%tab_basis(ib)%nb
            Basis2%nq = Basis2%nq*Basis2%tab_basis(ib)%nq
         END DO
      ELSE
         Basis2%Basis_name       = Basis1%Basis_name
         Basis2%nb_basis         = Basis1%nb_basis
         Basis2%nb               = Basis1%nb
         Basis2%nq               = Basis1%nq
         Basis2%Q0               = Basis1%Q0
         Basis2%SCALEQ           = Basis1%SCALEQ
         Basis2%Imp_k            = Basis1%Imp_k
         Basis2%A                = Basis1%A
         Basis2%B                = Basis1%B
         Basis2%S                = Basis1%S
      END IF

   END SUBROUTINE init_Basis1_TO_Basis2

   RECURSIVE SUBROUTINE Deallocate_Basis(Basis)
      USE QDUtil_m

      TYPE(Basis_t), intent(inout)     :: Basis
      integer                          :: i

      ! write(out_unit,*) '********************************************************************'
      !write(out_unit,*) 'Deallocate_Basis'
      !write(out_unit,*) "Basis_name=", Basis%Basis_name
      !write(out_unit,*) "n_basis=", Basis%nb_basis
      !write(out_unit,*) 'nb,nq',Basis%nb,Basis%nq
      !write(out_unit,*) "Q0=", Basis%Q0
      !write(out_unit,*) "scaleQ =", Basis%scaleQ

      IF (allocated(Basis%tab_basis)) THEN
         deallocate (Basis%NDindexq%Tab0)
         deallocate (Basis%NDindexb%Tab0)
         deallocate (Basis%tab_basis)
      END IF
      write (out_unit, *)

      IF (allocated(Basis%x)) THEN
         deallocate (Basis%x)
         ! write(out_unit,*)' Basis table x is now deallocated.'
      END IF
      write (out_unit, *)
      IF (allocated(Basis%W)) THEN
         deallocate (Basis%W)
         !write(out_unit,*)' Basis table w is now deallocated.'
      END IF
      write (out_unit, *)
      IF (allocated(Basis%d0gb)) THEN
         deallocate (Basis%d0gb)
         ! write(out_unit,*)' Basis table d0gb is now deallocated.'
      END IF
      write (out_unit, *)
      IF (allocated(Basis%d0bgw)) THEN
         deallocate (Basis%d0bgw)
         ! write(out_unit,*)' Basis table d0bgw is now deallocated.'
      END IF
      write (out_unit, *)
      IF (allocated(Basis%d1gb)) THEN
         deallocate (Basis%d1gb)
         ! write(out_unit,*)' Basis table d1gb is now deallocated.'
      END IF
      write (out_unit, *)
      IF (allocated(Basis%d1gg)) THEN
         deallocate (Basis%d1gg)
         ! write(out_unit,*)' Basis table d1gb is now deallocated.'
      END IF
      write (out_unit, *)
      IF (allocated(Basis%d2gb)) THEN
         deallocate (Basis%d2gb)
         !write(out_unit,*)' Basis table d1gb is now deallocated.'
      END IF
      write (out_unit, *)
      IF (allocated(Basis%d2gg)) THEN
         deallocate (Basis%d2gg)
         !  write(out_unit,*)' Basis table d2gg is now deallocated.'
      END IF
      IF (allocated(Basis%tab_basis)) THEN
         DO i = 1, size(Basis%tab_basis)
            CALL Deallocate_Basis(Basis%tab_basis(i))
         END DO
      END IF
      ! write(out_unit,*) '********************************************************************'

   END SUBROUTINE Deallocate_Basis

   RECURSIVE SUBROUTINE Read_Basis(Basis, nio)
      USE QDUtil_m
      logical, parameter                       :: debug = .true.
      !logical,             parameter          ::debug = .false.
      TYPE(Basis_t), intent(inout)             :: Basis
      integer, intent(in)                      :: nio
      integer                                  :: err_io, nb, nq, i, j, nb_basis, ib
      character(len=Name_len)                  :: name
      real(kind=Rkind)                            :: A, B, scaleQ, Q0, d0, d2, X1, W1,Imp_k

      NAMELIST /basis_nD/ name, nb_basis, nb, nq, A, B, scaleQ, Q0,Imp_k
      nb_basis = 0
      nb = 0
      nq = 0
      A = ZERO
      B = ZERO
      Q0 = ZERO
      scaleQ = ONE
      name = '0'
      Imp_k= ZERO

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
         CALL string_uppercase_TO_lowercase(Basis%Basis_name)
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
         Basis%Basis_name = trim(adjustl(name))
         allocate (Basis%S(nb, nb))
         Basis%S(:, :) = ZERO
         do ib = 1, Basis%nb
            Basis%S(ib, ib) = ONE
         end do

         CALL string_uppercase_TO_lowercase(Basis%Basis_name)
      END IF
   END SUBROUTINE Read_Basis

   RECURSIVE SUBROUTINE construct_primitive_basis0(Basis)
      USE QDUtil_m
      logical, parameter                     :: debug = .true.
      !logical,             parameter        ::debug = .false.
      TYPE(Basis_t), intent(inout)           :: Basis
      integer, allocatable                   :: NDend_q(:)
      integer, allocatable                   :: NDend_b(:)
      integer                                :: nb, nq, i, j
      character(len=Name_len)                :: name
      ! write(out_unit,*) ' Begin  construct primitive  Basis '
      IF (allocated(Basis%tab_basis)) THEN
         allocate (NDend_q(Basis%nb_basis - 1))
         allocate (NDend_b(Basis%nb_basis - 1))
         DO i = 1, Basis%nb_basis - 1
            NDend_q(i) = Basis%tab_basis(i)%nq
            NDend_b(i) = Basis%tab_basis(i)%nb
         END DO
         CALL Init_NDindex(Basis%NDindexq, NDend_q, Basis%nb_basis - 1)
         CALL Init_NDindex(Basis%NDindexb, NDend_b, Basis%nb_basis - 1)

         DO i = 1, Basis%nb_basis
            CALL construct_primitive_basis(Basis%tab_basis(i))
         END DO

      ELSE
         SELECT CASE (Basis%Basis_name)
         CASE ('el')
            write (6, *) 'Electronic basis. Electronic state number:', basis%nb
            basis%nq = 0
         CASE ('boxab')
            CALL Construct_Basis_Sin(Basis)
            Basis%Q0 = Basis%A
            Basis%scaleQ = pi/(Basis%B - Basis%A)
         CASE ('fourier')
            CALL Construct_Basis_Fourier(Basis)
         CASE ('herm', 'ho')
            CALL Construct_Basis_Ho(Basis)
         CASE default
            STOP 'ERROR  Noting to construct'
         END SELECT
         !  this part wil not have sens for 'el' basis
         CALL Scale_Basis(Basis, Basis%Q0, Basis%scaleQ)
         CALL Calc_tranpose_d0gb(Basis)
         CALL Calc_dngg_grid(Basis)
         CALL CheckOrtho_Basis(Basis, nderiv=2)
      END IF
      ! write(out_unit,*) ' End  construct  primitive Basis '
   END SUBROUTINE construct_primitive_basis0

   SUBROUTINE construct_primitive_basis1(Basis, x, p,sx)
      USE QDUtil_m
      logical, parameter                     :: debug = .true.
      real(kind=Rkind), intent(in)           :: x(:), sx(:),p(:)
      real(kind=Rkind)                       :: x0, sx0,p0
      !logical,             parameter        ::debug = .false.
      TYPE(Basis_t), intent(inout)           :: Basis
      integer, allocatable                   :: NDend_q(:)
      integer, allocatable                   :: NDend_b(:)
      integer                                :: nb, nq, i, j
      character(len=Name_len)                :: name

       !write(out_unit,*) ' Begin  construct primitive  Basis '
       !write(out_unit,*) 'ici',x,sx,p

      IF (allocated(Basis%tab_basis)) THEN
         !write(out_unit,*) 'ici',x,sx,p,size(Basis%tab_basis)
         DO i = 1, size(Basis%tab_basis)
            if (Basis%tab_basis(i)%Basis_name == 'herm' .or. Basis%tab_basis(i)%Basis_name == 'ho') then
               call Construct_Basis_Hagedorn(Basis%tab_basis(i), x(i),p(i),sx(i))
               call Scale_Basis(Basis%tab_basis(i), Basis%tab_basis(i)%Q0, Basis%tab_basis(i)%scaleQ)
               call Calc_tranpose_d0gb(Basis%tab_basis(i))
               call Calc_dngg_grid(Basis%tab_basis(i))
               call CheckOrtho_Basis(Basis%tab_basis(i), nderiv=2)
            end if
         END DO
      ELSE
         if (Basis%Basis_name == 'herm' .or. Basis%Basis_name == 'ho') then
            call Construct_Basis_Hagedorn(Basis, x(1),p(1),sx(1))
            call Scale_Basis(Basis, Basis%Q0, Basis%scaleQ)
            call Calc_tranpose_d0gb(Basis)
            call Calc_dngg_grid(Basis)
            call CheckOrtho_Basis(Basis, nderiv=2)
         end if  
      END IF
       write(out_unit,*) ' End  construct  primitive Basis '
   END SUBROUTINE construct_primitive_basis1

   RECURSIVE SUBROUTINE construct_primitive_basis(Basis, x, p,sx)
      USE QDUtil_m
      logical, parameter                   :: debug = .true.
      !logical,             parameter      ::debug = .false.
      TYPE(Basis_t), intent(inout)         :: Basis
      real(kind=Rkind), intent(in), optional  :: x(:), sx(:) ,p(:)

      if (present(x) .and. present(sx) .and. present(p)) then
         !write(out_unit,*) ' S will be constructed for Ho Basis'
         call construct_primitive_basis1(Basis, x,p, sx)
      else
         call construct_primitive_basis0(Basis)
      end if

   END SUBROUTINE construct_primitive_basis

    SUBROUTINE Construct_Basis_Sin(Basis) ! sin : boxAB with A=0 and B=pi
      USE QDUtil_m

      TYPE(Basis_t), intent(inout)        :: Basis
      real(kind=Rkind)                       :: dx
      integer                             :: ib, iq, nb, nq

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

   SUBROUTINE Construct_Basis_Fourier(Basis) !basis_name fourier[-pi,pi]
      USE QDUtil_m

      TYPE(Basis_t), intent(inout)        :: Basis
      real(kind=Rkind)                    :: dx
      integer                             :: ib, iq, nb, nq

      nb = Basis%nb
      nq = Basis%nq
      dx = TWO*pi/nq

      !>*************** grid and weight *************************************8
      Basis%x = [(iq*dx - dx/2 - pi, iq=1, nq)]
      Basis%w = [(dx, iq=1, nq)]
      !> **********************allocation **********************************
      allocate (Basis%d0gb(nq, nb))
      allocate (Basis%d1gb(nq, nb, 1))
      allocate (Basis%d2gb(nq, nb, 1, 1))

      DO ib = 1, nb
         DO iq = 1, nq
            call d0d1d2fourier(Basis%x(iq),Basis%d0gb(iq, ib),Basis%d1gb(iq, ib, 1)&
                           &,Basis%d2gb(iq, ib, 1, 1),ib)
         END DO                  
      END DO
      !**************************************************************************
      IF (Basis%nb == Basis%nq .AND. mod(Basis%nb, 2) == 0) THEN
         Basis%d0gb(:, nb) = Basis%d0gb(:, nb)/sqrt(TWO)
         Basis%d1gb(:, nb, :) = Basis%d1gb(:, nb, :)/sqrt(TWO)
         Basis%d2gb(:, nb, :, :) = Basis%d2gb(:, nb, :, :)/sqrt(TWO)
      END IF


   END SUBROUTINE Construct_Basis_Fourier

   SUBROUTINE Construct_Basis_el(Basis) ! 'el' :
      USE QDUtil_m

      TYPE(Basis_t), intent(inout)  :: Basis
      Basis%nq = 0
      RETURN

   END SUBROUTINE Construct_Basis_el

  SUBROUTINE Construct_Basis_Ho(Basis) ! HO :
      USE QDUtil_m

      TYPE(Basis_t), intent(inout)               :: Basis
      integer                                    :: iq, ib
      real(kind=Rkind)                           :: alpha,SQeq,Peq
      integer                                    :: nb, nq
     
      nb    = Basis%nb
      nq    = Basis%nq

      SQeq  = Basis%scaleQ
      Peq    = Basis%Imp_k/SQeq

      if (allocated(Basis%x)) deallocate (Basis%x)
      if (allocated(Basis%w)) deallocate (Basis%w)
      if (allocated(Basis%d0gb)) deallocate (Basis%d0gb)
      if (allocated(Basis%d1gb)) deallocate (Basis%d1gb)
      if (allocated(Basis%d2gb)) deallocate (Basis%d2gb)

      if (allocated(Basis%dagb)) deallocate (Basis%dagb)
      if (allocated(Basis%dp0gb)) deallocate (Basis%dp0gb)
      if (allocated(Basis%dq0gb)) deallocate (Basis%dq0gb)


      allocate (Basis%x(nq))
      allocate (Basis%w(nq))
      call hercom(nq, Basis%x(:), Basis%w(:))

      allocate (Basis%d0gb(nq, nb))
      allocate (Basis%d1gb(nq, nb, 1))
      allocate (Basis%d2gb(nq, nb, 1, 1))

      DO iq = 1, nq
         DO ib = 1, nb
           call  d0d1d2poly_Hermite_exp_cplx(Basis%x(iq),Peq, ib - 1, Basis%d0gb(iq, ib),&
             Basis%d1gb(iq, ib, 1),Basis%d2gb(iq, ib, 1, 1), .TRUE.)

         END DO
      END DO
   END SUBROUTINE Construct_Basis_Ho


   SUBROUTINE Construct_Basis_Hagedorn(Basis,x,p,sx) ! Hagedorn :
      USE QDUtil_m
      TYPE(Basis_t), intent(inout)               :: Basis
      real(kind=Rkind),intent(in)                :: x,sx,p
      !----------------local variables ----------------------------------------------
      complex(kind=Rkind),allocatable            :: d0gb(:,:),d0bgw(:,:)
      real(kind=Rkind),allocatable               :: w(:),Q(:)
      real(kind=Rkind)                           :: Qeq, SQeq
      integer                                    :: nb, nq,ib,iq,jb
      real(kind=Rkind)                           :: Q0,s0,p0,B1(3),B2(3)
      
      nb    = Basis%nb
      nq    = Basis%nq

      Q0    = Basis%Q0
      s0    =  Basis%SCALEQ
      p0    = Basis%Imp_k
      B1(1)=Q0;B1(2)=p0;B1(3)=s0
      B2(1)=x;B2(2)=p;B2(3)=sx
    

      Basis%Q0    = x
      Basis%Imp_k = p/sx
      Basis%SCALEQ = sx
     
      
      SQeq = sqrt(s0*s0 + sx*sx)/sqrt(TWO)
      Qeq = (s0*s0*Q0 + sx*sx*x)/(s0*s0 + sx*sx)


      if (allocated(Basis%x)) deallocate (Basis%x)
      if (allocated(Basis%w)) deallocate (Basis%w)

      if (allocated(Basis%d0gb)) deallocate (Basis%d0gb)
      if (allocated(Basis%d1gb)) deallocate (Basis%d1gb)
      if (allocated(Basis%d2gb)) deallocate (Basis%d2gb)
      if(allocated(Basis%S)) deallocate(Basis%S)


      allocate (Basis%x(nq))
      allocate (Basis%w(nq))
      allocate (w(nq))
      allocate (Q(nq))
      allocate(Basis%S(nb,nb))

      call hercom(nq, Basis%x(:), Basis%w(:))
      w(:) = Basis%w(:)/SQeq
      Q (:) = Qeq + Basis%x(:)/SQeq

      allocate (Basis%d0gb(nq, nb))
      allocate (Basis%d1gb(nq, nb, 1))
      allocate (Basis%d2gb(nq, nb, 1, 1))

      allocate (d0gb(nq, nb)) 
      allocate (d0bgw(nb, nq))

      DO iq = 1, nq
         DO ib = 1, nb

            call d0d1d2poly_Hermite_exp_cplx(Basis%x(iq),p/sx,ib-1,Basis%d0gb(iq, ib),&
             Basis%d1gb(iq, ib, 1),Basis%d2gb(iq, ib, 1, 1), .TRUE.)

             d0gb(iq, ib) =sqrt(s0)*poly_Hermite_exp_cplx(s0*(Q(iq)-Q0), p0/s0, ib-1)

             d0bgw(ib, iq)= conjg(sqrt(sx)*poly_Hermite_exp_cplx(sx*(Q(iq)-x), p/sx,ib -1)*w(iq))

         END DO
      END DO

      !Basis%S = matmul(d0bgw,d0gb)
      call Hagedorn_ovelp_mat(Basis%S,nb,nq,B1,B2)
      !call  Write_VecMat(Basis%S, out_unit, 5,  info='S')
      deallocate(d0gb,d0bgw,w,Q)

   END SUBROUTINE Construct_Basis_Hagedorn

   SUBROUTINE CheckOrtho_Basis(Basis, nderiv)
      USE QDUtil_m

      TYPE(Basis_t), intent(in)     :: Basis
      integer, intent(in)           :: nderiv

      integer                       :: ib, jb
      real(kind=Rkind), ALLOCATABLE    :: S(:, :)
      real(kind=Rkind)                 :: Sii, Sij

      IF (Basis%Basis_name == 'el') Then
         print *, 'This routine is .not. possible Basis el'
         RETURN
      END IF
      IF (Basis_IS_allocated(Basis)) THEN

         S = matmul(conjg(Basis%d0bgw), Basis%d0gb)
          do ib = 1,Basis%nb
          do jb = 1,Basis%nb
         write(115,*) ib,jb,S(ib,ib),S(ib,jb)
           end do
          end do
          IF (nderiv > -1) CALL   Write_VecMat(S,out_unit,5, info='S')
         Sii = ZERO
         Sij = ZERO

         DO ib = 1, Basis%nb
            IF (abs(S(ib, ib) - ONE) > Sii) Sii = abs(S(ib, ib) - ONE)
            S(ib, ib) = ZERO
              
         END DO
         Sij = maxval(S)
         write (out_unit, *) 'Sii-1,Sij', Sii, Sij

         IF (nderiv > 0) THEN
            write(out_unit,*)
            S = matmul(Basis%d0bgw, Basis%d1gb(:, :, 1))
            call   Write_VecMat(S,out_unit,5, info='<d0b|d1b>',Rformat='e13.4')
            call   Write_VecMat(S,out_unit,5, info='<d0b|d1b>')
         END IF

         IF (nderiv > 1) THEN
            write(out_unit,*)
            S = matmul(Basis%d0bgw, Basis%d2gb(:, :, 1, 1))
             call   Write_VecMat(S,out_unit,5, info='<d0b|d2b>',Rformat='e13.4')
             call   Write_VecMat(S,out_unit,5, info='<d0b|d1b>')
         END IF

      ELSE
         write (out_unit, *) ' WARNNING in CheckOrtho_Basis'
         write (out_unit, *) ' the basis is not allocated.'
      END IF

   END SUBROUTINE CheckOrtho_Basis

   SUBROUTINE Scale_Basis(Basis, x0, sx)
      USE QDUtil_m

      TYPE(Basis_t), intent(inout)  :: Basis
      real(kind=Rkind), intent(in)     :: x0, sx
      IF (Basis%nq == 0) RETURN
      IF (abs(sx) > ONETENTH**6 .AND. Basis_IS_allocated(Basis)) THEN

         Basis%x(:) = x0 + Basis%x(:)/sx
         Basis%w(:) = Basis%w(:)/sx

         Basis%d0gb(:, :) = Basis%d0gb(:, :)*sqrt(sx)
         Basis%d1gb(:, :, :) = Basis%d1gb(:, :, :)*sqrt(sx)*sx
         Basis%d2gb(:, :, :, :) = Basis%d2gb(:, :, :, :)*sqrt(sx)*sx*sx

      ELSE
         write (out_unit, *) ' ERROR in Scale_Basis'
         write (out_unit, *) ' sx is too small  or ...'
         write (out_unit, *) ' the basis is not allocated.'
         STOP 'ERROR in Scale_Basis'
      END IF

   END SUBROUTINE Scale_Basis

   SUBROUTINE Calc_tranpose_d0gb(Basis)
      USE QDUtil_m
      TYPE(Basis_t), INTENT(INOUT)     :: Basis
      !real(kind=Rkind),        INTENT(INOUT)  :: d0bgw(:,:)
      INTEGER                                 :: ib, iq,nb,nq

      nb= Basis%nb
      nq = Basis%nq

      if (allocated(Basis%d0bgw)) deallocate (Basis%d0bgw)
      if (Basis%Basis_name == 'el') then
         return
      end if

      allocate (Basis%d0bgw(nb, nq))

      Basis%d0bgw = transpose(Basis%d0gb)
      DO ib = 1, nb
         Basis%d0bgw(ib, :) = Basis%d0bgw(ib, :)*Basis%w(:)
      END DO

   END SUBROUTINE Calc_tranpose_d0gb

   SUBROUTINE Calc_dngg_grid(Basis)
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

      if (allocated(Basis%d1gg)) deallocate (Basis%d1gg)
      if (allocated(Basis%d2gg)) deallocate (Basis%d2gg)

      if (allocated(Basis%dagb)) deallocate (Basis%dagb)
      if (allocated(Basis%dq0gb)) deallocate (Basis%dq0gb)
      if (allocated(Basis%dp0gb)) deallocate (Basis%dp0gb)
       
      if (allocated(Basis%dagg)) deallocate (Basis%dagg)
      if (allocated(Basis%dq0gg)) deallocate (Basis%dq0gg)
      if (allocated(Basis%dp0gg)) deallocate (Basis%dp0gg)

      allocate (Basis%dagb(nq, nb))
      allocate (Basis%dq0gb(nq, nb))
      allocate (Basis%dp0gb(nq, nb))

      allocate (Basis%d1gg(nq, nq, 1))
      allocate (Basis%d2gg(nq,nq, 1, 1))


      allocate (Basis%dagg(nq, nq))
      allocate (Basis%dq0gg(nq, nq))
      allocate (Basis%dp0gg(nq, nq))

      do ib = 1,nb
         do iq = 1,nq

            Basis%dagb(iq, ib)  = (ONE/a)*(HALF*Basis%d0gb(iq, ib)+(Basis%x(iq)-Q0)*Basis%d1gb(iq, ib, 1))

            Basis%dq0gb(iq, ib) = -Basis%d1gb(iq, ib, 1)

            Basis%dp0gb(iq, ib) =EYE*a*(Basis%x(iq)-Q0)*Basis%d0gb(iq, ib)

         end do
      end do


      IF (debug) THEN
         CALL   Write_VecMat(Basis%d0bgw(:, :), out_unit, 5,  info='d0bgw')
         write (out_unit, *)
      END IF

      if (Basis%Basis_name == 'el') THEN
         RETURN
      end if

      Basis%d1gg(:, :, 1) = matmul(Basis%d1gb(:, :, 1), conjg(Basis%d0bgw))
      Basis%d2gg(:, :, 1, 1) = matmul(Basis%d2gb(:, :, 1, 1), conjg(Basis%d0bgw))

      Basis%dagg(:, :) = matmul(Basis%dagb, conjg(Basis%d0bgw))
      Basis%dq0gg(:, :) = matmul(Basis%dq0gb,conjg(Basis%d0bgw))
      Basis%dp0gg(:, :) = matmul(Basis%dp0gb,conjg(Basis%d0bgw))

      IF (debug) THEN
         CALL   Write_VecMat(Basis%d1gg(:, :, 1), out_unit, 5,  info='d1gg')
         write (out_unit, *)
         CALL   Write_VecMat(Basis%d2gg(:, :, 1, 1), out_unit, 5,  info='d2gg')
         write (out_unit, *)
         CALL   Write_VecMat(Basis%dp0gg(:, :), out_unit, 5,  info='dp0gg')
         write (out_unit, *)
         CALL   Write_VecMat(Basis%dq0gg(:, :), out_unit, 5,  info='dq0gg')
         write (out_unit, *)
         CALL   Write_VecMat(Basis%dagg(:, :), out_unit, 5,  info='dagg')
         write (out_unit, *) 'END Calc_dngg_grid'
         flush (out_unit)
      END IF

   END SUBROUTINE Calc_dngg_grid

   SUBROUTINE test_basitogridgridtobasis(Basis)
      USE QDUtil_m
      TYPE(Basis_t), intent(in)       :: Basis
      logical, parameter        :: debug = .true.
      COMPLEX(kind=Rkind), allocatable   :: G1(:), B1(:)!,G(:), Hpsi(:)
      COMPLEX(kind=Rkind), allocatable   :: G2(:), B2(:)
      COMPLEX(kind=Rkind), allocatable   :: B(:)
      REAL(kind=Rkind), allocatable       :: diff_g(:), diff_b(:)
      !REAL(KIND=Rkind)                      :: Norm0,Norm1,min_diff,max_diff
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
      print *, 'G in', G1
      Call GridTOBasis_nD_cplx(B, G1, Basis)
      Call BasisTOGrid_nD_cplx(G2, B, Basis)
      print *, 'G out', G2

      print *, '---------------------------------------------------------------------'

      diff_g(:) = ABS(G1(:) - G2(:))

      print *, '---------------------------------------------------------------------'
      Write (out_unit, *) 'maxval(diff_g(:))=', maxval(diff_g(:))
      Write (out_unit, *) 'MINVAL(diff_g(:))=', MINVAL(diff_g(:))
      print *, '---------------------------------------------------------------------'
      G2(:) = CZERO
      B2(:) = CZERO
      B1(:) = CONE
      print *, 'B in', B1
      Call BasisTOGrid_nD_cplx(G2, B1, Basis)
      Call GridTOBasis_nD_cplx(B2, G2, Basis)
      print *, 'B out', B2
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

     SUBROUTINE GridTOBasis_1D_cplx(BB, GG, Basis)
      USE QDUtil_m
      TYPE(Basis_t), intent(in), target              :: Basis
      complex(kind=Rkind), intent(inout)             :: BB(:, :, :)
      complex(kind=Rkind), intent(in)                :: GG(:, :, :)
      logical, parameter                             :: debug = .true.
      Integer                                        :: i1, i3


      IF (debug) THEN
         flush (out_unit)
      END IF

      BB = CZERO
      DO i3 = 1, ubound(GG, dim=3)
      DO i1 = 1, ubound(GG, dim=1)

         BB(i1, :, i3) = matmul(conjg(Basis%d0bgw),GG(i1, :, i3))

      END DO
      END DO

      IF (debug) THEN
         flush (out_unit)
      END IF
   END SUBROUTINE

   SUBROUTINE BasisTOGrid_1D_cplx(GB, BB, Basis)
      USE QDUtil_m
      TYPE(Basis_t), intent(in), target           :: Basis
      complex(kind=Rkind), intent(inout)             :: GB(:, :, :)
      complex(kind=Rkind), intent(in)                :: BB(:, :, :)
      logical, parameter                          :: debug = .true.
      integer                                     :: i1, i3, iq, ib

      IF (debug) THEN
         flush (out_unit)
      END IF

      GB = CZERO
      DO i3 = 1, ubound(BB, dim=3)
      DO i1 = 1, ubound(BB, dim=1)

         GB(i1, :, i3) = matmul(Basis%d0gb(:,:),BB(i1, :, i3))

      END DO
      END DO

      IF (debug) THEN
         flush (out_unit)
      END IF
   END SUBROUTINE

   SUBROUTINE Calc_index(Ib1, Ib2, Ib3, Iq1, Iq2, Iq3, Basis)
      TYPE(Basis_t), intent(in), target           :: Basis
      integer, intent(inout), allocatable, optional  :: Iq1(:), Iq2(:), Iq3(:)
      integer, intent(inout), allocatable, optional  :: Ib1(:), Ib2(:), Ib3(:)
      integer                                        :: Ndim
      integer                                        :: inb

      Ndim = size(Basis%tab_basis) - 1

      if (present(Ib3)) allocate (Ib3(Ndim))
      if (present(Ib2)) allocate (Ib2(Ndim))
      if (present(Ib1)) allocate (Ib1(Ndim))

      if (present(Iq3)) allocate (Iq3(Ndim))
      if (present(Iq2)) allocate (Iq2(Ndim))
      if (present(Iq1)) allocate (Iq1(Ndim))

      DO inb = 1, Ndim

         IF (inb == 1) THEN

            if (present(Iq1)) Iq1(1) = 1
            if (present(Ib1)) Ib1(1) = 1

            if (present(Iq2)) Iq2(1) = Basis%tab_basis(1)%nq
            if (present(Ib2)) Ib2(1) = Basis%tab_basis(1)%nb

            if (present(Iq3)) Iq3(1) = Product(Basis%tab_basis(2:Ndim)%nq)*Basis%tab_basis(Ndim + 1)%nb
            if (present(Ib3)) Ib3(1) = Product(Basis%tab_basis(2:Ndim + 1)%nb)

         ELSE IF (inb == Ndim) THEN

            if (present(Iq1)) Iq1(inb) = Product(Basis%tab_basis(1:Ndim - 1)%nq)
            if (present(Ib1)) Ib1(inb) = Product(Basis%tab_basis(1:Ndim - 1)%nb)

            if (present(Ib2)) Ib2(inb) = Basis%tab_basis(Ndim)%nb
            if (present(Iq2)) Iq2(inb) = Basis%tab_basis(Ndim)%nq

            if (present(Ib3)) Ib3(inb) = Basis%tab_basis(Ndim + 1)%nb
            if (present(Iq3)) Iq3(inb) = Basis%tab_basis(Ndim + 1)%nb
         ELSE

            if (present(Iq1)) Iq1(inb) = Product(Basis%tab_basis(1:inb - 1)%nq)
            if (present(Ib1)) Ib1(inb) = Product(Basis%tab_basis(1:inb - 1)%nb)

            if (present(Iq2)) Iq2(inb) = Basis%tab_basis(inb)%nq
            if (present(Ib2)) Ib2(inb) = Basis%tab_basis(inb)%nb

            if (present(Ib3)) Ib3(inb) = Product(Basis%tab_basis(inb + 1:Ndim + 1)%nb)
            if (present(Iq3)) Iq3(inb) = Product(Basis%tab_basis(inb + 1:Ndim)%nq)*Basis%tab_basis(Ndim + 1)%nb

         END IF
      END DO

   END SUBROUTINE

   SUBROUTINE BasisTOGrid_nD_cplx(G, B, Basis)
      USE QDUtil_m
      USE NDindex_m
      !Logical,           parameter                 :: debug = .true.
      Logical, parameter                            :: debug = .false.
      TYPE(Basis_t), intent(in)                     :: Basis
      complex(kind=Rkind), intent(in), target          :: B(:) !Vector on base,
      complex(kind=Rkind), intent(inout), target       :: G(:) !Vector on the grid, out
      complex(kind=Rkind), pointer                     :: BBG(:, :, :)
      complex(kind=Rkind), pointer                     :: BBB(:, :, :), GGB(:, :, :)
      integer, allocatable                             :: Ib3(:), Iq1(:), Iq2(:), ib1(:), ib2(:), iq3(:)
      complex(kind=Rkind), pointer                     :: GBB(:, :, :)
      complex(kind=Rkind), allocatable, target         :: GBB1(:), GGB2(:)

      Integer                                       :: ib, iq, nq, nb, inb, Ndim
      Integer                                       :: jb, jb1, jb2

      IF (debug) THEN
         write (out_unit, *) 'BEGINNING BasisTOGrid_nD_cplx'
         write (out_unit, *) 'intent(in) :: B(:)', B
         Call Write_Basis(Basis)
         flush (out_unit)
      END IF

      IF (.NOT. Basis_IS_Allocated(Basis)) THEN
         write (out_unit, *) ' ERROR in BasisTOGrid_Basis'
         write (out_unit, *) " the basis is not Allocated."
         STOP "ERROR BasisTOGrid_Basis: the basis is not Allocated."
      END IF

      IF (size(B) /= Basis%nb*Basis%tab_basis(size(Basis%tab_basis))%nb) THEN
         write (out_unit, *) ' ERROR in BasisTOGrid_Basis'
         write (out_unit, *) ' the size of B is different from nb.'
         write (out_unit, *) ' size(B), Basis%nb', size(B), Basis%nb*Basis%tab_basis(size(Basis%tab_basis))%nb
         STOP 'ERROR in BasisTOGrid_Basis: wrong B size.'
      END IF

      IF (size(G) /= Basis%nq*Basis%tab_basis(size(Basis%tab_basis))%nb) THEN
         write (out_unit, *) ' ERROR in GridTOBasis_Basis'
         write (out_unit, *) ' the size of G is different from nq.'
         write (out_unit, *) ' size(G), Basis%nq', size(G), Basis%nq*Basis%tab_basis(size(Basis%tab_basis))%nb
         STOP 'ERROR in BasisTOGrid_Basis: wrong G size..'
      END IF

      Ndim = size(Basis%tab_basis) - 1
      Call Calc_index(Ib1, Ib2, Ib3, Iq1, Iq2, Iq3, Basis)
      if (Ndim == 1) then
         BBB(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => B
         GGB(1:Iq1(1), 1:Iq2(1), 1:Iq3(1)) => G
         G = CZERO
         Call BasisTOGrid_1D_cplx(GGB, BBB, Basis%tab_basis(1))
      else

         Allocate (GBB1(iq1(1)*iq2(1)*ib3(1)))
         BBB(1:iq1(1), 1:ib2(1), 1:ib3(1)) => B
         GBB(1:iq1(1), 1:iq2(1), 1:ib3(1)) => GBB1
         GBB1 = CZERO

         Call BasisTOGrid_1D_cplx(GBB, BBB, Basis%tab_basis(1))

         DO inb = 2, Ndim - 1

            Allocate (GGB2(iq1(inb)*iq2(inb)*ib3(inb)))

            BBB(1:iq1(Inb), 1:ib2(inb), 1:ib3(inb)) => GBB1

            GBB(1:iq1(inb), 1:iq2(inb), 1:ib3(inb)) => GGB2

            Call BasisTOGrid_1D_cplx(GBB, BBB, Basis%tab_basis(inb))

            GBB1 = GGB2

            Deallocate (GGB2)

         END DO

         BBB(1:iq1(Ndim), 1:ib2(Ndim), 1:ib3(Ndim)) => GBB1

         GBB(1:iq1(Ndim), 1:iq2(Ndim), 1:ib3(Ndim)) => G

         Call BasisTOGrid_1D_cplx(GBB, BBB, Basis%tab_basis(Ndim))
         Deallocate (Ib1, Iq1, Iq2, Ib2, Ib3, Iq3)

      END IF

      IF (debug) THEN
         write (out_unit, *) 'intent(INOUT) :: G(:)', G
         write (out_unit, *) 'END BasisTOGrid_nD_cplx'
         flush (out_unit)
      END IF
   END SUBROUTINE

   SUBROUTINE GridTOBasis_nD_cplx(B, G, Basis)
      USE QDUtil_m
      !Logical,          parameter                     :: debug = .true.
      Logical, parameter                               :: debug = .false.
      TYPE(Basis_t), intent(in), target                :: Basis
      complex(kind=Rkind), intent(in), target             :: G(:)
      complex(kind=Rkind), intent(inout), target          :: B(:)
      complex(kind=Rkind), pointer                        :: BBB(:, :, :), GGB(:, :, :)
      complex(kind=Rkind), allocatable, target            :: BGG1(:), BGG2(:)
      complex(kind=Rkind), pointer                        :: GGG(:, :, :)
      Integer                                          :: ib, i1, i3, inb, Ndim, iq
      Integer, allocatable                             :: Ib1(:), Ib2(:), Iq3(:), Iq1(:), Iq2(:), Ib3(:)

      IF (debug) THEN
         write (out_unit, *) 'BEGINNING GridTOBasis_nD_cplx'
         write (out_unit, *) 'intent(in) :: G(:)', G
         !Call Write_Basis(Basis)
         flush (out_unit)
      END IF

      IF (.NOT. Basis_IS_Allocated(Basis)) THEN
         write (out_unit, *) ' ERROR in BasisTOGrid_nD_cplx'
         write (out_unit, *) " the basis is not Allocated."
         STOP "ERROR BasisTOGrid_Basis: the basis is not Allocated."
      END IF

      IF (size(B) /= Basis%nb*Basis%tab_basis(size(Basis%tab_basis))%nb) THEN
         write (out_unit, *) ' ERROR in BasisTOGrid_nD_cplx'
         write (out_unit, *) ' the size of G is different from nb.'
         write (out_unit, *) ' size(B), Basis%nb', size(B), Basis%nb*Basis%tab_basis(size(Basis%tab_basis))%nb
         STOP 'ERROR in GridTOBasis_Basis: wrong B size.'
      END IF

      IF (size(G) /= Basis%nq*Basis%tab_basis(size(Basis%tab_basis))%nb) THEN
         write (out_unit, *) ' ERROR in GridTOBasis_Basis'
         write (out_unit, *) ' the size of G is different from nq.'
         write (out_unit, *) ' size(G), Basis%nq', size(G), Basis%nq*Basis%tab_basis(size(Basis%tab_basis))%nb
         STOP 'ERROR in GridTOBasis_Basis: wrong G size'
      END IF

      Ndim = size(Basis%tab_basis) - 1

      Call Calc_index(Ib1, Ib2, Ib3, Iq1, Iq2, Iq3, Basis)
      If (Ndim == 1) then
         B = CZERO
         GGB(1:Iq1(1), 1:Iq2(1), 1:Iq3(1)) => G
         BBB(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => B
         Call GridTOBasis_1D_cplx(BBB, GGB, Basis%tab_basis(1))
      else

         Allocate (BGG1(Ib1(1)*Ib2(1)*Iq3(1)))
         BGG1 = CZERO
         GGG(1:Ib1(1), 1:Iq2(1), 1:Iq3(1)) => G
         GGB(1:Ib1(1), 1:Ib2(1), 1:Iq3(1)) => BGG1

         Call GridTOBasis_1D_cplx(GGB, GGG, Basis%tab_basis(1))

         DO inb = 2, Ndim - 1

            Allocate (BGG2(Ib1(inb)*Ib2(inb)*Iq3(inb)))

            BGG2 = CZERO

            GGG(1:Ib1(inb), 1:Iq2(inb), 1:Iq3(inb)) => BGG1
            GGB(1:Ib1(inb), 1:Ib2(inb), 1:Iq3(inb)) => BGG2

            Call GridTOBasis_1D_cplx(GGB, GGG, Basis%tab_basis(inb))

            BGG1 = BGG2

            deallocate (BGG2)
         END DO

         B(:) = CZERO

         GGG(1:Ib1(Ndim), 1:Iq2(Ndim), 1:Iq3(Ndim)) => BGG1
         GGB(1:Ib1(Ndim), 1:Ib2(Ndim), 1:Iq3(Ndim)) => B

         Call GridTOBasis_1D_cplx(GGB, GGG, Basis%tab_basis(Ndim))

         Deallocate (Iq1, Iq2, Iq3, Ib1, Ib2, Ib3)
      END IF
      IF (debug) THEN
         write (out_unit, *) 'END GridTOBasis_nD_cplx'
         flush (out_unit)
      END IF
   END SUBROUTINE

   SUBROUTINE Calc_Q_grid(Q, Basis, WnD)

      implicit none
      TYPE(Basis_t), intent(in)                                    :: Basis
      integer, ALLOCATABLE                                         :: Tab_iq(:), NDend(:)
      integer                                                      :: inb, ndim, iq
      real(Kind=Rkind), intent(inout), allocatable, optional          :: Q(:, :)
      real(Kind=Rkind), intent(inout), allocatable, optional          :: WnD(:)
      TYPE(NDindex_t)                                              :: NDindex
      logical                                                      :: Endloop
      ndim = SIZE(Basis%tab_basis) - 1
      allocate (Tab_iq(Ndim))
      allocate (NDend(Ndim))
      if (present(Q)) allocate (Q(Basis%nq, Ndim))
      if (present(WnD)) allocate (WnD(Basis%nq))
      do inb = 1, Ndim
         NDend(inb) = Basis%NDindexq%NDend(inb)
      end do
      CALL Init_NDindex(NDindex, NDend, Ndim)
      Call Init_tab_ind(Tab_iq, NDindex)
      Iq = 0
      DO
         Iq = Iq + 1
         CALL increase_NDindex(Tab_iq, NDindex, Endloop)
         IF (Endloop) exit
         if (present(WnD)) WnD(iq) = ONE
         do inb = 1, Ndim
            if (present(Q)) Q(iq, inb) = Basis%tab_basis(inb)%X(Tab_iq(inb))
            if (present(WnD)) WnD(iq) = WnD(iq)*Basis%tab_basis(inb)%w(Tab_iq(inb))
         end do
         !  if (present(Q))   print*,iq,Q(iq,:)
         ! if (present(WnD))  print*,iq,WnD(iq)
      END DO
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
      Hfi =  sqrt(B1(3))*poly_Hermite_exp_cplx(qi,B1(2)/B1(3), i- 1)

      ! Hfj = He(J,qj) = HermiteH[J,qj/Sqrt[2]] / Sqrt [ 2^J ]
       qj = (x - B2(1))*B2(3)
       Hfj= sqrt(B2(3))*poly_Hermite_exp_cplx(qj,B2(2)/B2(3), j- 1)
       

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
         call  Hermite_double_product_func(f(iq), Q(iq), B1,B2, j, i)
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
            call Calc_Hermitec_int(Hf, Q, B1,B2, jb, ib)
            Mat(jb,ib) = dot_product(Hf, w)
         End Do
      End Do

      deallocate (Hf,q,w)
      print*,"first basis parameters",B1
      print*,"second basis parameters",B2
      CALL  Write_VecMat(Mat, out_unit, 5,  info='Mat')

   END SUBROUTINE 
   


 
SUBROUTINE Calc_daG(daG, G, Basis)
   USE  QDUtil_m
   TYPE(Basis_t), intent(in)              :: Basis
   complex(kind=Rkind), intent(in)        :: G(:)
   complex(kind=Rkind), intent(inout)     :: daG(:)
   logical, parameter                     :: debug = .false.

   
   IF (debug) THEN
      write(out_unit,*) 'BEGINNING daG'
      flush (out_unit)
   END IF
          
            daG =  matmul(Basis%tab_basis(1)%dagg, G)
    
   IF (debug) THEN
              write(out_unit,*) 'END daG'
      flush (out_unit)
   END IF
END SUBROUTINE 


SUBROUTINE Calc_dq0G(dq0G, G, Basis)
   USE  QDUtil_m
   TYPE(Basis_t), intent(in)              :: Basis
   complex(kind=Rkind), intent(in)        :: G(:)
   complex(kind=Rkind), intent(inout)     :: dq0G(:)
   logical, parameter                     :: debug = .false.

   
   IF (debug) THEN
      write(out_unit,*) 'BEGINNING daG'
      flush (out_unit)
   END IF
          
            dq0G =  matmul(Basis%tab_basis(1)%dq0gg, G)
    
   IF (debug) THEN
              write(out_unit,*) 'END daG'
      flush (out_unit)
   END IF
END SUBROUTINE 



SUBROUTINE Calc_dp0G(dp0G, G, Basis)
   USE  QDUtil_m
   TYPE(Basis_t), intent(in)              :: Basis
   complex(kind=Rkind), intent(in)        :: G(:)
   complex(kind=Rkind), intent(inout)     :: dp0G(:)
   logical, parameter                     :: debug = .false.

   
   IF (debug) THEN
      write(out_unit,*) 'BEGINNING dp0G'
      flush (out_unit)
   END IF
          
            dp0G =  matmul(Basis%tab_basis(1)%dp0gg, G)
    
   IF (debug) THEN
              write(out_unit,*) 'END dp0G'
      flush (out_unit)
   END IF
END SUBROUTINE 


END MODULE Basis_m
