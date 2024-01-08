MODULE Basis_m
   USE QDUtil_m
   USE NDindex_m
   USE polyortho_m

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: Basis_t, Read_Basis, Write_Basis, Basis_IS_allocated, Deallocate_Basis, Basis_IS_allocatedtot
   PUBLIC :: Calc_dngg_grid, Calc_tranpose_d0gb, test_basitogridgridtobasis
   PUBLIC :: GridTOBasis_nD_cplx, BasisTOGrid_nD_cplx
   PUBLIC :: Calc_Q_grid, Calc_index,Rdensity_alloc
   PUBLIC :: Scale_Basis, init_Basis1_TO_Basis2,Calc_dngg_grid_0
   PUBLIC :: construct_primitive_basis,CheckOrtho_Basis
   PUBLIC :: pdv2psi_nD,Get_Basis_Parameters,Calc_Basis_dPtQtAt
   PUBLIC :: Calc_reduced_Density_surf,Calc_reduced_density,Rdensity_Writing, REDUCED_DENSIRY_t
   PUBLIC :: Hagedorn_construction,Construct_Basis_Ho,Construct_Hagedorn_Variational_Basis
   PUBLIC :: Change_Basis_Parameters,Complete_Hagedorn_none_variationnal_Basis,Calc_d0d1d2W
   PUBLIC :: Complete_Hagedorn_none_variationnal_Basis_temp,construct_primitive_basis_temp
   PUBLIC :: Weight_D_smol,Calc_n_smol,Calc_nterm_compact,Testsmol,Calc_nterm

   TYPE :: Basis_t
      integer                             :: nb_basis = ZERO
      integer                             :: nb ,nq
      integer                             :: A_smol, B_smol,LB,LG
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
      TYPE(NDindex_t)                     :: NDindexlq
      TYPE(NDindex_t)                     :: NDindexlb
      TYPE(Basis_t),      allocatable     :: tab_basis(:)     
      TYPE(Basis_t),      allocatable     :: tab_Smolyak(:)     

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
            call init_Basis1_TO_Basis2(Basis2%tab_basis(ib), Basis1%tab_basis(ib))
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
         Basis2%Alpha            = Basis1%Alpha
         Basis2%A                = Basis1%A
         Basis2%B                = Basis1%B
         Basis2%S                = Basis1%S
      END IF

   END SUBROUTINE init_Basis1_TO_Basis2

   RECURSIVE SUBROUTINE Deallocate_Basis(Basis)
      USE QDUtil_m
      TYPE(Basis_t), intent(inout)     :: Basis
      integer                          :: i

      ! write(out_unit,*) '-------------------------------------------------------------------------------'
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
          write(out_unit,*)' Basis table x is now deallocated.'
      END IF


      write (out_unit, *)
      IF (allocated(Basis%d0gb)) THEN
         deallocate (Basis%d0gb)
          write(out_unit,*)' Basis table d0gb is now deallocated.'
      END IF

      write (out_unit, *)
      IF (allocated(Basis%d0bgw)) THEN
         deallocate (Basis%d0bgw)
          write(out_unit,*)' Basis table d0bgw is now deallocated.'
      END IF
        
      write (out_unit, *)
      IF (allocated(Basis%dagb)) THEN
         deallocate (Basis%dagb)
         write(out_unit,*)' Basis table dagb is now deallocated.'
      END IF

      write (out_unit, *)
      IF (allocated(Basis%dagg)) THEN
         deallocate (Basis%dagg)
         write(out_unit,*)' Basis table dagg is now deallocated.'
      END IF


      write (out_unit, *)
      IF (allocated(Basis%dq0gb)) THEN
         deallocate (Basis%dq0gb)
         write(out_unit,*)' Basis table dq0gb is now deallocated.'
      END IF

      
      write (out_unit, *)
      IF (allocated(Basis%dq0gb)) THEN
         deallocate (Basis%dq0gb)
          write(out_unit,*)' Basis table dq0gb is now deallocated.'
      END IF

      write (out_unit, *)
      IF (allocated(Basis%dq0gg)) THEN
         deallocate (Basis%dq0gg)
         write(out_unit,*)' Basis table dq0gg is now deallocated.'
      END IF


      write (out_unit, *)
      IF (allocated(Basis%dp0gg)) THEN
         deallocate (Basis%dp0gg)
         write(out_unit,*)' Basis table dp0gg is now deallocated.'
      END IF


      write (out_unit, *)
      IF (allocated(Basis%d1gb)) THEN
         deallocate (Basis%d1gb)
          write(out_unit,*)' Basis table d1gb is now deallocated.'
      END IF
      write (out_unit, *)
      IF (allocated(Basis%d1gg)) THEN
         deallocate (Basis%d1gg)
          write(out_unit,*)' Basis table d1gb is now deallocated.'
      END IF

      write (out_unit, *)
      IF (allocated(Basis%d2gb)) THEN
         deallocate (Basis%d2gb)
         write(out_unit,*)' Basis table d1gb is now deallocated.'
      END IF
      write (out_unit, *)
      IF (allocated(Basis%d2gg)) THEN
         deallocate (Basis%d2gg)
           write(out_unit,*)' Basis table d2gg is now deallocated.'
      END IF
      IF (allocated(Basis%tab_basis)) THEN
         DO i = 1, size(Basis%tab_basis)
            CALL Deallocate_Basis(Basis%tab_basis(i))
         END DO
      END IF
      ! write(out_unit,*) '--------------------------------------------------------------------------------'

   END SUBROUTINE Deallocate_Basis

   RECURSIVE SUBROUTINE Read_Basis(Basis, nio)
      USE QDUtil_m
      logical, parameter                       :: debug = .true.
      !logical,             parameter          ::debug = .false.
      TYPE(Basis_t), intent(inout)             :: Basis
      integer, intent(in)                      :: nio
      integer                                  :: err_io,nb,nq,i,nb_basis,ib
      character(len=Name_len)                  :: name
      real(kind=Rkind)                         :: A, B, scaleQ, Q0,Imp_k
      integer                                  :: A_smol, B_smol,LB,LG
      complex(kind=Rkind)                      :: Alpha

      NAMELIST /basis_nD/ name, nb_basis, A_smol, B_smol,LB,LG ,nb, nq, A, B, scaleQ, Q0,Imp_k,Alpha
      nb_basis = 0
      nb = 0
      nq = 0
      A = ZERO
      B = ZERO
      A_smol=0
      B_smol=0
      LB=0
      LG=0
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
         Basis%Basis_name = name
         Basis%nb_basis = nb_basis
         call string_uppercase_TO_lowercase(Basis%Basis_name)
         If(Basis%Basis_name == 'dp') then
           allocate (Basis%tab_basis(nb_basis))
           DO i = 1, nb_basis
              CALL Read_Basis(Basis%tab_basis(i), nio)
           END DO
           Basis%NDindexb%sys_type = 'dp'
           Basis%NDindexq%sys_type = 'dp'
           Basis%nb = 1
           Basis%nq = 1

           DO i = 1, nb_basis
              if (Basis%tab_basis(i)%Basis_name == 'el') cycle
              Basis%nb = Basis%nb*Basis%tab_basis(i)%nb
              Basis%nq = Basis%nq*Basis%tab_basis(i)%nq
           END DO

         else if(Basis%Basis_name == 'smolyak') then
            Basis%NDindexlb%L = LB
            Basis%NDindexlb%Ndim = nb_basis-1
            Basis%NDindexlq%Ndim = nb_basis-1
            allocate(Basis%NDindexlq%NDend(nb_basis-1))
            allocate(Basis%NDindexlb%NDend(nb_basis-1))
            Basis%NDindexlq%NDend(:) =LG
            Basis%NDindexlb%NDend(:) =LB
            Basis%NDindexlb%sys_type = 'smolyak'
            Basis%NDindexlq%sys_type = 'smolyak'
           allocate (Basis%tab_Smolyak(nb_basis))
           DO i = 1, nb_basis
           CALL Read_Basis(Basis%tab_Smolyak(i), nio)
           END DO

         end if
      ELSE
         Basis%nb_basis = nb_basis
         Basis%nb = nb
         Basis%nq = nq
         Basis%Q0 = Q0
         Basis%SCALEQ = SCALEQ
         Basis%A = A
         Basis%B = B
         Basis%A_smol = A_smol
         Basis%B_smol = B_smol
         Basis%LB = LB
         Basis%LG = LG
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

  FUNCTION Calc_n_smol(A_smol,B_smol,l)
    Integer, intent(in)     :: A_smol
    Integer, intent(in)     :: B_smol
    Integer, intent(in)     :: l
    Integer                 :: Calc_n_smol
 
     Calc_n_smol = A_smol + B_smol*l

  END FUNCTION 

   SUBROUTINE Weight_D_smol ( D,L,Som_l,ndim)
     IMPLICIT NONE
     Integer,   intent(in)       :: ndim
     Integer,   intent(in)       :: Som_l
     Integer,   intent(in)       :: L
     Integer,   intent(inout)    :: D
 
     D = ((-1)**(L-Som_l))*binomial(L-Som_l,ndim-1)

   END SUBROUTINE 


   SUBROUTINE Calc_nterm_compact ( Ncomp,L,ndim)
    IMPLICIT NONE
    Integer,   intent(in)       :: ndim
    Integer,   intent(in)       :: L
    Integer,   intent(inout)    :: Ncomp

      Ncomp = Factorial(L+ndim)/(Factorial(L)*Factorial(ndim))

   END SUBROUTINE 


 SUBROUTINE Calc_nterm ( Nterm,Basis)
  IMPLICIT NONE
  TYPE(Basis_t), intent(inout)   :: Basis
  Integer,   intent(inout)       :: Nterm
  integer                        :: I,sum_l
  logical                        :: Endloop
  Integer,allocatable            :: Tab_lb(:),D_lb(:)


  Allocate(Tab_lb(Basis%NDindexlb%Ndim))
  call Init_NDindex(Basis%NDindexlb, Basis%NDindexlb%NDend, Basis%NDindexlb%Ndim)
  Call Init_tab_ind(Tab_lb,Basis%NDindexlb)
I=0
DO 
  I=I+1
  CALL increase_NDindex(Tab_lb,Basis%NDindexlb,Endloop)
   write(out_unit,*) i, Tab_lb(:)
    IF (Endloop) exit
 
END DO
Nterm = I
print*,'I',I

allocate(D_lb(Nterm))

Call Init_tab_ind(Tab_lb,Basis%NDindexlb)
I=0
DO 
  I=I+1
  CALL increase_NDindex(Tab_lb,Basis%NDindexlb,Endloop)
  sum_l = sum(Tab_lb)
  !call Weight_D_smol ( D_lb(I),Basis%NDindexlb%L,Sum_l,Basis%NDindexlb%Ndim)
   !write(out_unit,*) i, D_lb(:)
   write(out_unit,*) i,Basis%NDindexlb%L-sum_l ,Basis%NDindexlb%Ndim-1
    IF (Endloop) exit
 
END DO
 END SUBROUTINE 


   RECURSIVE SUBROUTINE construct_primitive_basis_temp(Basis)
   USE QDUtil_m
   logical, parameter                     :: debug = .true.
   !logical,             parameter        ::debug = .false.
   TYPE(Basis_t), intent(inout)           :: Basis
   integer                                :: ndim, ib
    ndim = Basis%nb_basis 
   IF (allocated(Basis%tab_basis)) THEN
      DO ib = 1, ndim
         CALL construct_primitive_basis_temp(Basis%tab_basis(ib))
      END DO
   ELSE
      SELECT CASE (Basis%Basis_name)
      CASE ('el')
         write (6, *) 'Electronic basis. Electronic state number:', basis%nb
         basis%nq = 0
      CASE ('boxab')
         call Construct_Basis_Sin(Basis)
         Basis%Q0 = Basis%A
         Basis%scaleQ = pi/(Basis%B - Basis%A)
         call Hagedorn_construction(Basis)
      CASE ('fourier')
         call Construct_Basis_Fourier(Basis)
         call Hagedorn_construction(Basis)
      CASE ('herm', 'ho')
         call Construct_Basis_Ho(Basis)
         call Complet_scaling_Basis(Basis)
      CASE default
         STOP 'ERROR  Noting to construct'
      END SELECT
      !  this part wil not have sens for 'el' basis
   END IF
   ! write(out_unit,*) ' End  construct  primitive Basis '
END SUBROUTINE 


 SUBROUTINE Testsmol(NDindex)
   USE QDUtil_m
   IMPLICIT NONE
   TYPE(NDindex_t),  intent(inout) :: NDindex
   integer, allocatable         :: Tab_ind(:)
   logical                      :: Endloop
 !  integer,     intent(in)      :: nl
   integer                      :: i
   !logical,    parameter      :: debug = .true.
   logical,     parameter      :: debug = .false.
   IF (debug) THEN
     write(out_unit,*) 'BEGINNING Testsmol'
     flush(out_unit)
   END IF
   Allocate(Tab_ind(NDindex%Ndim))
  ! write(out_unit,*) 'NDindex%NDim',NDindex%Ndim
   ! write(out_unit,*) 'NDindex%NDend',NDindex%NDend
   !STOP 'cc'
   call Init_NDindex(NDindex, NDindex%NDend, NDindex%Ndim)
   Call Init_tab_ind(Tab_ind,NDindex)
   i=0
  ! write(out_unit,*) i, tab_ind(:)
   !STOP 'cc'
   DO !i= 1,100
     i=i+1
     CALL increase_NDindex(Tab_ind,NDindex,Endloop)
      write(out_unit,*) i, tab_ind(:),sum(tab_ind(:))
       IF (Endloop) exit
    
   END DO
   
   IF (debug) THEN
     write(out_unit,*) 'END Testsmol'
     flush(out_unit)
   END IF
 END SUBROUTINE 

   RECURSIVE SUBROUTINE construct_primitive_basis(Basis)
      USE QDUtil_m
      logical, parameter                     :: debug = .true.
      !logical,             parameter        ::debug = .false.
      TYPE(Basis_t), intent(inout)           :: Basis
      integer, allocatable                   :: NDend_q(:)
      integer, allocatable                   :: NDend_b(:)
      integer, allocatable                   :: NDend_lb(:)
      integer, allocatable                   :: NDend_lq(:)
      integer                                :: ndim, ib

       ndim = Basis%nb_basis - 1

        SELECT CASE (Basis%Basis_name)
        CASE ('dp')
         allocate (NDend_q(ndim))
         allocate (NDend_b(ndim)) 
          DO ib = 1, ndim
            NDend_q(ib) = Basis%tab_basis(ib)%nq
            NDend_b(ib) = Basis%tab_basis(ib)%nb
          END DO 

         CASE ('smolyak')
         allocate (NDend_lb(ndim))
         allocate (NDend_lq(ndim))
          DO ib = 1, ndim
             NDend_lq(ib) = Basis%LG
             NDend_lb(ib) = Basis%LB
         END DO 

        END SELECT

      SELECT CASE (Basis%Basis_name)
        CASE ('dp')

         call Init_NDindex(Basis%NDindexq, NDend_q, ndim)
         call Init_NDindex(Basis%NDindexb, NDend_b, ndim)

        CASE('smolyak')

           Basis%NDindexlq%L  = Basis%LG
           Basis%NDindexlb%L  = Basis%LB
           call Init_NDindex(Basis%NDindexlq, NDend_lq, ndim)
           call Init_NDindex(Basis%NDindexlb, NDend_lb, ndim)

       End SELECT

       !call construct_primitive_basis_temp(Basis)
        
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

   SUBROUTINE Construct_Basis_Fourier(Basis) !basis_name fourier[-pi,pi]
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
      dx = TWO*pi/nq

      !>------------------grid and weight -----------------------------------------
      Basis%x = [(iq*dx - dx/2 - pi, iq=1, nq)]
      Basis%w = [(dx, iq=1, nq)]
      !> ---------------------allocation ---------------------------------------
      allocate (Basis%d0gb(nq, nb))
      allocate (Basis%d1gb(nq, nb, 1))
      allocate (Basis%d2gb(nq, nb, 1, 1))

      DO ib = 1, nb
         DO iq = 1, nq
            call d0d1d2fourier(Basis%x(iq),Basis%d0gb(iq, ib),Basis%d1gb(iq, ib, 1)&
                           &,Basis%d2gb(iq, ib, 1, 1),ib)
         END DO                  
      END DO
      !---------------------------------------------------------------------------
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
      integer                                    :: nb, nq
     
      nb    = Basis%nb
      nq    = Basis%nq   

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
           call d0d1d2poly_Hermite_exp(Basis%x(iq),ib - 1, Basis%d0gb(iq, ib),&
            Basis%d1gb(iq, ib, 1),Basis%d2gb(iq, ib, 1, 1), .true.) 
         END DO

      END DO
      
   END SUBROUTINE Construct_Basis_Ho


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
         !write(115,*) ib,jb,S(ib,ib),S(ib,jb)
           end do
          end do
          !IF (nderiv > -1) CALL   Write_VecMat(S,out_unit,5, info='S')
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
            S = matmul(conjg(Basis%d0bgw), Basis%d1gb(:, :, 1))
            !call   Write_VecMat(S,out_unit,5, info='<d0b|d1b>',Rformat='e13.4')
            !call   Write_VecMat(S,out_unit,5, info='<d0b|d1b>')
         END IF

         IF (nderiv > 1) THEN
            write(out_unit,*)
            S = matmul(conjg(Basis%d0bgw), Basis%d2gb(:, :, 1, 1))
             !call   Write_VecMat(S,out_unit,5, info='<d0b|d2b>',Rformat='e13.4')
             !call   Write_VecMat(S,out_unit,5, info='<d0b|d1b>')
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

   SUBROUTINE Hagedorn_construction(Basis)
       USE QDUtil_m
       TYPE(Basis_t), intent(inout)  :: Basis
        real(kind=Rkind):: Q0,SQ0

        SQ0 = Basis%scaleQ
        Q0 = Basis%Q0

        call Scale_Basis(Basis,Q0,SQ0)
        call Calc_Basis_dPtQtAt(Basis)
        Call Calc_tranpose_d0gb(Basis)
        Call Calc_dngg_grid(Basis)
        call CheckOrtho_Basis(Basis, nderiv=2)

   END SUBROUTINE

    SUBROUTINE Complet_scaling_Basis(Basis)
      USE QDUtil_m
      TYPE(Basis_t), intent(inout)  :: Basis
   
       
       call Scale_Basis(Basis,Basis%Q0,Basis%SCALEQ)
       call Complete_Hagedorn_none_variationnal_Basis(Basis)
       Call Calc_tranpose_d0gb(Basis)
       Call Calc_dngg_grid_0(Basis)
       call CheckOrtho_Basis(Basis, nderiv=2)
  END SUBROUTINE



   SUBROUTINE Calc_tranpose_d0gb(Basis)
      USE QDUtil_m
      TYPE(Basis_t), intent(inout)     :: Basis
      INTEGER                          :: ib,nb

      nb= Basis%nb

      if (allocated(Basis%d0bgw)) deallocate (Basis%d0bgw)
      if (Basis%Basis_name == 'el') then
         return
      end if

      Basis%d0bgw = transpose(Basis%d0gb)
      DO ib = 1, nb
         Basis%d0bgw(ib, :) = Basis%d0bgw(ib, :)*Basis%w(:)
      END DO

   END SUBROUTINE


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


   SUBROUTINE Calc_dngg_grid(Basis)
      USE QDUtil_m
      TYPE(Basis_t), intent(inout)    :: Basis
      integer                         :: ib,iq,nb,nq
      !logical, parameter             :: debug = .true.
      logical, parameter              :: debug = .false.

      nb = Basis%nb
      nq = Basis%nq

      IF (debug) THEN
         write (out_unit, *) 'BEGINNING Calc_dngg_grid'
         CALL Write_Basis(Basis)
         flush (out_unit)
      END IF

      if (allocated(Basis%d1gg)) deallocate (Basis%d1gg)
      if (allocated(Basis%d2gg)) deallocate (Basis%d2gg)
      if (allocated(Basis%dagg)) deallocate (Basis%dagg)
      if (allocated(Basis%dq0gg)) deallocate (Basis%dq0gg)
      if (allocated(Basis%dp0gg)) deallocate (Basis%dp0gg)

      allocate (Basis%d1gg(nq, nq, 1))
      allocate (Basis%d2gg(nq,nq, 1, 1))
      allocate (Basis%dagg(nq, nq))
      allocate (Basis%dq0gg(nq, nq))
      allocate (Basis%dp0gg(nq, nq))

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
         call   Write_VecMat(Basis%d1gg(:, :, 1), out_unit, 5,  info='d1gg')
         write (out_unit, *)
         call    Write_VecMat(Basis%d2gg(:, :, 1, 1), out_unit, 5,  info='d2gg')
         write (out_unit, *)
         call    Write_VecMat(Basis%dp0gg(:, :), out_unit, 5,  info='dp0gg')
         write (out_unit, *)
         call    Write_VecMat(Basis%dq0gg(:, :), out_unit, 5,  info='dq0gg')
         write (out_unit, *)
         call   Write_VecMat(Basis%dagg(:, :), out_unit, 5,  info='dagg')
         write (out_unit, *) 'END Calc_dngg_grid'
         flush (out_unit)
      END IF

   END SUBROUTINE 

   SUBROUTINE Calc_dngg_grid_0(Basis)
   USE QDUtil_m
   TYPE(Basis_t), intent(inout)    :: Basis
   integer                         :: ib,iq,nb,nq
   !logical, parameter             :: debug = .true.
   logical, parameter              :: debug = .false.
   nb = Basis%nb
   nq = Basis%nq
   IF (debug) THEN
      write (out_unit, *) 'BEGINNING Calc_dngg_grid'
      CALL Write_Basis(Basis)
      flush (out_unit)
   END IF
   if (allocated(Basis%d1gg)) deallocate (Basis%d1gg)
   if (allocated(Basis%d2gg)) deallocate (Basis%d2gg)
   allocate (Basis%d1gg(nq, nq, 1))
   allocate (Basis%d2gg(nq,nq, 1, 1))
   IF (debug) THEN
      CALL   Write_VecMat(Basis%d0bgw(:, :), out_unit, 5,  info='d0bgw')
      write (out_unit, *)
   END IF
   if (Basis%Basis_name == 'el') THEN
      RETURN
   end if
   Basis%d1gg(:, :, 1) = matmul(Basis%d1gb(:, :, 1), conjg(Basis%d0bgw))
   Basis%d2gg(:, :, 1, 1) = matmul(Basis%d2gb(:, :, 1, 1), conjg(Basis%d0bgw))
   IF (debug) THEN
      call   Write_VecMat(Basis%d1gg(:, :, 1), out_unit, 5,  info='d1gg')
      write (out_unit, *)
      call    Write_VecMat(Basis%d2gg(:, :, 1, 1), out_unit, 5,  info='d2gg')
      write (out_unit, *)
      write (out_unit, *) 'END Calc_dngg_grid'
      flush (out_unit)
   END IF
END SUBROUTINE 

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
      complex(kind=Rkind), intent(inout)          :: GB(:, :, :)
      complex(kind=Rkind), intent(in)             :: BB(:, :, :)
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
      TYPE(Basis_t), intent(in), target              :: Basis
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

   SUBROUTINE Calc_Q_grid(Q, Basis)

      implicit none
      TYPE(Basis_t), intent(in)                                    :: Basis
      integer, ALLOCATABLE                                         :: Tab_iq(:)
      integer                                                      :: inb, ndim, iq
      real(Kind=Rkind), intent(inout), allocatable                 :: Q(:, :)
      logical                                                      :: Endloop
      ndim = SIZE(Basis%tab_basis) - 1
       allocate (Q(Basis%nq, Ndim),Tab_iq(ndim))
      Call Init_tab_ind(Tab_iq, Basis%NDindexq)
      Iq = 0
      DO
         Iq = Iq + 1
         CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop)
         IF (Endloop) exit
         do inb = 1, Ndim
             Q(iq, inb) = Basis%tab_basis(inb)%x(Tab_iq(inb))
         end do
         !    print*,iq,Q(iq,:)
      END DO
      deallocate(Tab_iq)
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



SUBROUTINE Construct_Hagedorn_Variational_Basis(Basis,Qt,SQt,At,Pt)
  USE QDUtil_m
  IMPLICIT NONE
  TYPE(Basis_t),intent(inout)                     :: Basis
  real(kind=Rkind),intent(in)                     :: Qt(:), SQt(:),Pt(:)
  complex(kind=Rkind), intent(in)                 :: At(:)
  integer                                         :: ndim,ib

  ndim = size(Basis%tab_basis)-1
  call Change_Basis_Parameters(Basis,Qt,SQt,At,Pt)

  DO ib = 1,ndim
   if (Basis%tab_basis(ib)%Basis_name == 'herm' .or.Basis%tab_basis(ib)%Basis_name == 'ho') then
   call Construct_Basis_Ho(Basis%tab_basis(ib))
   call Scale_Basis(Basis%tab_basis(ib),Basis%tab_basis(ib)%Q0,Basis%tab_basis(ib)%SCALEQ)
   call Complete_Hagedorn_none_variationnal_Basis(Basis%tab_basis(ib))
   call Calc_tranpose_d0gb(Basis%tab_basis(ib))
   call Calc_Basis_dPtQtAt(Basis%tab_basis(ib))
   call Calc_dngg_grid(Basis%tab_basis(ib))
   call CheckOrtho_Basis(Basis%tab_basis(ib), nderiv=2)
   End If  
  End Do

End SUBROUTINE




SUBROUTINE Calc_d0d1d2W_temp(Q,Qt,SQt,Bt,Pt,d0W,d1W,d2W,nderiv)
  complex(kind=Rkind) , intent(inout)     :: d0W,d1W,d2W
  real(kind=Rkind)  , intent(in)          :: Q,Qt,SQt,Bt,Pt
  logical,intent(in)                      :: nderiv
  complex(kind=Rkind)                     :: cst1,cst2,cst
  real(kind=Rkind)                        ::DQ

    DQ = Q-Qt
    cst = -HALF*EYE*(Bt/SQt*SQt)*Q*Q+EYE*(Pt/SQt)*Q
    cst1 = -EYE*(Bt/SQt*SQt)*Q+EYE*(Pt/SQt)
    cst2 = -EYE*(Bt/SQt*SQt)+cst1*cst1
  If(nderiv) then
   d0W =exp(cst)
   d1W = cst1*d0W
   d2W =cst2*d0W
  Else
    d0W =exp(cst)
    d1W = CZERO
    d2W = CZERO
  End If
  !print*,'d0w, d1w,d2w',d0W,d1W,d2W
End SUBROUTINE


SUBROUTINE Calc_d0d1d2W(Q,Qt,SQt,Bt,Pt,d0W,d1W,d2W,nderiv)
  complex(kind=Rkind) , intent(inout)     :: d0W,d1W,d2W
  real(kind=Rkind)  , intent(in)          :: Q,Qt,SQt,Bt,Pt
  logical,intent(in)                      :: nderiv
  complex(kind=Rkind)                     :: cst1,cst2,cst
  real(kind=Rkind)                        ::DQ
   DQ = Q-Qt
  If(nderiv) then
   cst = -EYE*HALF*Bt*DQ*DQ+EYE*Pt*DQ
   cst1 = -EYE*Bt*DQ+EYE*Pt
   cst2 = -EYE*Bt+cst1*cst1
   d0W =exp(cst)
   d1W = cst1*d0W
   d2W =cst2*d0W
  Else
    cst = -EYE*HALF*Bt*DQ*DQ+EYE*Pt*DQ
    d0W =exp(cst)
    d1W = CZERO
    d2W = CZERO
  End If
  !print*,'d0w, d1w,d2w',d0W,d1W,d2W
End SUBROUTINE

SUBROUTINE Complete_Hagedorn_none_variationnal_Basis(Basis)
   TYPE(Basis_t), intent(inout)            :: Basis
   complex(kind=Rkind)                     :: d0W,d1W,d2W
   real(kind=Rkind)                        :: Qt,SQt,Bt,Pt
   integer                                 :: nb,nq,ib,iq
   nb  = Basis%nb
   nq  = Basis%nq
   Qt  = Basis%Q0
   SQt = Basis%SCALEQ
   Pt  = Basis%Imp_k
   Bt  = aimag(Basis%alpha)
  DO ib = 1,nb
     DO iq = 1,nq
         call Calc_d0d1d2W(Basis%x(iq),Qt,SQt,Bt,Pt,d0W,d1W,d2W,.true.)
         Basis%d2gb(iq, ib, 1, 1) =  Basis%d2gb(iq, ib, 1, 1)*d0W +TWO*Basis%d1gb(iq,ib,1)*d1W+Basis%d0gb(iq,ib)*d2W
         Basis%d1gb(iq,ib,1) = Basis%d1gb(iq,ib,1)*d0W+Basis%d0gb(iq,ib)*d1W
         Basis%d0gb(iq,ib) = Basis%d0gb(iq,ib)*d0W
     END DO
 END DO    
End SUBROUTINE




SUBROUTINE Complete_Hagedorn_none_variationnal_Basis_temp(Basis)
   TYPE(Basis_t), intent(inout)            :: Basis
   complex(kind=Rkind)                     :: d0W,d1W,d2W
   real(kind=Rkind)                        :: Qt,SQt,Bt,Pt
   integer                                 :: nb,nq,ib,iq
   nb  = Basis%nb
   nq  = Basis%nq
   Qt  = Basis%Q0
   SQt = Basis%SCALEQ
   Pt  = Basis%Imp_k
   Bt  = aimag(Basis%alpha)
  DO ib = 1,nb
     DO iq = 1,nq
         call Calc_d0d1d2W_temp(Basis%x(iq),Qt,SQt,Bt,Pt,d0W,d1W,d2W,.true.)
         Basis%d2gb(iq, ib, 1, 1) =  Basis%d2gb(iq, ib, 1, 1)*d0W +TWO*Basis%d1gb(iq,ib,1)*d1W+Basis%d0gb(iq,ib)*d2W
         Basis%d1gb(iq,ib,1) = Basis%d1gb(iq,ib,1)*d0W+Basis%d0gb(iq,ib)*d1W
         Basis%d0gb(iq,ib) = Basis%d0gb(iq,ib)*d0W
     END DO
 END DO    
End SUBROUTINE

SUBROUTINE Change_Basis_Parameters(Basis,Qt,SQt,At,Pt)

  TYPE(Basis_t),intent(inout)                      :: Basis
  real(kind=Rkind),intent(in)                      :: Qt(:), SQt(:),Pt(:)
  complex(kind=Rkind), intent(in)                  :: At(:)
  integer                                          :: ndim,ib

ndim = size(Basis%tab_basis)-1

Do ib = 1,ndim
    Basis%tab_basis(ib)%Q0  = Qt(ib)
    Basis%tab_basis(ib)%SCALEQ  = SQt(ib)
    Basis%tab_basis(ib)%Imp_k  = Pt(ib)
    Basis%tab_basis(ib)%alpha  = At(ib)
End DO

End SUBROUTINE


SUBROUTINE Get_Basis_Parameters(Basis,Qt,SQt,At,Pt)

  TYPE(Basis_t),intent(in)                      :: Basis
  real(kind=Rkind),intent(inout)                :: Qt(:), SQt(:),Pt(:)
  complex(kind=Rkind), intent(inout)            :: At(:)
  integer                                       :: ndim,ib

ndim = size(Basis%tab_basis)-1

Do ib = 1,ndim
  Qt(ib) = Basis%tab_basis(ib)%Q0
  SQt(ib) = Basis%tab_basis(ib)%SCALEQ
  Pt(ib) = Basis%tab_basis(ib)%Imp_k
  At(ib) = Basis%tab_basis(ib)%alpha
End DO

End SUBROUTINE

END MODULE Basis_m
