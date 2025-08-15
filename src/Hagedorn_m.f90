module Hagedorn_m
   USE QDUtil_m
   Use Basis_m
   Use Psi_m
   Use Op_m
   USE Ana_psi_m
   USE polyortho_m
   USE sub_propa_m
   implicit none
   PRIVATE
     !> Public interfaces for Hagedorn basis construction and projection operations
   PUBLIC :: Projection_temp                  ! Temporary projection routine
   PUBLIC :: test_psi_temp                    ! Temporary test of psi
   PUBLIC :: Calc_Basis_parameters, Calc_Basis_parameters_temp ! Compute basis parameters from a psi structure
   PUBLIC :: march_Hagedorn                   ! Marching algorithm using Hagedorn method
   PUBLIC :: march_Global                     ! Global propagation step
   
     !> Hagedorn propagation and testing routines
   PUBLIC :: Hagedorn_inv                     ! Hagedorn propagation with inversion
   PUBLIC :: Hagedorn_temp                    ! Temporary Hagedorn routine
   PUBLIC :: Test_Calc_S                      ! Test routine for computing matrix S
   PUBLIC :: Test_projection_hagedorn         ! Test of projection method using Hagedorn basis
   
     !> Utility routines for basis construction and manipulation
   PUBLIC :: build_identity_matrix            ! Build identity matrix for trivial bases
   PUBLIC :: build_overlap_matrices           ! Compute overlap matrices
   
     !> Hagedorn basis construction routines
   PUBLIC :: Construct_Hagedorn_none_Variational_Basis       ! Construct non-variational Hagedorn basis
   PUBLIC :: Construct_Hagedorn_none_Variational_Basis_temp  ! Temporary version of non-variational Hagedorn basis
   
CONTAINS

SUBROUTINE Construct_Hagedorn_none_Variational_Basis(Basis, Qt, SQt, At, Pt)
  USE QDUtil_m
  IMPLICIT NONE

  ! === Input/Output Arguments ===
  TYPE(Basis_t),            INTENT(INOUT)  :: Basis      ! Basis structure to update
  REAL(kind=Rkind),         INTENT(IN)     :: Qt(:)      ! <Q>
  REAL(kind=Rkind),         INTENT(IN)     :: SQt(:)     ! sqrt(<Q²> - <Q>²)
  COMPLEX(kind=Rkind),      INTENT(IN)     :: At(:)      ! Complex Gaussian widths
  REAL(kind=Rkind),         INTENT(IN)     :: Pt(:)      ! <P>

  ! === Local variables ===
  REAL(kind=Rkind), ALLOCATABLE :: Q0(:), SQ0(:), P0(:)
  COMPLEX(kind=Rkind), ALLOCATABLE :: A0(:)
  INTEGER :: ndim, ib, nb, nq, i

  ! === Determine number of dimensions (excluding electronic)
  ndim = SIZE(Basis%tab_basis) - 1

  ! === Allocate arrays to hold original basis parameters
  ALLOCATE(Q0(ndim), SQ0(ndim), P0(ndim), A0(ndim))

  ! === Store the current basis parameters (before change)
  CALL Get_Basis_Parameters(Basis, Q0, SQ0, A0, P0)

  ! === Update basis parameters with new Qt, SQt, At, Pt
  CALL Change_Basis_Parameters(Basis, Qt, SQt, At, Pt)

  ! === Rebuild primitive basis structure
  CALL construct_primitive_basis_temp(Basis)

  ! === For each dimension, calculate transformation matrix S
  DO ib = 1, ndim
    nb = Basis%tab_basis(ib)%nb
    nq = Basis%tab_basis(ib)%nq

    IF (Basis%tab_basis(ib)%Basis_name == 'herm' .OR. &
        Basis%tab_basis(ib)%Basis_name == 'ho') THEN

      ! Compute transformation matrix between old and new basis
      CALL Calc_S(Basis%tab_basis(ib)%S, nb, nq, &
                  Q0(ib), SQ0(ib), A0(ib), P0(ib), &
                  Qt(ib), SQt(ib), At(ib), Pt(ib))

    ELSE
      ! Set identity matrix for non-Hagedorn or non-HO basis types
      call build_identity_matrix(Basis%tab_basis(ib)%S, nb)
    END IF
  END DO

  ! === Clean up: deallocate local arrays
  DEALLOCATE(Q0, SQ0, P0, A0)

END SUBROUTINE Construct_Hagedorn_none_Variational_Basis

subroutine build_identity_matrix(Id, n)
  implicit none
  integer, intent(in)                            :: n
  complex(kind=Rkind), allocatable, intent(out)  :: Id(:,:)
  integer                                        :: ib

  ! Allocate the identity matrix
  allocate(Id(n,n))

  ! Initialize all elements to zero
  Id(:,:) = CZERO

  ! Set diagonal elements to 1
  do ib = 1, n
     Id(ib,ib) = CONE
  end do

end subroutine build_identity_matrix

SUBROUTINE Construct_Hagedorn_none_Variational_Basis_temp(Basis, Qt, SQt, At, Pt)
  USE QDUtil_m
  IMPLICIT NONE

  ! === Interface ===
  TYPE(Basis_t),            INTENT(INOUT)  :: Basis         ! Basis to be modified
  REAL(kind=Rkind),         INTENT(IN)     :: Qt(:)         ! <Q>
  REAL(kind=Rkind),         INTENT(IN)     :: SQt(:)        ! sqrt(<Q²> - <Q>²)
  COMPLEX(kind=Rkind),      INTENT(IN)     :: At(:)         ! Complex width A
  REAL(kind=Rkind),         INTENT(IN)     :: Pt(:)         ! <P>

  ! === Local variables ===
  REAL(kind=Rkind), ALLOCATABLE :: Q0(:), SQ0(:), P0(:)
  COMPLEX(kind=Rkind), ALLOCATABLE :: A0(:)
  INTEGER :: ndim, ib, nb, nq

  ! === Determine dimensionality
  ndim = SIZE(Basis%tab_basis) - 1

  ! === Allocate arrays to store current basis parameters
  ALLOCATE(Q0(ndim), SQ0(ndim), P0(ndim), A0(ndim))

  ! === Get original basis parameters
  CALL Get_Basis_Parameters(Basis, Q0, SQ0, A0, P0)

  ! === Change the parameters of the basis to new values
  CALL Change_Basis_Parameters(Basis, Qt, SQt, At, Pt)

  ! === Rebuild the primitive basis
  CALL construct_primitive_basis_temp(Basis)

  ! === Loop over all dimensions to compute S matrix
  DO ib = 1, ndim
    nb = Basis%tab_basis(ib)%nb
    nq = Basis%tab_basis(ib)%nq

    IF (Basis%tab_basis(ib)%Basis_name == 'herm' .OR. &
        Basis%tab_basis(ib)%Basis_name == 'ho') THEN

      CALL Calc_S(Basis%tab_basis(ib)%S, nb, nq, &
                  Q0(ib), SQ0(ib), A0(ib), P0(ib), &
                  Qt(ib), SQt(ib), At(ib), Pt(ib))
    ELSE
      CALL build_identity_matrix(Basis%tab_basis(ib)%S, nb)
    END IF
  END DO

  ! === Cleanup: deallocate local arrays
  DEALLOCATE(Q0, SQ0, P0, A0)

END SUBROUTINE Construct_Hagedorn_none_Variational_Basis_temp
SUBROUTINE Calc_Basis_parameters(psi, Qt, SQt, At, Pt)
   USE psi_m
   IMPLICIT NONE
  
   ! === Interfaces ===
   TYPE(psi_t),         INTENT(IN)     :: psi         ! Input wavefunction
   REAL(kind=Rkind),    INTENT(INOUT)  :: Qt(:)       ! Expectation value of position <Q>
   REAL(kind=Rkind),    INTENT(INOUT)  :: SQt(:)      ! Standard deviation of position
   COMPLEX(kind=Rkind), INTENT(INOUT)  :: At(:)       ! Complex width parameter A
   REAL(kind=Rkind),    INTENT(INOUT)  :: Pt(:)       ! Expectation value of momentum <P>
  
   ! === Local variables ===
   INTEGER, ALLOCATABLE                :: Tab_Iq(:, :)  ! Index table for integration
  
   ! === Build index table for integration (e.g. <ψ|Q|ψ>, <ψ|Q²|ψ>, etc.)
   CALL Calc_tab_Iq0(Tab_Iq, psi%Basis)
  
   ! === Compute expectation values of Q and its standard deviation
   CALL Calc_AVQ_SQ_nD(psi, Qt, SQt, Tab_Iq)
  
   ! === Compute expectation value of the momentum operator
   CALL Calc_Av_imp_k_nD(psi, Pt)
  
   ! === Compute complex width parameter A (related to beta/squeezing)
   CALL Calc_Avg_A_nD(psi, At)
  
   ! === Deallocate index table
   DEALLOCATE(Tab_Iq)
  
 END SUBROUTINE Calc_Basis_parameters
  
 SUBROUTINE Calc_Basis_parameters_temp(psi, Qt, SQt, At, Pt, propa)
   USE psi_m
   IMPLICIT NONE
  
   ! === Interfaces ===
   TYPE(psi_t),       INTENT(IN)     :: psi        ! Wavefunction input
   TYPE(propa_t),     INTENT(IN)     :: propa      ! Propagation settings
   REAL(kind=Rkind),  INTENT(INOUT)  :: Qt(:)      ! Position expectation value <Q>
   REAL(kind=Rkind),  INTENT(INOUT)  :: SQt(:)     ! Width sqrt(<Q^2> - <Q>^2)
   REAL(kind=Rkind),  INTENT(INOUT)  :: Pt(:)      ! Momentum expectation <P>
   COMPLEX(kind=Rkind), INTENT(INOUT):: At(:)      ! "A" parameter (related to Gaussian width)
  
   ! === Local ===
   INTEGER, ALLOCATABLE              :: Tab_Iq(:, :)
  
   ! === Compute the index table.
   CALL Calc_tab_Iq0(Tab_Iq, psi%Basis)
  
   ! === Default initialization
   Pt(:) = ZERO
  
   ! === Compute <Q> and Sqrt(<Q^2> - <Q>^2)
   CALL Calc_AVQ_SQ_nD(psi, Qt, SQt, Tab_Iq)
  
   ! === Set At = SQt^2 by default
   At(:) = SQt(:) * SQt(:)
  
   ! === Optionally compute momentum expectation value
   IF (propa%P) THEN
      CALL Calc_Av_imp_k_nD(psi, Pt)
   END IF
  
   ! === Optionally compute complex A parameters
   IF (propa%Beta) THEN
      CALL Calc_Avg_A_nD(psi, At)
   END IF
  
   ! === Deallocate index table
   DEALLOCATE(Tab_Iq)
  
 END SUBROUTINE Calc_Basis_parameters_temp


SUBROUTINE Calc_S(S,nb,nq,Q0,SQ0,A0,P0,Qt,SQt,At,Pt)

  IMPLICIT NONE
  ! Input/Output
  complex(kind=Rkind), intent(inout) :: S(:,:)

  ! Input parameters
  complex(kind=Rkind), intent(in)    :: A0, At
  real(kind=Rkind), intent(in)       :: Q0, SQ0, P0
  real(kind=Rkind), intent(in)       :: Qt, SQt, Pt
  integer, intent(in)                :: nb, nq

  ! Local variables
  integer :: ib, iq, jb
  real(kind=Rkind) :: SQeq, Qeq
  real(kind=Rkind), allocatable :: Q(:), W(:)
  real(kind=Rkind) :: B0, Bt
  complex(kind=Rkind), allocatable :: d0gb0(:,:),d0gb1(:,:), d0bgw(:,:)

  !--------------------------------------------------------------
  ! Compute effective scaling factor SQeq and center Qeq
  !--------------------------------------------------------------
  SQeq = sqrt(SQ0*SQ0 + SQt*SQt) / sqrt(TWO)
  Qeq  = (SQ0*SQ0*Q0 + SQt*SQt*Qt) / (SQ0*SQ0 + SQt*SQt)

  !--------------------------------------------------------------
  ! Extract imaginary parts of At and A0 (b_t and b_0)
  !--------------------------------------------------------------
  Bt = AIMAG(At)
  B0 = AIMAG(A0)

  !--------------------------------------------------------------
  ! Allocate grid points and weights for Hermite quadrature
  !--------------------------------------------------------------
  allocate(Q(nq), W(nq))
  call hercom(nq, Q(:), W(:))

  ! Scale quadrature weights and shift grid points
  W(:) = W(:) / SQeq
  Q(:) = Qeq + Q(:) / SQeq

  !--------------------------------------------------------------
  ! Allocate basis function arrays
  !--------------------------------------------------------------
  allocate(d0gb0(nq, nb))
  allocate(d0gb1(nq, nb))
  allocate(d0bgw(nb, nq))

  !--------------------------------------------------------------
  ! Debug prints for checking parameters
  !--------------------------------------------------------------
  !print*, 'Q0, A0, P0 =', Q0, A0, P0
  !print*, 'Qt, At, Pt =', Qt, At, Pt
  !print*, 'SQeq, Qeq =', SQeq, Qeq

  !--------------------------------------------------------------
  ! Compute Hagedorn basis functions at grid points Q(:)
  ! for both initial (Q0, A0, P0) and target (Qt, At, Pt)
  !--------------------------------------------------------------
  DO iq = 1, nq
      DO ib = 1, nb
          call d0poly_Hermite_exp_cplx(Q(iq), Q0, A0, P0, ib-1, d0gb0(iq, ib))
          call d0poly_Hermite_exp_cplx(Q(iq), Qt, At, Pt, ib-1, d0gb1(iq, ib))
      END DO
  END DO

  d0bgw = transpose(d0gb1)
   DO ib = 1, nb
      d0bgw(ib, :) = d0bgw(ib, :)*W(:)
   END DO

  !--------------------------------------------------------------
  ! Compute overlap matrix S = <d0bgw|d0gb>
  !--------------------------------------------------------------
  S = MATMUL(conjg(d0bgw), d0gb0)

  !--------------------------------------------------------------
  ! Print overlap matrix S for verification
  !--------------------------------------------------------------
  !call Write_VecMat(S, out_unit, 5, info='S')

  !--------------------------------------------------------------
  ! Deallocate arrays
  !--------------------------------------------------------------
  deallocate(d0gb0, d0gb1, d0bgw, W, Q)

End SUBROUTINE

SUBROUTINE Test_Calc_S()
  USE QDUtil_m
  IMPLICIT NONE

  integer, parameter :: nb = 3, nq = 25
  complex(kind=Rkind) :: S(nb, nb)
  real(kind=Rkind) :: Q0, SQ0, P0, Qt, SQt, Pt
  complex(kind=Rkind) :: A0, At

  !--------------------------------------------------------------
  ! Initialize Parameters for the Test
  !--------------------------------------------------------------

  A0  = CMPLX(1.0_Rkind, 0.1_Rkind, kind=Rkind)
  Q0  = 0.5_Rkind
  P0  = 0.25_Rkind
  SQ0 = 1.0_Rkind

  At  = CMPLX(1.0_Rkind, 0.2_Rkind, kind=Rkind)
  Qt  = 0.49_Rkind
  SQt = 1.0_Rkind
  Pt  = 0.24_Rkind
  print*, 'Q0, SQ0, A0, P0 =', Q0, SQ0, A0, P0
  print*, 'Qt, SQt, At, Pt =', Qt, SQt, At, Pt

  !--------------------------------------------------------------
  ! Call Calc_S to compute the overlap matrix S
  !--------------------------------------------------------------
  call Calc_S(S, nb, nq, Q0, SQ0, A0, P0, Qt, SQt, At, Pt)
    call  Write_VecMat(S, out_unit, nb,  info='S')
  !--------------------------------------------------------------
  ! End of test subroutine
  !--------------------------------------------------------------
  print*, 'Test_Calc_S finished.'

END SUBROUTINE Test_Calc_S


subroutine build_overlap_matrices(B1, B2)
  implicit none

  ! Input/output
  TYPE(Basis_t), intent(inout)         :: B1, B2

  ! Local variables
  integer                              :: ib, nb, nq, ndim
  character(len=10)                    :: bname
  real(kind=Rkind)                     :: Q0, SQ0, P0, Qt, SQt, Pt
  complex(kind=Rkind)                  :: A0, At

  ! Determine number of dimensions (suppose même taille pour B1 et B2)
  ndim = size(B1%tab_basis)

  ! Loop over each coordinate
  do ib = 1, ndim
     bname = B1%tab_basis(ib)%Basis_name

     ! Extract parameters of basis B1
     nb   = B1%tab_basis(ib)%nb
     nq   = B1%tab_basis(ib)%nq
     Q0   = B1%tab_basis(ib)%Q0
     SQ0  = B1%tab_basis(ib)%scaleQ
     A0   = B1%tab_basis(ib)%alpha
     P0   = B1%tab_basis(ib)%Imp_k

     ! Extract parameters of basis B2
     Qt   = B2%tab_basis(ib)%Q0
     SQt  = B2%tab_basis(ib)%scaleQ
     At   = B2%tab_basis(ib)%alpha
     Pt   = B2%tab_basis(ib)%Imp_k
      
     print*, 'Q0,Qt', Q0,Qt
     print*, 'SQ0,SQt', SQ0,SQt
     print*, 'A0,At', A0,At
     print*, 'P0,Pt', P0,Pt
     print*, 'nb,nq', nb,nq 
     ! Build S according to the basis type
     if (trim(bname) == 'hag' .or. trim(bname) == 'ho') then
        call Calc_S(B1%tab_basis(ib)%S, nb, nq, Q0, SQ0, A0, P0, Qt, SQt, At, Pt)
     else
        call build_identity_matrix(B1%tab_basis(ib)%S, nb)
     end if
     call Write_VecMat(B1%tab_basis(ib)%S, out_unit, nb, info='S')
  end do

end subroutine build_overlap_matrices



  SUBROUTINE projection_1D(BBB2, BBB1, Basis)
     USE QDUtil_m
     TYPE(Basis_t), intent(in), target              :: Basis
     complex(kind=Rkind), intent(inout)             :: BBB2(:, :, :)
     complex(kind=Rkind), intent(in)                :: BBB1(:, :, :)
     logical, parameter                             :: debug = .true.
     Integer                                        :: i1, i3
   
     IF (debug) THEN
        flush (out_unit)
     END IF
     BBB2 = CZERO
     DO i3 = 1, ubound(BBB1, dim=3)
     DO i1 = 1, ubound(BBB1, dim=1)
        BBB2(i1, :, i3) =  matmul(Basis%S,BBB1(i1, :, i3))
     END DO
     END DO
     IF (debug) THEN
        flush (out_unit)
     END IF
 END SUBROUTINE

  SUBROUTINE projection_1D_temp(BBB2, BBB1, S)
    USE QDUtil_m
    complex(kind=Rkind), intent(inout)             :: BBB2(:, :, :)
    complex(kind=Rkind), intent(in)                :: BBB1(:, :, :),S(:,:)
    logical, parameter                             :: debug = .true.
    Integer                                        :: i1, i3
  
    IF (debug) THEN
       flush (out_unit)
    END IF
    BBB2 = CZERO
    DO i3 = 1, ubound(BBB1, dim=3)
    DO i1 = 1, ubound(BBB1, dim=1)
       BBB2(i1, :, i3) =  matmul(S,BBB1(i1, :, i3))
    END DO
    END DO
    IF (debug) THEN
       flush (out_unit)
    END IF
END SUBROUTINE


   SUBROUTINE Projection_temp(psi, psi_dt)
   TYPE(psi_t), intent(in), target                  :: psi_dt
   TYPE(psi_t), intent(inout), target               :: psi
   complex(kind= Rkind), pointer                    :: BBB1(:, :, :), BBB2(:, :, :)
   complex(kind= Rkind), allocatable, target        :: B1(:), B2(:)
   real(kind= Rkind)                                :: Norm0,norm,E,E0
   logical, parameter                               :: debug = .true.
   integer                                          :: ib, ndim
   Integer, allocatable                             :: Ib1(:), Ib2(:), Ib3(:)

   call Calc_index(Ib1=Ib1, Ib2=Ib2, Ib3=Ib3, Basis=psi%Basis)
   ndim = size(psi%Basis%tab_basis) - 1
   call Calc_Norm_OF_psi(psi_dt,Norm0)
   call Calc_average_energy(psi_dt, E0)
   psi%CVec(:) = CZERO
   allocate (B1(Ib1(1)*Ib2(1)*Ib3(1)))
   allocate (B2(Ib1(1)*Ib2(1)*Ib3(1)))

    B1(:) = CZERO
    B2(:) = CZERO
    B1(:) = psi_dt%CVec(:)

    DO ib = 1,ndim
       BBB1(1:Ib1(ib), 1:Ib2(ib), 1:Ib3(ib)) => B1
       BBB2(1:Ib1(ib), 1:Ib2(ib), 1:Ib3(ib)) => B2
       call projection_1D(BBB2, BBB1, psi_dt%Basis%tab_basis(ib))
       B1(:) = B2(:)
    END DO
    psi%CVec(:) = B2(:)
    
   call  Calc_average_energy(psi, E)
     deallocate(B1)
     deallocate(B2)
    call Calc_Norm_OF_psi(psi,Norm)
   ! write(25,*) E0,E
   !  write(24,*) Norm0,Norm
    write (out_unit, *) 'Begin Hagedorn projection  E0,Norm0',E0,Norm0
    write (out_unit, *) 'END Hagedorn projection E,Norm ',E,Norm
END SUBROUTINE 

subroutine Test_projection_hagedorn(Basis,psi)
  implicit none

  TYPE(Basis_t), intent(in)                 :: Basis
  TYPE(psi_t), intent(in)                   :: psi

  ! Local variables
  TYPE(psi_t)                               :: psi2, psi1
  TYPE(Basis_t)                             :: B1, B2
  type(Op_t)                                :: H1,H2
  integer, allocatable                      :: Tab_iq1(:,:)
  integer, allocatable                      :: Tab_iq2(:,:)
  
  real(kind=Rkind)                          :: norm1  , norm2, E1, E2
  integer                                   :: ib, ndim, n_basis
  complex(kind=Rkind)                       :: alpha(2)
  real(kind=Rkind)                          :: q(2), p(2), sq(2)

  ! Initialization of parameters
  alpha = [1.0_Rkind + 0.1_Rkind*EYE, 1.0_Rkind + 0.1_Rkind*EYE]
  q     = [1.9_Rkind, 0.1_Rkind]
  p     = [0.1_Rkind, 0.1_Rkind]
  sq    = [1.0_Rkind, 1.0_Rkind]



  ! Copy Basis into B1 and B2
  call init_Basis1_TO_Basis2(B1, Basis)
  call init_Basis1_TO_Basis2(B2, Basis)

  ! Modify the parameters of B2
  B2%tab_basis(1)%alpha    = alpha(1)
  B2%tab_basis(2)%alpha    = alpha(2)
  B2%tab_basis(1)%q0       = q(1)
  B2%tab_basis(2)%q0       = q(2)
  B2%tab_basis(1)%Imp_k    = p(1)
  B2%tab_basis(2)%Imp_k    = p(2)
  B2%tab_basis(1)%scaleQ   = sq(1)
  B2%tab_basis(2)%scaleQ   = sq(2)

  ! Construct both primitive bases
  call construct_primitive_basis(B1)
  call construct_primitive_basis(B2)

  call Calc_tab_Iq0(Tab_Iq1,B1)
  call Calc_tab_Iq0(Tab_Iq2,B2)

  call Set_Op(H1, B1,Tab_Iq1)
  call Set_Op(H2, B2,Tab_Iq2)

  ! Initialize psi2 with B2 and psi1 with B1
  call init_psi(psi2, B2, cplx=.true., grid=.false.)
  call init_psi(psi1, B1, cplx=.true., grid=.false.)

  call build_overlap_matrices(B1,B2)

  ! Fill psi1 with CZERO except for the first coefficient
  psi1%CVec(:) = psi%CVec(:)


  ! Initialize psi2 to zero
  psi2%CVec(:) = CZERO

  ! Print the parameters of both bases for verification
  print *, 'B1'
  print *, "Basis is allocated", Basis_IS_allocated(B1)
  print *, "B1%tab_basis(1)%alpha(:)",    B1%tab_basis(1)%alpha
  print *, "B1%tab_basis(2)%alpha(:)",    B1%tab_basis(2)%alpha
  print *, "B1%tab_basis(1)%q0(:)",       B1%tab_basis(1)%q0
  print *, "B1%tab_basis(2)%q0(:)",       B1%tab_basis(2)%q0
  print *, "B1%tab_basis(1)%imp_k(:)",    B1%tab_basis(1)%Imp_k
  print *, "B1%tab_basis(2)%imp_k(:)",    B1%tab_basis(2)%Imp_k
  print *, "B1%tab_basis(1)%scaleQ(:)",  B1%tab_basis(1)%scaleQ
  print *, "B1%tab_basis(2)%scaleQ(:)",  B1%tab_basis(2)%scaleQ

  print *, 'B2'
  print *, "Basis is allocated", Basis_IS_allocated(B2)
  print *, "B2%tab_basis(1)%alpha(:)",    B2%tab_basis(1)%alpha
  print *, "B2%tab_basis(2)%alpha(:)",    B2%tab_basis(2)%alpha
  print *, "B2%tab_basis(1)%q0(:)",       B2%tab_basis(1)%q0
  print *, "B2%tab_basis(2)%q0(:)",       B2%tab_basis(2)%q0
  print *, "B2%tab_basis(1)%imp_k(:)",    B2%tab_basis(1)%Imp_k
  print *, "B2%tab_basis(2)%imp_k(:)",    B2%tab_basis(2)%Imp_k
  print *, "B2%tab_basis(1)%scaleQ(:)",  B2%tab_basis(1)%scaleQ
  print *, "B2%tab_basis(2)%scaleQ(:)",  B2%tab_basis(2)%scaleQ


   ! Compute the norm and energy of psi1
   call Calc_Av_E(E1, psi1, H1)
   call Calc_Norm_OF_psi(psi1, norm1)
   print *, 'Norm, E1', norm1, E1
  ! Display psi1 coefficients
  print *, '--------------------------------------------------------'
  do ib = 1, size(psi1%CVec)
     write(out_unit, *) ib, abs(psi1%CVec(ib))**2
  end do


  ! Project psi1 (constructed with B1) onto B2; result is stored in psi2
  call projection1(psi2, psi1)

  ! Display the coefficients of psi2
  do ib = 1, size(psi2%CVec)
     write(out_unit, *) ib, abs(psi2%CVec(ib))**2
  end do

  ! Compute the norm and energy of psi2
  call Calc_Av_E(E2, psi2, H2)
  call Calc_Norm_OF_psi(psi2, norm2)
  print *, 'Norm, E2', norm2, E2

  ! End of projection test
end subroutine test_projection_hagedorn

subroutine test_psi_temp(psi,propa)
   implicit none
   TYPE(psi_t), target, intent(inout)  :: psi
    type(propa_t),intent(in)           :: propa
    TYPE(Basis_t)                      ::  Basis1,Basis2
   TYPE(psi_t)                         :: psi2,psi1
   real(kind= Rkind)                   :: Norm0,norm,E,E0
   integer                             :: ib
  
   call init_Basis1_TO_Basis2(Basis1, psi%Basis)
   call construct_primitive_basis(Basis1)
   call init_Basis1_TO_Basis2(Basis2, psi%Basis)
   call construct_primitive_basis(Basis2)
   call init_psi(psi2, psi%Basis, cplx=.true., grid=.false.)
   call init_psi(psi1, Basis1, cplx=.true., grid=.false.)
   psi2%CVec(:) = CZERO
   print *, '--------------------------------------------------------'

       do ib = 1,size(psi%CVec)
         write(out_unit,*) ib,abs(psi%CVec(ib))**2
       end Do 

       call  Calc_average_energy(psi, E0) 
       call Hagedorn_temp(psi2, psi,propa)
       call  Calc_average_energy(psi2, E)


       do ib = 1,size(psi%CVec)
           write(out_unit,*) ib,abs(psi2%CVec(ib))**2
      end Do 
       write(out_unit,*)  'E0,E' ,E0,E 

       !call  Hagedorn_inv(psi1, psi2)

        do ib = 1,size(psi%CVec)
         write(out_unit,*) ib,abs(psi%CVec(ib))**2,abs(psi1%CVec(ib))**2,abs(psi%CVec(ib))**2-abs(psi1%CVec(ib))**2
       end Do 

    print *, '--------------------------------------------------------'

end subroutine 


  SUBROUTINE march_Hagedorn(psi, psi_dt, t, propa)
  USE psi_m
  type(propa_t),intent(in)                :: propa
  type(psi_t),  intent(inout)             :: psi
  type(psi_t),  intent(inout)             :: psi_dt
  real(kind=Rkind), intent(in)            :: t

   call  march(psi, psi_dt, t, propa)
   call Hagedorn_temp(psi, psi_dt,propa)

 END SUBROUTINE


   SUBROUTINE Hagedorn_temp(psi, psi_dt, propa)
    USE psi_m
    IMPLICIT NONE
  
    ! Declarations
    TYPE(psi_t),     INTENT(INOUT) :: psi         ! Wavefunction to project into new basis
    TYPE(psi_t),     INTENT(INOUT) :: psi_dt      ! Temporary wavefunction from time propagation
    TYPE(propa_t),   INTENT(IN)    :: propa       ! Propagation parameters
  
    COMPLEX(kind=Rkind), ALLOCATABLE :: At(:)     ! A(t) parameters
    REAL(kind=Rkind),    ALLOCATABLE :: Qt(:), SQt(:), Pt(:)  ! Q(t), sQt(t), and P(t)
    REAL(kind=Rkind)                  :: Norm0, norm
    INTEGER                           :: ndim
  
    ! Get number of degrees of freedom
    ndim = SIZE(psi%Basis%tab_basis) - 1
  
    ! Allocate arrays for Gaussian parameters
    ALLOCATE(Qt(ndim), SQt(ndim), Pt(ndim), At(ndim))
  
    ! Compute parameters from psi_dt
    CALL Calc_Basis_parameters_temp(psi_dt, Qt, SQt, At, Pt, propa)
  
    ! Reconstruct new non-variational Hagedorn basis using those parameters
    CALL Construct_Hagedorn_none_Variational_Basis_temp(psi_dt%Basis, Qt, SQt, At, Pt)
  
    ! Save norm of the unprojected psi_dt
    CALL Calc_Norm_OF_psi(psi_dt, Norm0)
  
    ! Project psi_dt into the new basis and update psi
    CALL Projection(psi, psi_dt)
  
    ! Optional renormalization if requested
    IF (propa%renorm) THEN
       CALL Calc_Norm_OF_psi(psi, norm)
       psi%CVec(:) = psi%CVec(:) / norm
    END IF
  
    ! Final norm check and output
    CALL Calc_Norm_OF_psi(psi, norm)
    WRITE(out_unit,*) ABS(Norm0 - norm), Norm0, norm
  
    ! Clean-up
    DEALLOCATE(Qt, SQt, At, Pt)
  
  END SUBROUTINE Hagedorn_temp



 SUBROUTINE Hagedorn_inv(psi0, psi_dt, renorm)
  USE psi_m
  IMPLICIT NONE

  TYPE(psi_t), INTENT(INOUT)             :: psi0      ! Wavefunction to be updated
  TYPE(psi_t), INTENT(IN)                :: psi_dt    ! Reference wavefunction
  LOGICAL, INTENT(IN)                    :: renorm    ! Whether to renormalize psi0

  ! Basis parameters
  COMPLEX(KIND=Rkind), ALLOCATABLE       :: A0(:)     ! Complex scaling parameter
  REAL(KIND=Rkind), ALLOCATABLE          :: Q0(:), SQ0(:), P0(:) ! Gaussian parameters

  ! Norm variables
  REAL(KIND=Rkind)                       :: Norm0, norm
  INTEGER                                :: ndim

  !-------------------------------- Initialization --------------------------------
  ndim = SIZE(psi0%Basis%tab_basis) - 1
  ALLOCATE(Q0(ndim), SQ0(ndim), P0(ndim), A0(ndim))

  ! Extract basis parameters from psi0
  CALL Get_Basis_Parameters(psi0%Basis, Q0, SQ0, A0, P0)

  ! Alternative: temporary calculation of parameters (commented)
  !CALL Calc_Basis_parameters_temp(psi0, Q0, SQ0, A0, P0, propa)

  ! Construct Hagedorn basis for psi_dt using parameters from psi0
  CALL Construct_Hagedorn_none_Variational_Basis_temp(psi_dt%Basis, Q0, SQ0, A0, P0)

  ! Compute norm of psi_dt before projection
  CALL Calc_Norm_OF_psi(psi_dt, Norm0)

  ! Project psi_dt onto psi0's basis
  CALL Projection(psi0, psi_dt)

  !-------------------------------- Optional Renormalization ----------------------
  IF (renorm) THEN
     CALL Calc_Norm_OF_psi(psi0, norm)
     psi0%CVec(:) = psi0%CVec(:) / norm
  END IF

  !-------------------------------- Final Norm Check ------------------------------
  CALL Calc_Norm_OF_psi(psi0, norm)
  WRITE(out_unit,*) ABS(Norm0 - norm), Norm0, norm

  !-------------------------------- Clean Up --------------------------------------
  DEALLOCATE(Q0, SQ0, A0, P0)

END SUBROUTINE

 SUBROUTINE march_Global(psi, psi_dt, t, propa, H)
  USE psi_m
  IMPLICIT NONE

  ! Input/output arguments
  type(propa_t), INTENT(IN)     :: propa     ! Propagation parameters
  type(psi_t),    INTENT(INOUT) :: psi       ! Current wavefunction (updated if needed)
  type(psi_t),    INTENT(INOUT) :: psi_dt    ! Temporary wavefunction (used for propagation)
  real(kind=Rkind), INTENT(IN)  :: t         ! Current time
  type(Op_t),     INTENT(INOUT) :: H         ! Hamiltonian operator (includes potential, kinetic, etc.)

  ! Local variables
  complex(kind=Rkind), ALLOCATABLE :: At(:)      ! Complex width parameters for each degree of freedom
  real(kind=Rkind),    ALLOCATABLE :: Qt(:)      ! Gaussian centers in position
  real(kind=Rkind),    ALLOCATABLE :: SQt(:)     ! Widths (real part of square root of At)
  real(kind=Rkind),    ALLOCATABLE :: Pt(:)      ! Gaussian centers in momentum
  integer, ALLOCATABLE             :: Tab_iq(:, :) ! Index mapping for operator construction

  real(kind=Rkind) :: Norm0, norm      ! Norms before and after propagation (for conservation check)
  real(kind=Rkind) :: E0, E            ! Energies before and after propagation
  integer          :: ndim             ! Number of degrees of freedom (basis size - 1)

  ! Determine number of degrees of freedom from the basis
  ndim = size(psi%Basis%tab_basis) - 1

  ! Allocate and initialize Gaussian parameters
  ALLOCATE(Qt(ndim), SQt(ndim), Pt(ndim), At(ndim))
  Qt(:)  = ZERO         ! Initial position centers
  SQt(:) = ONE          ! Initial widths (sqrt of Re(At))
  Pt(:)  = ZERO         ! Initial momenta
  At(:)  = ONE          ! Initial complex widths (simplified to 1 for default)

  ! If the propagator is Hagedorn-based, perform full Hagedorn propagation
  IF (propa%propa_name == 'hagedorn') THEN
     ! Step 1: Time propagation (temporary psi_dt)
     CALL march_temp(psi, psi_dt, t, propa, H)

     ! Step 2: Compute initial energy and norm (after temporary propagation)
     CALL Calc_Av_E(E0, psi_dt, H)
     CALL Calc_Norm_OF_psi(psi_dt, Norm0)

     ! Step 3: Project psi_dt into Hagedorn basis
     CALL Hagedorn_temp(psi, psi_dt, propa)

     ! Step 4: Reconstruct operator H with new basis
     DEALLOCATE(H%Scal_pot)
     CALL Calc_tab_Iq(Tab_iq, psi%Basis)
     CALL Set_Op(H, psi%Basis, Tab_iq)

     ! Step 5: Recalculate energy and norm for final psi
     CALL Calc_Av_E(E, psi, H)
     CALL Calc_Norm_OF_psi(psi, norm)

     ! Step 6: Output energy and norm differences for diagnostics
     WRITE(25,*) t, ABS(E0 - E), E0, E
     WRITE(24,*) t, ABS(Norm0 - norm), Norm0, norm

  ELSE
     ! Default propagation: simple time step without Hagedorn basis
     CALL march_temp(psi, psi_dt, t, propa, H)

     ! Directly update psi with temporary psi_dt
     psi%CVec(:) = psi_dt%CVec(:)
  END IF

END SUBROUTINE march_Global


 SUBROUTINE Poly_Hermite_modified_funct(Q,Qt,SQt,Bt,Pt,ib,d0)
     real(kind=Rkind)                               :: Q,Qt,SQt,Pt,Bt
     integer,intent(in)                             :: ib
     complex(kind=Rkind) , intent(inout)            :: d0
     complex(kind=Rkind)                            :: d0W,d1W,d2W
     real(kind=Rkind)                               :: DQ

       DQ = SQt*(Q-Qt)
       Pt = Pt/SQt
       Bt = Bt/(SQt*SQt)
       
      call Calc_d0d1d2W(Q,Qt,SQt,Bt,Pt,d0W,d1W,d2W,.false.)
      d0 = sqrt(SQt)*poly_Hermite(DQ,ib)*exp(-HALF*DQ*DQ)*d0W

 End SUBROUTINE

 SUBROUTINE Poly_Hermite_modified_funct_temp(Q,Qt,SQt,Bt,Pt,ib,d0)
    real(kind=Rkind) ,intent(in)                   :: Q,Qt,SQt,Pt,Bt
    integer,intent(in)                             :: ib
    complex(kind=Rkind) , intent(inout)            :: d0
    real(kind=Rkind)                               :: DQ,P
    complex(kind=Rkind)                            :: At

      DQ = SQt*(Q-Qt)
      At= complex(ONE,Bt/(SQt*SQt))
      P = Pt/SQt
     d0 = sqrt(SQt)*poly_Hermite(DQ,ib)*exp(-HALF*At*DQ*DQ+EYE*P*DQ)
End SUBROUTINE

 SUBROUTINE d0poly_Hermite_exp(Q,ib,d0)
    USE QDUtil_m
    IMPLICIT NONE
    complex(kind=Rkind),intent(inout) :: d0
    real(kind=Rkind),intent(in)       :: Q
    integer,intent(in)                :: ib
   
      d0 = poly_Hermite(Q,ib)*exp(-HALF*Q*Q)
    RETURN
 END SUBROUTINE


 SUBROUTINE d0poly_Hermite_exp_cplx(Q, Qt, At, Pt, ib, d0)
  USE QDUtil_m
  IMPLICIT NONE

  ! Input/output arguments
  complex(kind=Rkind), INTENT(INOUT) :: d0           ! Output: result of Hermite polynomial * exponential (complex)
  real(kind=Rkind),    INTENT(IN)    :: Q            ! Position variable
  real(kind=Rkind),    INTENT(IN)    :: Qt           ! Center of the Gaussian
  real(kind=Rkind),    INTENT(IN)    :: Pt           ! Momentum
  complex(kind=Rkind), INTENT(IN)    :: At           ! Complex width parameter 
  integer,            INTENT(IN)     :: ib           ! Hermite polynomial order

  ! Local variables
  complex(kind=Rkind) :: w                          ! Complex phase factor
  real(kind=Rkind)    :: SQ                         ! Square root of the real part of At
  real(kind=Rkind)    :: DQ                         ! Scaled displacement from Qt
  real(kind=Rkind)    :: Bt                         ! Imaginary part of At

  ! Compute square root of the real part of At (assumed to be positive)
  SQ = sqrt(real(At, kind=Rkind))

  ! Compute scaled coordinate (DQ = sqrt(Re(At)) * (Q - Qt))
  DQ = SQ * (Q - Qt)

  ! Extract the imaginary part of At
  Bt = aimag(At)

  ! Compute the complex exponential factor:
  !   exp[-(i/2) * Bt * (Q - Qt)^2 + i * Pt * (Q - Qt)]
  w = exp(-HALF * EYE * Bt * (Q - Qt)**2 + EYE * Pt * (Q - Qt))

  ! Call routine to compute the Hermite polynomial times Gaussian part
  CALL d0poly_Hermite_exp(DQ, ib, d0)

  ! Multiply the result by sqrt(SQ) and the complex exponential
  d0 = sqrt(SQ) * d0 * w

  ! Optional: write debug output
  ! WRITE(out_unit,*) 'l x,p,beta, d0 :', ib, Q, Pt, Bt, d0

  RETURN
END SUBROUTINE d0poly_Hermite_exp_cplx


end module Hagedorn_m