!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!
!> MODULE: lanczos_m
!
!> This module is used to do short iterative lanczos evolutions on a 
!> symetric tridiagonal matrix.
!----------------------------------------------------------------------------------
module lanczos_m
  USE UtilLib_m
  USE Op_m
  USE psi_m
  USE diago_m
  implicit none

contains

SUBROUTINE Krilov_Gram_Schmidt_basis(K,psi,num_max)
USE UtilLib_m
implicit none
 complex(kind=Rk),allocatable,intent(inout) :: K(:,:)
 TYPE(psi_t),intent(in)                     :: psi
 integer,intent(in)                         :: num_max

 ! locals variables ======================================================================================

 Complex(kind=Rk) ,allocatable              :: Q(:,:),V(:)
 TYPE(Op_t)                                 :: H
 TYPE(psi_t)                                :: psi_b1,psi_b2
 integer                                    :: i,j,nb
 real(kind=Rk) ,allocatable                 :: S(:,:)   ! Overlap matrix          


 nb = size(psi%CVec)
 allocate(Q(nb,num_max),K(nb,num_max),V(nb))
 CALL init_psi(psi_b1, psi%Basis, cplx=.TRUE., grid=.false.)
 CALL init_psi(psi_b2, psi%Basis, cplx=.TRUE., grid=.false.)

 ! loop for H|psi>, H^2|psi>, H^3|psi>...
 call calc_OpPsi(H, psi, psi_b1)
 Q(:,1) = psi_b1%CVec(:)

 ! construction de la base de krylov non normalisee========================================. 

 Do i=2,num_max

  call calc_OpPsi(H, psi_b1,psi_b2) 
  Q(:,i) = psi_b2%CVec(:)
  psi_b1%CVec(:) =psi_b2%CVec(:)
  psi_b2%CVec(:) =CZERO

 END Do

! FIRST Normalisation avec gram shmidt============================================================

 K(:,1) = Q(:,1)
 K(:,1)= K(:,1)/sqrt(dot_product(K(:,1),K(:,1)))

 Do i=2,num_max
  V(:) = Q(:,i)
   Do j=1,i-1
      V(:) = V(:)-dot_product(K(:,j),Q(:,i))*K(:,j)
   End do
   K(:,i) = V(:)/sqrt(dot_product(V,V))
   V= CZERO
 END do

! SECOND NORMALISATION =============================================================
K(:,1)= K(:,1)/sqrt(dot_product(K(:,1),K(:,1)))
Do i=2,num_max
 V(:) = K(:,i)
  Do j=1,i-1
     V(:) = V(:)-dot_product(K(:,j),K(:,i))*K(:,j)
  End do
  K(:,i) = V(:)/sqrt(dot_product(V,V))
  V= CZERO
END do

  !S = matmul(conjg(transpose(K)),K)
!
  !do i = 1,num_max
  !  write(out_unitp,*) (s(i,j),j=1,num_max)
  !end do

   call dealloc_psi(psi_b1)
   call dealloc_psi(psi_b2)
   deallocate(V,Q)
  END SUBROUTINE



SUBROUTINE Construct_triband_matrix(Mat,K,psi,num_max)
   real(kind=Rk),allocatable,intent(inout)                :: Mat(:,:)
   complex(kind=Rk),allocatable,intent(inout) ,optional   ::K(:,:)
   TYPE(psi_t),intent(in)                                 :: psi
   integer,intent(in)                                     :: num_max

   !locals variables =======================================================================

   complex(kind=Rk) ,allocatable                          :: K1(:,:)
   TYPE(Op_t)                                             :: H
   integer                                                :: i,j,nb
   TYPE(psi_t)                                            :: psi_b1,psi_b2,Hpsi_b2

  CALL init_psi(psi_b1, psi%Basis, cplx=.TRUE., grid=.false.)
  CALL init_psi(psi_b2, psi%Basis, cplx=.TRUE., grid=.false.)
  CALL init_psi(Hpsi_b2, psi%Basis, cplx=.TRUE., grid=.false.)
 
  nb= size(psi%CVec)
  if(present(K))  allocate (K(nb,num_max))
  allocate (Mat(num_max,num_max))

   call Krilov_Gram_Schmidt_basis(K1,psi,num_max)
   if(present(K)) K(:,:) = K1(:,:)

   Do i = 1,num_max
     psi_b1%CVec(:) = K1(:,i)

     Do j= 1,num_max
      psi_b2%CVec(:) = K1(:,j)
      call calc_OpPsi(H, psi_b2, Hpsi_b2)
      Mat(i,j) = dot_product(psi_b1%CVec,Hpsi_b2%CVec)
     END DO

   END do
 
   !do i = 1,num_max
   !  write(out_unitp,*) (Mat(i,j),j=1,num_max)
   !end do


END SUBROUTINE



 
 
 SUBROUTINE Lanczos_eign_syst_solve(EigenVal,Vec_Basis,psi,kmax)
      Use diago_m
      Use psi_m
      type(psi_t),intent(in)                        :: psi 
      real (kind=Rk), allocatable ,intent(inout)    :: EigenVal(:)
      complex (kind=Rk), allocatable ,intent(inout) :: Vec_Basis(:,:)
      integer,intent(in)                            :: kmax
   
   !locals variables===================================================================================
      real(kind=Rk) ,allocatable                    :: Matrix(:,:)
      real (kind=Rk), allocatable                   :: EigenVec(:,:)
      complex (kind=Rk), allocatable                :: K(:,:)
      !logical,          parameter                  :: debug = .true.
      logical,         parameter                    :: debug = .false.
      integer                                       :: ib,jb,nb
      
      
      IF (debug) THEN
          write(out_unitp,*) 'BEGINNING Eig_syst_solve'
          flush(out_unitp)
      END IF
      
       nb = size(psi%CVec)
      allocate(EigenVal(kmax ))
      allocate(EigenVec(kmax ,kmax ))
      allocate(Vec_Basis(nb ,kmax ))
      
      CALL  Construct_triband_matrix(Matrix,K,psi,kmax)
      CALL  diagonalization(matrix,EigenVal,EigenVec,kmax )
      DO ib= 1,kmax
         Vec_Basis(:,ib) =CZERO
        DO jb = 1,kmax
         Vec_Basis(:,ib)= Vec_Basis(:,ib)+EigenVec(jb,ib)*K(:,jb)
        END DO
      END DO
 
    ! Write(out_unitp,*)
    ! Write(out_unitp,*)
     
    !Write(out_unitp,*) 'eigenvalues = '
    !DO ib=1,kmax
    !   write(out_unitp,*) EigenVal(ib)
    !END DO
     
    !Write(out_unitp,*)
    !Write(out_unitp,*) 'EigenVec'
     
    ! DO ib=1,kmax
    !     write(*,*) (EigenVec(ib,jb),jb=1,kmax)
    !END DO

    !Write(out_unitp,*)
    !Write(out_unitp,*) 'Vec_Basis'
    !DO ib=1,kmax
    ! write(*,*) (Vec_Basis(ib,jb),jb=1,kmax)
    !END DO

    
      IF (debug) THEN
          write(out_unitp,*) 'END Eig_syst_solve'
          flush(out_unitp)
      END IF
      
  END SUBROUTINE 

SUBROUTINE TEST_Lonaczos_cplx(psi,kmax)
   integer,intent(in)                            :: kmax
   type(psi_t),intent(in)                        :: psi 
   real (kind=Rk), allocatable                   :: EigenVal(:)
   complex (kind=Rk), allocatable                :: Vec_Basis(:,:)


   CALL Lanczos_eign_syst_solve(EigenVal,Vec_Basis,psi,kmax)

  END SUBROUTINE



 SUBROUTINE  Calc_psi_step_cplx(psi_dt,psi,dt,Kmax)
  
  
     TYPE(psi_t), INTENT(INOUT)                    :: psi_dt
     TYPE(psi_t), INTENT(IN)                       :: psi
     integer,intent(IN)                            :: kmax
     real(kind=Rk), INTENT(IN)                     :: dt
     
     
     !locals variables===================================================================
     
     real (kind=Rk), allocatable                   :: EigenVal(:)
     complex (kind=Rk), allocatable                :: Vec_Basis(:,:)
     complex (kind=Rk), allocatable                :: C1(:),C2(:)
      
      allocate(C1(kmax),C2(kmax))
      
      CALL Lanczos_eign_syst_solve( EigenVal,Vec_Basis,psi,kmax)
      
      C1 = matmul(conjg(transpose(Vec_Basis)),psi%CVec)
      C2(:) = C1(:)*exp(-dt*EYE*EigenVal(:))
      psi_dt%CVec = matmul(Vec_Basis,C2)
       
       deallocate(C1,C2,Vec_Basis,EigenVal) 
  
  END SUBROUTINE
  

    end module lanczos_m
