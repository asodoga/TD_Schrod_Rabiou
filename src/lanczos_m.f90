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
  USE QDUtil_m
  USE Op_m
  USE psi_m
  implicit none

contains


SUBROUTINE Calc_Beta(Beta,K,Basis,nb_k)
  USE QDUtil_m
  implicit none
  complex(kind=Rkind),intent(in)             :: K(:,:)
  TYPE(Basis_t),intent(in)                   :: Basis
  real(kind=Rkind),allocatable,intent(inout) :: Beta(:)
  integer,intent(in)                         :: nb_k

  TYPE(psi_t)                                :: psi,psi_k,Hpsi_k
  real(kind=Rkind),allocatable               :: Mat(:,:)
  TYPE(Op_t)                                 :: H
  integer                                    :: I_k,J_k

  call init_psi(psi, Basis, cplx=.true., grid=.false.)
  call init_psi(psi_k, Basis, cplx=.true., grid=.false.)
  call init_psi(Hpsi_k,Basis, cplx=.true., grid=.false.)

  If(allocated(Beta)) deallocate(Beta)
  allocate (Mat(nb_k,nb_k),Beta(nb_k-1))


  Mat(:,:) = CONE

  Do I_k = 1,nb_k

  psi%CVec(:) = K(:,I_k)

  Do J_k= 1,nb_k

   psi_k%CVec(:) = K(:,J_k)
   call calc_OpPsi(H, psi_k, Hpsi_k)
   Mat(I_k,J_k) = dot_product(psi%CVec,Hpsi_k%CVec)

  END DO

END Do

   !CALL  Write_VecMat(Mat, out_unit, 5,  info='Mat')

   Beta(:) = ZERO

 Do I_k = 1,nb_k-1

   Beta(I_k) = Mat(I_k,I_k+1)

 End Do

call dealloc_psi(psi)
call dealloc_psi(psi_k)
call dealloc_psi(Hpsi_k)
deallocate(Mat)

END SUBROUTINE


   
SUBROUTINE Krilov_Gram_Schmidt_basis(K,psi,epsi_precision,deltat_t)
   USE QDUtil_m
   implicit none
    complex(kind=Rkind),allocatable,intent(inout) :: K(:,:)
    TYPE(psi_t),intent(in)                        :: psi
    real(kind=Rkind),intent(in)                   :: epsi_precision,deltat_t
   
    ! locals variables --------------------------------------------------------------
   
    integer                                       :: nb,Ib_k  
    real(kind=Rkind) ,allocatable                 :: Beta(:)
    complex(kind=Rkind),allocatable               :: Ktemp(:,:)
    real(kind=Rkind)                              :: epsi_precision_temp
   
    nb = size(psi%CVec)
    epsi_precision_temp = ZERO
            
     Ib_k = 0

      Do 
       
       Ib_k = Ib_k+1
       CALL  Krilov_Gram_Schmidt_basis_temp(Ktemp,psi,Ib_k)
       call Calc_Beta(Beta,Ktemp,psi%Basis,Ib_k)
       epsi_precision_temp = abs(((deltat_t**size(Beta))*product(Beta(:)))/Factorial(size(Beta)))
       write (out_unit, *) 'Ib_k,|psi_Ib_k+1-psi_Ib_k|=:',Ib_k,epsi_precision_temp
       If(epsi_precision_temp<=epsi_precision) then
         write (out_unit, *) 
         write (out_unit, *) '--------------------------------------------------------------------------------------------------'
         write (out_unit, *) 'SIL condition is fulfild after Ib_k =:', Ib_k, 'iteration'
         write (out_unit, *) 'deltat_t =:',deltat_t
         write (out_unit, *) '|psi_Ib_k+1-psi_Ib_k|=:',epsi_precision_temp
         write (out_unit, *) 'Krylov Basis will be constructed with size  nb_max =:',Ib_k 
         write (out_unit, *) '--------------------------------------------------------------------------------------------------'
         K = Ktemp
       exit
       elseif(epsi_precision_temp >= TEN**10)then
        write (out_unit, *) 
        write (out_unit, *) '--------------------------------------------------------------------------------------------------'
        write (out_unit, *) 'deltat_t =:',deltat_t
        write (out_unit, *) '|psi_Ib_k+1-psi_Ib_k|=:',epsi_precision_temp
        write (out_unit, *) 'Wrong  deltat_t,Check your Initial Lanczos data and try again'
        write (out_unit, *) '--------------------------------------------------------------------------------------------------'
        stop 'sorry check your result file for mor information'
        elseif(Ib_k==nb) then
        write (out_unit, *)
        write (out_unit, *)  'Krylov Basis size is greter than psi%Basis%nb '
        write (out_unit, *) 'psi%Basis%nb =:',psi%Basis%nb
        write (out_unit, *) 'Krylov nb_max =:',Ib_k
        write (out_unit, *) 'Check your data and try again'
          stop 'sorry check your result file for mor information'

       End If
     
      END Do
     deallocate(Ktemp,Beta)
  END SUBROUTINE


SUBROUTINE Krilov_Gram_Schmidt_basis_temp(K,psi,nb_k)
USE QDUtil_m
implicit none
 complex(kind=Rkind),allocatable,intent(inout) :: K(:,:)
 TYPE(psi_t),intent(in)                        :: psi
 integer,intent(in)                            :: nb_k

 ! locals variables --------------------------------------------------------------

 Complex(kind=Rkind) ,allocatable              :: Q(:,:),V(:)
 TYPE(Op_t)                                    :: H
 TYPE(psi_t)                                   :: psi0,Hpsi
 integer                                       :: I_k,J_k,nb      

 nb = size(psi%CVec)
 If(allocated(K)) deallocate(K)
 allocate(Q(nb,nb_k),K(nb,nb_k),V(nb))
 CALL init_psi(psi0, psi%Basis, cplx=.TRUE., grid=.false.)
 CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.false.)

   K = CZERO
   Q = CZERO
   Q(:,1) = psi%CVec(:)
   psi0%CVec(:) = psi%CVec(:)

 ! construction of krylov normalize basis  --------------------------------------. 
 ! loop for H|psi>, H^2|psi>, H^3|psi>......H^nb_k|psi>.
 Do I_k=2,nb_k

  call calc_OpPsi(H, psi0,Hpsi) 
  Q(:,I_k) = Hpsi%CVec(:)
  psi0%CVec(:) =Hpsi%CVec(:)
  Hpsi%CVec(:) =CZERO

 END Do

! 1st Normalisation With gram shmidt--------------------------------------------------.

 K(:,1) = Q(:,1)
 K(:,1)= K(:,1)/sqrt(dot_product(K(:,1),K(:,1)))

 Do I_k=2,nb_k
  V(:) = Q(:,I_k)
   Do J_k=1,I_k-1
      V(:) = V(:)-dot_product(K(:,J_k),Q(:,I_k))*K(:,J_k)
   End do
   K(:,I_k) = V(:)/sqrt(dot_product(V,V))
   V= CZERO
 END do

! 2cnd Normalisation ---------------------------------------------------------------.
K(:,1)= K(:,1)/sqrt(dot_product(K(:,1),K(:,1)))
Do I_k=2,nb_k
 V(:) = K(:,I_k)
  Do J_k=1,I_k-1
     V(:) = V(:)-dot_product(K(:,J_k),K(:,I_k))*K(:,J_k)
  End do
  K(:,I_k) = V(:)/sqrt(dot_product(V,V))
  V= CZERO
END do

   call dealloc_psi(psi0)
   call dealloc_psi(Hpsi)
   deallocate(V,Q)

  END SUBROUTINE



SUBROUTINE Construct_Krylov_matrix(Mat_k,K,psi,epsi_precision,deltat_t)
   USE QDUtil_m
   implicit none
   complex(kind=Rkind),allocatable,intent(inout)             :: Mat_k(:,:)
   complex(kind=Rkind),allocatable,intent(inout)             :: K(:,:)
   TYPE(psi_t),intent(in)                                    :: psi
   real(kind=Rkind),intent(in)                               :: epsi_precision,deltat_t

   !locals variables -----------------------------------------------------
   TYPE(Op_t)                                                :: H
   integer                                                   :: I_k,J_k,nb_k
   TYPE(psi_t)                                               :: psi0,psi_k,Hpsi_k

  CALL init_psi(psi0, psi%Basis,  cplx=.true., grid=.false.)
  CALL init_psi(psi_k, psi%Basis,  cplx=.true., grid=.false.)
  CALL init_psi(Hpsi_k, psi%Basis, cplx=.true., grid=.false.)

  call Krilov_Gram_Schmidt_basis(K,psi,epsi_precision,deltat_t)
  !CALL  Write_VecMat(K, out_unit, 5,  info='Ki')
   nb_k = size(K, dim = 2)
   allocate (Mat_k(nb_k,nb_k))

   psi0%CVec(:) = CZERO
   psi_k%CVec(:) = CZERO
   Hpsi_k%CVec(:) = CZERO

   Do I_k = 1,nb_k

     psi0%CVec(:) = K(:,I_k)

     Do J_k= 1,nb_k

      psi_k%CVec(:) = K(:,J_k)
      call calc_OpPsi(H, psi_k, Hpsi_k)
      Mat_k(I_k,J_k) = dot_product(psi0%CVec,Hpsi_k%CVec)

     END DO

   END do
 !CALL  Write_VecMat(Mat, out_unit, 5,  info='Mat')
 call dealloc_psi(psi0)
 call dealloc_psi(psi_k)
 call dealloc_psi(Hpsi_k)

END SUBROUTINE

SUBROUTINE TEST_Lonaczos_cplx(psi)
   USE QDUtil_m
  implicit none
   type(psi_t),intent(in)                        :: psi 
   real (kind=Rkind), allocatable                :: Eig(:)
   complex (kind=Rkind), allocatable             :: H(:,:),V(:,:)
   complex (kind=Rkind), allocatable             :: Vec_Basis(:,:)

   complex (kind=Rkind), allocatable             :: K(:,:)
   real(kind=Rkind),allocatable                  :: Beta(:)
   real(kind=Rkind)                              :: epsi_precision,deltat_t
   integer                                       :: n


   
  
   !CALL  Write_VecMat(K, out_unit, 5,  info='K')
  
  epsi_precision = ONETENTH**10
  deltat_t=ONETENTH
 !call  Krilov_Gram_Schmidt_basis0(K,psi,epsi_precision,deltat_t)
 !call Construct_triband_matrix(H,psi,epsi_precision,deltat_t)
 !CALL  Write_VecMat(H, out_unit, 5,  info='H')
   n = size(H,dim=1)
    write (out_unit, *) 'n',n
     allocate(Eig(n ))
     allocate(V(n ,n))
    ! CALL  diagonalization(H,Eig,V,n) 
     ! CALL  Write_VecMat(V, out_unit, 5,  info='V') 
       write (out_unit, *) 'Eig',Eig


  END SUBROUTINE



 SUBROUTINE  Calc_psi_step_cplx(psi_dt,psi,epsi_precision,deltat_t)
  
     USE QDUtil_m
     implicit none
     TYPE(psi_t), intent(in)                       :: psi
     real(kind=Rkind), intent(in)                  :: deltat_t,epsi_precision
     TYPE(psi_t), intent(inout)                    :: psi_dt
     
     
     !locals variables-----------------------------------------------------------------
     
     complex (kind=Rkind), allocatable             :: Eig_v(:)
     complex (kind=Rkind), allocatable             :: K(:,:)
     complex (kind=Rkind), allocatable             :: Eig_vec(:,:),H(:,:)
     complex (kind=Rkind), allocatable             :: Eig_vec_basis(:,:),C1(:),C2(:)
     integer                                       :: nb_k,nb
      
      call Construct_Krylov_matrix(H,K,psi,epsi_precision,deltat_t)      
      nb_k =size(H,dim = 1)
      nb =size(K,dim = 1)
      allocate(Eig_vec_basis(nb,nb_k),C1(nb_k),C2(nb_k))
      allocate(Eig_vec(nb_k,nb_k),Eig_v(nb_k))

      call  diagonalization(H,Eig_v,Eig_vec,nb_k) 
      !CALL  Write_VecMat(Eig_vec, out_unit, 5,  info='Vec') 
     ! write (out_unit, *) 'Eig',Eig_v
 
      Eig_vec_basis(:,:) =CZERO
      Eig_vec_basis = matmul(K, Eig_vec)
      C1 = matmul(conjg(transpose(Eig_vec_basis)),psi%CVec)
      C2(:) = CZERO
      C2(:) = C1(:)*exp(-EYE*Eig_v(:)*deltat_t)
      psi_dt%CVec =  matmul(Eig_vec_basis, c2)
   
      deallocate(Eig_vec_basis,K,Eig_vec,Eig_v,H,C1,C2)

  END SUBROUTINE

    end module lanczos_m
