module Hagedorn_m
   USE QDUtil_m
   Use Basis_m
   Use Psi_m
   USE Ana_psi_m
   USE polyortho_m
   USE sub_propa_m
   implicit none
   PRIVATE
   PUBLIC :: Projection_temp,test_psi_temp, Calc_Basis_parameters,march_Hagedorn,march_Global

CONTAINS


SUBROUTINE Construct_Hagedorn_none_Variational_Basis(Basis,Qt,SQt,At,Pt)
  USE QDUtil_m
  IMPLICIT NONE
  TYPE(Basis_t),intent(inout)                     :: Basis
  real(kind=Rkind),intent(in)                     :: Qt(:), SQt(:),Pt(:)
  complex(kind=Rkind), intent(in)                 :: At(:)
  real(kind=Rkind)           , allocatable        :: Q(:),W(:)
  real(kind=Rkind) ,allocatable                   :: Q0(:), SQ0(:),P0(:)
  complex(kind=Rkind),allocatable                 :: A0(:)
  integer                                         :: ndim,ib,nb,nq

  ndim = size(Basis%tab_basis)-1
  allocate(Q0(ndim), SQ0(ndim),P0(ndim),A0(ndim))
  call Get_Basis_Parameters(Basis,Q0,SQ0,A0,P0)
  call Change_Basis_Parameters(Basis,Qt,SQt,At,Pt)
   call construct_primitive_basis_temp(Basis)

  DO ib = 1,ndim
   if (Basis%tab_basis(ib)%Basis_name == 'herm' .or.&
   & Basis%tab_basis(ib)%Basis_name == 'ho') then
   nb = Basis%tab_basis(ib)%nb
   nq = Basis%tab_basis(ib)%nq
   call  Calc_S(Basis%tab_basis(ib)%S,nb,nq,Q0(ib),SQ0(ib),A0(ib),&
   &P0(ib),Qt(ib),SQt(ib),At(ib),Pt(ib)) 
  
   End If  
  End Do

End SUBROUTINE


SUBROUTINE Construct_Hagedorn_none_Variational_Basis_temp(Basis,Qt,SQt,At,Pt)
  USE QDUtil_m
  IMPLICIT NONE
  TYPE(Basis_t),intent(inout)                     :: Basis
  real(kind=Rkind),intent(in)                     :: Qt(:), SQt(:),Pt(:)
  complex(kind=Rkind), intent(in)                 :: At(:)
  real(kind=Rkind)           , allocatable        :: Q(:),W(:)
  real(kind=Rkind) ,allocatable                   :: Q0(:), SQ0(:),P0(:)
  complex(kind=Rkind),allocatable                 :: A0(:)
  integer                                         :: ndim,ib,nb,nq

  ndim = size(Basis%tab_basis)-1
  allocate(Q0(ndim), SQ0(ndim),P0(ndim),A0(ndim))
  call Get_Basis_Parameters(Basis,Q0,SQ0,A0,P0)
  call Change_Basis_Parameters(Basis,Qt,SQt,At,Pt)
  call construct_primitive_basis_temp(Basis)

  DO ib = 1,ndim
   if (Basis%tab_basis(ib)%Basis_name == 'herm' .or.&
   & Basis%tab_basis(ib)%Basis_name == 'ho') then
   nb = Basis%tab_basis(ib)%nb
   nq = Basis%tab_basis(ib)%nq
   call  Calc_S(Basis%tab_basis(ib)%S,nb,nq,Q0(ib),SQ0(ib),A0(ib),&
   &P0(ib),Qt(ib),SQt(ib),At(ib),Pt(ib))   
   End If  
  End Do

End SUBROUTINE

SUBROUTINE Calc_Basis_parameters(psi,Qt,SQt,At,Pt)

  TYPE(psi_t),intent(in)                   :: psi
  real(kind=Rkind),intent(inout)           :: Qt(:), SQt(:),Pt(:)
  complex(kind=Rkind), intent(inout)       :: At(:)
   call Calc_AVQ_nD(psi0=psi, AVQ=Qt, SQ=SQt)
   call Calc_Av_imp_k_nD(psi,Pt)
   call  Calc_Avg_A_nD(psi, At)
End SUBROUTINE

SUBROUTINE Calc_Basis_parameters_temp(psi,Qt,SQt,At,Pt,propa)

  TYPE(psi_t),intent(in)                   :: psi
   type(propa_t),intent(in)                :: propa
  real(kind=Rkind),intent(inout)           :: Qt(:), SQt(:),Pt(:)
  complex(kind=Rkind), intent(inout)       :: At(:)

   Pt(:) = ZERO
   call Calc_AVQ_nD(psi0=psi, AVQ=Qt, SQ=SQt)
   At(:) = SQt(:)*SQt(:)
  If(propa%P)  call Calc_Av_imp_k_nD(psi,Pt)
  If(propa%Beta) call  Calc_Avg_A_nD(psi, At)

End SUBROUTINE


SUBROUTINE Calc_S(S,nb,nq,Q0,SQ0,A0,P0,Qt,SQt,At,Pt)
IMPLICIT NONE
 complex(kind=Rkind), intent(inout)             :: S(:,:)
 complex(kind=Rkind), intent(in)                :: At,A0
 real(kind=Rkind),intent(in)                    :: Q0,SQ0,P0,Qt,SQt,Pt
 integer ,intent(in)                            :: nb,nq
 integer                                        :: ib,iq,jb
 real(kind=Rkind)                               :: SQeq,Qeq
 real(kind=Rkind)           , allocatable       :: Q(:),W(:)
  real(kind=Rkind)                              :: Bt,B0
 complex(kind=Rkind),allocatable                :: d0gb(:,:),d0bgw(:,:)

 SQeq = sqrt(SQ0*SQ0 + SQt*SQt)/sqrt(TWO)
 Qeq = (SQ0*SQ0*Q0 + SQt*SQt*Qt)/(SQ0*SQ0 + SQt*SQt)
 
 Bt  = aimag(At) 
 B0  = aimag(A0)
 
 allocate(Q(nq),W(nq))
 call hercom(nq, Q(:),W(:))
 w(:) = W(:)/SQeq
 Q (:) = Qeq + Q(:)/SQeq
 allocate (d0gb(nq, nb))
 allocate (d0bgw(nb, nq))

  print*,'Q0,A0,P0',Q0,A0,P0
  print*,'Qt,At,Pt',Qt,At,Pt
  !print*,'SQeq,Qeq',SQeq,Qeq

 DO iq = 1, nq
   DO ib = 1, nb
    call d0poly_Hermite_exp_cplx(Q(iq),Q0,A0,P0,ib-1,d0gb(iq, ib))
    call d0poly_Hermite_exp_cplx(Q(iq),Qt,At,Pt,ib-1,d0bgw(ib, iq))
     d0bgw(ib, iq) = d0bgw(ib, iq)*W(iq)
   END DO
END DO
  S = matmul(conjg(d0bgw),d0gb)

 call  Write_VecMat(S, out_unit, 5,  info='S')
 deallocate(d0gb,d0bgw,w,Q)

End SUBROUTINE

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
    write(25,*) E0,E
     write(24,*) Norm0,Norm
    write (out_unit, *) 'Begin Hagedorn projection  E0,Norm0',E0,Norm0
    write (out_unit, *) 'END Hagedorn projection E,Norm ',E,Norm
END SUBROUTINE 



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

       call  Hagedorn_Inv(psi1, psi2)

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

SUBROUTINE Hagedorn_temp(psi, psi_dt,propa)
USE psi_m
type(psi_t),  intent(inout)             :: psi
type(psi_t),  intent(inout)             :: psi_dt
TYPE(propa_t), intent(in)               :: propa
 complex(kind=Rkind),allocatable        :: At(:)
 real(kind=Rkind)   ,allocatable        :: Qt(:),SQt(:),Pt(:)
  real(kind= Rkind)                     :: Norm0,norm
  integer                               :: ndim
  ndim = size(psi%Basis%tab_basis) - 1
 allocate(Qt(ndim), SQt(ndim),Pt(ndim),At(ndim))
 call  Calc_Basis_parameters_temp(psi_dt,Qt,SQt,At,Pt,propa)
 call Construct_Hagedorn_none_Variational_Basis_temp(psi_dt%Basis,Qt,SQt,At,Pt)
 call Calc_Norm_OF_psi(psi_dt,Norm0)
 call Projection(psi,psi_dt)
 call Calc_Norm_OF_psi(psi,Norm)
 write(out_unit,*)  abs(Norm0-Norm), Norm0,Norm
    deallocate(Qt,SQt,At,Pt)
 end subroutine


SUBROUTINE Hagedorn_Inv(psi1, psi2)
  USE psi_m
  type(psi_t),  intent(inout)   ,target     :: psi1
  type(psi_t),  intent(in),target           :: psi2
  complex(kind=Rkind),allocatable           :: At(:),A0(:),S(:,:)
  real(kind=Rkind)   ,allocatable           :: Qt(:),SQt(:),Pt(:)
  real(kind=Rkind)   ,allocatable           :: Q0(:),SQ0(:),P0(:)
  Integer, allocatable                      :: Ib1(:), Ib2(:), Ib3(:)
  complex(kind= Rkind), pointer             :: BBB1(:, :, :), BBB2(:, :, :)
  complex(kind= Rkind), allocatable, target :: B1(:), B2(:)
  real(kind= Rkind)                         :: Norm0,norm
  integer                                   :: ndim,nb,nq,ib


  call Calc_index(Ib1=Ib1, Ib2=Ib2, Ib3=Ib3, Basis=psi1%Basis)
  allocate (B1(Ib1(1)*Ib2(1)*Ib3(1)))
  allocate (B2(Ib1(1)*Ib2(1)*Ib3(1)))
  ndim = size(psi1%Basis%tab_basis) - 1
  allocate(Qt(ndim), SQt(ndim),Pt(ndim),At(ndim))
  allocate(Q0(ndim), SQ0(ndim),P0(ndim),A0(ndim))
    call Get_Basis_Parameters(psi1%Basis,Q0,SQ0,A0,P0)
    call Get_Basis_Parameters(psi2%Basis,Qt,SQt,At,Pt)
     B2(:) = psi2%CVec(:)
     B1(:) = CZERO
     DO ib = 1,ndim
      nb = psi1%Basis%tab_basis(ib)%nb
      nq = psi1%Basis%tab_basis(ib)%nq
      allocate(S(nb,nb))
      call Calc_S(S,nb,nq,Qt(ib),SQt(ib),At(ib),Pt(ib),Q0(ib),SQ0(ib),A0(ib),P0(ib))
       BBB1(1:Ib1(ib), 1:Ib2(ib), 1:Ib3(ib)) => B1
       BBB2(1:Ib1(ib), 1:Ib2(ib), 1:Ib3(ib)) => B2
       call projection_1D_temp(BBB1, BBB2, S)
       B2(:) = B1(:)
       deallocate(S)
    END DO
   call Calc_Norm_OF_psi(psi2,Norm0)
   call Projection(psi1,psi2)
   call Calc_Norm_OF_psi(psi1,Norm)

   deallocate(Qt,SQt,At,Pt)
   deallocate(Q0,SQ0,A0,P0)
END SUBROUTINE

 
  SUBROUTINE march_Global(psi, psi_dt, t, propa)
  USE psi_m
  type(propa_t),intent(IN)                :: propa
  type(psi_t),  intent(inout)             :: psi
  type(psi_t),  intent(INOUT)             :: psi_dt
  real(kind=Rkind), intent(IN)            :: t

   complex(kind=Rkind),allocatable        :: At(:)
   real(kind=Rkind)   ,allocatable        :: Qt(:),SQt(:),Pt(:)
    real(kind= Rkind)                     :: Norm0,norm,E,E0
    integer                               :: ndim

    ndim = size(psi%Basis%tab_basis) - 1
    allocate(Qt(ndim), SQt(ndim),Pt(ndim),At(ndim))
    Qt(:)=ZERO; SQt(:)=ONE;Pt(:)=ZERO;At(:)=ONE

   If (propa%propa_name == 'hagedorn') Then
    call  march(psi, psi_dt, t, propa)
     call write_psi(psi=psi_dt, psi_cplx=.false., print_psi_grid=.false. &
     , print_basis=.false., t=t, int_print=23, real_part=.false.)
     write(23,*)
     call  Calc_average_energy(psi_dt, E0)
     call Calc_Norm_OF_psi(psi_dt,Norm0)
     call Hagedorn_temp(psi, psi_dt,propa)
     call  Calc_average_energy(psi, E)
     call Calc_Norm_OF_psi(psi,Norm)
     write(25,*) t, abs(E0-E), t,E0,E
    write(24,*) t, abs(Norm0-Norm),t, Norm0,Norm
   Else
     call  march(psi, psi_dt, t, propa)
     psi%CVec(:) = psi_dt%CVec(:)
   End If  
  
 END SUBROUTINE

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



 SUBROUTINE d0poly_Hermite_exp_cplx(Q,Qt,At,Pt,ib,d0)
    USE QDUtil_m
    IMPLICIT NONE
    complex(kind=Rkind),intent(inout) :: d0
    real(kind=Rkind),intent(in)       :: Q,Qt,Pt
    complex(kind=Rkind),intent(in)    :: At
    integer,intent(in)                :: ib
    real(kind=Rkind)                  :: SQ,DQ,Z_K,Bt

    SQ  =  sqrt(real(At,kind=Rkind)) 
    DQ  = SQ*(Q-Qt)
    Bt =  aimag(At)
    Z_K = mod(HALF*Bt*(Q-Qt)*(Q-Qt)-Pt*(Q-Qt),TWO*PI)

     call d0poly_Hermite_exp(DQ,ib,d0)
     d0 = sqrt(SQ)*d0*exp(-EYE*Z_K)
   !write(out_unit,*) 'l x,p,beta, d0 :',l,x,d0
    RETURN

 END SUBROUTINE

end module Hagedorn_m