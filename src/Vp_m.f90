!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!
!> MODULE: Vp_m
!
!> This module is used to perform Variatrional priciple methode in solving TDSE 
!> .
!----------------------------------------------------------------------------------
module Vp_m
  USE QDUtil_m
  USE Basis_m
  USE psi_m
  implicit none

contains


SUBROUTINE psiTOLambda(L,psi)
   Implicit none
   complex(Kind=Rkind) ,intent(inout) :: L(:)
   TYPE(psi_t),intent(in)             :: psi
   
   !> Locals variables ------------------------------------------------
   
   integer                            :: Ib,nb,ndim,nsurf
   !logical, parameter                 :: debug = .true.
   logical, parameter                 :: debug = .false.
   real(kind=Rkind), allocatable      :: Qt(:),Pt(:)
   complex(kind=Rkind), allocatable   :: At(:)
   
   
   
   nb = psi%Basis%nb
   ndim = size(psi%Basis%tab_basis)-1
   nsurf = psi%Basis%tab_basis(ndim+1)%nb
   allocate (Qt(ndim),Pt(ndim),At(ndim))

   IF (debug) THEN
    write (out_unit, *) 'psi in psiTOLambda:--------------------------------------------'
    Do Ib = 1,nb*nsurf
       write (out_unit, *) Ib,psi%CVec(Ib)
    End Do
  End IF
   
   L(:) = CZERO

   Do Ib = 1,ndim
    Qt(Ib) = psi%Basis%tab_basis(Ib)%Q0 
    Pt(Ib) = psi%Basis%tab_basis(Ib)%Imp_k 
    At(Ib) = psi%Basis%tab_basis(Ib)%alpha 
   End DO

   L(1:nb*nsurf) = psi%CVec(1:nb*nsurf)
   L(nb*nsurf+1:nb*nsurf+ndim) = At(1:ndim)
   L(nb*nsurf+ndim+1:nb*nsurf+2*ndim) = Qt(1:ndim)
   L(nb*nsurf+2*ndim+1:nb*nsurf+3*ndim) = Pt(1:ndim)
   
    IF (debug) THEN
      write (out_unit, *) 'Lambda in psiTOLambda:----------------------------------------------------------'
      Do Ib = 1,nb*nsurf+3*ndim
         write (out_unit, *) Ib,L(Ib)
      End Do
    End IF
   
   deallocate (Qt,Pt,At)

END SUBROUTINE

SUBROUTINE LambdaTOpsi(psi,L)
   Implicit none
   complex(Kind=Rkind) ,intent(in)    :: L(:)
   TYPE(psi_t),intent(inout)          :: psi
   
   !> Locals variables ------------------------------------------------
   
   integer                            :: Ib,nb,ndim,nsurf
   !logical, parameter                 :: debug = .true.
   logical, parameter                 :: debug = .false.
   real(kind=Rkind), allocatable      :: SQt(:), Qt(:),Pt(:)
   complex(kind=Rkind), allocatable   :: At(:)
   
   
   
   nb = psi%Basis%nb
   ndim = size(psi%Basis%tab_basis)-1
   nsurf = psi%Basis%tab_basis(ndim+1)%nb
   allocate (SQt(ndim),Qt(ndim),Pt(ndim),At(ndim))

   IF (debug) THEN
      write (out_unit, *) 'Lambda in LambdaTOpsi:-------------------------------------------------'
      Do Ib = 1,nb*nsurf+3*ndim
         write (out_unit, *) Ib,L(Ib)
      End Do
    End IF
   psi%CVec(:) = CZERO

    psi%CVec(1:nb*nsurf)= L(1:nb*nsurf) 
    At(1:ndim) = L(nb*nsurf+1:nb*nsurf+ndim)
    Qt(1:ndim) = L(nb*nsurf+ndim+1:nb*nsurf+2*ndim)
    Pt(1:ndim) = L(nb*nsurf+2*ndim+1:nb*nsurf+3*ndim)

   Do Ib = 1,ndim
       psi%Basis%tab_basis(Ib)%Q0 =Qt(Ib) 
       psi%Basis%tab_basis(Ib)%Imp_k =Pt(Ib)  
       psi%Basis%tab_basis(Ib)%alpha =At(Ib) 
   End DO

   SQt(:) = ONE
   call construct_primitive_basis(psi%Basis,Qt,Pt,At,SQt) 
   deallocate (SQt,Qt,Pt,At)



   IF (debug) THEN
      write (out_unit, *) 'psi in LambdaTOpsi:-----------------------------------------------------'
      Do Ib =1, nb*nsurf 
         write (out_unit, *) Ib,psi%CVec(Ib)
      End Do
      flush(out_unit)
    End IF  

END SUBROUTINE




 SUBROUTINE Vp_test(psi)
    TYPE(psi_t), intent(inout)          :: psi

    !> Locals variables ------------------------------------------------
    complex(Kind=Rkind) ,allocatable    :: L(:)
    integer                             :: Ib,nb,ndim,nsurf

    nb = psi%Basis%nb
    ndim = size(psi%Basis%tab_basis)-1
    nsurf = psi%Basis%tab_basis(ndim+1)%nb
    allocate (L(nb*nsurf+3*ndim))
  
    print*,'----------------Debut du test sur LambdaTOpsi psiTOLambda----------'
    
    call psiTOLambda(L,psi)
    call LambdaTOpsi(psi,L)
     print*,'---------------- Fin du test sur LambdaTOpsi psiTOLambda------------'
End SUBROUTINE



SUBROUTINE Find_I(Tab_Ib,I,Basis)
   implicit none
   integer,intent(in)                   ::I
   TYPE(Basis_t),intent(in)             :: Basis
   integer,intent(inout)                ::Tab_Ib(:)

   !> Locals variables ------------------------------------------------
   !logical, parameter                  :: debug = .true.
   logical, parameter                   :: debug = .false.
   integer                              :: Ib
   logical                              :: Endloop_b


   IF (debug) THEN
     write (out_unit, *) 'Tab_Ib :-------------------------------------------------'
   End IF

   Tab_Ib(:) = 0
   Call Init_tab_ind(Tab_Ib, Basis%NDindexb)
   DO Ib = 1,I
    CALL increase_NDindex(Tab_Ib, Basis%NDindexb, Endloop_b)
   End DO

   IF (debug) THEN
     write (out_unit, *) 'I , Tab_Ib',I,Tab_Ib
     write (out_unit, *) 'Tab_Ib :-------------------------------------------------'
   End IF
End SUBROUTINE



SUBROUTINE Calc_SIJ(S_IJ,I,J,Basis)
  implicit none
  
  TYPE(Basis_t),intent(in)             :: Basis
  real(kind=Rkind),intent(inout)       :: S_IJ
  integer ,intent(in)                  :: I,J
  
  !> Locals variables ------------------------------------------------
  !logical, parameter                   :: debug = .true.
  logical, parameter                  :: debug = .false.
  logical                              :: Endloop_q
  integer,allocatable                  :: Tab_Ib(:),Tab_Jb(:),Tab_Iq(:)
  integer                              :: Iq,ndim,Ib
  complex(kind=Rkind)                  :: C_I,C_J


  IF (debug) THEN
    write (out_unit, *) ''
    !flush(out_unit)   
  End IF
  
  ndim = size(Basis%tab_basis)-1 
  allocate(Tab_Ib(ndim),Tab_Jb(ndim),Tab_Iq(ndim))

  call Find_I(Tab_Ib,I,Basis)
  call Find_I(Tab_Jb,J,Basis)

  Call Init_tab_ind(Tab_Iq, Basis%NDindexq)
  Iq = 0
  S_IJ = ZERO

  DO
   Iq = Iq +1
   CALL increase_NDindex(Tab_Iq, Basis%NDindexq, Endloop_q)
   !print*,'Tab_Iq',Iq,Tab_Iq
   IF (Endloop_q) exit
   C_I = CONE
   C_J = CONE

   DO Ib = 1,ndim
      C_I = C_I*Basis%tab_basis(Ib)%d0bgw(Tab_Ib(Ib),Tab_Iq(Ib))
      C_J = C_J*Basis%tab_basis(Ib)%d0gb(Tab_Iq(Ib),Tab_Jb(Ib))
   End DO

   S_IJ = S_IJ+conjg(C_I)*C_J

  End DO


  IF (debug) THEN
   write (out_unit, *) 'I ,J ,S_IJ',I,J,S_IJ 
   flush(out_unit)   
   End IF

deallocate(Tab_Ib,tab_iq,Tab_Jb)

End SUBROUTINE


   end module Vp_m
