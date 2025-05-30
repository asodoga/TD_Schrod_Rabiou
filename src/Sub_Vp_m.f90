MODULE Sub_Vp_m
   USE QDUtil_m
   USE NDindex_m
   USE op_m
   USE Basis_m
   USE psi_m
   IMPLICIT NONE

   PRIVATE
   PUBLIC ::Vp_test_temp,Vp_step_psi


   contains

   SUBROUTINE psiTOLambda(L,psi)
   Implicit none
   complex(Kind=Rkind) ,intent(inout) :: L(:)
   TYPE(psi_t),intent(in)             :: psi
   
   !> Locals variables ------------------------------------------------
   
   integer                            :: Ib,nb,ndim,nsurf
   logical, parameter                 :: debug = .true.
   !logical, parameter                 :: debug = .false.
   real(kind=Rkind), allocatable      :: Qt(:),Pt(:),SQt(:)
   complex(kind=Rkind), allocatable   :: At(:)
   
   
   
  
   ndim = size(psi%Basis%tab_basis)-1
   nsurf = psi%Basis%tab_basis(ndim+1)%nb
    nb = psi%Basis%nb*nsurf
   allocate (Qt(ndim),Pt(ndim),At(ndim),SQt(ndim))

   IF (debug) THEN
    write (out_unit, *) 'psi in psiTOLambda:--------------------------------------------'
    Do Ib = 1,nb
       write (out_unit, *) Ib,psi%CVec(Ib)
    End Do
  End IF
   
   L(:) = CZERO
   call Get_Basis_Parameters(psi%Basis,Qt,SQt,At,Pt)
   L(1:nb) = psi%CVec(1:nb)
   L(nb+1:nb+ndim) = At(1:ndim)
   L(nb+ndim+1:nb+2*ndim) = Qt(1:ndim)
   L(nb+2*ndim+1:nb+3*ndim) = Pt(1:ndim)
   
    IF (debug) THEN
      write (out_unit, *) 'Lambda in psiTOLambda:----------------------------------------------------------'
      Do Ib = 1,nb+3*ndim
         write (out_unit, *) Ib,L(Ib)
      End Do
    End IF
   deallocate (Qt,Pt,At,SQt)

END SUBROUTINE

SUBROUTINE LambdaTOpsi(psi,L)
   Implicit none
   complex(Kind=Rkind) ,intent(in)    :: L(:)
   TYPE(psi_t),intent(inout)          :: psi
   
   !> Locals variables ------------------------------------------------
   
   integer                            :: Ib,nb,ndim,nsurf
   logical, parameter                 :: debug = .true.
   !logical, parameter                 :: debug = .false.
   real(kind=Rkind), allocatable      :: SQt(:), Qt(:),Pt(:)
   complex(kind=Rkind), allocatable   :: At(:)
   
   ndim = size(psi%Basis%tab_basis)-1
   nsurf = psi%Basis%tab_basis(ndim+1)%nb
   nb = psi%Basis%nb*nsurf
   allocate (SQt(ndim),Qt(ndim),Pt(ndim),At(ndim))

   IF (debug) THEN
      write (out_unit, *) 'Lambda in LambdaTOpsi:-------------------------------------------------'
      Do Ib = 1,nb+3*ndim
         write (out_unit, *) Ib,L(Ib)
      End Do
    End IF
   psi%CVec(:) = CZERO

    psi%CVec(1:nb)= L(1:nb) 
    At(1:ndim) = L(nb+1:nb+ndim)
    Qt(1:ndim) = L(nb+ndim+1:nb+2*ndim)
    Pt(1:ndim) = L(nb+2*ndim+1:nb+3*ndim)

   Do Ib = 1,ndim
       psi%Basis%tab_basis(Ib)%Q0 =Qt(Ib) 
       psi%Basis%tab_basis(Ib)%Imp_k =Pt(Ib)  
       psi%Basis%tab_basis(Ib)%alpha =At(Ib) 
   End DO

   SQt(:) = ONE
   call Construct_Hagedorn_Variational_Basis(psi%Basis,Qt,SQt,At,Pt)
   deallocate (SQt,Qt,Pt,At)



   IF (debug) THEN
      write (out_unit, *) 'psi in LambdaTOpsi:-----------------------------------------------------'
      Do Ib =1, nb 
         write (out_unit, *) Ib,psi%CVec(Ib)
      End Do
      flush(out_unit)
    End IF  

END SUBROUTINE




 SUBROUTINE Vp_test_temp(psi)
    TYPE(psi_t), intent(inout)          :: psi

    !> Locals variables ------------------------------------------------
    TYPE(psi_t)                         :: psif
    complex(Kind=Rkind) ,allocatable    :: L(:)
     complex(Kind=Rkind),allocatable    :: Mat(:,:),V(:)
     integer,allocatable                :: Tab_Ip(:)
    integer                             :: Ib,Jb,nb,ndim,nsurf

    nb = psi%Basis%nb
    ndim = size(psi%Basis%tab_basis)-1
    nsurf = psi%Basis%tab_basis(ndim+1)%nb
    allocate (L(nb*nsurf+3*ndim),Mat(nb*nsurf+3*ndim,nb*nsurf+3*ndim))
     call init_psi(psif, psi%Basis, cplx=.true., grid=.true.)
  
    print*,'----------------Debut du test sur LambdaTOpsi psiTOLambda----------'
    
   ! call psiTOLambda(L,psi)
   ! call LambdaTOpsi(psif,L)
   ! print*,'nb',nb,nsurf

    !  Do Ib =1,3
    !     Do Jb =1,3
    !     call  Calc_Saqp(Mat(Ib,Jb),psi,Ib,Jb)
    !     print*,'Ib,Jb,S ',Mat(Ib,Jb)
    !   END Do
    ! END Do
    
     !call Write_VecMat(Mat, out_unit, 3,  info='Mat')
     call Calc_GlobalOverlap_S(Mat,psi)
  ! Do Ib =1,nb*nsurf+3*ndim
  !  Do Jb =1,nb*nsurf+3*ndim
  !     print*,'Ib,Jb,S ',Ib,Jb,Mat(Ib,Jb)
  !  END Do
  !END Do

  DO ib =1,nb*nsurf+3*ndim
     print*,'Ib,S ',Ib,Mat(ib,:)
  End Do
      !call Write_VecMat(Mat, out_unit, 10,  info='Mat')
  ! call Tab(Tab_Ip,23-nb*nsurf,psi%Basis)
      !call Calc_V(V,psi)
      !print*,'V',V(:)
      !deallocate(Mat,V,L)
     print*,'---------------- Fin du test sur LambdaTOpsi psiTOLambda------------'
End SUBROUTINE



SUBROUTINE Find_I(Tab_Ib,I,Basis)
   implicit none
   integer,intent(in)                   :: I
   TYPE(Basis_t),intent(in)             :: Basis
   integer,intent(inout)                :: Tab_Ib(:)

   !> Locals variables ------------------------------------------------
   TYPE(NDindex_t)                       :: NDindex
   logical                               :: Endloop
   integer, allocatable                  :: NDend(:)
   !logical, parameter                   :: debug = .true.
   logical, parameter                    :: debug = .false.
   integer                               :: Ib,ndim


   IF (debug) THEN
     write (out_unit, *) 'Tab_Ib :-------------------------------------------------'
   End IF

   ndim = size(Basis%tab_basis)
   allocate (NDend(ndim))

   DO Ib = 1, ndim
     NDend(Ib) = Basis%tab_basis(Ib)%nb
   END DO

   call Init_NDindex(NDindex, NDend, ndim)
   call Init_tab_ind(Tab_Ib,NDindex)
   DO Ib = 1,I
    CALL increase_NDindex(Tab_Ib, NDindex, Endloop)
   ! write (out_unit, *) 'Ib , Tab_Ib',ib,Tab_Ib
   End DO

  call dealloc_NDindex(NDindex)
  deallocate(NDend)
   IF (debug) THEN
     write (out_unit, *) 'I , Tab_Ib',I,Tab_Ib
     write (out_unit, *) 'Tab_Ib :-------------------------------------------------'
   End IF
End SUBROUTINE



SUBROUTINE Calc_SIJ(S_IJ,I,J,Basis)
  implicit none
  
  TYPE(Basis_t),intent(in)             :: Basis
  complex(kind=Rkind),intent(inout)    :: S_IJ
  integer ,intent(in)                  :: I,J
  
  !> Locals variables ------------------------------------------------
  !logical, parameter                   :: debug = .true.
  logical, parameter                   :: debug = .false.
  logical                              :: Endloop_q
  integer,allocatable                  :: Tab_Ib(:),Tab_Jb(:),Tab_Iq(:)
  integer                              :: Iq,ndim,Ib
  complex(kind=Rkind)                  :: C_I,C_J
  real(kind=Rkind)                     :: delta
 

  IF (debug) THEN
    write (out_unit, *) ''
    !flush(out_unit)   
  End IF
  
  ndim = size(Basis%tab_basis)
  allocate(Tab_Ib(ndim),Tab_Jb(ndim),Tab_Iq(ndim-1))

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

   DO Ib = 1,ndim-1
      C_I = C_I*Basis%tab_basis(Ib)%d0bgw(Tab_Ib(Ib),Tab_Iq(Ib))
      C_J = C_J*Basis%tab_basis(Ib)%d0gb(Tab_Iq(Ib),Tab_Jb(Ib))
   End DO

   S_IJ = S_IJ+conjg(C_I)*C_J

  End DO
  call Kronecker_delta(Tab_Ib(ndim),Tab_Jb(ndim),delta)
  S_IJ = S_IJ*delta

  IF (debug) THEN
   write (out_unit, *) 'I ,J ,S_IJ',I,J,S_IJ 
   flush(out_unit)   
   End IF

deallocate(Tab_Ib,Tab_iq,Tab_Jb)

End SUBROUTINE




SUBROUTINE Calc_Overlap_S(Mat,Basis)

  implicit none
  complex(kind=Rkind) ,intent(inout) ,allocatable :: Mat(:,:)
  TYPE(Basis_t),intent(in)                        :: Basis
  integer                                         :: I1,I2,nb,nsurf,ndim

  ndim= size(Basis%tab_basis)
  nsurf = Basis%tab_basis(ndim)%nb
  nb = Basis%nb*nsurf
  If(allocated(Mat)) deallocate(Mat)
  allocate(Mat(nb,nb))

  DO I1= 1,nb
    DO I2= 1,nb
      call Calc_SIJ(Mat(I1,I2),I1,I2,Basis)
    End Do
  End Do

  !call Write_VecMat(Mat, out_unit, 9,  info='Mat')

End SUBROUTINE


SUBROUTINE Calc_dapsi(dapsi, psi, Ib)
   USE  QDUtil_m
   TYPE(psi_t), intent(in) ,target         :: psi
   TYPE(psi_t), intent(inout),target      :: dapsi
   integer,intent(in)                     :: Ib

   !> Locals variables ------------------------------------------------
   !logical, parameter                    :: debug = .true.
   logical, parameter                     :: debug = .false.
   complex(kind=Rkind), pointer           :: psi_ggb(:, :, :)
   complex(kind=Rkind), pointer           :: dapsi_ggb(:, :, :)
   complex(kind=Rkind), pointer           :: dagg(:, :)
   Integer, allocatable                   :: Iq1(:), Iq2(:), Iq3(:)
   !Integer, allocatable                   :: Ib1(:), Ib2(:), Ib3(:)
   integer                                :: i1,i3

   
   IF (debug) THEN
      write(out_unit,*) 'BEGINNING dapsi'
      flush (out_unit)
   END IF
         dapsi%CVec(:) =CZERO
         Call Calc_index(Iq1=Iq1, Iq2=Iq2,Iq3=Iq3, Basis=psi%Basis)

          psi_ggb(1:Iq1(Ib),1:Iq2(Ib),1:Iq3(Ib)) => psi%CVec
          dapsi_ggb(1:Iq1(Ib),1:Iq2(Ib),1:Iq3(Ib)) => dapsi%CVec
          dagg(1:Iq2(Ib),1:Iq2(Ib)) => psi%Basis%tab_basis(Ib)%dagg(:, :)


   DO i3 = 1, ubound(psi_ggb, dim=3)
     DO i1 = 1, ubound(psi_ggb, dim=1)
       dapsi_ggb(i1, :, i3) = dapsi_ggb(i1, :, i3)+matmul(dagg(:,:),psi_ggb(i1,:,i3))
     END DO
   END DO

   deallocate(Iq1, Iq2, Iq3)
       
   IF (debug) THEN
              write(out_unit,*) 'END dapsi'
      flush (out_unit)
   END IF
END SUBROUTINE 





SUBROUTINE Calc_dqpsi(dqpsi, psi, Ib)
   USE  QDUtil_m
   TYPE(psi_t), intent(in) ,target         :: psi
   TYPE(psi_t), intent(inout),target       :: dqpsi
   integer,intent(in)                      :: Ib


   !> Locals variables ------------------------------------------------
   !logical, parameter                    :: debug = .true.
   logical, parameter                     :: debug = .false.
   complex(kind=Rkind), pointer           :: psi_ggb(:, :, :)
   complex(kind=Rkind), pointer           :: dqpsi_ggb(:, :, :)
   complex(kind=Rkind), pointer           :: dqgg(:, :)
   Integer, allocatable                   :: Iq1(:), Iq2(:), Iq3(:)
   !Integer, allocatable                   :: Ib1(:), Ib2(:), Ib3(:)
   integer                                :: i1,i3

   
   IF (debug) THEN
      write(out_unit,*) 'BEGINNING dqpsi'
      flush (out_unit)
   END IF
         dqpsi%CVec(:) =CZERO
         Call Calc_index(Iq1=Iq1, Iq2=Iq2,Iq3=Iq3, Basis=psi%Basis)

          psi_ggb(1:Iq1(Ib),1:Iq2(Ib),1:Iq3(Ib)) => psi%CVec
          dqpsi_ggb(1:Iq1(Ib),1:Iq2(Ib),1:Iq3(Ib)) => dqpsi%CVec
          dqgg(1:Iq2(Ib),1:Iq2(Ib)) => psi%Basis%tab_basis(Ib)%dq0gg(:, :)


   DO i3 = 1, ubound(psi_ggb, dim=3)
     DO i1 = 1, ubound(psi_ggb, dim=1)
       dqpsi_ggb(i1, :, i3) = dqpsi_ggb(i1, :, i3)+matmul(dqgg(:,:),psi_ggb(i1,:,i3))
     END DO
   END DO

    deallocate(Iq1, Iq2, Iq3)
       
   IF (debug) THEN
              write(out_unit,*) 'END dqpsi'
      flush (out_unit)
   END IF
END SUBROUTINE 


SUBROUTINE Calc_dppsi(dppsi, psi, Ib)
   USE  QDUtil_m
   TYPE(psi_t), intent(in) ,target         :: psi
   TYPE(psi_t), intent(inout),target       :: dppsi
   integer,intent(in)                      :: Ib


   !> Locals variables ------------------------------------------------
   !logical, parameter                    :: debug = .true.
   logical, parameter                     :: debug = .false.
   complex(kind=Rkind), pointer           :: psi_ggb(:, :, :)
   complex(kind=Rkind), pointer           :: dppsi_ggb(:, :, :)
   complex(kind=Rkind), pointer           :: dpgg(:, :)
   Integer, allocatable                   :: Iq1(:), Iq2(:), Iq3(:)
   !Integer, allocatable                   :: Ib1(:), Ib2(:), Ib3(:)
   integer                                :: i1,i3

   
   IF (debug) THEN
      write(out_unit,*) 'BEGINNING dapsi'
      flush (out_unit)
   END IF
         dppsi%CVec(:) =CZERO
         Call Calc_index(Iq1=Iq1, Iq2=Iq2,Iq3=Iq3, Basis=psi%Basis)

          psi_ggb(1:Iq1(Ib),1:Iq2(Ib),1:Iq3(Ib)) => psi%CVec
          dppsi_ggb(1:Iq1(Ib),1:Iq2(Ib),1:Iq3(Ib)) => dppsi%CVec
          dpgg(1:Iq2(Ib),1:Iq2(Ib)) => psi%Basis%tab_basis(Ib)%dp0gg(:, :)


   DO i3 = 1, ubound(psi_ggb, dim=3)
     DO i1 = 1, ubound(psi_ggb, dim=1)
       dppsi_ggb(i1, :, i3) = dppsi_ggb(i1, :, i3)+matmul(dpgg(:,:),psi_ggb(i1,:,i3))
     END DO
   END DO

    deallocate(Iq1, Iq2, Iq3)
       
   IF (debug) THEN
              write(out_unit,*) 'END dapsi'
      flush (out_unit)
   END IF
END SUBROUTINE

SUBROUTINE Calc_daqppsi_IJ(dIpsi,dJpsi,psi0,I,J)
  implicit none
  TYPE(psi_t), intent(in) ,target      :: psi0
  TYPE(psi_t), intent(inout) ,target   :: dIpsi,dJpsi
  integer ,intent(in)                  :: I,J

  
  !> Locals variables ------------------------------------------------
  !logical, parameter                  :: debug = .true.
  logical, parameter                   ::  debug = .false.
  TYPE(psi_t)                          :: psi
  integer,allocatable                  :: Tab_Ip(:),Tab_Jp(:)
  integer                              ::   ndim
  
  IF (debug) THEN
    write (out_unit, *) ''
    !flush(out_unit)   
  End IF
  
  ndim =  size(psi0%Basis%tab_basis)
   call init_psi(psi, psi0%Basis,cplx=.true.,grid=.true.)
  IF (psi0%Grid) then
    psi%CVec(:) = psi0%CVec(:)
  ELSE
    call BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
  END IF

  call Tab(Tab_Ip,I,psi0%Basis)
  call Tab(Tab_Jp,J,psi0%Basis)
  dIpsi%CVec = CZERO
  dJpsi%CVec = CZERO

  If(Tab_Ip(2)==1) Then
    call Calc_dapsi(dIpsi, psi, Tab_Ip(1)) 
  Else If(Tab_Ip(2)==2) Then
   call Calc_dqpsi(dIpsi, psi, Tab_Ip(1))
  Else If(Tab_Ip(2)==3) Then
   call Calc_dppsi(dIpsi, psi, Tab_Ip(1))
  End If

 If(Tab_Jp(2)==1) Then
   call Calc_dapsi(dJpsi, psi, Tab_Jp(1)) 
 Else If(Tab_Jp(2)==2) Then
   call Calc_dqpsi(dJpsi, psi, Tab_Jp(1))
 Else If(Tab_Jp(2)==3) Then
   call Calc_dppsi(dJpsi, psi, Tab_Jp(1))
 End If

  IF (debug) THEN
   flush(out_unit)   
   End IF

deallocate(Tab_Ip,Tab_Jp)
call dealloc_psi(psi)

End SUBROUTINE


SUBROUTINE Calc_daqppsi_I(dIpsi,psi0,I)
  implicit none
  TYPE(psi_t), intent(in) ,target      :: psi0
  TYPE(psi_t), intent(inout) ,target   :: dIpsi
  integer ,intent(in)                  :: I

  
  !> Locals variables ------------------------------------------------
  !logical, parameter                  :: debug = .true.
  logical, parameter                   ::  debug = .false.
  TYPE(psi_t)                          :: psi
  integer,allocatable                  :: Tab_Ip(:)
  integer                              ::   ndim
  
  IF (debug) THEN
    write (out_unit, *) ''
    !flush(out_unit)   
  End IF
  
  ndim =  size(psi0%Basis%tab_basis)
   call init_psi(psi, psi0%Basis,cplx=.true.,grid=.true.)
  IF (psi0%Grid) then
    psi%CVec(:) = psi0%CVec(:)
  ELSE
    call BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
  END IF

  call Tab(Tab_Ip,I,psi0%Basis)
  dIpsi%CVec = CZERO

  If(Tab_Ip(2)==1) Then
    call Calc_dapsi(dIpsi, psi, Tab_Ip(1)) 
  Else If(Tab_Ip(2)==2) Then
   call Calc_dqpsi(dIpsi, psi, Tab_Ip(1))
  Else If(Tab_Ip(2)==3) Then
   call Calc_dppsi(dIpsi, psi, Tab_Ip(1))
  End If

  IF (debug) THEN
   flush(out_unit)   
   End IF

deallocate(Tab_Ip)
call dealloc_psi(psi)

End SUBROUTINE

SUBROUTINE Calc_Saqp(S_aqp,psi0,I,J)
  implicit none
  TYPE(psi_t), intent(in) ,target      :: psi0
  complex(kind=Rkind),intent(inout)    :: S_aqp
  integer ,intent(in)                  :: I,J

  
  !> Locals variables ------------------------------------------------
  !logical, parameter                   :: debug = .true.
  logical, parameter                   ::  debug = .false.
  TYPE(psi_t),target                   ::  psi,dJpsi,dIpsi
  complex(kind=Rkind), pointer         ::  psi_gb(:, :),psi_gbI(:, :),psi_gbJ(:, :)
  integer,allocatable                  ::  Tab_Iq(:)
  logical                              ::  Endloop_q
  integer                              ::  Iq,ndim,nsurf,Isurf,ib,nq
 real(kind=Rkind)                      ::  w
  complex(kind=Rkind) ,allocatable     ::  S_aqpEl(:)
  real(kind=Rkind) ,allocatable        ::  N (:)
 

  IF (debug) THEN
    write (out_unit, *) ''
    !flush(out_unit)   
  End IF
  
  ndim =  size(psi0%Basis%tab_basis)
  nq = psi0%Basis%nq
  nsurf = psi0%Basis%tab_basis(ndim)%nb
  allocate(Tab_Iq(ndim-1),S_aqpEl(nsurf),N(nsurf))
  

  call init_psi(psi, psi0%Basis,     cplx= .true., grid=.true.)
  call init_psi(dIpsi, psi0%Basis,   cplx= .true., grid=.true.)
  call init_psi(dJpsi, psi0%Basis,   cplx= .true., grid=.true.)

  IF (psi0%Grid) then
    psi%CVec(:) = psi0%CVec(:)
  ELSE
    call BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
  END IF
   call Calc_daqppsi_IJ(dIpsi,dJpsi,psi0,I,J)
   psi_gb(1:nq,1:nsurf) => psi%CVec
   psi_gbI(1:nq,1:nsurf) => dIpsi%CVec
   psi_gbJ(1:nq,1:nsurf) => dJpsi%CVec

   Iq = 0
   S_aqpEl(:) = CZERO
   N(:) = ZERO

   DO Isurf = 1,nsurf
   Call Init_tab_ind(Tab_Iq, psi0%Basis%NDindexq)
   DO
   Iq = Iq +1
   CALL increase_NDindex(Tab_Iq, psi0%Basis%NDindexq, Endloop_q)
   IF (Endloop_q) exit
   w = ONE
   DO Ib = 1, ndim-1
   w = w*psi%Basis%tab_basis(Ib)%w(Tab_iq(Ib))
   END DO

   N(Isurf) = N(Isurf) + conjg(psi_gb(Iq, Isurf))*psi_gb(Iq, Isurf)*w
   S_aqpEl(Isurf) = S_aqpEl(Isurf) + conjg(psi_gbI(Iq, Isurf))*psi_gbJ(Iq, Isurf)*w

 End DO
 End DO
S_aqp = CZERO
S_aqp = sum(S_aqpEl)/(sum(N)**2)


  IF (debug) THEN
   write (out_unit, *) 'I ,J ,S_aqp',I,J,S_aqp
   flush(out_unit)   
   End IF
   
 deallocate(Tab_Iq,S_aqpEl,N)
 call dealloc_psi(psi)
 call dealloc_psi(dJpsi)
 call dealloc_psi(dIpsi)

End SUBROUTINE




SUBROUTINE Calc_SIaqp(S_Iaqp,psi0,I,J)
  implicit none
  TYPE(psi_t), intent(in) ,target      :: psi0
  complex(kind=Rkind),intent(inout)    :: S_Iaqp
  integer ,intent(in)                  :: I,J

  
  !> Locals variables ------------------------------------------------
  !logical, parameter                   :: debug = .true.
  logical, parameter                   ::  debug = .false.
  TYPE(psi_t),target                   ::  psi,dJpsi
  complex(kind=Rkind), pointer         ::  psi_gb(:, :),psi_gbJ(:, :)
  integer,allocatable                  ::  Tab_Iq(:),Tab_Ib(:)
  logical                              ::  Endloop_q
  integer                              ::  Iq,ndim,nsurf,Isurf,ib,nq
  real(kind=Rkind)                     ::  w,N
  complex(kind=Rkind)                  ::  C_I  

  IF (debug) THEN
    write (out_unit, *) ''
    !flush(out_unit)   
  End IF
  
  ndim =  size(psi0%Basis%tab_basis)
  nq = psi0%Basis%nq
  nsurf = psi0%Basis%tab_basis(ndim)%nb
  allocate(Tab_Iq(ndim-1),Tab_Ib(ndim))
  

  call init_psi(psi, psi0%Basis,     cplx= .true., grid=.true.)
  call init_psi(dJpsi, psi0%Basis,   cplx= .true., grid=.true.)
  call Find_I(Tab_Ib,I,psi0%Basis)

  IF (psi0%Grid) then
    psi%CVec(:) = psi0%CVec(:)
  ELSE
    call BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
  END IF
   call Calc_daqppsi_I(dJpsi,psi0,J)
   psi_gb(1:nq,1:nsurf) => psi%CVec
   psi_gbJ(1:nq,1:nsurf) => dJpsi%CVec

   Call Init_tab_ind(Tab_Iq, psi0%Basis%NDindexq)
   Iq = 0
   S_Iaqp = CZERO
   N = ZERO
   Isurf = Tab_Ib(ndim)
   DO
   Iq = Iq +1
   CALL increase_NDindex(Tab_Iq, psi0%Basis%NDindexq, Endloop_q)
   IF (Endloop_q) exit
   w = ONE
   C_I = ONE
   DO Ib = 1, ndim- 1
   w = w*psi%Basis%tab_basis(Ib)%w(Tab_Iq(Ib))
   C_I = C_I*psi%Basis%tab_basis(Ib)%d0gb(Tab_Iq(Ib),Tab_Ib(Ib))
   END DO

   N = N + conjg(psi_gb(Iq, Isurf))*psi_gb(Iq, Isurf)*w
   S_Iaqp = S_Iaqp + conjg(C_I)*psi_gbJ(Iq, Isurf)*w

 End DO


call dealloc_psi(psi)
call dealloc_psi(dJpsi)
deallocate(Tab_Iq,Tab_Ib)

  IF (debug) THEN
   write (out_unit, *) 'I ,J ,S_Iaqp',I,J,S_Iaqp
   flush(out_unit)   
   End IF
 

End SUBROUTINE



SUBROUTINE Calc_SaqpJ(S_aqpJ,psi0,I,J)
  implicit none
  TYPE(psi_t), intent(in) ,target      :: psi0
  complex(kind=Rkind),intent(inout)    :: S_aqpJ
  integer ,intent(in)                  :: I,J

  
  !> Locals variables ------------------------------------------------
  !logical, parameter                   :: debug = .true.
  logical, parameter                   ::  debug = .false.
  TYPE(psi_t),target                   ::  psi,dIpsi
  complex(kind=Rkind), pointer         ::  psi_gb(:, :),psi_gbI(:, :)
  integer,allocatable                  ::  Tab_Iq(:),Tab_Jb(:)
  logical                              ::  Endloop_q
  integer                              ::  Iq,ndim,nsurf,Isurf,ib,nq
  real(kind=Rkind)                     ::  w,N
  complex(kind=Rkind)                  ::  C_J  

  IF (debug) THEN
    write (out_unit, *) ''
    !flush(out_unit)   
  End IF
  
  ndim =  size(psi0%Basis%tab_basis)
  nq = psi0%Basis%nq
  nsurf = psi0%Basis%tab_basis(ndim)%nb
  allocate(Tab_Iq(ndim-1),Tab_Jb(ndim))
  

  call init_psi(psi, psi0%Basis,     cplx= .true., grid=.true.)
  call init_psi(dIpsi, psi0%Basis,   cplx= .true., grid=.true.)
  call Find_I(Tab_Jb,J,psi0%Basis)

  IF (psi0%Grid) then
    psi%CVec(:) = psi0%CVec(:)
  ELSE
    call BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
  END IF
   call Calc_daqppsi_I(dIpsi,psi0,I)
   psi_gb(1:nq,1:nsurf) => psi%CVec
   psi_gbI(1:nq,1:nsurf) => dIpsi%CVec

   Call Init_tab_ind(Tab_Iq, psi0%Basis%NDindexq)
   Iq = 0
   S_aqpJ = CZERO
   N = ZERO
   Isurf = Tab_Jb(ndim)
   DO
   Iq = Iq+1
   CALL increase_NDindex(Tab_Iq, psi0%Basis%NDindexq, Endloop_q)
   IF (Endloop_q) exit
   w = ONE
   C_J = ONE
   DO Ib = 1, ndim-1
   w = w*psi%Basis%tab_basis(Ib)%w(Tab_Iq(Ib))
   C_J = C_J*psi%Basis%tab_basis(Ib)%d0gb(Tab_Iq(Ib),Tab_Jb(Ib))
   END DO

   N = N + conjg(psi_gb(Iq, Isurf))*psi_gb(Iq, Isurf)*w
   S_aqpJ = S_aqpJ + conjg(psi_gbI(Iq, Isurf))*C_J*w

 End DO

  deallocate(Tab_Iq,Tab_Jb)
  call dealloc_psi(psi)
  call dealloc_psi(dIpsi)


  IF (debug) THEN
   write (out_unit, *) 'I ,J ,S_aqpJ',I,J,S_aqpJ
   flush(out_unit)   
   End IF

End SUBROUTINE


SUBROUTINE Kronecker_delta(I,J,delta)
  implicit none
 integer,intent(in)              :: I,J
 real(kind=Rkind),intent(inout)  :: delta

 !> Locals variables ------------------------------------------------
 !logical, parameter                    :: debug = .true.
 logical, parameter                     :: debug = .false.

   IF (debug) THEN
           write(out_unit,*) 'Begin Kronecker'
   flush (out_unit)
  END IF

  If (I==J) Then
    delta = ONE
  Else
   delta = ZERO
  End If    

 ! write(out_unit,*)  'delta',delta

 IF (debug) THEN
   write(out_unit,*) 'END kronecker delta'
   flush (out_unit)
 END IF

End SUBROUTINE



SUBROUTINE Tab(Tab_Ip,I,Basis)
  integer ,intent(inout),allocatable    :: Tab_Ip(:)
  integer  ,intent(in)                  :: I
  TYPE(Basis_t),intent(in)              :: Basis
  integer                               :: ndim,II
  TYPE(NDindex_t)                       :: NDindex
  logical                               :: Endloop
  integer, allocatable                  :: NDend(:)

  ndim = size(Basis%tab_basis)-1
  If(allocated(Tab_Ip)) deallocate(Tab_Ip)
  allocate (NDend(2),Tab_Ip(2))
  NDend(1) = ndim
  NDend(2) = 3

  call Init_NDindex(NDindex, NDend, 2)
  call Init_tab_ind(Tab_Ip,NDindex)

  Do II = 1,I
  CALL increase_NDindex(Tab_Ip, NDindex, Endloop)
  !write(out_unit,*) 'I,Tab_Ip',II,Tab_Ip
  End DO

  call dealloc_NDindex(NDindex)
  deallocate(Ndend)
  flush(out_unit)
  !write(out_unit,*) 'Tab_Ip',Tab_Ip

  End SUBROUTINE




  SUBROUTINE Calc_SIHpsi(S_I,psi0,I)
  implicit none
  TYPE(psi_t), intent(in) ,target      :: psi0
  complex(kind=Rkind),intent(inout)    :: S_I
  integer ,intent(in)                  :: I

  
  !> Locals variables ------------------------------------------------
  !logical, parameter                   :: debug = .true.
  logical, parameter                   ::  debug = .false.
  TYPE(Op_t)                           :: H
  TYPE(psi_t),target                   ::  psi,Hpsi,Hpsi_g
  complex(kind=Rkind), pointer         ::  psi_gb(:, :),Hpsi_gb(:, :)
  integer,allocatable                  ::  Tab_Iq(:),Tab_Ib(:)
  logical                              ::  Endloop_q
  integer                              ::  Iq,ndim,nsurf,Isurf,ib,nq
  real(kind=Rkind)                     ::  w,N
  complex(kind=Rkind)                  ::  C_I  

  IF (debug) THEN
    write (out_unit, *) ''
    !flush(out_unit)   
  End IF
  
  ndim =  size(psi0%Basis%tab_basis)
  nq = psi0%Basis%nq
  nsurf = psi0%Basis%tab_basis(ndim)%nb
  allocate(Tab_Iq(ndim-1),Tab_Ib(ndim)) 

  call init_psi(psi, psi0%Basis,     cplx= .true., grid=.true.)
  call init_psi(Hpsi, psi0%Basis,   cplx= .true., grid=.false.)
  call init_psi(Hpsi_g, psi0%Basis,   cplx= .true., grid=.true.)
  call Find_I(Tab_Ib,I,psi0%Basis)
  call calc_OpPsi(H, psi0, Hpsi)
  call BasisTOGrid_nD_cplx(Hpsi_g%CVec, Hpsi%CVec, psi0%Basis)

  IF (psi0%Grid) then
    psi%CVec(:) = psi0%CVec(:)
  ELSE
    call BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
  END IF
  
   psi_gb(1:nq,1:nsurf) => psi%CVec
   Hpsi_gb(1:nq,1:nsurf) => Hpsi_g%CVec

   Call Init_tab_ind(Tab_Iq, psi0%Basis%NDindexq)
   Iq = 0
   S_I = CZERO
   N = ZERO
   Isurf = Tab_Ib(ndim)
   DO
   Iq = Iq +1
   CALL increase_NDindex(Tab_Iq, psi0%Basis%NDindexq, Endloop_q)
   IF (Endloop_q) exit
   w = ONE
   C_I = ONE
   DO Ib = 1, ndim- 1
   w = w*psi%Basis%tab_basis(Ib)%w(Tab_Iq(Ib))
   C_I = C_I*psi%Basis%tab_basis(Ib)%d0gb(Tab_Iq(Ib),Tab_Ib(Ib))
   END DO

   N = N + conjg(psi_gb(Iq, Isurf))*psi_gb(Iq, Isurf)*w
   S_I = S_I + conjg(C_I)*Hpsi_gb(Iq, Isurf)*w

 End DO

  IF (debug) THEN
   write (out_unit, *) 'I ,S_I',I,S_I
   flush(out_unit)   
   End IF
  call dealloc_psi(psi)
  call dealloc_psi(Hpsi)
  call dealloc_psi(Hpsi_g)
 deallocate(Tab_Iq,Tab_Ib)

End SUBROUTINE



 SUBROUTINE Calc_SaqpHpsi(S_I,psi0,I)
 implicit none
 TYPE(psi_t), intent(in) ,target      :: psi0
 complex(kind=Rkind),intent(inout)    :: S_I
 integer ,intent(in)                  :: I
 
 !> Locals variables ------------------------------------------------
 !logical, parameter                   :: debug = .true.
 logical, parameter                   ::  debug = .false.
 TYPE(Op_t)                           :: H
 TYPE(psi_t),target                   ::  psi,Hpsi,Hpsi_g,dpsi
 complex(kind=Rkind), pointer         ::  psi_gb(:, :),Hpsi_gb(:, :),dpsi_gb(:,:)
 integer,allocatable                  ::  Tab_Iq(:),Tab_Ip(:)
 logical                              ::  Endloop_q
 integer                              ::  Iq,ndim,nsurf,Isurf,ib,nq
  complex(kind=Rkind) ,allocatable    ::  S_IEl(:)
 real(kind=Rkind) ,allocatable        ::  N (:)
 real(kind=Rkind)                     ::  w

 IF (debug) THEN
   write (out_unit, *) ''
   !flush(out_unit)   
 End IF
 
 ndim =  size(psi0%Basis%tab_basis)
 nq = psi0%Basis%nq
 nsurf = psi0%Basis%tab_basis(ndim)%nb
 allocate(Tab_Iq(ndim-1),S_IEl(nsurf),N(nsurf)) 
 call init_psi(psi, psi0%Basis,   cplx= .true., grid=.true.)
 call init_psi(dpsi, psi0%Basis,  cplx= .true., grid=.true.)
 call init_psi(Hpsi, psi0%Basis,  cplx= .true., grid=.false.)
 call init_psi(Hpsi_g, psi0%Basis,cplx= .true., grid=.true.)
 call Tab(Tab_Ip,I,psi0%Basis) 
 call calc_OpPsi(H, psi0, Hpsi)
 call BasisTOGrid_nD_cplx(Hpsi_g%CVec, Hpsi%CVec, psi0%Basis)
 call Calc_daqppsi_I(dpsi,psi0,I)

 IF (psi0%Grid) then
   psi%CVec(:) = psi0%CVec(:)
 ELSE
   call BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
 END IF
 
  psi_gb(1:nq,1:nsurf) =>  psi%CVec
  Hpsi_gb(1:nq,1:nsurf) => Hpsi_g%CVec
  dpsi_gb(1:nq,1:nsurf) => dpsi%CVec

  N(:) = ZERO
  S_IEl(:) = CZERO
  S_I = CZERO
 DO Isurf = 1,nsurf
  Tab_Iq = 0
  Call Init_tab_ind(Tab_Iq, psi0%Basis%NDindexq)
  Iq = 0
  DO
  Iq = Iq +1
  CALL increase_NDindex(Tab_Iq, psi0%Basis%NDindexq, Endloop_q)
  IF (Endloop_q) exit
  w = ONE
  DO Ib = 1, ndim- 1
  w = w*psi%Basis%tab_basis(Ib)%w(Tab_Iq(Ib))
  END DO
  N(Isurf) = N(Isurf) + conjg(psi_gb(Iq, Isurf))*psi_gb(Iq, Isurf)*w
  S_IEl(Isurf) = S_IEl(Isurf) + conjg(dpsi_gb(Iq, Isurf))*Hpsi_gb(Iq, Isurf)*w
End DO
End DO

S_I = sum(S_IEl)/(sum(N)**2)

call dealloc_psi(psi)
call dealloc_psi(Hpsi)
call dealloc_psi(Hpsi_g)
deallocate(Tab_Iq,Tab_Ip)
call dealloc_psi(dpsi)
deallocate(S_IEl,N)

 IF (debug) THEN
  write (out_unit, *) 'I ,S_I',I,S_I
  flush(out_unit)   
  End IF

End SUBROUTINE


SUBROUTINE Calc_V(V,psi0)

  TYPE(psi_t), intent(in) ,target                  :: psi0
  complex(kind=Rkind),intent(inout),allocatable    :: V(:)

  integer                                          ::  Ib,ndim,nsurf,nb,n


   ndim =  size(psi0%Basis%tab_basis)-1
   nsurf = psi0%Basis%tab_basis(ndim+1)%nb
   n =psi0%Basis%nb*nsurf
   nb =n +3*ndim

   If(allocated(V)) deallocate(V)
   allocate(V(nb))

   Do Ib = 1,nb
     If(Ib<= n) Then
       call Calc_SIHpsi(V(Ib),psi0,Ib)
     Else
       call Calc_SaqpHpsi(V(Ib),psi0,Ib-n)
        !write (out_unit, *) 'Ib , nb,nsurf,Ib-nb*nsurf',Ib,psi0%Basis%nb,nsurf,Ib-psi0%Basis%nb*nsurf
     End If   
   End Do
   V(:) = -EYE*V(:)

End SUBROUTINE


SUBROUTINE Calc_GlobalOverlap_S(Mat,psi0)

  implicit none
  complex(kind=Rkind) ,intent(inout) ,allocatable    :: Mat(:,:)
  TYPE(psi_t), intent(in) ,target                    :: psi0

  !> ---------------local variables--------------------------------------------
 ! logical, parameter                                 :: debug = .true.
 logical, parameter                                 ::  debug = .false.
  integer                                            :: I,J,nb,nsurf
  integer                                            :: ndim,n

 ndim =  size(psi0%Basis%tab_basis)-1
 nsurf = psi0%Basis%tab_basis(ndim+1)%nb
 n = psi0%Basis%nb*nsurf
 nb = n+3*ndim


  If(allocated(Mat)) deallocate(Mat)
  allocate(Mat(nb,nb))

  Mat(:,:) = CZERO

DO I=1,nb
   If(I<=n) Then
   DO J=1,nb
     If(J<= n) Then
       call Calc_SIJ(Mat(I,J),I,J,psi0%Basis)
     Else
      call Calc_SIaqp(Mat(I,J),psi0,I,J-n) 
     End If
   End DO
   Else
   DO J=1,nb
     If(J<=n) Then
       call Calc_SaqpJ(Mat(I,J),psi0,I-n,J)
      Else
       call Calc_Saqp(Mat(I,J),psi0,I-n,J-n)
        !write (out_unit, *) 'Ib , nb,nsurf,Ib-nb*nsurf',i,psi0%Basis%nb,nsurf,i-psi0%Basis%nb*nsurf
        !write (out_unit, *) 'Jb , nb,nsurf,Ib-nb*nsurf',j,psi0%Basis%nb,nsurf,j-psi0%Basis%nb*nsurf
     End If
   End DO
   END If
 End DO

  IF (debug) THEN
    call Write_VecMat(Mat, out_unit, 16,  info='Mat')
    flush(out_unit)   
  End IF
  

End SUBROUTINE




 SUBROUTINE Vp_step_psi(psi, psi_dt,delta_t)
    USE  QDUtil_m
    USE op_m
    USE psi_m
    TYPE(psi_t), intent(inout)       :: psi_dt
    TYPE(psi_t), intent(in)          :: psi
    real(kind=Rkind),intent(in)      :: delta_t

    !> ---------------local variables------------------------
   !logical, parameter                :: debug = .true.
   logical, parameter               ::  debug = .false.
     TYPE(psi_t)                     :: psi_temp
    complex(kind=Rkind),allocatable  :: Vec1(:), Vec2(:),Vec3(:),Vec4(:)
    complex(kind=Rkind),allocatable  :: CL_t(:),CL_tdt(:),CL_temp(:)
    integer                          :: ndim,nb,nsurf
    

    !  variables locales
    ndim = size(psi%Basis%tab_basis)-1
    nsurf = psi%Basis%tab_basis(ndim+1)%nb
    nb =psi%Basis%nb*nsurf+3*ndim 
    psi_dt%CVec = CZERO

    allocate(CL_t(nb),CL_tdt(nb),CL_temp(nb))
    call init_psi(psi_temp, psi%basis, cplx=.true., grid=.false.)

    call psiTOLambda(CL_tdt,psi) !> labda init
    call Runge_Kutta_Vp_Func_cplx(Vec1,psi)
    CL_temp(:) = CL_tdt(:) + (delta_t*HALF)*Vec1(:)
    call LambdaTOpsi(psi_temp,CL_temp)
    call Runge_Kutta_Vp_Func_cplx(Vec2,psi_temp)
    CL_temp(:) = CL_tdt(:) + (delta_t*HALF)*Vec2(:)
    call LambdaTOpsi(psi_temp,CL_temp)
    call Runge_Kutta_Vp_Func_cplx(Vec3,psi_temp)
    CL_temp(:) = CL_tdt(:) + delta_t*Vec3(:)
    call LambdaTOpsi(psi_temp,CL_temp)
    call Runge_Kutta_Vp_Func_cplx(Vec4,psi_temp)

    CL_tdt(:) =CL_tdt(:) + (delta_t*SIXTH)*(Vec1(:) + TWO*Vec2(:) + TWO*Vec3(:) + Vec4(:))
    call LambdaTOpsi(psi_dt,CL_tdt)

    IF (debug) THEN
      print*,'CL_tdt',CL_tdt
      flush(out_unit)   
    End IF

    call dealloc_psi(psi_temp) 
    deallocate(Vec1,Vec2,Vec3,Vec4) 
    deallocate(CL_t,CL_tdt,CL_temp)  
  
 END SUBROUTINE 


 SUBROUTINE Runge_Kutta_Vp_Func_cplx(CL,psi)
   USE  QDUtil_m
    TYPE(psi_t),intent(in)                          :: psi
    complex(kind=Rkind),allocatable ,intent(inout)  :: CL(:)

     !> ---------------local variables------------------------
     logical, parameter                    :: debug = .true.
     !logical, parameter                   ::  debug = .false.
     complex(kind=Rkind),allocatable       :: CMat(:,:),CVec(:),CMat_inv(:,:)
     integer                               :: ndim,nb,nsurf,nb1
    
     ndim = size(psi%Basis%tab_basis)-1
     nb1 = psi%Basis%nb
     nsurf = psi%Basis%tab_basis(ndim+1)%nb
     nb =nb1*nsurf+3*ndim 
     If(allocated(CL)) deallocate(CL)
     allocate(CL(nb))!,CMat_inv(nb*nsurf+3*ndim,nb*nsurf+3*ndim))
     call Calc_GlobalOverlap_S(CMat,psi)
     call Calc_V(CVec,psi)
     !CMat_inv = inv_OF_Mat_TO(Cmat)
     !call Write_VecMat(CMat_inv, out_unit, 8,  info='CMat_inv')
     CL = LinearSys_Solve(CMat,CVec)
     !CL = matmul(CMat_inv,CVec)
     !print*,'CL',CL
     deallocate(CMat,CVec)

END SUBROUTINE
  

END MODULE 
