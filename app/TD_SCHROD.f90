PROGRAM TD_SCHROD
  USE NumParameters_m
  USE GWPnD_m
  USE Basis_m
  USE psi_m
  Use Molec_m
  !USE psi0_m
  USE Propa_m
  IMPLICIT NONE
  TYPE (Basis_t), target      :: Basis
  TYPE (psi_t)                ::Hpsi
  TYPE (psi_t)                :: psi0,psif,psi,B,G
  TYPE(propa_t)               :: propa
  real(Kind = Rk)             :: E
  !complex(kind=Rk),allocatable            :: G1,B1
  REAL(kind=Rk),allocatable   :: Mat_V(:,:)
    real(Kind = Rk) ,ALLOCATABLE        ::Q(:,:)
  !====================================================================
  ! for QML
  integer :: ndim,nsurf,option
  logical :: adiabatic
  character (len=16)                  :: pot_name

  ndim      = 0
  nsurf     = 0
  pot_name  = 'read_model'
  adiabatic = .FALSE.
  option    = 1
  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)
   write(out_unitp,*)'ndim,nsurf',ndim,nsurf
  write(out_unitp,*) 'pot_name',pot_name

  allocate(Mat_V(nsurf,nsurf))
  call sub_pot(Mat_V,[1._Rk])
  write(out_unitp,*) 'pot_name',pot_name
  !====================================================================
  ! read some informations (basis set/grid) : numbers of basis functions, grid points ...
  ! the basis/grid informations have to be put in a module
  CALL Read_Basis(Basis,nio=in_unitp)
!====================================================================
!print*,"Basis is allocated",Basis_IS_allocated(Basis)
write(out_unitp,*) 'Initialization of a complex psi'
CALL init_psi(psi0,   Basis,    cplx=.TRUE.   ,grid =.true.)
  CALL init_psi(psi,   Basis,    cplx=.TRUE.   ,grid =.false.)
CALL init_psi(psif,   Basis,    cplx=.TRUE.   ,grid =.true.)
CALL init_psi(Hpsi   ,Basis,    cplx=.TRUE.   ,grid =.false.)

!=====================test================================
  CALL init_psi(B,   Basis,    cplx=.TRUE.   ,grid =.false.)
  CALL init_psi(G,   Basis,    cplx=.TRUE.   ,grid =.true.)


  !G1%CVec(:)= CZERO
  !call write_psi(B1)
  !print*,"========================"
  !print*,"shape(B1)=",shape(B1%CVec)
  !call bg2D(G1%CVec,B1%CVec,Basis)

  !call write_psi(G1)
  !print*,"========================"
  ! B1%CVec(:)= CZERO
   !call gb2D(B1%CVec,G1%CVec,Basis)
  !call GridTOBasis_nD_cplx(B1%CVec,G1%CVec,Basis)
  !call write_psi(B1)
 ! print*,"autre"
  !B1%CVec(:) = CZERO
  !G1%CVec(:)= CONE
  !call write_psi(G1)
  !call gb2D(B1%CVec,G1%CVec,Basis)
  print*,"========================"
  B%CVec(:) = CONE
  G%CVec(:)= CZERO
  CALL bg2D(G%CVec,B%CVec,Basis)
  !call BasisTOGrid_nD_cplx(G%CVec,B%CVec,Basis)
   call write_psi(G)
   B%CVec=CZERO
    print*,"========================"
    print*,"========================"
   !call GridTOBasis_nD_cplx(B%CVec,G%CVec,Basis)
    call gb2D(B%CVec,G%CVec,Basis)
   call write_psi(B)

 !call Init_GWP_parameters(Mparagwp,Basis,nb_GWP=2,nio=in_unitp)
 !call  init_psi0_nD(psi0%CVec,Basis,Basis%NDindexq%NDend,nio=in_unitp)
   ! call write_psi1(psi0,1)
 !call Calc_average_energy(psi0,Basis,E)
 !CALL test_basitogridgridtobasis(Basis)
 STOP 'calcul de Hpsi est fait'

  CALL read_propa(propa)
 CALL propagation(psif,psi0,propa,Basis)
  CALL Write_psi(psif)


  write(out_unitp,*) 'deallocation'
  CALL dealloc_psi(psi0)
  CALL dealloc_psi(psif)
  CALL dealloc_psi(Hpsi)

END PROGRAM TD_SCHROD
