PROGRAM TD_SCHROD
  USE NumParameters_m
  USE GWPnD_m
  USE Basis_m
  USE psi_m
  Use Molec_m
  USE Propa_m
  IMPLICIT NONE
  TYPE (Basis_t), target      :: Basis
  TYPE (psi_t)                ::Hpsi
  TYPE (psi_t)                :: psi0,psif,psi
  TYPE(propa_t)               :: propa
  real(Kind = Rk)             :: E,Norm
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

 !call Init_GWP_parameters(Mparagwp,Basis,nb_GWP=1,nio=in_unitp)
 call  init_psi0_nD(psi0%CVec,Basis,nio=in_unitp)
  call  Calc_Norm(psi0, Norm,Basis)
  psi0%CVec(:) = psi0%CVec(:)/Norm
 call Calc_average_energy(psi0,Basis,E)
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
