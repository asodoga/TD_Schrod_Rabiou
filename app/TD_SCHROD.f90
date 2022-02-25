PROGRAM TD_SCHROD
  USE NumParameters_m
  USE Basis_m
  USE psi_m
  USE op_m
  USE Propa_m
  IMPLICIT NONE

  TYPE (Basis_t), target :: Basis
  TYPE (psi_t)           :: psi,Hpsi
  TYPE (psi_t)           :: psi0,psif
  TYPE (op_t)            :: H
  TYPE(propa_t)            :: propa
  real(kind=Rk)          :: Norm
  !real(kind=Rk)          :: t0,tf,delta_t,Norm,eps

  !====================================================================
  ! read some informations (basis set/grid) : numbers of basis functions, grid points ...
  ! the basis/grid informations have to be put in a module
  CALL Read_Basis(Basis,nio=in_unitp)
  
  !====================================================================

  write(out_unitp,*) 'Initialization of a real psi'

  CALL init_psi(psi,Basis,cplx=.FALSE.) ! to be changed
  psi%RVec(:) = ONE
  CALL Write_psi(psi)

  write(out_unitp,*) 'Initialization of a complex psi'
  CALL init_psi(psi,Basis,cplx=.TRUE.) ! to be changed
  psi%CVec(:) = CONE
  CALL Write_psi(psi)

  write(out_unitp,*) ' | H | Psi > calculation'
  CALL Set_op(H,Basis) ! to be change
  CALL calc_OpPsi(H,psi,Hpsi)
  CALL Write_psi(Hpsi)


  CALL init_psi(psi0,Basis,cplx=.TRUE.) ! to be changed
  CALL init_psi(psif,Basis,cplx=.TRUE.) ! to be changed
  psi0%CVec(:) = ZERO
  psi0%CVec(1) = ONE
  CALL Calc_Norm(psi0, Norm)
  !Norm = sqrt(real(dot_product(psi0%CVec,psi0%CVec), kind=Rk))
  write(out_unitp,*) 'norm,psi0',Norm
  CALL write_propa(propa)
  CALL propagation(psif,psi0,H,propa)
  CALL Write_psi(psif)


  write(out_unitp,*) 'deallocation'
  CALL dealloc_Op(H)
  CALL dealloc_psi(psi0)
  CALL dealloc_psi(psif)
  CALL dealloc_psi(psi)
  CALL dealloc_psi(Hpsi)

END PROGRAM TD_SCHROD
