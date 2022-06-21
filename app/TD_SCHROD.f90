PROGRAM TD_SCHROD
  USE NumParameters_m
  USE Basis_m
  USE psi0_m
  USE psi_m
  USE op_m
  USE Propa_m
  IMPLICIT NONE
  TYPE (Basis_t), target      :: Basis
  TYPE (psi_t)                :: psi,Hpsi
  TYPE (psi_t)                :: psi0,psif
  TYPE(propa_t)               :: propa
  TYPE(psi_t)                 :: G
  TYPE(psi_t)                 :: B
  real(Kind = Rk)             :: E,dot_prdct
  TYPE(para_psi0_t)           :: paragwp0
  integer                     :: nio
 !====================================================================
  ! read some informations (basis set/grid) : numbers of basis functions, grid points ...
  ! the basis/grid informations have to be put in a module

  CALL Read_Basis(Basis,nio=in_unitp)

!====================================================================
!print*,"Basis is allocated",Basis_IS_allocated(Basis)

write(out_unitp,*) 'Initialization of a complex psi'
CALL init_psi(psi ,   Basis,    cplx=.TRUE.   ,grid =.true. ) ! to be changed
CALL init_psi(psi0,   Basis,    cplx=.TRUE.   ,grid =.true.) ! to be changed
CALL init_psi(psif,   Basis,    cplx=.TRUE.   ,grid =.true.) ! to be changed
CALL init_psi(G   ,   Basis,    cplx=.TRUE.   ,grid =.true.) ! to be changed
CALL init_psi(B   ,   Basis,    cplx=.TRUE.   ,grid =.false.) ! to be changed
CALL init_psi(Hpsi   ,Basis,    cplx=.TRUE.   ,grid =.false.) ! to be changed
   call Read_psi0(paragwp0,nio=in_unitp)
  call write_psi0(paragwp0, in_unitp)

   !call  GWP0(psi0,paragwp0,Basis)

! call GWP0(psi0,para_psi0,Basis)
  !CALL Calc_average_energy(psi0,Basis,E)
STOP 'calcul de Hpsi est fait'

  CALL read_propa(propa)
 CALL propagation(psif,psi0,propa,Basis)
  CALL Write_psi(psif)


  write(out_unitp,*) 'deallocation'
  !CALL dealloc_Op(H)
  CALL dealloc_psi(psi0)
  CALL dealloc_psi(psif)
  CALL dealloc_psi(psi)
  CALL dealloc_psi(Hpsi)

END PROGRAM TD_SCHROD
