PROGRAM TD_SCHROD
  USE NumParameters_m
  USE Basis_m
  USE GWP1D_m
  USE GWPnD_m
  USE GWPmnD_m
  USE psi_m
  USE psi0_m
  USE op_m
  USE Propa_m
  IMPLICIT NONE
  TYPE (Basis_t), target      :: Basis
  TYPE (psi_t)                :: psi,Hpsi
  TYPE (psi_t)                :: psi0,psif
  TYPE(propa_t)               :: propa
  TYPE(psi_t)                 :: G
  TYPE(psi_t)                 :: B
  real(Kind = Rk)             :: E

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
call init_psi0(psi0,Basis,I_ElecChannel=3,ndim=1, nb_GWP= 3)
  CALL Calc_average_energy(psi0,Basis,E)
   !CALL Write_psi1(psi0,no=1)
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
