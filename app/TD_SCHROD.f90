PROGRAM TD_SCHROD
  USE NumParameters_m
  USE GWPnD_m
  USE Basis_m
  USE psi_m
  !USE psi0_m
  USE Propa_m
  IMPLICIT NONE
  TYPE (Basis_t), target      :: Basis
  TYPE (psi_t)                ::Hpsi
  TYPE (psi_t)                :: psi0,psif
  TYPE(propa_t)               :: propa
  real(Kind = Rk)             :: E
 ! type(GWPnD_t)               paragwp

  !====================================================================
  ! read some informations (basis set/grid) : numbers of basis functions, grid points ...
  ! the basis/grid informations have to be put in a module
  CALL Read_Basis(Basis,nio=in_unitp)
!====================================================================
!print*,"Basis is allocated",Basis_IS_allocated(Basis)
write(out_unitp,*) 'Initialization of a complex psi'
CALL init_psi(psi0,   Basis,    cplx=.TRUE.   ,grid =.true.)
CALL init_psi(psif,   Basis,    cplx=.TRUE.   ,grid =.true.)
CALL init_psi(Hpsi   ,Basis,    cplx=.TRUE.   ,grid =.false.)
 !call Init_GWP_parameters(Mparagwp,Basis,nb_GWP=2,nio=in_unitp)
 call  init_psi0_nD(psi0%CVec,Basis,NDend=[100],nio=in_unitp)
    call write_psi1(psi0,1)
 call Calc_average_energy(psi0,Basis,E)
STOP 'calcul de Hpsi est fait'

  CALL read_propa(propa)
 CALL propagation(psif,psi0,propa,Basis)
  CALL Write_psi(psif)


  write(out_unitp,*) 'deallocation'
  CALL dealloc_psi(psi0)
  CALL dealloc_psi(psif)
  CALL dealloc_psi(Hpsi)

END PROGRAM TD_SCHROD
