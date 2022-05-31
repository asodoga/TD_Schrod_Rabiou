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
  TYPE(propa_t)          :: propa
  REAL(kind=Rk)          :: Norm ,E !,phase,sigma,k0,Q0,sig0
  TYPE(psi_t)                   :: G
  TYPE(psi_t)                   :: B ,diff
  integer                       :: iq
 !====================================================================
  ! read some informations (basis set/grid) : numbers of basis functions, grid points ...
  ! the basis/grid informations have to be put in a module

  CALL Read_Basis(Basis,nio=in_unitp)

!====================================================================
!print*,"Basis is allocated",Basis_IS_allocated(Basis)


    OPEN(unit=14,file = 'ENERGY' )
!  write(out_unitp,*) 'Initialization of a real psi'

  !CALL init_psi(psi,Basis,cplx=.FALSE.) ! to be changed
!  psi%RVec(:) = ONE

write(out_unitp,*) 'Initialization of a complex psi'
CALL init_psi(psi ,Basis,  cplx=.TRUE.   ,grid =.true. ) ! to be changed
CALL init_psi(psi0,Basis,  cplx=.TRUE.   ,grid =.true.) ! to be changed
CALL init_psi(psif,Basis,  cplx=.TRUE.   ,grid =.true.) ! to be changed
CALL init_psi(G   ,Basis,  cplx=.TRUE.   ,grid =.true.) ! to be changed
CALL init_psi(B   ,Basis,  cplx=.TRUE.   ,grid =.false.) ! to be changed
CALL init_psi(Hpsi   ,Basis,  cplx=.TRUE.,grid =.true.) ! to be changed
  CALL init_psi(diff   ,Basis,  cplx=.TRUE.,grid =.true.) ! to be changed
 CALL initial_wp(B,psi0,G,Basis)

 CALL Write_psi(G)




 write(out_unitp,*) ' | H | Psi > calculation'
  CALL Calc_Hpsi(G%CVec,Hpsi%CVec,Basis)
  CALL Write_psi(Hpsi)
  CALL Calc_average_energy(G,Basis,E)
  !CALL Set_op(H,Basis) ! to be change
  !CALL calc_OpPsi(H,psi,Hpsi)
  write(out_unitp,*) ' diff  calculation'
  do iq = 1, Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nb
      diff%CVec(iq) = Hpsi%CVec(iq)- E*G%CVec(iq)
      write(*,*) iq , diff%CVec(iq)
  end do

!STOP 'calcul de Hpsi est fait'

 ! WRITE(14,*) E
  psi0%CVec(:) = G%CVec(:)

 ! Norm = sqrt(real(dot_product(psi0%CVec,psi0%CVec), kind=Rk))

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
