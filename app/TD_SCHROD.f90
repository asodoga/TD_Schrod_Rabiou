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
  REAL(kind=Rk)          :: Norm,phase,sigma,k0,Q0,sig0,Norm1
  !COMPLEX(kind=Rk),    ALLOCATABLE   :: G(:)
  !COMPLEX(kind=Rk),    ALLOCATABLE   :: B(:)


  TYPE(psi_t)                   :: G
  TYPE(psi_t)                   :: B
  INTEGER                       IQ , IB

  !ALLOCATE(G(Basis%NQ), B(Basis%NB))

  !====================================================================
  ! read some informations (basis set/grid) : numbers of basis functions, grid points ...
  ! the basis/grid informations have to be put in a module

  CALL Read_Basis(Basis,nio=in_unitp)

  !====================================================================

!  write(out_unitp,*) 'Initialization of a real psi'

  !CALL init_psi(psi,Basis,cplx=.FALSE.) ! to be changed
!  psi%RVec(:) = ONE
!  CALL Write_psi(psi)

OPEN(unit=11,file = 'norm' )
OPEN(unit=12,file = 'G' )
OPEN(unit=13,file = 'B' )


  CALL init_psi(G,Basis,cplx=.TRUE.)
  CALL init_psi(B,Basis,cplx=.TRUE.)
  sigma = 0.5
  sig0 = TWO
  k0 = 1
  phase = 0
  Q0 = 0
  G%CVec(:)  = EXP(-((Basis%x(:)-Q0)/(2d0*sig0))**2)* EXP(EYE*k0*Basis%x(:))
  !G%CVec(:)  = EXP(-(ONETENTH**3)*((Basis%x(:)-Q0)/sigma)**2)*EXP(EYE*k0*(Basis%x(:)-Q0)+ EYE*phase)

   CALL Calc_Norm_Grid(G, Norm1)
   G%CVec(:) = G%CVec(:)/Norm1
   CALL Calc_Norm_Grid(G, Norm1)
   DO  IQ = 1, G%Basis%nq
     write(12,*) G%Basis%x(IQ), ABS(G%CVec(IQ))**2
   ENDDO
 CALL GridTOBasis_Basis_cplx(B%CVec,G%CVec,Basis)
 CALL Calc_Norm(B, Norm)
 !B%CVec(:) = B%CVec(:)/Norm
 write(11,*) Norm1,Norm
 DO  IB = 1, G%Basis%nb
     write(13,*) G%Basis%x(IB), ABS(G%CVec(IB))**2
   ENDDO
  write(out_unitp,*) 'Initialization of a complex psi'
  CALL init_psi(psi,Basis,cplx=.TRUE.) ! to be changed
  psi%CVec(:) = B%CVec(:)
  CALL Write_psi(psi)

  write(out_unitp,*) ' | H | Psi > calculation'
  CALL Set_op(H,Basis) ! to be change
  CALL calc_OpPsi(H,psi,Hpsi)
  CALL Write_psi(Hpsi)


  CALL init_psi(psi0,Basis,cplx=.TRUE.) ! to be changed
  CALL init_psi(psif,Basis,cplx=.TRUE.) ! to be changed
  psi0%CVec(:) = B%CVec(:)
!  psi0%CVec(1) = ONE
  CALL Calc_Norm(psi0, Norm)
  !Norm = sqrt(real(dot_product(psi0%CVec,psi0%CVec), kind=Rk))
  write(out_unitp,*) 'norm,psi0',Norm
  CALL read_propa(propa)
  CALL propagation(psif,psi0,H,propa)
  CALL Write_psi(psif)


  write(out_unitp,*) 'deallocation'
  CALL dealloc_Op(H)
  CALL dealloc_psi(psi0)
  CALL dealloc_psi(psif)
  CALL dealloc_psi(psi)
  CALL dealloc_psi(Hpsi)

END PROGRAM TD_SCHROD
