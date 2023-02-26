PROGRAM TD_SCHROD
  USE GWPnD_m
  USE Basis_m
  USE psi_m
  USE op_m
  USE Propa_m
  USE Ana_psi_m
  USE lanczos_m
  IMPLICIT NONE
  TYPE (Basis_t), target       :: Basis_i,Basis_f
    TYPE (Op_t)                  :: H
  TYPE (psi_t)                 :: psi0,psif,psi
  TYPE(propa_t)                :: propa
  real(Kind = Rk)              :: E,Norm,AVQ,DQ 

!====================================================================
! for QML
integer :: ndim,nsurf,option,iq
logical :: adiabatic
character (len=16)                  :: pot_name

ndim      = 1
nsurf     = 1
pot_name  = 'read_model'
adiabatic = .false.
option    = 1
CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)
 write(out_unitp,*)'ndim,nsurf',ndim,nsurf
write(out_unitp,*) 'pot_name'
  !====================================================================
  ! read some informations (basis set/grid) : numbers of basis functions, grid points ...
  ! the basis/grid informations have to be put in a module
  CALL Read_Basis(Basis_i,nio=in_unitp)
  CALL construct_primitive_basis(Basis_i)
  CALL init_Basis1_TO_Basis2 (Basis_i,Basis_f)
  CALL construct_primitive_basis(Basis_f)
  !Call Write_Basis(Basis_f)
!====================================================================
!print*,"Basis is allocated",Basis_IS_allocated(Basis)
  write(out_unitp,*) 'Initialization of a complex psi'
  CALL init_psi(psi0,   Basis_i,    cplx=.TRUE.   ,grid =.false.)
  CALL init_psi(psif,   Basis_f,    cplx=.TRUE.   ,grid =.false.)
  CALL init_psi(psi,   Basis_i,    cplx=.TRUE.   ,grid =.false.)
  !CALL construct_primitive_basis(Basis_f)
 ! Call Write_Basis(Basis_f)
   CALL  GWP_init(psi0,1,nio=in_unitp)

   !psi0%CVec(:) = CZERO
   !psi0%CVec(1) =ONE
 ! call  Projection(psi0,psi,q10=ZERO,q20=ZERO,sci=ONE,scj=ONE )
  !call   init_basis(Basis1, psi%Basis,ONE,ONE)
  !CALL  init_psiHG(psi,Basis1)
  !CALL Write_psi(psi0)
  !CALL Hagedorn_psi1_TO_psi2(psif,psi0,0.0000000001_Rk,ONE)
 !CALL Test_Hagedorn(psi0,psif,q10=ZERO,q20=HALF,sci=ONE,scj=ONE)
  !call Calc_average_energy(psi0,E)
  !call Calc_std_dev_AVQ_1D(psi0,1,AVQ,DQ)
  !call Set_Op(H,Basis)
  ! CALL Make_Mat_OP(H)
  !call  write_Op(H)
 !STOP 'calcul de H|psi> est fait'

  CALL read_propa(propa)
  !CALL propagation_DP(Psif,Psi0,propa)
  CALL propagation_Hagedorn(Psif,Psi0,Basis_f,propa)
 ! CALL Write_psi(psif)


  write(out_unitp,*) 'deallocation'
  CALL dealloc_psi(psi0)






  CALL dealloc_psi(psif)


END PROGRAM TD_SCHROD
