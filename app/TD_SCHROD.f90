PROGRAM TD_SCHROD
   ! ======================================================
   ! PROGRAM: TD_SCHROD
   ! Quantum Dynamics Propagation using Wavepacket Methods
   ! Solves time-dependent Schrödinger equation
   ! ======================================================

   ! --------------------------
   ! Module declarations
   ! --------------------------
   USE QDUtil_m           ! Utility functions and constants
   USE Basis_m            ! Basis set definitions
   USE psi_m              ! Wavefunction operations
   USE op_m               ! Quantum operators (Hamiltonian)
   USE Param_WP0_m        ! Initial wavepacket parameters
   USE Propa_m            ! Propagation methods
   USE Ana_psi_m          ! Wavefunction analysis
   USE lanczos_m          ! Lanczos algorithm
   USE Sub_Vp_m           ! Variational principle subroutines

   IMPLICIT NONE

   ! --------------------------
   ! Main quantum system variables
   ! --------------------------
   TYPE(Basis_t), target  :: Basis    ! Basis set definition
   TYPE(Op_t)            :: H         ! Hamiltonian operator
   TYPE(psi_t)           :: psi0      ! Initial wavefunction
   TYPE(psi_t)           :: psif      ! Final wavefunction
   TYPE(propa_t)         :: propa     ! Propagation parameters
   
   ! --------------------------
   ! Wavepacket specific variables
   ! --------------------------
   TYPE(GWP_t), allocatable :: tab_GWP(:)  
   integer, allocatable    :: Tab_iq(:,:)  

   ! --------------------------
   ! QML (Quantum Model Library) interface
   ! --------------------------
   integer :: ndim = -1               ! Default Number of dimensions
   integer :: nsurf = -1              ! Default Number of surfaces
   logical :: adiabatic = .false.     ! Adiabatic representation flag
   integer :: option = 1              ! QML option
   character(len=16) :: pot_name = 'read_model'  ! Potential name

   ! ======================================================
   ! Initialization Section
   ! ======================================================

   ! Initialize QML quantum model
   call sub_Init_Qmodel(ndim, nsurf, pot_name, adiabatic, option)
   write (out_unit, *) 'ndim,nsurf', ndim, nsurf
   write (out_unit, *) 'pot_name: ', pot_name
   
   ! --------------------------
   ! Basis set initialization
   ! --------------------------
   call Read_Basis(Basis, nio=in_unit)       ! Read basis parameters
   call construct_primitive_basis(Basis)      ! Construct basis set

   ! --------------------------
   ! Wavefunction initialization
   ! --------------------------
   write (out_unit, *) "----- Initialization of psi0 -----"
   call init_psi(psi0, Basis, cplx=.true., grid=.false.)
   call init_psi(psif, Basis, cplx=.true., grid=.false.)
   
   ! Initialize Gaussian Wavepacket(s)
   call Read_tab_GWP(tab_GWP=tab_GWP, nb_GWP=1, nio=in_unit)
   
   ! --------------------------
   ! Operator initialization
   ! --------------------------
   call Calc_tab_Iq0(Tab_Iq, psi0%Basis)  
   call Set_Op(H, psi0%Basis, Tab_Iq)     ! Set up Hamiltonian

   ! Set initial wavefunction from GWP
   call psi_init_GWP(psi=psi0, Tab_GWP=tab_GWP)

   ! ======================================================
   ! Propagation Section
   ! ======================================================
   
   ! Read propagation parameters
   call read_propa(propa)

   ! Main propagation routine
   call propagation(psif, psi0, propa)
   print*, 'Propagation completed'

   ! ======================================================
   ! Cleanup Section
   ! ======================================================
   write (out_unit, *) 'Deallocating memory...'
   call dealloc_psi(psi0)
   call dealloc_psi(psif)

END PROGRAM TD_SCHROD