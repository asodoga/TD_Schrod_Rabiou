module sub_propa_m
   USE QDUtil_m
   USE psi_m
   USE Op_m
   USE Ana_psi_m
   Use lanczos_m
   USE Sub_Vp_m

   implicit none

   TYPE propa_t
      real(kind=Rkind)               :: t0
      real(kind=Rkind)               :: tf
      real(kind=Rkind)               :: delta_t
      real(kind=Rkind)               :: eps
      integer                        :: max_iter
      character(len=:), allocatable  :: propa_name
      character(len=:), allocatable  :: propa_name2
      logical                        :: Beta 
      logical                        :: P
      logical                        :: renorm
   END TYPE propa_t

   public :: march_taylor, marh_RK4th,March,march_ITP,march_SIL,march_VP
   public :: mEyeHPsi, write_propa, Analyse, creat_file_unit, read_propa,diff2
   public :: test_taylor,march_taylor_nolim,Calc_average_energy

contains


SUBROUTINE march_temp(psi, psi_dt, t, propa,H)
   USE psi_m
   TYPE(propa_t), intent(in)                :: propa
   TYPE(psi_t), intent(in)                  :: psi
   TYPE(psi_t), intent(inout)               :: psi_dt
   TYPE(Op_t), intent(in)                   :: H
   real(kind=Rkind), intent(in)             :: t

   character(len=Name_len)                  :: name

    name = propa%propa_name2 
   call string_uppercase_TO_lowercase(name)
   select case (name)
   case ('rk4') ! rk4 : Runge-kutta time propagation
      CALL marh_RK4th(psi, psi_dt, t, propa)
   case ('taylor') ! taylor : Taylor propagation
      CALL march_taylor_temp(psi, psi_dt, t, propa,H)
   case ('sil') ! SIL: short iterative lanczos
     CALL march_SIL(psi, psi_dt, t, propa)
   case ('itp') ! ITP : imaginary times propagation
      call Imaginary_time_propagation(psi, psi_dt, propa)
   case ('vp') ! VP : Variational principle times propagation
      call march_VP(psi, psi_dt, t, propa)   
   case default
      write (out_unit, *) ' March name is not in the list'
   end select      
END SUBROUTINE


!=====================================================================
!> SUBROUTINE march
!! Performs time propagation of the wavefunction `psi`
!! according to the propagation method specified in `propa`.
!!
!! This routine acts as a dispatcher:
!! - Converts the propagation method name to lowercase.
!! - Calls the corresponding propagation routine:
!!     * RK4     : 4th-order Runge-Kutta
!!     * Taylor  : Taylor series expansion
!!     * SIL     : Short Iterative Lanczos
!!     * ITP     : Imaginary Time Propagation
!!     * VP      : Variational Principle
!!
!! Arguments:
!!   - psi     : (IN)     initial wavefunction
!!   - psi_dt  : (INOUT)  wavefunction after propagation
!!   - t       : (IN)     current time
!!   - propa   : (IN)     propagation parameters
!=====================================================================
SUBROUTINE march(psi, psi_dt, t, propa)
   USE psi_m

   TYPE(propa_t), INTENT(IN)     :: propa
   TYPE(psi_t),    INTENT(IN)    :: psi
   TYPE(psi_t),    INTENT(INOUT) :: psi_dt
   real(kind=Rkind), INTENT(IN)  :: t

   real(kind=Rkind)              :: Qt, sQt, Norm, Norm0
   character(len=Name_len)       :: name

   ! Retrieve and convert the propagation method name to lowercase
   name = propa%propa_name2
   call string_uppercase_TO_lowercase(name)

   ! Select the propagation method based on the method name
   select case (name)

   case ('rk4')        ! RK4: 4th-order Runge-Kutta time propagation
      CALL marh_RK4th(psi, psi_dt, t, propa)

   case ('taylor')     ! Taylor: Taylor series time propagation
      CALL march_taylor(psi, psi_dt, t, propa)

   case ('sil')        ! SIL: Short Iterative Lanczos method
      CALL march_SIL(psi, psi_dt, t, propa)

   case ('itp')        ! ITP: Imaginary Time Propagation
      CALL Imaginary_time_propagation(psi, psi_dt, propa)

   case ('vp')         ! VP: Variational Principle time propagation
      CALL march_VP(psi, psi_dt, t, propa)

   case default        ! Unknown propagation method
      write (out_unit, *) ' March name is not in the list'

   end select

END SUBROUTINE march


   SUBROUTINE CopyPsi(psi2,psi1)
     implicit none
    TYPE(psi_t), intent(in)                    :: psi1
    TYPE(psi_t), intent(inout)                 :: psi2
    real(Kind=Rkind), allocatable             ::Pt(:),Qt(:),SQt(:)
    complex(Kind=Rkind), allocatable          :: At(:)
    integer                                   :: Ndim,i1

    Ndim = size(psi1%Basis%tab_basis) - 1
    allocate(At(Ndim),SQt(Ndim),Pt(Ndim),Qt(Ndim))
     Qt(:) = ZERO; SQt(:) = ONE; Pt(:) = ZERO;At = CONE
       stop 'att copie imp sub-propa ligne 102'
      !call Calc_AVQ_SQ_nD(psi1, Qt, SQt)
      call Calc_Av_imp_k_nD(psi1,Pt)
      call  Calc_Avg_A_nD(psi1, At)
   
     !call construct_primitive_basis(psi2%Basis, Qt,Pt,At,SQt)
      psi2%CVec(:) = psi1%CVec(:)

 END SUBROUTINE


   SUBROUTINE Analyse(psi, t)
      implicit none
      TYPE(psi_t), INTENT(IN)              :: psi
      real(Kind=Rkind), allocatable        :: pop(:), Qm(:), Qp(:)
      real(kind=Rkind), intent(in)         :: t
      real(kind=Rkind)                     :: Norm, E
      integer                              :: Ndim
      Ndim = size(psi%Basis%tab_basis)
      allocate (Pop(psi%Basis%tab_basis(Ndim)%nb))
      allocate (Qp(psi%Basis%tab_basis(Ndim)%nb))
      allocate (Qm(Ndim - 1))
      pop(:) = ZERO
      Qm(:) = ZERO
      Qp(:) = ZERO
      E = ZERO; Norm = ZERO
      !-------------------------------------- beging Anapsi------------------------
      call Population(psi, pop)
      call Calc_average_energy(psi, E)
      call Calc_Norm_OF_psi(psi, Norm)
      call Qpop(psi, Qp)
      write (3, *) t, E, Norm, pop
      write (4, *) t, Qm
      write (5, *) t, Qp
      deallocate (pop, Qm, Qp)
      !-----------------------------------------And Anapsi---------------------------------

   END SUBROUTINE Analyse

!=====================================================================
!> SUBROUTINE march_taylor
!! Performs time propagation of the wavefunction `psi`
!! using the Taylor series expansion method.
!!
!! The algorithm:
!!   1. Initialize working wavefunctions (`Hpsi`, `psi0`).
!!   2. Copy the initial wavefunction into `psi_dt` and `psi0`.
!!   3. Iteratively apply the Hamiltonian using mEyeHPsi,
!!      scaling each term by Δt / k, and summing into `psi_dt`.
!!   4. Stop when the norm of the current term is smaller than
!!      the convergence threshold (`propa%eps`) or when `max_iter`
!!      is reached.
!!
!! Arguments:
!!   - psi     : (IN)     Initial wavefunction
!!   - psi_dt  : (INOUT)  Wavefunction after propagation
!!   - t       : (IN)     Current time
!!   - propa   : (IN)     Propagation parameters
!=====================================================================
   SUBROUTINE march_taylor(psi, psi_dt, t, propa)
      USE op_m
      USE psi_m
      USE Basis_m
   
      TYPE(psi_t),    INTENT(INOUT) :: psi_dt
      TYPE(psi_t),    INTENT(IN)    :: psi
      TYPE(propa_t),  INTENT(IN)    :: propa
      real(kind=Rkind), INTENT(IN)  :: t
   
      TYPE(psi_t)                   :: psi0      ! Temporary wavefunction
      TYPE(psi_t)                   :: Hpsi      ! Hamiltonian applied to psi
      real(kind=Rkind)              :: alpha     ! Stability threshold
   
      ! Local variables -------------------------------------------------
      real(kind=Rkind)              :: Rkk, Norm, Norm0
      integer                       :: kk
   
      ! Initialize work arrays for Hamiltonian application
      CALL init_psi(Hpsi, psi%basis, cplx=.TRUE., grid=.FALSE.)
      CALL init_psi(psi0, psi%basis, cplx=.TRUE., grid=.FALSE.)
   
      write (out_unit, *) 'BEGINNING march_taylor', t, propa%delta_t
   
      Rkk   = ONE
      alpha = TEN**10
   
      ! Initialize psi_dt and psi0 with the initial wavefunction
      psi_dt%CVec = psi%CVec
      psi0%CVec   = psi%CVec
   
      ! Taylor expansion loop
      DO kk = 1, propa%max_iter, 1
   
         ! Apply Hamiltonian: H|psi0> → Hpsi
         CALL mEyeHPsi(psi0, Hpsi)
   
         ! Scale by (Δt / kk) and accumulate into psi_dt
         psi0%CVec(:)   = Hpsi%CVec(:) * (propa%delta_t / kk)
         psi_dt%CVec(:) = psi_dt%CVec(:) + psi0%CVec(:)
   
         ! Reset Hpsi to zero for next iteration
         Hpsi%CVec(:) = CZERO
   
         ! Check the norm of the current term
         CALL Calc_Norm_OF_Psi(psi0, Norm)
         write (out_unit, *) 'sqrt(<Hpsi|Hpsi>) = ', kk, Norm
   
         ! Stop if the term is too large (instability)
         IF (Norm >= alpha) THEN
            STOP "wrong choice of delta_t"
   
         ! Convergence reached
         ELSEIF (Norm <= propa%eps) THEN
            write (out_unit, *) 'Taylor condition is fulfilled after', kk, 'iterations'
            EXIT
         END IF
   
      END DO
   
      ! Final norm check
      CALL Calc_Norm_OF_Psi(psi, Norm0)
      CALL Calc_Norm_OF_Psi(psi_dt, Norm)
      write (out_unit, *) '<psi_dt|psi_dt> = ', Norm, 'abs(<psi_dt|psi_dt> - <psi0|psi0>) =', abs(Norm0 - Norm)
   
      write (out_unit, *) 'END march_taylor'
   
      ! Free allocated memory
      CALL dealloc_psi(psi0)
      CALL dealloc_psi(Hpsi)
   
   END SUBROUTINE march_taylor


   SUBROUTINE march_taylor_temp(psi, psi_dt, t, propa,H)
      USE op_m
      USE psi_m
      USE Basis_m
      TYPE(psi_t), intent(inout)       :: psi_dt
      TYPE(psi_t), intent(in)          :: psi
      TYPE(propa_t), intent(in)        :: propa
      TYPE(Op_t), intent(in)           :: H
      real(kind=Rkind), intent(in)     :: t

      TYPE(psi_t)                      :: psi0
      TYPE(psi_t)                      :: Hpsi
      real(kind=Rkind)                 :: alpha

      ! variables locales-------------------------------------------------------------------------------

      real(kind=Rkind)                 :: Rkk, Norm, Norm0
      integer                          :: kk

      CALL init_psi(Hpsi, psi%basis, cplx=.true., grid=.false.) 
      CALL init_psi(Psi0, psi%basis, cplx=.true., grid=.false.) 

      write (out_unit, *) 'BEGINNIG march_taylor  ', t, propa%delta_t

      Rkk = ONE
      alpha = TEN**10
      Psi_dt%CVec = Psi%CVec
      Psi0%CVec = Psi%CVec

      Do kk = 1, propa%max_iter, 1

         call  mEyeHPsi_temp(psi0, Hpsi,H)
          psi0%CVec(:) = Hpsi%CVec(:)*(propa%delta_t/kk)
          psi_dt%CVec(:) = psi_dt%CVec(:) + psi0%CVec(:)
          Hpsi%CVec(:) = CZERO
         call Calc_Norm_OF_Psi(psi0, Norm)
         write (out_unit, *) 'sqrt(<Hpsi|Hpsi>) = ', kk, Norm

         if (Norm >= alpha) then

            stop "wrong choice of delta_t"

         elseif (Norm <= propa%eps) Then

            write (out_unit, *) 'Taylor condition is fulfild after', kk, 'iteration'
            exit

         End if

      End do

      CALL Calc_Norm_OF_Psi(Psi, Norm0)
      CALL Calc_Norm_OF_Psi(Psi_dt, Norm)
      write (out_unit, *) '<psi_dt|psi_dt> = ', Norm, 'abs(<psi_dt|psi_dt> - <psi0|psi0>)  =', abs(Norm0 - Norm)
      write (out_unit, *) 'END march_taylor'
      CALL dealloc_psi(psi0)
      CALL dealloc_psi(Hpsi)

   END SUBROUTINE 


 SUBROUTINE march_SIL(psi, psi_dt, t, propa)
     USE lanczos_m
     USE psi_m
     USE Basis_m
     TYPE(psi_t),  intent(inout) :: psi_dt
     TYPE(psi_t),  intent(in)    :: psi
     TYPE(propa_t),intent(in)    :: propa
     real(kind=Rkind),intent(in) :: t

       ! variables locales--------------------------------------------------------------------

     real(kind=Rkind)                   :: Norm, Norm0
     logical, parameter                 :: debug=.false.

         IF (debug) THEN
            !write(out_unit,*) 'psi_t',psi%CVec
           flush(out_unit)
         END IF
    
          write (out_unit, *) 'BEGINNIG march_SIL ', t, propa%delta_t
          call Calc_psi_step_cplx(psi_dt,psi,propa%eps,propa%delta_t)
    
     call Calc_Norm_OF_Psi(psi, Norm0)
     call Calc_Norm_OF_Psi(psi_dt, Norm)
     write (out_unit, *) '<psi_dt|psi_dt> = ', Norm, 'abs(<psi_dt|psi_dt> - <psi|psi>)  =', abs(Norm0 - Norm)
     write (out_unit, *) 'END march_SIL'

      IF (debug) THEN
        flush(out_unit)
     END IF

  END SUBROUTINE 


   SUBROUTINE marh_RK4th(psi, psi_dt, t, propa)
      USE op_m
      USE psi_m

      TYPE(psi_t), intent(inout)       :: psi_dt
      TYPE(psi_t), intent(in)          :: psi
      TYPE(psi_t)                      :: K1, K2, K3, K4, psi_inter
      TYPE(propa_t), intent(in)        :: propa

      real(kind=Rkind), intent(in)     :: t
      real(kind=Rkind)                 ::  Norm, Norm0
      integer                          :: iq

      !  variables locales

      call init_psi(K1, psi%basis, cplx=.true., grid=.false.)
      call init_psi(K2, psi%basis, cplx=.true., grid=.false.)
      call init_psi(K3, psi%basis, cplx=.true., grid=.false.)
      call init_psi(K4, psi%basis, cplx=.true., grid=.false.)
      call init_psi(psi_inter, psi%basis, cplx=.true., grid=.false.)

      write (out_unit, *) 'BEGINNIG march_RK4th', t, propa%delta_t
      psi_dt%CVec(:) = psi%CVec(:)
      CALL mEyeHPsi(psi, K1)

      psi_inter%CVec = psi%CVec + (propa%delta_t*HALF)*K1%CVec
      CALL mEyeHPsi(psi_inter, K2)
      psi_inter%CVec = psi%CVec + (propa%delta_t*HALF)*K2%CVec
      CALL mEyeHPsi(psi_inter, K3)
      psi_inter%CVec = psi%CVec + propa%delta_t*K3%CVec
      CALL mEyeHPsi(psi_inter, K4)
      psi_dt%CVec(:) = psi_dt%CVec(:) + (propa%delta_t*SIXTH)*(K1%CVec(:) + TWO*K2%CVec(:) + TWO*K3%CVec(:) + K4%CVec(:))
      CALL Calc_Norm_OF_Psi(Psi_dt, Norm)
      CALL Calc_Norm_OF_Psi(Psi, Norm0)

      write (out_unit, *) '<psi|psi> = ', Norm0, '<psi_dt|psi_dt> = ', Norm, 'abs(<psi|psi> - <psi_dt|psi_dt))=', ABS(Norm0 - Norm)
      write (out_unit, *) 'END marh_RK4th'
      call dealloc_psi(K1)
       call dealloc_psi(K2)
       call dealloc_psi(K3)
       call dealloc_psi(K4)
       call dealloc_psi(psi_inter)
   END SUBROUTINE marh_RK4th

   SUBROUTINE march_ITP(psi, psi_dt, propa,plus)
      USE op_m
      USE psi_m
      TYPE(psi_t), INTENT(INOUT)       :: psi_dt
      TYPE(psi_t), INTENT(IN)          :: psi
      logical,INTENT(IN)               :: plus

      !variables locales -------------------------------------------------------
      TYPE(psi_t)                      :: psi0
      TYPE(psi_t)                      :: Hpsi
      TYPE(propa_t), INTENT(IN)        :: propa
      real(kind=Rkind)                 ::  Norm
      integer                          :: kk,max_iter
      complex(kind=Rkind)              :: idelta_t,Rkk

      call  init_psi(Hpsi, psi%basis, cplx=.TRUE., grid=.false.) 
      call  init_psi(Psi0, psi%basis, cplx=.TRUE., grid=.false.) 

      write (out_unit, *) 'BEGINNIG march_ITP  '
      Rkk         = ONE
      max_iter = 25
      Psi_dt%CVec = Psi%CVec
      Psi0%CVec   = Psi%CVec

      if(plus .eqv. .true.) then
         idelta_t = EYE*propa%delta_t
      else
         idelta_t = -EYE*propa%delta_t
      end if   
      write (out_unit, *) 'Imaginary time step',idelta_t
      Do kk = 1,max_iter, 1

         CALL mEyeHPsi(psi0, Hpsi)
         Rkk            = Rkk*(idelta_t/kk)
         psi_dt%CVec(:) = psi_dt%CVec(:) + Rkk*Hpsi%CVec(:)
         psi0%CVec(:)   = Hpsi%CVec(:)
         Hpsi%CVec(:)   = CZERO

         call Calc_Norm_OF_Psi(psi0, Norm)
         Norm = abs(Rkk)*Norm
         write (out_unit, *) 'sqrt(<Hpsi|Hpsi>) = ', kk, Norm
      End do
      call Calc_Norm_OF_Psi(psi_dt, Norm)
      write (out_unit, *) 'sqrt(<psi_dt|psi_dt>) = ', Norm
      psi_dt%CVec(:) = psi_dt%CVec(:)/Norm 

      write (out_unit, *) 'END march_ITP'
      call dealloc_psi(psi0)
      call  dealloc_psi(Hpsi)
   END SUBROUTINE 


   SUBROUTINE Imaginary_time_propagation(psi, psi_dt, propa)
      USE op_m
      USE psi_m
      TYPE(psi_t), intent(inout)       :: psi_dt
      TYPE(psi_t), intent(in)          :: psi
      TYPE(propa_t), intent(in)        :: propa
     
      !  variables locales --------------------------------------------------------
      TYPE(psi_t)                      :: psi0
      real(kind=Rkind)                 :: E_old , E_new,Rkk,delta_E
      integer                          :: iq,kk,it,nt
      
      write (out_unit, *) 'BEGINNIG imaginary time propagation'

      call init_psi(psi0, psi%basis, cplx=.true., grid=.false.) 

      psi0%CVec(:) = psi%CVec(:)
      nt = 2000 

   Do it = 1,nt
    call march_ITP(psi=psi0, psi_dt=psi_dt, propa=propa,plus=.false.) 

    call  Calc_average_energy(psi0, E_old)
    call  Calc_average_energy(psi_dt, E_new)
    delta_E = abs(E_new-E_old)
    psi0%CVec(:) = psi_dt%CVec(:)
    psi_dt%CVec  = CZERO

        
    write (out_unit, *) '----------------------------------------------------------------'
    write (out_unit, *) '--E_old--',E_old,'--delta_E--',delta_E,'--it--',it
    write (out_unit, *) '--E_new--',E_new,'--delta_E--',delta_E,'--it--',it
    write (out_unit, *) '---------------------------------------------------------'

    if(delta_E <= ONETENTH**15) then
      write (out_unit, *) '---the relaxation is fulfild after---',it, '---iteration----'
      exit
   end if

   End Do     
     
  write(out_unit,*) ' End imaginary time propagation'
 
  call dealloc_psi(psi0) 
      
   END SUBROUTINE 

   !----------------------------------------------------------------------
! SUBROUTINE: read_propa
!
! DESCRIPTION:
!   Reads propagation parameters (time interval, time step, tolerance,
!   method names, and various options) from the "prop" namelist and
!   stores them in the "propa" structure.
!   Method names are converted to lowercase.
!----------------------------------------------------------------------

SUBROUTINE read_propa(propa)
   USE psi_m
   IMPLICIT NONE

   TYPE(propa_t), INTENT(INOUT) :: propa
   REAL(KIND=Rkind)             :: t0, tf, delta_t, eps
   CHARACTER(LEN=40)            :: propa_name, propa_name2
   INTEGER                      :: max_iter
   LOGICAL                      :: Beta, P, renorm

   NAMELIST /prop/ t0, tf, delta_t, max_iter, eps, propa_name, propa_name2, &
                   Beta, P, renorm

   ! Default values
   t0          = ZERO
   tf          = TEN
   delta_t     = ONETENTH**3
   eps         = ONETENTH**10
   max_iter    = 500
   propa_name  = 'non_hagedorn'
   propa_name2 = 'rk4'
   Beta        = .TRUE.
   P           = .TRUE.
   renorm      = .TRUE.

   ! Read the namelist
   READ (*, NML=prop)

   ! Assign to structure
   propa%t0          = t0
   propa%tf          = tf
   propa%delta_t     = delta_t
   propa%eps         = eps
   propa%max_iter    = max_iter
   propa%propa_name  = propa_name
   propa%propa_name2 = propa_name2
   propa%Beta        = Beta
   propa%P           = P 
   propa%renorm      = renorm

   ! Convert method names to lowercase
   CALL string_uppercase_TO_lowercase(propa%propa_name)
   CALL string_uppercase_TO_lowercase(propa%propa_name2)

END SUBROUTINE read_propa

   SUBROUTINE mEyeHPsi(psi, Hpsi) !calcul de -iHpsi
      USE op_m
      USE psi_m

      TYPE(psi_t), intent(in)       :: psi
      TYPE(psi_t), intent(inout)    :: Hpsi
      TYPE(Op_t)                    :: H
      stop 'verifie mEyeHPsi '
      !CALL calc_OpPsi(H, psi, Hpsi)

      Hpsi%CVec(:) = -EYE*Hpsi%CVec(:)
   END SUBROUTINE 

   !----------------------------------------------------------------------
! SUBROUTINE: mEyeHPsi_temp
!
! DESCRIPTION:
!   Computes the action of the Hamiltonian operator H on a wavefunction
!   psi, multiplies the result by -i, and stores it in Hpsi.
!
!   Mathematically:
!      Hpsi = -i * ( H * psi )
!
! ARGUMENTS:
!   psi   : Input wavefunction.
!   Hpsi  : Output wavefunction after applying -i*H.
!   H     : Hamiltonian operator.
!----------------------------------------------------------------------

SUBROUTINE mEyeHPsi_temp(psi, Hpsi, H)
   USE op_m
   USE psi_m
   TYPE(psi_t), INTENT(IN)    :: psi
   TYPE(psi_t), INTENT(INOUT) :: Hpsi
   TYPE(Op_t),  INTENT(IN)    :: H

   ! Apply the operator H to psi
   CALL calc_OpPsi(H, psi, Hpsi)

   ! Multiply by -i
   Hpsi%CVec(:) = -EYE * Hpsi%CVec(:)
END SUBROUTINE mEyeHPsi_temp

!----------------------------------------------------------------------
! SUBROUTINE: write_propa
!
! DESCRIPTION:
!   Writes the current propagation parameters stored in the "propa"
!   structure to the output unit (out_unit) in a clean, aligned format.
!   This is useful for logging or verifying the propagation settings
!   before or during the simulation.
!
! ARGUMENTS:
!   propa : Propagation parameters structure (read-only here).
!----------------------------------------------------------------------

SUBROUTINE write_propa(propa)
   USE psi_m
   IMPLICIT NONE

   TYPE(propa_t), INTENT(INOUT) :: propa

   WRITE (out_unit, '(A15, " : ", ES15.8)') 't0',         propa%t0
   WRITE (out_unit, '(A15, " : ", ES15.8)') 'tf',         propa%tf
   WRITE (out_unit, '(A15, " : ", ES15.8)') 'delta_t',    propa%delta_t
   WRITE (out_unit, '(A15, " : ", ES15.8)') 'eps',        propa%eps
   WRITE (out_unit, '(A15, " : ", I0)')     'max_iter',   propa%max_iter
   WRITE (out_unit, '(A15, " : ", A)')      'propa_name',  TRIM(propa%propa_name)
   WRITE (out_unit, '(A15, " : ", A)')      'propa_name2', TRIM(propa%propa_name2)
   WRITE (out_unit, '(A15, " : ", L1)')     'Beta',       propa%Beta
   WRITE (out_unit, '(A15, " : ", L1)')     'P',          propa%P
   WRITE (out_unit, '(A15, " : ", L1)')     'renorm',     propa%renorm
END SUBROUTINE write_propa

!----------------------------------------------------------------------
! SUBROUTINE: Calc_average_energy
!
! DESCRIPTION:
!   Calculates the average energy of a wavefunction Psi using the
!   Hamiltonian operator H. The energy is computed as <Psi | H | Psi>
!   and normalized by <Psi | Psi> to ensure proper scaling.
!
! ARGUMENTS:
!   Psi   : Input wavefunction.
!   E     : Output average energy.
!----------------------------------------------------------------------

   SUBROUTINE Calc_average_energy(Psi, E)
      !>-------------------------------------------------------
      !>     E = <Psi | H | Psi>
      !>--------------------------------------------------------
      USE QDUtil_m
      USE psi_m
      TYPE(psi_t), intent(in)                        :: psi
      REAL(kind=Rkind), intent(inout)                :: E
      TYPE(psi_t)                                    :: Hpsi, psi_b
      TYPE(Op_t)                                     :: H
      REAL(KIND=Rkind)                               :: Norm
      if (Psi%Grid) then

         Print*,"psi  is on Grid"
         CALL init_psi(psi_b, psi%Basis, cplx=.TRUE., grid=.false.)
         call GridTOBasis_nD_cplx(psi_b%CVec, psi%CVec, psi%Basis)
         CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.false.)
         call calc_OpPsi(H, psi_b, Hpsi)
         !stop 'verifie le calcul d E'
         write(out_unit,*) 'rappele toi , sur la calcule de l energie'
         E = real(dot_product(Hpsi%CVec, psi_b%CVec), kind=Rkind)

      else
         Print*,"psi is on basis"
         CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.false.)
         call calc_OpPsi(H, psi, Hpsi)
         !stop 'verifie le calcul d E'
         write(out_unit,*) 'rappele toi , sur la calcule de l energie'
         E = real(dot_product(Hpsi%CVec, psi%CVec), kind=Rkind)

      end if
      call Calc_Norm_OF_Psi(psi, Norm)
      E = E/Norm**2
      !print *, "<psi|H|psi> = ", E, "<psi|psi> =", Norm

      CALL dealloc_psi(Hpsi)
      CALL dealloc_psi(psi_b)
   End SUBROUTINE Calc_average_energy

   !----------------------------------------------------------------------
   ! SUBROUTINE: Calc_Av_E
   !
   ! DESCRIPTION:
   !   Computes the expectation value E = <Psi | H | Psi>
   !   Normalizes by <Psi|Psi> to ensure proper scaling.
   !
   ! ARGUMENTS:
   !   E     : Output average energy.
   !   psi   : Input wavefunction.
   !   H     : Hamiltonian operator.
   !----------------------------------------------------------------------

   SUBROUTINE Calc_Av_E(E, psi, H)
      !> -----------------------------------------------------------
      !> Computes the expectation value E = <Psi | H | Psi>
      !> Normalizes by <Psi|Psi> to ensure proper scaling.
      !> Includes a local debug flag to optionally print E.
      !> -----------------------------------------------------------
   
      USE QDUtil_m
      USE psi_m
      IMPLICIT NONE
   
      ! ===== Arguments =====
      REAL(KIND=Rkind), INTENT(INOUT) :: E
      TYPE(psi_t), INTENT(IN)         :: psi
      TYPE(Op_t),  INTENT(IN)         :: H
   
      ! ===== Local variables =====
      TYPE(psi_t)      :: Hpsi, psi_b
      REAL(KIND=Rkind) :: Norm
      LOGICAL          :: debug  ! Debug flag
   
      ! ===== Debug control =====
      debug = .FALSE.  ! Set to .TRUE. to print E
   
      ! ===== Expectation value calculation =====
      IF (psi%Grid) THEN
         ! Convert wavefunction from grid to basis representation
         CALL init_psi(psi_b, psi%Basis, cplx=.TRUE., grid=.FALSE.)
         CALL GridTOBasis_nD_cplx(psi_b%CVec, psi%CVec, psi%Basis)
   
         ! Apply operator H to wavefunction in basis
         CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.FALSE.)
         CALL calc_OpPsi(H, psi_b, Hpsi)
   
         ! Compute <psi_b | H | psi_b>
         E = REAL(DOT_PRODUCT(Hpsi%CVec, psi_b%CVec), KIND=Rkind)
      ELSE
         ! Directly apply H to wavefunction
         CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.FALSE.)
         CALL calc_OpPsi(H, psi, Hpsi)
   
         ! Compute <psi | H | psi>
         E = REAL(DOT_PRODUCT(Hpsi%CVec, psi%CVec), KIND=Rkind)
      END IF
   
      ! ===== Normalize by the norm of psi =====
      CALL Calc_Norm_OF_Psi(psi, Norm)
      E = E / Norm**2
   
      ! ===== Debug output =====
      IF (debug) THEN
         WRITE(*,'(A,ES16.8)') 'DEBUG: <psi|H|psi> normalized = ', E
      END IF
   
      ! ===== Cleanup =====
      CALL dealloc_psi(Hpsi)
      CALL dealloc_psi(psi_b)
   END SUBROUTINE Calc_Av_E

   SUBROUTINE Calc_Av_E_temp(E,psi,H)
      !>-------------------------------------------------------
      !>     E = <Psi | H | Psi>
      !>--------------------------------------------------------
      USE QDUtil_m
      USE psi_m
      REAL(kind=Rkind), intent(inout)                :: E
      TYPE(psi_t), intent(in)                        :: psi
      TYPE(Op_t) ,intent(in)                         :: H

      TYPE(psi_t)                                    :: Hpsi, psi_b
      REAL(kind=Rkind)                               :: Norm

      if (psi%Grid) then

         CALL init_psi(psi_b, psi%Basis, cplx=.true., grid=.false.)
         call GridTOBasis_nD_cplx(psi_b%CVec, psi%CVec, psi%Basis)
         CALL init_psi(Hpsi, psi%Basis, cplx=.true., grid=.false.)
         call calc_OpPsi(H, psi_b, Hpsi)
         E = real(dot_product(Hpsi%CVec, psi_b%CVec), kind=Rkind)

      else
         call init_psi(Hpsi, psi%Basis, cplx=.true., grid=.false.)
         call calc_OpPsi(H, psi, Hpsi)
         E = real(dot_product(Hpsi%CVec, psi%CVec), kind=Rkind)

      end if

      call Calc_Norm_OF_Psi(psi, Norm)
      E = E/Norm**2
      print *, "<psi|H|psi> = ", E, "<psi|psi> =", Norm
      call dealloc_psi(Hpsi)
      call dealloc_psi(psi_b)
   End SUBROUTINE 


   SUBROUTINE diff2()
    complex(kind=Rkind), allocatable         :: df(:, :)
    complex(kind=Rkind), allocatable         :: f1(:, :)
    complex(kind=Rkind), allocatable         :: f2(:, :)
    integer                             :: iostat, iq = 1, n = 6001, m = 784
    open (200, file='psi_dt_on_basis0_non_hagedorn_taylor.txt', status="old")
    open (201, file='psi_dt_on_basis0_hagedorn_taylor.txt', status="old")
    open (202, file='psi_dt_on_basis_diff.txt')
    allocate (f1(n, m), f2(n, m), df(n, m))
    do while (iq < n)
       read (200, *, IOSTAT=iostat) f1(iq, :)
       read (201, *, IOSTAT=iostat) f2(iq, :)
       iq = iq + 1
    end do
    df(:, :) = ZERO
    !df(:, :) = f1(:, :) - f2(:, :)
    do iq = 1, n
       df(iq,:) = f1(iq,:)-f2(iq,:)
       write (202, *)  iq,sqrt(real(dot_product(df(iq, :), df(iq, :)), kind= Rkind))
    end do
  END SUBROUTINE

  SUBROUTINE diff()
      real(kind=Rkind), allocatable         :: df(:, :)
      real(kind=Rkind), allocatable         :: f1(:, :)
      real(kind=Rkind), allocatable         :: f2(:, :)
      integer                             :: iostat, iq = 1, n = 1000, m = 3

      open (22, file='psi.dat', status="old")
      open (23, file='psih.dat', status="old")
      open (24, file='diff.dat')
      allocate (f1(n, m), f2(n, m), df(n, m))

      do while (iq < n)
         read (22, *, IOSTAT=iostat) f1(iq, :)
         read (23, *, IOSTAT=iostat) f2(iq, :)
         iq = iq + 1
      end do

      do iq = 1, n
         if (iq >= 550) then
            f1(iq, :) = ZERO
            f2(iq, :) = ZERO
         end if
         df(iq, 1) = f1(iq, 1)

         ! print*,iq,f(iq,:)
      end do
      df(:, :) = ZERO
      df(:, :) = f1(:, :)
      df(:, 3) = abs(f1(:, 3) - f2(:, 3))
      do iq = 1, 549
         !df(iq,3) = f1(iq,3)-f2(iq,3)
         !print*,iq,df(iq,:)   ,   abs(f1(iq,3)-f2(iq,3))
         write (24, *) df(iq, :)
      end do
    END SUBROUTINE

    !===============================================================
!> @brief   Crée et ouvre un fichier texte pour l'enregistrement
!>          des résultats d'une propagation.
!>
!> @param[in] nio     Numéro d'unité logique (unit) pour le fichier.
!> @param[in] name    Nom de base à utiliser pour créer le fichier.
!> @param[in] propa   Structure contenant les informations sur la
!>                    propagation (nom, type, pas de temps, etc.).
!>
!> @details
!> Le nom final du fichier est construit comme :
!>   name_propaName_propaName2.txt
!> Exemple :
!>   si name = "Qt", propa%propa_name = "STD", propa%propa_name2 = "nb10"
!>   alors le fichier sera : "Qt_STD_nb10.txt"
!>
!> Le fichier est ouvert avec le buffering désactivé (`BUFFERED='NO'`)
!> pour que chaque écriture soit immédiatement visible sur disque.
!===============================================================
SUBROUTINE creat_file_unit(nio, name, propa)

   !==== Arguments ====
   INTEGER,          INTENT(IN) :: nio       ! Numéro d'unité du fichier
   CHARACTER(*),     INTENT(IN) :: name      ! Nom de base du fichier
   TYPE(propa_t),    INTENT(IN) :: propa     ! Données de la propagation

   !==== Variables locales ====
   CHARACTER(100) :: name_tot   ! Nom complet du fichier
   CHARACTER(20)  :: dt         ! Pas de temps converti en chaîne
   CHARACTER(8)   :: fmt        ! Format pour afficher le pas de temps

   !==== Construction du format pour écrire delta_t ====
   fmt = "(E0.1)"
   WRITE(dt, fmt) propa%delta_t

   !==== Construction du nom complet du fichier ====
   name_tot = TRIM(name) // '_' // TRIM(propa%propa_name) // '_' // &
              TRIM(propa%propa_name2) // '.txt'

   !==== Ouverture du fichier ====
   ! BUFFERED='NO' → désactive le buffering pour écriture immédiate
   OPEN(UNIT=nio, FILE=name_tot, STATUS='UNKNOWN', ACTION='WRITE', &
        FORM='FORMATTED')

END SUBROUTINE creat_file_unit

     !==============================================================
     !> @brief
     !>   Write a complex autocorrelation value at a given time 
     !>   to a file unit with different output options.
     !>
     !> @details
     !>   Controlled by the character string 'op':
     !>     - "t_z"     : write (t, z)   where z is complex
     !>     - "t_norm"  : write (t, |z|)
     !>     - "t_re_im" : write (t, Re(z), Im(z))
     !>     - "t_re"    : write (t, Re(z))
     !>     - "t_im"    : write (t, Im(z))
     !>     - "z"       : write z only
     !>     - "norm"    : write |z| only
     !>     - "re_im"   : write Re(z), Im(z)
     !>
     !> @param[in] z    Complex autocorrelation value
     !> @param[in] nio  File unit number (must be open before calling)
     !> @param[in] op   Output format option (see list above)
     !> @param[in] t    Time value (required if op begins with "t_")
     !==============================================================
   SUBROUTINE Write_Complex(z, nio, op, t)
      IMPLICIT NONE
      COMPLEX(KIND=Rkind), INTENT(in) :: z(:)
      INTEGER, INTENT(in)              :: nio
      CHARACTER(*), INTENT(in)         :: op
      REAL(KIND=Rkind), INTENT(in)     :: t

      ! Test if z is a scalar (size=1) or vector (size>1)
      IF (SIZE(z) == 1) THEN
         CALL Write_Complex_scalar(z(1), nio, op, t)
      ELSE
         CALL Write_Complex_vector(z, nio, op, t)
      END IF

   END SUBROUTINE Write_Complex

   !=============================
   ! Cas scalaire
   !=============================
   SUBROUTINE Write_Complex_scalar(z, nio, op, t)
      IMPLICIT NONE
      COMPLEX(KIND=Rkind), INTENT(in) :: z
      INTEGER,             INTENT(in) :: nio
      CHARACTER(*),        INTENT(in) :: op
      REAL(KIND=Rkind),    INTENT(in) :: t

      REAL(KIND=Rkind) :: zr, zi, znorm
      zr    = REAL(z, KIND=Rkind)
      zi    = AIMAG(z)
      znorm = ABS(z)

      SELECT CASE (TRIM(ADJUSTL(op)))
      CASE ("t_z","t_re_im")
         WRITE(nio,*) t, zr, zi
      CASE ("t_norm")
         WRITE(nio,*) t, znorm
      CASE ("t_re")
         WRITE(nio,*) t, zr
      CASE ("t_im")
         WRITE(nio,*) t, zi
      CASE ("z","re_im")
         WRITE(nio,*) zr, zi
      CASE ("norm")
         WRITE(nio,*) znorm
      CASE DEFAULT
         WRITE(*,*) "Error: Unknown option in Write_Complex_scalar -> ", op
      END SELECT
   END SUBROUTINE Write_Complex_scalar

   !=============================
   ! Cas vecteur
   !=============================
   SUBROUTINE Write_Complex_vector(z, nio, op, t)
      IMPLICIT NONE
      COMPLEX(KIND=Rkind), DIMENSION(:), INTENT(in) :: z
      INTEGER,             INTENT(in) :: nio
      CHARACTER(*),        INTENT(in) :: op
      REAL(KIND=Rkind),    INTENT(in) :: t

      INTEGER :: i, n
      n = SIZE(z)

      SELECT CASE (TRIM(ADJUSTL(op)))
      CASE ("t_z","t_re_im")
         WRITE(nio,*) t, REAL(z(:),KIND=Rkind), AIMAG(z(:))
      CASE ("t_norm")
         WRITE(nio,*) t, ABS(z(:))
      CASE ("t_re")
         WRITE(nio,*) t, REAL(z(:),KIND=Rkind)
      CASE ("t_im")
         WRITE(nio,*) t, AIMAG(z(:))
      CASE ("z","re_im")
         WRITE(nio,*) REAL(z(:),KIND=Rkind), AIMAG(z(:))
      CASE ("norm")
         WRITE(nio,*) ABS(z(:))
      CASE DEFAULT
         WRITE(*,*) "Error: Unknown option in Write_Complex_vector -> ", op
      END SELECT
   END SUBROUTINE Write_Complex_vector


    SUBROUTINE march_VP(psi, psi_dt, t, propa)
      USE lanczos_m
      USE psi_m
      type(psi_t),   intent(INOUT)   :: psi_dt
      type(psi_t),   intent(IN)      :: psi
      type(propa_t), intent(IN)      :: propa
      real(kind=Rkind), INTENT(IN)   :: t

        ! variables locales -------------------------------------------------------------------------------

      real(kind=Rkind)                    :: Norm, Norm0,E 
      logical, parameter                  :: debug=.false.      

      IF (debug) THEN
         write(out_unit,*) 'psi_t',psi%CVec
        flush(out_unit)
      END IF
     
      write (out_unit, *) 'BEGINNIG march VP ', t, propa%delta_t
      call Vp_step_psi(psi, psi_dt,propa%delta_t)
      call Calc_Norm_OF_Psi(psi, Norm0)
      call Calc_Norm_OF_Psi(psi_dt, Norm)
      write (out_unit, *) '<psi_dt|psi_dt> = ', Norm, 'abs(<psi_dt|psi_dt> - <psi|psi>)  =', abs(Norm0 - Norm)
      write (out_unit, *) 'END march VP'
      
       IF (debug) THEN
         flush(out_unit)
      END IF

   END SUBROUTINE 



   SUBROUTINE GWP0(G0, Basis)
      USE  QDUtil_m
      USE Basis_m
      TYPE(Basis_t), intent(in), target                :: Basis
      complex(kind=Rkind) ,intent(inout)               :: G0(:)
   
      integer, allocatable                             :: Tab_iq(:)
      integer                                          :: inb, ndim, iq,nq,i
      real(Kind=Rkind)               , allocatable     :: Q(:)
      real(Kind=Rkind)                                 :: Q1,Q2,S1,S2
      logical                                          :: Endloop
   
      ndim = size(Basis%tab_basis) - 1
      nq =Basis%nq
      Q1=TWO; Q2=ZERO
      S1=1.29; S2=sqrt(TWO)
   
       allocate (Q(ndim),Tab_iq(ndim))
   
      Call Init_tab_ind(Tab_iq, Basis%NDindexq)
      Iq = 0
      DO
         Iq = Iq + 1
         CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop)
         IF (Endloop) exit
         do inb = 1, ndim
             Q(inb) = Basis%tab_basis(inb)%x(Tab_iq(inb))
         end do
        G0(Iq) = exp(-((Q(1)-Q1)/S1)**2+((Q(2)-Q2)/S2)**2) /(sqrt(sqrt(pi/TWO)*S1)*sqrt(sqrt(pi/TWO)*S2))
      END DO
       G0 = G0/(sqrt(sum(abs(G0(:))**2)))
       print*,"NG=",sqrt(sum(abs(G0(:))**2))
      deallocate(Tab_iq,Q)
   END SUBROUTINE 



   
   SUBROUTINE march_taylor_nolim(psi, psi_dt, propa,H)
      USE op_m
      USE psi_m
      USE Basis_m
      TYPE(psi_t), intent(inout)       :: psi_dt
      TYPE(psi_t), intent(in)          :: psi
      TYPE(propa_t), intent(in)        :: propa
      TYPE(Op_t), intent(in)           :: H

      TYPE(psi_t)                      :: psi0
      TYPE(psi_t)                      :: Hpsi
      real(kind=Rkind)                 :: alpha

      ! variables locales-------------------------------------------------------------------------------

      real(kind=Rkind)                 :: Rkk, Norm, Norm0
      integer                          :: kk

      CALL init_psi(Hpsi, psi%basis, cplx=.true., grid=.false.) 
      CALL init_psi(Psi0, psi%basis, cplx=.true., grid=.false.) 


      Rkk = ONE
      Psi_dt%CVec = Psi%CVec
      Psi0%CVec = Psi%CVec

      Do kk = 1, propa%max_iter, 1

         call  mEyeHPsi_temp(psi0, Hpsi,H)
          psi0%CVec(:) = Hpsi%CVec(:)*(propa%delta_t/kk)
          psi_dt%CVec(:) = psi_dt%CVec(:) + psi0%CVec(:)
          Hpsi%CVec(:) = CZERO
         call Calc_Norm_OF_Psi(psi0, Norm)
         write (5000, *)  kk, Norm

      End do

      CALL Calc_Norm_OF_Psi(Psi, Norm0)
      CALL Calc_Norm_OF_Psi(Psi_dt, Norm)
      write (out_unit, *) '<psi_dt|psi_dt> = ', Norm, 'abs(<psi_dt|psi_dt> - <psi0|psi0>)  =', abs(Norm0 - Norm)
      write (out_unit, *) 'END march_taylor'
      CALL dealloc_psi(psi0)
      CALL dealloc_psi(Hpsi)

   END SUBROUTINE 

   SUBROUTINE test_taylor(psi, psi_dt, propa, H)
      USE op_m
      USE psi_m
      USE Basis_m
      
      TYPE(psi_t), intent(inout)       :: psi_dt
      TYPE(psi_t), intent(in)          :: psi
      TYPE(propa_t), intent(in)        :: propa
      TYPE(Op_t), intent(in)           :: H
      
      CHARACTER(len=100)               :: delta_t_string
      CHARACTER(len=100)               :: name_tot
  
      ! Convertir delta_t en chaîne de caractères
      WRITE(delta_t_string, '(F6.2)') propa%delta_t 
  
      ! Construire le nom du fichier
     name_tot = 'ty_norm_' // TRIM(ADJUSTL(delta_t_string)) // '.txt'
  
      ! Ouvrir le fichier
      OPEN(unit=5000, file=name_tot)
  
      ! Appeler la fonction de propagation
      CALL march_taylor_nolim(psi, psi_dt, propa, H)
  END SUBROUTINE test_taylor


      
end module sub_propa_m
