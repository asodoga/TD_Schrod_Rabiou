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


   SUBROUTINE march(psi, psi_dt, t, propa)
      USE psi_m
      TYPE(propa_t), INTENT(IN)                :: propa
      TYPE(psi_t), INTENT(IN)                  :: psi
      TYPE(psi_t), INTENT(INOUT)               :: psi_dt
      real(kind=Rkind), INTENT(IN)             :: t

      real(kind=Rkind)                          :: Qt, sQt, Norm, Norm0
      character(len=Name_len)                   :: name


       name = propa%propa_name2 
      call string_uppercase_TO_lowercase(name)

      select case (name)
      case ('rk4') ! rk4 : Runge-kutta time propagation
         CALL marh_RK4th(psi, psi_dt, t, propa)
      case ('taylor') ! taylor : Taylor propagation
         CALL march_taylor(psi, psi_dt, t, propa)
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

   SUBROUTINE march_taylor(psi, psi_dt, t, propa)
      USE op_m
      USE psi_m
      USE Basis_m

      TYPE(psi_t), INTENT(INOUT)       :: psi_dt
      TYPE(psi_t), INTENT(IN)          :: psi
      TYPE(propa_t), INTENT(IN)        :: propa
      real(kind=Rkind), INTENT(IN)     :: t
      TYPE(psi_t)                      :: psi0
      TYPE(psi_t)                      :: Hpsi
      real(kind=Rkind)                 :: alpha
 
      ! variables locales-------------------------------------------------------------------------------

      real(kind=Rkind)                 :: Rkk, Norm, Norm0
      integer                          :: kk
      CALL init_psi(Hpsi, psi%basis, cplx=.TRUE., grid=.false.) ! to be changed
      CALL init_psi(Psi0, psi%basis, cplx=.TRUE., grid=.false.) ! to be changed

      write (out_unit, *) 'BEGINNIG march_taylor  ', t, propa%delta_t
      Rkk = ONE
      alpha = TEN**10
      Psi_dt%CVec = Psi%CVec
      Psi0%CVec = Psi%CVec
      Do kk = 1, propa%max_iter, 1
         CALL mEyeHPsi(psi0, Hpsi)

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


   SUBROUTINE read_propa(propa)
      USE psi_m
      implicit none
      TYPE(propa_t), intent(inout)   :: propa
      real(kind=Rkind)               :: t0, tf, delta_t, eps
      character(len=40)              :: propa_name, propa_name2
      integer                        :: max_iter
      logical                        :: Beta,P,renorm

      namelist /prop/ t0, tf, delta_t,max_iter,eps,propa_name, propa_name2 , Beta,P,renorm
      t0 = ZERO
      tf = TEN
      delta_t = ONETENTH**3
      eps = ONETENTH**10
      max_iter = 500
      propa_name = 'non_hagedorn'
      propa_name2 = 'rk4'
      Beta        = .true.
      P           = .true.
      renorm      = .true.

      read (*, nml=prop)

      propa%t0 = t0
      propa%tf = tf
      propa%delta_t = delta_t
      propa%eps = eps
       propa%max_iter = max_iter
      propa%propa_name = propa_name
      propa%propa_name2 = propa_name2
      propa%Beta =   Beta
      propa%P = P 
      propa%renorm =renorm
      call string_uppercase_TO_lowercase( propa%propa_name)
      call string_uppercase_TO_lowercase( propa%propa_name2)

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

   SUBROUTINE mEyeHPsi_temp(psi, Hpsi,H) !calcul de -iHpsi
      USE op_m
      USE psi_m
      TYPE(psi_t), intent(in)       :: psi
      TYPE(psi_t), intent(inout)    :: Hpsi
      TYPE(Op_t) ,intent(in)        :: H

      call calc_OpPsi(H, psi, Hpsi)
      Hpsi%CVec(:) = -EYE*Hpsi%CVec(:)

   END SUBROUTINE 

   SUBROUTINE write_propa(propa)
      USE psi_m
      implicit none
      TYPE(propa_t), intent(inout) :: propa

      write (out_unit, *) 't0 = ', propa%t0
      write (out_unit, *) 'tf = ', propa%tf
      write (out_unit, *) 'deltat_t = ', propa%delta_t
      write (out_unit, *) 'eps = ', propa%eps
      write (out_unit, *) 'max_iter = ', propa%max_iter
      write (out_unit, *) 'propa_name = ', propa%propa_name
      write (out_unit, *) 'propa_name2 = ', propa%propa_name2
       write (out_unit, *) 'Beta = ', propa%Beta
      write (out_unit, *) 'P = ', propa%P
      write (out_unit, *) 'renorm = ', propa%renorm

   END SUBROUTINE write_propa

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

         !Print*,"psi  is on Grid"
         CALL init_psi(psi_b, psi%Basis, cplx=.TRUE., grid=.false.)
         call GridTOBasis_nD_cplx(psi_b%CVec, psi%CVec, psi%Basis)
         CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.false.)
         !call calc_OpPsi(H, psi_b, Hpsi)
         stop 'verifie le calcul d E'
         E = real(dot_product(Hpsi%CVec, psi_b%CVec), kind=Rkind)

      else
         !Print*,"psi is on basis"
         CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.false.)
         !call calc_OpPsi(H, psi, Hpsi)
         stop 'verifie le calcul d E'
         E = real(dot_product(Hpsi%CVec, psi%CVec), kind=Rkind)

      end if
      call Calc_Norm_OF_Psi(psi, Norm)
      E = E/Norm**2
      !print *, "<psi|H|psi> = ", E, "<psi|psi> =", Norm

      CALL dealloc_psi(Hpsi)
      CALL dealloc_psi(psi_b)
   End SUBROUTINE Calc_average_energy


   SUBROUTINE Calc_Av_E(E,psi,H)
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

  SUBROUTINE  creat_file_unit(nio, name, propa)
      character(*), intent(in)    :: name
      type(propa_t), intent(in)   :: propa
      character(100)              :: name_tot
      integer, intent(in)         :: nio
      character(len=20)           :: dt
      character(len=8)            :: fmt

      fmt = "(E0.1)"

      write (dt, fmt) propa%delta_t
      name_tot = trim(name)//'_'//trim(propa%propa_name)//'_'//trim(propa%propa_name2)//'.txt'
     ! name_tot = 'file'//'_'//trim(name)!trim(name_tot)

      open (unit=nio, file=name_tot)

    END SUBROUTINE

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

      
end module sub_propa_m
