module Propa_m
   USE QDUtil_m
   USE psi_m
   USE Basis_m
   USE Op_m
   USE Ana_psi_m
   Use lanczos_m
   USE Auto_corr_m

   implicit none

   TYPE propa_t
      real(kind=Rkind)      :: t0
      real(kind=Rkind)      :: tf
      real(kind=Rkind)      :: delta_t
      real(kind=Rkind)      :: eps
      integer             :: max_iter
      integer             :: Kmax
      character(len=:), allocatable  :: propa_name
      character(len=:), allocatable  :: propa_name2
   END TYPE propa_t

   public :: march_taylor, marh_RK4th, read_propa
   public :: propagation, Hagedorn,H_test
   public :: mEyeHPsi, write_propa, Analyse, creat_file_unit

contains

   SUBROUTINE march(psi, psi_dt, t, propa)
      USE psi_m
      TYPE(propa_t), INTENT(IN)                :: propa
      TYPE(psi_t), INTENT(IN)                  :: psi
      TYPE(psi_t), INTENT(INOUT)               :: psi_dt
      real(kind=Rkind), INTENT(IN)             :: t

      real(kind=Rkind)                          :: Qt, sQt, Norm, Norm0

      select case (propa%propa_name2)
      case ('rk4')
         CALL marh_RK4th(psi, psi_dt, t, propa)
      case ('taylor')
         CALL march_taylor(psi, psi_dt, t, propa)
      case ('SIL')
        CALL march_SIL(psi, psi_dt, t, propa)
      case ('ITP') ! ITP : imaginary times propagation
         call Imaginary_time_propagation(psi, psi_dt, propa)
      case ('VP') ! VP : Variational principle times propagation
         call march_VP(psi, psi_dt, t, propa)   
      case default
         write (out_unit, *) 'name is not in the list'
      end select      

   END SUBROUTINE

   SUBROUTINE propagation(psif, psi0, propa)
      USE psi_m
      USE Basis_m

      TYPE(psi_t), intent(inout)       :: psif
      TYPE(psi_t), intent(in)          :: psi0
      TYPE(propa_t), intent(inout)     :: propa
      logical, parameter               :: debug = .true.

      ! variables locales------------------------------------------------------------------
      TYPE(Basis_t)                       :: Basis1,Basis2
      REAL(kind=Rkind)                    :: t, t_deltat, Norm, E, y
      REAL(kind=Rkind), allocatable       :: Qt(:), SQt(:), Auto_corr_function(:), populat(:),Pt(:)
      complex(kind=Rkind)                 ::  x
      complex(kind=Rkind) ,allocatable    :: Alpha(:)
      integer                             :: Ndim

      INTEGER                          :: i, nt, Iq, nf
      TYPE(psi_t)                      :: psi, psi_dt, psi00
      if (debug) then

         write (out_unit, *) 'BEGINNIG propagation', propa%t0, propa%tf, propa%delta_t
         ! write(out_unit,*) ''

         write (out_unit, *) '-------------propagation parameters---------------'
         Call write_propa(propa)
      else
         STOP ' check your data!'
         flush (out_unit)

      end if
      call creat_file_unit(nio=10, name='psi', propa=propa)
      call creat_file_unit(nio=11, name='Qt', propa=propa)
      call creat_file_unit(nio=12, name='E', propa=propa)
      call creat_file_unit(nio=13, name='SQt', propa=propa)
      call creat_file_unit(nio=14, name='Norm', propa=propa)
      call creat_file_unit(nio=15, name='Auto_corr_func', propa=propa)
      !call creat_file_unit(nio=16, name='psi_Ha', propa=propa)
      !call creat_file_unit(nio=17, name='psi_NHa', propa=propa)
      call creat_file_unit(nio=18, name='pop', propa=propa)
      call creat_file_unit(nio=19, name='Imp_k', propa=propa)
      call creat_file_unit(nio=20, name='alpha', propa=propa)

      Ndim = size(psi0%Basis%tab_basis) - 1
      allocate (Qt(Ndim), SQt(Ndim),Pt(Ndim),Alpha(Ndim))
      allocate (populat(psi0%Basis%tab_basis(Ndim + 1)%nb))
      Qt(:) = ZERO; SQt(:) = ONE
      nt = int((propa%tf - propa%t0)/propa%delta_t)
      E = ZERO;Alpha= CZERO


      call init_Basis1_TO_Basis2(Basis1, psi0%Basis)
      call init_Basis1_TO_Basis2(Basis2, psi0%Basis)
      call construct_primitive_basis(Basis1)
      call construct_primitive_basis(Basis2)


      call init_psi(psi, Basis1, cplx=.TRUE., grid=.false.)
      call init_psi(psi_dt,Basis2, cplx=.TRUE., grid=.false.)
      call init_psi(psi00, psi0%Basis, cplx=.TRUE., grid=.false.)

      psi%CVec(:) = psi0%CVec(:)
      !call write_psi(psi=psi, psi_cplx=.false., print_psi_grid=.true. &
      !               , print_basis=.false., t=ZERO, int_print=20, real_part=.false.)
      !call Calc_average_energy(Psi0, E)
      !call Calc_average_energy(psi, E)
      ! ---------------------------------- Beging  propagation----------------------------------------------------------
      DO i = 0, nt
         t = i*propa%delta_t
         t_deltat = t + propa%delta_t
         write (out_unit, *) propa%propa_name2, i, t, t_deltat
         call Calc_AVQ_nD0(Psi0=psi, AVQ=Qt, SQ=SQt)
         call Calc_average_energy(psi, E)
         call Calc_Norm_OF_Psi(psi, Norm)
         call Calc_Av_imp_k_nD(psi,Pt)
         !call Calc_Integral_cplx(psi, Alpha, 2)
         !call Population(psi, populat)

          write (11,*) t, Qt
          write (12,*) t, E
          write (13,*) t, SQt
          write (14,*) t, Norm
          write (18,*) t, populat
          write (19,*) t, Pt
         !write (20,*) t, Alpha
         !write (11, '(F18.6,2X,F18.6,F18.6,2X,F18.6)') t, Qt
         !write (12, '(F18.6,2X,F18.6,F18.6,2X,F18.6)') t, E
         !write (13, '(F18.6,2X,F18.6,F18.6,2X,F18.6)') t, SQt
         !write (14, '(F18.6,2X,F18.6,F18.6,2X,F18.6)') t, Norm
         !write (18, '(F18.6,2X,F18.6,F18.6,2X,F18.6)') t, populat(:)
         !write (19, '(F18.6,2X,F18.6,F18.6,2X,F18.6)') t, Pt(:)
         !write (20, '(F18.6,2X,F18.6,F18.6,2X,F18.6)') t, Alpha

         if (mod(i, 1) == 0) then
            call write_psi(psi=psi, psi_cplx=.true., print_psi_grid=.false. &
                           , print_basis=.false., t=t, int_print=10, real_part=.false.)
             write(10,*)
         end if

         CALL march(psi, psi_dt, t, propa)
         if (propa%propa_name == 'hagedorn') Then
            call Hagedorn(psi, psi_dt)
         else
            psi%CVec(:) = psi_dt%CVec(:)
         end if
         if(propa%propa_name2 == 'ITP') exit
       ! if (propa%propa_name == 'hagedorn') Then
       !    call Calc_Auto_corr(psi0, psi, x, y, propa%propa_name)
       !    write (15, '(F10.6,2X,F10.6,F10.6,2X,F10.6)') t, abs(x), y
       !    psi00%CVec = CZERO
       !    call Hagedorn0(psi00, psi)
       !    call write_psi(psi=psi00, psi_cplx=.true., print_psi_grid=.false. &
       !                   , print_basis=.false., t=t, int_print=16, real_part=.true.)
       !    write (16, *), ''
       ! else
       !    call Calc_Auto_corr(psi0, psi_dt, x, y, propa%propa_name)
       !    write (15, '(F18.6,2X,F18.6,F18.6,2X,F18.6)') t, abs(x), y
       !    call write_psi(psi=psi, psi_cplx=.true., print_psi_grid=.false. &
       !                   , print_basis=.false., t=t, int_print=17, real_part=.true.)
       !    write (17, *), ''
       ! end if

      END DO
      psif = psi_dt
      CALL Calc_Norm_OF_Psi(psif, Norm)
      IF (debug) THEN
         write (out_unit, *) 'END propagation'
         write (out_unit, *) 'norm,psi_dt', Norm

         ! if (propa%propa_name == 'hagedorn') Then
         !    psi00%CVec = CZERO
         !    call Hagedorn0(psi00, psi)
         !    call write_psi(psi=psi00, psi_cplx=.true., print_psi_grid=.false. &
         !                   , print_basis=.false., t=t, int_print=16, real_part=.true.)
         ! else
         !    call write_psi(psi=psif, psi_cplx=.true., print_psi_grid=.false. &
         !                   , print_basis=.false., t=t, int_print=17, real_part=.true.)
         ! end if
!
         flush (out_unit)
      END IF

   END SUBROUTINE

   SUBROUTINE Hagedorn(psi, psi_dt)
      USE psi_m
      USE Basis_m

      TYPE(psi_t), intent(inout)                      :: psi, psi_dt
      logical, parameter                              :: debug = .true.

      ! variables locales--------------------------------------------------------------
      REAL(kind=Rkind), allocatable                   :: Qt(:), SQt(:),Pt(:)
      REAL(kind=Rkind)                                :: Norm, E,E0
      integer                                         :: Ndim
     ! write (out_unit, *) 'Beging Hagedorn'
      !call Write_Basis(psi%Basis)

      call Calc_Norm_OF_Psi(psi_dt,Norm)
       call Calc_average_energy(psi_dt, E0)
       write(out_unit,*) 'HAgedorn in <psi|psi> =',Norm,'E0=',E0

      Ndim = size(psi_dt%Basis%tab_basis) - 1
      allocate (Qt(Ndim), SQt(Ndim),Pt(Ndim))
      Qt(:) = ZERO; SQt(:) = ONE; Pt(:) = ZERO

      call Calc_AVQ_nD0(psi0=psi_dt, AVQ=Qt, SQ=SQt)
      call Calc_Av_imp_k_nD(psi_dt,Pt)

      call construct_primitive_basis(psi_dt%Basis, Qt,Pt,SQt)
      call projection(psi, psi_dt)
      call construct_primitive_basis(psi%Basis, Qt,Pt,SQt)
     
      !write (out_unit, *) 'End Hagedorn'

      call Calc_Norm_OF_Psi(psi,Norm)
       call Calc_average_energy(psi, E)
      write(out_unit,*) 'HAgedorn out <psi_dt|psi_dt> =',Norm,'E=',E
      IF (debug) THEN
         flush (out_unit)
      END IF

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
      !===================================== beging Anapsi==================
      call Population(psi, pop)
      call Calc_average_energy(psi, E)
      call Calc_Norm_OF_psi(psi, Norm)
      call Qpop(psi, Qp)
      write (3, *) t, E, Norm, pop
      write (4, *) t, Qm
      write (5, *) t, Qp
      deallocate (pop, Qm, Qp)
      !====================================And Anapsi========================

   END SUBROUTINE Analyse

   SUBROUTINE march_taylor(psi, psi_dt, t, propa)
      USE op_m
      USE psi_m
      USE Basis_m

      TYPE(psi_t), INTENT(INOUT)       :: psi_dt
      TYPE(psi_t), INTENT(IN)          :: psi
      TYPE(psi_t)                      :: psi0
      TYPE(psi_t)                      :: Hpsi
      TYPE(propa_t), INTENT(IN)        :: propa
      real(kind=Rkind), INTENT(IN)     :: t
      real(kind=Rkind)                 :: alpha
 
      ! variables locales

      real(kind=Rkind)                 :: Rkk, Norm, Norm0
      integer                          :: kk
      CALL init_psi(Hpsi, psi%basis, cplx=.TRUE., grid=.false.) ! to be changed
      CALL init_psi(Psi0, psi%basis, cplx=.TRUE., grid=.false.) ! to be changed

      write (out_unit, *) 'BEGINNIG march_taylor  ', t, propa%delta_t
      !write(out_unit,*) 'psi',psi%CVec
      Rkk = ONE
      alpha = TEN**10
        !!======================debut ordre 1==========================
      !psi_dt%CVec    = psi%CVec
      !CALL calc_OpPsi(H,psi,Hpsi)
      !Rkindk = Rkindk*delta_t
      !Hpsi%CVec(:)    = - EYE*Hpsi%CVec(:)
      !psi_dt%CVec(:) = psi_dt%CVec(:) +Rkindk*Hpsi%CVec(:)
      !CALL Calc_Norm_OF_Psi(HPsi,Norm)
      !CALL Calc_Norm_OF_Psi(Psi_dt,Norm)

      !write(out_unit,*) 'norm,Hpsi',Rkindk*Norm
      !write(out_unit,*) 'norm,psi_dt',Norm

      !write(out_unit,*) 'psi_dt',psi_dt%CVec
      !write(out_unit,*) 'END march_taylor'

        !!===========================Ordre deux etplus=======================
      Psi_dt%CVec = Psi%CVec
      Psi0%CVec = Psi%CVec
      Do kk = 1, propa%max_iter, 1
         CALL mEyeHPsi(psi0, Hpsi)
         Rkk = Rkk*(propa%delta_t/kk)
         psi_dt%CVec(:) = psi_dt%CVec(:) + Rkk*Hpsi%CVec(:)
         psi0%CVec(:) = Hpsi%CVec(:)
         Hpsi%CVec(:) = CZERO
         call Calc_Norm_OF_Psi(psi0, Norm)
         Norm = Rkk*Norm
         write (out_unit, *) 'sqrt(<Hpsi|Hpsi>) = ', kk, Norm
         if (Norm >= alpha) then
            stop "wrong choice of delta_t"
         elseif (Norm <= propa%eps) Then

            print *, 'Taylor condition is fulfild after', kk, 'iteration'
            exit
         End if
      End do
      CALL Calc_Norm_OF_Psi(Psi, Norm0)
      CALL Calc_Norm_OF_Psi(Psi_dt, Norm)
      write (out_unit, *) '<psi_dt|psi_dt> = ', Norm, 'abs(<psi_dt|psi_dt> - <psi0|psi0>)  =', abs(Norm0 - Norm)
      write (out_unit, *) 'END march_taylor'
      !CALL dealloc_psi(psi0)
      !CALL dealloc_psi(Hpsi)
   END SUBROUTINE march_taylor


 SUBROUTINE march_SIL(psi, psi_dt, t, propa)
     USE lanczos_m
     USE psi_m
     USE Basis_m
     TYPE(psi_t), INTENT(INOUT)   :: psi_dt
     TYPE(psi_t), INTENT(IN)      :: psi
     TYPE(propa_t), INTENT(IN)    :: propa
     real(kind=Rkind), INTENT(IN)    :: t

       ! variables locales ============================================================================

     real(kind=Rkind)                    :: Norm, Norm0,E
      complex (kind=Rkind), allocatable  :: Vec_Basis(:,:)
      real (kind=Rkind), allocatable     :: EigenVal(:)
       TYPE(psi_t)        :: Psi0
     logical, parameter               :: debug=.false.

         IF (debug) THEN
            !write(out_unit,*) 'psi_t',psi%CVec
           flush(out_unit)
         END IF
    
          write (out_unit, *) 'BEGINNIG march_SIL ', t, propa%delta_t
          write(out_unit,*) 'Krylov Basis size',propa%Kmax
          call init_psi(psi0, psi%Basis, cplx=.TRUE., grid=.false.)
          CALL Lanczos_eign_syst_solve( EigenVal,Vec_Basis,psi,propa%kmax)
          psi0%CVec = Vec_Basis(:,1)
           call Calc_average_energy(psi0, E)
          print*,'==E1==',EigenVal(1:5),E,'=='
          CALL Calc_psi_step_cplx(psi_dt,psi,propa%delta_t,propa%Kmax)
    
     CALL Calc_Norm_OF_Psi(psi, Norm0)
     CALL Calc_Norm_OF_Psi(psi_dt, Norm)
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
      ! write(out_unit,*) 'psi'
      psi_dt%CVec(:) = psi%CVec(:)
      CALL mEyeHPsi(psi, K1)

      psi_inter%CVec = psi%CVec + (propa%delta_t*HALF)*K1%CVec
      !CALL Write_psi_basis(psi_inter,t,10)
      CALL mEyeHPsi(psi_inter, K2)
      psi_inter%CVec = psi%CVec + (propa%delta_t*HALF)*K2%CVec
      !CALL Write_psi_basis(psi_inter,t,11)
      print *, 'cc'
      CALL mEyeHPsi(psi_inter, K3)
      psi_inter%CVec = psi%CVec + propa%delta_t*K3%CVec
      !CALL Write_psi_basis(psi_inter,t,12)
      CALL mEyeHPsi(psi_inter, K4)
      psi_dt%CVec(:) = psi_dt%CVec(:) + (propa%delta_t*SIXTH)*(K1%CVec(:) + TWO*K2%CVec(:) + TWO*K3%CVec(:) + K4%CVec(:))
      CALL Calc_Norm_OF_Psi(Psi_dt, Norm)
      CALL Calc_Norm_OF_Psi(Psi, Norm0)

      write (out_unit, *) '<psi|psi> = ', Norm0, '<psi_dt|psi_dt> = ', Norm, 'abs(<psi|psi> - <psi_dt|psi_dt))=', ABS(Norm0 - Norm)
      write (out_unit, *) 'END marh_RK4th'
      deallocate (k1%CVec)
      deallocate (k2%CVec)
      deallocate (k3%CVec)
      deallocate (k4%CVec)
      deallocate (psi_inter%CVec)
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
     
      !  variables locales ==============================================================================
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

        
    write (out_unit, *) '================================================================='
    write (out_unit, *) '--E_old--',E_old,'--delta_E--',delta_E,'--it--',it
    write (out_unit, *) '--E_new--',E_new,'--delta_E--',delta_E,'--it--',it
    write (out_unit, *) '================================================================'

    if(delta_E <= ONETENTH**15) then
      write (out_unit, *) 'the relaxation is fulfild after',it, 'iteration'
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
      integer                        ::  max_iter,kmax

      namelist /prop/ t0, tf, delta_t, eps, max_iter,Kmax,propa_name, propa_name2
      t0 = ZERO
      tf = 10._Rkind
      delta_t = 0.001
      eps = ONETENTH**10
      max_iter = 5000
      Kmax     = 10
      propa_name = 'non_hagedorn'
      propa_name2 = 'rk4'

      read (*, nml=prop)

      propa%t0 = t0
      propa%tf = tf
      propa%delta_t = delta_t
      propa%eps = eps
      propa%max_iter = max_iter
      propa%Kmax = Kmax
      propa%propa_name = propa_name
      propa%propa_name2 = propa_name2

   END SUBROUTINE read_propa
   SUBROUTINE mEyeHPsi(psi, Hpsi) !calcul de -iHpsi
      USE op_m
      USE psi_m

      TYPE(psi_t), intent(in)       :: psi
      TYPE(psi_t), intent(inout)    :: Hpsi
      TYPE(Op_t)                    :: H
      CALL calc_OpPsi(H, psi, Hpsi)

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
      write (out_unit, *) 'Kmax= ', propa%Kmax
      write (out_unit, *) 'propa_name = ', propa%propa_name
      write (out_unit, *) 'propa_name2 = ', propa%propa_name2

   END SUBROUTINE write_propa

   SUBROUTINE Calc_average_energy(Psi, E)
      !>======================================================
      !>     E = <Psi | H | Psi>
      !>======================================================
      USE QDUtil_m
      USE psi_m
      USE Basis_m

      TYPE(psi_t), intent(in)                        :: psi
      TYPE(psi_t)                                    :: Hpsi, psi_b
      REAL(KIND=Rkind), intent(inout)                :: E
      TYPE(Op_t)                                     :: H
      REAL(KIND=Rkind)                               :: Norm
      if (Psi%Grid) then

         !Print*,"psi  is on Grid"
         CALL init_psi(psi_b, psi%Basis, cplx=.TRUE., grid=.false.)
         call GridTOBasis_nD_cplx(psi_b%CVec, psi%CVec, psi%Basis)
         CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.false.)
         call calc_OpPsi(H, psi_b, Hpsi)
         E = real(dot_product(Hpsi%CVec, psi_b%CVec), kind=Rkind)

      else
         !Print*,"psi is on basis"
         CALL init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.false.)
         call calc_OpPsi(H, psi, Hpsi)
         E = real(dot_product(Hpsi%CVec, psi%CVec), kind=Rkind)

      end if
      call Calc_Norm_OF_Psi(psi, Norm)
      E = E/Norm**2
      print *, "<psi|H|psi> = ", E, "<psi|psi> =", Norm

      CALL dealloc_psi(Hpsi)
      CALL dealloc_psi(psi_b)
   End SUBROUTINE Calc_average_energy

   subroutine diff()
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
   end subroutine diff

   subroutine creat_file_unit(nio, name, propa)
      character(*), intent(in)    :: name
      type(propa_t), intent(in)   :: propa
      character(100)              :: name_tot
      integer, intent(in)         :: nio
      character(len=20)           :: dt
      character(len=8)            :: fmt

      fmt = "(E0.1)"

      write (dt, fmt) propa%delta_t
      name_tot = trim(name)//'_'//trim(propa%propa_name)//'_'//trim(propa%propa_name2)//'.dat'
      name_tot = trim(name_tot)

      open (unit=nio, file=name_tot)

   end subroutine


   SUBROUTINE Calc_varia_princinpe_overlap_s(S,V,psi)
      type(psi_t), intent(in)           :: psi
      complex(kind=Rkind),allocatable,intent(inout) :: V(:),S(:,:)
      !------------------locals variables----------------------------------------------------
      integer                            :: nb,nq,ib,jb,n
      type(psi_t)                        :: Hpsi,psi_g,dapsi,dq0psi,dp0psi,dag,dpg,dqg
      complex(kind=Rkind),allocatable    :: Q(:,:),QT(:,:)

      TYPE(Op_t)                         :: H

      nb = psi%Basis%nb
      nq = psi%Basis%nq
      n = nb+3

      call init_psi(Hpsi, psi%Basis, cplx=.TRUE., grid=.false.)  
      call init_psi(dapsi, psi%Basis, cplx=.TRUE., grid=.false.)
      call init_psi(dq0psi, psi%Basis, cplx=.TRUE., grid=.false.)
      call init_psi(dp0psi, psi%Basis, cplx=.TRUE., grid=.false.)

      call init_psi(psi_g, psi%Basis, cplx=.TRUE., grid=.true.)
      call init_psi(dag, psi%Basis, cplx=.TRUE., grid=.true.)
      call init_psi(dpg, psi%Basis, cplx=.TRUE., grid=.true.)
      call init_psi(dqg, psi%Basis, cplx=.TRUE., grid=.true.)

      allocate(V(n),S(n,n),Q(nq,n),QT(n,nq))

      call  calc_OpPsi(H, psi, Hpsi)
      call  CalcDapqpsi(dapsi,dq0psi,dp0psi,psi)
      call  BasisTOGrid_nD_cplx(psi_g%CVec, Hpsi%CVec, psi%Basis)
      call  BasisTOGrid_nD_cplx(dag%CVec, dapsi%CVec, psi%Basis)
      call  BasisTOGrid_nD_cplx(dqg%CVec, dq0psi%CVec, psi%Basis)
      call  BasisTOGrid_nD_cplx(dpg%CVec, dp0psi%CVec, psi%Basis)
       !stop 'ici'
      do ib = 1,nb
         QT(ib,:) = psi%Basis%tab_basis(1)%d0bgw(ib,:)
      end do

      QT(nb+1,:) = dag%CVec(:)*psi%Basis%tab_basis(1)%w(:)
      QT(nb+2,:) = dqg%CVec(:)*psi%Basis%tab_basis(1)%w(:)
      QT(nb+3,:) = dpg%CVec(:)*psi%Basis%tab_basis(1)%w(:)

      do ib = 1,nb
         Q(:,ib) = psi%Basis%tab_basis(1)%d0gb(:,ib)
      end do

      Q(:,nb+1) = dag%CVec(:)
      Q(:,nb+2) = dqg%CVec(:)
      Q(:,nb+3) = dpg%CVec(:)


      do ib = 1,n
         V(ib) = -EYE*dot_product(conjg(QT(ib,:)),psi_g%CVec(:))
      end do
       S = matmul(conjg(QT),Q)

       CALL   Write_VecMat(S, out_unit, 5,  info='S')
       CALL   Write_VecMat(V, out_unit, 5,  info='V')

       !deallocation ------------------------------------------------------
       call dealloc_psi(Hpsi)
       call dealloc_psi(psi_g)
       call dealloc_psi(dapsi)
       call dealloc_psi(dq0psi)
       call dealloc_psi(dp0psi)

       call dealloc_psi(dag)
       call dealloc_psi(dpg)
       call dealloc_psi(dqg)
       deallocate(QT,Q)
      
      END SUBROUTINE

      SUBROUTINE Calc_vp_func(CVec,psi)

         type(psi_t), intent(in)           :: psi
         complex(kind=Rkind),allocatable,intent(inout) :: CVec(:)

         !locals variables-----------------------------------------------------------
         complex(kind=Rkind),allocatable  :: C1Vec(:),CSMat(:,:)!,S(:,:),V(:),res(:)
         integer                          :: n,nb

         nb = psi%Basis%nb
         n = nb+3
         allocate(CVec(n))
         !allocate(V(2),S(2,2),res(2))
         !S(1,1) =TWO;S(1,2) =ONE;S(2,1) =ONE;S(2,2) =-ONE;V(1) =SIX;V(2) =TWO
         !res = LinearSys_Solve(S,V,LS_type=0)
         !print*,'-----------------res----------------',res
         call Calc_varia_princinpe_overlap_s(CSMat,C1Vec,psi)
 
         CVec = LinearSys_Solve(CSMat,C1Vec,LS_type=0)
         deallocate(C1Vec,CSMat)
     END SUBROUTINE

     SUBROUTINE Lambda_rk4(Clambda,psi,propa)

      type(psi_t), intent(in)           :: psi
      type(propa_t), intent(in)         :: propa
      complex(kind=Rkind),allocatable,intent(inout) :: Clambda(:)

        !locals variables-----------------------------------------------------------
      type(psi_t)                        :: Vec,Vec1,Vec2,Vec3,Vec4
      complex(kind=Rkind),allocatable    :: lambda1(:),lambda2(:),lambda3(:),lambda4(:),Clambda0(:)
      real(kind=Rkind)                   :: dt
      integer                            :: n,nb

      nb = psi%Basis%nb
      n = nb+3

      call init_psi(Vec, psi%Basis,  cplx=.TRUE., grid=.false.)  
      call init_psi(Vec1, psi%Basis, cplx=.TRUE., grid=.false.)  
      call init_psi(Vec2, psi%Basis, cplx=.TRUE., grid=.false.)  
      call init_psi(Vec3, psi%Basis, cplx=.TRUE., grid=.false.)  
      call init_psi(Vec4, psi%Basis, cplx=.TRUE., grid=.false.) 

      allocate(Clambda0(n))!(,lambda2(n),lambda3(n),lambda4(n),lambda(n))
      
             dt = propa%delta_t
             !call Calc_vp_func(lambda,Vec)

             Clambda0(1:nb) = psi%CVec(1:nb)
             Clambda0(nb+1) = psi%Basis%tab_basis(1)%scaleQ
             Clambda0(nb+1) = psi%Basis%tab_basis(1)%Q0
             Clambda0(nb+1) = psi%Basis%tab_basis(1)%Imp_k
              
             Vec%CVec(:) = psi%CVec(:)
             call mEyeHPsi(Vec, Vec1)
             call Calc_vp_func(lambda1,Vec1)
             Vec%CVec = psi%CVec + (dt*HALF)*Vec1%CVec
   
             call mEyeHPsi(Vec, Vec2)
             call Calc_vp_func(lambda2,Vec2)
             Vec%CVec = psi%CVec + (dt*HALF)*Vec2%CVec
          
             call mEyeHPsi(Vec, Vec3)
             Vec%CVec = psi%CVec + dt*Vec3%CVec
             call Calc_vp_func(lambda3,Vec3)
            
             call mEyeHPsi(Vec, Vec4)
             call Calc_vp_func(lambda4,Vec4)

             Clambda(:) = Clambda0(:)+ dt*SIXTH*(lambda1(:)+TWO*lambda2(:)+TWO*lambda3(:)+lambda4(:))

             
      call dealloc_psi(Vec)
      call dealloc_psi(Vec1)
      call dealloc_psi(Vec2)
      call dealloc_psi(Vec3)
      call dealloc_psi(Vec4)

      deallocate(lambda1,lambda2,lambda3,lambda4,Clambda0)

     END SUBROUTINE


     SUBROUTINE march_VP(psi, psi_dt, t, propa)
      USE lanczos_m
      USE psi_m
      USE Basis_m
      TYPE(psi_t), INTENT(INOUT)   :: psi_dt
      TYPE(psi_t), INTENT(IN)      :: psi
      TYPE(propa_t), INTENT(IN)    :: propa
      real(kind=Rkind), INTENT(IN)    :: t

        ! variables locales -------------------------------------------------------------------------------

      real(kind=Rkind)                    :: Norm, Norm0,E
       complex (kind=Rkind), allocatable  :: Clambda(:)
      logical, parameter                  :: debug=.false.
      real(kind=Rkind)                    :: SQt(1),Qt(1),Pt(1)
      integer                             :: n,nb

          nb = psi%Basis%nb
          n = nb+3

          IF (debug) THEN
             !write(out_unit,*) 'psi_t',psi%CVec
            flush(out_unit)
          END IF
     
           write (out_unit, *) 'BEGINNIG march VP ', t, propa%delta_t
           allocate(Clambda(n))
           call Lambda_rk4(Clambda,psi,propa)

            psi_dt%CVec(1:nb) = Clambda(1:nb)
            SQt = abs(Clambda(nb+1))
            Qt =  abs(Clambda(nb+2))
            Pt =  abs(Clambda(nb+1))
            call construct_primitive_basis(psi%Basis, x=Qt, sx=SQt,p=Pt) 
            call construct_primitive_basis(psi_dt%Basis, x=Qt, sx=SQt,p=Pt)

            write (out_unit, *) 'lambda',Clambda
     
      CALL Calc_Norm_OF_Psi(psi, Norm0)
      CALL Calc_Norm_OF_Psi(psi_dt, Norm)
      write (out_unit, *) '<psi_dt|psi_dt> = ', Norm, 'abs(<psi_dt|psi_dt> - <psi|psi>)  =', abs(Norm0 - Norm)
      write (out_unit, *) 'END march VP'
      
      
       IF (debug) THEN
         flush(out_unit)
      END IF
   END SUBROUTINE 


   SUBROUTINE H_test(psi)
      TYPE(psi_t), INTENT(INOUT)      :: psi
      TYPE(psi_t)                     :: psi0
      REAL(kind=Rkind), allocatable   :: Qt(:), SQt(:),Pt(:)
      TYPE(Basis_t)                   :: Basis0
      integer                         :: Ndim

   
      call init_Basis1_TO_Basis2(Basis0, psi%Basis)
      call construct_primitive_basis(Basis0)
      call init_psi(psi0, Basis0, cplx=.TRUE., grid=.false.)
     ! allocate (Qt(Ndim), SQt(Ndim),Pt(Ndim))
  

      !Qt(:) = ONE; SQt(:) = ONE; Pt(:) = ONE

      print*,'-------------------------debut du test-------------------------------'
      !print*,'psi in',psi%CVec
      
      call Hagedorn(psi0, psi)
     ! print*,'psi out',psi0%CVec

       print*,'-------------------------fin du test-------------------------------'


      end SUBROUTINE




end module Propa_m
