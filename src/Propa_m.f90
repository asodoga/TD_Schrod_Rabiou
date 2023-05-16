module Propa_m
   USE NumParameters_m
   USE psi_m
   USE Basis_m
   USE Op_m
   USE Ana_psi_m
   Use lanczos_m
   USE Auto_corr_m

   implicit none

   TYPE propa_t
      real(kind=Rk)      :: t0
      real(kind=Rk)      :: tf
      real(kind=Rk)      :: delta_t
      real(kind=Rk)      :: eps
      integer             :: max_iter
      character(len=:), allocatable  :: propa_name
      character(len=:), allocatable  :: propa_name2
   END TYPE propa_t

   public :: march_taylor, marh_RK4th, read_propa
   public :: propagation, Hagedorn,TEST_D,Calc_psi_tild_coeff
   public :: mEyeHPsi, write_propa, Analyse, creat_file_unit

contains

   SUBROUTINE march(psi, psi_dt, t, propa)
      USE psi_m
      TYPE(propa_t), INTENT(IN)               :: propa
      TYPE(psi_t), INTENT(IN)               :: psi
      TYPE(psi_t), INTENT(INOUT)            :: psi_dt
      real(kind=Rk), INTENT(IN)                :: t

      real(kind=Rk)               :: Qt, sQt, Norm, Norm0

      select case (propa%propa_name2)
      case ('rk4')
         CALL marh_RK4th(psi, psi_dt, t, propa)
      case ('taylor')
         CALL march_taylor(psi, psi_dt, t, propa)
      case default
         write (out_unitp, *) 'name is not in the list'
      end select
      CALL Calc_Norm_OF_Psi(psi, Norm0)
      CALL Calc_Norm_OF_Psi(psi_dt, Norm)
      ! write(out_unitp,*) '<psi_dt|psi_dt> = ',Norm , 'abs(<psi_dt|psi_dt> - <psi0|psi0>)  =',abs(Norm0-Norm)

   END SUBROUTINE

   SUBROUTINE propagation(psif, psi0, propa)
      USE psi_m
      USE Basis_m

      TYPE(psi_t), intent(inout)       :: psif
      TYPE(psi_t), intent(in)          :: psi0
      TYPE(propa_t), intent(inout)     :: propa
      logical, parameter               :: debug = .true.
      TYPE(Basis_t), target            :: Basis_1, Basis_2

      ! variables locales
      REAL(kind=Rk)                    :: t, t_deltat, Norm, E, y
      REAL(kind=Rk), allocatable       :: Qt(:), SQt(:), Auto_corr_function(:), populat(:), K(:)
      complex(kind=Rk)                 ::  x
      integer                          :: Ndim

      INTEGER                          :: i, nt, Iq, nf
      TYPE(psi_t)                      :: psi, psi_dt, psi00
      if (debug) then

         write (out_unitp, *) 'BEGINNIG propagation', propa%t0, propa%tf, propa%delta_t
         ! write(out_unitp,*) ''

         write (out_unitp, *) '-------------propagation parameters---------------'
         Call write_propa(propa)
      else
         STOP ' check your data!'
         flush (out_unitp)

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
      call creat_file_unit(nio=19, name='imp_k', propa=propa)

      Ndim = size(psi0%Basis%tab_basis) - 1
      allocate (Qt(Ndim), SQt(Ndim))
      allocate (K(Ndim))
      allocate (populat(psi0%Basis%tab_basis(Ndim + 1)%nb))
      Qt(:) = ZERO; SQt(:) = ONE
      nt = int((propa%tf - propa%t0)/propa%delta_t)
      E = ZERO

      CALL init_Basis1_TO_Basis2(Basis_1, psi0%Basis)
      CALL init_Basis1_TO_Basis2(Basis_2, psi0%Basis)
      CALL construct_primitive_basis(Basis_1)
      CALL construct_primitive_basis(Basis_2)

      CALL init_psi(psi, Basis_1, cplx=.TRUE., grid=.false.)
      CALL init_psi(psi_dt, Basis_2, cplx=.TRUE., grid=.false.)
      CALL init_psi(psi00, psi0%Basis, cplx=.TRUE., grid=.false.)

      psi%CVec(:) = psi0%CVec(:)
      !call write_psi(psi=psi, psi_cplx=.false., print_psi_grid=.true. &
      !               , print_basis=.false., t=ZERO, int_print=20, real_part=.false.)

      ! ---------------------------------- Beging  propagation----------------------------------------------------------
      DO i = 0, nt
         t = i*propa%delta_t
         t_deltat = t + propa%delta_t
         write (out_unitp, *) propa%propa_name2, i, t, t_deltat
         CALL Calc_AVQ_nD0(Psi0=psi, AVQ=Qt, SQ=SQt)
         call Calc_average_energy(psi, E)
         call Calc_Norm_OF_Psi(psi, Norm)
         call Population(psi, populat)
         call Calc_Av_imp_k_nD(psi, K)
         write (11, '(F18.6,2X,F18.6,F18.6,2X,F18.6)') t, Qt
         write (12, '(F18.6,2X,F18.6,F18.6,2X,F18.6)') t, E
         write (13, '(F18.6,2X,F18.6,F18.6,2X,F18.6)') t, SQt
         write (14, '(F18.6,2X,F18.6,F18.6,2X,F18.6)') t, Norm
         write (18, '(F18.6,2X,F18.6,F18.6,2X,F18.6)') t, populat(:)
         write (19, '(F18.6,2X,F18.6,F18.6,2X,F18.6)') t, K
         if (mod(i, 1) == 0) then
            call write_psi(psi=psi, psi_cplx=.false., print_psi_grid=.true. &
                           , print_basis=.false., t=t, int_print=10, real_part=.false.)
            ! write(10,*)
         end if

         CALL march(psi, psi_dt, t, propa)
         if (propa%propa_name == 'hagedorn') Then
            call Hagedorn(psi, psi_dt)
         else
            psi%CVec(:) = psi_dt%CVec(:)
         end if

         if (propa%propa_name == 'hagedorn') Then
            call Calc_Auto_corr(psi0, psi, x, y, propa%propa_name)
            write (15, '(F10.6,2X,F10.6,F10.6,2X,F10.6)') t, abs(x), y
            psi00%CVec = CZERO
            call Hagedorn0(psi00, psi)
            call write_psi(psi=psi00, psi_cplx=.true., print_psi_grid=.false. &
                           , print_basis=.false., t=t, int_print=16, real_part=.true.)
            write (16, *) ''
         else
            call Calc_Auto_corr(psi0, psi_dt, x, y, propa%propa_name)
            write (15, '(F18.6,2X,F18.6,F18.6,2X,F18.6)') t, abs(x), y
            call write_psi(psi=psi, psi_cplx=.true., print_psi_grid=.false. &
                           , print_basis=.false., t=t, int_print=17, real_part=.true.)
            write (17, *) ''
         end if

      END DO
      psif = psi_dt
      CALL Calc_Norm_OF_Psi(psif, Norm)
      IF (debug) THEN
         write (out_unitp, *) 'END propagation'
         write (out_unitp, *) 'norm,psi_dt', Norm

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
         flush (out_unitp)
      END IF

   END SUBROUTINE


   SUBROUTINE Calc_psi_tild_coeff(psi_2,psi_1)

      TYPE(psi_t), intent(inout) ,target:: psi_2
      TYPE(psi_t), intent(in)  ,target  :: psi_1
      TYPE(psi_t),target                :: psi_g1,psi_g2
      logical, parameter                :: debug = .true.
      complex(kind=Rk),pointer          :: GB1(:,:,:),GB2(:,:,:) 
      real(kind=Rk) ,allocatable        :: Pt(:) ,Qt(:),SQt(:),x(:)
      complex(kind=Rk) ,allocatable     :: W(:)
      Integer, allocatable              :: Ib1(:), Ib2(:), Ib3(:),Iq1(:), Iq2(:), Iq3(:)
      integer                           :: nq,inb,Ndim,inbe,i1,i3
      real(kind=Rk)                     :: E,Norm
      
        call Calc_average_energy(psi_1, E)
        Ndim = size(psi_1%Basis%tab_basis)-1
        call init_psi(psi_g1, psi_1%Basis, cplx=.TRUE., grid=.true.)
        call init_psi(psi_g2, psi_1%Basis, cplx=.TRUE., grid=.true.)
        Call Calc_index(Ib1=Ib1, Ib2=Ib2, Ib3=Ib3,Iq1=Iq1, Iq2=Iq2,&
        & Iq3=Iq3, Basis=psi_1%Basis)
        allocate(Pt(Ndim),Qt(Ndim),SQt(Ndim))
        
        CALL Calc_AVQ_nD0(Psi0=psi_1, AVQ=Qt, SQ=SQt)
        call Calc_Av_imp_k_nD(psi_1, Pt)
  
        Call BasisTOGrid_nD_cplx(psi_g1%CVec,psi_1%CVec,psi_1%Basis)

        
        DO inb = 1,Ndim

            allocate(x(psi_1%Basis%tab_basis(inb)%nq),W(psi_1%Basis%tab_basis(inb)%nq))
            x(:) = psi_1%Basis%tab_basis(inb)%x(:)-Qt(inb)
            w(:) = exp(-EYE*Pt(inb)*x(:))
            GB1(1:Iq1(inb),1:Iq2(inb),1:Iq3(inb)) => psi_g1%CVec
            GB2(1:Iq1(inb),1:Iq2(inb),1:Iq3(inb))=> psi_g2%CVec

           DO i3 = 1, ubound(GB1, dim=3)
              DO i1 = 1, ubound(GB1, dim=1)
                 GB2(i1, :, i3) = GB1(i1, :, i3)*W(:)
              END DO
           END DO
           
           psi_g1%CVec(:) = psi_g2%CVec(:)
           deallocate(w,x)
           psi_2%Basis%tab_basis(inb)%imp_k = Pt(inb)
           print*,'pt',psi_2%Basis%tab_basis(Inb)%imp_k 
           If(inb == Ndim) then
            call GridTOBasis_nD_cplx(psi_2%CVec,psi_g2%CVec,psi_1%Basis)
           End if
         
         END DO
         call Calc_average_energy(psi_2, E)
         call Calc_Av_imp_k_nD(psi_2, Pt)
          deallocate(Ib1, Ib2, Ib3,Iq1, Iq2, Iq3)
           call dealloc_psi(psi_g1)
            call dealloc_psi(psi_g2)

   END SUBROUTINE


    SUBROUTINE TEST_D(psi2,psi1)
      TYPE(psi_t), intent(inout)   :: psi2
      TYPE(psi_t), intent(in)      :: psi1
      TYPE(psi_t)                  :: psi_g1,psi_g2
      real(kind=Rk) ,allocatable   :: Pt(:) ,Qt(:),SQt(:),x(:)
      real(kind=Rk)                :: E,Norm
      integer                      :: iq

      print*,'Begin tes'
      print*,'<psi|H|psi>'
      call Calc_average_energy(psi1, E)
      call Calc_Norm_OF_Psi(psi1,Norm)
      print*,'<psi|psi>',Norm
      call init_psi(psi_g1, psi1%Basis, cplx=.TRUE., grid=.true.)
      call init_psi(psi_g2, psi1%Basis, cplx=.TRUE., grid=.true.)
       
      allocate(Pt(1),x(psi1%Basis%nq),Qt(1),SQt(1))
      CALL Calc_AVQ_nD0(Psi0=psi1, AVQ=Qt, SQ=SQt)
      Qt(:)=ZERO;SQt(:)=ONE;x(:)=ZERO
      call Calc_Av_imp_k_nD(psi1, Pt)
      call BasisTOGrid_nD_cplx(psi_g1%CVec,psi1%CVec,psi1%Basis)

      x(:) =psi1%Basis%tab_basis(1)%x(:)-Qt(1)

      psi_g2%CVec(:) = psi_g1%CVec(:)*exp(-EYE*Pt(1)*x(:))

      call GridTOBasis_nD_cplx(psi2%CVec,psi_g2%CVec,psi1%Basis)

      psi2%Basis%tab_basis(1)%imp_k = Pt(1)
      print*,'pt',psi2%Basis%tab_basis(1)%imp_k 

      print*,'<psi_tild|H|psi_tild>'

      call Calc_average_energy(psi2, E)
      call Calc_Norm_OF_Psi(psi2,Norm)
      print*,'<psi_tild|psi_tild>',Norm
      print*,'End tes'

    END SUBROUTINE


   SUBROUTINE Hagedorn(psi, psi_dt)
      USE psi_m
      USE Basis_m

      TYPE(psi_t), intent(inout)                  :: psi, psi_dt
      logical, parameter                          :: debug = .true.

      ! variables locales
      REAL(kind=Rk), allocatable                  :: Qt(:), SQt(:), Pt(:)
      REAL(kind=Rk)                               :: Norm, E
      integer                                     :: Ndim
      write (out_unitp, *) 'Beging Hagedorn'

      call Calc_Norm_OF_Psi(psi_dt,Norm)
       write(out_unitp,*) '<psi_dt|psi_dt> =',Norm

      Ndim = size(psi_dt%Basis%tab_basis) - 1
      allocate (Qt(Ndim), SQt(Ndim), Pt(Ndim))
      Qt(:) = ZERO; SQt(:) = ONE; Pt(:) = ZERO

      call Calc_AVQ_nD0(Psi0=psi_dt, AVQ=Qt, SQ=SQt)
      call Calc_Av_imp_k_nD(psi_dt, Pt)
      call construct_primitive_basis(psi_dt%Basis, Qt, SQt, Pt)
      call projection(psi, psi_dt)
      call construct_primitive_basis(psi%Basis, Qt, SQt, Pt)
      write (out_unitp, *) 'End Hagedorn'

      call Calc_Norm_OF_Psi(psi,Norm)
      write(out_unitp,*) '<psi|psi> =',Norm
      IF (debug) THEN
         flush (out_unitp)
      END IF

   END SUBROUTINE

   SUBROUTINE Analyse(psi, t)
      implicit none
      TYPE(psi_t), INTENT(IN)       :: psi
      real(Kind=Rk), allocatable     :: pop(:), Qm(:), Qp(:)
      real(kind=Rk), intent(in)        :: t
      real(kind=Rk)                   :: Norm, E
      integer                         :: Ndim
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

      TYPE(psi_t), INTENT(INOUT):: psi_dt
      TYPE(psi_t), INTENT(IN)    :: psi
      TYPE(psi_t)                 :: psi0
      TYPE(psi_t)                 :: Hpsi
      TYPE(propa_t), INTENT(IN)    :: propa
      real(kind=Rk), INTENT(IN)   :: t
      real(kind=Rk)               :: alpha

      ! variables locales

      real(kind=Rk)                :: Rkk, Norm, Norm0
      integer                      :: kk
      CALL init_psi(Hpsi, psi%basis, cplx=.TRUE., grid=.false.) ! to be changed
      CALL init_psi(Psi0, psi%basis, cplx=.TRUE., grid=.false.) ! to be changed

      write (out_unitp, *) 'BEGINNIG march_taylor  ', t, propa%delta_t
      !write(out_unitp,*) 'psi',psi%CVec
      Rkk = ONE
      alpha = TEN**10
        !!======================debut ordre 1==========================
      !psi_dt%CVec    = psi%CVec
      !CALL calc_OpPsi(H,psi,Hpsi)
      !Rkk = Rkk*delta_t
      !Hpsi%CVec(:)    = - EYE*Hpsi%CVec(:)
      !psi_dt%CVec(:) = psi_dt%CVec(:) +Rkk*Hpsi%CVec(:)
      !CALL Calc_Norm_OF_Psi(HPsi,Norm)
      !CALL Calc_Norm_OF_Psi(Psi_dt,Norm)

      !write(out_unitp,*) 'norm,Hpsi',Rkk*Norm
      !write(out_unitp,*) 'norm,psi_dt',Norm

      !write(out_unitp,*) 'psi_dt',psi_dt%CVec
      !write(out_unitp,*) 'END march_taylor'

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
         write (out_unitp, *) 'sqrt(<Hpsi|Hpsi>) = ', kk, Norm
         if (Norm >= alpha) then
            stop "wrong choice of delta_t"
         elseif (Norm <= propa%eps) Then

            print *, 'Taylor condition is fulfild after', kk, 'iteration'
            exit
         End if
      End do
      CALL Calc_Norm_OF_Psi(Psi, Norm0)
      CALL Calc_Norm_OF_Psi(Psi_dt, Norm)
      write (out_unitp, *) '<psi_dt|psi_dt> = ', Norm, 'abs(<psi_dt|psi_dt> - <psi0|psi0>)  =', abs(Norm0 - Norm)
      write (out_unitp, *) 'END march_taylor'
      !CALL dealloc_psi(psi0)
      !CALL dealloc_psi(Hpsi)
   END SUBROUTINE march_taylor

   SUBROUTINE marh_RK4th(psi, psi_dt, t, propa)
      USE op_m
      USE psi_m
      USE Basis_m

      TYPE(psi_t), intent(inout)     :: psi_dt
      TYPE(psi_t), intent(in)        :: psi
      TYPE(psi_t)                     :: K1, K2, K3, K4, psi_inter
      TYPE(propa_t), intent(in)        :: propa

      real(kind=Rk), intent(in)        :: t
      real(kind=Rk)                    ::  Norm, Norm0
      integer                          :: iq
      !  variables locales

      call init_psi(K1, psi%basis, cplx=.true., grid=.false.)
      call init_psi(K2, psi%basis, cplx=.true., grid=.false.)
      call init_psi(K3, psi%basis, cplx=.true., grid=.false.)
      call init_psi(K4, psi%basis, cplx=.true., grid=.false.)
      call init_psi(psi_inter, psi%basis, cplx=.true., grid=.false.)
      write (out_unitp, *) 'BEGINNIG march_RK4th', t, propa%delta_t
      ! write(out_unitp,*) 'psi'
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

      write (out_unitp, *) '<psi|psi> = ', Norm0, '<psi_dt|psi_dt> = ', Norm, 'abs(<psi|psi> - <psi_dt|psi_dt))=', ABS(Norm0 - Norm)
      write (out_unitp, *) 'END marh_RK4th'
      deallocate (k1%CVec)
      deallocate (k2%CVec)
      deallocate (k3%CVec)
      deallocate (k4%CVec)
      deallocate (psi_inter%CVec)
   END SUBROUTINE marh_RK4th

   SUBROUTINE read_propa(propa)
      USE psi_m
      implicit none
      TYPE(propa_t), intent(inout) :: propa
      real(kind=Rk)                 :: t0, tf, delta_t, eps
      character(len=40)             :: propa_name, propa_name2
      integer                        ::  max_iter

      namelist /prop/ t0, tf, delta_t, eps, max_iter, propa_name, propa_name2
      t0 = ZERO
      tf = 10._Rk
      delta_t = 0.001
      eps = ONETENTH**10
      max_iter = 5000
      propa_name = 'non_hagedorn'
      propa_name2 = 'rk4'

      read (*, nml=prop)

      propa%t0 = t0
      propa%tf = tf
      propa%delta_t = delta_t
      propa%eps = eps
      propa%max_iter = max_iter
      propa%propa_name = propa_name
      propa%propa_name2 = propa_name2

   END SUBROUTINE read_propa
   SUBROUTINE mEyeHPsi(Psi, HPsi) !calcul de -iHpsi
      USE op_m
      USE psi_m

      TYPE(psi_t), intent(in)      :: Psi
      TYPE(psi_t), intent(inout)   :: HPsi
      TYPE(Op_t)                    :: H
      CALL calc_OpPsi(H, Psi, HPsi)

      HPsi%CVec(:) = -EYE*HPsi%CVec(:)
   END SUBROUTINE mEyeHPsi

   SUBROUTINE write_propa(propa)
      USE psi_m
      implicit none
      TYPE(propa_t), intent(inout) :: propa

      write (out_unitp, *) 't0 = ', propa%t0
      write (out_unitp, *) 'tf = ', propa%tf
      write (out_unitp, *) 'deltat_t = ', propa%delta_t
      write (out_unitp, *) 'eps = ', propa%eps
      write (out_unitp, *) 'max_iter = ', propa%max_iter
      write (out_unitp, *) 'propa_name = ', propa%propa_name
      write (out_unitp, *) 'propa_name2 = ', propa%propa_name2

   END SUBROUTINE write_propa

   SUBROUTINE Calc_average_energy(Psi, E)
      !>======================================================
      !>     E = <Psi | H | Psi>
      !>======================================================
      USE UtilLib_m
      USE psi_m
      USE Basis_m

      TYPE(psi_t), intent(in)                      :: Psi
      TYPE(psi_t)                                  :: HPsi, Psi_b
      REAL(KIND=Rk), intent(inout)                 :: E
      TYPE(Op_t)                                   :: H
      REAL(KIND=Rk)                                :: Norm

      
      if (Psi%Grid) then

         !Print*,"psi  is on Grid"
         CALL init_psi(Psi_b, Psi%Basis, cplx=.TRUE., grid=.false.)
         call GridTOBasis_nD_cplx(Psi_b%CVec, Psi%CVec, Psi%Basis)
         CALL init_psi(Hpsi, Psi%Basis, cplx=.TRUE., grid=.false.)
         call calc_OpPsi(H, Psi_b, Hpsi)
         E = real(dot_product(Hpsi%CVec, Psi_b%CVec), kind=Rk)

      else
         !Print*,"psi is on basis"
         CALL init_psi(Hpsi, Psi%Basis, cplx=.TRUE., grid=.false.)

         call calc_OpPsi(H, Psi, Hpsi)
         E = real(dot_product(Hpsi%CVec, Psi%CVec), kind=Rk)

      end if
      call Calc_Norm_OF_Psi(Psi, Norm)
      E = E/Norm**2
      print *, "<Psi|H|Psi> = ", E, "Norm=", Norm
      !stop 'cc energie'
      CALL dealloc_psi(HPsi)
      CALL dealloc_psi(Psi_b)
   End SUBROUTINE Calc_average_energy

   subroutine diff()
      real(kind=Rk), allocatable         :: df(:, :)
      real(kind=Rk), allocatable         :: f1(:, :)
      real(kind=Rk), allocatable         :: f2(:, :)
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

end module Propa_m
