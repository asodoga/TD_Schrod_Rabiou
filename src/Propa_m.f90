module Propa_m
     USE NumParameters_m
     USE psi_m
     USE Basis_m
     USE Op_m
     USE Ana_psi_m
     Use lanczos_m

     implicit none

     TYPE propa_t
        real (kind=Rk)      :: t0
        real (kind=Rk)      :: tf
        real (kind=Rk)      :: delta_t
        real (kind=Rk)      :: eps
        integer             :: max_iter
        integer             :: kmax
        character(len=:), allocatable  :: propa_name
    END TYPE propa_t

    public :: march_taylor,marh_RK4th,read_propa,autocor_func
    public :: propagation_Hagedorn,propagation_DP
    public :: mEyeHPsi,write_propa,fft_autocor_func,Analyse,Hagedorn_psi1_TO_psi2

contains
    SUBROUTINE propagation_Hagedorn(psif,psi0,Basis_f,propa)
        USE psi_m
        USE Basis_m

        TYPE (psi_t),  intent(inout)     :: psif
        TYPE (psi_t),  intent(in)        :: psi0
        TYPE(propa_t), intent(inout)     :: propa
        TYPE(Basis_t) ,intent(inout)     :: Basis_f
        TYPE(Basis_t) ,target            :: Basis_1,Basis_2
        logical, parameter               :: debug = .true.
        REAL(kind=Rk)                    :: Qt,SQt,Q0,SQ0

        ! variables locales
        REAL(kind=Rk)                    :: t ,t_deltat, Norm,E

        INTEGER                          :: i,nt,Iq,nf
        TYPE (psi_t)                     :: psi,psi_dt
        if (debug) then

            write(out_unitp,*) 'BEGINNIG propagation', propa%t0,propa%tf,propa%delta_t
           ! write(out_unitp,*) ''

            write(out_unitp,*) '-------------propagation parameters---------------'
            Call write_propa(propa)
        else
            STOP ' check your data!'
            flush(out_unitp)

        endif
        nt = int((propa%tf-propa%t0)/propa%delta_t)
        open(10, file = '<xHG>.dat')
        !@@@@@@@@@@@@@@@ initialization of psi and psi_dt and construction of Basis_loc@@@@@@@@@@@@@@@
        CALL init_Basis1_TO_Basis2 (psi0%Basis,Basis_1)
        CALL init_Basis1_TO_Basis2 (psi0%Basis,Basis_2)
        CALL construct_primitive_basis(Basis_1)
        CALL construct_primitive_basis(Basis_2)
        CALL init_psi(psi,Basis_1,cplx=.TRUE.,grid =.false.)
        CALL init_psi(psi_dt,Basis_2,cplx=.TRUE.,grid =.false.  )

        ! ************************************** initialisation of psi and <Q> ,DQ *************************************
       ! Q0 = ZERO ; SQ0 = ZERO
        psi%CVec(:) = psi0%CVec(:)
        ! ******************************* Beging  propagation ********************************************************
        DO i=0,nt
            t = i*propa%delta_t
            t_deltat = t + propa%delta_t
            write(out_unitp,*) propa%propa_name,i,t,t_deltat
           ! call Analyse(psi,t)
            Write(10,*)   t, Qt
            if(mod(i,25)== 0)  call Write_p(psi,i)
            CALL  march(psi,psi_dt,t,propa)
            CALL  Calc_std_dev_AVQ_1D(psi_dt,1,Qt,SQt)
           ! Qt = Qt+Q0
            CALL  Calc_basis(Basis_1, psi0%Basis,Qt,ONE)
            CALL Hagedorn_transfo(psi_dt,psi)
            CALL  Calc_basis(Basis_2, psi0%Basis,Qt,ONE)
            ! Q0 = Qt

        END DO
        psif%CVec(:) = psi_dt%CVec(:)
        CALL  Calc_basis(Basis_f, psi0%Basis,Qt,ONE)
        call Calc_Norm_OF_Psi(psif, Norm)

        IF (debug) THEN
            write(out_unitp,*) 'END propagation'
            write(out_unitp,*) 'norm,psi_dt',Norm
            call write_psi(psif)

            flush(out_unitp)
        END IF

    END SUBROUTINE propagation_Hagedorn


    SUBROUTINE propagation_DP(psif,psi0,propa)
        USE psi_m
        USE Basis_m

        TYPE (psi_t),  intent(inout)     :: psif
        TYPE (psi_t),  intent(in)        :: psi0
        TYPE(propa_t), intent(inout)     :: propa
        logical, parameter               :: debug = .true.

        ! variables locales
        REAL(kind=Rk)                    :: t ,t_deltat, Norm,E,Qt,SQt

        INTEGER                          :: i,nt,Iq,nf
        TYPE (psi_t)                     :: psi,psi_dt
        if (debug) then

            write(out_unitp,*) 'BEGINNIG propagation', propa%t0,propa%tf,propa%delta_t
            ! write(out_unitp,*) ''

            write(out_unitp,*) '-------------propagation parameters---------------'
            Call write_propa(propa)
        else
            STOP ' check your data!'
            flush(out_unitp)

        endif
        nt = int((propa%tf-propa%t0)/propa%delta_t)
        open(12, file = '<x>.dat')

        !@@@@@@@@@@@@@@@ initialization of psi and psi_dt and construction of Basis_loc@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        CALL init_psi(psi,psi0%Basis,cplx=.TRUE.,grid =.false.)
        CALL init_psi(psi_dt,psi0%Basis,cplx=.TRUE.,grid =.false.  )
        ! ************************************** initialisation of psi and <Q> ,DQ *************************************
        psi%CVec(:) = psi0%CVec(:)
        ! ******************************* Beging  propagation ********************************************************
        DO i=0,nt
            t = i*propa%delta_t
            t_deltat = t + propa%delta_t
            write(out_unitp,*) propa%propa_name,i,t,t_deltat
            ! call Analyse(psi,t)
            CALL  Calc_std_dev_AVQ_1D(psi,1,Qt,SQt)
            Write(12,*)   t, Qt
            if(mod(i,25)== 0)  call Write_p(psi,i)
            CALL  march(psi,psi_dt,t,propa)
            psi%CVec(:) = psi_dt%CVec(:)

        END DO
        psif%CVec(:) = psi_dt%CVec(:)
        call Calc_Norm_OF_Psi(psif, Norm)
        IF (debug) THEN
            write(out_unitp,*) 'END propagation'
            write(out_unitp,*) 'norm,psi_dt',Norm
            call write_psi(psif)

            flush(out_unitp)
        END IF

    END SUBROUTINE propagation_DP

    SUBROUTINE Analyse(psi,t)
        implicit none
        TYPE (psi_t)  ,INTENT(IN)       :: psi
        real(Kind = Rk),allocatable     :: pop(:),Qm(:),Qp(:)
        real(kind=Rk),intent(in)        :: t
        real(kind=Rk)                   :: Norm,E
        integer                         :: Ndim
        Ndim = size(psi%Basis%tab_basis)
        allocate(Pop(psi%Basis%tab_basis(Ndim)%nb))
        allocate(Qp(psi%Basis%tab_basis(Ndim)%nb))
        allocate(Qm(Ndim-1))
        pop(:) = ZERO
        Qm(:)=  ZERO
        Qp(:) = ZERO
        E = ZERO;Norm = ZERO
        !===================================== beging Anapsi==================
        call Population(psi,pop)
        call Calc_average_energy(psi,E)
        call Average_Q(psi,Qm)
        call Calc_Norm_OF_psi(psi, Norm)
        call Qpop(psi,Qp)
        write(3,*) t,E,Norm,pop
        write(4,*) t,Qm
        write(5,*) t,Qp
        deallocate(pop,Qm,Qp)
        !====================================And Anapsi========================


     END SUBROUTINE Analyse


    SUBROUTINE Hagedorn_psi1_TO_psi2(psi_dt,Basis)
        USE psi_m
        USE Basis_m

        TYPE (psi_t),  intent(inout)     :: psi_dt
        TYPE (psi_t)                     :: psi
        TYPE(Basis_t)                    :: Basis

        CALL init_psi(psi,Basis,cplx=.TRUE.,grid =.false.  )
        CALL Hagedorn_transfo(psi_dt,psi)
             psi_dt%CVec(:) = psi%CVec(:)
    END SUBROUTINE Hagedorn_psi1_TO_psi2


    SUBROUTINE march_taylor(psi,psi_dt,t,propa)
       USE op_m
       USE psi_m
       USE Basis_m

       TYPE (psi_t)  , INTENT(INOUT):: psi_dt
       TYPE (psi_t)  ,INTENT(IN)    :: psi
       TYPE (psi_t)                 :: psi0
       TYPE (psi_t)                 :: Hpsi
       TYPE(propa_t) ,INTENT(IN)    :: propa
       real(kind=Rk) , INTENT(IN)   :: t
        real(kind=Rk)               :: alpha

       ! variables locales

       real(kind=Rk)                :: Rkk, Norm,Norm0
       integer                      :: kk
     CALL init_psi(Hpsi ,psi%basis,  cplx=.TRUE.   ,grid =.false. ) ! to be changed
     CALL init_psi(Psi0 ,psi%basis,  cplx=.TRUE.   ,grid =.false. ) ! to be changed


       write(out_unitp,*) 'BEGINNIG march_taylor  ',t,propa%delta_t
       !write(out_unitp,*) 'psi',psi%CVec
       Rkk = ONE
       alpha = 10._Rk**15
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
        Psi_dt%CVec    = Psi%CVec
       Psi0%CVec    = Psi%CVec
        Do  kk = 1,propa%max_iter,1
            CALL mEyeHPsi(psi0,Hpsi)
            Rkk = Rkk*(propa%delta_t/kk)
            psi_dt%CVec(:) = psi_dt%CVec(:) + Rkk*Hpsi%CVec(:)
            psi0%CVec(:) =Hpsi%CVec(:)
            Hpsi%CVec(:) = CZERO
           call Calc_Norm_OF_Psi(psi0,Norm)
            Norm =   Rkk*Norm
            write(out_unitp,*) '<Hpsi|Hpsi> = ',kk,Norm
            if (Norm >= alpha ) then
                stop "wrong choice of delta_t"
            elseif(Norm <= propa%eps) Then

                print*,'Taylor condition is fulfild after',kk,'iteration'
                exit
            Endif
        Enddo
       CALL Calc_Norm_OF_Psi(Psi,Norm0)
       CALL Calc_Norm_OF_Psi(Psi_dt,Norm)
       write(out_unitp,*) '<psi_dt|psi_dt> = ',Norm , 'abs(<psi_dt|psi_dt> - <psi0|psi0>)  =',abs(Norm0-Norm)
       write(out_unitp,*) 'END march_taylor'
       !CALL dealloc_psi(psi0)
       !CALL dealloc_psi(Hpsi)
    END SUBROUTINE march_taylor


    SUBROUTINE marh_RK4th(psi,psi_dt,t,propa)
       USE op_m
       USE psi_m
       USE Basis_m

       TYPE (psi_t),  intent(inout)     :: psi_dt
       TYPE (psi_t),  intent(in)        :: psi
       TYPE (psi_t)                     :: K1,K2,K3,K4 ,psi_inter
       TYPE(propa_t), intent(in)        :: propa

       real(kind=Rk), intent(in)        :: t
       real(kind=Rk)                    ::  Norm, Norm0
       integer                          :: iq
       !  variables locales

       call init_psi(K1,psi%basis,cplx = .true.,grid = .false.)
       call init_psi(K2,psi%basis,cplx = .true.,grid = .false.)
       call init_psi(K3,psi%basis,cplx = .true.,grid = .false.)
       call init_psi(K4,psi%basis,cplx = .true.,grid = .false.)
       call init_psi(psi_inter,psi%basis,cplx = .true.,grid = .false.)
       write(out_unitp,*) 'BEGINNIG march_RK4th',t,propa%delta_t
      ! write(out_unitp,*) 'psi'
       psi_dt%CVec(:)    = psi%CVec(:)
       CALL mEyeHPsi(psi,K1)

       psi_inter%CVec = psi%CVec+(propa%delta_t*HALF)*K1%CVec
       !CALL Write_psi_basis(psi_inter,t,10)
       CALL mEyeHPsi(psi_inter,K2)
       psi_inter%CVec = psi%CVec+(propa%delta_t*HALF)*K2%CVec
       !CALL Write_psi_basis(psi_inter,t,11)
       CALL mEyeHPsi(psi_inter,K3)
       psi_inter%CVec = psi%CVec+propa%delta_t*K3%CVec
       !CALL Write_psi_basis(psi_inter,t,12)
       CALL mEyeHPsi(psi_inter,K4)
       psi_dt%CVec(:) = psi_dt%CVec(:)+(propa%delta_t*SIXTH)*(K1%CVec(:)+ TWO*K2%CVec(:)+ TWO*K3%CVec(:)+K4%CVec(:))
       CALL Calc_Norm_OF_Psi(Psi_dt,Norm)
       CALL Calc_Norm_OF_Psi(Psi,Norm0)

       write(out_unitp,*) '<psi|psi> = ',Norm0 , '<psi_dt|psi_dt> = ',Norm ,'abs(<psi|psi> - <psi_dt|psi_dt))=',ABS(Norm0-Norm)
       write(out_unitp,*) 'END marh_RK4th'
      deallocate(k1%CVec)
      deallocate(k2%CVec)
      deallocate(k3%CVec)
      deallocate(k4%CVec)
      deallocate(psi_inter%CVec)
    END SUBROUTINE marh_RK4th
    subroutine march_SIL(psi,psi_dt,propa)
        implicit none
        TYPE (psi_t),  intent(inout)     :: psi_dt
        TYPE (psi_t),  intent(in)        :: psi
        TYPE(propa_t), intent(in)        :: propa
        real(Kind = Rk)                  :: Norm,Norm0
         write(out_unitp,*) 'BEGING march_SIL'
        call Lan_eign_syst_solve(psi,psi_dt,propa%kmax ,propa%delta_t)
        CALL Calc_Norm_OF_Psi(Psi,Norm0)
        CALL Calc_Norm_OF_Psi(Psi_dt,Norm)
        write(out_unitp,*) '<psi_dt|psi_dt> = ',Norm , 'abs(<psi_dt|psi_dt> - <psi0|psi0>)  =',abs(Norm0-Norm)
        write(out_unitp,*) 'END march_SIL'

    end subroutine march_SIL




    SUBROUTINE read_propa( propa)
      USE psi_m
      implicit none
      TYPE (propa_t),  intent(inout) :: propa
      real(kind= Rk)                 :: t0,tf,delta_t ,eps
      character(len= 40)             :: propa_name
      integer                        ::  max_iter,kmax

      namelist /prop/ t0,tf,delta_t ,eps,max_iter, propa_name,kmax
      t0  = ZERO
      tf  = 10._Rk
      delta_t = 0.001
      eps= ONETENTH**10
      max_iter = 5000
      propa_name = 'rk4'
      kmax       = 50

      read(*,nml=prop)

      propa%t0 = t0
      propa%tf = tf
      propa%delta_t = delta_t
      propa%eps = eps
      propa%max_iter = max_iter
      propa%propa_name = propa_name
      propa%kmax = kmax

    END SUBROUTINE read_propa
    SUBROUTINE mEyeHPsi (Psi,HPsi) !calcul de -iHpsi
       USE op_m
       USE psi_m

       TYPE (psi_t),  intent(in)      :: Psi
       TYPE (psi_t),  intent(inout)   :: HPsi
       TYPE (Op_t)                    :: H
       CALL calc_OpPsi(H,Psi,HPsi)

       HPsi%CVec(:)    = - EYE*HPsi%CVec(:)
    END SUBROUTINE mEyeHPsi

    SUBROUTINE march(psi,psi_dt,t,propa)
       USE psi_m
       TYPE (propa_t) , INTENT(IN)      :: propa
       TYPE (psi_t)   , INTENT(IN)      :: psi
       TYPE (psi_t)   , INTENT(INOUT)   :: psi_dt
       real(kind=Rk)  ,INTENT(IN)       :: t

       select case (propa%propa_name)
            case ('rk4')
            CALL marh_RK4th(psi,psi_dt,t,propa)
            case ('taylor')
            CALL  march_taylor(psi,psi_dt,t,propa)
            case ('SIL')
           CALL  march_SIL(psi,psi_dt,propa)
            case default
            write(out_unitp,*) 'name is not in the list'
       end select

    END SUBROUTINE march




    SUBROUTINE write_propa( propa)
      USE psi_m
      implicit none
      TYPE (propa_t),  intent(inout) :: propa

      write(out_unitp,*) 't0 = ',  propa%t0
      write(out_unitp,*) 'tf = ',propa%tf
      write(out_unitp,*) 'deltat_t = ',propa%delta_t
      write(out_unitp,*) 'eps = ',propa%eps
      write(out_unitp,*) 'max_iter = ',propa%max_iter
      write(out_unitp,*) 'propa_name = ',propa%propa_name
      write(out_unitp,*) 'kmax = ',      propa%kmax

    END SUBROUTINE write_propa

    SUBROUTINE Calc_average_energy(Psi,E)
        !>======================================================
        !>     E = <Psi | H | Psi>
        !>======================================================
        USE UtilLib_m
        USE psi_m
        USE Basis_m

        TYPE (psi_t),  intent(in)                    :: Psi
        TYPE (psi_t)                                 :: HPsi ,Psi_b
        REAL(KIND= Rk) ,intent(inout)                :: E
        TYPE(Op_t)                                   :: H
        REAL(KIND= Rk)                               :: Norm
       if(Psi%Grid) then

           !Print*,"psi  is on Grid"
           CALL init_psi(Psi_b,Psi%Basis,  cplx=.TRUE.,grid =.false.)
           call GridTOBasis_nD_cplx(Psi_b%CVec,Psi%CVec,Psi%Basis)
           CALL init_psi(Hpsi,Psi%Basis,  cplx=.TRUE.,grid =.false.)
           call calc_OpPsi(H,Psi_b,Hpsi)
           E = real(dot_product(Hpsi%CVec,Psi_b%CVec),kind=Rk)

       else
           !Print*,"psi is on basis"
           CALL init_psi(Hpsi,Psi%Basis,  cplx=.TRUE.,grid =.false.)
           call calc_OpPsi(H,Psi,Hpsi)
           E = real(dot_product(Hpsi%CVec,Psi%CVec),kind=Rk)

       end if
        call Calc_Norm_OF_Psi(Psi,Norm)
        E = E/Norm**2
        print*,"<Psi|H|Psi> = ",E,"Norm=",Norm

        CALL dealloc_psi(HPsi)
        CALL dealloc_psi(Psi_b)
    End SUBROUTINE Calc_average_energy

      SUBROUTINE autocor_func(psi,psi_dt,corre_coeff,arg_corre_coeff)
        USE UtilLib_m
        USE op_m
        USE psi_m

        TYPE (psi_t),  intent(in)                :: psi,psi_dt
        COMPLEX(KIND=Rk)  ,  INTENT(INOUT)       :: corre_coeff
        REAL(KIND=Rk)     ,INTENT(INOUT)         ::arg_corre_coeff

        corre_coeff     = dot_product(psi%CVec,psi_dt%CVec)
        arg_corre_coeff =ATAN(AIMAG(corre_coeff )/REAL(corre_coeff ))
      End SUBROUTINE autocor_func

      SUBROUTINE fft_autocor_func(autocor_function,time,fft_autocor_function,delta_t,N)
        USE NumParameters_m
        USE UtilLib_m
        COMPLEX(KIND=Rk),INTENT(IN) ,allocatable, DIMENSION(:)        :: autocor_function(:)
        COMPLEX(KIND=Rk),INTENT(INOUT) , ALLOCATABLE                  :: fft_autocor_function(:)
        REAL(KIND=Rk),INTENT(IN) , DIMENSION(:)                       :: time
        REAL(KIND=Rk),ALLOCATABLE, DIMENSION(:)                       :: w
        INTEGER ,intent(in)                                           :: N
         REAL(KIND=Rk)                                                :: delta_t
        REAL(KIND=Rk)                                                 :: wm,wmax,dw
        INTEGER                                                       :: Iw,nw,I
        OPEN(UNIT=100,FILE="fft_autocor_function.dat")
        fft_autocor_function(:)= CZERO
        wm = -0.001_Rk ; wmax = 5._Rk ;dw = 0.005_Rk
        nw = int((wmax-wm)/dw)
         !allocate( fft_autocor_function(0:nw))
        ALLOCATE(w(0:nw-1))
        DO Iw = 0,nw-1
           w(Iw) = wm + float(Iw)*dw
           fft_autocor_function(Iw)= ZERO
           DO I = 1,N
               fft_autocor_function(Iw) = fft_autocor_function(Iw)+ autocor_function(Iw)*EXP(EYE*w(Iw)*time(I))*delta_t
           ENDDO !itime
           WRITE(100,*)  w(Iw), ABS(fft_autocor_function(Iw))
        ENDDO !omega
      END SUBROUTINE fft_autocor_func

end module Propa_m
