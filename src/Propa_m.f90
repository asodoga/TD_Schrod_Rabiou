module Propa_m
     USE NumParameters_m
     USE psi_m
     USE Basis_m
     USE Op_m

     implicit none

     TYPE propa_t
        real (kind=Rk)      :: t0
        real (kind=Rk)      :: tf
        real (kind=Rk)      :: delta_t
        real (kind=Rk)      :: eps
        integer             :: max_iter
        character(len=:), allocatable  :: propa_name
    END TYPE propa_t

    public :: propagation,march_taylor,marh_RK4th,read_propa,autocor_func
    public :: mEyeHPsi,write_propa,ana_wp,spectrum

contains
    SUBROUTINE propagation(psif,psi0,propa,Basis)
       USE op_m
       USE psi_m
       USE Basis_m

       TYPE (psi_t),  intent(inout) :: psif
       TYPE (psi_t),  intent(in)    :: psi0
       TYPE(Basis_t),    intent(in) :: Basis
       TYPE (psi_t)                 :: G,rho_num ,rho



       TYPE(propa_t), intent(inout) :: propa
       logical, parameter           :: debug = .true.

       ! variables locales
       REAL(kind=Rk)                :: t ,t_deltat, Norm

       INTEGER                      :: i,nt,IQ,nf
       TYPE (psi_t)                 :: psi,psi_dt,rho_ana
       if (debug) then

                 write(out_unitp,*) 'BEGINNIG propagation', propa%t0,propa%tf,propa%delta_t
                 write(out_unitp,*) ''

                 write(out_unitp,*) '-------------propagation parameters---------------'
                 Call write_propa(propa)
       else
                 STOP ' check your data!'
                 flush(out_unitp)

       endif

       nt = int((propa%tf-propa%t0)/propa%delta_t)

       CALL init_psi(psi,psi0%Basis,cplx=.TRUE.,grid =.true.) ! to be changed
       CALL init_psi(psi_dt,psi0%Basis,cplx=.TRUE.,grid =.true.) ! to be changed
       CALL init_psi(G,psi0%Basis,cplx=.TRUE.,grid =.true.) ! to be changed
       CALL init_psi(rho_num,psi0%Basis,cplx=.TRUE.,grid =.false.) ! to be changed
       CALL init_psi(rho,psi0%Basis,cplx=.TRUE.,grid =.true.) ! to be changed
       CALL init_psi(rho_ana,psi0%Basis,cplx=.TRUE.,grid =.true.) ! to be changed


        psi%CVec = psi0%CVec

       DO i=0,nt-1

            t = i*propa%delta_t
            t_deltat = t + propa%delta_t
      write(out_unitp,*) propa%propa_name,i,t,t_deltat

            CALL march(psi,psi_dt,t,propa,Basis)

            psi%CVec(:) = psi_dt%CVec(:)

               nf = int(nt/5)
             IF( nf == 0)then
             nf = 1

             ENDIF

            IF(   MOD(i,nf) == 0  )Then
             OPEN(unit=i+10)
             !CALL BasisTOGrid_cplx(G%CVec,psi_dt%CVec,psi_dt%Basis)
             !rho_num%CVec(:)= G%CVec(:)
             !CALL ana_wp(rho_ana,t)
             !CALL calc_rho(rho_ana,rho_num,rho)
             DO  IQ = 1,Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nb
            WRITE(i+10,*) iQ, ABS(psi_dt%CVec(IQ))**2
             ENDDO
              CLOSE(UNIT=i+10)
            ENDIF


       END DO
       psif%CVec(:) = psi%CVec(:)
       CALL Calc_Norm_Grid(psi_dt%CVec, Norm,Basis)

       CALL dealloc_psi(psi)
       CALL dealloc_psi(psi_dt)
       IF (debug) THEN
           write(out_unitp,*) 'END propagation'
           write(out_unitp,*) 'norm,psi_dt',Norm
           write(out_unitp,*) 'psi_dt',psi_dt%CVec
           flush(out_unitp)
       END IF




    END SUBROUTINE propagation

    SUBROUTINE march_taylor(psi,psi_dt,t,propa,Basis)
       USE op_m
       USE psi_m
       USE Basis_m

       TYPE (psi_t)  , INTENT(INOUT):: psi_dt
       TYPE (psi_t)  ,INTENT(INOUT) :: psi
       TYPE(Basis_t)    ,INTENT(IN) :: Basis
       TYPE (psi_t)                 :: Hpsi
       TYPE(propa_t) ,INTENT(IN)    :: propa
       real(kind=Rk) , INTENT(IN)   :: t

       ! variables locales


       real(kind=Rk)                :: Rkk, Norm
       integer                      :: kk
     CALL init_psi(Hpsi ,Basis,  cplx=.TRUE.   ,grid =.true. ) ! to be changed



       write(out_unitp,*) 'BEGINNIG march_taylor  ',t,propa%delta_t
       !write(out_unitp,*) 'psi',psi%CVec
       Rkk = ONE
       !!======================debut ordre 1==========================
       !psi_dt%CVec    = psi%CVec
       !CALL calc_OpPsi(H,psi,Hpsi)
       !Rkk = Rkk*delta_t
       !Hpsi%CVec(:)    = - EYE*Hpsi%CVec(:)
       !psi_dt%CVec(:) = psi_dt%CVec(:) +Rkk*Hpsi%CVec(:)
       !CALL Calc_Norm(Hpsi, Norm)
       !CALL Calc_Norm(psi_dt, Norm)

       !write(out_unitp,*) 'norm,Hpsi',Rkk*Norm
       !write(out_unitp,*) 'norm,psi_dt',Norm



       !write(out_unitp,*) 'psi_dt',psi_dt%CVec
       !write(out_unitp,*) 'END march_taylor'

       !!===========================Ordre deux etplus=======================

        psi_dt%CVec    = psi%CVec
        Do  kk = 1,propa%max_iter,1
            CALL mEyeHPsi(psi,Hpsi,Basis)
            Rkk = Rkk*(propa%delta_t/kk)
            psi_dt%CVec(:) = psi_dt%CVec(:) +Rkk*Hpsi%CVec(:)

            !CALL Calc_Norm(psi_dt, Norm)
            !write(out_unitp,*) 'norm,psi_dt',Norm
            CALL Calc_Norm_Grid(Hpsi%CVec,Norm,Basis)
            Norm =   Rkk*Norm
            write(out_unitp,*) 'norm,Hpsi',kk,Norm

            If(Norm <= propa%eps) Then
                print*,'Taylor condition is fulfild after',kk,'iteration'
                exit
            else
                psi%CVec(:)    = Hpsi%CVec(:)
            Endif
              !CALL BasisTOGrid_cplx(G%CVec,psi_dt%CVec,psi_dt%Basis)
        Enddo

       CALL Calc_Norm_Grid(psi_dt%CVec, Norm,Basis)
        write(out_unitp,*) 'norm,psi_dt',Norm , 'Norm precision =',abs(Norm-ONE)
       write(out_unitp,*) 'psi_dt',psi_dt%CVec
        write(out_unitp,*) 'END march_taylor'

    END SUBROUTINE march_taylor


    SUBROUTINE marh_RK4th(psi,psi_dt,t,propa,Basis)
       USE op_m
       USE psi_m
       USE Basis_m

       TYPE (psi_t),  intent(inout)    :: psi_dt
       TYPE (psi_t),  intent(in)       :: psi
       TYPE(Basis_t),    intent(in)    :: basis
       TYPE (psi_t)       :: K1,K2,K3,K4 ,psi_inter
       TYPE(propa_t), intent(in)        :: propa

       real(kind=Rk), intent(in)    :: t
       real(kind=Rk)                ::  Norm, Norm0
       !  variables locales

       call init_psi(K1,basis,cplx = .true.,grid = .true.)
         call init_psi(K2,basis,cplx = .true.,grid = .true.)
           call init_psi(K3,basis,cplx = .true.,grid = .true.)
             call init_psi(K4,basis,cplx = .true.,grid = .true.)
               call init_psi(psi_inter,basis,cplx = .true.,grid = .true.)
       write(out_unitp,*) 'BEGINNIG march_RK4th',t,propa%delta_t
       write(out_unitp,*) 'psi',psi%CVec

       psi_dt%CVec    = psi%CVec
       CALL mEyeHPsi(psi,K1,basis)

       psi_inter%CVec = psi%CVec+(propa%delta_t/2._Rk)*K1%CVec
       CALL mEyeHPsi(psi_inter,K2,basis)
       psi_inter%CVec = psi%CVec+(propa%delta_t/2._Rk)*K2%CVec
       CALL mEyeHPsi(psi_inter,K3,basis)
       psi_inter%CVec = psi%CVec+propa%delta_t*K3%CVec
       CALL mEyeHPsi(psi_inter,K4,basis)
       psi_dt%CVec(:) = psi_dt%CVec(:)+(propa%delta_t/6._Rk)*(K1%CVec(:)+2*K2%CVec(:)+2*K3%CVec(:)+K4%CVec(:))
       CALL Calc_Norm_Grid(psi_dt%CVec, Norm,basis)
       CALL Calc_Norm_Grid(psi%CVec,Norm0,basis)
       write(out_unitp,*) 'norm,psi_dt',Norm , 'Norm precision =',ABS(ONE-Norm)
       write(out_unitp,*) 'psi_dt',psi_dt%CVec
       write(out_unitp,*) 'END marh_RK4th'

    END SUBROUTINE marh_RK4th




    SUBROUTINE read_propa( propa)
      USE psi_m
      implicit none
      TYPE (propa_t),  intent(inout) :: propa
      real(kind= Rk)                 :: t0,tf,delta_t ,eps
      character(len= 40)             :: propa_name
      !character(len=:),ALLOCATABLE   :: propa_name
      integer                        ::  max_iter

      namelist /prop/ t0,tf,delta_t ,eps,max_iter, propa_name
      t0  = ZERO
      tf  = 10._Rk
      delta_t = 0.001
      eps= ONETENTH**10._Rk
      max_iter = 5000
      propa_name = 'rk4'

      read(*,nml=prop)

      propa%t0 = t0
      propa%tf = tf
      propa%delta_t = delta_t
      propa%eps = eps
      propa%max_iter = max_iter
      propa%propa_name = propa_name

    END SUBROUTINE read_propa
    SUBROUTINE mEyeHPsi (psi_in,psi_out,Basis) !calcul de -iHpsi
       USE op_m
       USE psi_m

       TYPE (psi_t),  intent(in)   :: psi_in
       TYPE (psi_t),  intent(inout):: psi_out
       TYPE(Basis_t)  ,  intent(in)   :: Basis



       !CALL calc_OpPsi(H,psi_in,psi_out)
       CALL Calc_Hpsi(psi_in%CVec,psi_out%CVec,basis)

       psi_out%CVec(:)    = - EYE*psi_out%CVec(:)


    END SUBROUTINE mEyeHPsi

    SUBROUTINE march(psi,psi_dt,t,propa,Basis)
       USE op_m
       USE psi_m
       TYPE(Basis_t)     , INTENT(IN)   :: Basis
       TYPE (propa_t) , INTENT(IN)     :: propa
       TYPE (psi_t)   , INTENT(INOUT)  :: psi
       TYPE (psi_t)   , INTENT(INOUT)  :: psi_dt
       real(kind=Rk)  ,INTENT(IN)     :: t



       !ALLOCATE(psi_dt%CVec(Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nb  ))
       select case (propa%propa_name)
            case ('rk4')
            CALL marh_RK4th(psi,psi_dt,t,propa,Basis)
            case ('taylor')
            CALL  march_taylor(psi,psi_dt,t,propa,Basis)
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

    END SUBROUTINE write_propa


      SUBROUTINE ana_wp(rho_ana,t)
        USE NumParameters_m
          USE op_m
          USE psi_m

          TYPE (psi_t),  intent(inout) :: rho_ana
          real(kind=Rk), intent(in)    :: t

          REAL(kind=Rk)                :: alpha,hbar,v,k0,mass,w0,Norm1
          complex (kind=Rk)            :: c0,c3
          TYPE (psi_t)                 ::c1,c2,c4

          hbar  = ONE
          mass  = ONE
          alpha = TWO
          k0    = ONE
          w0    = (hbar*k0*k0)/(TWO*mass)
          v     = (hbar*k0)/mass

          CALL init_psi(rho_ana,rho_ana%Basis,cplx=.TRUE.,grid =.true.) ! to be changed
          CALL init_psi(c1,c1%Basis,cplx=.TRUE.,grid =.true.) ! to be changed
          CALL init_psi(c2,c2%Basis,cplx=.TRUE. ,grid =.true.) ! to be changed
          CALL init_psi(c4,c4%Basis,cplx=.TRUE.,grid =.true.) ! to be changed

          c0 = SQRT(PI/(alpha*alpha + EYE*(hbar*t)/(TWO*mass)))
          c1%CVec(:) = EYE*(k0*rho_ana%Basis%x(:)-w0*t)
          c2%CVec(:) = -(rho_ana%Basis%x(:)-v*t)**TWO
          c3 = FOUR*(alpha*alpha +EYE*(hbar*t)/(TWO*mass))
          c4%CVec(:) = c2%CVec(:)/c3
          rho_ana%CVec(:) = c0*EXP(c1%CVec(:))*EXP(c4%CVec(:))

          CALL Calc_Norm_Grid(rho_ana%CVec, Norm1,rho_ana%Basis)
          rho_ana%CVec(:) = rho_ana%CVec(:)/Norm1

        END SUBROUTINE ana_wp


      SUBROUTINE Calc_average_energy(psi_in,Basis,E)
      !======================================================
      !
      !     E = <Psi | H | Psi>
      !
      !======================================================
      USE UtilLib_m
      USE psi_m
      USE Basis_m

      TYPE (psi_t),  intent(in)                    :: psi_in
      REAL(KIND= Rk) ,intent(inout)                :: E
      TYPE (psi_t)                                 :: Hpsi,B,B1
      Type(Basis_t), intent(in)                    :: Basis
      REAL(KIND= Rk)                               :: dott

      CALL init_psi(Hpsi,Basis,  cplx=.TRUE.,grid =.true.)
      CALL init_psi(B ,Basis,    cplx=.TRUE.,grid   =.false.)
      CALL init_psi(B1,Basis,    cplx=.TRUE.,grid   =.false.)

      CALL Calc_Hpsi(psi_in%CVec,Hpsi%CVec,Basis)
      CALL GridTOBasis_cplx(B%CVec,Hpsi%CVec,Basis)
      CALL GridTOBasis_cplx(B1%CVec,psi_in%CVec,Basis)
           dott =  real( dot_product(B1%CVec,B1%CVec),kind=Rk)

           E = real( dot_product(B%CVec,B1%CVec),kind=Rk)
           E = E/dott

      print*,' E = <Psi | H | Psi> =',E, '<Psi | Psi> =',dott
      End SUBROUTINE Calc_average_energy

      SUBROUTINE autocor_func(psi_in,psi_out,ki,argki)
        USE UtilLib_m
        USE op_m
        USE psi_m

        TYPE (psi_t),  intent(in)                :: psi_in,psi_out
        COMPLEX(KIND=Rk)  ,  INTENT(INOUT)       :: ki
        REAL(KIND=Rk)     ,INTENT(INOUT)         ::argki

            ki= dot_product(psi_in%CVec,psi_out%CVec)
            argki =ATAN(AIMAG(ki)/REAL(ki))
      End SUBROUTINE autocor_func

      SUBROUTINE spectrum(f,time,sf,delta_t,N)
        USE NumParameters_m
        USE UtilLib_m
        COMPLEX(KIND=Rk),INTENT(IN) , DIMENSION(:)    :: f
        COMPLEX(KIND=Rk),INTENT(INOUT) , DIMENSION(0:999) :: sf
        REAL(KIND=Rk),INTENT(IN) , DIMENSION(:)       :: time
        REAL(KIND=Rk),ALLOCATABLE, DIMENSION(:)       :: w
        REAL(KIND=Rk),INTENT(IN)                      :: delta_t
        INTEGER,INTENT(IN)                            :: N
        REAL(KIND=Rk)                                 :: wm,wmax,dw
        INTEGER                                       :: Iw,nw,I
        OPEN(UNIT=20,FILE="spectrum")
        sf(:) = CZERO
        wm = -0.5;wmax = 9.5;dw = 0.01
        nw = int((wmax-wm)/dw)
        ALLOCATE(w(0:nw-1))
        DO Iw = 0,nw-1
           w(Iw) = wm + float(Iw)*dw
           Sf(Iw)= ZERO
           DO I = 1,N
               sf(Iw) = sf(Iw)+ f(I)*EXP(EYE*w(Iw)*time(I))*delta_t
           ENDDO !itime
           WRITE(20,*)  w(Iw), ABS(sf(Iw))
        ENDDO !omega
      END SUBROUTINE spectrum

end module Propa_m
