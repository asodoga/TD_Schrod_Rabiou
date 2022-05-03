module Propa_m
     USE NumParameters_m
     USE psi_m

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
    public ::mEyeHPsi,write_propa,initial_wp,ana_wp,spectrum
    public :: ENERGY
contains
    SUBROUTINE propagation(psif,psi0,H,propa)
       USE op_m
       USE psi_m
       USE Basis_m

       TYPE (psi_t),  intent(inout) :: psif
       TYPE (psi_t),  intent(in)    :: psi0
       TYPE(Op_t),    intent(inout) :: H
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

       CALL init_psi(psi,psi0%Basis,cplx=.TRUE.) ! to be changed
       CALL init_psi(psi_dt,psi0%Basis,cplx=.TRUE.) ! to be changed
       CALL init_psi(G,psi0%Basis,cplx=.TRUE.) ! to be changed
       CALL init_psi(rho_num,psi0%Basis,cplx=.TRUE.) ! to be changed
       CALL init_psi(rho,psi0%Basis,cplx=.TRUE.) ! to be changed
       CALL init_psi(rho_ana,psi0%Basis,cplx=.TRUE.) ! to be changed


        psi%CVec = psi0%CVec

       DO i=0,nt-1

            t = i*propa%delta_t
            t_deltat = t + propa%delta_t
      write(out_unitp,*) propa%propa_name,i,t,t_deltat

            CALL march(psi,psi_dt,H,t,propa)

            psi%CVec(:,1) = psi_dt%CVec(:,1)

               nf = int(nt/5)
             IF( nf == 0)then
             nf = 1

             ENDIF

            IF(   MOD(i,nf) == 0  )Then
             OPEN(unit=i+10)
             CALL BasisTOGrid_Basis_cplx(G%CVec,psi_dt%CVec,psi_dt%Basis)
             rho_num%CVec(:,1)= G%CVec(:,1)
             CALL ana_wp(rho_ana,t)
             !CALL calc_rho(rho_ana,rho_num,rho)
             DO  IQ = 1, psi_dt%Basis%nq
             WRITE(i+10,*) G%Basis%x(IQ), ABS(rho_num%CVec(IQ,1))**2, ABS(rho%CVec(IQ,1)),ABS(rho_ana%CVec(IQ,1))**2
             ENDDO
              CLOSE(UNIT=i+10)
            ENDIF


       END DO
       psif%CVec(:,1) = psi%CVec(:,1)
       CALL Calc_Norm(psi_dt, Norm)

       CALL dealloc_psi(psi)
       CALL dealloc_psi(psi_dt)
       IF (debug) THEN
           write(out_unitp,*) 'END propagation'
           write(out_unitp,*) 'norm,psi_dt',Norm
           write(out_unitp,*) 'psi_dt',psi_dt%CVec
           flush(out_unitp)
       END IF




    END SUBROUTINE propagation

    SUBROUTINE march_taylor(psi,psi_dt,H,t,propa)
       USE op_m
       USE psi_m
       USE Basis_m

       TYPE (psi_t)  , INTENT(INOUT):: psi_dt
       TYPE (psi_t)  ,INTENT(INOUT) :: psi
       TYPE(Op_t)    ,INTENT(IN)    :: H
       TYPE (psi_t)                 :: Hpsi,G
       TYPE(propa_t) ,INTENT(IN)    :: propa
       real(kind=Rk) , INTENT(IN)   :: t

       ! variables locales


       real(kind=Rk)                :: Rkk, Norm
       integer                      :: kk



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
            CALL mEyeHPsi(H,psi,Hpsi)
            Rkk = Rkk*(propa%delta_t/kk)
            psi_dt%CVec(:,1) = psi_dt%CVec(:,1) +Rkk*Hpsi%CVec(:,1)

            !CALL Calc_Norm(psi_dt, Norm)
            !write(out_unitp,*) 'norm,psi_dt',Norm
            CALL Calc_Norm(Hpsi, Norm)
            Norm =   Rkk*Norm
            write(out_unitp,*) 'norm,Hpsi',kk,Norm

            If(Norm <= propa%eps) Then
                print*,'Taylor condition is fulfild after',kk,'iteration'
                exit
            else
                psi%CVec(:,1)    = Hpsi%CVec(:,1)
            Endif
              CALL BasisTOGrid_Basis_cplx(G%CVec,psi_dt%CVec,psi_dt%Basis)
        Enddo

       CALL Calc_Norm(psi_dt, Norm)
        write(out_unitp,*) 'norm,psi_dt',Norm , 'Norm precision =',abs(Norm-ONE)
       write(out_unitp,*) 'psi_dt',psi_dt%CVec
        write(out_unitp,*) 'END march_taylor'

    END SUBROUTINE march_taylor


    SUBROUTINE marh_RK4th(psi,psi_dt,H,t,propa)
       USE op_m
       USE psi_m
       USE Basis_m

       TYPE (psi_t),  intent(inout) :: psi_dt
       TYPE (psi_t),  intent(in)    :: psi
       TYPE(Op_t),    intent(in)    :: H
       TYPE (psi_t)                 :: K1,K2,K3,K4 ,psi_inter
       TYPE(propa_t), intent(in)    :: propa

       real(kind=Rk), intent(in)    :: t
       real(kind=Rk)                ::  Norm, Norm0
       !  variables locales

       write(out_unitp,*) 'BEGINNIG march_RK4th',t,propa%delta_t
       write(out_unitp,*) 'psi',psi%CVec

       psi_dt%CVec    = psi%CVec
       CALL mEyeHPsi(H,psi,K1)

       psi_inter%CVec = psi%CVec+(propa%delta_t/2._Rk)*K1%CVec
       CALL mEyeHPsi(H,psi_inter,K2)
       psi_inter%CVec = psi%CVec+(propa%delta_t/2._Rk)*K2%CVec
       CALL mEyeHPsi(H,psi_inter,K3)
       psi_inter%CVec = psi%CVec+propa%delta_t*K3%CVec
       CALL mEyeHPsi(H,psi_inter,K4)
       psi_dt%CVec(:,1) = psi_dt%CVec(:,1)+(propa%delta_t/6._Rk)*(K1%CVec(:,1)+2*K2%CVec(:,1)+2*K3%CVec(:,1)+K4%CVec(:,1))
       CALL Calc_Norm(psi_dt, Norm)
       CALL Calc_Norm(psi,Norm0)
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
    SUBROUTINE mEyeHPsi (H,psi_in,psi_out) !calcul de -iHpsi
       USE op_m
       USE psi_m

       TYPE (psi_t),  intent(in)   :: psi_in
       TYPE (psi_t),  intent(inout):: psi_out
       TYPE(Op_t)  ,  intent(in)   :: H



       CALL calc_OpPsi(H,psi_in,psi_out)

       psi_out%CVec(:,1)    = - EYE*psi_out%CVec(:,1)


    END SUBROUTINE mEyeHPsi

    SUBROUTINE march(psi,psi_dt,H,t,propa)
       USE op_m
       USE psi_m
       TYPE(Op_t)     , INTENT(IN)     :: H
       TYPE (propa_t) , INTENT(IN)     :: propa
       TYPE (psi_t)   , INTENT(INOUT)  :: psi
       TYPE (psi_t)   , INTENT(INOUT)  :: psi_dt
       real(kind=Rk)  ,INTENT(IN)     :: t




       select case (propa%propa_name)
            case ('rk4')
            CALL marh_RK4th(psi,psi_dt,H,t,propa)
            case ('taylor')
            CALL  march_taylor(psi,psi_dt,H,t,propa)
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

          CALL init_psi(rho_ana,rho_ana%Basis,cplx=.TRUE.) ! to be changed
          CALL init_psi(c1,c1%Basis,cplx=.TRUE.) ! to be changed
          CALL init_psi(c2,c2%Basis,cplx=.TRUE.) ! to be changed
          CALL init_psi(c4,c4%Basis,cplx=.TRUE.) ! to be changed

          c0 = SQRT(PI/(alpha*alpha + EYE*(hbar*t)/(TWO*mass)))
          c1%CVec(:,1) = EYE*(k0*rho_ana%Basis%x(:)-w0*t)
          c2%CVec(:,1) = -(rho_ana%Basis%x(:)-v*t)**TWO
          c3 = FOUR*(alpha*alpha +EYE*(hbar*t)/(TWO*mass))
          c4%CVec(:,1) = c2%CVec(:,1)/c3
          rho_ana%CVec(:,1) = c0*EXP(c1%CVec(:,1))*EXP(c4%CVec(:,1))

          CALL Calc_Norm_Grid(rho_ana, Norm1)
          rho_ana%CVec(:,1) = rho_ana%CVec(:,1)/Norm1

        END SUBROUTINE ana_wp







    SUBROUTINE initial_wp(B,psi0,G)
      USE NumParameters_m
      USE Basis_m
      USE psi_m

      TYPE(psi_t),INTENT(INOUT)     :: B,G
      TYPE(psi_t),INTENT(IN)        :: psi0
      INTEGER                       :: IQ , IB

       COMPLEX(KIND=Rk)             :: alpha0,gamma0
       REAL(KIND= Rk)               :: k,mass, omega,Norm,Norm1
       REAL(kind=Rk)                 ::alpha,k0,phase,Q0,sig0,sigma

      OPEN(unit=11,file = 'norm' )
      OPEN(unit=12,file = 'G' )
      OPEN(unit=13,file = 'B' )

        CALL init_psi(G,psi0%Basis,cplx=.TRUE.)
        CALL init_psi(B,psi0%Basis,cplx=.TRUE.)
        !sigma = HALF
      !  sig0 = TWO
        k0 = ONE
      !  phase = ZERO
        Q0 = ZERO
        alpha= TWO
       !omega = SQRT(k/mass)
       !alpha0 = HALF*EYE*mass*SQRT(k*mass)


      !  gamma0     = -EYE*LOG((SQRT(mass*omega)/PI)**0.25)
        !G%CVec(:)  = EXP(EYE*alpha0*(G%Basis%x(:)-ONE)**2+EYE*gamma0)
        !G%CVec(:)  =SQRT(PI/alpha**2)*EXP(EYE*k0*(G%Basis%x(:)-Q0))*EXP(-(G%Basis%x(:)-Q0)**2/(FOUR*alpha**2))
        !G%CVec(:)  = EXP(-((Basis%x(:)-Q0)/(2d0*sig0))**2)* EXP(EYE*k0*Basis%x(:))
        !G%CVec(:)  = EXP(-(ONETENTH**3)*((Basis%x(:)-Q0)/sigma)**2)*EXP(EYE*k0*(Basis%x(:)-Q0)+ EYE*phase)
        G%CVec(:,1) = CONE


         CALL Calc_Norm_Grid(G, Norm1)
         G%CVec(:,1) = G%CVec(:,1)/Norm1
        CALL Calc_Norm_Grid(G, Norm1)
         DO  IQ = 1, G%Basis%nq
           write(12,*) G%Basis%x(IQ), ABS(G%CVec(IQ,1))**2
         ENDDO
       CALL GridTOBasis_Basis_cplx(B%CVec,G%CVec,G%Basis)
       CALL Calc_Norm(B, Norm)
       !B%CVec(:) = B%CVec(:)/Norm
       write(11,*) Norm1,Norm
       DO  IB = 1, G%Basis%nb
           write(13,*) G%Basis%x(IB), ABS(B%CVec(IB,1))**2
         ENDDO
      END SUBROUTINE initial_wp


      SUBROUTINE ENERGY(psi,H,E)
        USE UtilLib_m
        USE op_m
        USE psi_m

        TYPE (psi_t),  intent(in)       :: psi
        TYPE(Op_t)  ,  intent(in)       :: H
        REAL(KIND= Rk)                  :: E
        TYPE (psi_t)                    :: Hpsi

            CALL calc_OpPsi(H,psi,Hpsi)
            E= REAL(dot_product(psi%CVec(:,1),Hpsi%CVec(:,1)),KIND= Rk)

      End SUBROUTINE ENERGY
      SUBROUTINE autocor_func(psi_in,psi_out,ki,argki)
        USE UtilLib_m
        USE op_m
        USE psi_m

        TYPE (psi_t),  intent(in)                :: psi_in,psi_out
        COMPLEX(KIND=Rk)  ,  INTENT(INOUT)       :: ki
        REAL(KIND=Rk)     ,INTENT(INOUT)         ::argki

            ki= dot_product(psi_in%CVec(:,1),psi_out%CVec(:,1))
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
