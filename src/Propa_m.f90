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


contains
    SUBROUTINE propagation(psif,psi0,H,propa)
       USE op_m
       USE psi_m
       USE Basis_m

       TYPE (psi_t),  intent(inout) :: psif
       TYPE (psi_t),  intent(in)    :: psi0
       TYPE(Op_t),    intent(inout) :: H
       TYPE (psi_t)                 :: G,rho_exact,rho_calc   !,rho



       TYPE(propa_t), intent(inout) :: propa
       logical, parameter           :: debug = .true.

       ! variables locales
       REAL(kind=Rk)                :: t ,t_deltat, Norm,Norm1
       REAL(kind=Rk)                :: alpha,hbar,v,k0,mass,w0
       complex (kind=Rk)            :: c0,c3
       TYPE (psi_t)                 ::c1,c2,c4
       INTEGER                      :: i,nt,IQ,nf
       TYPE (psi_t)                 :: psi,psi_dt
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
       CALL init_psi(rho_exact,psi0%Basis,cplx=.TRUE.) ! to be changed
       CALL init_psi(rho_calc,psi0%Basis,cplx=.TRUE.) ! to be changed
       CALL init_psi(c1,psi0%Basis,cplx=.TRUE.) ! to be changed
       CALL init_psi(c2,psi0%Basis,cplx=.TRUE.) ! to be changed
       CALL init_psi(c4,psi0%Basis,cplx=.TRUE.) ! to be changed
       !CALL init_psi(rho,psi0%Basis,cplx=.TRUE.) ! to be changed

              psi%CVec(:) = psi0%CVec

        hbar = ONE
        mass = ONE
        alpha = TWO
        k0 = ONE
        w0 = (hbar*k0*k0)/(TWO*mass)
         v= (hbar*k0)/mass

       DO i=0,nt-1

            t = i*propa%delta_t
            t_deltat = t + propa%delta_t

             c0 = SQRT(PI/(alpha*alpha + EYE*(hbar*t)/(TWO*mass)))
             c1%CVec(:) = EYE*(k0*rho_exact%Basis%x(:)-w0*t)
             c2%CVec(:) = -(rho_exact%Basis%x(:)-v*t)**TWO
             c3 = FOUR*(alpha*alpha +EYE*(hbar*t)/(TWO*mass))
             c4%CVec(:) = c2%CVec(:)/c3
            rho_exact%CVec(:) = c0*EXP(c1%CVec(:))*EXP(c4%CVec(:))

            CALL Calc_Norm_Grid(rho_exact, Norm1)
            rho_exact%CVec(:) = rho_exact%CVec(:)/Norm1
            !  CALL Calc_Norm_Grid(rho_exact, Norm1)



            write(out_unitp,*) propa%propa_name,i,t,t_deltat

            CALL march(psi,psi_dt,H,t,propa)

            psi%CVec(:) = psi_dt%CVec
            rho_calc%CVec(:)= psi_dt%CVec
            !rho%CVec(:)  = rho_exact%CVec(:)-rho_calc%CVec(:)
               nf = int(nt/5)
             IF( nf == 0)then
             nf = 1

             ENDIF

            IF(   MOD(i,nf) == 0  )Then
             OPEN(unit=i+10)
             CALL BasisTOGrid_Basis_cplx(G%CVec,psi_dt%CVec,psi_dt%Basis)
             DO  IQ = 1, psi_dt%Basis%nq
             WRITE(i+10,*) G%Basis%x(IQ), ABS(G%CVec(IQ))**2, ABS(rho_exact%CVec(IQ))**2-ABS(G%CVec(IQ))**2,&
             ABS(rho_exact%CVec(IQ))**2
             ENDDO
              CLOSE(UNIT=i+10)
            ENDIF


       END DO
       psif%CVec(:) = psi%CVec
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
            psi_dt%CVec(:) = psi_dt%CVec(:) +Rkk*Hpsi%CVec(:)

            !CALL Calc_Norm(psi_dt, Norm)
            !write(out_unitp,*) 'norm,psi_dt',Norm
            CALL Calc_Norm(Hpsi, Norm)
            Norm =   Rkk*Norm
            write(out_unitp,*) 'norm,Hpsi',kk,Norm

            If(Norm <= propa%eps) Then
                print*,'Taylor condition is fulfild after',kk,'iteration'
                exit
            else
                psi%CVec(:)    = Hpsi%CVec(:)
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
       psi_dt%CVec(:) = psi_dt%CVec(:)+(propa%delta_t/6._Rk)*(K1%CVec(:)+2*K2%CVec(:)+2*K3%CVec(:)+K4%CVec(:))
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

       psi_out%CVec(:)    = - EYE*psi_out%CVec(:)


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




    !SUBROUTINE out_propa(psi0,psi_dt,propa)
        !USE op_m
        !USE psi_m

        !TYPE (psi_t),  intent(inout) :: psi_dt
        !TYPE (psi_t),  intent(inout) :: psi0
        !TYPE (psi_t)                 :: cor_fonct
        !TYPE(propa_t), intent(in)    :: propa
        !real(kind= Rk)               :: energy
       ! TYPE (psi_t)                 :: Hpsi
        !real(kind=Rk), intent(in)    :: t
        !real(kind=Rk)                ::  Norm


       ! CALL  energy(H,psi)

       !END SUBROUTINE out_propa














end module Propa_m
