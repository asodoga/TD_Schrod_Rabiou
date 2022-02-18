module Propa_m
  USE NumParameters_m
   USE psi_m
  implicit none
   
   TYPE propa_t
   real (kind=Rk)      :: t0
   real (kind=Rk)      :: tf
   real (kind=Rk)      :: delta_t
   real (kind=Rk)      :: eps 
   integer (kind=Rk)   :: max_iter
  END TYPE propa_t

  
contains
  SUBROUTINE propagation(psif,psi0,H,propa)
  USE op_m
  USE psi_m

  TYPE (psi_t),  intent(inout) :: psif
  TYPE (psi_t),  intent(in)    :: psi0
  TYPE(Op_t),    intent(in)    :: H

 TYPE(propa_t), intent(inout)    :: propa

  ! variables locales
  real(kind=Rk) :: t ,t_deltat
  integer       :: i,nt
  TYPE (psi_t)  :: psi,psi_dt
  CALL read_propa(propa)
  write(out_unitp,*) 'propa= ',propa
  

  write(out_unitp,*) 'BEGINNIG propagation', propa%t0,propa%tf,propa%delta_t

  nt = int((propa%tf-propa%t0)/propa%delta_t)

  CALL init_psi(psi,psi0%Basis,cplx=.TRUE.) ! to be changed
  CALL init_psi(psi_dt,psi0%Basis,cplx=.TRUE.) ! to be changed

  psi%CVec(:) = psi0%CVec

  DO i=0,nt-1

    t = i*propa%delta_t
    t_deltat = t + propa%delta_t

    write(out_unitp,*) 'marh_RK4th',i,t,t_deltat
    CALL marh_RK4th(psi,psi_dt,H,t,propa)

    psi%CVec(:) = psi_dt%CVec

  END DO

  psif%CVec(:) = psi%CVec

  CALL dealloc_psi(psi)
  CALL dealloc_psi(psi_dt)


  write(out_unitp,*) 'END propagation'

  END SUBROUTINE propagation

  SUBROUTINE marh_RK4th(psi,psi_dt,H,t,propa)
  USE op_m
  USE psi_m

  TYPE (psi_t),  intent(inout) :: psi_dt
  TYPE (psi_t),  intent(inout) :: psi
  TYPE(Op_t),    intent(in)    :: H
  TYPE (psi_t)                 :: K1,K2,K3,K4 ,psi_inter
  TYPE(propa_t), intent(inout) :: propa

  real(kind=Rk), intent(in)    :: t
  real(kind=Rk)                ::  Norm
  ! variables locales
   CALL read_propa(propa)

  write(out_unitp,*) 'BEGINNIG march_RK4th',t,propa%delta_t
  write(out_unitp,*) 'psi',psi%CVec
  
    psi_dt%CVec    = psi%CVec
    CALL calc_OpPsi(H,psi,K1)
    K1%CVec(:)     = - EYE*K1%CVec(:)
    psi_inter%CVec = psi%CVec+(propa%delta_t/2._Rk)*K1%CVec
    CALL calc_OpPsi(H,psi_inter,K2)
    psi_inter%CVec = 0._Rk 
    K2%CVec(:)     = - EYE*K2%CVec(:)
    psi_inter%CVec = psi%CVec+(propa%delta_t/2._Rk)*K2%CVec
    CALL calc_OpPsi(H,psi_inter,K3)
    psi_inter%CVec = 0._Rk
    K3%CVec(:)     = - EYE*K3%CVec(:)
    psi_inter%CVec = psi%CVec+propa%delta_t*K3%CVec
    CALL calc_OpPsi(H,psi_inter,K4)
    K4%CVec(:)     = - EYE*K4%CVec(:)
    psi_dt%CVec(:) = psi_dt%CVec(:)+(propa%delta_t/6._Rk)*(K1%CVec(:)+2*K2%CVec(:)+2*K3%CVec(:)+K4%CVec(:))
    CALL Calc_Norm(psi_dt, Norm)
    write(out_unitp,*) 'norm,psi_dt',Norm , 'Norm precision =',ABS(Norm-ONE)
    write(out_unitp,*) 'psi_dt',psi_dt%CVec
    write(out_unitp,*) 'END marh_RK4th'
  
  END SUBROUTINE marh_RK4th
  
  
  !=====================march rk4th============================================
  
  
  
  
   SUBROUTINE march_taylor(psi,psi_dt,H,t,propa)
  USE op_m
  USE psi_m

  TYPE (psi_t),  intent(inout) :: psi_dt
  TYPE (psi_t),  intent(inout) :: psi
  TYPE(Op_t),    intent(in)    :: H
  TYPE (psi_t)                 :: Hpsi 
  TYPE(propa_t), intent(inout) :: propa

  real(kind=Rk), intent(in)    :: t
  real(kind=Rk)                :: Rkk, Norm
  integer(kind = Rk)           :: kk
  ! variables locales
   CALL read_propa(propa)

  write(out_unitp,*) 'BEGINNIG march_taylor',t,propa%delta_t
  write(out_unitp,*) 'psi',psi%CVec
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
  !=====================fin ordre deux================================
    
    
    
    
  !!===========================Ordre deux etplus=======================
   
   psi_dt%CVec    = psi%CVec
   Do kk = 1,propa%max_iter,1
   CALL calc_OpPsi(H,psi,Hpsi)
    Rkk = Rkk*(propa%delta_t/kk)
    Hpsi%CVec(:)    = - EYE*Hpsi%CVec(:)
    psi_dt%CVec(:) = psi_dt%CVec(:) +Rkk*Hpsi%CVec(:)
      
       !CALL Calc_Norm(psi_dt, Norm)
        !write(out_unitp,*) 'norm,psi_dt',Norm
        CALL Calc_Norm(Hpsi,Norm)
         Norm =   Rkk*Norm
  write(out_unitp,*) 'norm,Hpsi',Norm
 
  If(Norm .le. propa%eps) Then
   print*,'Taylor condition is fulfild after',kk,'iteration'
   exit
   else 
   psi%CVec(:)    = Hpsi%CVec(:)
   
  
  
  Endif
  
  Enddo
  
  CALL Calc_Norm(psi_dt, Norm)
        write(out_unitp,*) 'norm,psi_dt',Norm , 'Norm precision =',ABS(Norm-ONE)
        write(out_unitp,*) 'psi_dt',psi_dt%CVec
  write(out_unitp,*) 'END march_taylor'
  
  END SUBROUTINE march_taylor
  
  
  
  
  
  !=============================subroutine read_propa================================
  SUBROUTINE read_propa(propa)
  USE psi_m
  USE NumParameters_m
  implicit none
  TYPE (propa_t),  intent(inout) :: propa
  real(kind= Rk)                 :: t0, tf, delta_t, eps
  integer(kind= Rk)              :: max_iter
  
   
   
   NAMELIST / prop / t0, tf, delta_t, eps, max_iter
  
   t0       = ZERO
   tf       = 10._Rk  
   delta_t  = (ONETENTH)**3 
   eps      = ONETENTH**10._Rk 
   max_iter = 5000
   
  !read(in_unitp,nml=prop)
  propa%t0        = t0
  propa%tf        = tf
  propa%delta_t   = delta_t
  propa%eps       = eps
  propa%max_iter  = max_iter
     
 !write(out_unitp,*) 'propa= ',propa
   
   END SUBROUTINE read_propa
  
  
  
  
  
  

end module Propa_m
