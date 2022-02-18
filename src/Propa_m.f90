module Propa_m
  USE NumParameters_m
   USE psi_m
  implicit none
   
   TYPE propa_t
   real (kind=Rk)      :: t0
   real (kind=Rk)      :: tf
   real (kind=Rk)      :: delta_t
   real (kind=Rk)      :: eps 
   integer(kind=Rk)    :: max_iter
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
  CALL read_propa( propa)
  write(out_unitp,*) 'propa_t= ',propa

  write(out_unitp,*) 'BEGINNIG propagation', propa%t0,propa%tf,propa%delta_t

  nt = int((propa%tf-propa%t0)/propa%delta_t)

  CALL init_psi(psi,psi0%Basis,cplx=.TRUE.) ! to be changed
  CALL init_psi(psi_dt,psi0%Basis,cplx=.TRUE.) ! to be changed

  psi%CVec(:) = psi0%CVec

  DO i=0,nt-1

    t = i*propa%delta_t
    t_deltat = t + propa%delta_t

    write(out_unitp,*) 'march taylor',i,t,t_deltat
    CALL march_taylor(psi,psi_dt,H,t,propa)

    psi%CVec(:) = psi_dt%CVec

  END DO

  psif%CVec(:) = psi%CVec

  CALL dealloc_psi(psi)
  CALL dealloc_psi(psi_dt)


  write(out_unitp,*) 'END propagation'

  END SUBROUTINE propagation

  SUBROUTINE march_taylor(psi,psi_dt,H,t,propa)
  USE op_m
  USE psi_m

  TYPE (psi_t),  intent(inout) :: psi_dt
  TYPE (psi_t),  intent(inout) :: psi
  TYPE(Op_t),    intent(in)    :: H
  TYPE (psi_t)                 :: Hpsi 
  TYPE(propa_t) ,intent(inout)    :: propa

  real(kind=Rk), intent(in)    :: t
  real(kind=Rk)                :: Rkk, Norm
  integer(kind = Rk)           :: kk
  ! variables locales
   CALL read_propa( propa)

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
        CALL Calc_Norm(Hpsi, Norm)
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
        write(out_unitp,*) 'norm,psi_dt',Norm , 'Norm precision =',abs(Norm-ONE)
        write(out_unitp,*) 'psi_dt',psi_dt%CVec
  write(out_unitp,*) 'END march_taylor'
  
  END SUBROUTINE march_taylor
  
  
  
  
  SUBROUTINE read_propa( propa)
  USE psi_m
implicit none
  TYPE (propa_t),  intent(inout) :: propa
  real(kind= Rk) :: t0,tf,delta_t ,eps
  integer(kind= Rk) ::  max_iter
 
  namelist /prop/ t0,tf,delta_t ,eps,max_iter
   t0  = ZERO 
   tf  = 10._Rk  
   delta_t = 1._Rk  
   eps= ONETENTH**10  
   max_iter = 5000

  read(*,nml=prop)
  
  propa%t0 = t0
  propa%tf = tf
  propa%delta_t = delta_t
  propa%eps = eps
  propa%max_iter = max_iter
     
 !write(out_unitp,*) 'propa_t= ',propa_t 
   
   END SUBROUTINE read_propa
  
  
  
  
  
  

end module Propa_m
