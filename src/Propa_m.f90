module Propa_m
  USE NumParameters_m
   USE psi_m
  implicit none


contains
  SUBROUTINE propagation(psif,psi0,H,t0,tf,delta_t)
  USE op_m
  USE psi_m

  TYPE (psi_t),  intent(inout) :: psif
  TYPE (psi_t),  intent(in)    :: psi0
  TYPE(Op_t),    intent(in)    :: H

  real(kind=Rk), intent(in)    :: t0,tf,delta_t

  ! variables locales
  real(kind=Rk) :: t,t_deltat
  integer       :: i,nt
  TYPE (psi_t)  :: psi,psi_dt

  write(out_unitp,*) 'BEGINNIG propagation', t0,tf,delta_t

  nt = int((tf-t0)/delta_t)

  CALL init_psi(psi,psi0%Basis,cplx=.TRUE.) ! to be changed
  CALL init_psi(psi_dt,psi0%Basis,cplx=.TRUE.) ! to be changed

  psi%CVec(:) = psi0%CVec

  DO i=0,nt-1

    t = i*delta_t
    t_deltat = t + delta_t

    write(out_unitp,*) 'march taylor',i,t,t_deltat
    CALL march_taylor(psi,psi_dt,H,t,delta_t)

    psi%CVec(:) = psi_dt%CVec

  END DO

  psif%CVec(:) = psi%CVec

  CALL dealloc_psi(psi)
  CALL dealloc_psi(psi_dt)


  write(out_unitp,*) 'END propagation'

  END SUBROUTINE propagation

  SUBROUTINE march_taylor(psi,psi_dt,H,t,delta_t)
  USE op_m
  USE psi_m

  TYPE (psi_t),  intent(inout) :: psi_dt
  TYPE (psi_t),  intent(inout)    :: psi
  TYPE(Op_t),    intent(in)    :: H
  TYPE (psi_t)                 :: Hpsi 

  real(kind=Rk), intent(in)    :: t,delta_t
  real(kind=Rk)                :: Rkk, Norm,eps
  integer(kind = Rk)           :: kk
  ! variables locales


  write(out_unitp,*) 'BEGINNIG march_taylor',t,delta_t
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
   Do kk = 1,5000,1
   CALL calc_OpPsi(H,psi,Hpsi)
    Rkk = Rkk*(delta_t/kk)
    Hpsi%CVec(:)    = - EYE*Hpsi%CVec(:)
    psi_dt%CVec(:) = psi_dt%CVec(:) +Rkk*Hpsi%CVec(:)
      
       !CALL Calc_Norm(psi_dt, Norm)
        !write(out_unitp,*) 'norm,psi_dt',Norm
        CALL Calc_Norm(Hpsi, Norm)
         Norm =   Rkk*Norm
  write(out_unitp,*) 'norm,Hpsi',Norm
 
  If(Norm .le. eps) Then
   print*,'Taylor condition is fulfild after',kk,'iteration'
   exit
   else 
   psi%CVec(:)    = Hpsi%CVec(:)
   
  
  
  Endif
  
  Enddo
  
  CALL Calc_Norm(psi_dt, Norm)
        write(out_unitp,*) 'norm,psi_dt',Norm , 'Norm precision =',Norm-ONE
        write(out_unitp,*) 'psi_dt',psi_dt%CVec
  write(out_unitp,*) 'END march_taylor'
  
  END SUBROUTINE march_taylor
  
  
  
  
  
  
  
  

end module Propa_m
