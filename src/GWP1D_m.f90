module GWP1D_m
USE NumParameters_m
USE Basis_m
USE psi_m
implicit none


TYPE  :: GWP1D_t ! gaussian 1D
  real(KIND=Rk)      ::Q0
  real(KIND=Rk)      ::K
  real(KIND=Rk)      ::phase
  real(KIND=Rk)      ::DQ
END TYPE    GWP1D_t
public::  Read_GWP1D,Write_GWP1D
contains
SUBROUTINE Read_GWP1D(paragwp1D,nio)
   USE NumParameters_m
    implicit none
    integer           , intent(in)                  :: nio
    logical           ,parameter                    :: debug = .true.
    type( GWP1D_t)    ,intent(inout)                :: paragwp1D
    real(kind= Rk)                                  :: DQ, Q0, k, phase
   integer                                          :: err_io

    namelist/defWP0/ DQ, Q0, k, phase
      DQ=0.2;
      Q0=0.0;
      k=0.0;
      phase=0.;
 read(nio,nml=defWP0,IOSTAT=err_io)

 IF (err_io /= 0) THEN
      write(out_unitp,*) ' ERROR in Read_GWP1D *** defWP0'
      stop 'problem while reading Read_GWP1D  *** defWP0'
   END IF

       paragwp1D%Q0    = Q0
       paragwp1D%K     = K
       paragwp1D%phase = phase
       paragwp1D%DQ    = DQ

END SUBROUTINE Read_GWP1D




    SUBROUTINE Write_GWP1D(paragwp1D)
   USE NumParameters_m
   Use psi_m
    implicit none
    type( GWP1D_t)    ,intent(inout)                :: paragwp1D

   write(out_unitp,*) 'Q0 =',  paragwp1D%Q0
   write(out_unitp,*) 'K= ',  paragwp1D%K
   write(out_unitp,*) 'phase=',paragwp1D%phase
   write(out_unitp,*) 'DQ=', paragwp1D%DQ

END SUBROUTINE Write_GWP1D
    SUBROUTINE GWP_1D(paragwp1D,psi0_1D,Basis,I_ElecChannel)
        implicit none
        type(GWP1D_t), intent(in)               :: paragwp1D
        type(Basis_t), intent(in)               :: Basis
        type(psi_t)   ,intent(inout)            :: psi0_1D
        complex(kind= Rk) ,allocatable          ::  psi0_1D_el(:,:)
        integer           , intent(in)          :: I_ElecChannel
        integer                                 :: i
        real(kind= Rk)                          :: dot_prdct
        ALLOCATE(psi0_1D_el(Basis%tab_basis(1)%nq, Basis%tab_basis(2)%nb))
        do i = 1,I_ElecChannel
            psi0_1D_el(:,i) = exp(-((basis%tab_basis(1)%x(:)-paragwp1D%Q0)/paragwp1D%DQ)**2  &
                         &+ EYE*paragwp1D%K*(basis%tab_basis(1)%x(:)-paragwp1D%Q0) +EYE*paragwp1D%phase )
        end do
        psi0_1D%CVec(:)  = reshape( psi0_1D_el,[Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nb])
            !CALL Calc_dot_product(psi0_1D%CVec,dot_prdct,Basis,grid=.true.,yes=.false.)
            !psi0_1D%CVec(:)  =  psi0_1D%CVec(:)/SQRT(dot_prdct)
           ! CALL Calc_dot_product(psi0_1D%CVec,dot_prdct,Basis,grid=.true.,yes=.true.)

    END SUBROUTINE GWP_1D


end module GWP1D_m