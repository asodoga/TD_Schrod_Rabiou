module Molec_m
  USE NumParameters_m
  implicit none
  private

  real(kind=Rk) :: mass = ONE

  public :: Calc_pot,mass,sub_pot

contains
  FUNCTION Calc_pot(Q)
    real(kind=Rk)             :: Calc_pot

    real(kind=Rk), intent(in) :: Q

  Calc_pot = mass*HALF * Q*Q

  END FUNCTION Calc_pot


  SUBROUTINE  sub_pot(Mat_V,Q)
       USE NumParameters_m
       REAL(kind=Rk), intent(inout)   :: Mat_V(:,:)
       REAL(kind=Rk), intent(in)      :: Q(:)
       INTEGER                        :: i,j
       !REAL(kind=Rk), intent(in)      :: Q

        !IF (size(Q) /= 1) STOP 'wrong dimension'
        !IF (size(Mat_V,dim=1) /= 2) STOP 'wrong number of electronic state'
                Mat_V(:,:) = 0
            do i = 1, size(Mat_V(1,:))

               Mat_V(i,i) = (Q(1)+i)**2
                do j = 1, size(Mat_V(:,1))
                    if (abs(i-j)== 1) then

                         Mat_V(i,j) = 0.001*Q(1)
                    end if
                 end do
            end do

               !Mat_V(1,1) = (Q(1)-1)**2

               !Mat_V(2,2) = (Q(1)+1)**2

              ! Mat_V(1,2) = 0.001*Q(1)
               !Mat_V(2,1) = 0.001*Q(1)


  END SUBROUTINE sub_pot



end module Molec_m
