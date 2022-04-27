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
  ! for dibatic one need multidimenssional pot
  SUBROUTINE  sub_pot(Mat_V,Q)
    REAL(kind=Rk), intent(inout)   :: Mat_V(:,:,:)
    REAL(kind=Rk), intent(in)      :: Q(:)
    INTEGER                        :: iq

    Do iq = 1, SIZE(Q)
    Mat_V(Q(iq),1,1) = (Q(iq)+1)**2
    Mat_V(Q(iq),2,2) = (Q(iq)-1)**2
    Mat_V(Q(iq),1,2) = 0.01*Q(iq)
    Mat_V(Q(iq),2,1) = 0.01*Q(iq)
   END DO


  END SUBROUTINE sub_pot

end module Molec_m
