module Molec_m
  USE NumParameters_m
  implicit none
  private

  real(kind=Rk) :: mass = ONE
  integer :: potential_type = 1 !0 internal,1 QLM,2 options

  public :: Calc_pot,mass,sub_pot

contains
  FUNCTION Calc_pot(Q)
    real(kind=Rk)             :: Calc_pot

    real(kind=Rk), intent(in) :: Q

  Calc_pot = mass*HALF * Q*Q

  END FUNCTION Calc_pot

  SUBROUTINE  sub_pot(Mat_V,Q,potential_type)
       USE NumParameters_m
       REAL(kind=Rk), intent(inout)   :: Mat_V(:,:)
       REAL(kind=Rk), intent(in)      :: Q(:)
       INTEGER                        :: i,j
       integer , intent(in)           :: potential_type


        !IF (size(Q) /= 1) STOP 'wrong dimension'
       ! IF (size(Mat_V,dim=1) /= 2) STOP 'wrong number of electronic state'

       SELECT CASE (potential_type)
   CASE (0)

                Mat_V(:,:) = 0
            do i = 1, size(Mat_V(1,:))

               Mat_V(i,i) = HALF*dot_product(Q,Q)
                do j = 1, size(Mat_V(:,1))
                    if (abs(i-j)== 1) then

                         Mat_V(i,j) = ZERO
                         !0.001_Rk*Q(1)
                    end if
                 end do
            end do

   CASE (1) !QML
     CALL sub_Qmodel_V(Mat_V,Q)
   CASE (2)!QML
      stop "ERROR potential_type=2 not define"
      CASE DEFAULT
      stop "no default in sub_pot"
END SELECT


  END SUBROUTINE sub_pot



end module Molec_m
