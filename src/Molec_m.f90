module Molec_m
  USE QDUtil_m
  implicit none
  private

  real(kind=Rkind) :: mass =  2000._Rkind
  integer :: potential_type = 1 !0 internal,1 QLM,2 options
  public :: Calc_pot,mass,sub_pot

contains
  FUNCTION Calc_pot(Q)
    real(kind=Rkind)             :: Calc_pot

    real(kind=Rkind), intent(in) :: Q

  Calc_pot = mass*HALF * Q*Q

  END FUNCTION Calc_pot

  SUBROUTINE  sub_pot(Mat_V,Q,potential_type)
       USE QDUtil_m
       REAL(kind=Rkind), intent(inout)   :: Mat_V(:,:)
       REAL(kind=Rkind), intent(in)      :: Q(:)
       INTEGER                        :: i,j,iq
       integer , intent(in)           :: potential_type


        !IF (size(Q) /= 1) STOP 'wrong dimension'
       ! IF (size(Mat_V,dim=1) /= 2) STOP 'wrong number of electronic state'

       SELECT CASE (potential_type)
   CASE (0)

                Mat_V(:,:) = 0
            do i = 1, size(Mat_V(1,:))
               do iq = 1, size(Q)
                Mat_V(i,i) = Mat_V(i,i)+ HALF*(Q(iq))**2
               end do
                !do j = 1, size(Mat_V(:,1))
                !    if (abs(i-j)== 1) then
                !         Mat_V(i,j) = ZERO !0.111803_Rkind*Q(1)
                !    end if
                ! end do
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
