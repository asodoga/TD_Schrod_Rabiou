module GWPmnD_m
    USE NumParameters_m
    USE Basis_m
    USE psi_m
    USE GWP1D_m
    USE GWPnD_m
    implicit none


TYPE  ::  GWPmnD_t !  multidim gaussian mnD (nb_GWP)
     complex(kind= Rk),allocatable   :: Coef(:)
     type(GWPnD_t),allocatable       :: mnD(:) ! nb_GWP
END TYPE   GWPmnD_t
public::  Read_GWPmnD,Write_GWPmnD
contains


SUBROUTINE Read_GWPmnD(paragwp,nio,ndim,nb_GWP)
   USE NumParameters_m
    implicit none
    integer           , intent(in)                  :: nio,nb_GWP,ndim
    logical           ,parameter                    :: debug = .true.
    type( GWPmnD_t)    ,intent(inout)               :: paragwp
   integer                                          :: err_io,i
   complex(kind= Rk)                                :: Coef
   integer                                          :: I_ElecChannel


   namelist/defGWP /I_ElecChannel ,Coef
    I_ElecChannel=2 ;
    Coef=CONE;
        allocate(paragwp%mnD(nb_GWP ))
        allocate(paragwp%Coef(nb_GWP ))
   do i = 1,nb_GWP,1
        read(nio,nml=defGWP ,IOSTAT=err_io)
        paragwp%Coef(i) = Coef
       CALL Read_GWPnD(paragwp%mnD(i),nio,ndim)
   end do

END SUBROUTINE Read_GWPmnD

    SUBROUTINE Write_GWPmnD(paragwp,nb_GWP)
   USE NumParameters_m
    implicit none
    type( GWPmnD_t)    ,intent(inout)                  :: paragwp
   integer                                             :: nb_GWP
     integer                                           :: i
  do i = 1, nb_GWP,1
       CALL Write_GWPnD(paragwp%mnD(i))
        write(*,*)'Coef = ',paragwp%Coef(i)
  end do

   END SUBROUTINE Write_GWPmnD


    SUBROUTINE GWP_mnD(paragwp,psi0_mnD,Basis,I_ElecChannel,nb_GWP)
        implicit none
        type(GWPmnD_t), intent(in)               :: paragwp
        type(Basis_t), intent(inout)             :: Basis
        type(psi_t)   ,intent(inout)             :: psi0_mnD
        type(psi_t),allocatable                  :: GWmnD(:)
        integer           , intent(in)           :: I_ElecChannel,nb_GWP
        integer                                  :: i

        allocate(GWmnD(nb_GWP))
        psi0_mnD%CVec(:)=CZERO
        do i = 1, nb_GWP,1
             CALL init_psi(GWmnD(i)  ,Basis,    cplx=.TRUE.   ,grid =.true.)
            CALL GWP_nD(paragwp%mnD(i), GWmnD(i),Basis,I_ElecChannel)
            psi0_mnD%CVec(:) = psi0_mnD%CVec(:)+paragwp%Coef(i)*GWmnD(i)%CVec(:)
        end do

    END SUBROUTINE GWP_mnD
    

    end module GWPmnD_m

