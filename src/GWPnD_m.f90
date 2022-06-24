module GWPnD_m
    USE NumParameters_m
    USE Basis_m
    USE psi_m
    USE GWP1D_m
    implicit none


TYPE  ::  GWPnD_t ! gaussian nD (ndim)
     integer                         :: ndim
     type(GWP1D_t),allocatable       :: nD(:) ! ndim
END TYPE   GWPnD_t
public::  Read_GWPnD,Write_GWPnD
contains
SUBROUTINE Read_GWPnD(paragwp,nio,ndim)
   USE NumParameters_m
    implicit none
    integer           , intent(in)                  :: nio,ndim
    logical           ,parameter                    :: debug = .true.
    type( GWPnD_t)    ,intent(inout)                :: paragwp
   integer                                          :: err_io,i

       paragwp%ndim = ndim
       allocate(paragwp%nD( paragwp%ndim))
   do i = 1, paragwp%ndim,1
       CALL Read_GWP1D(paragwp%nD(i),nio)
   end do

END SUBROUTINE Read_GWPnD

    SUBROUTINE Write_GWPnD(paragwp)
   USE NumParameters_m
    implicit none
    type( GWPnD_t)    ,intent(inout)                  :: paragwp
     integer                                          :: i
  do i = 1, paragwp%ndim,1
       CALL Write_GWP1D(paragwp%nD(i))
  end do
         write(out_unitp,*) 'ndim =',  paragwp%ndim

   END SUBROUTINE Write_GWPnD


    SUBROUTINE GWP_nD(paragwp,psi0_nD,Basis,I_ElecChannel)
        implicit none
        type(GWPnD_t), intent(in)               :: paragwp
        type(Basis_t), intent(inout)            :: Basis
        type(psi_t)   ,intent(inout)            :: psi0_nD
        type(psi_t),allocatable                 :: GWnD(:)
        integer           , intent(in)          :: I_ElecChannel
        integer                                 :: i
        real(kind= Rk)                          :: dot_prdct

        allocate(GWnD(paragwp%ndim))
        psi0_nD%CVec(:)=CONE
        do i = 1, paragwp%ndim,1
            CALL init_psi(GWnD(i)  ,Basis,    cplx=.TRUE.   ,grid =.true.)
            CALL GWP_1D(paragwp%nD(i), GWnD(i),Basis,I_ElecChannel)
            psi0_nD%CVec(:) = psi0_nD%CVec(:)*GWnD(i)%CVec(:)
        end do

    END SUBROUTINE GWP_nD



    end module GWPnD_m
