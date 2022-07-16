module GWPnD_m
    USE NumParameters_m
    USE Basis_m
    USE GWP1D_m
    USE NDindex_m
    USE psi_m
    implicit none


TYPE  ::  GWPnD_t ! gaussian nD (ndim)
     integer                         :: ndim
     type(GWP1D_t),allocatable       :: nD(:) ! ndim
END TYPE   GWPnD_t
public::  Read_GWPnD,Write_GWPnD ,Tab_GWPnD,init_psi0_nD


contains
SUBROUTINE Read_GWPnD(paragwp,nio,Basis)
   USE NumParameters_m
    implicit none
    integer           , intent(in)                  :: nio
    logical           ,parameter                    :: debug = .true.
    type(Basis_t), intent(in)                       :: Basis
    type( GWPnD_t)    ,intent(inout)                :: paragwp
   integer                                          :: i

       paragwp%ndim = size(Basis%tab_basis)-1
       allocate(paragwp%nD(paragwp%ndim))
   do i = 1, paragwp%ndim,1
       CALL Read_GWP1D(paragwp%nD(i),nio)
   end do

END SUBROUTINE Read_GWPnD

    SUBROUTINE Write_GWPnD(paragwp)
   USE NumParameters_m
    implicit none
    type( GWPnD_t)    ,intent(inout)                  :: paragwp
     integer                                          :: i
     write(out_unitp,*) '*************************************************'
     write(out_unitp,*) '*************************************************'
         do i = 1, paragwp%ndim,1
              CALL Write_GWP1D(paragwp%nD(i))

         end do
         write(out_unitp,*) 'ndim =',  paragwp%ndim
     write(out_unitp,*) '*************************************************'
     write(out_unitp,*) '*************************************************'
   END SUBROUTINE Write_GWPnD


    SUBROUTINE Tab_GWPnD(paragwp,Q,psi0nD)
        implicit none
        type(GWPnD_t), intent(in)                              :: paragwp
        complex(kind=Rk) ,allocatable                          :: TabGWPnD(:)
        real(kind=Rk) ,intent(in)                              :: Q(:)
        complex(kind=Rk), intent(inout)                        :: psi0nD
        integer                                                :: inb
        allocate(TabGWPnD(paragwp%ndim))
        do inb = 1, paragwp%ndim,1
            CALL GWP01D(paragwp%nD(inb),Q(inb), TabGWPnD(inb))
        end do
        psi0nD = product(TabGWPnD)
        END SUBROUTINE Tab_GWPnD

    SUBROUTINE init_psi0_nD(psi0_nD,Basis,nio)
        implicit none
        TYPE (Basis_t)  ,intent(in)                          :: Basis
        complex(Kind = Rk), intent(inout)                    ::psi0_nD(:)
         real(Kind= Rk), allocatable                         ::Q(:,:)
        real(Kind= Rk)                                       ::Norm
         type(GWPnD_t)                                       :: paragwp
        complex(kind=Rk)                                     :: psi0nD
        integer           , intent(in)                       :: nio
        integer                                              ::iq
        call calc_Q_grid(Q,Basis)
        call    Read_GWPnD(paragwp,nio,Basis)
        call    Write_GWPnD(paragwp)
          DO   iq = 1, Basis%nq,1
          call Tab_GWPnD(paragwp,Q(iq,:),psi0_nD(iq))
          END DO

    END SUBROUTINE init_psi0_nD


    end module GWPnD_m
