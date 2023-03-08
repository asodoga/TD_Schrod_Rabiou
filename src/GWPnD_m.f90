module GWPnD_m
    USE NumParameters_m
    USE Basis_m
    USE GWP1D_m
   ! USE NDindex_m
    USE psi_m
    !Use propa_m
    implicit none


TYPE  ::  GWPnD_t ! gaussian nD (ndim)
     integer                         :: ndim
     type(GWP1D_t),allocatable       :: nD(:) ! ndim
END TYPE   GWPnD_t
public::  Read_GWPnD,Write_GWPnD ,Tab_GWPnD,psi0_init1,psi0_init


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
        psi0nD = CONE
        do inb = 1, paragwp%ndim,1
            CALL GWP01D(paragwp%nD(inb),Q(inb), TabGWPnD(inb))
            psi0nD = psi0nD *TabGWPnD(inb)
        end do

        END SUBROUTINE Tab_GWPnD

    SUBROUTINE GWP_init(psi0,Istate,nio)
        type(psi_t), intent(inout)        ::psi0
        type(psi_t)  ,target              ::psi
        complex(Kind = Rk), allocatable   ::g(:)
        type(GWPnD_t)                     :: paragwp
        integer, intent(in)               :: nio,Istate
        real(Kind= Rk), allocatable       ::Q(:,:)
        complex(Kind= Rk), allocatable    ::gb(:,:)
        real(Kind= Rk)                    ::NormG,NormB
        integer                           ::Iq,inb,Ndim

        Ndim = size(Psi0%Basis%tab_basis)-1
        allocate(g(psi0%Basis%nq))
        allocate(gb(psi0%Basis%nq,Psi0%Basis%tab_basis(Ndim+1)%nb))
        CALL init_psi(psi,   psi0%Basis,    cplx=.TRUE.   ,grid  =.true.)

        call    Read_GWPnD(paragwp,nio,Psi0%Basis)
        call    Write_GWPnD(paragwp)

        call calc_Q_grid(Q,Psi0%Basis)
        DO iq = 1,psi0%Basis%nq
            call Tab_GWPnD(paragwp,Q(Iq,:),g(Iq))
        END DO
        gb(:,:) = CZERO
        gb(:,Istate) =  g(:)
        psi%CVec(:) = reshape(gb,shape=[psi0%Basis%nq*Psi0%Basis%tab_basis(Ndim+1)%nb])
        call Calc_Norm_OF_Psi(psi,NormG)
        psi%CVec(:) = psi%CVec(:)/NormG
        DO iq = 1,psi0%Basis%nq
            write(22,*) Q(Iq,:),abs(psi%CVec(Iq))**2
        END DO
        call Calc_Norm_OF_Psi(psi,NormG)
        print*,'NormG = ',NormG
        call GridTOBasis_nD_cplx(psi0%CVec,psi%CVec,psi0%Basis)
        call Calc_Norm_OF_Psi(psi0,NormB)
        print*,'NormB = ',NormB
        call BasisTOGrid_nD_cplx(psi%CVec,psi0%CVec,psi0%Basis)

        deallocate(g)
        deallocate(Q)
        CALL dealloc_psi(psi)
    END SUBROUTINE GWP_init




    SUBROUTINE psi0_init1(psi0,nio)
        type(psi_t), intent(inout)        ::psi0
        type(psi_t)  ,target              ::psi
        complex(Kind = Rk), allocatable   ::g(:)
        type(GWPnD_t)                     :: paragwp
        integer, intent(in)               :: nio
        real(Kind= Rk), allocatable       ::Q(:,:)
        complex(Kind= Rk), pointer        ::gb(:,:)
        real(Kind= Rk)                    ::Dq,Dphi,NormG,NormB,Q0,DQ1
        integer                           ::iq,inb,Ndim
        call    Read_GWPnD(paragwp,nio,Psi0%Basis)
        call    Write_GWPnD(paragwp)
        Ndim = size(Psi0%Basis%tab_basis)
        Q0 = ZERO
        !open (unit = 19, file = "psi0_Retinal")
        call calc_Q_grid(Q,Psi0%Basis)
        allocate(g(psi0%Basis%nq))
        CALL init_psi(psi,   psi0%Basis,    cplx=.TRUE.   ,grid  =.true.)
        DO iq = 1,psi0%Basis%nq
            DO inb = 1,Ndim-1
                g(iq) = exp(-((Q(iq,inb)-Q0)/paragwp%nD(inb)%DQ)**2 ) / sqrt(sqrt(pi/TWO)*paragwp%nD(inb)%DQ)
            END DO
        END DO
        psi%CVec(:) = CZERO
        gb( 1:psi0%Basis%nq, 1:psi0%Basis%tab_basis(size(psi0%Basis%tab_basis))%nb)   =>    psi%CVec
        gb(:,1) = CZERO
        gb(:,2) = g(:)
        !CALL Write_psi_basis(psi,0._Rk,2)
        call Calc_Norm_OF_Psi(psi,NormG)
        print*,'NormG = ',NormG
        call GridTOBasis_nD_cplx(psi0%CVec,psi%CVec,psi0%Basis)
        !call GTB_nDcplx(Psi0%CVec,psi%CVec,psi0%Basis)
        !CALL Write_psi_basis(psi0,0._Rk,2)
        call Calc_Norm_OF_Psi(psi0,NormB)
        print*,'NormB = ',NormB
        deallocate(g)
        deallocate(Q)
        CALL dealloc_psi(psi)
    END SUBROUTINE psi0_init1

    SUBROUTINE psi0_init(psi0,nio)
        type(psi_t), intent(inout)        ::psi0
        type(psi_t)  ,target              ::psi
        complex(Kind = Rk), allocatable   ::g(:)
        type(GWPnD_t)                     :: paragwp
        integer, intent(in)               :: nio
        real(Kind= Rk), allocatable       ::Q(:,:)
        complex(Kind= Rk), pointer        ::gb(:,:)
        real(Kind= Rk)                    ::Dq,Dphi,NormG,NormB,Q0,DQ1
        integer                           ::iq
        call    Read_GWPnD(paragwp,nio,Psi0%Basis)
        call    Write_GWPnD(paragwp)
        DQ  = 1.474256632_Rk
        DQ1  = sqrt(2.0_Rk)
        Dphi =  0.181095798_Rk
        Q0 = ZERO
        !open (unit = 19, file = "psi0_Retinal")
        call calc_Q_grid(Q,Psi0%Basis)
        allocate(g(psi0%Basis%nq))
        CALL init_psi(psi,   psi0%Basis,    cplx=.TRUE.   ,grid  =.true.)
        DO iq = 1,psi0%Basis%nq
            g(iq) = exp(-   ((Q(iq,1)-Q0)/Dphi)**2 ) / sqrt(sqrt(pi/TWO)*Dphi)
            g(iq) = g(iq) *  exp(-   ((Q(iq,2)-Q0)/DQ)**2 ) / sqrt(sqrt(pi/TWO)*DQ)
            g(iq) = g(iq) *  exp(-   ((Q(iq,3)-Q0)/DQ1)**2 ) / sqrt(sqrt(pi/TWO)*DQ1)
            g(iq) = g(iq) *  exp(-   ((Q(iq,4)-Q0)/DQ1)**2 ) / sqrt(sqrt(pi/TWO)*DQ1)

        END DO
        psi%CVec(:) = CZERO
        gb( 1:psi0%Basis%nq, 1:psi0%Basis%tab_basis(size(psi0%Basis%tab_basis))%nb)   =>    psi%CVec
        gb(:,1) = CZERO
        gb(:,2) = g(:)
        !CALL Write_psi_basis(psi,0._Rk,2)
        call Calc_Norm_OF_Psi(psi,NormG)
        psi%CVec(:) = psi%CVec(:)/NormG
        call Calc_Norm_OF_Psi(psi,NormG)
        print*,'NormG = ',NormG
       call GridTOBasis_nD_cplx(psi0%CVec,psi%CVec,psi0%Basis)
       !call GTB_nDcplx(Psi0%CVec,psi%CVec,psi0%Basis)
        !CALL Write_psi_basis(psi0,0._Rk,2)
        call Calc_Norm_OF_Psi(psi0,NormB)
        print*,'NormB = ',NormB
        deallocate(g)
        deallocate(Q)
       CALL dealloc_psi(psi)
    END SUBROUTINE psi0_init
    end module GWPnD_m
