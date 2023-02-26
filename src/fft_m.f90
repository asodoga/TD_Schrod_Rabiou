module fft_m
    USE NumParameters_m
    USE Basis_m
    USE NDindex_m
    USE psi_m
    implicit none
    private
    public:: fft

contains
    SUBROUTINE fft_d(psi,psi_chp,sign0)
        implicit none
   !################################variables declaration############################################
        type (Psi_t) ,intent(in)              ,target           :: Psi
        complex (kind=Rk), pointer                              :: Psi_ggb(:,:)  , Psi_ggb_chp(:,:)
        type (Psi_t) ,intent(inout)              ,target        :: Psi_chp
        integer                                                 :: inb ,i,k
        integer ,intent(in)                                     :: sign0
        logical,          parameter                             :: debug = .true.
        integer                                                 :: iq,i1,i3,inb,ndim
        integer , allocatable                                   :: Iq1(:),Iq2(:),Iq3(:),Ib1(:),Ib2(:),Ib3(:)

  !================================debuging==================================================================
        IF (debug) THEN
            !write(out_unitp,*) 'BEGINNING Average_Q'
            flush(out_unitp)
        END IF
  !=======================================initialisation=========================================================
        CALL init_psi(psi_chp, psi%Basis,    cplx=.TRUE.   ,grid =.true.)
        call Calc_iqib( Ib1,Ib2,Ib3,Iq1,Iq2,Iq3,Ndim,Psi%Basis)
   ! #################################################evaluation of S^[i]############################################

        DO inb = 1,ndim-1
            Psi_ggb(1:Iq1(inb),1:Iq2(inb),1:iq3(inb))  => Psi%CVec
            Psi_ggb_chp(1:Iq1(inb),1:Iq2(inb),1:iq3(inb)) => Psi_chp%CVec

            DO J=1,psi%Basis%tab_basis(inb)%nq
                DO k=1,psi%Basis%tab_basis(inb)%nq
                    Psi_ggb_chp(i1,J,i3) = Psi_ggb_chp(i1,J,i3)+&
                    &Psi(i1,:,i3)*(exp(sign0*EYE*TWO*J*k*(psi%Basis%tab_basis(inb)%nq)**-1))
                END DO
            END DO
        END DO


    END SUBROUTINE fft_d











end module fft_m