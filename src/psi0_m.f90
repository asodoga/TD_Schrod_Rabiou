module psi0_m
    USE NumParameters_m
    USE Basis_m
    USE GWP1D_m
    USE GWPnD_m
    USE GWPmnD_m
   !implicit none
   !integer                   :: nb_GWP, ndim, I_ElecChannel
    contains


    SUBROUTINE init_psi0(psi0,Basis,I_ElecChannel,ndim, nb_GWP)

        implicit none
        TYPE (Basis_t), intent(inout)       :: Basis
        integer , intent(in)             :: I_ElecChannel,ndim, nb_GWP
        type(psi_t)  , intent(inout)     :: psi0
        type( GWPmnD_t)                  :: paragwp
        real(kind=RK)                    :: dot_prdct

        call Read_GWPmnD(paragwp,in_unitp,ndim,nb_GWP)
        call Write_GWPmnD(paragwp,nb_GWP)
        call  GWP_mnD(paragwp,psi0,Basis,I_ElecChannel,nb_GWP)
        CALL Calc_dot_product(psi0%CVec,dot_prdct,Basis,grid=.true.,yes=.false.)
              psi0%CVec(:)  =  psi0%CVec(:)/SQRT(dot_prdct)
        CALL Calc_dot_product(psi0%CVec,dot_prdct,Basis,grid=.true.,yes=.true.)


    END SUBROUTINE init_psi0


    SUBROUTINE write_psi1(psi,no)
    TYPE(psi_t), intent(in)  :: psi
    integer,   intent(in)    :: no
    integer                   :: i

    IF (associated(psi%Basis)) THEN
      write(out_unitp,*) ' The basis is linked to psi.'
    END IF

    IF (allocated(psi%RVec)) THEN
      write(n0,*) 'Writing psi (real):'
      write(no,*) psi%RVec
      write(no,*) 'END Writing psi'
    END IF
    IF (allocated(psi%CVec)) THEN
      write(no,*) 'Writing psi (complex):'
      do i=1, size(psi%CVec)
      write(no,*) i, psi%CVec(i)
      end do
      write(no,*) 'END Writting psi'
    END IF

  END SUBROUTINE write_psi1












end module psi0_m