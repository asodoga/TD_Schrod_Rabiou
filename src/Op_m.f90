module Op_m
  USE NumParameters_m
  USE Basis_m, only : Basis_t
  USE Molec_m
  Use Psi_m
  implicit none
  private

  TYPE :: Op_t

    TYPE (Basis_t),    pointer     :: Basis

    real (kind=Rk),    allocatable :: RMat(:,:)
  END TYPE Op_t

  public :: Op_t,write_Op,Set_Op,dealloc_Op,calc_OpPsi,Calc_Hpsi  ,Kpsi_nD,Make_Mat_OP

contains
  SUBROUTINE alloc_Op(Op,nb)
    TYPE(Op_t),  intent(inout) :: Op
    integer,     intent(in)    :: nb

    CALL dealloc_Op(Op)

    IF (nb < 1) STOP 'ERROR in init_Op: nb < 1!'

    allocate(Op%RMat(nb,nb))

  END SUBROUTINE alloc_Op

  SUBROUTINE dealloc_Op(Op)
    TYPE(Op_t), intent(inout) :: Op


    nullify(Op%Basis)

    IF (allocated(Op%RMat)) THEN
      deallocate(Op%RMat)
    END IF

  END SUBROUTINE dealloc_Op

  SUBROUTINE write_Op(Op)
    USE UtilLib_m, only : Write_RMat

    TYPE(Op_t), intent(in) :: Op

   ! integer :: i

    IF (associated(Op%Basis)) THEN
      write(out_unitp,*) ' The basis is linked to Op.'
    END IF

    IF (allocated(Op%RMat)) THEN
      write(out_unitp,*) 'Writing Op (real):'
      write(out_unitp,*)
      CALL Write_RMat(Op%RMat,out_unitp,5,name_info='Op%Rmat')
      write(out_unitp,*) 'END Writing Op'
    END IF

  END SUBROUTINE write_Op


    SUBROUTINE Make_Mat_OP(Op)
        USE Basis_m

        USE Psi_m
        TYPE (Op_t),     intent(inout)      :: Op
        logical,          parameter        :: debug = .true.
        !logical,         parameter          :: debug = .false.
        integer                             :: ib,iq,jb
        TYPE(psi_t)                         :: Psi_g
        TYPE(psi_t)                        :: Psi_b
        TYPE(psi_t)                        :: OpPsi_g, OpPsi_b

        IF (debug) THEN
            write(out_unitp,*) 'BEGINNING Make_Mat_OP'
            call write_Op(Op)
            call write_basis(Op%Basis)
            flush(out_unitp)
        END IF
        CALL init_psi(psi_g   ,Op%Basis,    cplx=.TRUE.   ,grid =.true.)
        CALL init_psi(psi_b   ,Op%Basis,    cplx=.TRUE.   ,grid =.false.)
        CALL init_psi(OpPsi_g  ,Op%Basis,    cplx=.TRUE.   ,grid =.true.)
        CALL init_psi(OpPsi_b ,Op%Basis,    cplx=.TRUE.   ,grid =.false.)

        DO ib=1,Op%Basis%nb
            Psi_b%CVec(:)  = ZERO
            Psi_b%CVec(ib) = ONE
           ! stop " coucou"
            CALL BasisTOGrid_nD_cplx(Psi_g%CVec, Psi_b%CVec,Op%Basis)
            CALL Calc_Hpsi(Psi_g%CVec,OpPsi_g%CVec,Op%Basis)
            CALL  GridTOBasis_nD_cplx(OpPsi_b%CVec, OpPsi_g%CVec,Op%Basis)
            Op%RMat(:,ib) = real(OpPsi_b%CVec(:),kind= Rk)
        END DO
        !deallocate(Psi_b)
        !deallocate(Psi_g)
        !deallocate(OpPsi_g)

        IF (debug) THEN
            write(out_unitp,*) 'END Make_Mat_OP'
            flush(out_unitp)
        END IF

    END SUBROUTINE Make_Mat_OP



    SUBROUTINE Set_Op(Op,Basis)
   USE Basis_m

   TYPE(Op_t),     intent(inout)       :: Op
   TYPE (Basis_t), intent(in),  target :: Basis
   REAL(kind=Rk)  ,ALLOCATABLE      :: d0bgw(:,:)

   IF (.NOT. Basis_IS_allocated(Basis)) THEN
     STOP 'ERROR in Set_Op: the Basis is not initialized'
   END IF
   CALL alloc_Op(Op,Basis%nb*Basis%tab_basis(size(Basis%tab_basis))%nb)

   Op%Basis => Basis

 END SUBROUTINE Set_Op




  SUBROUTINE calc_OpPsi(Op,Psi,HPsi)
    USE psi_m

    TYPE(Op_t),  intent(in)          :: Op
    TYPE(psi_t), intent(in)          :: Psi 
    TYPE(psi_t), intent(inout)       :: HPsi
    TYPE(psi_t)                      :: Psi_g  ,HPsi_g


    IF (allocated(Op%RMat)) THEN
        print*,"The Hamilitonian is allocated"
        !> This part is used when the hamilitonian is constructed
        !> else the the subroutine will use Calc_HPsi either psi in on grid or not
          HPsi%CVec = matmul(Op%RMat,Psi%CVec)

        !> This part evaluate H|Psi> using Cal_HPsi
        !> When Psi is on Grid  H|Psi> will be evaluate directly using  Cal_HPsi
        !> Elsif Psi is on Basis, the evaluation will be in thre part
        !> (I)   Psi_b => Psi_g
        !> (II)  Hpsi_g = H.psi_g
        !> (III) Hpsi_g => Hpsi_b
    ELSE
         !print*,"The Hamilitonian is not allocated"
        IF (psi%Grid) THEN
            CALL init_psi(HPsi_g,   Psi%Basis,.true.,.true.)
            CALL Calc_Hpsi(Psi%CVec,Hpsi_g%CVec,Psi%Basis)
            CALL  GridTOBasis_nD_cplx(HPsi%CVec,HPsi_g%CVec,psi%Basis)
           ! CALL GTB_nDcplx(HPsi%CVec,HPsi_g%CVec,psi%Basis)
            CALL dealloc_psi(HPsi_g)
        ELSE
           CALL init_psi(Psi_g,Psi%Basis,.true.,.true.)
           CALL init_psi(HPsi_g,   Psi%Basis,.true.,.true.)
           CALL BasisTOGrid_nD_cplx(Psi_g%CVec,Psi%CVec,Psi%Basis)
           CALL Calc_Hpsi(Psi_g%CVec,Hpsi_g%CVec,Psi%Basis)
           CALL GridTOBasis_nD_cplx(HPsi%CVec,HPsi_g%CVec,psi%Basis)
           CALL dealloc_psi(psi_g)
           CALL dealloc_psi(HPsi_g)
        END IF
    END IF

  END SUBROUTINE calc_OpPsi
  SUBROUTINE Calc_Hpsi(psi_g,HPsi_g,Basis)
        USE Basis_m
        USE Psi_m
        USE Molec_m

        complex(kind=Rk),intent(in) ,target           :: Psi_g(:)
        complex(kind=Rk),intent(inout)                :: HPsi_g(:)
        type (Basis_t), intent(in),  target           :: Basis
        complex(kind=Rk),allocatable,target           ::VPsi_g(:),KPsi_g(:)
        complex(kind=Rk), pointer                     :: VPsi_gb(:,:),Psi_gb(:,:)
        real(kind=Rk), allocatable                    :: V(:,:,:), Q(:,:)
        INTEGER                                       ::iq,i2,j2,i1 ,Ndim
         !open(11, file = 'Pot_Retinal11.dat')
         !open(12, file = 'Pot_Retinal22.dat')
         !open(13, file = 'Pot_Retinal12.dat')
         !open(14, file = 'Pot_Retinal21.dat')
        IF(.not. allocated(Basis%tab_basis)) THEN
            STOP 'ERROR in Set_Op: the Basis%tab_bais is not initialized'
        END IF
        ! action potential V|Psi_g>
          Ndim = size(Basis%tab_basis)
          allocate(VPsi_g(Basis%nq*Basis%tab_basis(Ndim)%nb))
          allocate(KPsi_g(Basis%nq*Basis%tab_basis(Ndim)%nb ))
          allocate(V(Basis%nq,Basis%tab_basis(Ndim)%nb,Basis%tab_basis(Ndim)%nb))
          V(:,:,:)= ZERO
          VPsi_g(:)= CZERO
          KPsi_g(:)= CZERO
          HPsi_g(:)= CZERO
          call Calc_Q_grid(Q,Basis)
          Do iq=1,Basis%nq

              CALL sub_pot(V(iq,:,:),Q(iq,:),1)
             ! if(mod(iq,1000)==0)then
             !     Write(11,*) "  "
             !     Write(12,*) "  "
             !     Write(13,*) "  "
             !     Write(14,*) "  "
             ! else
                !Write(110,*) Q(iq,:),V(iq,1,1)
                !if(mod(iq,25)==0) Write(110,*)
                
             !     Write(12,*) Q(iq,:),V(iq,2,2)
             !     Write(13,*) Q(iq,:),V(iq,1,2)
             !     Write(14,*) Q(iq,:), V(iq,2,1)
             ! end if



          END DO
           VPsi_gb(1:Basis%nq,1:Basis%tab_basis(Ndim)%nb)  =>    VPsi_g
           Psi_gb( 1:Basis%nq, 1:Basis%tab_basis(Ndim)%nb)   =>    Psi_g


          DO i2=1,Basis%tab_basis(Ndim)%nb
             DO j2=1,Basis%tab_basis(Ndim)%nb
                 VPsi_gb(:,i2) = VPsi_gb(:,i2) + V(:,i2,j2)*Psi_gb(:,j2)
             END DO

          END DO
          ! action potential K|Psi_g>
          call Kpsi_nD(KPsi_g,Psi_g,Basis)
          HPsi_g(:) =  VPsi_g(:)+KPsi_g(:)

          DEALLOCATE(VPsi_g)
          DEALLOCATE(KPsi_g)
          DEALLOCATE(Q)
        DEALLOCATE(V)
  END SUBROUTINE Calc_Hpsi

  SUBROUTINE Kpsi_nD(KPsi_g,Psi_g,Basis)
      USE Basis_m
      USE UtilLib_m
      USE Molec_m
      TYPE(Basis_t),        intent(in),target      :: Basis
      complex (kind=Rk), intent(in) ,target        :: Psi_g(:)
      complex (kind=Rk), intent(inout),target      :: KPsi_g(:)
      complex (kind=Rk), pointer                   :: Psi_ggb(:,:,:)
      real (kind=Rk),    pointer                   :: d2gg(:,:)
      complex (kind=Rk), pointer                   :: KPsi_ggb(:,:,:)
      real(kind=Rk), allocatable                   :: GGdef(:,:)
      real(Kind= Rk)                               ::Mass_qphi(2)
      logical,           parameter                 :: debug = .true.
      integer                                      :: iq,i1,i3,inb,Ndim
      integer , allocatable                        :: Iq1(:),Iq2(:),Iq3(:),Ib1(:),Ib2(:),Ib3(:)

      IF (debug) THEN
        !write(out_unitp,*) 'BEGINNING Kpsi'
        flush(out_unitp)
      END IF

      Ndim = size(Basis%tab_basis)
      allocate(GGdef(Ndim-1,Ndim-1))
      CALL get_Qmodel_GGdef(GGdef)
      call Calc_index( Ib1,Ib2,Ib3,Iq1,Iq2,Iq3,Basis)
       KPsi_g(:) = CZERO
      DO inb = 1,Ndim-1
        KPsi_ggb(1:Iq1(inb),1:Iq2(inb),1:Iq3(inb))=> KPsi_g
        Psi_ggb(1:Iq1(inb),1:Iq2(inb),1:Iq3(inb))  => Psi_g
        d2gg(1:Iq2(inb),1:Iq2(inb)) => Basis%tab_basis(inb)%d2gg
        DO i3=1,ubound(Psi_ggb,dim=3)
           DO i1=1,ubound(Psi_ggb,dim=1)
              KPsi_ggb(i1,:,i3) =KPsi_ggb(i1,:,i3)-HALF*GGdef(inb,inb)*matmul(d2gg,Psi_ggb(i1,:,i3))
              ! KPsi_ggb(i1,:,i3) =KPsi_ggb(i1,:,i3)-HALF*matmul(d2gg,Psi_ggb(i1,:,i3))
           END DO
        END DO
      END DO
      Deallocate(Ib1,Ib2,Ib3,Iq1,Iq2,Iq3)
      IF (debug) THEN
       !	write(out_unitp,*) 'END KPsi_nD'
       	flush(out_unitp)
      END IF
  END SUBROUTINE  Kpsi_nD

end module Op_m
