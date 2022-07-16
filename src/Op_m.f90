module Op_m
  USE NumParameters_m
  USE Basis_m, only : Basis_t
  USE Molec_m
  implicit none
  private

  TYPE :: Op_t

    TYPE (Basis_t),    pointer     :: Basis

    real (kind=Rk),    allocatable :: RMat(:,:)
  END TYPE Op_t

  public :: Op_t,write_Op,Set_Op,dealloc_Op,calc_OpPsi,Calc_Hpsi  ,Kpsi_nD

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



  SUBROUTINE Set_Op(Op,Basis)
   USE Basis_m
   USE Molec_m

   TYPE(Op_t),     intent(inout)       :: Op
   TYPE (Basis_t), intent(in),  target :: Basis
   REAL(kind=Rk)  ,ALLOCATABLE      :: d0bgw(:,:)


   integer :: ib,iq,i2,j2

   real (kind=Rk), allocatable :: mat_pot_grid(:,:,:),OpPsi_g(:),OpPsi_ge(:,:),Q(:),V(:,:)
    REAL(kind=Rk)              :: Q1

   IF (.NOT. Basis_IS_allocated(Basis)) THEN
     STOP 'ERROR in Set_Op: the Basis is not initialized'
   END IF

   CALL alloc_Op(Op,Basis%nb)

   Op%Basis => Basis

   ! calculation of a potential on the grid
       allocate(d0bgw(Basis%nq,Basis%nb))
       allocate(Q(1))
       allocate(mat_pot_grid(Basis%tab_basis(1)%nq , Basis%tab_basis(2)%nb , Basis%tab_basis(2)%nb))
       allocate( V( Basis%tab_basis(2)%nb , Basis%tab_basis(2)%nb)  )

         Do iq=1,Basis%tab_basis(1)%nq
            Q=Basis%tab_basis(1)%x(iq)
            CALL sub_pot(mat_pot_grid(iq,:,:),Q,1)
         END DO

   ! ==============================END PRINT V(:,1,1)
   do iq = 1,Basis%tab_basis(1)%nq
           write(25,*) Basis%tab_basis(1)%x(iq),mat_pot_grid(iq,1,1)
   end do
      allocate(  OpPsi_g(Basis%nb) ,  OpPsi_ge( Basis%tab_basis(1)%nq,Basis%tab_basis(2)%nb ) )
      OpPsi_g(:)   = 0
      OpPsi_ge( :,: ) = 0
      ! calculation of Op|b_i>
    DO ib=1,Basis%nb
        DO i2=1,Basis%tab_basis(2)%nb
        DO j2=1,Basis%tab_basis(2)%nb
            DO iq=1,Basis%nq
             OpPsi_ge(iq,i2) = OpPsi_ge(iq,i2)+mat_pot_grid(iq,i2,j2)*Basis%tab_basis(1)%d0gb(iq,ib)
             !OpPsi_ge(iq,i2) = OpPsi_ge(iq,i2) + V(i2,j2) * Basis%tab_basis(1)%d0gb(iq,ib)
            END DO

        END DO

        END DO

      OpPsi_g(:)= reshape( OpPsi_ge,[Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nb])
          DO iq=1,Basis%nq
             OpPsi_g(iq) = OpPsi_g(iq) - HALF/mass * Basis%tab_basis(1)%d2gb(iq,ib,1,1) ! -1/2mass d2./dx2 part
          END DO
      !  CALL TRANSPOS(d0bgw,Basis%tab_basis(1))
      Op%RMat(:,ib) = matmul(Basis%d0bgw,OpPsi_g)
    END DO

   CALL write_Op(Op)


 END SUBROUTINE Set_Op




  SUBROUTINE calc_OpPsi(Op,Psi,OpPsi)
    USE psi_m, ONLY : psi_t

    TYPE(Op_t),  intent(in)     :: Op
    TYPE(psi_t), intent(in)     :: Psi
    TYPE(psi_t), intent(inout)  :: OpPsi


    IF (.NOT. allocated(Op%RMat)) THEN
      STOP 'ERROR in calc_OpPsi: Op is not initialized!'
    END IF

    IF (allocated(Psi%RVec)) THEN
      OpPsi%RVec = matmul(Op%RMat,Psi%RVec)
    ELSE IF (allocated(Psi%CVec)) THEN
      OpPsi%CVec = matmul(Op%RMat,Psi%CVec)
    ELSE
      STOP 'ERROR in calc_OpPsi: Psi is not initialized!'
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
        IF(.not. allocated(Basis%tab_basis)) THEN
            STOP 'ERROR in Set_Op: the Basis%tab_bais is not initialized'
        END IF
        ! action potential V|Psi_g>
          Ndim = size(Basis%tab_basis)
          allocate(VPsi_g(Basis%nq*Basis%tab_basis(Ndim)%nb))
          allocate(KPsi_g(Basis%nq*Basis%tab_basis(Ndim)%nb ))
          allocate(V(Basis%nq,Basis%tab_basis(Ndim)%nb,Basis%tab_basis(Ndim)%nb))
          call Calc_Q_grid(Q,Basis)
          Do iq=1,Basis%nq
              CALL sub_pot(V(iq,:,:),Q(iq,:),0)
          END DO
           VPsi_gb(1:Basis%nq,1:Basis%tab_basis(Ndim)%nb)  =>    VPsi_g
           Psi_gb(1:Basis%nq,1:Basis%tab_basis(Ndim)%nb)   =>    Psi_g

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
      logical,           parameter                 :: debug = .true.
      integer                                      :: iq,i1,i3,inb,ndim
      integer , allocatable                        :: Iq1(:),Iq2(:),Iq3(:),Ib1(:),Ib2(:),Ib3(:)

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING Kpsi'
        flush(out_unitp)
      END IF
      Ndim = size(Basis%tab_basis)
      call Calc_iqib( Ib1,Ib2,Ib3,Iq1,Iq2,Iq3,Ndim,Basis)
       KPsi_g(:) = CZERO
      DO inb = 1,ndim-1
        KPsi_ggb(1:Iq1(inb),1:Iq2(inb),1:iq3(inb))=> KPsi_g
        Psi_ggb(1:Iq1(inb),1:Iq2(inb),1:iq3(inb))  => Psi_g
        d2gg(1:Iq2(inb),1:Iq2(inb)) => Basis%tab_basis(inb)%d2gg
        DO i3=1,ubound(Psi_ggb,dim=3)
           DO i1=1,ubound(Psi_ggb,dim=1)
              KPsi_ggb(i1,:,i3) =   -(HALF/mass)*matmul(d2gg,Psi_ggb(i1,:,i3))
           END DO
        END DO
      END DO
      IF (debug) THEN
       	write(out_unitp,*) 'END KPsi_nD'
       	flush(out_unitp)
      END IF
  END SUBROUTINE  Kpsi_nD

end module Op_m
