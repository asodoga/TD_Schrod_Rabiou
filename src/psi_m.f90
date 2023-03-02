module psi_m
  USE UtilLib_m
  USE NumParameters_m
  USE Basis_m
  USE NDindex_m
  implicit none

  TYPE :: psi_t
    type(Basis_t),    pointer      :: Basis
    real (kind=Rk),    allocatable :: RVec(:)
    complex (kind=Rk), allocatable :: CVec(:)
    logical                        :: Grid = .true.
      CONTAINS
          PRIVATE
      PROCEDURE, PASS         :: Copy_psi    ! Copy content from other psi1 to psi2 instance,
      GENERIC,PUBLIC          :: ASSIGNMENT(=)  => Copy_psi

  END TYPE psi_t

   public :: psi_t,write_psi,init_psi,dealloc_psi,write_psi_Grid, Write_p
   public :: write_psi_basis,Calc_Norm_OF_PsiBasis,Calc_Norm_OF_PsiGrid,Calc_Norm_OF_Psi
   public:: Projection,Test_Hagedorn,init_psiHG
contains

     SUBROUTINE copy_psi(psi_out,psi_in)

        CLASS(psi_t), intent(in)     :: psi_in
        CLASS(psi_t), intent(inout)  :: psi_out

        IF (allocated(psi_in%RVec)) THEN
            write(out_unitp,*) 'Coping psi_in in psi_out (real):'
            psi_out%RVec(:) = psi_in%RVec(:)
            !CALL init_Basis1_TO_Basis2 (psi_in%Basis,psi_out%Basis)
            !CALL  construct_primitive_basis(psi_out%Basis)
             psi_out%Basis => psi_in%Basis
            write(out_unitp,*) 'END Coping psi_in in psi_out'
        END IF
        IF (allocated(psi_in%CVec)) THEN
            write(out_unitp,*) 'Coping psi_in in psi_out (complex):'
            psi_out%CVec(:) = psi_in%CVec(:)
            !CALL init_Basis1_TO_Basis2 (psi_in%Basis,psi_out%Basis)
            !CALL  construct_primitive_basis(psi_out%Basis)
            psi_out%Basis => psi_in%Basis
            write(out_unitp,*) 'END Coping psi_in in psi_out'
        END IF

     END SUBROUTINE copy_psi




  SUBROUTINE init_psi(psi,Basis,cplx,grid)
  USE Basis_m

    TYPE(psi_t),    intent(inout)      :: psi
    TYPE (Basis_t), intent(in), target :: Basis
    logical,        intent(in)         :: cplx,grid

    CALL dealloc_psi(psi)

    IF (.NOT. Basis_IS_allocated(Basis)) STOP 'ERROR in init_psi: the Basis is not initialized'


    IF (Basis%nb < 1) STOP 'ERROR in init_psi: Basis%nb < 1!'

    Psi%Basis => Basis
    Psi%Grid = Grid
  If( Grid)THEN !allocation on grid
      !print*,"psi is on Grid"
    IF (cplx) THEN
     IF(allocated(Basis%tab_basis))THEN
      allocate(psi%CVec(Basis%nq*Basis%tab_basis(size(Basis%tab_basis))%nb))
     else
      allocate(psi%CVec(Basis%nq))
     END IF
    ELSE
      IF(allocated(Basis%tab_basis))THEN
        allocate(psi%RVec(Basis%nq*Basis%tab_basis(size(Basis%tab_basis))%nb))
      ELSE
        allocate(psi%RVec(Basis%nq))
      END IF
    END IF
  ELSE ! allocation on basis
      !print*,"psi is on Basis"
    IF (cplx) THEN
     IF(allocated(Basis%tab_basis))THEN
      allocate(psi%CVec(Basis%nb*Basis%tab_basis(size(Basis%tab_basis))%nb))
     else
      allocate(psi%CVec(Basis%nb))
     END IF
    ELSE
      IF(allocated(Basis%tab_basis))THEN
        allocate(psi%RVec(Basis%nb*Basis%tab_basis(size(Basis%tab_basis))%nb))
      ELSE
        allocate(psi%RVec(Basis%nb))
      END IF
    END IF
  END IF
  END SUBROUTINE init_psi

  SUBROUTINE dealloc_psi(psi)
    TYPE(psi_t), intent(inout) :: psi

    nullify(psi%Basis)

    IF (allocated(psi%RVec)) THEN
      deallocate(psi%RVec)
    END IF
    IF (allocated(psi%CVec)) THEN
      deallocate(psi%CVec)
    END IF

  END SUBROUTINE dealloc_psi

  SUBROUTINE write_psi(psi)
    TYPE(psi_t), intent(in) :: psi
    integer                 :: i

    IF (associated(psi%Basis)) THEN
     ! write(out_unitp,*) ' The basis is linked to psi.'
    END IF
      !CALL  Write_Basis(psi%Basis)
    IF (allocated(psi%RVec)) THEN
     ! write(out_unitp,*) 'Writing psi (real):'
      write(out_unitp,*) psi%RVec
     ! write(out_unitp,*) 'END Writing psi'
    END IF
    IF (allocated(psi%CVec)) THEN
       Write(out_unitp,*) 'Writing psi (complex):'
       do i=1, size(psi%CVec)
       write(out_unitp,*) i, psi%CVec(i)
       end do
      ! write(out_unitp,*) psi%CVec
       write(out_unitp,*) 'END Writting psi'
    END IF

  END SUBROUTINE write_psi



    SUBROUTINE Projection(psi_dt_2,psi_dt_1)
        TYPE(psi_t), intent(in)        :: psi_dt_1
        TYPE(psi_t), intent(inout)     :: psi_dt_2

        logical,  parameter            :: debug = .true.
        integer                        :: ib1,ib2

        psi_dt_2%CVec(:) =CZERO
        !write(out_unitp,*) 'writing psi_in'
        !CALL Write_psi(psi_dt_1)
        write(out_unitp,*) 'Begin Hagedorn projection'
        write(out_unitp,*) 'b1'
            CALL Write_RMat(psi_dt_1%Basis%tab_basis(1)%S,out_unitp,5)
            write(out_unitp,*) 'b2'
            CALL Write_RMat(psi_dt_2%Basis%tab_basis(1)%S,out_unitp,5)
           
            DO ib2=1,psi_dt_1%Basis%tab_basis(1)%nb
                DO ib1=1,psi_dt_1%Basis%tab_basis(1)%nb
                    psi_dt_2%CVec(ib2) = psi_dt_2%CVec(ib2) + psi_dt_1%Basis%tab_basis(1)%S(ib1,ib2)*psi_dt_1%CVec(ib1)
                END DO
                
            END DO
            !write(out_unitp,*) 'writing psi_out'
            !CALL Write_psi(psi_dt_2)

            write(out_unitp,*) 'END Hagedorn projection'
            
    END SUBROUTINE Projection


    SUBROUTINE init_psiHG(psi,Basis)
        USE UtilLib_m
        TYPE(psi_t)   ,        INTENT(INOUT)               :: psi
        TYPE(Basis_t)   ,        INTENT(IN),target          :: Basis
        write(out_unitp,*) 'Beging init_psiHG'
        call write_psi(psi)
        psi%Basis  => Basis
        call write_psi(psi)
        write(out_unitp,*) 'END init_psiHG'
    END SUBROUTINE
     SUBROUTINE Write_p(psi,i)
         TYPE(psi_t), intent(in) :: psi
         TYPE(psi_t)             :: psi1
         integer,intent(in)      :: i
         integer                 :: iq ,ib
         if(psi%grid)then
             print*,'pis is on grid'
             CALL init_psi(psi1,psi%Basis,cplx=.TRUE.,grid =.true.)
             call BasisTOGrid_nD_cplx(psi1%CVec,psi%CVec,psi%Basis)
             do iq = 1,psi%Basis%nq
                write(100+i,*)   psi1%Basis%tab_basis(1)%X(iq), abs(psi1%CVec(iq))**2
             end do
         else
             print*,'pis is on basis'
             do ib = 1,psi%Basis%nb
                 write(i,*)  abs(psi%CVec(ib))**2
             end do
         end if
         END SUBROUTINE


    SUBROUTINE Test_Hagedorn(psi1,psi2,q10,q20,sci,scj)
        USE UtilLib_m
        TYPE(psi_t), intent(in)        :: psi1
        TYPE(psi_t), intent(inout)     :: psi2
        real(Kind = Rk),allocatable    :: S(:,:)
        real(kind=Rk)      ,intent(in) :: q10,q20,sci,scj
        integer                        :: iq,jq
        real(kind=Rk)                  :: Norm1,Norm2

        allocate(S(size(psi1%CVec),size(psi1%CVec)))
        print*,''
        Do iq = 1,size(psi1%CVec)
            Do jq =1,size(psi1%CVec)
                CALL  Hermite_product_integral ( S(iq,jq), psi1%Basis%tab_basis(1)%X ,&
                & psi1%Basis%tab_basis(1)%W , iq-1, jq-1,q10,q20,sci,scj )
            End Do
        End Do
        Do iq = 1,size(psi1%CVec)
            print*, S(iq,:)
        End Do
        print*,''

        print*,'psiB1',psi1%CVec
        CALL Calc_Norm_OF_Psi(psi1,Norm1)

        print*, '<psis1|psi1> = ',Norm1
        print*,''
       CALL Projection(psi2,psi1)
        print*,''

        CALL Calc_Norm_OF_Psi(psi2,Norm2)
        print*, '<psis2|psi2> = ',Norm2
        print*,'psiB2',psi2%CVec

        deallocate(S)

    END SUBROUTINE Test_Hagedorn

    SUBROUTINE Hagedorn_transfo(psi1,psi2)
        USE UtilLib_m
        TYPE(psi_t), intent(in)        :: psi1
        TYPE(psi_t), intent(inout)     :: psi2
        real(Kind = Rk),allocatable    :: S(:,:)
        integer                        :: iq,jq
        real(kind=Rk)                  :: Norm1,Norm2

        allocate(S(psi1%Basis%nb,psi1%Basis%nb))
        print*,' Beging Hagedorn projection'
        S(:,:) = ZERO
        psi2%CVec(:) = CZERO
        Do iq = 1,psi1%Basis%nb
            Do jq =1,psi1%Basis%nb
                CALL  Hermite_product_integral ( psi1%Basis%tab_basis(1)%S(iq,jq), psi1%Basis%tab_basis(1)%X ,&
                        &psi1%Basis%tab_basis(1)%W , iq-1, jq-1,psi1%Basis%tab_basis(1)%Q0,psi2%Basis%tab_basis(1)%Q0,&
                        psi1%Basis%tab_basis(1)%SCALEQ,psi2%Basis%tab_basis(1)%SCALEQ )
            End Do
        End Do
        CALL Projection(psi2,psi2)
        deallocate(S)
          print*,'END Hagedorn projection'
    END SUBROUTINE Hagedorn_transfo




    SUBROUTINE Calc_Norm_OF_Psi(Psi,Norm)
        implicit none
        type  (Psi_t),   intent(in)      :: Psi
        real (kind=Rk), intent(inout)    :: Norm
        IF (psi%Grid) THEN
            CALL Calc_Norm_OF_PsiGrid(Psi,Norm)
        ELSE
            CALL Calc_Norm_OF_PsiBasis(psi,Norm)
        END IF
    END SUBROUTINE Calc_Norm_OF_Psi

    SUBROUTINE Calc_Norm_OF_PsiGrid(Psi_g,Norm)

        USE UtilLib_m
        logical,         parameter      :: debug = .false.
        TYPE(Psi_t), intent(in)         :: Psi_g
        complex (kind=Rk),allocatable   :: Psi_gb(:,:)
        logical                         :: Endloop_q
        real(kind=Rk),intent(inout)     :: Norm
        real(kind=Rk),allocatable       :: Norme(:)
        real(kind=Rk)                   :: WnD
        integer,        allocatable     :: Tab_iq(:)
        integer                         :: iq,inb,inbe

        IF (debug) THEN
            write(out_unitp,*) 'Beging NormGrid'
            flush(out_unitp)
        END IF

        Allocate(Psi_gb(Psi_g%Basis%nq,Psi_g%Basis%tab_basis(size(Psi_g%Basis%tab_basis))%nb))
        Allocate(Tab_iq(size(Psi_g%Basis%tab_basis)-1))
        Allocate(Norme(size(Psi_g%Basis%tab_basis)-1))
        Psi_gb(:,:) = reshape(Psi_g%CVec,shape= [psi_g%Basis%nq,Psi_g%Basis%tab_basis(size(Psi_g%Basis%tab_basis))%nb])
        DO inbe = 1,Psi_g%Basis%tab_basis(size(Psi_g%Basis%tab_basis))%nb !electronic state
            Norme(inbe) = ZERO
            Call Init_tab_ind(Tab_iq,Psi_g%Basis%NDindexq)
            Iq = 0
            DO
                Iq = Iq+1
                CALL increase_NDindex(Tab_iq,Psi_g%Basis%NDindexq,Endloop_q)
                !if(inbe == 1)  print*,Tab_iq
                IF (Endloop_q) exit
                WnD= ONE
                DO inb = 1,size(psi_g%Basis%tab_basis)-1
                    WnD  =WnD* Psi_g%Basis%tab_basis(inb)%w(tab_iq(inb))
                END DO
                Norme(inbe) = Norme(inbe) + conjg(Psi_gb(iq,inbe))*Psi_gb(iq,inbe)*WnD
            END DO
            Norm = sqrt(sum(Norme))
        END DO
        Deallocate(Tab_iq)
        Deallocate(Psi_gb)
        IF (debug) THEN
            write(out_unitp,*) 'END NormGrid'
            flush(out_unitp)
        END IF



    END SUBROUTINE Calc_Norm_OF_PsiGrid

    SUBROUTINE Calc_Norm_OF_PsiBasis(Psi,Norm)
        TYPE (psi_t),  intent(in)     :: Psi
        real(kind = Rk),intent(inout) :: Norm
        Norm = sqrt(real(dot_product(Psi%CVec,Psi%CVec), kind=Rk))
        ! write(out_unitp,*) 'norm,psi',Norm
    END SUBROUTINE Calc_Norm_OF_PsiBasis

     SUBROUTINE write_psi_basis(psi,t,nio)
    TYPE(psi_t), intent(in) ,target:: psi
    real(kind=Rk),   intent(in)    :: t
    real(kind=Rk)                  :: c1mod1, c1mod2
    real(Kind= Rk), allocatable    ::Q(:,:)
    integer,   intent(in)          :: nio
    complex(Kind= Rk), pointer     ::gb0(:,:),gb1(:,:)
    integer                        :: i,Ndim
    Ndim = size(Psi%Basis%tab_basis)
      if(psi%Grid) then
          !> *************************wrinting psi on grid***********************************
          gb0( 1:psi%Basis%nq, 1:psi%Basis%tab_basis(Ndim)%nb)   =>    psi%CVec
          call calc_Q_grid(Q,psi%Basis)
          if(nio== 0) then
              do i=1, psi%Basis%nq
                 ! c1mod1 = abs(gb(i,1))**2
                  !if(c1mod1 <= 1.d-10) gb(i,1) =CZERO
                  !c1mod2 = abs(gb(i,2))**2
                  !if(c1mod2 <= 1.d-10) gb(i,2) =CZERO
                  write(*,*) t,Q(i,:), real(gb0(i,1)),aimag(gb0(i,1)),real(gb0(i,2)),aimag(gb0(i,2))
              end do
          else
              do i=1, psi%Basis%nq
                  !c1mod1 = abs(gb(i,1))**2
                  !if(c1mod1 <= 1.d-10) gb(i,1) =CZERO
                  !c1mod2 = abs(gb(i,2))**2
                  !if(c1mod2 <= 1.d-10) gb(i,2) =CZERO
                  write(nio,*) t,Q(i,:), real(gb0(i,1)),aimag(gb0(i,1)) ,real(gb0(i,2)),aimag(gb0(i,2))
              end do
          endif
          else
          !> *************************wrinting psi on Basis***********************************
          gb1( 1:psi%Basis%nb, 1:psi%Basis%tab_basis(Ndim)%nb)   =>    psi%CVec
          if(nio== 0) then
          do i=1, psi%Basis%nb
              !c1mod1 = abs(gb(i,1))**2
              !if(c1mod1 <= 1.d-10) gb(i,1) =CZERO
             ! c1mod2 = abs(gb(i,2))**2
             ! if(c1mod2 <= 1.d-10) gb(i,2) =CZERO
          write(*,*) t, real(gb1(i,1)),aimag(gb1(i,1)),real(gb1(i,2)),aimag(gb1(i,2))
          end do
          else
          do i=1, psi%Basis%nb
              !c1mod1 = abs(gb(i,1))**2
              !if(c1mod1 <= 1.d-10) gb(i,1) =CZERO
             ! c1mod2 = abs(gb(i,2))**2
             ! if(c1mod2 <= 1.d-10) gb(i,2) =CZERO
              write(nio,*) t, real(gb1(i,1)),aimag(gb1(i,1)),real(gb1(i,2)),aimag(gb1(i,2))
          end do
          endif

      end if
  END SUBROUTINE write_psi_basis

    SUBROUTINE write_psi_Grid(psi,nio)
        TYPE(psi_t), intent(in)        :: psi
        TYPE(psi_t)                    :: psi_g
        integer,   intent(in)          :: nio
        integer                         :: i
            !> *************************wrinting psi on grid***********************************

        CALL init_psi(psi_g,   psi%Basis,    cplx=.TRUE.   ,grid =.true.)
        call BasisTOGrid_nD_cplx(psi_g%CVec,psi%CVec,psi%Basis)

           do i=1, psi%Basis%tab_basis(1)%nq
                Write(nio,*) psi_g%Basis%tab_basis(1)%x(i) , ABS(psi_g%CVec(i))**2,&
                 & ABS(psi_g%CVec(psi%Basis%tab_basis(1)%nq+i))**2
           end do
           call dealloc_psi(psi_g)
    END SUBROUTINE write_psi_Grid

    subroutine test(Basis)
        implicit none
        TYPE (Basis_t), target  ,intent(in)     :: Basis
        TYPE (psi_t)                            :: psiB,psiG
        real(kind=Rk)                           :: NormG,NormB
        CALL init_psi(psiB,   Basis,    cplx=.TRUE.   ,grid =.false.)
        CALL init_psi(psiG,   Basis,    cplx=.TRUE.   ,grid =.true.)

        psiB%CVec(:) = CZERO
        psiB%CVec(1) = CONE
        call Calc_Norm_OF_Psi(PsiB,NormB)
        psiB%CVec(:) = psiB%CVec(:)/NormB
           call Calc_Norm_OF_Psi(psiB,NormB)
        call  BasisTOGrid_nD_cplx(PsiG%CVec,psiB%CVec,Basis)
        !call BasisTOGrid_nE_cplx(PsiG%CVec,psiB%CVec,Basis)
        call Calc_Norm_OF_Psi(psiG,NormG)
          print*,'NormB = ',NormB
          print*,'NormG = ',NormG
        print*,'=================================================='
        psiG%CVec(:) = CONE
       ! psiG%CVec(1) = CONE
        call Calc_Norm_OF_Psi(PsiG,NormG)
        psiG%CVec(:) = psiG%CVec(:)/NormG
        call GridTOBasis_nD_cplx(PsiB%CVec,psiG%CVec,Basis)
        call Calc_Norm_OF_Psi(psiB,NormB)
        call Calc_Norm_OF_Psi(psiG,NormG)
        print*,'NormG = ',NormG
        print*,'NormB = ',NormB
    end subroutine test

end module psi_m
