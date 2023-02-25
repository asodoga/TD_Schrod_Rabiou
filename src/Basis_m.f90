MODULE Basis_m
  USE NumParameters_m
  USE NDindex_m


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Basis_t,Read_Basis,Write_Basis,Basis_IS_allocated
  PUBLIC :: Calc_dngg_grid,Basis_IS_allocatedtot,Calc_tranpose_d0gb,test_basitogridgridtobasis
  PUBLIC :: GridTOBasis_nD_cplx,BasisTOGrid_nD_cplx
  PUBLIC ::  BasisTOGrid_1D_cplx, GridTOBasis_1D_cplx
  PUBLIC ::   Calc_Q_grid,construct_primitive_basis,init_Basis1_TO_Basis2
  PUBLIC ::Calc_iqib,Scale_Basis,Construct_poly_Hermite , Hermite_product_integral
  PUBLIC ::Deallocate_Basis  ,Hermite_double_product_func,Calc_basis

  TYPE :: Basis_t
    integer                      :: nb_basis   = 0
    integer                      :: nb         = 0
    integer                      :: nq         = 0
    real(kind=Rk)                :: Q0         = 0_Rk
    real(kind=Rk)                :: scaleQ     = 0_Rk
    real(kind=Rk)                :: A          = 0_Rk
    real(kind=Rk)                :: B          = 0_Rk
    character(len=:),allocatable :: Basis_name
    real(kind=Rk),   allocatable :: x(:)
    real(kind=Rk),   allocatable :: w(:)
    real(kind=Rk),   allocatable :: d0gb(:,:)      ! basis functions d0gb(nq,nb)
    real(kind=Rk),   allocatable :: d1gb(:,:,:)    ! basis functions d2gb(nq,nb,1)
    real(kind=Rk),   allocatable :: d1gg(:,:,:)    ! basis functions d2gg(nq,nq,1)
    real(kind=Rk),   allocatable :: d2gb(:,:,:,:)  ! basis functions d2gb(nq,nb,1,1)
    real(kind=Rk),   allocatable :: d2gg(:,:,:,:)  ! basis functions d2gg(nq,nq,1,1)
    real(kind=Rk),   allocatable :: d0bgw(:,:)     ! transpose of basis functions d0gb(nb,nq)
    TYPE(NDindex_t)              :: NDindexq
    TYPE(NDindex_t)              :: NDindexb
    TYPE (Basis_t),  allocatable :: tab_basis(:)   !  for more than one Basis.
    CONTAINS
    PRIVATE
    PROCEDURE, PASS         :: Copy_basis
    GENERIC,PUBLIC          :: ASSIGNMENT(=)             => Copy_basis

  END TYPE Basis_t

CONTAINS

    SUBROUTINE copy_basis(basis_out,basis_in)

        CLASS(Basis_t), intent(in)     :: basis_in
        CLASS(Basis_t), intent(inout)  :: basis_out

            CALL init_Basis1_TO_Basis2 (Basis_in,Basis_out)
            CALL  construct_primitive_basis(Basis_out)
    END SUBROUTINE copy_basis


    RECURSIVE FUNCTION Basis_IS_allocated(Basis) RESULT(alloc)

      TYPE(Basis_t),   intent(in)  :: Basis
      logical                      :: alloc
      integer                      :: i


              IF( Basis%Basis_name == 'el'.AND. Basis%nb >0 )Then
                alloc = .TRUE.
                RETURN
              END IF

      alloc = allocated(Basis%tab_basis)
      IF ( allocated(Basis%tab_basis)) THEN
        Do i=1,size(Basis%tab_basis)
          alloc  = alloc .and. Basis_IS_allocated(Basis%tab_basis(i))
        END DO
      ELSE
        alloc =             allocated(Basis%x)
        alloc = alloc .AND. allocated(Basis%w)
        alloc = alloc .AND. allocated(Basis%d0gb)
        alloc = alloc .AND. allocated(Basis%d1gb)
        alloc = alloc .AND. allocated(Basis%d2gb)
        !alloc = alloc .AND. allocated(Basis%d0bgw)
      END IF
    END FUNCTION Basis_IS_allocated

    RECURSIVE FUNCTION Basis_IS_allocatedtot(Basis) RESULT(alloc)

        TYPE(Basis_t),   intent(in)  :: Basis
        logical                      :: alloc
        integer                      :: i



        alloc = allocated(Basis%tab_basis)
        IF ( allocated(Basis%tab_basis)) THEN
          Do i=1,size(Basis%tab_basis)
            alloc  = alloc .and. Basis_IS_allocated(Basis%tab_basis(i))
          END DO
        ELSE
          alloc =             allocated(Basis%x)
          alloc = alloc .AND. allocated(Basis%w)
          alloc = alloc .AND. allocated(Basis%d0gb)
          alloc = alloc .AND. allocated(Basis%d1gb)
          alloc = alloc .AND. allocated(Basis%d2gb)
          alloc = alloc .AND. allocated(Basis%d1gg)
          alloc = alloc .AND. allocated(Basis%d2gg)

        END IF

      END FUNCTION Basis_IS_allocatedtot
      RECURSIVE SUBROUTINE Write_Basis(Basis)
       USE UtilLib_m

         TYPE(Basis_t),       intent(in)  :: Basis
         integer                          :: i

       !write(out_unitp,*) '---------------------------------------------------------------------'
       !write(out_unitp,*) 'Write_Basis'
       !write(out_unitp,*) "Basis_name=", Basis%Basis_name
       !write(out_unitp,*) "n_basis=", Basis%nb_basis
       !write(out_unitp,*) 'nb,nq',Basis%nb,Basis%nq
       !write(out_unitp,*) "Q0,scaleQ ", Basis%Q0,Basis%scaleQ
       !write(out_unitp,*)  'A,B',Basis%A,Basis%B

       IF (.NOT.allocated(Basis%x)) THEN
            !write(out_unitp,*)' Basis table x is not allocated.'
           ELSE
             CALL Write_RVec(Basis%x,out_unitp,5,name_info='x')
           END IF
           write(out_unitp,*)
           IF (.NOT.allocated(Basis%W)) THEN
             !write(out_unitp,*)' Basis table w is not allocated.'
           ELSE
            CALL Write_RVec(Basis%w,out_unitp,5,name_info='w')
           END IF
         !! write(out_unitp,*)
         ! IF (.NOT.allocated(Basis%d0gb)) THEN
         !   !write(out_unitp,*)' Basis table d0gb is not allocated.'
         ! ELSE
         ! !  CALL Write_RMat(Basis%d0gb,out_unitp,5,name_info='d0gb')
         ! END IF
         !
         ! !write(out_unitp,*)
         !IF (.NOT.allocated(Basis%d0bgw)) THEN
         !   !write(out_unitp,*)' Basis table d0bgw is not allocated.'
         ! ELSE
         !  !CALL Write_RMat(Basis%d0bgw,out_unitp,5,name_info='d0gbw')
         !END IF
         !! write(out_unitp,*)
         ! IF (.NOT.allocated(Basis%d1gb)) THEN
         !   !write(out_unitp,*)' Basis table d1gb is not allocated.'
         ! ELSE
         !  ! CALL Write_RMat(Basis%d1gb(:,:,1),out_unitp,5,name_info='d1gb')
         ! END IF
         ! !write(out_unitp,*)
         ! IF (.NOT.allocated(Basis%d1gg)) THEN
         !  ! write(out_unitp,*)' Basis table d1gb is not allocated.'
         ! ELSE
         !  ! CALL Write_RMat(Basis%d1gg(:,:,1),out_unitp,5,name_info='d1gg')
         ! END IF
         !! write(out_unitp,*)
         ! IF (.NOT.allocated(Basis%d2gb)) THEN
         !  ! write(out_unitp,*)' Basis table d1gb is not allocated.'
         ! ELSE
         !   !CALL Write_RMat(Basis%d2gb(:,:,1,1),out_unitp,5,name_info='d2gb')
         ! END IF
         !! write(out_unitp,*)
         ! IF (.NOT.allocated(Basis%d2gg)) THEN
         !   !write(out_unitp,*)' Basis table d2gg is not allocated.'
         ! ELSE
         !   !CALL Write_RMat(Basis%d2gg(:,:,1,1),out_unitp,5,name_info='d2gg')
         ! END IF
         !
         !  write(out_unitp,*) 'nb_basis',Basis%nb_basis
         IF (allocated(Basis%tab_basis)) THEN
           DO i=1,size(Basis%tab_basis)
               if (Basis%tab_basis(i)%Basis_name /= 'el') CALL Write_Basis(Basis%tab_basis(i))
           END DO
         END IF
         write(out_unitp,*) '--------------------------------------------------------------------------'

       END SUBROUTINE Write_Basis

    RECURSIVE SUBROUTINE Deallocate_Basis(Basis)
        USE UtilLib_m

        TYPE(Basis_t),       intent(inout)  :: Basis
        integer                          :: i

       ! write(out_unitp,*) '********************************************************************'
        !write(out_unitp,*) 'Deallocate_Basis'
        !write(out_unitp,*) "Basis_name=", Basis%Basis_name
        !write(out_unitp,*) "n_basis=", Basis%nb_basis
        !write(out_unitp,*) 'nb,nq',Basis%nb,Basis%nq
        !write(out_unitp,*) "Q0=", Basis%Q0
        !write(out_unitp,*) "scaleQ =", Basis%scaleQ

        IF(allocated(Basis%tab_basis))THEN
            deallocate(Basis%NDindexq%Tab0)
            deallocate(Basis%NDindexb%Tab0)
            deallocate(Basis%tab_basis)
       END IF
        write(out_unitp,*)

        IF (allocated(Basis%x)) THEN
            deallocate(Basis%x)
           ! write(out_unitp,*)' Basis table x is now deallocated.'
        END IF
        write(out_unitp,*)
        IF (allocated(Basis%W)) THEN
            deallocate(Basis%W)
            !write(out_unitp,*)' Basis table w is now deallocated.'
        END IF
        write(out_unitp,*)
        IF (allocated(Basis%d0gb)) THEN
            deallocate(Basis%d0gb)
           ! write(out_unitp,*)' Basis table d0gb is now deallocated.'
        END IF
        write(out_unitp,*)
         IF (allocated(Basis%d0bgw)) THEN
             deallocate(Basis%d0bgw)
           ! write(out_unitp,*)' Basis table d0bgw is now deallocated.'
         END IF
        write(out_unitp,*)
        IF (allocated(Basis%d1gb)) THEN
            deallocate(Basis%d1gb)
           ! write(out_unitp,*)' Basis table d1gb is now deallocated.'
        END IF
        write(out_unitp,*)
        IF (allocated(Basis%d1gg)) THEN
            deallocate(Basis%d1gg)
           ! write(out_unitp,*)' Basis table d1gb is now deallocated.'
        END IF
        write(out_unitp,*)
        IF (allocated(Basis%d2gb)) THEN
            deallocate(Basis%d2gb)
            !write(out_unitp,*)' Basis table d1gb is now deallocated.'
        END IF
        write(out_unitp,*)
        IF (allocated(Basis%d2gg)) THEN
            deallocate(Basis%d2gg)
          !  write(out_unitp,*)' Basis table d2gg is now deallocated.'
        END IF
        IF (allocated(Basis%tab_basis)) THEN
            DO i=1,size(Basis%tab_basis)
                CALL Deallocate_Basis(Basis%tab_basis(i))
            END DO
        END IF
       ! write(out_unitp,*) '********************************************************************'

    END SUBROUTINE Deallocate_Basis

    RECURSIVE SUBROUTINE Read_Basis(Basis,nio)
        USE UtilLib_m
        logical,             parameter      :: debug = .true.
        !logical,             parameter     ::debug = .false.
        TYPE(Basis_t),       intent(inout)  :: Basis
        integer,             intent(in)     :: nio
        integer                             :: err_io,nb,nq,i,j,nb_basis
        character (len=Name_len)            :: name
        real(kind=Rk)                       :: A,B,scaleQ,Q0,d0,d2,X1,W1

        NAMELIST /basis_nD/ name,nb_basis,nb,nq,A,B,scaleQ,Q0
        nb_basis  = 0
        nb        = 0
        nq        = 0
        A         = ZERO
        B         = ZERO
        Q0        = ZERO
        scaleQ    = ONE
        name      = '0'

        read(nio,nml=basis_nD,IOSTAT=err_io)
        write(out_unitp,nml=basis_nD)
        IF (err_io < 0) THEN
            write(out_unitp,basis_nD)
            write(out_unitp,*) ' ERROR in Read_Basis'
            write(out_unitp,*) ' while reading the namelist "basis_nD"'
            write(out_unitp,*) ' end of file or end of record'
            write(out_unitp,*) ' Probably, you forget a basis set ...'
            write(out_unitp,*) ' Check your data !!'
            STOP ' ERROR in Read_Basis: problems with the namelist.'
        END IF
        IF (err_io > 0) THEN
            write(out_unitp,basis_nD)
            write(out_unitp,*) ' ERROR in Read_Basis'
            write(out_unitp,*) ' while reading the namelist "basis_nD"'
            write(out_unitp,*) ' Probably, some arguments of namelist are wrong.'
            write(out_unitp,*) ' Check your data !!'
            STOP ' ERROR in Read_Basis: problems with the namelist.'
        END IF

        IF (nb_basis > 1) THEN
            Basis%Basis_name     = 'Dp'
            Basis%nb_basis  = nb_basis
            CALL string_uppercase_TO_lowercase(Basis%Basis_name)
            allocate(Basis%tab_basis(nb_basis))
            DO i=1,nb_basis
                CALL Read_Basis(Basis%tab_basis(i),nio)
            END DO
            Basis%nb =1
            Basis%nq =1
            DO i=1,nb_basis
                if (Basis%tab_basis(i)%Basis_name == 'el') cycle
                Basis%nb = Basis%nb * Basis%tab_basis(i)%nb
                Basis%nq = Basis%nq * Basis%tab_basis(i)%nq
            END DO
        ELSE
            Basis%nb_basis             =   nb_basis
            Basis%nb                   =   nb
            Basis%nq                   =   nq
            Basis%Q0                   =   Q0
            Basis%SCALEQ               =   SCALEQ
            Basis%A                    =    A
            Basis%B                    =    B
            Basis%Basis_name           =    trim(adjustl(name))
            CALL string_uppercase_TO_lowercase(Basis%Basis_name)
        END IF
    END SUBROUTINE Read_Basis

    RECURSIVE SUBROUTINE construct_primitive_basis(Basis)
        USE UtilLib_m
        logical,             parameter      :: debug = .true.
        !logical,             parameter      ::debug = .false.
        TYPE(Basis_t),       intent(inout)  :: Basis
        integer, allocatable                :: NDend_q(:)
        integer, allocatable                :: NDend_b(:)
        integer                             :: nb,nq,i,j
        character (len=Name_len)            :: name
           ! write(out_unitp,*) ' Begin  construct primitive  Basis '
        IF(allocated(Basis%tab_basis))THEN
            allocate(NDend_q(Basis%nb_basis-1))
            allocate(NDend_b(Basis%nb_basis-1))
            DO i=1,Basis%nb_basis-1
                NDend_q(i)=Basis%tab_basis(i)%nq
                NDend_b(i)=Basis%tab_basis(i)%nb
            END DO
            CALL Init_NDindex(Basis%NDindexq,NDend_q,Basis%nb_basis-1)
            CALL Init_NDindex(Basis%NDindexb,NDend_b,Basis%nb_basis-1)

            DO i=1,Basis%nb_basis
                CALL construct_primitive_basis(Basis%tab_basis(i))
            END DO

        ELSE
            SELECT CASE (Basis%Basis_name)
            CASE('el')
                write(6,*) 'Electronic basis. Electronic state number:',basis%nb
                basis%nq = 0
            CASE ('boxab')
                CALL Construct_Basis_Sin(Basis)
                Basis%Q0      = Basis%A
                Basis%scaleQ  = pi/(Basis%B-Basis%A)
            CASE ('fourier')
                CALL Construct_Basis_Fourier(Basis)
            CASE ('herm','ho')
                CALL Construct_Basis_Ho(Basis)
            CASE default
                STOP 'ERROR  Noting to construct'
            END SELECT
            !  this part wil not have sens for 'el' basis
            CALL Scale_Basis(Basis,Basis%Q0,Basis%scaleQ)
            CALL   Calc_tranpose_d0gb(Basis)
            CALL Calc_dngg_grid(Basis)
            CALL CheckOrtho_Basis(Basis,nderiv=2)
        END IF
       ! write(out_unitp,*) ' End  construct  primitive Basis '
    END SUBROUTINE construct_primitive_basis

  SUBROUTINE Construct_Basis_Sin(Basis) ! sin : boxAB with A=0 and B=pi
 USE UtilLib_m

    TYPE(Basis_t),       intent(inout)  :: Basis
    real(kind=Rk)                       :: dx
    integer                             :: ib,iq,nb,nq

    nb = Basis%nb
    nq = Basis%nq
    dx = pi/nq

    ! grid and weight
    Basis%x = [(dx*(iq-HALF),iq=1,nq)]
    Basis%w = [(dx,iq=1,nq)]

    allocate(Basis%d0gb(nq,nb))
    allocate(Basis%d1gb(nq,nb,1))
    allocate(Basis%d2gb(nq,nb,1,1))


    DO ib=1,nb
      Basis%d0gb(:,ib)     =          sin(Basis%x(:)*ib) / sqrt(pi*HALF)
      Basis%d1gb(:,ib,1)   =  ib    * cos(Basis%x(:)*ib) / sqrt(pi*HALF)
      Basis%d2gb(:,ib,1,1) = -ib**2 * Basis%d0gb(:,ib)

    END DO

    IF (nb == nq) THEN
      Basis%d0gb(:,nb)      = Basis%d0gb(:,nb)      / sqrt(TWO)
      Basis%d1gb(:,nb,:)    = Basis%d1gb(:,nb,:)    / sqrt(TWO)
      Basis%d2gb(:,nb,:,:)  = Basis%d2gb(:,nb,:,:)  / sqrt(TWO)
    END IF

  END SUBROUTINE Construct_Basis_Sin



  SUBROUTINE Construct_Basis_Fourier(Basis) !basis_name fourier[-pi,pi]
    USE UtilLib_m

    TYPE(Basis_t),       intent(inout)  :: Basis
    real(kind=Rk)                       :: dx
    integer                             :: ib,iq,nb,nq,k

    nb = Basis%nb
    nq = Basis%nq
    dx = TWO*pi/nq

  !>*************** grid and weight *************************************8
    Basis%x = [ (iq*dx - dx/2 -pi ,iq=1,nq)]
    Basis%w = [(dx,iq=1,nq)]
  !> **********************allocation **********************************
    allocate(Basis%d0gb(nq,nb))
    allocate(Basis%d1gb(nq,nb,1))
    allocate(Basis%d2gb(nq,nb,1,1))
 !>***************************************************************************
    DO ib=1,nb
        k = int(ib/2)
      if(mod(ib,2)==0)then
          Basis%d0gb(:,ib)     =  sin(Basis%x(:)*k) / sqrt(pi)
          Basis%d1gb(:,ib,1)   =  k * cos(Basis%x(:)*k) / sqrt(pi)
          Basis%d2gb(:,ib,1,1) = -k**2 * Basis%d0gb(:,ib)
 !>*****************************************************************************
      elseif(ib==1)then
          Basis%d0gb(:,ib)     =  ONE/sqrt(TWO*pi)
          Basis%d1gb(:,ib,1)   =  ZERO
          Basis%d2gb(:,ib,1,1) =  ZERO
 !>*****************************************************************************
      else
          Basis%d0gb(:,ib)     =  cos(Basis%x(:)*k) / sqrt(pi)
          Basis%d1gb(:,ib,1)   =  -k   * sin(Basis%x(:)*k) / sqrt(pi)
          Basis%d2gb(:,ib,1,1) = -k**2 * Basis%d0gb(:,ib)
    end if
    END DO
      !**************************************************************************
    IF (Basis%nb == Basis%nq .AND. mod(Basis%nb,2) == 0) THEN
        Basis%d0gb(:,nb)      = Basis%d0gb(:,nb)      / sqrt(TWO)
        Basis%d1gb(:,nb,:)    = Basis%d1gb(:,nb,:)    / sqrt(TWO)
        Basis%d2gb(:,nb,:,:)  = Basis%d2gb(:,nb,:,:)  / sqrt(TWO)
    END IF

  END SUBROUTINE Construct_Basis_Fourier




  SUBROUTINE Construct_Basis_el(Basis) ! 'el' :
  USE UtilLib_m

     TYPE(Basis_t),       intent(inout)  :: Basis
      Basis%nq = 0
      RETURN



  END SUBROUTINE Construct_Basis_el



 SUBROUTINE Construct_Basis_Ho(Basis) ! HO :
 USE UtilLib_m

    TYPE(Basis_t),       intent(inout)  :: Basis
    integer                :: iq,ib


    allocate(Basis%x(Basis%nq))
    allocate(Basis%w(Basis%nq))

    call hercom(Basis%nq, Basis%x(:), Basis%w(:))

    allocate(Basis%d0gb(Basis%nq,Basis%nb))
    allocate(Basis%d1gb(Basis%nq,Basis%nb,1))
    allocate(Basis%d2gb(Basis%nq,Basis%nb,1,1))


    DO iq = 1, Basis%nq
        DO ib = 1, Basis%nb
          CALL Construct_Basis_poly_Hermite_exp(Basis%x(iq),Basis%d0gb(iq,ib),&
           Basis%d1gb(iq,ib,1),Basis%d2gb(iq,ib,1,1), ib-1,.TRUE.)
        END DO
    END DO

 END SUBROUTINE Construct_Basis_Ho

    !*******************************************************************************************************************


     SUBROUTINE Construct_poly_Hermite(poly_Hermite,x,l)
        Implicit none
        real(kind = Rk),intent(inout)        :: poly_Hermite
        real(kind = Rk)                      :: pl0,pl1,pl2,norme
        real(kind = Rk),intent(in)           :: x
        integer,intent(in)                   :: l
        integer                              :: i

        IF ( l < 0 ) THEN
            Write(out_unitp,*) 'Bad arguments in poly_hermite :'
            Write(out_unitp,*) ' l < 0 : ',l
            STOP
        END IF
        norme  =  sqrt(PI)

        IF (l == 0) THEN
            poly_Hermite  = exp ( - 0.5*x*x  )
            poly_Hermite  = poly_Hermite/sqrt(norme)
        ELSE IF (l == 1) THEN
            norme          = norme * TWO
            poly_Hermite   = TWO*x*exp ( - 0.5*x*x  )
            poly_Hermite   =  poly_Hermite /sqrt(norme)
        ELSE

            pl2 = ONE
            pl1 = TWO * x
            norme = norme * TWO

            DO i=2,l
                norme = norme * TWO * i
                pl0 = TWO*( x*pl1 - (i-1)*pl2 )
                pl2 = pl1
                pl1 = pl0
            END DO
             pl0          = pl0*exp ( - 0.5*x*x  )
             poly_Hermite = pl0/sqrt(norme)


        END IF

    END SUBROUTINE Construct_poly_Hermite

   SUBROUTINE Hermite_double_product_func( Hf,x,i,j,q0i,q0j,sci,scj)
       !    Hf(i,x) is the normalized probabilist's Hermite polynomial of degree I.
       !    Hj(j,x) is the normalized probabilist's Hermite polynomial of degree J.
       !    Hf(X) = Hen(I,X) * Hen(J,X)*exp ( -x^2-X*(q0i+q0j)-0.5*(q0i^2+q0j^2) )
       !    Hen(I,X) = Hfi
       !    Hen(J,X)  = Hfj

       USE UtilLib_m
       real(kind=Rk),       intent(in)     :: x,q0i,q0j,sci,scj
       real(kind=Rk),    intent(inout)     :: Hf !f(q)
       real(kind=Rk)                       :: qi,qj
       integer,intent(in)                  :: i,j
       real(kind=Rk)                       ::  Hfi, Hfj


       ! Hfi = He(I,qi) = HermiteH[I,qi/Sqrt[2]] / Sqrt [ 2^I ]
               qi = (x-q0i)/sci
        CALL Construct_poly_Hermite(Hfi,qi,I)

       ! Hfj = He(J,qj) = HermiteH[J,qj/Sqrt[2]] / Sqrt [ 2^J ]
          qj = (x-q0j)/scj
       CALL  Construct_poly_Hermite(Hfj,qj,J)

       !  Hf = Hfi *Hfj = He(J,x)*He(J,x)

        Hf = Hfi *Hfj


   END SUBROUTINE Hermite_double_product_func


    SUBROUTINE Hermite_product_integral ( intfq, x_table , w_table ,i,j,q0i,q0j,sci,scj )
        !    VALUE = integral ( -oo < x < +oo ) f(x)dx
        !    Parameters:
        !    Input, integer ( kind = Rk ) I, J, the polynomial indices.
        !    Output, real ( kind = Rk ) HERMITE_PRODUCT_INTEGRAL, the value of
        !    the integral.
        USE UtilLib_m
        real(kind=Rk)       ,intent(in)     :: q0i,q0j,sci,scj
        real(kind=Rk)       ,intent(inout)  :: intfq !f(q)
        integer,intent(in)                  :: i,j
        real(kind=Rk)      ,intent(in)      :: x_table(:), w_table(:)
        real(kind=Rk)      ,allocatable     :: Hermite_func_table(:)
        integer                             :: iq

        allocate(Hermite_func_table(size(x_table)))
        Do iq= 1,size(x_table)
            CALL  Hermite_double_product_func( Hermite_func_table(iq),x_table(iq),i,j,q0i,q0j,sci,scj)
        End Do

        intfq = dot_product( Hermite_func_table , w_table )
        !print*,'<b|b> =' ,i,j,intfq
        deallocate( Hermite_func_table )

    END SUBROUTINE Hermite_product_integral

     RECURSIVE SUBROUTINE  init_Basis1_TO_Basis2 (Basis1,Basis2)
        USE UtilLib_m
        TYPE(Basis_t), intent(in)           :: Basis1
        TYPE(Basis_t), intent(inout)        :: Basis2
        integer                             :: ib



        IF(allocated(Basis1%tab_basis))THEN
            call Deallocate_Basis(Basis2)
            Basis2%Basis_name     = Basis1%Basis_name
            Basis2%nb_basis  = Basis1%nb_basis
            allocate(Basis2%tab_basis(Basis2%nb_basis))
            DO ib=1,Basis1%nb_basis
              CALL  init_Basis1_TO_Basis2(Basis1%tab_basis(ib),Basis2%tab_basis(ib))
            END DO
            Basis2%nb =1
            Basis2%nq =1
            DO ib=1,Basis1%nb_basis
                if (Basis2%tab_basis(ib)%Basis_name == 'el') cycle
                Basis2%nb = Basis2%nb * Basis2%tab_basis(ib)%nb
                Basis2%nq = Basis2%nq * Basis2%tab_basis(ib)%nq
            END DO
        ELSE
            Basis2%Basis_name           = Basis1%Basis_name
            Basis2%nb_basis             =  Basis1%nb_basis
            Basis2%nb                   =  Basis1%nb
            Basis2%nq                   =  Basis1%nq
            Basis2%Q0                   =  Basis1%Q0
            Basis2%SCALEQ               =  Basis1%SCALEQ
            Basis2%A                    =  Basis1%A
            Basis2%B                    =  Basis1%B
        END IF

    END SUBROUTINE init_Basis1_TO_Basis2 


    !******************************************************************************************************************




 FUNCTION poly_Hermite(x,l)
 Implicit none
    real(kind = Rk):: poly_Hermite
    real(kind = Rk):: pl0,pl1,pl2,norme,x
    integer        :: i,l

    IF ( l < 0 ) THEN
       Write(out_unitp,*) 'Bad arguments in poly_hermite :'
       Write(out_unitp,*) ' l < 0 : ',l
       STOP
    END IF
       norme  =  sqrt(PI)

    IF (l == 0) THEN
       poly_Hermite = ONE/sqrt(norme)
    ELSE IF (l == 1) THEN
       norme = norme * TWO
       poly_Hermite = TWO * x/sqrt(norme)
    ELSE

       pl2 = ONE
       pl1 = TWO * x
       norme = norme * TWO

     DO i=2,l
       norme = norme * TWO * i
       pl0 = TWO*( x*pl1 - (i-1)*pl2 )
       pl2 = pl1
       pl1 = pl0
     END DO
       poly_Hermite = pl0/sqrt(norme)
     END IF

 END FUNCTION poly_Hermite

 FUNCTION gamma_perso(n)
 Implicit none
    real(kind = Rk)  :: gamma_perso
    real(kind = Rk)  :: a
    integer          :: i,n

    IF (n <= 0) THEN
       write(out_unitp,*) 'ERROR: gamma( n<=0)',n
       STOP
    END IF
      a = ONE
    DO i = 1,n-1
      a = a * dble (i)
    END DO
      gamma_perso = a

 END FUNCTION gamma_perso

 SUBROUTINE herrec ( p2, dp2, p1, x, nq )
 Implicit none
    integer       ::i
    integer       :: nq
    real(kind = Rk):: dp0,dp1,dp2,p0,p1,p2,x

      p1  = ONE
      dp1 = ZERO

      p2  = x
      dp2 = ONE

    DO i = 2, nq

        p0  = p1
        dp0 = dp1

        p1  = p2
        dp1 = dp2

        p2  = x * p1 - HALF * ( dble ( i ) - ONE ) * p0
        dp2 = x * dp1 + p1 - HALF * ( dble ( i ) - ONE ) * dp0

   END DO
 END SUBROUTINE herrec

  SUBROUTINE herroot ( x, nq, dp2, p1 )
  Implicit none
     integer          :: i
     integer          :: nq
     real(kind = Rk),parameter  :: eps = TEN**(-TWELVE) ! 1.0d-12
     real(kind = Rk)  :: d,dp2,p1,p2,x

    DO i = 1, 10
        CALL herrec ( p2, dp2, p1, x, nq )
        d = p2 / dp2
        x = x - d
        IF ( ABS ( d ) <= eps * ( ABS ( x ) + ONE ) ) THEN
          RETURN
        END IF

    END DO

 END SUBROUTINE herroot

 SUBROUTINE hercom (nq,xp,w)
 Implicit none
   integer        :: i,nq
   real(kind = Rk):: cc,dp2,p1,s,temp,x
   real(kind = Rk):: w(nq),xp(nq)

   CC = 1.7724538509_Rk * gamma_perso(nq ) / ( TWO**( nq-1) )

   S = ( TWO * dble (real(nq,Kind=Rk) ) + ONE )**( SIXTH )

   DO i = 1, ( nq + 1 ) / 2
     IF ( i == 1 ) THEN
       x = s**3 - 1.85575_Rk / s
     ELSE IF ( i == 2 ) THEN
       x = x - 1.14_Rk * ( ( dble ( nq ) )**0.426_Rk ) / x
     ELSE IF ( i == 3 ) THEN
       x = 1.86_Rk * x - 0.86_Rk * xp(1)
     ELSE IF ( i == 4 ) THEN
       x = 1.91_Rk * x - 0.91_Rk * xp(2)
     ELSE
       x = TWO * x - xp(i-2)
     END IF
     CALL herroot ( x,  nq, dp2, p1 )
     xp(i) = x
     W(i) = cc / dp2 / p1
     xp( nq-i+1) = - x
     w( nq-i+1) = w(i)
   END DO
   DO i = 1,  nq/2
     temp = xp(i)
     xp(i) = xp( nq+1-i)
     xp( nq+1-i) = temp
   END DO
   DO i = 1, nq
      w(i) = w(i)*exp(xp(i)*xp(i))
   END DO

 END SUBROUTINE hercom

 SUBROUTINE Construct_Basis_poly_Hermite_exp(x,d0gb,d1gb,d2gb,l,deriv)

   logical        :: deriv
   integer        :: l
   real(kind = RK):: pexp,x,d0gb,d1gb,d2gb

   IF (deriv) THEN
      d0gb = poly_Hermite( x,l)
      IF (l == 0) THEN
        d1gb     = ZERO
        d2gb     = ZERO
      ELSE IF (l == 1) THEN
        d1gb = sqrt(TWO)*poly_Hermite( x,0)
        d2gb = ZERO
      ELSE IF (l == 2) THEN
        d1gb = sqrt(TWO*l) * poly_Hermite( x,l-1)
        d2gb = TWO*( x*d1gb-d0gb *l)
      ELSE
        d1gb = sqrt(TWO*l) * poly_Hermite( x,l-1)
        d2gb = TWO*( x* d1gb-d0gb*l)
      END IF
      pexp = exp(- HALF* x* x)
      d2gb = (d2gb-TWO*x*d1gb+( x* x-ONE)*d0gb)*pexp
      d1gb = (d1gb- x*d0gb)*pexp
      d0gb = d0gb*pexp
    ELSE
     d0gb = poly_Hermite(x ,l)*exp(-HALF* x* x)
     d1gb = ZERO
     d2gb = ZERO
    END IF

  END SUBROUTINE Construct_Basis_poly_Hermite_exp


  SUBROUTINE CheckOrtho_Basis(Basis,nderiv)
  USE UtilLib_m

    TYPE(Basis_t),           intent(in)     :: Basis
    integer,                 intent(in)     :: nderiv

    integer                      :: ib,jb
    real(kind=Rk), ALLOCATABLE   :: S(:,:)
    real(kind=Rk), ALLOCATABLE   :: d0bgw(:,:)
    real(kind=Rk)                :: Sii,Sij


    IF( Basis%Basis_name == 'el')Then
      print*,'This routine is .not. possible Basis el'
      RETURN
    END IF
    IF (Basis_IS_allocated(Basis)) THEN
      d0bgw = transpose(Basis%d0gb)
      DO ib=1,Basis%nb
        d0bgw(ib,:) = d0bgw(ib,:) * Basis%w(:)
      END DO

      S = matmul(d0bgw,Basis%d0gb)
      ! do ib = 1,Basis%nb
          ! do jb = 1,Basis%nb
       !write(115,*) S(ib,ib),S(ib,jb)

             !  end do

          ! end do
      ! IF (nderiv > -1) CALL Write_RMat(S,out_unitp,5,name_info='S')
      Sii = ZERO
      Sij = ZERO
      DO ib=1,Basis%nb
        IF (abs(S(ib,ib)-ONE) > Sii) Sii = abs(S(ib,ib)-ONE)
        S(ib,ib) = ZERO
      END DO
      Sij = maxval(S)
      write(out_unitp,*) 'Sii-1,Sij',Sii,Sij

      IF (nderiv > 0) THEN
        !write(out_unitp,*)
        S = matmul(d0bgw,Basis%d1gb(:,:,1))
        !CALL Write_RMat(S,out_unitp,5,name_info='<d0b|d1b>',Rformat='e13.4')
        !CALL Write_RMat(S,out_unitp,5,name_info='<d0b|d1b>')
      END IF

      IF (nderiv > 1) THEN
        !write(out_unitp,*)
        S = matmul(d0bgw,Basis%d2gb(:,:,1,1))
       ! CALL Write_RMat(S,out_unitp,5,name_info='<d0b|d2b>',Rformat='e13.4')
       ! CALL Write_RMat(S,out_unitp,5,name_info='<d0b|d1b>')
      END IF

    ELSE
      write(out_unitp,*) ' WARNNING in CheckOrtho_Basis'
      write(out_unitp,*) ' the basis is not allocated.'
    END IF

  END SUBROUTINE CheckOrtho_Basis
  SUBROUTINE Scale_Basis(Basis,x0,sx)
  USE UtilLib_m

    TYPE(Basis_t),       intent(inout)  :: Basis
    real(kind=Rk),       intent(in)     :: x0,sx

    IF (Basis%nq == 0) RETURN
    IF (abs(sx) > ONETENTH**6 .AND. Basis_IS_allocated(Basis)) THEN

      Basis%x(:) = x0 + Basis%x(:) / sx
      Basis%w(:) =      Basis%w(:) / sx

      Basis%d0gb(:,:)     = Basis%d0gb(:,:)     * sqrt(sx)
      Basis%d1gb(:,:,:)   = Basis%d1gb(:,:,:)   * sqrt(sx)*sx
      Basis%d2gb(:,:,:,:) = Basis%d2gb(:,:,:,:) * sqrt(sx)*sx*sx
    ELSE
      write(out_unitp,*) ' ERROR in Scale_Basis'
      write(out_unitp,*) ' sx is too small  or ...'
      write(out_unitp,*) ' the basis is not allocated.'
      STOP 'ERROR in Scale_Basis'
    END IF

  END SUBROUTINE Scale_Basis

     SUBROUTINE Calc_basis(Basis2,Basis1,x0,sx)
         USE UtilLib_m
         TYPE(Basis_t)   ,        INTENT(IN)        :: Basis1
         TYPE(Basis_t)   ,        INTENT(INOUT)     :: Basis2
         REAL(kind=Rk),              INTENT(IN)     :: X0,SX
         write(out_unitp,*) 'Beging init_Basis'
        ! CALL Write_Basis(Basis1)
         Call init_Basis1_TO_Basis2 (Basis1,Basis2)
         Basis2%tab_basis(1)%Q0      = X0
         Basis2%tab_basis(1)%SCALEQ  = SX
         Call construct_primitive_basis(Basis2)
         !CALL Write_Basis(Basis2)
         write(out_unitp,*) 'END init_Basis'
         END SUBROUTINE

SUBROUTINE Calc_tranpose_d0gb(Basis)
  USE UtilLib_m
  TYPE(Basis_t)   ,        INTENT(INOUT)     :: Basis
  !real(kind=Rk),        INTENT(INOUT)  :: d0bgw(:,:)
  INTEGER                                 :: IB,IQ
  if(Basis%Basis_name == 'el')THEN
    RETURN
  end if
  ALLOCATE(Basis%d0bgw(Basis%nb,Basis%nq))
    Basis%d0bgw = TRANSPOSE(Basis%d0gb)
    DO IB=1,Basis%nb
      Basis%d0bgw(IB,:) = Basis%d0bgw(IB,:) * Basis%w(:)
    END DO

  END  SUBROUTINE Calc_tranpose_d0gb

  SUBROUTINE Calc_dngg_grid(Basis)
USE UtilLib_m
  TYPE(Basis_t), intent(inout)    :: Basis
  integer                         :: ib
  logical,         parameter     ::debug = .false.


  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING Calc_dngg_grid'
    CALL Write_Basis(Basis)
    flush(out_unitp)
  END IF

  allocate(Basis%d1gg(Basis%nq,Basis%nq,1))
  allocate(Basis%d2gg(Basis%nq,Basis%nq,1,1))

  IF (debug) THEN
    CALL Write_RMat(Basis%d0bgw(:,:),out_unitp,5,name_info='d0bgw')
    write(out_unitp,*)
  END IF

  if(Basis%Basis_name == 'el')THEN
    RETURN
  end if

  Basis%d1gg(:,:,1)   = matmul(Basis%d1gb(:,:,1),Basis%d0bgw)
  Basis%d2gg(:,:,1,1) = matmul(Basis%d2gb(:,:,1,1),Basis%d0bgw)

  IF (debug) THEN
    CALL Write_RMat(Basis%d1gg(:,:,1),out_unitp,5,name_info='d1gg')
    write(out_unitp,*)
    CALL Write_RMat(Basis%d2gg(:,:,1,1),out_unitp,5,name_info='d2gg')
    CALL Write_Basis(Basis)
    write(out_unitp,*) 'END Calc_dngg_grid'
    flush(out_unitp)
  END IF

END SUBROUTINE Calc_dngg_grid


  SUBROUTINE test_basitogridgridtobasis(Basis)
  USE NumParameters_m
  USE UtilLib_m
    TYPE(Basis_t),    intent(in)       :: Basis
    logical,          parameter        :: debug = .true.
    COMPLEX(kind=Rk),    allocatable   :: G1(:),B1(:)!,G(:), Hpsi(:)
    COMPLEX(kind=Rk),    allocatable   :: G2(:),B2(:)
    COMPLEX(kind=Rk),    allocatable   :: B(:)
    REAL(kind=Rk),    allocatable       :: diff_g(:),diff_b(:)
    !REAL(KIND=Rk)                      :: Norm0,Norm1,min_diff,max_diff
    integer                             :: iq,ndim
    ndim = size(Basis%tab_basis)


    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Test'
      flush(out_unitp)
    END IF

    allocate(B(Basis%nb*Basis%tab_basis(ndim)%nb))
    allocate(G1(Basis%nq*Basis%tab_basis(ndim)%nb))
    allocate(diff_g(Basis%nq*Basis%tab_basis(ndim)%nb))
    allocate(diff_b(Basis%nb*Basis%tab_basis(ndim)%nb))
    allocate(G2(Basis%nq*Basis%tab_basis(ndim)%nb))
    allocate(B1(Basis%nb*Basis%tab_basis(ndim)%nb))
    allocate(B2(Basis%nb*Basis%tab_basis(ndim)%nb))

    B(:)=CZERO
    B1(:)=CONE
    G1(:)=CONE
    G2(:)=CZERO
    Call GridTOBasis_nD_cplx(B,G1,Basis)
   Call BasisTOGrid_nD_cplx(G2,B,Basis)
   print*,'##############################################################'

    diff_g(:) = ABS(G1(:)-G2(:))

      print*,'##############################################################'
     Write(out_unitp,*) 'maxval(diff_g(:))=',maxval(diff_g(:))
     Write(out_unitp,*) 'MINVAL(diff_g(:))=',MINVAL(diff_g(:))
       print*,'##############################################################'
      G2(:)=CZERO
      B2(:)=CZERO
      B1(:)=CONE
      Call BasisTOGrid_nD_cplx(G2,B1,Basis)
      Call GridTOBasis_nD_cplx(B2,G2,Basis)
      diff_b(:) = ABS(B2(:)-B1(:))
      print*,'##############################################################'
     Write(out_unitp,*) 'maxval(diff_b(:))=',maxval(diff_b(:))
     Write(out_unitp,*) 'MINVAL(diff_b(:))=',MINVAL(diff_b(:))
       print*,'##############################################################'



    IF (debug) THEN
      write(out_unitp,*) 'END Test'
      flush(out_unitp)
    END IF

END SUBROUTINE test_basitogridgridtobasis





    SUBROUTINE BasisTOGrid_1D_cplx(Psi_ggb,Psi_bbb,Basis)
            USE UtilLib_m
            TYPE(Basis_t)    , intent(in),target     :: Basis
            complex (kind=Rk), intent(inout)         :: Psi_ggb(:,:,:)
            complex (kind=Rk), intent(in)            :: Psi_bbb(:,:,:)
           !complex (kind=Rk), target,allocatable    :: Psi_b(:)
           !complex (kind=Rk), target  ,allocatable  :: Psi_b(:)
            logical          , parameter             :: debug = .true.
            integer                                  ::i1,i3


            IF (debug) THEN
              !write(out_unitp,*) 'BEGINNING GridTOBasis_1D_cplx'
             flush(out_unitp)
            END IF
                !Psi_ggb(i1,:,i3)= 0
              DO i3=1,ubound(Psi_bbb,dim=3)
              DO i1=1,ubound(Psi_bbb,dim=1)

                 Psi_ggb(i1,:,i3) =  matmul( Basis%d0gb  ,Psi_bbb(i1,:,i3))
                  !write(*,*) i1,i3 ,   Psi_ggb(i1,:,i3)

              END DO
              END DO

            IF (debug) THEN
            	!write(out_unitp,*) 'END GridTOBasis_1D_cplx'
            	flush(out_unitp)
            END IF

    END SUBROUTINE  BasisTOGrid_1D_cplx

    SUBROUTINE  GridTOBasis_1D_cplx(Psi_bbb,Psi_ggb,Basis)
           USE UtilLib_m
          TYPE(Basis_t)    , intent(in),target     :: Basis
          complex (kind=Rk), intent(inout)         :: Psi_bbb(:,:,:)
          complex (kind=Rk), intent(in)            :: Psi_ggb(:,:,:)
          logical          , parameter             :: debug = .true.
          integer                                  ::i1,i3


          IF (debug) THEN
           ! write(out_unitp,*) 'BEGINNING  GridTOBasis_1D_cplx'
            flush(out_unitp)
          END IF
           Psi_bbb(:,:,:) = CZERO
            DO i3=1,ubound(Psi_ggb,dim=3)
               DO i1=1,ubound(Psi_ggb,dim=1)

                   Psi_bbb(i1,:,i3) =  matmul( Basis%d0bgw  ,Psi_ggb(i1,:,i3))

               END DO
            END DO

          IF (debug) THEN
          	!write(out_unitp,*) 'END  GridTOBasis_1D_cplx'
          	flush(out_unitp)
          END IF

    END SUBROUTINE  GridTOBasis_1D_cplx


    !++++++++++++++++++++++++++++++++++++++++nD++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE Psi_permut(Psi_2,Psi_1)
        complex(kind=Rk),  intent(inout) ,target       :: Psi_1(:)
        complex(kind=Rk),  intent(inout),target        :: Psi_2(:)
        Psi_1(:) = Psi_2(:)
        Psi_2(:) = CZERO
    END SUBROUTINE Psi_permut
    SUBROUTINE Calc_iqib( Ib1,Ib2,Ib3,Iq1,Iq2,Iq3,Ndim,Basis)
         TYPE(Basis_t),     intent(in),target ::Basis
        integer,intent(inout) , allocatable   :: Ib1(:),Ib2(:),Iq3(:),Iq1(:),Iq2(:),Ib3(:)
        integer,intent(in)                    :: Ndim
         integer                              :: inb
         allocate(Ib3(Ndim-1))
         allocate(Ib2(Ndim-1))
         allocate(Ib1(Ndim-1))

         allocate(Iq3(Ndim-1))
         allocate(Iq2(Ndim-1))
         allocate(Iq1(Ndim-1))

        DO inb = 1,Ndim-1
            if (inb == 1)then
                 Ib3(1) = Product(Basis%tab_basis(2:Ndim-1)%nb)*Basis%tab_basis(Ndim)%nb
                 Iq3(1) =  Product(Basis%tab_basis(2:Ndim-1)%nq)*Basis%tab_basis(Ndim)%nb

                 Iq2(1) =  Basis%tab_basis(1)%nq
                 Ib2(1) =  Basis%tab_basis(1)%nb

                 Iq1(1) = 1
                 Ib1(1) = 1

                elseif(inb == Ndim-1)then

                 Ib3(inb) =  Basis%tab_basis(Ndim)%nb
                 Iq3(inb) =  Basis%tab_basis(Ndim)%nb

                 Ib2(inb) =  Basis%tab_basis(Ndim-1)%nb
                 Iq2(inb) =  Basis%tab_basis(Ndim-1)%nq

                 Iq1(inb) =  Product(Basis%tab_basis(1:Ndim-2)%nq)
                 Ib1(inb) =  Product(Basis%tab_basis(1:Ndim-2)%nb)

                 else
         Ib3(inb) = Product(Basis%tab_basis(inb+1:Ndim-1)%nb)*Basis%tab_basis(Ndim)%nb
         Iq3(inb) = Product(Basis%tab_basis(inb+1:Ndim-1)%nq)*Basis%tab_basis(Ndim)%nb
         Iq2(inb) = Basis%tab_basis(inb)%nq
         Ib2(inb) = Basis%tab_basis(inb)%nb
         Iq1(inb) = Product(Basis%tab_basis(1:inb-1)%nq)
         Ib1(inb) = Product(Basis%tab_basis(1:inb-1)%nb)
                 end if
         END DO
    END SUBROUTINE Calc_iqib

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++






         SUBROUTINE GridTOBasis_nD_cplx(Psi_b,Psi_g,Basis)
        USE UtilLib_m
        TYPE(Basis_t),   intent(in),target                      :: Basis
        complex (kind=Rk), intent(in) ,target                   :: Psi_g(:)
        complex (kind=Rk), intent(inout),target                 :: Psi_b(:)
        complex (kind=Rk), pointer                              :: Psi_ggb(:,:,:)
        complex(kind=Rk) , allocatable  ,target                 :: Psi_x1(:),Psi_x2(:)
        complex (kind=Rk), pointer                              :: Psi_bbb(:,:,:)
        logical,         parameter                              :: debug = .true.
        integer                                                 :: ib,i1,i3,inb,Ndim,iq
        integer , allocatable                                   :: Ib1(:),Ib2(:),Iq3(:),Iq1(:),Iq2(:),Ib3(:)

        IF (debug) THEN
         !write(out_unitp,*) 'BEGINNING GridTOBasis_dnD_cplx'
         !write(*,*) "Psi_g"
         !do iq = 1, size(Psi_g)
         !write(*,*) iq, Psi_g(iq)
         !end do
          flush(out_unitp)
        END IF
           Ndim = size(Basis%tab_basis)
           call Calc_iqib( Ib1,Ib2,Ib3,Iq1,Iq2,Iq3,Ndim,Basis)
        If (Ndim<= 2)then
            Psi_b(:) = CZERO
            Psi_ggb(1:Iq1(1),1:Iq2(1),1:Iq3(1)) => Psi_g
            Psi_bbb(1:Ib1(1),1:Ib2(1),1:Ib3(1)) => Psi_b
            call GridTOBasis_1D_cplx(Psi_bbb,Psi_ggb,Basis%tab_basis(1))

        else
            allocate(Psi_x1(Ib1(1)*Ib2(1)*Iq3(1)))
            Psi_x1(:) = CZERO
            Psi_bbb( 1:Ib1(1),1:Ib2(1),1:Iq3(1))   => Psi_x1
            Psi_ggb(1:Iq1(1),1:Iq2(1),1:Iq3(1))    => Psi_g
            call GridTOBasis_1D_cplx(Psi_bbb,Psi_ggb,Basis%tab_basis(1))
            DO inb = 2,Ndim-2
                allocate(Psi_x2(Ib1(inb)*Ib2(inb)*Iq3(inb)))
                Psi_x2(:) = CZERO
                Psi_bbb( 1:Ib1(inb),1:Ib2(inb),1:Iq3(inb))    => Psi_x2
                Psi_ggb( 1:Ib1(inb),1:Ib2(inb),1:Iq3(inb))    => Psi_x1
                call GridTOBasis_1D_cplx(Psi_bbb,Psi_ggb,Basis%tab_basis(inb))
                deallocate(Psi_x1)
                allocate(Psi_x1(Ib1(inb)*Ib2(inb)*Iq3(inb)))
                call Psi_permut(Psi_x2,Psi_x1)
                deallocate(Psi_x2)
            END DO
            Psi_b(:) = CZERO
            Psi_ggb(1:Iq1(Ndim-1),1:Iq2(Ndim-1),1:Iq3(Ndim-1)) => Psi_x1
            Psi_bbb(1:Ib1(Ndim-1),1:Ib2(Ndim-1),1:Ib3(Ndim-1)) => Psi_b
            call GridTOBasis_1D_cplx(Psi_bbb,Psi_ggb,Basis%tab_basis(Ndim-1))
       End If
       IF (debug) THEN
          !write(out_unitp,*) 'END GridTOBasis_nD_cplx'
          !write(*,*) "Psi_b"
          !do iq = 1, size(Psi_b)
          !  write(*,*) iq, Psi_b(iq)
          !end do
          	flush(out_unitp)
       END IF
       deallocate (Iq1,Iq2,Iq3,Ib1,Ib2,Ib3)
    END SUBROUTINE GridTOBasis_nD_cplx

    SUBROUTINE BasisTOGrid_nD_cplx(Psi_g,Psi_b,Basis)
        USE UtilLib_m
        TYPE(Basis_t),   intent(in),target                   :: Basis
        complex (kind=Rk), intent(in) ,target                :: Psi_b(:)
        complex (kind=Rk), intent(inout),target              :: Psi_g(:)
        complex (kind=Rk), pointer                           :: Psi_bbb(:,:,:)
        complex (kind=Rk) ,allocatable,target                 ::Psi_x1(:),Psi_x2(:)
        complex (kind=Rk), pointer                           :: Psi_ggb(:,:,:)
        logical,         parameter                           :: debug = .true.
        integer                                              :: ib,i1,i3,inb,Ndim,iq
        integer , allocatable                                ::Ib3(:),Iq1(:),Iq2(:),ib1(:),ib2(:),iq3(:)

        IF (debug) THEN
         ! write(out_unitp,*) 'BEGINNING BasisTOGrid_dnD_cplx'
        !write(*,*) "Psi_b"
        !do iq = 1, size(Psi_b)
        !   write(*,*) iq, Psi_b(iq)
        !end do
          flush(out_unitp)
        END IF
       Ndim =size(Basis%tab_basis)
       call Calc_iqib( Ib1,Ib2,Ib3,Iq1,Iq2,Iq3,Ndim,Basis)
        If (Ndim<= 2)then
          Psi_bbb(1:ib1(1),1:ib2(1),1:ib3(1)) => Psi_b
          Psi_ggb(1:iq1(1),1:iq2(1),1:iq3(1)) => Psi_g
          call BasisTOGrid_1D_cplx(Psi_ggb,Psi_bbb,Basis%tab_basis(1))
        else
          allocate(Psi_x1(iq1(1)*iq2(1)*ib3(1)))
          Psi_bbb(1:ib1(1),1:ib2(1),1:ib3(1)) => Psi_b
          Psi_ggb(1:iq1(1),1:iq2(1),1:ib3(1)) => Psi_x1
          call BasisTOGrid_1D_cplx(Psi_ggb,Psi_bbb,Basis%tab_basis(1))
             DO inb = 2,Ndim-2
                allocate(Psi_x2(iq1(inb)*iq2(inb)*ib3(inb)))
                Psi_bbb(1:ib1(Inb),1:ib2(inb),1:ib3(inb)) => Psi_x1
                Psi_ggb(1:iq1(inb),1:iq2(inb),1:ib3(inb)) => Psi_x2
                call BasisTOGrid_1D_cplx(Psi_ggb,Psi_bbb,Basis%tab_basis(inb))
                deallocate(Psi_x1)
                allocate(Psi_x1(iq1(inb)*iq2(inb)*ib3(inb)))
                call Psi_permut(Psi_x2,Psi_x1)
                deallocate(Psi_x2)
             END DO
          Psi_ggb(1:iq1(Ndim-1),1:iq2(Ndim-1),1:iq3(Ndim-1)) => Psi_g
          Psi_bbb(1:ib1(Ndim-1),1:ib2(Ndim-1),1:ib3(Ndim-1)) => Psi_x1
          call BasisTOGrid_1D_cplx(Psi_ggb,Psi_bbb,Basis%tab_basis(Ndim-1))
       End If
       IF (debug) THEN
       	  !write(out_unitp,*) 'END BasisTOGrid_dnD_cplx'
          !  write(*,*) "Psi_g"
          !  do iq = 1, size(Psi_g)
          !  write(*,*) iq, Psi_g(iq)
          ! end do
       	  flush(out_unitp)
       END IF
       deallocate (Iq1,Iq2,Iq3,Ib1,Ib2,Ib3)
    END SUBROUTINE BasisTOGrid_nD_cplx

    SUBROUTINE Calc_Q_grid(Q,Basis)

       implicit none
       TYPE (Basis_t)  ,intent(in)                         :: Basis
       integer ,ALLOCATABLE                                :: Tab_iq(:),NDend(:)
       integer                                             :: inb,ndim,iq
       real(Kind = Rk), intent(inout),ALLOCATABLE          ::Q(:,:)
       TYPE (NDindex_t)                                    :: NDindex
       logical                                             ::Endloop
       ndim = SIZE(Basis%tab_basis)-1
       allocate(Tab_iq(Ndim))
       allocate(NDend(Ndim))
       allocate(Q(Basis%nq,Ndim))

       do inb= 1, Ndim
           NDend(inb) = Basis%NDindexq%NDend(inb)
       end do
      CALL Init_NDindex(NDindex,NDend,Ndim)
      Call Init_tab_ind(Tab_iq,NDindex)
      Iq=0
      DO
        Iq=Iq+1
        CALL increase_NDindex(Tab_iq,NDindex,Endloop)
         IF (Endloop) exit
          do inb = 1,Ndim
            Q(iq,inb) = Basis%tab_basis(inb)%X(Tab_iq(inb))
          end do
      END DO
    END SUBROUTINE Calc_Q_grid

END MODULE Basis_m
