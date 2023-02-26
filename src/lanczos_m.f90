!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!
!> MODULE: lanczos_m
!
!> This module is used to do short iterative lanczos evolutions on a 
!> symetric tridiagonal matrix.
!----------------------------------------------------------------------------------
module lanczos_m
  use NumParameters_m
  use Basis_m
  use psi_m
  use Op_m
  Use diago_m
  implicit none
private
public::find_alpha,fill_lan_Matrix,construct_Diag_Offdiag  , Lan_eign_syst_solve
contains
  subroutine find_alpha(psi, alpha,Basis)
    implicit none

    TYPE(psi_t), intent(in)                            :: psi
    TYPE (Basis_t), intent(in),  target                :: Basis
    real(kind=Rk),intent(inout)                        :: alpha
    TYPE(psi_t)                                        :: HPsi
    TYPE(Op_t)                                         :: Op

                  !!! alpha(k) = <q(:,k)|H|q(:,k)>!!!!
    CALL init_psi(Hpsi,  Basis,    cplx=.TRUE.   ,grid =.false.)
    call calc_OpPsi(Op,psi,Hpsi)
    alpha = real(dot_product(psi%CVec,Hpsi%CVec),kind=Rk)

 end subroutine find_alpha
  SUBROUTINE fill_lan_Matrix(psi0,tdm,kmax,Alpha,Beta)
    implicit none
    real(kind=Rk),intent(inout),allocatable       :: tdm(:,:)
    type(psi_t),intent(in)                        :: psi0
    integer,intent(in)                            ::kmax
    integer                                       ::i,j
    real(kind=Rk),allocatable,intent(inout)       :: Alpha(:),Beta(:)
    !****************************************************************************************
    call construct_Diag_Offdiag(psi0 ,Alpha,Beta,kmax)
    allocate(tdm(kmax,kmax))
    tdm(:,:) = ZERO
    do i = 1,kmax,1
         tdm(i,i)= Alpha(i)
         tdm(i,i+1)= Beta(i+1)
         tdm(i,i-1)= Beta(i)

    end do
   ! do i = 1,kmax
       ! write(*,*) Alpha(i),Beta(i)
    !end do
   ! print*,"*********************************************************"
    !do i = 1,kmax
     !   write(*,*) tdm(:,i)
   ! end do

  END SUBROUTINE fill_lan_Matrix

  SUBROUTINE construct_Diag_Offdiag(psi0,Alpha,Beta,kmax)
    implicit none
    type(psi_t),intent(in)                       :: psi0
    TYPE(Op_t)                                   :: H
    real(kind=Rk),intent(inout),allocatable      :: Alpha(:),Beta(:)
    real(kind=Rk)                                :: N,error_threshold
    complex(kind=Rk),allocatable                 ::  q(:,:),rkk(:,:)
    TYPE(psi_t)                                 :: qk ,Hq,Hqk
    integer,intent(in)                           ::kmax
    integer                                      ::k,Ndim

    !*******************************************************************
    Ndim = size(psi0%Basis%tab_basis)
     allocate(q(psi0%Basis%nq*psi0%Basis%tab_basis(Ndim)%nb,0:kmax))
     allocate(rkk(psi0%Basis%nq*psi0%Basis%tab_basis(Ndim)%nb,0:kmax))
     allocate(Alpha(0:kmax))
     allocate(Beta(0:kmax))
    CALL init_psi(qk,   psi0%Basis,    cplx=.TRUE.   ,grid =.false.)
    CALL init_psi(Hq,   psi0%Basis,    cplx=.TRUE.   ,grid =.false.)
    CALL init_psi(Hqk,  psi0%Basis,    cplx=.TRUE.   ,grid =.false.)
     Beta(0) =ZERO
     Beta(1) =ZERO
     Alpha(:)  =  ZERO
     rkk(:,:) = CZERO
     q(:,0)   = CZERO
     error_threshold = 0.0000000000000001_Rk
     call Calc_Norm_OF_Psi(psi0,N)
     q(:,1) = psi0%CVec(:)/N
     k =1
    do while (k<kmax)
        !****************<qk|H|qk>******************
        qk%CVec(:) =q(:,k)
        call find_alpha(qk,Alpha(k),psi0%Basis)
        call calc_OpPsi(H,qk,Hqk)
        rkk(:,k) = Hqk%CVec(:)-Alpha(k)*qk%CVec(:)-Beta(k-1)*q(:,k-1)
        Beta(k+1) = sqrt(real(dot_product(rkk(:,k),rkk(:,k)),kind=Rk))
        q(:,k+1)=  rkk(:,k)/Beta(k+1)
        k = k+1

    end do

  END SUBROUTINE construct_Diag_Offdiag



    SUBROUTINE Lan_eign_syst_solve(psi0,psi_dt,kmax,delta_t)
        Use diago_m
        Use psi_m
        type(psi_t),intent(in)                     ::psi0
        type(psi_t),intent(inout)                  ::psi_dt
        real(kind=Rk) ,allocatable                 ::Matrix(:,:)
        complex(kind=Rk) ,allocatable              :: ak(:)
        real(kind=Rk) ,intent(in)                  ::delta_t
         integer,intent(in)                        :: kmax
        !logical,          parameter               :: debug = .true.
        logical,         parameter                 ::debug = .false.
        integer                                    :: ib,jb
        real (kind=Rk), allocatable                :: EigenVal(:),EigenVec(:,:)
        real(kind=Rk),allocatable                  :: Alpha(:),Beta(:)

        IF (debug) THEN
            write(out_unitp,*) 'BEGINNING Eig_syst_solve'
            flush(out_unitp)
        END IF
        psi_dt%CVec(:)=CZERO
        call fill_lan_Matrix(psi0,Matrix,kmax ,Alpha,Beta)
        allocate(EigenVal(kmax ))
        allocate(EigenVec(kmax ,kmax ))
        allocate(ak(kmax ))
        CALL  diagonalization(matrix,EigenVal,EigenVec,kmax )
        DO ib=1,kmax
            ak(ib) = real(dot_product(EigenVec(:,ib),psi0%CVec),kind=Rk)
            ak(ib) = ak(ib)*exp(-EYE*EigenVal(ib)*delta_t)
            psi_dt%CVec(:) =   psi_dt%CVec(:)+ ak(ib)*EigenVec(:,ib) 
        END DO
       ! call Write_psi(psi_dt)
     !  Write(out_unitp,*)
     !  Write(out_unitp,*)
     !  Write(out_unitp,*) 'eigenvalues = '
     !  DO ib=1,kmax
     !      write(out_unitp,*) EigenVal(ib)
     !  END DO
     !  Write(out_unitp,*)
     !  Write(out_unitp,*)
     !  DO ib=1,kmax
     !      write(*,*) (EigenVec(ib,jb),jb=1,kmax)
     !  END DO

        IF (debug) THEN
            write(out_unitp,*) 'END Eig_syst_solve'
            flush(out_unitp)
        END IF
    END SUBROUTINE Lan_eign_syst_solve
    end module lanczos_m
