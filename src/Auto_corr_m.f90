module Auto_corr_m
   USE QDUtil_m
   USE psi_m
   USE Basis_m
   Use Hagedorn_m

   implicit none
   public :: Calc_Auto_corr, Calc_fft_Auto_corr

contains

   SUBROUTINE Calc_Auto_corr(psi0, psi_dt, corre_coeff, arg_corre_coeff, propa_name,t,it)
      USE QDUtil_m
      USE psi_m

      TYPE(psi_t), intent(in)                    :: psi0, psi_dt
      real(kind=Rkind),intent(in),optional       :: t
      integer,intent(in),optional                :: it
      complex(kind=Rkind), intent(inout)         :: corre_coeff
      real(kind=Rkind), intent(inout)            :: arg_corre_coeff
      character(*), intent(in)                   :: propa_name

      !local variables---------------------------------------------------
      TYPE(psi_t)                                :: psi,psi_t0
      TYPE(Basis_t)                              ::  Basis0,Basis_dt
      real(kind=Rkind)                           :: X, Y
      integer                                    :: ib

      write (out_unit, *) 'Beging Calc_Auto_corr'

      if (propa_name == 'hagedorn') then

         call init_Basis1_TO_Basis2(Basis0, psi0%Basis)
         call construct_primitive_basis(Basis0)
         call init_Basis1_TO_Basis2(Basis_dt, psi_dt%Basis)
         call construct_primitive_basis(Basis_dt)
         call init_psi(psi, Basis0, cplx=.TRUE., grid=.false.)
         call init_psi(psi_t0, Basis_dt, cplx=.TRUE., grid=.false.)

         psi%CVec = CZERO
         psi_t0%CVec(:) =psi_dt%CVec(:)
         call  Hagedorn_Inv2(psi, psi_t0)
         corre_coeff = dot_product(psi0%CVec, psi%CVec)/dot_product(psi%CVec, psi%CVec)**2
         X = real(corre_coeff, kind=Rkind)
         Y = aimag(corre_coeff)
         arg_corre_coeff = atan2(Y, X)
         if(present(t)) then
           call write_psi(psi=psi, psi_cplx=.false., print_psi_grid=.false. &
            , print_basis=.false., t=t, int_print=27, real_part=.false.)
             write(27,*)
         end if

          if(present(it)) call test_write(psi,ib=27)    
           if(present(it)) write(27,*)     

         call dealloc_psi(psi)
         call dealloc_psi(psi_t0)

      else

         corre_coeff = dot_product(psi0%CVec, psi_dt%CVec)/dot_product(psi_dt%CVec, psi_dt%CVec)**2
         X = real(corre_coeff, kind=Rkind)
         Y = aimag(corre_coeff)
         arg_corre_coeff = atan2(Y, X)
         if(present(t)) then
           call write_psi(psi=psi_dt, psi_cplx=.false., print_psi_grid=.false. &
           , print_basis=.false., t=t, int_print=27, real_part=.false.)
            write(27,*)
         end if
         if(present(it)) call test_write(psi_dt,ib=27)
         if(present(it)) write(27,*)
      end if
      !write (out_unit, *) 'corre_coeff =', corre_coeff, 'arg_corre_coeff=', arg_corre_coeff

      write (out_unit, *) 'End Calc_Auto_corr'

   End SUBROUTINE

   SUBROUTINE Calc_fft_Auto_corr(autocor_function, time, fft_autocor_function, delta_t, N)
      USE QDUtil_m
      COMPLEX(KIND=Rkind), INTENT(IN), allocatable, DIMENSION(:)        :: autocor_function(:)
      COMPLEX(KIND=Rkind), INTENT(INOUT), ALLOCATABLE                   :: fft_autocor_function(:)
      REAL(KIND=Rkind), INTENT(IN), DIMENSION(:)                        :: time
      REAL(KIND=Rkind), ALLOCATABLE, DIMENSION(:)                       :: w
      INTEGER, intent(in)                                               :: N
      REAL(KIND=Rkind)                                                  :: delta_t
      REAL(KIND=Rkind)                                                  :: wm, wmax, dw
      INTEGER                                                           :: Iw, nw, I
      OPEN (UNIT=100, FILE="fft_autocor_function.dat")
      fft_autocor_function(:) = CZERO
      wm = -0.001_Rkind; wmax = 5._Rkind; dw = 0.005_Rkind
      nw = int((wmax - wm)/dw)
      !allocate( fft_autocor_function(0:nw))
      ALLOCATE (w(0:nw - 1))
      DO Iw = 0, nw - 1
         w(Iw) = wm + float(Iw)*dw
         fft_autocor_function(Iw) = ZERO
         DO I = 1, N
            fft_autocor_function(Iw) = fft_autocor_function(Iw) + autocor_function(Iw)*EXP(EYE*w(Iw)*time(I))*delta_t
         END DO !itime
         WRITE (100, *) w(Iw), ABS(fft_autocor_function(Iw))
      END DO !omega
   END SUBROUTINE



   SUBROUTINE test_write(psi,ib)

   implicit none

   TYPE(psi_t), intent(in)                    :: psi
   integer,intent(in)                          :: ib
  integer                                     :: i

  !do i =1,size(psi%CVec)

  write(ib,*) psi%CVec(:)
 ! end do

   END SUBROUTINE
end module Auto_corr_m

