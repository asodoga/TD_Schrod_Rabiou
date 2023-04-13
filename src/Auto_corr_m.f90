module Auto_corr_m
   USE UtilLib_m
   USE Ana_psi_m
   USE psi_m

   implicit none
   public :: Calc_Auto_corr, Calc_fft_Auto_corr

contains

   SUBROUTINE Calc_Auto_corr(psi0, psi_dt, corre_coeff, arg_corre_coeff, propa_name)
      USE UtilLib_m
      USE psi_m

      TYPE(psi_t), intent(in)                 :: psi0, psi_dt
      complex(kind=Rk), intent(inout)         :: corre_coeff
      real(kind=Rk), intent(inout)            :: arg_corre_coeff
      !real(kind=Rk)     ,intent(in)          :: T
      character(*), intent(in)                :: propa_name

      !local variables---------------------------------------------------
      TYPE(psi_t)                             :: psi
      real(kind=Rk)                           :: X, Y
      real(kind=Rk), allocatable              :: Qt(:), SQt(:)
      integer                                 :: Ndim

      write (out_unitp, *) 'Beging Calc_Auto_corr'

      if (propa_name == 'hagedorn') then

         call init_psi(psi, psi0%Basis, cplx=.TRUE., grid=.false.)
         psi%CVec = CZERO
         Ndim = size(psi0%Basis%tab_basis) - 1
         allocate (Qt(Ndim), SQt(Ndim))
         Qt(:) = ZERO; SQt(:) = ONE

         call Calc_AVQ_nD(psi0=psi0, AVQ=Qt, SQ=SQt)
         call construct_primitive_basis(psi_dt%Basis, Qt, SQt)
         call projection(psi, psi_dt)

         corre_coeff = dot_product(psi0%CVec, psi%CVec)
         X = real(corre_coeff, kind=RK)
         Y = aimag(corre_coeff)
         arg_corre_coeff = atan2(Y, X)
         call dealloc_psi(psi)

      else
         corre_coeff = dot_product(psi0%CVec, psi_dt%CVec)
         X = real(corre_coeff, kind=RK)
         Y = aimag(corre_coeff)
         arg_corre_coeff = atan2(Y, X)
      end if
      write (out_unitp, *) 'corre_coeff =', corre_coeff, 'arg_corre_coeff=', arg_corre_coeff

      write (out_unitp, *) 'End Calc_Auto_corr'

   End SUBROUTINE

   SUBROUTINE Calc_fft_Auto_corr(autocor_function, time, fft_autocor_function, delta_t, N)
      USE NumParameters_m
      USE UtilLib_m
      COMPLEX(KIND=Rk), INTENT(IN), allocatable, DIMENSION(:)        :: autocor_function(:)
      COMPLEX(KIND=Rk), INTENT(INOUT), ALLOCATABLE                  :: fft_autocor_function(:)
      REAL(KIND=Rk), INTENT(IN), DIMENSION(:)                       :: time
      REAL(KIND=Rk), ALLOCATABLE, DIMENSION(:)                       :: w
      INTEGER, intent(in)                                           :: N
      REAL(KIND=Rk)                                                :: delta_t
      REAL(KIND=Rk)                                                 :: wm, wmax, dw
      INTEGER                                                       :: Iw, nw, I
      OPEN (UNIT=100, FILE="fft_autocor_function.dat")
      fft_autocor_function(:) = CZERO
      wm = -0.001_Rk; wmax = 5._Rk; dw = 0.005_Rk
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

end module Auto_corr_m

