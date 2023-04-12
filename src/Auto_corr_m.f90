module Auto_corr_m
   USE UtilLib_m
   USE psi_m

   implicit none
   public :: Calc_Auto_corr, Calc_fft_Auto_corr

contains

   SUBROUTINE Calc_Auto_corr(psi0, psi_dt, corre_coeff, arg_corre_coeff, propa_name)
      USE UtilLib_m
      USE psi_m

      TYPE(psi_t), intent(inout)         :: psi0, psi_dt
      complex(kind=Rk), intent(inout)         :: corre_coeff
      real(kind=Rk), intent(inout)         :: arg_corre_coeff
      !real(kind=Rk)     ,intent(in)            ::T
      character(*), intent(in)            :: propa_name

      !local variables---------------------------------------------------
      TYPE(psi_t)                             :: psi, psi1
      real(kind=Rk)                            :: X, Y

      CALL init_psi(psi, psi0%Basis, cplx=.TRUE., grid=.false.)
      CALL init_psi(psi1, psi0%Basis, cplx=.TRUE., grid=.false.)

      IF (psi0%Grid) then
         call GridTOBasis_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
      ELSE
         psi%CVec(:) = psi0%CVec(:)
      END IF

      if (propa_name == 'hagedorn') then
         psi1%CVec = CZERO
         psi_dt%Basis => psi0%Basis
         call projection(psi1, psi_dt)
         corre_coeff = dot_product(psi1%CVec, psi_dt%CVec)
         X = real(corre_coeff, kind=RK)
         Y = aimag(corre_coeff)
         arg_corre_coeff = atan2(Y, X)

      else
         corre_coeff = dot_product(psi%CVec, psi_dt%CVec)
         X = real(corre_coeff, kind=RK)
         Y = aimag(corre_coeff)
         arg_corre_coeff = atan2(Y, X)
      end if
      write (out_unitp, *) 'corre_coeff =', corre_coeff, 'arg_corre_coeff=', arg_corre_coeff
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
