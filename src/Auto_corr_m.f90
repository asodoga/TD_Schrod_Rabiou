module Auto_corr_m
   USE UtilLib_m
   USE Ana_psi_m
   USE psi_m
   USE Basis_m

   implicit none
   public :: Calc_Auto_corr, Calc_fft_Auto_corr, Test_calc_S, Calc_Num_S

contains

   SUBROUTINE Calc_Auto_corr(psi0, psi_dt, corre_coeff, arg_corre_coeff, propa_name)
      USE UtilLib_m
      USE psi_m

      TYPE(psi_t), intent(in)                 :: psi0, psi_dt
      complex(kind=Rk), intent(inout)         :: corre_coeff
      real(kind=Rk), intent(inout)            :: arg_corre_coeff
      character(*), intent(in)                :: propa_name

      !local variables---------------------------------------------------
      TYPE(psi_t)                             :: psi
      real(kind=Rk)                           :: X, Y

      write (out_unitp, *) 'Beging Calc_Auto_corr'

      if (propa_name == 'hagedorn') then

         call init_psi(psi, psi0%Basis, cplx=.TRUE., grid=.false.)
         psi%CVec = CZERO
         call Hagedorn0(psi, psi_dt)
         corre_coeff = dot_product(Psi0%CVec, psi%CVec)/dot_product(Psi%CVec, psi%CVec)
         X = real(corre_coeff, kind=RK)
         Y = aimag(corre_coeff)
         arg_corre_coeff = atan2(Y, X)
         call dealloc_psi(psi)

      else
         corre_coeff = dot_product(psi0%CVec, psi_dt%CVec)/dot_product(psi_dt%CVec, psi_dt%CVec)
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

   SUBROUTINE Hagedorn1(psi_dt_2, psi_dt_1)
      TYPE(psi_t), intent(in), target                  :: psi_dt_1
      TYPE(psi_t), intent(inout), target               :: psi_dt_2
      complex(kind=Rk), pointer                        :: BBB1(:, :, :), BBB2(:, :, :)
      complex(kind=Rk), allocatable, target            :: B1(:), B2(:)
      !logical, parameter                              :: debug = .true.
      integer                                          :: inb, i1, i3, Ndim, iq, jq
      Integer, allocatable                             :: Ib1(:), Ib2(:), Ib3(:)
      complex(Kind=Rk), allocatable                    :: S(:, :)
      real(Kind=Rk)                                    :: p1, p2, s1, s2, x1, x2

      Call Calc_index(Ib1=Ib1, Ib2=Ib2, Ib3=Ib3, Basis=psi_dt_1%Basis)
      Ndim = size(psi_dt_1%Basis%tab_basis) - 1
      write (out_unitp, *) 'Begin Hagedorn projection'

      If (Ndim == 1) then
         BBB1(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_1%CVec
         BBB2(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_2%CVec

         !constuction of ovelap  S(:,:)
         !----------------------------------------------------------------------------
         allocate (S(psi_dt_1%Basis%tab_basis(1)%nb, psi_dt_1%Basis%tab_basis(1)%nb))
         s1 = psi_dt_1%Basis%tab_basis(1)%scaleQ
         s2 = psi_dt_2%Basis%tab_basis(1)%scaleQ
         p1 = psi_dt_1%Basis%tab_basis(1)%Imp_k
         p2 = psi_dt_2%Basis%tab_basis(1)%Imp_k
         x1 = psi_dt_1%Basis%tab_basis(1)%Q0
         x2 = psi_dt_2%Basis%tab_basis(1)%Q0

         Do iq = 1, psi_dt_1%Basis%tab_basis(1)%nb
            Do jq = 1, psi_dt_1%Basis%tab_basis(1)%nb

               CALL Hermite_product_integral(S(iq, jq), psi_dt_1%Basis%tab_basis(1)%X, psi_dt_1%Basis%tab_basis(1)%W, iq, jq, x1,&
               & x2, s1, s2, p1, p2)

            End Do
         End Do
         !----------------------------------------------------------------------------------

         psi_dt_2%CVec(:) = CZERO
         DO i3 = 1, ubound(BBB1, dim=3)
         DO i1 = 1, ubound(BBB1, dim=1)
            BBB2(i1, :, i3) = matmul(BBB1(i1, :, i3), S)
         END DO
         END DO

         deallocate (S)
      else

         psi_dt_2%CVec = CZERO
         allocate (B1(Ib1(1)*Ib2(1)*Ib3(1)))
         BBB1(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_1%CVec
         BBB2(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => B1

         !constuction of ovelap  S(:,:)
         !----------------------------------------------------------------------------
         allocate (S(psi_dt_1%Basis%tab_basis(1)%nb, psi_dt_1%Basis%tab_basis(1)%nb))
         s1 = psi_dt_1%Basis%tab_basis(1)%scaleQ
         s2 = psi_dt_2%Basis%tab_basis(1)%scaleQ
         p1 = psi_dt_1%Basis%tab_basis(1)%Imp_k
         p2 = psi_dt_2%Basis%tab_basis(1)%Imp_k
         x1 = psi_dt_1%Basis%tab_basis(1)%Q0
         x2 = psi_dt_2%Basis%tab_basis(1)%Q0

         Do iq = 1, psi_dt_1%Basis%tab_basis(1)%nb
            Do jq = 1, psi_dt_1%Basis%tab_basis(1)%nb
               CALL Hermite_product_integral(S(iq, jq), psi_dt_1%Basis%tab_basis(1)%X, psi_dt_1%Basis%tab_basis(1)%W, iq, jq, x1,&
               & x2, s1, s2, p1, p2)

            End Do

         End Do
         !----------------------------------------------------------------------------------

         DO i3 = 1, ubound(BBB1, dim=3)
         DO i1 = 1, ubound(BBB1, dim=1)

            BBB2(i1, :, i3) = matmul(BBB1(i1, :, i3), S)

         END DO
         END DO
         deallocate (S)
         Do inb = 2, Ndim

            allocate (B2(Ib1(inb)*Ib2(inb)*Ib3(inb)))
            BBB1(1:Ib1(inb), 1:Ib2(inb), 1:Ib3(inb)) => B1
            BBB2(1:Ib1(inb), 1:Ib2(inb), 1:Ib3(inb)) => B2

            !constuction of ovelap  S(:,:)
            !----------------------------------------------------------------------------
            allocate (S(psi_dt_1%Basis%tab_basis(Inb)%nb, psi_dt_1%Basis%tab_basis(Inb)%nb))
            s1 = psi_dt_1%Basis%tab_basis(Inb)%scaleQ
            s2 = psi_dt_2%Basis%tab_basis(Inb)%scaleQ
            p1 = psi_dt_1%Basis%tab_basis(Inb)%Imp_k
            p2 = psi_dt_2%Basis%tab_basis(Inb)%Imp_k
            x1 = psi_dt_1%Basis%tab_basis(Inb)%Q0
            x2 = psi_dt_2%Basis%tab_basis(Inb)%Q0

            Do iq = 1, psi_dt_1%Basis%tab_basis(Inb)%nb
               Do jq = 1, psi_dt_1%Basis%tab_basis(Inb)%nb
             CALL Hermite_product_integral(S(iq, jq), psi_dt_1%Basis%tab_basis(Inb)%X, psi_dt_1%Basis%tab_basis(Inb)%W, iq, jq, x1,&
                                                                                                               & x2, s1, s2, p1, p2)
               End Do
            End Do
            !----------------------------------------------------------------------------------

            DO i3 = 1, ubound(BBB1, dim=3)
            DO i1 = 1, ubound(BBB1, dim=1)
               BBB2(i1, :, i3) = matmul(BBB1(i1, :, i3), S)
            END DO
            END DO
            B1 = B2
            B2 = CZERO
            if (inb == Ndim) psi_dt_2%CVec = B1
            deallocate (S)
         END DO

      END IF
      write (out_unitp, *) 'END Hagedorn projection'

   END SUBROUTINE

   SUBROUTINE Hagedorn0(psi_dt_2, psi_dt_1)
      TYPE(psi_t), intent(in), target                  :: psi_dt_1
      TYPE(psi_t), intent(inout), target               :: psi_dt_2
      complex(kind=Rk), pointer                        :: BBB1(:, :, :), BBB2(:, :, :)
      complex(kind=Rk), allocatable, target            :: B1(:), B2(:)
      !logical, parameter                              :: debug = .true.
      integer                                          :: inb, i1, i3, Ndim, iq, jq
      Integer, allocatable                             :: Ib1(:), Ib2(:), Ib3(:), nb, nq
      real(Kind=Rk), allocatable                       :: x(:), w(:)
      complex(Kind=Rk), allocatable                    :: S(:, :)
      real(Kind=Rk)                                    :: s1, s2, s3, x1, x2, x3, p1, p2

      Call Calc_index(Ib1=Ib1, Ib2=Ib2, Ib3=Ib3, Basis=psi_dt_1%Basis)
      Ndim = size(psi_dt_1%Basis%tab_basis) - 1
      write (out_unitp, *) 'Begin Hagedorn projection'

      If (Ndim == 1) then
         BBB1(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_1%CVec
         BBB2(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_2%CVec
         !constuction of ovelap  S(:,:)
         !----------------------------------------------------------------------------
         allocate (S(psi_dt_1%Basis%tab_basis(1)%nb, psi_dt_1%Basis%tab_basis(1)%nb))
         allocate (x(psi_dt_1%Basis%tab_basis(1)%nq))
         allocate (w(psi_dt_1%Basis%tab_basis(1)%nq))
         call hercom(psi_dt_1%Basis%tab_basis(1)%nq, x(:), w(:))
         x1 = psi_dt_1%Basis%tab_basis(1)%Q0
         x2 = psi_dt_2%Basis%tab_basis(1)%Q0
         s1 = psi_dt_1%Basis%tab_basis(1)%scaleQ
         s2 = psi_dt_2%Basis%tab_basis(1)%scaleQ
         p1 = psi_dt_1%Basis%tab_basis(1)%Imp_k
         p2 = psi_dt_2%Basis%tab_basis(1)%Imp_k
         nb = psi_dt_1%Basis%tab_basis(1)%nb
         nq = psi_dt_1%Basis%tab_basis(1)%nq

         call Calc_Num_S(S=S, nb=nb, nq=nq, x1=x1, x2=x2, s1=s1, s2=s1, p1=p1, p2=p2)
         !----------------------------------------------------------------------------------
         psi_dt_2%CVec(:) = CZERO
         call project_1D(BBB2, BBB1, S)
         deallocate (S, x, w)
      else
         psi_dt_2%CVec = CZERO
         allocate (B1(Ib1(1)*Ib2(1)*Ib3(1)))
         BBB1(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_1%CVec
         BBB2(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => B1
         !constuction of ovelap  S(:,:)
         !----------------------------------------------------------------------------
         allocate (S(psi_dt_1%Basis%tab_basis(1)%nb, psi_dt_1%Basis%tab_basis(1)%nb))
         allocate (x(psi_dt_1%Basis%tab_basis(1)%nq))
         allocate (w(psi_dt_1%Basis%tab_basis(1)%nq))
         call hercom(psi_dt_1%Basis%tab_basis(1)%nq, x(:), w(:))
         x1 = psi_dt_1%Basis%tab_basis(1)%Q0
         x2 = psi_dt_2%Basis%tab_basis(1)%Q0
         s1 = psi_dt_1%Basis%tab_basis(1)%scaleQ
         s2 = psi_dt_2%Basis%tab_basis(1)%scaleQ
         p1 = psi_dt_1%Basis%tab_basis(1)%Imp_k
         p2 = psi_dt_2%Basis%tab_basis(1)%Imp_k

         if (psi_dt_1%Basis%tab_basis(1)%Basis_name == 'herm' .or. psi_dt_1%Basis%tab_basis(1)%Basis_name == 'ho') then

            call Calc_Num_S(S=S, nb=nb, nq=nq, x1=x1, x2=x2, s1=s1, s2=s1, p1=p1, p2=p2)

         else
            S(:, :) = ZERO
            Do iq = 1, psi_dt_1%Basis%tab_basis(1)%nb
               S(iq, iq) = ONE
            End Do
         end if
         !----------------------------------------------------------------------------------
         call project_1D(BBB2, BBB1, S)
         deallocate (S, x, w)
         Do inb = 2, Ndim
            allocate (B2(Ib1(inb)*Ib2(inb)*Ib3(inb)))
            BBB1(1:Ib1(inb), 1:Ib2(inb), 1:Ib3(inb)) => B1
            BBB2(1:Ib1(inb), 1:Ib2(inb), 1:Ib3(inb)) => B2
            !constuction of ovelap  S(:,:)
            !----------------------------------------------------------------------------
            allocate (S(psi_dt_1%Basis%tab_basis(Inb)%nb, psi_dt_1%Basis%tab_basis(Inb)%nb))
            allocate (x(psi_dt_1%Basis%tab_basis(Inb)%nq))
            allocate (w(psi_dt_1%Basis%tab_basis(Inb)%nq))
            call hercom(psi_dt_1%Basis%tab_basis(Inb)%nq, x(:), w(:))
            x1 = psi_dt_1%Basis%tab_basis(Inb)%Q0
            x2 = psi_dt_2%Basis%tab_basis(Inb)%Q0
            s1 = psi_dt_1%Basis%tab_basis(Inb)%scaleQ
            s2 = psi_dt_2%Basis%tab_basis(Inb)%scaleQ
            p1 = psi_dt_1%Basis%tab_basis(Inb)%Imp_k
            p2 = psi_dt_2%Basis%tab_basis(Inb)%Imp_k

            if (psi_dt_1%Basis%tab_basis(Inb)%Basis_name == 'herm' .or. psi_dt_1%Basis%tab_basis(Inb)%Basis_name == 'ho') then

               Do iq = 1, psi_dt_1%Basis%tab_basis(Inb)%nb
                  call Calc_Num_S(S=S, nb=nb, nq=nq, x1=x1, x2=x2, s1=s1, s2=s1, p1=p1, p2=p2)
               End Do

            else
               S(:, :) = ZERO
               Do iq = 1, psi_dt_1%Basis%tab_basis(Inb)%nb
                  S(iq, iq) = ONE
               End Do
            end if
            ! ---------------------------------------------------------------------------------
            call project_1D(BBB2, BBB1, S)
            B1 = B2
            B2 = CZERO
            deallocate (B2)
            if (inb == Ndim) psi_dt_2%CVec = B1
            deallocate (S, x, w)

         END DO

      END IF
      write (out_unitp, *) 'END Hagedorn projection'
   END SUBROUTINE

   SUBROUTINE project_1D(BBB2, BBB1, S)
      complex(kind=Rk), intent(in)         :: BBB1(:, :, :)
      complex(kind=Rk), intent(inout)      ::  BBB2(:, :, :)
      integer                              ::  i1, i3
      complex(Kind=Rk), intent(in)         :: S(:, :)

      DO i3 = 1, ubound(BBB1, dim=3)
         DO i1 = 1, ubound(BBB1, dim=1)
            BBB2(i1, :, i3) = matmul(BBB1(i1, :, i3), S)
         END DO
      END DO

   END SUBROUTINE

   SUBROUTINE Calc_Num_S(S, nb, nq, x1, x2, s1, s2, p1, p2)
      real(Kind=Rk)                    :: s1, s2, x1, x2, p1, p2
      integer, intent(in)              :: nb, nq
      complex(Kind=Rk), intent(inout)  :: S(:, :)
      real(Kind=Rk), allocatable       :: x(:), w(:)
      real(Kind=Rk)                    ::  s3, x3
      integer                          :: iq, jq
      allocate (x(nq))
      allocate (w(nq))
      call hercom(nq, x(:), w(:))
      s3 = sqrt(s1*s1 + s2*s2)/sqrt(TWO)
      x3 = (s1*s1*x1 + s2*s2*x2)/(s1*s1 + s2*s2)
      w(:) = w(:)/s3
      x(:) = x3 + x(:)/s3
      s(:, :) = ZERO

      Do iq = 1, nb
         Do jq = 1, nb
            CALL Hermite_product_integral(S(iq, jq), x, w, iq, jq, x1, x2, s1, s2, p1, p2)
         End Do
      End Do
      !print *, ' Beging print Overlap matrix'
      !Do jq = 1, nb
      !   write (*, *) S(jq, :)
      !End Do
      !print *, ' End print Overlap matrix'
   END SUBROUTINE

   SUBROUTINE Test_calc_S(S, nb, nq, x1, x2, s1, s2, p1, p2)
      real(Kind=Rk)            :: s1, s2, x1, x2, p1, p2
      integer, intent(in)      :: nb, nq
      complex(Kind=Rk), intent(inout) :: S(:, :)
      real(Kind=Rk), allocatable     :: x(:), w(:)
      real(Kind=Rk)                  ::  s3, x3
      integer                        :: iq, jq

      allocate (x(nq))
      allocate (w(nq))
      call hercom(nq, x(:), w(:))
      s3 = sqrt(s1*s1 + s2*s2)/sqrt(TWO)
      x3 = (s1*s1*x1 + s2*s2*x2)/(s1*s1 + s2*s2)
      w(:) = w(:)/s3
      x(:) = x3 + x(:)/s3
      s(:, :) = ZERO

      Do iq = 1, nb
         Do jq = 1, nb
            CALL Hermite_product_integral(S(iq, jq), x, w, iq, jq, x1, x2, s1, s2, p1, p2)
         End Do
      End Do
      print *, ' Beging print Overlap matrix'
      Do jq = 1, nb
         write (*, *) S(jq, :)
      End Do
      print *, ' End print Overlap matrix'
   END SUBROUTINE

end module Auto_corr_m

