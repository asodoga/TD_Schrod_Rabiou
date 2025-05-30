module psi_m

   USE param_WP0_m
   USE Basis_m
   USE NDindex_m
   implicit none

   TYPE :: psi_t
   
      type(Basis_t), pointer      :: Basis
      real(kind= Rkind), allocatable :: RVec(:)
      complex(kind= Rkind), allocatable :: CVec(:)
      logical                        :: Grid = .true.

   CONTAINS

      PRIVATE
      PROCEDURE, PASS         :: Copy_psi    ! Copy content from other psi1 to psi2 instance,
      GENERIC, PUBLIC         :: ASSIGNMENT(=) => Copy_psi

   END TYPE psi_t

   public :: psi_t, write_psi, init_psi, dealloc_psi, write_psi_grid
   public :: write_psi_basis, Calc_Norm_OF_PsiBasis, Calc_Norm_OF_PsiGrid, Calc_Norm_OF_Psi
   public:: Projection,Projection1,psi_init_GWP,Ecrire_psi
contains


   SUBROUTINE copy_psi(psi_out, psi_in)

      CLASS(psi_t), intent(in)     :: psi_in
      CLASS(psi_t), intent(inout)  :: psi_out

      IF (allocated(psi_in%RVec)) THEN
         !write(out_unit,*) 'Coping psi_in in psi_out (real):'
         psi_out%RVec(:) = psi_in%RVec(:)
         !CALL init_Basis1_TO_Basis2 (psi_in%Basis,psi_out%Basis)
         !CALL  construct_primitive_basis(psi_out%Basis)
         psi_out%Basis => psi_in%Basis
         !write(out_unit,*) 'END Coping psi_in in psi_out'
      END IF
      IF (allocated(psi_in%CVec)) THEN
         !write(out_unit,*) 'Coping psi_in in psi_out (complex):'
         psi_out%CVec(:) = psi_in%CVec(:)
         !CALL init_Basis1_TO_Basis2 (psi_in%Basis,psi_out%Basis)
         !CALL  construct_primitive_basis(psi_out%Basis)
         psi_out%Basis => psi_in%Basis
         !write(out_unit,*) 'END Coping psi_in in psi_out'
      END IF

   END SUBROUTINE copy_psi

   SUBROUTINE init_psi(psi, Basis, cplx, grid)
      USE Basis_m

      TYPE(psi_t), intent(inout)          :: psi
      TYPE(Basis_t), intent(in), target   :: Basis
      logical, intent(in)                 :: cplx, grid

      CALL dealloc_psi(psi)

      IF (.NOT. Basis_IS_allocated(Basis)) STOP 'ERROR in init_psi: the Basis is not initialized'

      IF (Basis%nb < 1) STOP 'ERROR in init_psi: Basis%nb < 1!'

      Psi%Basis => Basis
      Psi%Grid = Grid

      If (Grid) THEN !allocation on grid
         !print*,"psi is on Grid"
         IF (cplx) THEN

            IF (allocated(Basis%tab_basis)) THEN

               allocate (psi%CVec(Basis%nq*Basis%tab_basis(size(Basis%tab_basis))%nb))

            else

               allocate (psi%CVec(Basis%nq))

            END IF

         ELSE

            IF (allocated(Basis%tab_basis)) THEN

               allocate (psi%RVec(Basis%nq*Basis%tab_basis(size(Basis%tab_basis))%nb))

            ELSE

               allocate (psi%RVec(Basis%nq))

            END IF

         END IF

      ELSE ! allocation on basis
         !print*,"psi is on Basis"

         IF (cplx) THEN

            IF (allocated(Basis%tab_basis)) THEN

               allocate (psi%CVec(Basis%nb*Basis%tab_basis(size(Basis%tab_basis))%nb))

            else

               allocate (psi%CVec(Basis%nb))

            END IF

         ELSE

            IF (allocated(Basis%tab_basis)) THEN

               allocate (psi%RVec(Basis%nb*Basis%tab_basis(size(Basis%tab_basis))%nb))

            ELSE

               allocate (psi%RVec(Basis%nb))

            END IF

         END IF

      END IF

   END SUBROUTINE init_psi


   SUBROUTINE dealloc_psi(psi)

      TYPE(psi_t), intent(inout) :: psi

      nullify (psi%Basis)

      IF (allocated(psi%RVec)) THEN

         deallocate (psi%RVec)

      END IF

      IF (allocated(psi%CVec)) THEN

         deallocate (psi%CVec)

      END IF


   END SUBROUTINE dealloc_psi

   SUBROUTINE write_psi(psi, psi_cplx, print_psi_grid, print_basis, t, int_print, real_part)
      TYPE(psi_t), intent(in)                                 :: psi
      logical, intent(in)                                     :: print_psi_grid, print_basis, psi_cplx, real_part
      integer, intent(in), optional                           :: int_print
      real(kind= Rkind), intent(in), optional                 :: t
      logical, parameter                                      :: debug = .true.

      !local variable
      TYPE(psi_t)                                              :: psi_g, psi_b
      integer                                                  :: i

      if (print_psi_grid) then

         if (psi%Grid) then

            if (present(t) .and. present(int_print)) then

               call write_psi_grid(psi, print_cplx=psi_cplx, t=t, nio=int_print, real_part=real_part)

            elseif (present(t)) then

               call write_psi_grid(psi, print_cplx=psi_cplx, t=t, real_part=real_part)

            elseif (present(int_print)) then

               call write_psi_grid(psi, print_cplx=psi_cplx, nio=int_print, real_part=real_part)

            elseif (.not. present(t) .and. .not. present(int_print)) then

               call write_psi_grid(psi, print_cplx=psi_cplx, real_part=real_part)

            end if

         else

            call init_psi(psi_g, psi%Basis, cplx=.TRUE., grid=.true.)
            call BasisTOGrid_nD_cplx(psi_g%CVec, psi%CVec, psi%Basis)

            if (present(t) .and. present(int_print)) then

               call write_psi_grid(psi_g, print_cplx=psi_cplx, t=t, nio=int_print, real_part=real_part)

            elseif (present(t)) then

               call write_psi_grid(psi_g, print_cplx=psi_cplx, t=t, real_part=real_part)

            elseif (present(int_print)) then

               call write_psi_grid(psi_g, print_cplx=psi_cplx, nio=int_print, real_part=real_part)

            elseif (.not. present(t) .and. .not. present(int_print)) then

               call write_psi_grid(psi_g, print_cplx=psi_cplx, real_part=real_part)

            end if

         end if

         call dealloc_psi(psi_g)

      else

         if (psi%Grid) then

            call init_psi(psi_b, psi%Basis, cplx=.TRUE., grid=.false.)
            call GridTOBasis_nD_cplx(psi_b%CVec, psi%CVec, psi%Basis)

            if (present(t) .and. present(int_print)) then

               call write_psi_basis(psi_b, print_cplx=psi_cplx, t=t, nio=int_print, real_part=real_part)

            elseif (present(t)) then

               call write_psi_basis(psi_b, print_cplx=psi_cplx, t=t, real_part=real_part)

            elseif (present(int_print)) then

               call write_psi_basis(psi_b, print_cplx=psi_cplx, nio=int_print, real_part=real_part)

            elseif (.not. present(t) .and. .not. present(int_print)) then

               call write_psi_basis(psi_b, print_cplx=psi_cplx, real_part=real_part)

            end if

            call dealloc_psi(psi_b)

         else

            if (present(t) .and. present(int_print)) then

               call write_psi_basis(psi, print_cplx=psi_cplx, t=t, nio=int_print, real_part=real_part)

            elseif (present(t)) then

               call write_psi_basis(psi, print_cplx=psi_cplx, t=t, real_part=real_part)

            elseif (present(int_print)) then

               call write_psi_basis(psi, print_cplx=psi_cplx, nio=int_print, real_part=real_part)

            elseif (.not. present(t) .and. .not. present(int_print)) then

               call write_psi_basis(psi, print_cplx=psi_cplx, real_part=real_part)

            end if

         end if

      end if

      if (print_basis) then

         call Write_Basis(psi%Basis)

      end if

      if (debug) then

         flush (out_unit)

      end if


   END SUBROUTINE write_psi

   SUBROUTINE Projection(psi_dt_2, psi_dt_1)

      TYPE(psi_t), intent(in), target                  :: psi_dt_1
      TYPE(psi_t), intent(inout), target               :: psi_dt_2

      complex(kind= Rkind), pointer                    :: BBB1(:, :, :), BBB2(:, :, :)
      complex(kind= Rkind), allocatable, target        :: B1(:), B2(:)
      real(kind= Rkind)                                :: Norm0

      logical, parameter                               :: debug = .true.
      integer                                          :: inb, i1, i3, Ndim
      Integer, allocatable                             :: Ib1(:), Ib2(:), Ib3(:)


      call Calc_index(Ib1=Ib1, Ib2=Ib2, Ib3=Ib3, Basis=psi_dt_1%Basis)
      Ndim = size(psi_dt_1%Basis%tab_basis) - 1
      call Calc_Norm_OF_psi(psi_dt_1,Norm0)

      write (out_unit, *) 'Begin Hagedorn projection <psi|psi>=',Norm0
      !write (out_unit, *) 'out',psi_dt_2%CVec
      !call  Write_VecMat(psi_dt_1%Basis%tab_basis(1)%S, out_unit, 5,  info='S')

      If (Ndim == 1) then

         BBB1(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_1%CVec
         BBB2(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_2%CVec

         psi_dt_2%CVec(:) = CZERO

         DO i3 = 1, ubound(BBB1, dim=3)

         DO i1 = 1, ubound(BBB1, dim=1)

            BBB2(i1, :, i3) = matmul(psi_dt_1%Basis%tab_basis(1)%S,BBB1(i1, :, i3))

         END DO

         END DO

      else

         allocate (B1(Ib1(1)*Ib2(1)*Ib3(1)))
         BBB1(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_1%CVec
         BBB2(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => B1

         DO i3 = 1, ubound(BBB1, dim=3)

         DO i1 = 1, ubound(BBB1, dim=1)

            BBB2(i1, :, i3) = matmul(psi_dt_1%Basis%tab_basis(1)%S,BBB1(i1, :, i3))

         END DO

         END DO

         Do inb = 2, Ndim

            allocate (B2(Ib1(inb)*Ib2(inb)*Ib3(inb)))
            BBB1(1:Ib1(inb), 1:Ib2(inb), 1:Ib3(inb)) => B1
            BBB2(1:Ib1(inb), 1:Ib2(inb), 1:Ib3(inb)) => B2

            DO i3 = 1, ubound(BBB1, dim=3)

            DO i1 = 1, ubound(BBB1, dim=1)

               BBB2(i1, :, i3) = matmul( psi_dt_1%Basis%tab_basis(inb)%S,BBB1(i1, :, i3))

            END DO

            END DO

            B1 = B2
            B2 = CZERO
            deallocate (B2)
            if (inb == Ndim) psi_dt_2%CVec = B1

         END DO

      END IF

        call Calc_Norm_OF_psi(psi_dt_2,Norm0)

        deallocate(Ib1, Ib2, Ib3)
        if(allocated(B1)) deallocate(B1)
        if(allocated(B2)) deallocate(B2)

      write (out_unit, *) 'END Hagedorn projection <psi|psi>=',Norm0
      !write (out_unit, *) 'out',psi_dt_2%CVec

   END SUBROUTINE Projection


   SUBROUTINE Projection1(psi_dt_2, psi_dt_1)
   TYPE(psi_t), intent(in), target                  :: psi_dt_1
   TYPE(psi_t), intent(inout), target               :: psi_dt_2
   complex(kind= Rkind), pointer                    :: BBB1(:, :, :), BBB2(:, :, :)
   complex(kind= Rkind), allocatable, target        :: B1(:), B2(:)
   real(kind= Rkind)                                :: Norm0
   logical, parameter                               :: debug = .true.
   integer                                          :: inb, i1, i3, Ndim
   Integer, allocatable                             :: Ib1(:), Ib2(:), Ib3(:)
   call Calc_index(Ib1=Ib1, Ib2=Ib2, Ib3=Ib3, Basis=psi_dt_1%Basis)
   Ndim = size(psi_dt_1%Basis%tab_basis) - 1
   call Calc_Norm_OF_psi(psi_dt_1,Norm0)
   write (out_unit, *) 'Begin Hagedorn projection <psi|psi>=',Norm0
   !write (out_unit, *) 'out',psi_dt_2%CVec
   !call  Write_VecMat(psi_dt_1%Basis%tab_basis(1)%S, out_unit, 5,  info='S')
   If (Ndim == 1) then
      BBB1(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_1%CVec
      BBB2(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_2%CVec
      psi_dt_2%CVec(:) = CZERO
      DO i3 = 1, ubound(BBB1, dim=3)
      DO i1 = 1, ubound(BBB1, dim=1)
         BBB2(i1, :, i3) = matmul(psi_dt_1%Basis%tab_basis(1)%S,BBB1(i1, :, i3))
      END DO
      END DO
   else
      allocate (B1(Ib1(1)*Ib2(1)*Ib3(1)))
      BBB1(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => psi_dt_1%CVec
      BBB2(1:Ib1(1), 1:Ib2(1), 1:Ib3(1)) => B1
      DO i3 = 1, ubound(BBB1, dim=3)
      DO i1 = 1, ubound(BBB1, dim=1)
         BBB2(i1, :, i3) = matmul(psi_dt_1%Basis%tab_basis(1)%S,BBB1(i1, :, i3))
      END DO
      END DO
      Do inb = 2, Ndim
         allocate (B2(Ib1(inb)*Ib2(inb)*Ib3(inb)))
         BBB1(1:Ib1(inb), 1:Ib2(inb), 1:Ib3(inb)) => B1
         BBB2(1:Ib1(inb), 1:Ib2(inb), 1:Ib3(inb)) => B2
         DO i3 = 1, ubound(BBB1, dim=3)
         DO i1 = 1, ubound(BBB1, dim=1)
            BBB2(i1, :, i3) = matmul( psi_dt_1%Basis%tab_basis(inb)%S,BBB1(i1, :, i3))
         END DO
         END DO
         B1 = B2
         B2 = CZERO
         deallocate (B2)
         if (inb == Ndim) psi_dt_2%CVec = B1
      END DO
   END IF
     call Calc_Norm_OF_psi(psi_dt_2,Norm0)
     deallocate(Ib1, Ib2, Ib3)
     if(allocated(B1)) deallocate(B1)
     if(allocated(B2)) deallocate(B2)
   write (out_unit, *) 'END Hagedorn projection <psi|psi>=',Norm0
   !write (out_unit, *) 'out',psi_dt_2%CVec
END SUBROUTINE 


   SUBROUTINE Calc_Norm_OF_Psi(psi, Norm)

      implicit none
      type(Psi_t), intent(in)             :: psi
      real(kind= Rkind), intent(inout)    :: Norm

      IF (psi%Grid) THEN

         CALL Calc_Norm_OF_PsiGrid(Psi, Norm)

      ELSE

         CALL Calc_Norm_OF_PsiBasis(psi, Norm)

      END IF

      !write(out_unit,*) '<psi|psi> =',Norm

   END SUBROUTINE Calc_Norm_OF_Psi

   SUBROUTINE Calc_Norm_OF_PsiGrid(psi_g, Norm)

      USE QDUtil_m
      logical, parameter                   :: debug = .false.
      TYPE(Psi_t), intent(in)              :: psi_g
      real(kind= Rkind), intent(inout)     :: Norm


      complex(kind= Rkind), allocatable    :: psi_gb(:, :)
      logical                              :: Endloop_q    
      real(kind= Rkind), allocatable       :: Norme(:)
      real(kind= Rkind)                    :: WnD
      integer, allocatable                 :: Tab_iq(:)
      integer                              :: iq, inb, inbe

      IF (debug) THEN

         write (out_unit, *) 'Beging NormGrid'
         flush (out_unit)

      END IF

      Allocate (Psi_gb(Psi_g%Basis%nq, Psi_g%Basis%tab_basis(size(Psi_g%Basis%tab_basis))%nb))
      Allocate (Tab_iq(size(Psi_g%Basis%tab_basis) - 1))
      Allocate (Norme( psi_g%Basis%tab_basis(size(psi_g%Basis%tab_basis))%nb))

      Psi_gb(:, :) = reshape(Psi_g%CVec, shape=[psi_g%Basis%nq, Psi_g%Basis%tab_basis(size(Psi_g%Basis%tab_basis))%nb])

      DO inbe = 1, psi_g%Basis%tab_basis(size(psi_g%Basis%tab_basis))%nb !electronic state

         Norme(inbe) = ZERO

         Call Init_tab_ind(Tab_iq, psi_g%Basis%NDindexq)

         Iq = 0

         DO

            Iq = Iq + 1

            CALL increase_NDindex(Tab_iq, psi_g%Basis%NDindexq, Endloop_q)

            !if(inbe == 1)  print*,Tab_iq

            IF (Endloop_q) exit

            WnD = ONE

            DO inb = 1, size(psi_g%Basis%tab_basis) - 1

               WnD = WnD*Psi_g%Basis%tab_basis(inb)%w(Tab_iq(inb))

            END DO

            Norme(inbe) = Norme(inbe) + conjg(Psi_gb(iq, inbe))*Psi_gb(iq, inbe)*WnD

         END DO

         Norm = sqrt(sum(Norme))

      END DO

      Deallocate (Tab_iq,Norme)
      Deallocate (Psi_gb)

      IF (debug) THEN

         write (out_unit, *) 'END NormGrid'
         flush (out_unit)

      END IF

   END SUBROUTINE Calc_Norm_OF_PsiGrid


   SUBROUTINE Calc_Norm_OF_PsiBasis(psi, Norm)

      TYPE(psi_t), intent(in)          :: Psi
      real(kind= Rkind), intent(inout) :: Norm

      Norm = sqrt(real(dot_product(Psi%CVec, Psi%CVec), kind= Rkind))

      ! write(out_unit,*) 'norm,psi',Norm
   END SUBROUTINE Calc_Norm_OF_PsiBasis


   SUBROUTINE write_psi_basis(psi, print_cplx, t, nio, real_part)

      TYPE(psi_t), intent(in), target                :: psi
      real(kind= Rkind), intent(in), optional        :: t
      integer, intent(in), optional                  :: nio
      logical, intent(in)                            :: print_cplx, real_part

      complex(Kind= Rkind), pointer                  ::psibe(:, :)
      integer                                        :: ib, Ndim

      Ndim = size(psi%Basis%tab_basis)
      psibe(1:psi%Basis%nb, 1:psi%Basis%tab_basis(Ndim)%nb) => psi%CVec

      if (print_cplx) then 

         if (real_part) then

           ! write (*, *) '****************** Beging wrinting psi on basis****************************'

            if (present(nio) .and. present(t)) then

               do ib = 1, psi%Basis%nb

                  write (nio, *) ib, real(psibe(ib, :), kind= Rkind), aimag(psibe(ib, :)), t

               end do

            elseif (present(t)) then

               do ib = 1, psi%Basis%nb

                  write (*, *) ib, real(psibe(ib, :), kind= Rkind), aimag(psibe(ib, :)), t

               end do

            elseif (present(nio)) then

               do ib = 1, psi%Basis%nb

                  write (nio, *) ib, real(psibe(ib, :), kind= Rkind), aimag(psibe(ib, :))

               end do

            elseif (.not. present(nio) .and. .not. present(t)) then

               do ib = 1, psi%Basis%nb

                  write (out_unit, *) ib, real(psibe(ib, :), kind= Rkind), aimag(psibe(ib, :))

               end do

            else

               print *, 'no default case'

            end if

            !write (*, *) '****************** End wrinting psi on basis********************************'

         else

           ! write (*, *) '****************** Beging wrinting psi on basis****************************'

            if (present(nio) .and. present(t)) then

               do ib = 1, psi%Basis%nb

                  write (nio, *) t, ib, psibe(ib, :)

               end do

            elseif (present(t)) then

               do ib = 1, psi%Basis%nb

                  write (*, *) t, ib, psibe(ib, 1)

               end do

            elseif (present(nio)) then

               do ib = 1, psi%Basis%nb

                  write (nio, *) ib, psibe(ib, :)

               end do

            elseif (.not. present(nio) .and. .not. present(t)) then

               do ib = 1, psi%Basis%nb

                  write (out_unit, *) ib, psibe(ib, :)

               end do

            else

               print *, 'no default case'

            end if

            !write (*, *) '****************** End wrinting psi on basis********************************'

         end if

      else

         !write (*, *) '****************** Beging wrinting <psi/psi> on basis****************************'

         if (present(nio) .and. present(t)) then

            do ib = 1, psi%Basis%nb

               write (nio, *) t, ib, abs(psibe(ib, :))**2

            end do

         elseif (present(t)) then

            do ib = 1, psi%Basis%nb

               write (*, *) t, ib, abs(psibe(ib, :))**2

            end do

         elseif (present(nio)) then

            do ib = 1, psi%Basis%nb

               write (nio, *) ib, abs(psibe(ib, :))**2

            end do

         elseif (.not. present(nio) .and. .not. present(t)) then

            do ib = 1, psi%Basis%nb

               write (out_unit, *) ib, abs(psibe(ib, :))**2

            end do

         else

            print *, 'no default case'

         end if

      end if

   END SUBROUTINE write_psi_basis


   SUBROUTINE write_psi_grid(psi, print_cplx, t, nio, real_part)

      TYPE(psi_t), intent(in), target                 :: psi
      real(kind= Rkind), intent(in), optional         :: t
      integer, intent(in), optional                   :: nio
      logical, intent(in)                             :: print_cplx, real_part

      complex(Kind= Rkind), pointer                   :: psige(:, :)
      integer                                         :: iq, Ndim
      real(Kind= Rkind), allocatable                  :: Q(:, :)

      Ndim = size(psi%Basis%tab_basis)
      psige(1:psi%Basis%nq, 1:psi%Basis%tab_basis(Ndim)%nb) => psi%CVec

      call calc_Q_grid(Q, psi%Basis)

      if (print_cplx) then

         if (real_part) then

           ! write (*, *) '****************** Beging wrinting psi on Grid****************************'

            if (present(nio) .and. present(t)) then

               do iq = 1, psi%Basis%nq

                  write (nio, *) Q(iq, :), real(psige(iq, :), kind= Rkind), aimag(psige(iq, :)), t
                  if (mod(iq, psi%Basis%tab_basis(1)%nq) == 0) write (nio, *)

               end do

            elseif (present(t)) then

               do iq = 1, psi%Basis%nq

                  write (*, *) Q(iq, :), real(psige(iq, :), kind= Rkind), aimag(psige(iq, :)), t
                  if (mod(iq, psi%Basis%tab_basis(1)%nq) == 0) write (*, *)

               end do

            elseif (present(nio)) then

               do iq = 1, psi%Basis%nq

                  write (nio, *) Q(iq, :), real(psige(iq, :), kind= Rkind), aimag(psige(iq, :))
                  if (mod(iq, psi%Basis%tab_basis(1)%nq) == 0) write (nio, *)

               end do

            elseif (.not. present(nio) .and. .not. present(t)) then

               do iq = 1, psi%Basis%nq

                  write (out_unit, *) Q(iq, :), real(psige(iq, :), kind= Rkind), aimag(psige(iq, :))
                  if (mod(iq, psi%Basis%tab_basis(1)%nq) == 0) write (out_unit, *)

               end do

            else

               print *, 'no default case'

            end if

          !  write (*, *) '****************** End wrinting psi on Grid ********************************'

         else

           ! write (*, *) '****************** Beging wrinting psi on Grid****************************'

            if (present(nio) .and. present(t)) then

               do iq = 1, psi%Basis%nq

                  write (nio, *) Q(iq, :), psige(iq, :), t

               end do
            elseif (present(t)) then

               do iq = 1, psi%Basis%nq

                  write (*, *) Q(iq, :), psige(iq, :), t

               end do

            elseif (present(nio)) then

               do iq = 1, psi%Basis%nq

                  write (nio, *) Q(iq, :), psige(iq, :)

               end do

            elseif (.not. present(nio) .and. .not. present(t)) then

               do iq = 1, psi%Basis%nq

                  write (out_unit, *) Q(iq, :), psige(iq, :)

               end do

            else

               print *, 'no default case'

            end if

            !write (*, *) '****************** End wrinting psi on Grid ********************************'


         end if

      else

        ! write (*, *) '****************** Beging wrinting <psi/psi> on Grid ****************************'

         if (present(nio) .and. present(t)) then

            do iq = 1, psi%Basis%nq

               write (nio, *) t, Q(iq, :), abs(psige(iq, :))**2
               if (mod(iq, psi%Basis%tab_basis(1)%nq) == 0) write (nio, *)

            end do

         elseif (present(t)) then

            do iq = 1, psi%Basis%nq

               write (*, *) t, Q(iq, :), abs(psige(iq, :))**2
               if (mod(iq, psi%Basis%tab_basis(1)%nq) == 0) write (*, *)

            end do

         elseif (present(nio)) then

            do iq = 1, psi%Basis%nq

               write (nio, *) Q(iq, :), abs(psige(iq, :))**2

               if (mod(iq, psi%Basis%tab_basis(1)%nq) == 0) write (nio, *)

            end do

         elseif (.not. present(nio) .and. .not. present(t)) then

            do iq = 1, psi%Basis%nq

               write (out_unit, *) Q(iq, :), abs(psige(iq, :))**2
               if (mod(iq, psi%Basis%tab_basis(1)%nq) == 0) write (out_unit, *)

            end do

         else

            print *, 'no default case'

         end if

      end if
      deallocate(Q)

   END SUBROUTINE



   SUBROUTINE psi_init_GWP(psi, Tab_GWP)
      USE QDUtil_m
      TYPE(psi_t), intent(inout), target         :: psi
      TYPE(GWP_t)                                :: Tab_GWP(:)

      !------------------locale variables--------------------------------------------
      TYPE(psi_t), target                        :: psi0
      real(kind= Rkind), allocatable             :: Q(:, :)
      integer                                    :: iq, iGWO
      complex(Kind= Rkind), pointer              :: gb(:, :)
      real(Kind= Rkind)                          :: NormG, NormB
      integer                                    :: ndim,n_surf,nq,i_surf

      !initialisation-----------------------------------------------------------------
      CALL init_psi(psi0, psi%Basis, cplx=.true., grid=.true.)
      call calc_Q_grid(Q, psi%Basis)

      psi%CVec(:) = CZERO
      ndim = size(psi%Basis%tab_basis)-1
      n_surf = psi%Basis%tab_basis(ndim+1)%nb
      nq = psi%Basis%nq
      i_surf =Tab_GWP(1)%Elecindex
      gb(1:nq, 1:n_surf) => psi0%CVec
      gb(:, :) = CZERO
      DO iq = 1,nq

         gb(iq, i_surf) = calc_GWP(Tab_GWP(1), Q(iq, :))
         ! write(25,*)  Q(iq,:),   abs( gb(iq,Tab_GWP(1)%Elecindex))**2
         !if(mod(iq,10)==0) write(25,*)
      END DO
      !------------------transformation Grid to Basis-------------------------
      call Calc_Norm_OF_Psi(psi0, NormG)
      !psi0%CVec(:) = psi0%CVec(:)/NormG
      !call Calc_Norm_OF_Psi(psi0, NormG)
      write (out_unit, *)'------------------- Begining Norm Checking---------------------------------------'
      write (out_unit, *)
      write (out_unit, *) 'NormG = ', NormG

      call GridTOBasis_nD_cplx(psi%CVec, psi0%CVec, psi%Basis)

      call Calc_Norm_OF_Psi(psi, NormB)
      !psi%CVec(:) = psi%CVec(:)/NormB
      !call Calc_Norm_OF_Psi(psi, NormB)

      print *, 'NormB = ', NormB

      !psi%CVec(:) = psi%CVec(:)/NormB
      !call Calc_Norm_OF_Psi(psi, NormB)

      print *, 'NormB (renormed)= ', NormB
      write (out_unit, *)'------------------------------End Norm Checking---------------------------------'
      write (out_unit, *)
      deallocate (Q)
      CALL dealloc_psi(psi0)

   END SUBROUTINE

   subroutine test(Basis)
      implicit none
      TYPE(Basis_t), target, intent(in)     :: Basis
      TYPE(psi_t)                            :: psiB, psiG
      real(kind= Rkind)                           :: NormG, NormB
      CALL init_psi(psiB, Basis, cplx=.TRUE., grid=.false.)
      CALL init_psi(psiG, Basis, cplx=.TRUE., grid=.true.)

      psiB%CVec(:) = CZERO
      psiB%CVec(1) = CONE
      call Calc_Norm_OF_Psi(PsiB, NormB)
      psiB%CVec(:) = psiB%CVec(:)/NormB
      call Calc_Norm_OF_Psi(psiB, NormB)
      call BasisTOGrid_nD_cplx(PsiG%CVec, psiB%CVec, Basis)
      !call BasisTOGrid_nE_cplx(PsiG%CVec,psiB%CVec,Basis)
      call Calc_Norm_OF_Psi(psiG, NormG)
      print *, 'NormB = ', NormB
      print *, 'NormG = ', NormG
      print *, '--------------------------------------------------------'
      psiG%CVec(:) = CONE
      ! psiG%CVec(1) = CONE
      call Calc_Norm_OF_Psi(PsiG, NormG)
      psiG%CVec(:) = psiG%CVec(:)/NormG
      call GridTOBasis_nD_cplx(PsiB%CVec, psiG%CVec, Basis)
      call Calc_Norm_OF_Psi(psiB, NormB)
      call Calc_Norm_OF_Psi(psiG, NormG)
      print *, 'NormG = ', NormG
      print *, 'NormB = ', NormB
   end subroutine test



   subroutine Ecrire_psi(psi,nio,t)
      implicit none
      TYPE(psi_t)           , intent(in)     :: psi
      TYPE(psi_t)    ,target                        :: psiG
      real(kind= Rkind), allocatable         :: Q(:, :)
      complex(Kind= Rkind), pointer          :: gb(:, :)
      real(kind= Rkind),intent(in)           :: t
      integer,intent(in)                     :: nio
      integer                                :: iq,ib,nq,nsurf,ndim,nio_2,nq1
      CHARACTER(len=100)                     :: t_string
      CHARACTER(len=100)                     :: name_1,name_2

      ndim = size(psi%Basis%tab_basis)-1
      nsurf = psi%Basis%tab_basis(ndim+1)%nb
      nq = psi%Basis%nq
      nio_2 =nio+1
      nq1 = psi%Basis%tab_basis(1)%nq

      ! Convertir delta_t en chaîne de caractères
      WRITE(t_string, '(F12.6)') t 
      name_1 = 'density_1_t=' // TRIM(ADJUSTL(t_string)) // '.txt'
      name_2 = 'density_2_t=' // TRIM(ADJUSTL(t_string)) // '.txt'
      
      open (unit=nio, file=name_1)
      open (unit=nio_2, file=name_2)

      call calc_Q_grid(Q, psi%Basis)
      CALL init_psi(psiG, psi%Basis, cplx=.TRUE., grid=.true.)
      call BasisTOGrid_nD_cplx(psiG%CVec, psi%CVec, psi%Basis)

      gb(1:nq, 1:nsurf) => psiG%CVec
      DO iq = 1,nq
       write(nio,*)  Q(iq,:),    abs( gb(iq,1))**2
       write(nio_2,*)  Q(iq,:),   abs( gb(iq,2))**2
       if (mod(iq, nq1 ) == 0) write(nio,*)
       if (mod(iq, nq1 ) == 0)  write(nio_2,*)
      END DO

      
      deallocate (Q)
      CALL dealloc_psi(psiG)
   end subroutine 





end module psi_m
