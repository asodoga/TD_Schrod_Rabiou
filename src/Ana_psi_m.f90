! Module for wavefunction analysis
module Ana_psi_m
   USE Basis_m
   USE NDindex_m
   USE psi_m
   implicit none
   private
   public :: Population, Qpop
   public :: Calc_Av_imp_k_nD
   public :: Calc_Avg_A_nD, test_analyse_psi, Calc_AVQ_SQ_nD

contains

   ! Test routine for wavefunction analysis
   SUBROUTINE test_analyse_psi(psi)
      USE QDUtil_m
      USE Basis_m
      type(psi_t), intent(in) :: psi

      integer, allocatable :: Tab_Iq(:, :)
      real(kind=Rkind), allocatable :: SQt(:), Qt(:)
      complex(kind=Rkind), allocatable :: At(:)
      integer :: ndim
     
      ! Get number of dimensions (excluding electronic states)
      ndim = size(psi%Basis%tab_basis) - 1
      allocate(Qt(ndim), SQt(ndim), At(ndim))
   
      ! Calculate index table and average values
      call Calc_tab_Iq0(Tab_Iq, psi%Basis)
      call Calc_Avg_A_nD(psi, At)

      ! Output results
      write(out_unit,*) 'At', At
      write(out_unit,*) '=============================================================='
   
   END SUBROUTINE test_analyse_psi
   
   ! Calculate population in each electronic state
   subroutine Population(Psi, Pop)
      implicit none
      type(Psi_t), intent(in), target :: Psi
      complex(kind=Rkind), pointer :: Psi_bb(:, :)
      real(Kind=Rkind), intent(inout), allocatable :: Pop(:)
      integer :: inb, nsurf, ndim, nb
      real(Kind=Rkind) :: Norm

      ! Get dimensions
      ndim = size(Psi%Basis%tab_basis)-1
      nb = Psi%Basis%nb
      nsurf = Psi%Basis%tab_basis(ndim+1)%nb
      
      ! Set pointer to wavefunction coefficients
      Psi_bb(1:nb, 1:nsurf) => Psi%CVec
      
      ! Calculate norm
      call Calc_Norm_OF_Psi(Psi, Norm)

      ! Calculate population for each electronic state
      do inb = 1, nsurf
         Pop(inb) = real(dot_product(Psi_bb(:, inb), Psi_bb(:, inb)), kind=Rkind)/Norm
      end do

   end subroutine Population

   ! Calculate average position for each electronic state
   SUBROUTINE Qpop(Psi, Qp)
      USE Basis_m
      USE QDUtil_m
      type(Psi_t), intent(in), target :: Psi
      type(Psi_t), target :: Psi_g
      complex(kind=Rkind), pointer :: psi_gb(:, :)
      real(kind=Rkind), intent(inout) :: Qp(:)
      real(kind=Rkind) :: Norm(2)
      logical, parameter :: debug = .true.
      integer :: inb, ndim

      IF (debug) THEN
         flush(out_unit)
      END IF
      
      ! Get number of dimensions
      Ndim = size(Psi%Basis%tab_basis)
      
      ! Initialize grid wavefunction
      call init_psi(psi_g, psi%Basis, cplx=.TRUE., grid=.true.)
      Psi_g%CVec(:) = CZERO
      
      ! Transform from basis to grid representation
      call BasisTOGrid_nD_cplx(Psi_g%CVec, Psi%CVec, Psi%Basis)
      
      ! Calculate average position for each electronic state
      do inb = 1, Psi%Basis%tab_basis(2)%nb
         psi_gb(1:Psi%Basis%tab_basis(1)%nq, 1:Psi%Basis%tab_basis(2)%nb) => psi_g%CVec
         
         ! Calculate <x> and norm
         Qp(inb) = real(dot_product(psi_gb(:, inb), Psi%Basis%tab_basis(1)%w*Psi%Basis%tab_basis(1)%x*psi_gb(:, inb)), kind=Rkind)
         Norm(inb) = real(dot_product(psi_gb(:, inb), Psi%Basis%tab_basis(1)%w*psi_gb(:, inb)), kind=Rkind)
         
         ! Normalize if norm is non-zero
         if (abs(Norm(inb)) > tiny(1.0_Rkind)) then
            Qp(inb) = Qp(inb)/Norm(inb)
         end if
      end do

      ! Output results
      print *, Qp, Norm
      IF (debug) THEN
         flush(out_unit)
      END IF
   END SUBROUTINE Qpop

   ! Calculate average values for analysis
   SUBROUTINE Calc_Avg_A_nD(psi, At)
      USE QDUtil_m
      type(psi_t), intent(in) :: psi
      complex(kind=Rkind), intent(inout) :: At(:)
      
      ! Local variables
      real(kind=Rkind), allocatable :: VQ(:), VP(:), VQQ(:), VQP(:), SQt(:)
      real(kind=Rkind), allocatable :: CA(:), CB(:)
      integer, allocatable :: Tab_Iq(:, :)
      logical, parameter :: debug = .true.
      integer :: Ndim, Ib
      
      IF (debug) THEN
         flush(out_unit)
      END IF
      
      ! Get number of dimensions
      Ndim = size(psi%Basis%tab_basis) - 1 
      call Calc_tab_Iq0(Tab_Iq, psi%Basis)
      
      ! Allocate and initialize arrays
      allocate(VQ(Ndim), VQQ(Ndim), VP(Ndim), VQP(Ndim), SQt(Ndim))
      allocate(CA(Ndim), CB(Ndim))
      VQ(:) = ZERO; VP(:) = ZERO; VQQ(:) = ZERO; VQP(:) = ZERO; SQt(:) = ONE

      ! Calculate various averages
      call Calc_AVQ_SQ_nD_scd(psi, VQ, SQt, VQQ, Tab_Iq)
      call Calc_Av_imp_k_nD(psi, VP)
      call Calc_VQP_nD(VQP, psi)
      write(out_unit, *) 'VQQ = ', VQQ
      
      ! Calculate complex coefficients
      At(:) = CZERO 
      CB(:) = ZERO
      CA(:) = ZERO
      
      Do Ib = 1, Ndim
         CA(Ib) = SQt(Ib)**2
         CB(Ib) = (VQP(Ib)-TWO*VP(Ib)*VQ(Ib))/(TWO*(VQQ(Ib)-VQ(Ib)*VQ(Ib))) 
         At(Ib) = complex(CA(Ib), -CB(Ib))   
      End do
      
      ! Clean up
      deallocate(VQQ, VQP, CB, CA, VQ, VP, SQt)   
   END SUBROUTINE Calc_Avg_A_nD

   ! Calculate average momentum for 1D
   SUBROUTINE Calc_Av_imp_k_1D(psi0, K, nio)
      USE QDUtil_m
      type(psi_t), intent(in), target :: psi0
      real(kind=Rkind), intent(inout) :: K
      integer, intent(in) :: nio
      
      ! Local variables
      type(psi_t), target :: psi, ikpsi
      type(psi_t), target :: psi_b, ikpsi_b
      logical, parameter :: debug = .true.
      complex(kind=Rkind), pointer :: GB(:, :, :)
      complex(kind=Rkind), pointer :: d1gg(:, :, :)
      complex(kind=Rkind), pointer :: ikpsi0(:, :, :)
      Integer, allocatable :: Iq1(:), Iq2(:), Iq3(:)
      Integer, allocatable :: Ib1(:), Ib2(:), Ib3(:)
      integer :: i1, i3
      
      IF (debug) THEN
         flush(out_unit)
      END IF
      
      ! Initialize wavefunctions
      call init_psi(psi, psi0%Basis, cplx=.true., grid=.true.)
      call init_psi(ikpsi, psi0%Basis, cplx=.true., grid=.true.)
      call init_psi(psi_b, psi0%Basis, cplx=.true., grid=.false.)
      call init_psi(ikpsi_b, psi0%Basis, cplx=.true., grid=.false.)
      
      ! Calculate indices
      Call Calc_index(Ib1=Ib1, Ib2=Ib2, Ib3=Ib3, Iq1=Iq1, Iq2=Iq2, Iq3=Iq3, Basis=psi0%Basis)
      
      ! Initialize wavefunctions
      psi%CVec = CZERO
      ikpsi%CVec = CZERO
      ikpsi_b%CVec = CZERO
      psi_b%CVec = CZERO
      
      ! Transform to grid if needed
      IF (psi0%Grid) then
         psi%CVec(:) = psi0%CVec(:)
      ELSE
         call BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
      END IF
       
      ! Set pointers
      GB(1:Iq1(nio), 1:Iq2(nio), 1:Iq3(nio)) => psi%CVec
      ikpsi0(1:Iq1(nio), 1:Iq2(nio), 1:Iq3(nio)) => ikpsi%CVec
      d1gg(1:Iq2(nio), 1:Iq2(nio), 1:1) => psi%Basis%tab_basis(nio)%d1gg(:, :, 1)
      
      ! Calculate momentum operator acting on wavefunction
      DO i3 = 1, ubound(GB, dim=3)
         DO i1 = 1, ubound(GB, dim=1)
            ikpsi0(i1, :, i3) = ikpsi0(i1, :, i3)-EYE*matmul(d1gg(:,:,1), GB(i1,:,i3))
         END DO
      END DO
      
      ! Transform back to basis representation
      call GridTOBasis_nD_cplx(psi_b%CVec, psi%CVec, psi0%Basis)
      call GridTOBasis_nD_cplx(ikpsi_b%CVec, ikpsi%CVec, psi0%Basis)
      
      ! Calculate expectation value
      K = real(dot_product(psi_b%CVec, ikpsi_b%CVec), kind=Rkind)
      
      ! Clean up
      call dealloc_psi(psi)
      call dealloc_psi(ikpsi)
      call dealloc_psi(ikpsi_b)
      call dealloc_psi(psi_b)
      Deallocate(Ib1, Ib2, Ib3, Iq1, Iq2, Iq3)
   END SUBROUTINE Calc_Av_imp_k_1D

   ! Calculate average momentum for all dimensions
   SUBROUTINE Calc_Av_imp_k_nD(psi0, K)
      USE QDUtil_m
      type(psi_t), intent(in) :: psi0
      real(kind=Rkind), intent(inout) :: K(:)
      
      ! Local variables
      type(psi_t), target :: psi
      logical, parameter :: debug = .true.
      integer :: Ndim, Inb

      IF (debug) THEN
         flush(out_unit)
      END IF
      
      ! Get number of dimensions
      Ndim = size(psi0%Basis%tab_basis) - 1
      
      ! Initialize wavefunction
      call init_psi(psi, psi0%Basis, cplx=.TRUE., grid=.true.)
      psi%CVec(:) = CZERO
      K(:) = ZERO

      ! Transform to grid if needed
      IF (psi0%Grid) then
         psi%CVec(:) = psi0%CVec(:)
      ELSE
         call BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
      END IF

      ! Calculate momentum for each dimension
      Do Inb = 1, Ndim
         call Calc_Av_imp_k_1D(psi, K(inb), inb)
      End do

      ! Output results
      write(out_unit, *) '<psi/-id_xi/psi> =', K

      IF (debug) THEN
         flush(out_unit)
      END IF
      CALL dealloc_psi(psi)
   END SUBROUTINE Calc_Av_imp_k_nD

   ! Calculate VQP for 1D
   SUBROUTINE Calc_VQP_1D(psi0, VQP, Ib)
      USE QDUtil_m
      type(psi_t), intent(in), target :: psi0
      real(kind=Rkind), intent(inout) :: VQP
      integer, intent(in) :: Ib
      
      ! Local variables
      type(psi_t), target :: psi, d1psi
      complex(kind=Rkind), allocatable :: psi_gb(:, :), d1psi_gb(:, :)
      real(kind=Rkind), allocatable :: VQPEl(:)
      real(kind=Rkind), allocatable :: N(:)
      integer, allocatable :: Tab_iq(:)
      logical :: Endloop_q
      real(kind=Rkind) :: W, X
      complex(kind=Rkind), pointer :: GB(:, :, :)
      complex(kind=Rkind), pointer :: d1gg(:, :, :)
      complex(kind=Rkind), pointer :: d1psi0(:, :, :)
      Integer, allocatable :: Iq1(:), Iq2(:), Iq3(:)
      Integer, allocatable :: Ib1(:), Ib2(:), Ib3(:)
      integer :: i1, i3, inbe,inb,Iq

      ! Initialize wavefunctions
      call init_psi(psi, psi0%Basis, cplx=.true., grid=.true.)
      call init_psi(d1psi, psi0%Basis, cplx=.true., grid=.true.)
      
      ! Calculate indices
      Call Calc_index(Ib1=Ib1, Ib2=Ib2, Ib3=Ib3, Iq1=Iq1, Iq2=Iq2, Iq3=Iq3, Basis=psi0%Basis)

      ! Initialize wavefunctions
      psi%CVec = CZERO
      d1psi%CVec = CZERO

      ! Transform to grid if needed
      IF (psi0%Grid) then
         psi%CVec(:) = psi0%CVec(:)
      ELSE
         call BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
      END IF
     
      ! Set pointers
      GB(1:Iq1(Ib), 1:Iq2(Ib), 1:Iq3(Ib)) => psi%CVec
      d1psi0(1:Iq1(Ib), 1:Iq2(Ib), 1:Iq3(Ib)) => d1psi%CVec
      d1gg(1:Iq2(Ib), 1:Iq2(Ib), 1:1) => psi%Basis%tab_basis(Ib)%d1gg(:, :, 1)
      
      ! Calculate first derivative
      DO i3 = 1, ubound(GB, dim=3)
         DO i1 = 1, ubound(GB, dim=1)
            d1psi0(i1, :, i3) = d1psi0(i1, :, i3) + matmul(d1gg(:,:,1), GB(i1,:,i3))
         END DO
      END DO
     
      ! Reshape wavefunctions
      Allocate(psi_gb(psi%Basis%nq, psi%Basis%tab_basis(size(psi%Basis%tab_basis))%nb))
      Allocate(d1psi_gb(d1psi%Basis%nq, d1psi%Basis%tab_basis(size(d1psi%Basis%tab_basis))%nb))
      Allocate(VQPEl(psi%Basis%tab_basis(size(psi%Basis%tab_basis))%nb))
      Allocate(Tab_iq(size(Psi%Basis%tab_basis) - 1), N(size(Psi%Basis%tab_basis) - 1))

      psi_gb(:, :) = reshape(psi%CVec, shape=[psi%Basis%nq, psi%Basis%tab_basis(size(psi%Basis%tab_basis))%nb])
      d1psi_gb(:, :) = reshape(d1psi%CVec, shape=[d1psi%Basis%nq, d1psi%Basis%tab_basis(size(d1psi%Basis%tab_basis))%nb])
      N(:) = ZERO

      ! Calculate VQP for each electronic state
      DO inbe = 1, psi%Basis%tab_basis(size(psi%Basis%tab_basis))%nb
         VQPel(inbe) = ZERO
         Call Init_tab_ind(Tab_iq, psi%Basis%NDindexq)
         Iq = 0
         DO
            Iq = Iq + 1
            CALL increase_NDindex(Tab_iq, psi%Basis%NDindexq, Endloop_q)
            IF (Endloop_q) exit
            W = ONE
            DO inb = 1, size(psi%Basis%tab_basis) - 1
               W = W*psi%Basis%tab_basis(inb)%w(tab_iq(inb))
            END DO

            X = psi%Basis%tab_basis(ib)%x(tab_iq(ib))
            N(inbe) = N(inbe) + real(conjg(psi_gb(iq, inbe))*psi_gb(iq, inbe)*W, kind=Rkind)
            VQPel(inbe) = VQPel(inbe)-EYE*conjg(psi_gb(iq, inbe))*W*(psi_gb(iq, inbe)+TWO*X*d1psi_gb(iq, inbe))
         END DO
      END DO
      
      ! Calculate total VQP
      VQP = sum(real(VQPel, kind=Rkind))/(Sum(N)**2)
    
      ! Clean up
      call dealloc_psi(psi)
      call dealloc_psi(d1psi)
      deallocate(Ib1, Ib2, Ib3, Iq1, Iq2, Iq3)
      deallocate(N, VQPEl, psi_gb, d1psi_gb, Tab_iq)
   END SUBROUTINE Calc_VQP_1D

   ! Calculate VQP for all dimensions
   SUBROUTINE Calc_VQP_nD(VQP, psi)
      TYPE(psi_t), intent(in) :: psi
      real(kind=Rkind), intent(inout) :: VQP(:)
      integer :: Ib, ndim

      ! Get number of dimensions
      ndim = size(psi%Basis%tab_basis) - 1

      ! Calculate VQP for each dimension
      Do Ib = 1, ndim
         call Calc_VQP_1D(psi, VQP(Ib), Ib)
      End do
   END SUBROUTINE Calc_VQP_nD

   ! Calculate average position and spread
   SUBROUTINE Calc_AVQ_SQ_nD(psi_in, AVQ, SQ, Tab_Iq)
      USE QDUtil_m
      TYPE(Psi_t), intent(in) :: psi_in
      integer, intent(in) :: Tab_Iq(:, :)
      real(kind=Rkind), intent(inout) :: AVQ(:), SQ(:)

      ! Local variables
      logical, parameter :: debug = .false.
      TYPE(Psi_t) :: psi
      complex(kind=Rkind), allocatable :: psi_gb(:, :)
      real(kind=Rkind), allocatable :: AVQel(:, :), SQel(:, :), Q(:)
      real(kind=Rkind), allocatable :: N(:)
      integer :: iq, inbe, inb, nsurf, nq, ndim
      real(kind=Rkind) :: W

      IF (debug) THEN
         write(out_unit, *) 'Beging AVQ'
         flush(out_unit)
      END IF

      ! Get dimensions
      ndim = size(psi_in%Basis%tab_basis) - 1
      nq = psi_in%Basis%nq
      nsurf = psi_in%Basis%tab_basis(ndim+1)%nb
      
      ! Allocate arrays
      allocate(N(nsurf), AVQel(ndim, nsurf), SQel(ndim, nsurf), Q(ndim), psi_gb(nq, nsurf))
      
      ! Initialize grid wavefunction
      CALL init_psi(psi, psi_in%Basis, cplx=.true., grid=.true.)

      ! Transform to grid if needed
      IF (psi_in%Grid) then
         psi%CVec(:) = psi_in%CVec(:)
      ELSE
         CALL BasisTOGrid_nD_cplx(psi%CVec, psi_in%CVec, psi_in%Basis)
      END IF

      ! Reshape wavefunction
      psi_gb(:, :) = reshape(psi%CVec, shape=[nq, nsurf])
      N(:) = ZERO
      AVQel(:, :) = ZERO
      SQel(:, :) = ZERO

      ! Calculate averages for each electronic state
      DO inbe = 1, nsurf 
         DO Iq = 1, nq
            w = ONE
            Q(:) = ZERO

            ! Calculate weight and position
            DO inb = 1, ndim
               w = w*psi%Basis%tab_basis(inb)%w(Tab_Iq(inb,Iq))
               Q(inb) = psi%Basis%tab_basis(inb)%x(Tab_Iq(inb,Iq))
            END DO
         
            ! Accumulate values
            N(inbe) = N(inbe) + real(conjg(psi_gb(iq, inbe))*psi_gb(iq, inbe)*W, kind=Rkind)
            AVQel(:, inbe) = AVQel(:, inbe) + real(conjg(psi_gb(iq, inbe))*(Q(:)*w)*psi_gb(iq, inbe), kind=Rkind)
            SQel(:, inbe) = SQel(:, inbe) + real(conjg(psi_gb(iq, inbe))*(Q(:)*Q(:)*w)*psi_gb(iq, inbe), kind=Rkind)
         END DO
      END DO

      ! Calculate final averages
      DO inb = 1, ndim
         AVQ(inb) = sum(AVQel(inb,:))/(Sum(N)**2)
         SQ(inb) = sum(SQel(inb, :))/(Sum(N)**2)
         SQ(inb) = TWO*(SQ(inb) - AVQ(inb)*AVQ(inb))
         SQ(inb) = sqrt(ONE/SQ(inb))
      END DO   

      ! Clean up
      Deallocate(Psi_gb, N, AVQel, SQel, Q)
      CALL dealloc_psi(psi)
      
      IF (debug) THEN
         write(out_unit, *) 'END AVQ'
         flush(out_unit)
      END IF
   END SUBROUTINE Calc_AVQ_SQ_nD

   ! Calculate average position and spread (second version)
   SUBROUTINE Calc_AVQ_SQ_nD_scd(psi_in, AVQ, SQ, VQQ, Tab_Iq)
      USE QDUtil_m
      TYPE(Psi_t), intent(in) :: psi_in
      integer, intent(in) :: Tab_Iq(:, :)
      real(kind=Rkind), intent(inout) :: AVQ(:), SQ(:), VQQ(:)

      ! Local variables
      logical, parameter :: debug = .false.
      TYPE(Psi_t) :: psi
      complex(kind=Rkind), allocatable :: psi_gb(:, :)
      real(kind=Rkind), allocatable :: AVQel(:, :), SQel(:, :), Q(:)
      real(kind=Rkind), allocatable :: N(:)
      integer :: iq, inbe, inb, nsurf, nq, ndim
      real(kind=Rkind) :: W

      IF (debug) THEN
         write(out_unit, *) 'Beging AVQ'
         flush(out_unit)
      END IF

      ! Get dimensions
      ndim = size(psi_in%Basis%tab_basis) - 1
      nq = psi_in%Basis%nq
      nsurf = psi_in%Basis%tab_basis(ndim+1)%nb
      
      ! Allocate arrays
      allocate(N(nsurf), AVQel(ndim, nsurf), SQel(ndim, nsurf), Q(ndim), psi_gb(nq, nsurf))
      
      ! Initialize grid wavefunction
      CALL init_psi(psi, psi_in%Basis, cplx=.true., grid=.true.)

      ! Transform to grid if needed
      IF (psi_in%Grid) then
         psi%CVec(:) = psi_in%CVec(:)
      ELSE
         CALL BasisTOGrid_nD_cplx(psi%CVec, psi_in%CVec, psi_in%Basis)
      END IF

      ! Reshape wavefunction
      psi_gb(:, :) = reshape(psi%CVec, shape=[nq, nsurf])
      N(:) = ZERO
      AVQel(:, :) = ZERO
      SQel(:, :) = ZERO

      ! Calculate averages for each electronic state
      DO inbe = 1, nsurf 
         DO Iq = 1, nq
            w = ONE
            Q(:) = ZERO

            ! Calculate weight and position
            DO inb = 1, ndim
               w = w*psi%Basis%tab_basis(inb)%w(Tab_Iq(inb,Iq))
               Q(inb) = psi%Basis%tab_basis(inb)%x(Tab_Iq(inb,Iq))
            END DO
         
            ! Accumulate values
            N(inbe) = N(inbe) + real(conjg(psi_gb(iq, inbe))*psi_gb(iq, inbe)*W, kind=Rkind)
            AVQel(:, inbe) = AVQel(:, inbe) + real(conjg(psi_gb(iq, inbe))*(Q(:)*w)*psi_gb(iq, inbe), kind=Rkind)
            SQel(:, inbe) = SQel(:, inbe) + real(conjg(psi_gb(iq, inbe))*(Q(:)*Q(:)*w)*psi_gb(iq, inbe), kind=Rkind)
         END DO
      END DO

      ! Calculate final averages
      DO inb = 1, ndim
         AVQ(inb) = sum(AVQel(inb,:))/(Sum(N)**2)
         VQQ(inb) = sum(SQel(inb,:))/(Sum(N)**2)
         SQ(inb) = sum(SQel(inb, :))/(Sum(N)**2)
         SQ(inb) = sqrt(SQ(inb) - AVQ(inb)*AVQ(inb))
         SQ(inb) = ONE/(SQ(inb)*sqrt(TWO))
      END DO   

      ! Clean up
      Deallocate(Psi_gb, N, AVQel, SQel, Q)
      CALL dealloc_psi(psi)
      
      IF (debug) THEN
         write(out_unit, *) 'END AVQ'
         flush(out_unit)
      END IF
   END SUBROUTINE Calc_AVQ_SQ_nD_scd

end module Ana_psi_m