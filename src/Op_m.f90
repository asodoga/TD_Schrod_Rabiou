!==============================================================================
! Module Op_m - Opérateurs quantiques et fonctions associées
!
! Ce module définit les opérateurs quantiques et fournit des routines pour leur
! manipulation dans le cadre de calculs de chimie quantique.
!==============================================================================
module Op_m
   !$ USE omp_lib                     ! Pour le parallélisme OpenMP
   USE QDUtil_m                      ! Utilitaires généraux
   USE Basis_m, only: Basis_t        !  base q
   USE Molec_m                       ! Données moléculaires
   Use Psi_m                         ! Fonctions d'onde
   implicit none
   private

   !---------------------------------------------------------------------------
   ! TYPE Op_t - Structure pour représenter un opérateur 
   !
   ! Contient :
   ! - Basis : pointeur vers la base associée
   ! - RMat : matrice réelle de l'opérateur
   ! - CMat : matrice complexe de l'opérateur
   ! - Scal_pot : potentiel scalaire (pour les surfaces d'énergie potentielle)
   !---------------------------------------------------------------------------
   TYPE :: Op_t
      TYPE(Basis_t), pointer           :: Basis      ! Base associée
      real(kind=Rkind), allocatable    :: RMat(:, :) ! Matrice réelle
      complex(kind=Rkind), allocatable :: CMat(:, :) ! Matrice complexe
      real(kind=Rkind), allocatable    :: Scal_pot(:, :, :) ! Potentiel scalaire
   END TYPE Op_t

   ! Liste des procédures publiques
   public :: Op_t, write_Op, Set_Op, dealloc_Op, calc_OpPsi, Calc_Hpsi, Kpsi_nD, Make_Mat_H
   public :: test_op, Calc_Scalar_Pot, test_openmp_op, calc_tab_Iq, calc_nac, calc_VV, Popu

contains
   !===========================================================================
   ! SUBROUTINE alloc_Op - Allocation mémoire pour un opérateur
   !
   ! Paramètres :
   ! - Op : opérateur à allouer (inout)
   ! - nb : taille de la matrice à allouer (in)
   !===========================================================================
   SUBROUTINE alloc_Op(Op, nb)
      TYPE(Op_t), intent(inout) :: Op  ! Opérateur à allouer
      integer, intent(in)    :: nb     ! Taille de la matrice

      ! Désallocation préalable pour éviter les fuites mémoire
      CALL dealloc_Op(Op)

      ! Vérification de la taille
      IF (nb < 1) STOP 'ERROR in init_Op: nb < 1!'

      ! Allocation de la matrice réelle
      allocate (Op%RMat(nb, nb))
   END SUBROUTINE alloc_Op

   !===========================================================================
   ! SUBROUTINE dealloc_Op - Désallocation mémoire d'un opérateur
   !
   ! Paramètres :
   ! - Op : opérateur à désallouer (inout)
   !===========================================================================
   SUBROUTINE dealloc_Op(Op)
      TYPE(Op_t), intent(inout) :: Op  ! Opérateur à désallouer

      ! Déréférencement de la base
      nullify (Op%Basis)

      ! Désallocation de la matrice réelle si allouée
      IF (allocated(Op%RMat)) THEN
         deallocate (Op%RMat)
      END IF
   END SUBROUTINE dealloc_Op

   !===========================================================================
   ! SUBROUTINE write_Op - Écriture des informations d'un opérateur
   !
   ! Paramètres :
   ! - Op : opérateur à afficher (in)
   !===========================================================================
   SUBROUTINE write_Op(Op)
      USE QDUtil_m
      TYPE(Op_t), intent(in) :: Op  ! Opérateur à afficher

      ! Vérification de l'association à une base
      IF (associated(Op%Basis)) THEN
         write (out_unit, *) ' The basis is linked to Op.'
      END IF

      ! Affichage de la matrice réelle si allouée
      IF (allocated(Op%RMat)) THEN
         write (out_unit, *) 'Writing Op (real):'
         write (out_unit, *)
         CALL Write_VecMat(Op%RMat, out_unit, 5, info='Op%Rmat')
         write (out_unit, *) 'END Writing Op'
      END IF
   END SUBROUTINE write_Op

   !===========================================================================
   ! SUBROUTINE Make_Mat_H - Construction de la matrice Hamiltonienne
   !
   ! Paramètres :
   ! - Basis : base quantique (in)
   ! - H : Hamiltonien (inout)
   !===========================================================================
   SUBROUTINE Make_Mat_H(Basis,H)
      Use QDUtil_m
      USE Basis_m
      USE Psi_m

      TYPE(Basis_t),intent(in)            :: Basis  ! Base quantique
      TYPE(Op_t), intent(inout)           :: H      ! Hamiltonien

      ! Variables locales
      logical, parameter                  :: debug = .false.  ! Mode debug
      integer, allocatable                :: Tab_iq(:, :)     ! Table d'indices
      TYPE(psi_t)                         :: psi, Hpsi        ! Fonctions d'onde
      integer                             :: nb, nsurf, ib, jb, ndim

      IF (debug) THEN
         write (out_unit, *) 'BEGINNING Make_Mat_H'
         call write_basis(Basis)
         flush (out_unit)
      END IF

      ! Initialisation des fonctions d'onde
      call init_psi(psi, Basis, cplx=.true., grid=.false.)
      call init_psi(Hpsi, Basis, cplx=.true., grid=.false.)

      ! Détermination des dimensions
      ndim = size(Basis%tab_basis) - 1
      nsurf = Basis%tab_basis(ndim+1)%nb
      nb = (Basis%nb)*nsurf
      
      ! Allocation et initialisation
      allocate(H%CMat(nb,nb))
      call Calc_tab_Iq0(Tab_Iq, Basis)
      call Set_Op(H, Basis, Tab_Iq)

      ! Construction de la matrice Hamiltonienne
      DO jb = 1, nb
         psi%CVec(:) = CZERO
         psi%CVec(jb) = CONE
         call calc_Oppsi(H, psi, Hpsi)

         DO ib = 1, nb
            psi%CVec(:) = CZERO
            psi%CVec(ib) = CONE
            H%CMat(ib, jb) = dot_product(psi%CVec, Hpsi%CVec)
         END DO
      END DO

      ! Nettoyage mémoire
      call dealloc_psi(psi)
      call dealloc_psi(Hpsi)
      deallocate(Tab_Iq)

      IF (debug) THEN
         write (out_unit, *) 'END Make_Mat_H'
         flush (out_unit)
      END IF
   END SUBROUTINE Make_Mat_H

   !===========================================================================
   ! SUBROUTINE Set_Op - Initialisation d'un opérateur
   !
   ! Paramètres :
   ! - Op : opérateur à initialiser (inout)
   ! - Basis : base  (in)
   ! - Tab_Iq : table d'indices (in)
   !===========================================================================
   SUBROUTINE Set_Op(Op, Basis, Tab_Iq)
      USE Basis_m
      TYPE(Op_t), intent(inout)         :: Op      ! Opérateur
      TYPE(Basis_t), intent(in), target :: Basis   ! Base quantique
      integer, intent(in)               :: Tab_iq(:, :)  ! Table d'indices

      Integer :: ndim, nsurf, nq

      ! Vérification de l'allocation de la base
      IF (.NOT. Basis_IS_allocated(Basis)) THEN
         STOP 'ERROR in Set_Op: the Basis is not initialized'
      END IF

      ! Détermination des dimensions
      ndim = size(Basis%tab_basis) - 1
      nq = Basis%nq
      nsurf = Basis%tab_basis(ndim+1)%nb

      ! Calcul du potentiel scalaire avec OpenMP
      call Calc_Scalar_Pot_openmp(Op%Scal_pot, Basis, Tab_Iq)
   END SUBROUTINE Set_Op

   !===========================================================================
   ! SUBROUTINE Kpsi_nD - Application de l'opérateur énergie cinétique
   !
   ! Paramètres :
   ! - KPsi_g : résultat de l'application (inout)
   ! - Psi_g : fonction d'onde d'entrée (in)
   ! - Basis : base  (in)
   !===========================================================================
   SUBROUTINE Kpsi_nD(KPsi_g, Psi_g, Basis)
      USE Basis_m
      USE QDUtil_m
      USE Molec_m
      TYPE(Basis_t), intent(in), target               :: Basis
      complex(kind=Rkind), intent(in), target         :: psi_g(:)
      complex(kind=Rkind), intent(inout), target      :: Kpsi_g(:)
      
      ! Pointeurs et variables locales
      complex(kind=Rkind), pointer :: psi_ggb(:, :, :), d2gg(:, :), Kpsi_ggb(:, :, :)
      real(kind=Rkind), allocatable :: GGdef(:, :)
      logical, parameter :: debug = .true.
      integer :: iq, i1, i3, inb, ndim
      integer, allocatable :: Iq1, Iq2, Iq3

      IF (debug) THEN
         flush (out_unit)
      END IF
      
      ! Initialisation
      ndim = size(Basis%tab_basis)-1
      allocate(GGdef(ndim, ndim))
      CALL get_Qmodel_GGdef(GGdef)

      Kpsi_g(:) = CZERO
      
      ! Application de l'opérateur énergie cinétique pour chaque dimension
      DO inb = 1, ndim 
         ! Détermination des tailles selon la dimension courante
         IF (inb == 1) THEN
            Iq1 = 1
            Iq2 = Basis%tab_basis(1)%nq
            Iq3 = Product(Basis%tab_basis(2:ndim)%nq)*Basis%tab_basis(ndim + 1)%nb
         ELSE IF (inb == ndim) THEN
            Iq1 = Product(Basis%tab_basis(1:ndim - 1)%nq)
            Iq2 = Basis%tab_basis(ndim)%nq
            Iq3 = Basis%tab_basis(ndim + 1)%nb
         ELSE
            Iq1 = Product(Basis%tab_basis(1:inb - 1)%nq)
            Iq2 = Basis%tab_basis(inb)%nq
            Iq3 = Product(Basis%tab_basis(inb + 1:ndim)%nq)*Basis%tab_basis(ndim + 1)%nb
         END IF

         ! Association des pointeurs
         Kpsi_ggb(1:Iq1, 1:Iq2, 1:Iq3) => Kpsi_g
         psi_ggb(1:Iq1, 1:Iq2, 1:Iq3) => psi_g
         d2gg(1:Iq2, 1:Iq2) => Basis%tab_basis(inb)%d2gg

         ! Calcul pour chaque bloc
         DO i3 = 1, ubound(psi_ggb, dim=3)
            DO i1 = 1, ubound(psi_ggb, dim=1)
               KPsi_ggb(i1, :, i3) = KPsi_ggb(i1, :, i3) - HALF*GGdef(inb, inb)*matmul(d2gg, psi_ggb(i1, :, i3))
            END DO
         END DO
      END DO
      
      ! Nettoyage
      deallocate(Iq1, Iq2, Iq3)
      IF (debug) THEN
         flush (out_unit)
      END IF
   END SUBROUTINE Kpsi_nD

   !===========================================================================
   ! SUBROUTINE Calc_Scalar_Pot - Calcul du potentiel scalaire
   !
   ! Paramètres :
   ! - V : potentiel à calculer (inout)
   ! - Basis : base  (in)
   !===========================================================================
   SUBROUTINE Calc_Scalar_Pot(V, Basis)
      USE QDUtil_m
      USE Basis_m
      TYPE(Basis_t), intent(in), target :: Basis
      real(kind=Rkind), allocatable, intent(inout) :: V(:, :, :)
      
      ! Variables locales
      integer, allocatable :: Tab_iq(:)
      integer :: inb, ndim, iq, nq, nsurf, i
      real(Kind=Rkind), allocatable :: Q(:)
      logical :: Endloop

      ! Initialisation
      ndim = size(Basis%tab_basis) - 1
      nq = Basis%nq
      nsurf = Basis%tab_basis(ndim+1)%nb

      ! Ouverture des fichiers de sortie
      open (unit = 200, file = "Vm.txt")
      open (unit = 205, file = "x.txt")
      open (unit = 206, file = "y.txt")
      open (unit = 201, file = "V22-ad=f.txt")
      open (unit = 202, file = "V12-ad=f.txt")

      ! Allocation et initialisation
      allocate(Q(ndim), Tab_iq(ndim))
      allocate(V(nq, nsurf, nsurf))

      ! Écriture des coordonnées
      do i = 1, Basis%tab_basis(1)%nq
         write(205,*) Basis%tab_basis(1)%x(i) 
         write(206,*) Basis%tab_basis(2)%x(i) 
      end do

      ! Calcul du potentiel pour chaque point de la grille
      Call Init_tab_ind(Tab_iq, Basis%NDindexq)
      Iq = 0
      DO
         Iq = Iq + 1
         i = Tab_iq(2)
         CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop)
         IF (Endloop) exit
         
         ! Récupération des coordonnées
         do inb = 1, ndim
            Q(inb) = Basis%tab_basis(inb)%x(Tab_iq(inb))
         end do
         
         ! Calcul du potentiel
         CALL sub_Qmodel_V(V(iq, :, :), Q(:))
         
         ! Gestion des sauts de ligne pour la visualisation
         if ((Tab_iq(2)/= i)) then 
            write(200,*)
            write(201,*)
            write(202,*)
         end if    
         
         ! Écriture des résultats
         write(200,*) Q(:), V(iq, 1, 1)
         write(201,*) Q(:), V(iq, 2, 2)
         write(202,*) Q(:), V(iq, 2, 1)
      END DO

      ! Nettoyage
      deallocate(Tab_iq, Q)
   END SUBROUTINE Calc_Scalar_Pot

   !===========================================================================
   ! SUBROUTINE test_op - Tests des opérateurs
   !
   ! Paramètres :
   ! - Basis : base  (in)
   ! - psi0 : fonction d'onde de test (in)
   !===========================================================================
   SUBROUTINE test_op(Basis, psi0)
      USE QDUtil_m
      USE Basis_m
      TYPE(Basis_t), intent(in), target :: Basis
      TYPE(psi_t), intent(in) :: psi0
      
      ! Variables locales
      TYPE(Op_t) :: H
      TYPE(psi_t) :: psi
      real(kind=Rkind), allocatable :: V(:, :, :)
      complex(kind=Rkind), allocatable :: CEigVal(:), CEigVec(:, :)
      real(kind=Rkind), allocatable :: prob(:), vec(:), Pop(:)
      integer :: nb, ndim, ib, nsurf, jb

      ! Ouverture des fichiers de sortie
      open(unit=50, file='proj.txt')
      open(unit=51, file='EingVec.txt')
      open(unit=52, file='population.txt')

      ! Construction de la matrice Hamiltonienne
      call Make_Mat_H(Basis, H)
      ndim = size(Basis%tab_basis) - 1
      nb = Basis%nb*Basis%tab_basis(ndim+1)%nb
      nsurf = Basis%tab_basis(ndim+1)%nb
      
      ! Allocation mémoire
      allocate(CEigVal(nb), CEigVec(nb, nb), prob(nb), vec(nb), Pop(nsurf))
      call init_psi(psi, psi0%Basis, cplx=.true., grid=.false.)
      
      ! Diagonalisation
      call diagonalization(H%CMat, CEigVal, CEigVec)

      ! Calcul des projections et populations
      Do ib = 1, nb
         psi%CVec(:) = CEigVec(:, ib)
         prob(ib) = sum(CEigVec(:, ib)*psi0%CVec(:))
         call Popu(psi, Pop)
         write(51,*) ib, CEigVal(ib)%re
         write(50,*) ib, prob(ib)
         write(52,*) ib, Pop(:)
      End Do

      ! Nettoyage
      call dealloc_psi(psi)
   END SUBROUTINE test_op

   !===========================================================================
   ! SUBROUTINE Popu - Calcul des populations
   !
   ! Paramètres :
   ! - Psi : fonction d'onde (in)
   ! - Pop : populations calculées (inout)
   !===========================================================================
   subroutine Popu(Psi, Pop)
      implicit none
      type(Psi_t), intent(in), target :: Psi
      real(Kind=Rkind), intent(inout), allocatable :: Pop(:)
      
      ! Variables locales
      complex(kind=Rkind), pointer :: Psi_bb(:, :)
      integer :: inb, nsurf, ndim, nb
      real(Kind=Rkind) :: Norm

      ! Initialisation
      ndim = size(Psi%Basis%tab_basis)-1
      nb = Psi%Basis%nb
      nsurf = Psi%Basis%tab_basis(ndim+1)%nb
      Psi_bb(1:nb, 1:nsurf) => Psi%CVec
      
      ! Calcul de la norme
      call Calc_Norm_OF_Psi(Psi, Norm)
      
      ! Calcul des populations
      do inb = 1, nsurf
         Pop(inb) = dot_product(Psi_bb(:, inb), Psi_bb(:, inb))/Norm
      end do
   end subroutine Popu

   !===========================================================================
   ! SUBROUTINE Calc_Hpsi - Application de l'Hamiltonien
   !
   ! Paramètres :
   ! - psi_g : fonction d'onde en représentation grille (in)
   ! - HPsi_g : résultat de l'application (inout)
   ! - Basis : (in)
   ! - V : potentiel (in)
   !===========================================================================
   SUBROUTINE Calc_Hpsi(psi_g, HPsi_g, Basis, V)
      USE Basis_m
      USE Psi_m
      USE Molec_m
      complex(kind=Rkind), intent(in), target :: psi_g(:)
      complex(kind=Rkind), intent(inout) :: HPsi_g(:)
      type(Basis_t), intent(in), target :: Basis
      real(kind=Rkind), intent(in) :: V(:, :, :)
      
      ! Variables locales
      complex(kind=Rkind), allocatable, target :: VPsi_g(:), KPsi_g(:)
      complex(kind=Rkind), pointer :: VPsi_gb(:, :), Psi_gb(:, :)
      Integer :: i, j, ndim, nsurf, nq

      ! Vérification de l'allocation
      IF (.not. allocated(Basis%tab_basis)) THEN
         STOP 'ERROR in Set_Op: the Basis%tab_bais is not initialized'
      END IF

      ! Initialisation
      ndim = SIZE(Basis%tab_basis) - 1
      nq = Basis%nq
      nsurf = Basis%tab_basis(ndim+1)%nb

      ! Allocation mémoire
      allocate(VPsi_g(nq*nsurf))
      allocate(KPsi_g(nq*nsurf))
      VPsi_g(:) = CZERO
      KPsi_g(:) = CZERO
      HPsi_g(:) = CZERO

      ! Application du potentiel
      VPsi_gb(1:nq, 1:nsurf) => VPsi_g
      Psi_gb(1:nq, 1:nsurf) => Psi_g

      DO i = 1, nsurf
         DO j = 1, nsurf
            VPsi_gb(:, i) = VPsi_gb(:, i) + V(:, i, j)*Psi_gb(:, j)
         END DO
      END DO

      ! Application de l'énergie cinétique
      call Kpsi_nD(KPsi_g, Psi_g, Basis)
      HPsi_g(:) = VPsi_g(:) + KPsi_g(:)

      ! Nettoyage
      DEALLOCATE(VPsi_g)
      DEALLOCATE(KPsi_g)
   END SUBROUTINE Calc_Hpsi

   !===========================================================================
   ! SUBROUTINE calc_OpPsi - Application d'un opérateur à une fonction d'onde
   !
   ! Paramètres :
   ! - Op : opérateur à appliquer (in)
   ! - psi : fonction d'onde d'entrée (in)
   ! - Oppsi : résultat de l'application (inout)
   !===========================================================================
   SUBROUTINE calc_OpPsi(Op, psi, Oppsi)
      USE psi_m
      TYPE(Op_t), intent(in) :: Op
      TYPE(psi_t), intent(in) :: psi
      TYPE(psi_t), intent(inout) :: Oppsi
      
      ! Variables locales
      TYPE(psi_t) :: psi_g, Oppsi_g

      ! Cas d'une matrice réelle
      IF (allocated(Op%RMat)) THEN
         Oppsi%CVec = matmul(Op%RMat, psi%CVec)
      ELSE
         ! Cas général avec transformation grille/base
         IF (psi%Grid) THEN
            call init_psi(Oppsi_g, Psi%Basis, .true., .true.)
            call Calc_Hpsi(psi%CVec, Oppsi_g%CVec, psi%Basis, Op%Scal_pot)
            call GridTOBasis_nD_cplx(Oppsi%CVec, Oppsi_g%CVec, psi%Basis)
            call dealloc_psi(Oppsi_g)
         ELSE
            call init_psi(psi_g, Psi%Basis, .true., .true.)
            call init_psi(Oppsi_g, psi%Basis, .true., .true.)
            call BasisTOGrid_nD_cplx(psi_g%CVec, psi%CVec, psi%Basis)
            call Calc_Hpsi(psi_g%CVec, Oppsi_g%CVec, psi%Basis, Op%Scal_pot)
            call GridTOBasis_nD_cplx(Oppsi%CVec, Oppsi_g%CVec, psi%Basis)
            call dealloc_psi(psi_g)
            call dealloc_psi(Oppsi_g)
         END IF
      END IF
   END SUBROUTINE calc_OpPsi

   !===========================================================================
   ! SUBROUTINE Calc_tab_Iq - Calcul de la table d'indices
   !
   ! Paramètres :
   ! - Tab_Iq : table à calculer (inout)
   ! - Basis : base  (in)
   !===========================================================================
   SUBROUTINE Calc_tab_Iq(Tab_Iq, Basis)
      USE QDUtil_m
      USE Basis_m
      TYPE(Basis_t), intent(in), target :: Basis
      integer, allocatable, intent(inout) :: Tab_iq(:, :)
      
      ! Variables locales
      integer, allocatable :: Tab_iq0(:)
      integer :: ndim, iq, nq
      logical :: Endloop

      ! Initialisation
      ndim = size(Basis%tab_basis) - 1
      nq = Basis%nq

      ! Allocation et initialisation
      allocate(Tab_iq(ndim, nq), Tab_iq0(ndim))
      Call Init_tab_ind(Tab_iq0, Basis%NDindexq)
      
      ! Remplissage de la table
      Iq = 0
      DO
         Iq = Iq + 1
         CALL increase_NDindex(Tab_Iq0, Basis%NDindexq, Endloop)
         IF (Endloop) exit
         Tab_iq(:, Iq) = Tab_Iq0
      END DO
      
      ! Nettoyage
      deallocate(Tab_Iq0)
   END SUBROUTINE Calc_tab_Iq

   !===========================================================================
   ! SUBROUTINE Calc_Scalar_Pot_openmp - Calcul parallèle du potentiel scalaire
   !
   ! Paramètres :
   ! - V : potentiel à calculer (inout)
   ! - Basis : base (in)
   ! - Tab_Iq : table d'indices (in)
   !===========================================================================
   SUBROUTINE Calc_Scalar_Pot_openmp(V, Basis, Tab_Iq)
      !$ USE omp_lib
      USE QDUtil_m
      USE Basis_m
      TYPE(Basis_t), intent(in), target :: Basis
      real(kind=Rkind), allocatable, intent(inout) :: V(:, :, :)
      integer, intent(in) :: Tab_iq(:, :)
      
      ! Variables locales
      integer :: Ib, ndim, iq, nq, nsurf, maxth
      real(kind=Rkind), allocatable :: Q(:)
      logical :: Endloop

      ! Initialisation
      ndim = size(Basis%tab_basis) - 1
      nq = Basis%nq
      nsurf = Basis%tab_basis(ndim+1)%nb

      ! Allocation et initialisation
      allocate(Q(ndim))
      allocate(V(nq, nsurf, nsurf))
      V(:, :, :) = ZERO
      Q(:) = ZERO

      ! Détermination du nombre de threads
      maxth = 1
      !$ maxth = omp_get_max_threads()
      write(*,*) 'nbr de coeurs :', maxth
      
      ! Parallélisation avec OpenMP
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP SHARED(Basis, Tab_Iq, maxth, V, nq, ndim, Ib) &
      !$OMP PRIVATE(Iq, Q) &
      !$OMP NUM_THREADS(maxth)

      !$OMP BARRIER
      !$OMP DO SCHEDULE(STATIC)
      DO Iq = 1, nq
         ! Récupération des coordonnées
         DO Ib = 1, ndim
            Q(Ib) = Basis%tab_basis(Ib)%x(Tab_Iq(Ib, Iq))
         ENDDO

         ! Calcul du potentiel
         CALL sub_Qmodel_V(V(iq, :, :), Q(:))
      END DO
      !$OMP END DO
      !$OMP BARRIER
      !$OMP END PARALLEL

      ! Nettoyage
      deallocate(Q)
   END SUBROUTINE Calc_Scalar_Pot_openmp

   !===========================================================================
   ! SUBROUTINE test_openmp_op - Tests de parallélisation OpenMP
   !
   ! Paramètres :
   ! - Basis : base quantique (in)
   !===========================================================================
   SUBROUTINE test_openmp_op(Basis)
      USE QDUtil_m
      USE Basis_m
      TYPE(Basis_t), intent(in), target :: Basis
      
      ! Variables locales
      real(kind=Rkind), allocatable :: V(:, :, :)
      integer, allocatable :: Tab_Iq(:, :)
      real(kind=Rkind) :: t1, t2, tps, tpsopenmp
      integer :: iq, nq

      nq = Basis%nq
      print*, "nq=", nq
      
      ! Calcul de la table d'indices et du potentiel
      call Calc_tab_Iq(Tab_Iq, Basis)
      call Calc_Scalar_Pot(V, Basis)
   END SUBROUTINE test_openmp_op

   !===========================================================================
   ! SUBROUTINE calc_nac - Calcul des couplages non-adiabatiques
   !
   ! Paramètres :
   ! - Basis : base (in)
   !===========================================================================
   SUBROUTINE calc_nac(Basis)
      USE QDUtil_m
      USE Basis_m
      TYPE(Basis_t), intent(in), target :: Basis
      
      ! Variables locales
      real(kind=Rkind), allocatable :: NAC(:, :, :), G(:, :, :), V(:, :)
      integer, allocatable :: Tab_iq(:)
      integer :: inb, ndim, iq, nq, nsurf, i, nq1, nq2
      real(Kind=Rkind), allocatable :: Q(:)
      logical :: Endloop

      ! Initialisation
      ndim = size(Basis%tab_basis) - 1
      nq = Basis%nq
      nsurf = Basis%tab_basis(ndim+1)%nb
      nq1 = Basis%tab_basis(1)%nq
      nq2 = Basis%tab_basis(2)%nq

      ! Ouverture des fichiers de sortie
      open (unit = 300, file = "N121-ad=t.txt")
      open (unit = 301, file = "N122-ad=t.txt")
      open (unit = 302, file = "N211-ad=t.txt")
      open (unit = 303, file = "N212-ad=t.txt")
      open (unit = 500, file = "x.txt")
      open (unit = 501, file = "y.txt")

      ! Allocation et initialisation
      allocate(Q(ndim), Tab_iq(ndim))
      allocate(V(nsurf, nsurf), NAC(nsurf, nsurf, ndim), G(nsurf, nsurf, ndim))

      ! Calcul pour chaque point de la grille
      Call Init_tab_ind(Tab_iq, Basis%NDindexq)
      Iq = 0
      DO
         Iq = Iq + 1
         i = Tab_iq(2)
         CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop)
         IF (Endloop) exit
         
         ! Récupération des coordonnées
         do inb = 1, ndim
            Q(inb) = Basis%tab_basis(inb)%x(Tab_iq(inb))
         end do
         
         ! Calcul des couplages
         CALL sub_Qmodel_VG_NAC(V, G, NAC, Q)
         
         ! Gestion des sauts de ligne
         if ((Tab_iq(2)/= i)) then 
            write(300,*)
            write(301,*)
            write(302,*)
            write(302,*)
         end if  

         ! Écriture des résultats
         write(300,*) Q(:), NAC(1, 2, 1)
         write(301,*) Q(:), NAC(1, 2, 2)
         write(302,*) Q(:), NAC(2, 1, 1)
         write(303,*) Q(:), NAC(2, 1, 2)
      END DO
      
      ! Écriture des coordonnées
      do i = 1, nq1
         write(500,*) Basis%tab_basis(1)%x(i)
      end do
      do i = 1, nq2
         write(501,*) Basis%tab_basis(2)%x(i)
      end do

      ! Nettoyage
      deallocate(Tab_iq, Q, NAC, G, V)
   END SUBROUTINE calc_nac

   !===========================================================================
   ! SUBROUTINE calc_VV - Calcul du potentiel
   !
   ! Paramètres :
   ! - Basis : base quantique (in)
   !===========================================================================
   SUBROUTINE calc_VV(Basis)
      USE QDUtil_m
      USE Basis_m
      TYPE(Basis_t), intent(in), target :: Basis
      
      ! Variables locales
      integer, allocatable :: Tab_iq(:)
      integer :: inb, ndim, iq, nq, nsurf, i, nq1, nq2
      real(Kind=Rkind), allocatable :: Q(:), V(:, :)
      logical :: Endloop

      ! Initialisation
      ndim = size(Basis%tab_basis) - 1
      nq = Basis%nq
      nsurf = Basis%tab_basis(ndim+1)%nb
      nq1 = Basis%tab_basis(1)%nq
      nq2 = Basis%tab_basis(2)%nq

      ! Allocation et initialisation
      allocate(Q(ndim), Tab_iq(ndim))
      allocate(V(nsurf, nsurf))

      ! Calcul pour chaque point de la grille
      Call Init_tab_ind(Tab_iq, Basis%NDindexq)
      Iq = 0
      DO
         Iq = Iq + 1
         i = Tab_iq(2)
         CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop)
         IF (Endloop) exit
         
         ! Récupération des coordonnées
         do inb = 1, ndim
            Q(inb) = Basis%tab_basis(inb)%x(Tab_iq(inb))
         end do
         
         ! Calcul du potentiel
         CALL sub_Qmodel_V(V, Q)
      END DO

      ! Nettoyage
      deallocate(Tab_iq, Q, V)
   END SUBROUTINE calc_VV

end module Op_m