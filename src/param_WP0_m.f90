
      MODULE param_WP0_m
      USE UtilLib_m
      IMPLICIT NONE

        TYPE GWP1D_t

        real (kind=Rkind) :: sigma = ONETENTH  ! width of WP0
        real (kind=Rkind) :: Q0    = ZERO      ! position of WP0
        real (kind=Rkind) :: imp_k = ZERO      ! impultion for WP0
        real (kind=Rkind) :: phase = ZERO      ! phase for WP0

        END TYPE GWP1D_t

        TYPE GWP_t
          complex (kind=Rkind)        :: Coef           = CZERO

          TYPE (GWP1D_t), allocatable :: tab_GWP1D(:)

        END TYPE GWP_t

        TYPE param_WP0


        integer             :: nb_WP0             = 1       ! default: 1
        !for each variable Qi : exp[-((Q-Qeq)/sigma)2+i*imp_k*(Q-Qeq)]
        real (kind=Rkind), allocatable :: WP0sigma(:)       ! WP0sigma(nb_act1) : sigma for WP0
        real (kind=Rkind), allocatable :: WP0Qeq(:)         ! WP0Qeq(nb_act1)   : position of WP0
        real (kind=Rkind), allocatable :: WP0imp_k(:)       ! WP0imp_k(nb_act1) : impultion for WP0
        real (kind=Rkind), allocatable :: WP0phase(:)       ! WP0imp_k(nb_act1) : phase for WP0
        TYPE (GWP_t), allocatable :: tab_GWP0(:)

        END TYPE param_WP0

        CONTAINS

  SUBROUTINE Read_GWP1D(GWP1D)
    TYPE (GWP1D_t), intent(inout) :: GWP1D

    !------ initial WP definition -----------------------------
    !     GWP(Q)=exp[-((Q-Qeq)/sigma)2+i*imp_k*(Q-Qeq)+i*phase]
    real (kind=Rkind) :: sigma,imp_k,Qeq,phase
    integer           :: Rerr

    NAMELIST /defWP0/sigma,Qeq,imp_k,phase

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Read_GWP1D'
    logical, parameter :: debug =.FALSE.
    !logical, parameter :: debug =.TRUE.
    !-----------------------------------------------------------


    sigma    = ONETENTH
    Qeq      = ZERO
    imp_k    = ZERO
    phase    = ZERO

    read(in_unitp,defWP0,iostat=Rerr)
    IF (Rerr /= 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' problem while reading the namelist "defWP0"'
      write(out_unitp,defWP0)
      STOP 'ERROR in Read_GWP1D: problem while reading the namelist "defWP0"'
    END IF

    GWP1D = GWP1D_t(sigma,Qeq,imp_k,phase)

  END SUBROUTINE Read_GWP1D

  SUBROUTINE Write_GWP1D(GWP1D)
    TYPE (GWP1D_t), intent(in) :: GWP1D

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Write_GWP1D'
    logical, parameter :: debug =.FALSE.
    !logical, parameter :: debug =.TRUE.
    !-----------------------------------------------------------

    IF (debug) write(out_unitp,*) 'Q0,sigma,imp_k,phase',           &
           GWP1D%Q0,GWP1D%sigma,GWP1D%imp_k,GWP1D%phase

    write(out_unitp,'(2i6,4f12.6)')                                             &
           GWP1D%iQact,GWP1D%iQdyn,GWP1D%Q0,GWP1D%sigma,GWP1D%imp_k,GWP1D%phase

  END SUBROUTINE Write_GWP1D
  FUNCTION calc_GWP1D(GWP1D,Q)
    complex (kind=Rkind)          :: calc_GWP1D
    real (kind=Rkind), intent(in) :: Q
    TYPE (GWP1D_t),    intent(in) :: GWP1D

    real (kind=Rkind) :: ze,zk,DQ
    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='calc_GWP1D'
    logical, parameter :: debug =.FALSE.
    !logical, parameter :: debug =.TRUE.
    !-----------------------------------------------------------

    DQ = Q-GWP1D%Q0
    ze = (DQ/GWP1D%sigma)**2
    zk = mod(DQ*GWP1D%imp_k + GWP1D%phase,TWO*pi)

    calc_GWP1D = exp(-ze + EYE*zk) / sqrt(sqrt(pi/TWO)*GWP1D%sigma)

  END FUNCTION calc_GWP1D

  SUBROUTINE Read_GWP(GWP,ndim,nb_bi,nb_elec)
    TYPE (GWP_t), intent(inout) :: GWP
    integer,      intent(in)    :: ndim,nb_elec,nb_bi


    !------ initial WP definition -----------------------------
    integer                     :: i,Rerr
    integer                     :: I_ElecChannel
    integer                     :: I_HAChannel
    integer                     :: tab_iQact(ndim)
    integer                     :: tab_iQdyn(ndim)
    integer                     :: nb_Qactm1,nb_Qdynm1

    complex (kind=Rkind)        :: Coef

    NAMELIST /defGWP/ Coef,I_ElecChannel,I_HAChannel


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Read_GWP'
    logical, parameter :: debug =.FALSE.
    !logical, parameter :: debug =.TRUE.
    !-----------------------------------------------------------

    IF (print_level > 0)                                                        &
       write(out_unitp,*) 'WP0(Q)=exp[-((Q-Q0)/sigma)^2+i*imp_k*(Q-Q0)+i*phase]'

    I_ElecChannel  = -1
    I_HAChannel    = -1
    Coef           = CZERO

    read(in_unitp,defGWP,iostat=Rerr)
    IF (Rerr /= 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' problem while reading the namelist "defGWP"'
      write(out_unitp,defGWP)
      STOP 'ERROR in Read_GWP: problem while reading the namelist "defGWP"'
    END IF

    IF (I_ElecChannel  == -1) I_ElecChannel  = 1
    IF (I_ElecChannel < 1 .OR. I_ElecChannel > nb_elec) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '  I_ElecChannel is out-of-range.'
      write(out_unitp,*) ' I_ElecChannel',I_ElecChannel
      write(out_unitp,*) ' The number of electronic states is: ',nb_elec
      STOP 'ERROR in Read_GWP: I_ElecChannel is out-of-range'
    END IF

    IF (I_HAChannel    == -1) I_HAChannel    = 1
    IF (I_HAChannel < 1 .OR. I_HAChannel > nb_bi) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '  I_HAChannel is out-of-range.'
      write(out_unitp,*) ' I_HAChannel',I_HAChannel
      write(out_unitp,*) ' The number of adiabatic (HAC) channel is: ',nb_bi
      STOP 'ERROR in Read_GWP: I_HAChannel is out-of-range'
    END IF

    GWP%I_ElecChannel = I_ElecChannel
    GWP%I_HAChannel   = I_HAChannel
    GWP%Coef          = Coef

    allocate(GWP%tab_GWP1D(ndim))
    DO i=1,ndim
      CALL Read_GWP1D(GWP%tab_GWP1D(i))
    END DO

    tab_iQact(:) = GWP%tab_GWP1D(:)%iQact
    tab_iQdyn(:) = GWP%tab_GWP1D(:)%iQdyn

    nb_Qactm1 = count(tab_iQact == -1)
    nb_Qdynm1 = count(tab_iQdyn == -1)

    IF ((nb_Qactm1 > 0 .AND. nb_Qactm1 < ndim) .OR. (nb_Qdynm1 > 0 .AND. nb_Qdynm1 < ndim)) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '  All or none iQact must be defined'
      write(out_unitp,*) '    or'
      write(out_unitp,*) '  all or none iQdyn must be defined'
      write(out_unitp,*) ' tab_iQact',tab_iQact
      write(out_unitp,*) ' tab_iQdyn',tab_iQdyn
      STOP 'ERROR in Read_GWP: wrong iQact and/or iQdyn'
    END IF
    IF (nb_Qactm1 == 0 .AND. nb_Qdynm1 == 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '  All iQact and all iQact are defined'
      write(out_unitp,*) '  You MUST define only iQact or iQact'
      write(out_unitp,*) ' tab_iQact',tab_iQact
      write(out_unitp,*) ' tab_iQdyn',tab_iQdyn
      STOP 'ERROR in Read_GWP: wrong iQact and/or iQdyn'
    END IF

    IF ( nb_Qactm1 == ndim .AND. nb_Qdynm1 == ndim) THEN
      GWP%tab_GWP1D(:)%iQact = [(i,i=1,ndim)]
    END IF

  END SUBROUTINE Read_GWP
  FUNCTION calc_GWP(GWP,Q)
    complex (kind=Rkind)          :: calc_GWP
    real (kind=Rkind), intent(in) :: Q(:)
    TYPE (GWP_t),      intent(in) :: GWP

    integer :: i
    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='calc_GWP'
    logical, parameter :: debug =.FALSE.
    !logical, parameter :: debug =.TRUE.
    !-----------------------------------------------------------

    calc_GWP = CONE
    DO i=1,size(GWP%tab_GWP1D)
      calc_GWP = calc_GWP * calc_GWP1D(GWP%tab_GWP1D(i),Q(i))
    END DO
    calc_GWP = calc_GWP * GWP%Coef

  END FUNCTION calc_GWP
  SUBROUTINE Write_GWP(GWP)
    TYPE (GWP_t), intent(in) :: GWP

    integer :: i
    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Write_GWP'
    logical, parameter :: debug =.FALSE.
    !logical, parameter :: debug =.TRUE.
    !-----------------------------------------------------------

    write(out_unitp,*) 'I_HAChannel,I_ElecChannel,Coef',GWP%I_HAChannel,GWP%I_ElecChannel,GWP%Coef
    write(out_unitp,'(a)') ' iQact iQdyn Qeq         sigma       imp_k       phase'
    DO i=1,size(GWP%tab_GWP1D)
      CALL Write_GWP1D(GWP%tab_GWP1D(i))
    END DO

  END SUBROUTINE Write_GWP

  SUBROUTINE Read_tab_GWP(tab_GWP,nb_GWP,mole,nb_bi,nb_elec)
    USE mod_Coord_KEO

    TYPE (GWP_t), allocatable, intent(inout) :: tab_GWP(:)
    integer,                   intent(in)    :: nb_GWP,nb_elec,nb_bi

    TYPE (CoordType),          intent(in)    :: mole


    !------ initial WP definition -----------------------------
      integer                :: iGWP,iQdyn,i

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Read_tab_GWP'
    logical, parameter :: debug =.FALSE.
    !logical, parameter :: debug =.TRUE.
    !-----------------------------------------------------------

    allocate(tab_GWP(nb_GWP))
    DO i=1,nb_GWP
      CALL Read_GWP(tab_GWP(i),mole%nb_act1,nb_bi,nb_elec)
    END DO

    ! check if iQact is defined
    DO iGWP=1,size(tab_GWP)
      IF (count(tab_GWP(iGWP)%tab_GWP1D(:)%iQact == -1) /= 0) THEN ! iQact is not defined => transfert iQdyn to iQact
        DO i=1,mole%nb_act1
          iQdyn = tab_GWP(iGWP)%tab_GWP1D(i)%iQdyn
          tab_GWP(iGWP)%tab_GWP1D(i)%iQact = mole%liste_QdynTOQact(iQdyn)
        END DO
      END IF
    END DO

    CALL Write_Tab_GWP(tab_GWP)


  END SUBROUTINE Read_tab_GWP
  SUBROUTINE Read_tab_GWP_v0(tab_GWP,nb_GWP,ndim,nb_bi,nb_elec)
    TYPE (GWP_t), allocatable, intent(inout) :: tab_GWP(:)
    integer,                   intent(in)    :: nb_GWP,ndim,nb_elec,nb_bi


    !------ initial WP definition -----------------------------
    integer                     :: i

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Read_tab_GWP_v0'
    logical, parameter :: debug =.FALSE.
    !logical, parameter :: debug =.TRUE.
    !-----------------------------------------------------------

    allocate(tab_GWP(nb_GWP))
    DO i=1,nb_GWP
      CALL Read_GWP(tab_GWP(i),ndim,nb_bi,nb_elec)
    END DO

  END SUBROUTINE Read_tab_GWP_v0
  SUBROUTINE Write_Tab_GWP(tab_GWP)
    TYPE (GWP_t), intent(in) :: tab_GWP(:)

    integer :: i
    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Write_Tab_GWP'
    logical, parameter :: debug =.FALSE.
    !logical, parameter :: debug =.TRUE.
    !-----------------------------------------------------------

    DO i=1,size(tab_GWP)
      CALL Write_GWP(tab_GWP(i))
    END DO

  END SUBROUTINE Write_Tab_GWP
!================================================================
!
!    alloc / dealloc param_WP0
!
!================================================================
      !!@description: alloc param_WP0
      !!@param: para_poly
      SUBROUTINE alloc_param_WP0(para_WP0,WP0Grid_Gaussian,WP0_CleanChannel)
      IMPLICIT NONE
      TYPE (param_WP0), intent(inout) :: para_WP0
      logical, intent(in) :: WP0Grid_Gaussian,WP0_CleanChannel

        IF (para_WP0%nb_act1 > 0 .AND. WP0Grid_Gaussian) THEN
          IF (.NOT. allocated(para_WP0%WP0sigma)) THEN
            CALL alloc_NParray(para_WP0%WP0sigma,[para_WP0%nb_act1],    &
                              "para_WP0%WP0sigma","alloc_param_WP0")
          END IF

          IF (.NOT. allocated(para_WP0%WP0Qeq)) THEN
            CALL alloc_NParray(para_WP0%WP0Qeq,[para_WP0%nb_act1],      &
                              "para_WP0%WP0Qeq","alloc_param_WP0")
          END IF

          IF (.NOT. allocated(para_WP0%WP0imp_k)) THEN
            CALL alloc_NParray(para_WP0%WP0imp_k,[para_WP0%nb_act1],    &
                              "para_WP0%WP0imp_k","alloc_param_WP0")
          END IF

          IF (.NOT. allocated(para_WP0%WP0phase)) THEN
            CALL alloc_NParray(para_WP0%WP0phase,[para_WP0%nb_act1],    &
                              "para_WP0%WP0phase","alloc_param_WP0")
          END IF

        END IF

        IF (WP0_CleanChannel .AND.                                      &
            .NOT. allocated(para_WP0%WP0_CleanChannellist)) THEN
          CALL alloc_NParray(para_WP0%WP0_CleanChannellist,               &
                                     [para_WP0%WP0_nb_CleanChannel],  &
                          "para_WP0%WP0_CleanChannellist","alloc_param_WP0")
        END IF

      END SUBROUTINE alloc_param_WP0

      !!@description: dealloc param_WP0
      !!@param: para_poly
      SUBROUTINE dealloc_param_WP0(para_WP0)
      IMPLICIT NONE
      TYPE (param_WP0), intent(inout) :: para_WP0

        IF (allocated(para_WP0%WP0sigma)) THEN
          CALL dealloc_NParray(para_WP0%WP0sigma,                         &
                              "para_WP0%WP0sigma","dealloc_param_WP0")
        END IF

        IF (allocated(para_WP0%WP0Qeq)) THEN
          CALL dealloc_NParray(para_WP0%WP0Qeq,                           &
                              "para_WP0%WP0Qeq","dealloc_param_WP0")
        END IF

        IF (allocated(para_WP0%WP0imp_k)) THEN
          CALL dealloc_NParray(para_WP0%WP0imp_k,                         &
                              "para_WP0%WP0imp_k","dealloc_param_WP0")
        END IF

        IF (allocated(para_WP0%WP0phase)) THEN
          CALL dealloc_NParray(para_WP0%WP0phase,                         &
                              "para_WP0%WP0phase","dealloc_param_WP0")
        END IF

        IF (allocated(para_WP0%WP0_CleanChannellist)) THEN
          CALL dealloc_NParray(para_WP0%WP0_CleanChannellist,             &
                            "para_WP0%WP0_CleanChannellist","dealloc_param_WP0")
        END IF
      END SUBROUTINE dealloc_param_WP0

      END MODULE param_WP0_m
