MODULE NDindex_m
  USE QDUtil_m, out_unit => out_unit
  IMPLICIT NONE

   TYPE :: NDindex_t
    Character(len=:),allocatable  :: sys_type
    integer                       :: Ndim  = 0
    integer                       :: Nterm = 0
    integer                       :: Nterm_compact = 0
    integer                       :: L        = 0
    integer, allocatable          :: Li_smol(:)
    integer, allocatable          :: Tab0(:)      ! = [0,1]
    integer, allocatable          :: NDsize(:)
    integer, allocatable          :: NDend(:)
    Integer, allocatable          :: NDinit(:)
  END TYPE NDindex_t


  PRIVATE
  PUBLIC :: increase_NDindex,Init_tab_ind,Testindex,NDindex_t,Init_NDindex,dealloc_NDindex!,alloc_Tab

CONTAINS

 SUBROUTINE Init_NDindex(NDindex,NDend,Ndim)
    TYPE(NDindex_t),intent(inout)   :: NDindex
    Integer ,      intent(in)       :: NDend(:)
    Integer ,      intent(in)       :: Ndim

    !Logical,    parameter          :: debug = .true.
    Logical,     parameter          :: debug = .false.

    IF (debug) THEN
       write(out_unit,*) 'BEGINNING Init_NDindex'
      flush(out_unit)
    END IF

    NDindex%Ndim  = Ndim
    NDindex%NDend = NDend

    Allocate(NDindex%Tab0(NDindex%Ndim))
    SELECT CASE (NDindex%sys_type)
    CASE ('smolyak')
      NDindex%Tab0(:)   =  0
      NDindex%Tab0(1)   = -1
    CASE ('dp')
      NDindex%Tab0(:)   =  1
      NDindex%Tab0(1)   =  0
    CASE default
    STOP 'ERROR in Read_Basis: no default sym_type.'
    END SELECT

    IF (debug) THEN
      write(out_unit,*)  NDindex%Ndim
      write(out_unit,*) 'NDindex%Ndend', NDindex%Ndend
      write(out_unit,*) 'NDindex%Tab0', NDindex%Tab0(:)
      write(out_unit,*)  'END Init_NDindex'
      flush(out_unit)
    END IF

  END SUBROUTINE Init_NDindex

  SUBROUTINE Init_tab_ind(Tab_ind,NDindex)
    USE QDUtil_m
    IMPLICIT NONE
    TYPE(NDindex_t),  intent(in) :: NDindex
    integer,     intent(inout)   :: Tab_ind(:)

    !logical,    parameter       :: debug = .true.
    logical,     parameter       :: debug = .false.

    IF (debug) THEN
      write(out_unit,*) 'BEGINNING Init_tab_ind'
      flush(out_unit)
    END IF

     Tab_ind(:)= NDindex%Tab0(:)

    IF (debug) THEN
      write(out_unit,*) 'END Init_tab_ind'
      flush(out_unit)
    END IF

  END SUBROUTINE Init_tab_ind


  SUBROUTINE increase_NDindex(Tab_ind,NDindex,Endloop)
    USE QDUtil_m
    IMPLICIT NONE
    TYPE(NDindex_t), intent(in)     :: NDindex
    Integer, intent(inout)          :: Tab_ind(:)
    !Logical, parameter             :: debug = .true.
    Logical, parameter              :: debug = .false.
    Logical, intent(inout)          :: Endloop
    character(len=Name_len)         :: name
    Integer                         :: i


    IF (debug) THEN
      Write(out_unit,*)'BEGINNING Tab_ind'
      Write(out_unit,*)'NDindex%Ndend', NDindex%Ndend
      Write(out_unit,*)'Tab_ind1', Tab_ind
      flush(out_unit)
    END IF

    IF (debug)  Write(out_unit,*)'Tab_ind', Tab_ind

      name=NDindex%sys_type
      CALL string_uppercase_TO_lowercase(name)
    Tab_ind(1)=Tab_ind(1)+1

    SELECT CASE (NDindex%sys_type)

    CASE ('smolyak')
     DO i=1,NDindex%Ndim-1
       IF(Sum(Tab_ind)>NDindex%L .or.Tab_ind(i) >NDindex%Ndend(i) ) THEN
          Tab_ind(i+1)  = Tab_ind(i+1)+1
          Tab_ind(i)    = 0
       END IF
     END DO
     Endloop = (Tab_ind(NDindex%Ndim) == NDindex%Ndend(NDindex%Ndim))
    CASE ('dp')
     DO i=1,NDindex%Ndim-1
       IF(Tab_ind(i) > NDindex%Ndend(i)) THEN
          Tab_ind(i+1) = Tab_ind(i+1)+1
          Tab_ind(i)   = 1
       END IF
     END DO
     IF (debug) Write(out_unit,*)'Tab_indfin', Tab_ind(NDindex%Ndim)
     Endloop = (Tab_ind(NDindex%Ndim) == NDindex%Ndend(NDindex%Ndim)+1)
    CASE default
      STOP 'ERROR in Read_Basis: no default sym_type.'
    END SELECT

    IF (debug) THEN
      Write(out_unit,*) 'END Tab_ind'
      flush(out_unit)
    END IF

  END SUBROUTINE increase_NDindex

    SUBROUTINE Write_NDindex(NDindex)
      USE QDUtil_m
      IMPLICIT NONE

      TYPE(NDindex_t),  intent(in) :: NDindex
      !logical,    parameter      :: debug = .true.
      logical,     parameter      :: debug = .false.

      IF (debug) THEN
        write(out_unit,*) 'BEGINNING Write_NDindex'
        flush(out_unit)
      END IF

      write(out_unit,*) 'NDindex%Ndim =',NDindex%Ndim
      write(out_unit,*)'NDindex%Ndend =', NDindex%Ndend
      write(out_unit,*)'NDindex%Tab0 =', NDindex%Tab0(:)

      IF (debug) THEN
        write(out_unit,*) ' END Write_NDindex'
        flush(out_unit)
      END IF

    END SUBROUTINE Write_NDindex


  SUBROUTINE dealloc_NDindex(NDindex)
   USE QDUtil_m
   IMPLICIT NONE
   TYPE(NDindex_t),  intent(inout) :: NDindex
   !logical,    parameter      :: debug = .true.
   logical,     parameter      :: debug = .false.
   IF (debug) THEN
     write(out_unit,*) 'BEGINNING dealloc_NDindex'
     flush(out_unit)
   END IF
   If(allocated(NDindex%Ndend)) deallocate(NDindex%Ndend)
   If(allocated(NDindex%Tab0)) deallocate(NDindex%Tab0)
 
   IF (debug) THEN
     write(out_unit,*) ' END dealloc_NDindex'
     flush(out_unit)
   END IF
 END SUBROUTINE

  SUBROUTINE Testindex(NDindex)
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(NDindex_t),  intent(inout) :: NDindex
    integer, allocatable         :: Tab_ind(:)
    logical                      :: Endloop
  !  integer,     intent(in)      :: nl
    integer                      :: i
    !logical,    parameter      :: debug = .true.
    logical,     parameter      :: debug = .false.

    IF (debug) THEN
      write(out_unit,*) 'BEGINNING Testindex'
      flush(out_unit)
    END IF
    Allocate(Tab_ind(NDindex%Ndim))
   ! write(out_unit,*) 'NDindex%NDim',NDindex%Ndim
    ! write(out_unit,*) 'NDindex%NDend',NDindex%NDend
    !STOP 'cc'
    call Init_NDindex(NDindex, NDindex%NDend, NDindex%Ndim)
    Call Init_tab_ind(Tab_ind,NDindex)

    i=0
   ! write(out_unit,*) i, tab_ind(:)
    !STOP 'cc'
    DO !i= 1,100
      i=i+1
      CALL increase_NDindex(Tab_ind,NDindex,Endloop)
       write(out_unit,*) i, tab_ind(:),sum(tab_ind(:))
        IF (Endloop) exit
     
    END DO

    IF (debug) THEN
      write(out_unit,*) 'END Testindex'
      flush(out_unit)
    END IF

  END SUBROUTINE Testindex
END MODULE NDindex_m
