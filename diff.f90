PROGRAM diff_WP
  implicit none
  
  integer, parameter :: Rkind=8
  character (len=:), allocatable :: nameA, nameB, diff_file
  
  
  
  integer :: nioA,nioB,niodiff, it,ib1,ib2,idum
  real(kind=Rkind) :: t,normdiff,maxnormdiff
  complex(kind=Rkind) :: psiA,psiB
  
  integer :: nbA,nbB, maxit,nbbA,nbbB
  logical :: stA,stB
  
  namelist /res / nbA,nbB,stA,stB,maxit
  
  
  
  !----------------------------------------------------------------------
  nbA   = 20
  nbB   = 30
  maxit = 600
  stA   = .TRUE.
  stB   = .TRUE.
  read(*,res)
  
  nbbA = nbA**2
  nbbB = nbB**2
  
  IF (stA) THEN 
    nameA = 'results_std_nb' // int_TO_string(nbbA) // '/psi_dt_on_basis0_non_hagedorn_taylor.txt'
  ELSE
    nameA = 'results_H_nb' // int_TO_string(nbbA) // '/psi_dt_on_basis0_hagedorn_taylor.txt'
  END IF
  IF (stB) THEN 
    nameB = 'results_std_nb' // int_TO_string(nbbB) // '/psi_dt_on_basis0_non_hagedorn_taylor.txt'
  ELSE
    nameB = 'results_H_nb' // int_TO_string(nbbB) // '/psi_dt_on_basis0_hagedorn_taylor.txt'
  END IF
  
  IF (stA) THEN 
    diff_file = 'diff_std-nb' // int_TO_string(nbbA)
  ELSE
    diff_file = 'diff_H-nb' // int_TO_string(nbbA)
  END IF
  IF (stB) THEN 
    diff_file = diff_file // '_std-nb' // int_TO_string(nbbB)
  ELSE
    diff_file = diff_file // '_H-nb' // int_TO_string(nbbB)
  END IF
  write(*,*) 'nameA:     ',nameA
  write(*,*) 'nameB:     ',nameB
  write(*,*) 'diff_file: ',diff_file
  !----------------------------------------------------------------------
  
  !----------------------------------------------------------------------
  open(newunit=nioA,file=nameA)
  open(newunit=nioB,file=nameB)
  open(newunit=niodiff,file=diff_file)
  
  maxnormdiff = 0._Rkind
  DO it=0,maxit
    normdiff = 0._Rkind
    print*,'it', it ; flush(6)
    DO ib1=1,max(nbA,nbB)
    DO ib2=1,max(nbA,nbB)
      psiA = 0._Rkind
      psiB = 0._Rkind
      !print*,'ib1,ib2, t',ib1,ib2, it ; flush(6)
      IF (ib1 <= nbA .AND. ib2 <= nbA)then
        !print*, 'lecture du fichier A' ; flush(6)
         read(nioA,*) t,idum,psiA
      end if
      IF (ib1 <= nbB .AND. ib2 <= nbB) then 
        !print*, 'lecture du fichier B' ; flush(6)
        read(nioB,*) t,idum,psiB
      end if
  
      normdiff = normdiff + abs(psiA-psiB)**2
    END DO
    END DO
    read(nioA,*)
    read(nioB,*)
  
    normdiff = sqrt(normdiff)
    maxnormdiff = max(maxnormdiff,normdiff)
    write(niodiff,*) t,normdiff
  END DO
  write(*,*) 'largest ',diff_file,' :',maxnormdiff
  
  close(niodiff)
  close(nioA)
  close(nioB)
  !----------------------------------------------------------------------
  
  CONTAINS
  FUNCTION int_TO_string(i) RESULT(string)
  IMPLICIT NONE
  
  character (len=:),  allocatable             :: string
  integer,                         intent(in) :: i
  
  
  character (len=:), allocatable  :: name_int
  integer :: clen
  
  ! first approximated size of name_int
  IF (i == 0) THEN
    clen = 1
  ELSE IF (i < 0) THEN
    clen = int(log10(abs(real(i,kind=8))))+2
  ELSE
    clen = int(log10(real(i,kind=8)))+1
  END IF
  
  ! allocate name_int
  allocate(character(len=clen) :: name_int)
  
  ! write i in name_int
  write(name_int,'(i0)') i
  
  ! transfert name_int in QDUtil_int_TO_char
  string = trim(adjustl(name_int))
  
  ! deallocate name_int
  deallocate(name_int)
  
  END FUNCTION int_TO_string
  
  END PROGRAM diff_WP
  