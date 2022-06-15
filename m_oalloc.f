MODULE oalloc

IMPLICIT NONE

INTEGER::info1
INTEGER,PARAMETER :: MKS = KIND(1.0E0)
INTEGER,PARAMETER :: MKD = KIND(1.0D0)
INTEGER,PARAMETER :: MK = MKS
REAL(MKS),DIMENSION(:),ALLOCATABLE :: swork
REAL(MKD),DIMENSION(:),ALLOCATABLE :: dwork
REAL(MK),DIMENSION(:),ALLOCATABLE :: Mold,Mnew

  INTERFACE alloc
      MODULE PROCEDURE salloc,dalloc
  END INTERFACE alloc

  CONTAINS
      SUBROUTINE salloc(array,Dimz)
      INTEGER :: Dimz    
      REAL(MKS),DIMENSION(:),ALLOCATABLE :: array

      IF(ALLOCATED(array))THEN
         swork = array
      ELSE
         ALLOCATE(array(Dimz),STAT=info1)
         !PRINT*,info1
      ENDIF
      END SUBROUTINE salloc

      SUBROUTINE dalloc(array,Dimz)
      INTEGER :: Dimz  
      REAL(MKD),DIMENSION(:),ALLOCATABLE :: array

      IF(ALLOCATED(array))THEN
         dwork = array
      ELSE
         ALLOCATE(array(Dimz),STAT=info1)
         !PRINT*,info1
      ENDIF
      END SUBROUTINE dalloc

END MODULE oalloc
