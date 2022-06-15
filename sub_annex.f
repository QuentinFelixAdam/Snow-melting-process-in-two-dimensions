SUBROUTINE lire()
USE m_global
IMPLICIT NONE

! FILES FOR CASE

nlinescase = 11*10000

ALLOCATE(V_air(nlinescase))

ALLOCATE(f_s(nlinescase))

ALLOCATE(T_air(nlinescase))

ALLOCATE(fD(nlinescase))


END SUBROUTINE lire

