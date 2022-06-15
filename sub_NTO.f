!Place the copying of the new-to-old data field in a subroutine

ELEMENTAL SUBROUTINE NTO(M1,M2)
USE m_global
IMPLICIT NONE
REAL,INTENT(INOUT)::M1,M2

M1 = M2

END SUBROUTINE NTO
