!PLOT Results | FILE OUTPUT

SUBROUTINE output()
USE m_global
USE oalloc
IMPLICIT NONE

OPEN(15,FILE='diff.dat')
DO j=1,Nz
    WRITE(15,'(3E12.4)')REAL(j-1)*dz,Mnew(j)
    WRITE(15,'(A)')
ENDDO
CLOSE(15)
END SUBROUTINE output
