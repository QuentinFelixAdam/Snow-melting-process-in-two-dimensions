SUBROUTINE lire()
USE m_global
IMPLICIT NONE

PRINT*,'-----------------------------------'
PRINT*,'------SIMULATION HAS STARTED-------'
PRINT*,'-----------------------------------'
PRINT*,'-----------LOADING DATA------------'

! FILES FOR CASE

nlinescase = 0

OPEN(36,FILE='Windspeed.dat',STATUS='old')

DO
   READ(36,*,END=22)
   nlinescase = nlinescase + 1
END DO
22 CLOSE(36)

ALLOCATE(V_air(nlinescase))

OPEN(36,FILE='Windspeed.dat',STATUS='old')

DO i = 1,nlinescase
   READ(36,*) V_air(i)
END DO

CLOSE(36)

PRINT*,'...'

nlinescase = 0

OPEN(36,FILE='Shortwave.dat',STATUS='old')

DO
   READ(36,*,END=23)
   nlinescase = nlinescase + 1
END DO
23 CLOSE(36)

ALLOCATE(f_s(nlinescase))

OPEN(36,FILE='Shortwave.dat',STATUS='old')

DO i = 1,nlinescase
   READ(36,*) f_s(i)
END DO

CLOSE(36)


PRINT*,'...'

nlinescase = 0

OPEN(36,FILE='Airtemp.dat',STATUS='old')

DO
   READ(36,*,END=24)
   nlinescase = nlinescase + 1
END DO
24 CLOSE(36)

ALLOCATE(T_air(nlinescase))

OPEN(36,FILE='Airtemp.dat',STATUS='old')

DO i = 1,nlinescase
   READ(36,*) T_air(i)
END DO

CLOSE(36)

PRINT*,'...'

nlinescase = 0

OPEN(36,FILE='Longwave.dat',STATUS='old')

DO
   READ(36,*,END=25)
   nlinescase = nlinescase + 1
END DO
25 CLOSE(36)

ALLOCATE(fD(nlinescase))

OPEN(36,FILE='Longwave.dat',STATUS='old')

DO i = 1,nlinescase
   READ(36,*) fD(i)
END DO

CLOSE(36)

PRINT*,'...'
PRINT*,'-----------DATA LOADED-------------'

END SUBROUTINE lire

