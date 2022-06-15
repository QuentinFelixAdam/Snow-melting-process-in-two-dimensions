!Define boundaries | INITIALIZATION

SUBROUTINE init()
USE m_global
USE oalloc
IMPLICIT NONE

alpha = k/(rho*c)
alpha1 = k1/(rho1*c1)
alpha2 = k2/(rho2*c2)
alpha_snow = k_snow/(rho_snow*c_snow)

indice_R = INT(z_R*Nz)

avgrho1 = 2.0*rho*rho1/(rho+rho1)
avgrho2 = 2.0*rho1*rho2/(rho1+rho2)

dz = L/(Nz-1.0)
dy = Ly/(Ny-1.0)

indexsnowinit = INT(REAL(depthsnowinit)/REAL(dz))

IF(indexsnowinit.EQ.0) THEN

allocationdepth = Nz

ALLOCATE(Told(allocationdepth,Ny))
ALLOCATE(Tnew(allocationdepth,Ny))
ALLOCATE(c_eff(allocationdepth,Ny))
ALLOCATE(k_eff(allocationdepth,Ny))
ALLOCATE(rho_eff(allocationdepth,Ny))
ALLOCATE(d_T_theta(allocationdepth,Ny))
ALLOCATE(theta(allocationdepth,Ny))
ALLOCATE(Tplot(allocationdepth,Ny))
ALLOCATE(Taux_value(allocationdepth,Ny))
ALLOCATE(Taux_time(allocationdepth,Ny))
ALLOCATE(Taux_value_right(allocationdepth,Ny))
ALLOCATE(Taux_time_right(allocationdepth,Ny))
ALLOCATE(Taux_new(allocationdepth,Ny))
ALLOCATE(Taux_old(allocationdepth,Ny))

ELSE

allocationdepth = Nz + indexsnowinit + 2

ALLOCATE(Told(allocationdepth,Ny))
ALLOCATE(Tnew(allocationdepth,Ny))
ALLOCATE(c_eff(allocationdepth,Ny))
ALLOCATE(k_eff(allocationdepth,Ny))
ALLOCATE(d_T_theta(allocationdepth,Ny))
ALLOCATE(theta(allocationdepth,Ny))
ALLOCATE(Tplot(allocationdepth,Ny))
ALLOCATE(Taux_value(allocationdepth,Ny))
ALLOCATE(Taux_time(allocationdepth,Ny))
ALLOCATE(Taux_value_right(allocationdepth,Ny))
ALLOCATE(Taux_time_right(allocationdepth,Ny))
ALLOCATE(Taux_new(allocationdepth,Ny))
ALLOCATE(Taux_old(allocationdepth,Ny))

END IF

ALLOCATE(depthsnow(Ny))
ALLOCATE(snowhelp(Ny))
ALLOCATE(indexsnowtot(Ny))
ALLOCATE(indexmeltup(Ny))
ALLOCATE(indexmeltdown(Ny))
ALLOCATE(indexsnowdead(Ny))
ALLOCATE(indexmeltuptot(Ny))
ALLOCATE(indexmeltdowntot(Ny))
ALLOCATE(indexmeltedold(Ny))
ALLOCATE(cumul_melted(Ny))
ALLOCATE(cumul_vertical(Ny))





DO m = 1,Ny
   depthsnow(m) = indexsnowinit
END DO

DO i = indexsnowinit,allocationdepth
DO j = 1,Ny

 Told(i,j) = 5.3*(i-indexsnowinit)*dz-5.8
 c_eff(i,j) = c_snow_f

END DO
END DO

DO i = indexsnowinit,allocationdepth
DO j = 1,Ny

 Told(i,j) = 5.0*(i-indexsnowinit)*dz-5.0
 c_eff(i,j) = c_snow_f

END DO
END DO

DO i = 1,indexsnowinit
DO j = 1,Ny

 Told(i,j) = -5.0
 c_eff(i,j) = c_snow_f

END DO
END DO

Tnew = Told

ALLOCATE(w_R(allocationdepth,Ny))

DO i = 1,allocationdepth
DO j = 1,Ny

w_R(i,j) = 0.0

END DO
END DO

DO j = 1,Ny

 Told(0,j) = -1.0

END DO

IF(indexsnowinit.EQ.0) THEN

indice_R = (INT(z_R/dz)) + 1
indexribbony1z = indice_R - INT(ribbonz*0.5/dz)
indexribbony2z = indice_R + INT(ribbonz*0.5/dz)

ELSE

indice_R = (INT(z_R/dz)) + 1 + indexsnowinit + 1
indexribbony1z = indice_R - INT(ribbonz*0.5/dz)
indexribbony2z = indice_R + INT(ribbonz*0.5/dz)

END IF

IF(indice_R.EQ.1) THEN
indice_R = 2
END IF


IF(MOD(nribbon,2).NE.0) THEN
   DO j = -FLOOR(nribbon/2.0),FLOOR(nribbon/2.0)
      indexribbony1 = INT((j*Delta_R + REAL(Ly)/2.0 - ribbony/2.0)/dy)
      indexribbony2 = INT((j*Delta_R + REAL(Ly)/2.0 + ribbony/2.0)/dy)
      DO i = indexribbony1+2,indexribbony2+1
         DO m = indexribbony1z,indexribbony2z
         w_R(m,i) = 1.0
         !w_R(indice_R,i) = 1.0
         END DO
      END DO
   END DO
ELSE
   DO j = -INT(nribbon/2.0)+1,INT(nribbon/2.0)
      indexribbony1 = INT((Delta_R*(j-0.5) + Ly/2.0 - ribbony/2.0)/dy)
      indexribbony2 = INT((Delta_R*(j-0.5) + Ly/2.0 + ribbony/2.0)/dy)
      IF (j.EQ.0) THEN
         DO i = indexribbony1+1,indexribbony2+1
            DO m = indexribbony1z,indexribbony2z
            w_R(m,i) = 1.0
            END DO
         END DO
      ELSE
         DO i = indexribbony1+1,indexribbony2+1
            DO m = indexribbony1z,indexribbony2z
            w_R(m,i) = 1.0
            END DO
         END DO
      END IF
   END DO
END IF

IF(indexsnowinit.EQ.0) THEN

indexl1 = INT(L1/dz) + 1
indexl2 = INT((L1+L2)/dz) + 1

ELSE

indexl1 = INT(L1/dz) + 1 + indexsnowinit + 1
indexl2 = INT((L1+L2)/dz) + 1 + indexsnowinit + 1

END IF

DO i = 1,allocationdepth
DO j = 1,Ny
   numberindexribbon = numberindexribbon + w_R(i,j)
END DO
END DO

indexNyleft = INT((Ly*0.5 - Delta_R*0.5)/dy) + 1
indexNyright = INT((Ly*0.5)/dy) + 1

!indexNyleft = INT((Ly*0.5 - Delta_R)/dy) + 0
!indexNyright = INT((Ly*0.5 + Delta_R)/dy) + 1

IF (1.EQ.0) THEN

OPEN(78, IOSTAT=stat, FILE='Temperature2D.dat', STATUS='old')
IF (stat.EQ.0) THEN
CLOSE(78, STATUS='delete')
END IF

OPEN(26,FILE='Temperature2D.dat',POSITION='append')

DO j = 1,INT(allocationdepth/3.0)
DO m = indexNyleft,indexNyright
WRITE(26,*) (j-1)*dz,(m-1)*dy,w_R(j,m)
END DO
WRITE(26,'(A)')
END DO

PRINT*,numberindexribbon,(indexsnowinit-1)*dz

!DO i=1,Ny
!WRITE(26,*) Tnew(i,:)
!WRITE (26, '("X")' )
!END DO

CLOSE(26)

STOP

END IF

END SUBROUTINE init
