PROGRAM temp
USE m_global
USE oalloc
IMPLICIT NONE

! MEASURING CPU TIME
CALL CPU_TIME(t1)

! READING THE SIMULATION PARAMETERS
CALL readpara()

CALL lire()

! CALLING SUBROUTINES TO GENERATE SOLAR RADIATION + AIR TEMPERATURE


! DEFINING THE SPATIAL INCREMENT
dz = L/(Nz-1.0)

dy = Ly/(Ny-1.0)

! DEFINING INVERSE OF COMMONLY USED VARIABLES TO OPTIMIZE SERIAL PERFORMANCE
rdz2 = 1.0/(dz**2.0)

rdy2 = 1.0/(dy**2.0)

rkdz = 1.0/(k*dz)

rkdrdz = 1.0/(k*Delta_R*dz)

rk1 = 1.0/(k1)

rk2 = 1.0/(k2)

avgrho0 = 2.0*rho*rho_snow/(rho_snow+rho)
avgrho1 = 2.0*rho*rho1/(rho+rho1)
avgrho2 = 2.0*rho2*rho1/(rho2+rho1)

avgc1 = 2.0*c*c1/(c+c1)
avgc2 = 2.0*c2*c1/(c2+c1)

avgk1 = 2.0*k*k1/(k+k1)
avgk2 = 2.0*k2*k1/(k2+k1)

! END OF DEFINING INVERSE

! SUBROUTINE INIT INITIALIZE THE MATRIX Tnew and Told + BOUNDARY CONDITIONS

nlinescase = 11*10000

CALL init()

avgalpha1 = 2.0*alpha*alpha1/(alpha+alpha1)
avgalpha2 = 2.0*alpha2*alpha1/(alpha2+alpha1)

rribbonyribbonz = REAL(nribbon)/REAL(dy*dz*numberindexribbon)

! ALLOCATE MEMORY FOR Maxsurf AND Tdepth
! Maxsurf is used to track the surface temperature and plot the min and max values at the end
! Tdepth is used to track the temperature at one time step along the depth
ALLOCATE(Maxsurf(nstep))
ALLOCATE(Tdepth(Nz))

listalpha(1) = alpha
listalpha(2) = alpha1
listalpha(3) = alpha2
minalpha = MAXVAL(listalpha)

countdays = 0

!nlinescase = 1000000

!nlinescase = INT(1000000.0/dt)

PRINT*,'dy = ',dy,'dz = ',dz, nlinescase

IF(0.EQ.0) THEN
DO i = 1,1500
   !IF((MOD(i,1)).EQ.0)THEN
      IF(countdays.LT.10) THEN
         WRITE(FileName,'(A,I1.1,A)') 'Tdepth.',countdays,'.dat'
         OPEN(26, IOSTAT=stat, FILE=FileName, STATUS='old')
         IF (stat.EQ.0) THEN
            CLOSE(26, STATUS='delete')
         END IF
      ELSE IF(countdays.LT.100) THEN
         WRITE(FileName,'(A,I2.2,A)') 'Tdepth.',countdays,'.dat'
         OPEN(26, IOSTAT=stat, FILE=FileName, STATUS='old')
         IF (stat.EQ.0) THEN
            CLOSE(26, STATUS='delete')
         END IF
      ELSE IF(countdays.LT.1000) THEN
         WRITE(FileName,'(A,I3.3,A)') 'Tdepth.',countdays,'.dat'
         OPEN(26, IOSTAT=stat, FILE=FileName, STATUS='old')
         IF (stat.EQ.0) THEN
            CLOSE(26, STATUS='delete')
         END IF
      ELSE
         WRITE(FileName,'(A,I4.4,A)') 'Tdepth.',countdays,'.dat'
         OPEN(26, IOSTAT=stat, FILE=FileName, STATUS='old')
         IF (stat.EQ.0) THEN
            CLOSE(26, STATUS='delete')
         END IF
      END IF
      countdays = countdays + 1
   !END IF
END DO
END IF

countdays = 0

!STOP

!OPEN(77, IOSTAT=stat, FILE='Temperature.dat', STATUS='old')
!IF (stat.EQ.0) THEN
!CLOSE(77, STATUS='delete')
!END IF

OPEN(78, IOSTAT=stat, FILE='Temperature2D.dat', STATUS='old')
IF (stat.EQ.0) THEN
CLOSE(78, STATUS='delete')
END IF

parite = 2

IF(MOD(parite+1,2).EQ.0) THEN
OPEN(51,FILE='Tmap.dat',STATUS='old')
DO i=1,Nz
READ(51,*) Tnew(i,:)
END DO
CLOSE(51)
Told = Tnew

OPEN(52,FILE='Taux_time.dat',STATUS='old')
DO i=1,Nz
READ(52,*) Taux_time(i,:)
END DO
CLOSE(52)

OPEN(53,FILE='Taux_value.dat',STATUS='old')
DO i=1,Nz
READ(53,*) Taux_value(i,:)
END DO
CLOSE(53)

OPEN(54,FILE='cumul_melted.dat',STATUS='old')
DO i=1,Ny
READ(54,*) cumul_melted(i)
END DO
CLOSE(54)

END IF



OPEN(16,FILE='Temperature.dat',POSITION='append')


! OPEN PARALLELIZATION REGION
!!$OMP Parallel DEFAULT(__auto)
!!$OMP DO

!power = 10.0

PRINT*,'-----------------------------------'
PRINT*,'------START OF CALCULATIONS--------'
PRINT*,'-----------------------------------'

! LOOP OVER THE NUMBER OF STEPS
DO i = 1,nlinescase

V_air(i) = 2.0
T_air(i) = -5.0
f_s(i) = 0.0
!fD(i) = 0.0

   Tplot = Tnew
   DO j = 1,indexsnowinit + indexl1
      DO m = indexNyleft-1,indexNyright+1
         IF (Taux_value(j,m).EQ.1.0) THEN
            Tplot(j,m) = 1000.0
         END IF
      END DO
   END DO

   hc = 1.78*V_air(i) + 2.21

!!$OMP DO

   cumul_lateral = 0

      DO WHILE(Tnew(indexsnowinit,cumul_lateral+indexNyleft+1).GE.&
      temp_melted.AND.cumul_lateral.LE.indexNyright-indexNyleft) 

         cumul_lateral = cumul_lateral + 1

      END DO

   IF (cumul_lateral.EQ.indexNyright-indexNyleft) THEN

      shift_right = 1

      !PRINT*,'HEY'
      !READ(*,*)

   ELSE

      shift_right = 0

   END IF


DO m = indexNyleft+1,indexNyright-1

! DO LOOP FOR FIRST LAYER

   indexsnowtot(m) = 0
   indexmeltup(m) = 0
   indexmeltdown(m) = 0
   indexsnowdead(m) = 0
   indexmeltedold(m) = 0
   cumul_vertical(m) = 0

IF(shift_right.GT.0) THEN

   DO j = 1,indexsnowinit

      IF(Tnew(j,m).GE.temp_melted) THEN
         indexsnowtot(m) = indexsnowtot(m) + 1
         IF((Taux_value(j,m).EQ.1.0).AND.(Taux_time(j,m).GT.0.0)&
         .AND.(Taux_time(j,m).LT.i + nlinescase*0)) THEN 
            indexmeltedold(m) = indexmeltedold(m) + 1
         ELSE IF(Taux_value(j,m).EQ.0.0) THEN
            Taux_value(j,m) = 1.0
            Taux_time(j,m) = i + nlinescase*0
         END IF
      END IF

   END DO

   !DO j = 1,indexsnowinit

   !   IF(Tnew(j,m).GE.temp_melted) THEN
   !      Taux_new(j,m) = 1.0
   !   END IF

   !END DO

   !cumul_vertical(m) = SUM(Taux_new(1:indexsnowinit,m))&
   !- SUM(Taux_old(1:indexsnowinit,m))

   !IF(cumul_vertical(m).LT.0) THEN

   !   cumul_vertical(m) = 0

   !END IF

   indexsnowtot(indexNyleft) = indexsnowtot(indexNyleft+1)
   indexsnowtot(indexNyright) = indexsnowtot(indexNyright-1)

   !cumul_melted(m) = cumul_melted(m) + indexsnowtot(m)

   !cumul_melted(indexNyleft) = cumul_melted(indexNyleft+1)
   !cumul_melted(indexNyright) = cumul_melted(indexNyright-1)

END IF

   DO j = 1,indexsnowinit

      IF(Tnew(j,m).GE.temp_melted) THEN
         Taux_new(j,m) = 1.0
      END IF

   END DO

   cumul_vertical(m) = SUM(Taux_new(1:indexsnowinit,m))&
   - SUM(Taux_old(1:indexsnowinit,m))

   IF(cumul_vertical(m).LT.0) THEN

      cumul_vertical(m) = 0

   END IF

   cumul_melted(m) = cumul_melted(m) + cumul_vertical(m)

   cumul_melted(indexNyleft) = cumul_melted(indexNyleft+1)
   cumul_melted(indexNyright) = cumul_melted(indexNyright-1)

END DO ! Shortening of END DO for "m" to fix the values of indexsnowtot(m)

DO m = indexNyleft+1,indexNyright-1 ! RESTART THINGS FOR "m"

   DO j = 1,indexsnowinit ! COUNT HOW MANY NODES ARE MELTED WITHIN SNOW PACK

      !PRINT*,indexsnowtot(indexNyleft),indexsnowtot(indexNyleft+1)

      indexsnowtot(indexNyleft) = indexsnowtot(indexNyleft+1)
      indexsnowtot(indexNyright) = indexsnowtot(indexNyright-1)

      Taux_value(j,indexNyleft) = Taux_value(j,indexNyleft+1)
      Taux_value(j,indexNyright) = Taux_value(j,indexNyright-1)

      Taux_time(j,indexNyleft) = Taux_time(j,indexNyleft+1)
      Taux_time(j,indexNyright) = Taux_time(j,indexNyright-1)

      !PRINT*,indexsnowtot(indexNyleft),indexsnowtot(indexNyleft+1)
      !READ(*,*)

   END DO ! END COUNT HOW MANY NODES ARE MELTED WITHIN SNOW PACK

   !IF(indexsnowtot(m).GE.2) THEN
   !   PRINT*,indexsnowtot(m),Taux_value(1:indexsnowinit,m)
   !   READ(*,*)
   !END IF

   IF(indexsnowtot(m).EQ.indexsnowinit) THEN ! CHECK HOW MANY NODES ARE MELTED
      
      DO j = 1, indexsnowinit + 1 ! EVERYTHING MELTED INCLUDING INTERFACE
         Tnew(j,m) = T_air(i)
         Told(j,m) = T_air(i)
         Tnew(j,indexNyleft) = T_air(i)
         Tnew(j,indexNyright) = T_air(i)
      END DO
      GO TO 87 ! GO TO ASPHALT LAYER

   ELSE ! NOT EVERYTHING IS MELTED

      DO j = 1,indexsnowinit

         IF((Tnew(j,m).GE.temp_melted).AND.(Tnew(j+1,m).LT.temp_melted))&
         THEN ! HOW MANY NODES MELTED UP?

            indexmeltup(m) = j
            indexmeltuptot(m) = j

         ELSE IF((Tnew(j+1,m).GE.temp_melted).AND.(Tnew(j,m).LT.temp_melted))&
         THEN ! HOW MANY NODES MELTED DOWN?

            indexmeltdown(m) = j + 1
            indexmeltdowntot(m) = j + 1

         END IF

      END DO

      indexmeltup(indexNyleft) = indexmeltup(indexNyleft+1)
      indexmeltup(indexNyright) = indexmeltup(indexNyright-1)

      indexmeltuptot(indexNyleft) = indexmeltuptot(indexNyleft+1)
      indexmeltuptot(indexNyright) = indexmeltuptot(indexNyright-1)

      indexmeltdown(indexNyleft) = indexmeltdown(indexNyleft+1)
      indexmeltdown(indexNyright) = indexmeltdown(indexNyright-1)

      indexmeltdowntot(indexNyleft) = indexmeltdown(indexNyleft+1)
      indexmeltdowntot(indexNyright) = indexmeltdown(indexNyright-1)

      DO j = 1,indexsnowinit ! COUNT HOW MANY NODES ARE INACTIVE
         
         IF(Tnew(j,m).EQ.T_air(i)) THEN
            indexsnowdead(m) = indexsnowdead(m) + 1
         END IF

      END DO ! END COUNT HOW MANY NODES ARE INACTIVE

      indexsnowdead(indexNyleft) = indexsnowdead(indexNyleft+1)
      indexsnowdead(indexNyright) = indexsnowdead(indexNyright-1)

      indexmeltup(m) = indexmeltup(m) - indexsnowdead(m)
      indexmeltdown(m) = indexsnowinit - indexmeltdown(m) + 1
      !indexmeltdown(m) = indexsnowtot(m) - indexsnowdead(m) - indexmeltup(m)

      indexmeltup(indexNyleft) = indexmeltup(indexNyleft+1)
      indexmeltup(indexNyright) = indexmeltup(indexNyright-1)

      indexmeltdown(indexNyleft) = indexmeltdown(indexNyleft+1)
      indexmeltdown(indexNyright) = indexmeltdown(indexNyright-1)

      indexmeltdowntot(indexNyleft) = indexmeltdown(indexNyleft+1)
      indexmeltdowntot(indexNyright) = indexmeltdown(indexNyright-1)

   !indexsnowtot(m) = 0
   !indexmeltup(m) = 0
   !indexmeltdown(m) = 0
   !indexsnowdead(m) = 0

      !IF(indexsnowtot(m).GT.1) THEN
      !PRINT*,'HERE2',Tnew(1:indexsnowinit,m)
      !PRINT*,'HERE2','indexsnowtot: ',indexsnowtot(m), 'indexmeltdown:',indexmeltdown(m)
      !READ(*,*)
      !END IF

      !DO j = indexsnowtot(m) + 2, indexsnowinit ! REASSIGN NODES

      !   Tnew(j,indexNyleft) = Tnew(j,indexNyleft+1)
      !   Told(j,indexNyleft) = Told(j,indexNyleft+1)

      !  Tnew(j,m) = Told(j-indexmeltdown(m),m)
         !Told(j,m) = Told(j-indexmeltdown(m),m)

      !   Tnew(j,indexNyright) = Tnew(j,indexNyright-1)
      !   Told(j,indexNyright) = Told(j,indexNyright-1)

      !END DO ! END REASSIGN NODES

      !Tnew(1 + indexsnowtot(m),m) = Told(1 + indexsnowdead(m),m) ! REASSIGN TOP NODE FOR HEAT FLUX
      !Told(1 + indexsnowtot(m),m) = Told(1 + indexmeltuptot(m),m) ! REASSIGN TOP NODE FOR HEAT FLUX

      !Tnew(1 + indexsnowtot(m),indexNyleft) =&
      !Tnew(1 + indexmeltuptot(m),indexNyleft+1) ! REASSIGN TOP NODE FOR HEAT FLUX

      !Tnew(1 + indexsnowtot(m),indexNyright) =&
      !Tnew(1 + indexmeltuptot(m),indexNyright-1) ! REASSIGN TOP NODE FOR HEAT FLUX

      IF(shift_right.GT.0) THEN
      DO j = 1, indexsnowtot(m) ! DEACTIVATE DEAD NODES

         Tnew(j,m) = T_air(i)
         Told(j,m) = T_air(i)
         Tnew(j,indexNyleft) = T_air(i)
         Tnew(j,indexNyright) = T_air(i)

      END DO ! END DEACTIVATE DEAD NODES
      END IF

      !IF(indexsnowtot_right.GE.1) THEN
      !PRINT*,'HERE3',Tnew(1:indexsnowinit,indexNyright-1)
      !PRINT*,'HERE3','indexsnowtot: ',indexsnowtot(m), 'indexmeltdown:',indexmeltdown(m)
      !READ(*,*)
      !END IF

      IF(shift_right.EQ.0) THEN
         indexsnowtot(m) = 0
      !ELSE
      !   PRINT*,indexsnowtot(indexNyleft:indexNyright),m,shift_right,i
      !   READ(*,*)
      END IF

      IF((Told(1,m).GE.T_l).AND.(Told(1,m).LE.temp_ref)) THEN
         theta(1,m) = -theta_snow*Told(1,m)/(T_l) + theta_snow
         d_T_theta(1,m) = -theta_snow/(T_l)
      ELSE IF (Told(1,m).LT.T_l) THEN
         theta(1,m) = 0.0
         d_T_theta(1,m) = 0.0
      ELSE
         theta(1,m) = theta_snow
         d_T_theta(1,m) = 0.0
      END IF

      exponentlatent = theta(1,m)/theta_snow

      c_eff(1,m) = c_snow_f*(1.0-exponentlatent) + &
      c_snow_u*exponentlatent
      
      k_eff(1,m) = k_snow_f**(1.0-exponentlatent) * &
      k_snow_u**exponentlatent

      Tnew(1,m) = Told(1,m) + ( 2.0*&
      k_eff(1,m)*dt*rdz2 * (Told(2,m)&
      - Told(1,m) + dz/k_eff(1,m) * &
      (a*f_s(i) + hc*(T_air(i) - Told(1,m)) + 0.0 * &
      (fD(i) - sigma*(Told(1,m) + 273.15)**4.0)  )  ) &
      + k_eff(1,m)*dt*rdy2 * (Told(1,m+1)&
      - 2.0*Told(1,m) + Told(1,m-1) ) ) /&
      ( c_eff(1,m) * rho_snow + latentheat * &
      d_T_theta(1,m))

      Tnew(1,indexNyright) = &
      Tnew(1,indexNyright-1) 

      Tnew(1,indexNyleft) = & ! END START HEAT FLUX SNOW LAYER
      Tnew(1,indexNyleft+1)

      !indexsnowtot(m) = 0

      DO j = 2, indexsnowinit ! SNOW LAYER

         IF(j-indexsnowtot(m).LE.0) THEN
         Tnew(j-indexsnowtot(m),m) = T_air(i)
         Told(j-indexsnowtot(m),m) = T_air(i)
         Tnew(j-indexsnowtot(m),indexNyleft) = T_air(i)
         Tnew(j-indexsnowtot(m),indexNyright) = T_air(i)
         END IF

         Tnew(j,indexNyleft) = Tnew(j,indexNyleft+1)

         IF((Told(j,m).GE.T_l).AND.(Told(j,m).LE.temp_ref)) THEN
            theta(j,m) = -theta_snow*Told(j,m)/(T_l) + theta_snow
            d_T_theta(j,m) = -theta_snow/(T_l)
         ELSE IF (Told(j,m).LT.T_l) THEN
            theta(j,m) = 0.0
            d_T_theta(j,m) = 0.0
         ELSE
            theta(j,m) = theta_snow
            d_T_theta(j,m) = 0.0
         END IF

         exponentlatent = theta(j,m)/theta_snow

         c_eff(j,m) = c_snow_f*(1.0-exponentlatent) + &
         c_snow_u*exponentlatent
      
         k_eff(j,m) = k_snow_f**(1.0-exponentlatent) * &
         k_snow_u**exponentlatent

         Tnew(j,m) = Told(j-indexsnowtot(m),m) + ( k_eff(j,m)*dt*rdz2*( &
         Told(j-1-indexsnowtot(m),m) + Told(j+1-indexsnowtot(m),m) - 2.0*&
         Told(j-indexsnowtot(m),m) ) + k_eff(j,m)*dt*rdy2*( &
         Told(j-indexsnowtot(m-1),m-1) + Told(j-indexsnowtot(m+1),m+1) - 2.0*&
         Told(j-indexsnowtot(m),m) ) )/( c_eff(j,m) * rho_snow + latentheat*d_T_theta(j,m) )

         Tnew(j,indexNyright) = Tnew(j,indexNyright-1)

         !PRINT*,Tnew(j,m),(j-1)*dz,(m-1)*dy

      END DO ! END SNOW LAYER

      Tnew(indexsnowinit + 1,indexNyleft) = & ! INTERFACE SNOW / ASPHALT
      Tnew(indexsnowinit + 1,indexNyleft+1)

      avgc1 = 2.0*(c_eff(indexsnowinit + 2,m)*&
      c_eff(indexsnowinit - 0,m))/(c_eff(indexsnowinit &
      + 2,m)+c_eff(indexsnowinit - 0,m))

      avgk1 = 2.0*(k_eff(indexsnowinit + 2,m)*&
      k_eff(indexsnowinit - 0,m))/(k_eff(indexsnowinit &
      + 2,m)+k_eff(indexsnowinit - 0,m))

      IF((d_T_theta(indexsnowinit + 2,m).EQ.0.0).AND.&
      (d_T_theta(indexsnowinit - 0,m).EQ.0.0)) THEN
            avgd_T_theta_snow = 0.0
      ELSE
         avgd_T_theta_snow = 2.0*(d_T_theta(indexsnowinit + 2,m)*&
         d_T_theta(indexsnowinit - 0,m))/(d_T_theta(indexsnowinit&
         + 2,m)+d_T_theta(indexsnowinit - 0,m))
      END IF

      Tnew(indexsnowinit + 1,m) = Told(indexsnowinit + 1,m) &
      + ( dt*rdz2 * ( k_eff(indexsnowinit + 2,m)*Told(indexsnowinit + 2 &
      ,m)-(k_eff(indexsnowinit + 2,m) + k_eff(indexsnowinit + 0 &
      ,m))*Told(indexsnowinit + 1,m) + k_eff(indexsnowinit + 0  &
      ,m) *Told(indexsnowinit + 0,m) ) + avgk1 * dt * rdy2 * &
      ( Told(indexsnowinit + 1,m+1) + Told(indexsnowinit + 1,m-1) &
      -2.0*Told(indexsnowinit + 1,m) ) ) / (avgc1*avgrho0 + latentheat*avgd_T_theta_snow)

      Tnew(indexsnowinit + 1,indexNyright) = &
      Tnew(indexsnowinit + 1,indexNyright-1) ! END INTERFACE SNOW / ASPHALT

      DO j = indexsnowinit + 2, indexsnowinit + 2 + indexl1 - 1 ! ASPHALT LAYER

         Tnew(j,indexNyleft) = Tnew(j,indexNyleft+1)

         k_eff(j,m) = k

         c_eff(j,m) = c

         Tnew(j,m) = Told(j,m) + ( k_eff(j,m)*dt*rdz2*( Told(j-1,m) + Told(j+1,m) -&
         2.0*Told(j,m) ) + k_eff(j,m)*dt*rdy2*( Told(j,m-1) + Told(j,m+1) -&
         2.0*Told(j,m) ) + dt*power*w_R(j,m)*rribbonyribbonz)/( c_eff(j,m) * rho )

         Tnew(j,indexNyright) = Tnew(j,indexNyright-1)

      END DO ! END ASPHALT LAYER

      GO TO 89

   END IF

   87 k_eff(indexsnowinit + 1 + 1,m) = k
   c_eff(indexsnowinit + 1 + 1,m) = c
   d_T_theta(indexsnowinit + 1 + 1,m) = 0.0

   Tnew(indexsnowinit + 1 + 1,m) = Told(indexsnowinit + 1 + 1,m) + ( 2.0*&
   k_eff(indexsnowinit + 1 + 1,m)*dt*rdz2 * (Told(indexsnowinit + 1 + 2,m)&
   - Told(indexsnowinit + 1 + 1,m) +  dz/k_eff(indexsnowinit + 1 + 1,m) * &
   (a*f_s(i) + hc*(T_air(i) - Told(indexsnowinit + 1 + 1,m)) + eps_s * &
   (fD(i) - sigma*(Told(indexsnowinit + 1 + 1,m) + 273.15)**4.0)  )  ) &
   + k_eff(indexsnowinit + 1 + 1,m)*dt*rdy2 * (Told(indexsnowinit + 1 + 1,m+1)&
   - 2.0*Told(indexsnowinit + 1 + 1,m) + Told(indexsnowinit + 1 + 1,m-1) ) ) /&
   ( c_eff(indexsnowinit + 1 + 1,m) * rho + latentheat * &
   d_T_theta(indexsnowinit + 1 + 1,m))

   !PRINT*,'HERE3',indexsnowtot(m)

   Tnew(indexsnowinit + 1 + 1,indexNyleft) =&
   Tnew(indexsnowinit + 1 + 1,indexNyleft+1)

   Tnew(indexsnowinit + 1 + 1,indexNyright) =&
   Tnew(indexsnowinit + 1 + 1,indexNyright-1)

   DO j = indexsnowinit + 3, indexsnowinit + 2 + indexl1 - 1

      Tnew(j,indexNyleft) = Tnew(j,indexNyleft+1)

      c_eff(j,m) = c
      
      k_eff(j,m) = k

      Tnew(j,m) = Told(j,m) + ( k_eff(j,m)*dt*rdz2*( Told(j-1,m) + Told(j+1,m) -&
      2.0*Told(j,m) ) + k_eff(j,m)*dt*rdy2*( Told(j,m-1) + Told(j,m+1) -&
      2.0*Told(j,m) ) + dt*power*w_R(j,m)*rribbonyribbonz)/( c_eff(j,m) * rho )

      Tnew(j,indexNyright) = Tnew(j,indexNyright-1)

   END DO
!   !$OMP END DO

!   !$OMP END DO
! CALCULATIONS FOR FIRST INTERFACE

   Tnew(indexsnowinit + 2 + indexl1,indexNyleft) = &
   Tnew(indexsnowinit + 2 + indexl1,indexNyleft+1)

   89 avgc1 = 2.0*(c_eff(indexsnowinit + 2 + indexl1-1,m)*&
   c_eff(indexsnowinit + 2 + indexl1+1,m))/(c_eff(indexsnowinit &
   + 2 + indexl1 - 1,m)+c_eff(indexsnowinit + 2 + indexl1 + 1,m))

   avgk1 = 2.0*(k_eff(indexsnowinit + 2 + indexl1-1,m)*&
   k_eff(indexsnowinit + 2 + indexl1+1,m))/(k_eff(indexsnowinit &
   + 2 + indexl1 - 1,m)+k_eff(indexsnowinit + 2 + indexl1 + 1,m))

   IF((d_T_theta(indexsnowinit + 2 + indexl1-1,m).EQ.0.0).AND.&
   (d_T_theta(indexsnowinit + 2 + indexl1+1,m).EQ.0.0)) THEN
      avgd_T_theta1 = 0.0
   ELSE
      avgd_T_theta1 = 2.0*(d_T_theta(indexsnowinit + 2 + indexl1-1,m)*&
      d_T_theta(indexsnowinit + 2 + indexl1+1,m))/(d_T_theta(indexsnowinit&
      + 2 + indexl1-1,m)+d_T_theta(indexsnowinit + 2 + indexl1+1,m))
   END IF

   Tnew(indexsnowinit + 2 + indexl1,m) = Told(indexsnowinit + 2 + indexl1,m) &
   + ( dt*rdz2 * ( k_eff(indexsnowinit + 2 + indexl1+1,m)*Told(indexsnowinit + 2 &
   + indexl1+1,m)-(k_eff(indexsnowinit + 2 + indexl1+1,m) + k_eff(indexsnowinit + 2 &
   + indexl1-1,m))*Told(indexsnowinit + 2 + indexl1,m) + k_eff(indexsnowinit + 2 + &
   indexl1-1,m) *Told(indexsnowinit + 2 + indexl1-1,m) ) + avgk1 * dt * rdy2 * &
   ( Told(indexsnowinit + 2 + indexl1,m+1) + Told(indexsnowinit + 2 + indexl1,m-1) &
   -2.0*Told(indexsnowinit + 2 + indexl1,m) ) ) / (avgc1*avgrho1 + latentheat*avgd_T_theta1)

   Tnew(indexsnowinit + 2 + indexl1,indexNyright) = &
   Tnew(indexsnowinit + 2 + indexl1,indexNyright-1)

! DO LOOP FOR SECOND LAYER

   DO j = indexsnowinit + 2 + indexl1 + 1,indexsnowinit + 2 + indexl2 - 1

      Tnew(j,indexNyleft) = Tnew(j,indexNyleft+1)

      IF(Told(j,m).GE.temp_ref) THEN
         theta(j,m) = theta_0_1
         d_T_theta(j,m) = 0.0
      ELSE
         theta(j,m) = ((temp_ref-Told(j,m))/(temp_ref+273.15))**beta
         theta(j,m) = theta_0_1*(1.0 - theta(j,m))
         d_T_theta(j,m) = ((temp_ref-Told(j,m))/(temp_ref+273.15))**beta
         d_T_theta(j,m) = theta_0_1*beta*d_T_theta(j,m)/&
                        (temp_ref - Told(j,m))
      END IF

      exponentlatent = theta(j,m)/theta_0_1

      c_eff(j,m) = c1frozen*(1.0-exponentlatent) + &
      c1unfrozen*exponentlatent
      
      k_eff(j,m) = k1frozen**(1.0-exponentlatent) * &
      k1unfrozen**exponentlatent

      Tnew(j,m) = Told(j,m) + ( k_eff(j,m)*dt*rdz2*( Told(j-1,m) + Told(j+1,m) -&
      2.0*Told(j,m) ) + k_eff(j,m)*dt*rdy2*( Told(j,m-1) + Told(j,m+1) -&
      2.0*Told(j,m) ) )/( c_eff(j,m) * rho1 + latentheat*d_T_theta(j,m) )

      Tnew(j,indexNyright) = Tnew(j,indexNyright-1)

   END DO

! CALCULATIONS FOR SECOND INTERFACE

   Tnew(indexsnowinit + 2 + indexl2,indexNyleft) = &
   Tnew(indexsnowinit + 2 + indexl2,indexNyleft+1)

   avgc2 = 2.0*(c_eff(indexsnowinit + 2 + indexl2-1,m)*&
   c_eff(indexsnowinit + 2 + indexl2+1,m))/(c_eff(indexsnowinit + 2 &
   + indexl2-1,m)+c_eff(indexsnowinit + 2 + indexl2+1,m))

   avgk2 = 2.0*(k_eff(indexsnowinit + 2 + indexl2-1,m)*&
   k_eff(indexsnowinit + 2 + indexl2+1,m))/(k_eff(indexsnowinit + 2 &
   + indexl2-1,m)+k_eff(indexsnowinit + 2 + indexl2+1,m))

   IF((d_T_theta(indexsnowinit + 2 + indexl2-1,m).EQ.0.0).AND.&
   (d_T_theta(indexsnowinit + 2 + indexl2+1,m).EQ.0.0)) THEN
      avgd_T_theta2 = 0.0
   ELSE
      avgd_T_theta2 = 2.0*(d_T_theta(indexsnowinit + 2 + indexl2-1,m)*&
      d_T_theta(indexsnowinit + 2 + indexl2+1,m))/(d_T_theta(indexsnowinit + 2 &
      + indexl2-1,m)+d_T_theta(indexsnowinit + 2 + indexl2+1,m))
   END IF

   Tnew(indexsnowinit + 2 + indexl2,m) = Told(indexsnowinit + 2 + indexl2,m) &
   + ( dt*rdz2 * ( k_eff(indexsnowinit + 2 + indexl2+1,m)*Told(indexsnowinit + 2 &
   + indexl2+1,m)-  (k_eff(indexsnowinit + 2 + indexl2+1,m) + k_eff(indexsnowinit + 2 &
   + indexl2-1,m))*Told(indexsnowinit + 2 + indexl2,m) + k_eff(indexsnowinit + 2 + &
   indexl2-1,m) * Told(indexsnowinit + 2 + indexl2-1,m) ) + avgk2 * dt * rdy2 * &
   ( Told(indexsnowinit + 2 + indexl2,m+1) + Told(indexsnowinit + 2 + indexl2,m-1) &
   -2.0*Told(indexsnowinit + 2 + indexl2,m) ) ) / (avgc2*avgrho2 + latentheat*avgd_T_theta2)

   Tnew(indexsnowinit + 2 + indexl2,indexNyright) = &
   Tnew(indexsnowinit + 2 + indexl2,indexNyright-1)


! DO LOOP FOR THIRD LAYER

   DO j = indexsnowinit + 2 + indexl2+1, allocationdepth-1

      Tnew(j,indexNyleft) = Tnew(j,indexNyleft+1)

      IF(Told(j,m).GE.temp_ref) THEN
         theta(j,m) = theta_0_2
         d_T_theta(j,m) = 0.0
      ELSE
         theta(j,m) = ((temp_ref-Told(j,m))/(temp_ref+273.15))**beta
         theta(j,m) = theta_0_2*(1.0 - theta(j,m))
         d_T_theta(j,m) = ((temp_ref-Told(j,m))/(temp_ref+273.15))**beta
         d_T_theta(j,m) = theta_0_2*beta*d_T_theta(j,m)/&
                        (temp_ref - Told(j,m))
      END IF

      exponentlatent = theta(j,m)/theta_0_2

      c_eff(j,m) = c2frozen*(1.0-exponentlatent) + &
      c2unfrozen*exponentlatent
      
      k_eff(j,m) = k2frozen**(1.0-exponentlatent) * &
      k2unfrozen**exponentlatent

      Tnew(j,m) = Told(j,m) + ( k_eff(j,m)*dt*rdz2*( Told(j-1,m) + Told(j+1,m) -&
      2.0*Told(j,m) ) + k_eff(j,m)*dt*rdy2*( Told(j,m-1) + Told(j,m+1) -&
      2.0*Told(j,m) ) )/( c_eff(j,m) * rho2 + latentheat*d_T_theta(j,m) )

      Tnew(j,indexNyright) = Tnew(j,indexNyright-1)

   END DO

! UPDATE ARRAY 
!   !$OMP WORKSHARE
   
!   !$OMP END WORKSHARE

   Tnew(allocationdepth,m) =&
   Tnew(allocationdepth-1,m)

   Tnew(allocationdepth,indexNyleft) =&
   Tnew(allocationdepth-1,indexNyleft)

   Tnew(allocationdepth,indexNyright) =&
   Tnew(allocationdepth-1,indexNyright)

END DO ! END FOR m (Ny)

DO j = 1,allocationdepth
   Tnew(j,indexNyleft) = Tnew(j,indexNyleft+1)
   Told(j,indexNyleft) = Told(j,indexNyleft+1)
END DO

DO j = 1,allocationdepth
   Tnew(j,indexNyright) = Tnew(j,indexNyright-1)
   Told(j,indexNyright) = Told(j,indexNyright-1)
END DO

!PRINT*,REAL(i)/REAL(nstep)*100, '%',Tnew(1,INT(Ny/2.0))

!DO j = 1,Nz
!   Tnew(j,Ny) = Tnew(j,Ny-1)
!   Tnew(j,1) = Tnew(j,2)
!END DO

Told = Tnew
Taux_old = Taux_new
Taux_new = 0

   IF((MOD(i-1,1000).EQ.0)) THEN !.AND.(i.GT.1150*60).AND.(i.LT.1160*60))THEN
      IF(countdays.LT.10) THEN
        WRITE(FileName,'(A,I1.1,A)') 'Tdepth.',countdays,'.dat'
      ELSE IF(countdays.LT.100) THEN
         WRITE(FileName,'(A,I2.2,A)') 'Tdepth.',countdays,'.dat'
      ELSE IF(countdays.LT.1000) THEN
         WRITE(FileName,'(A,I3.3,A)') 'Tdepth.',countdays,'.dat'
      ELSE 
         WRITE(FileName,'(A,I4.4,A)') 'Tdepth.',countdays,'.dat'
      END IF
      OPEN(36,FILE=Filename,POSITION='append')
      !DO j = 1,indexsnowinit + indice_R + 10
         DO m = indexNyleft,indexNyright
            WRITE(36,*) cumul_melted(m)*dz
         END DO
         !WRITE(36,'(A)')
      !END DO
      countdays = countdays + 1
   END IF

!(m-1)*dy, (indexsnowinit - cumul_melted(m))*dz 
!(j-1)*dz,(m-1)*dy,indexsnowtot(m),Tnew(j,m)

END DO

PRINT*,'-----------------------------------'
PRINT*,'--------END OF CALCULATIONS--------'
PRINT*,'-----------------------------------'

!!$OMP END DO
!!$OMP END Parallel

!WRITE(16,*) Tnew(INT(Nz/2.0)+1+10,INT(Ny/2.0)+1+10),Tnew(INT(Nz/2.0)+1-10,INT(Ny/2.0)+1-10)

CLOSE(16)

OPEN(28,FILE='Temperature2D.dat',POSITION='append')

DO j = 1,indexsnowinit + indice_R + 10
   DO m = indexNyleft,indexNyright
      WRITE(28,*) (j-1)*dz,(m-1)*dy,Tnew(j,m)
   END DO
   WRITE(28,'(A)')
END DO

CLOSE(28)

IF(MOD(parite,2).EQ.0) THEN

OPEN(79, IOSTAT=stat, FILE='Tmap.dat', STATUS='old')
IF (stat.EQ.0) THEN
CLOSE(79, STATUS='delete')
END IF

OPEN(82, IOSTAT=stat, FILE='Taux_value_right.dat', STATUS='old')
IF (stat.EQ.0) THEN
CLOSE(82, STATUS='delete')
END IF

OPEN(84, IOSTAT=stat, FILE='Taux_time_right.dat', STATUS='old')
IF (stat.EQ.0) THEN
CLOSE(84, STATUS='delete')
END IF

OPEN(86, IOSTAT=stat, FILE='cumul_melted.dat', STATUS='old')
IF (stat.EQ.0) THEN
CLOSE(86, STATUS='delete')
END IF

OPEN(26,FILE='Tmap.dat',POSITION='append')

DO i=1,Nz
WRITE(26,*) Tnew(i,:)
!WRITE (26, '("X")' )
END DO

CLOSE(26)

OPEN(28,FILE='Taux_value.dat',POSITION='append')

DO i=1,Nz
WRITE(28,*) Taux_value(i,:)
!WRITE (28, '("X")' )
END DO

CLOSE(28)

OPEN(30,FILE='Taux_time.dat',POSITION='append')

DO i=1,Nz
WRITE(30,*) Taux_time(i,:)
!WRITE (30, '("X")' )
END DO

CLOSE(30)

OPEN(32,FILE='cumul_melted.dat',POSITION='append')

DO i=1,Ny
WRITE(32,*) cumul_melted(i)
!WRITE (32, '("X")' )
END DO

CLOSE(32)

END IF


! END MAIN PROGRAM
END PROGRAM temp
