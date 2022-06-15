MODULE m_global
USE oalloc
IMPLICIT NONE
INTEGER :: Nz, counter, n,Tbase_air, Tamp_air,stat
INTEGER :: i,j, info, nstep, Pswr_n, N_winter, T_e
REAL :: dt, dz, t1, t2,alpha,k,z_R,Delta_R, Delta_z,z_e
REAL :: rho, c, k, Deltat_SS, sigma, eps_s, a,L,somme,hc,V_air
REAL:: Tol,power,alpha1,c1,rho1,k1,rk1,rk2
REAL::c2,rho2,k2,alpha2,L1,L2,rdz2,rkdz,rkdrdz,rk1drdz,rk2drdz
REAL:: Ly,dy,rdy2,minalpha,mindiscr,ribbony,ribbonz
REAL:: rribbonyribbonz
INTEGER:: indexl1,indexl2,Ny,m,nribbon, indexribbony1, indexribbony2
INTEGER:: indice_R,Middleribbon,nlinescase,indexribbony1z,indexribbony2z
REAL(KIND(1.0E0)),DIMENSION(:,:),ALLOCATABLE:: w_R,Told,Tnew,Tplot
REAL,DIMENSION(3):: listalpha
REAL,DIMENSION(2):: listdiscr
REAL,DIMENSION(:),ALLOCATABLE:: f_s,T_air,Maxsurf,Tdepth,V_air, fD
INTEGER :: countdays,nstepbegin,ndays,m,indexfak1,indexfak2,indexfak3
CHARACTER(LEN=20)::Filename
REAL:: avgrho1,avgrho2,avgc1,avgc2,avgk1,avgk2,avgalpha1,avgalpha2
REAL(KIND(1.0D0)),DIMENSION(:,:),ALLOCATABLE:: c_eff,k_eff,theta,d_T_theta
REAL::latentheat,theta_0_1,theta_0_2,avgd_T_theta1,avgd_T_theta2
REAL(KIND(1.0D0))::c1unfrozen, c1frozen, k1unfrozen, k1frozen,beta
REAL(KIND(1.0D0))::c2unfrozen, c2frozen, k2unfrozen, k2frozen
REAL::exponentlatent,temp_ref, numberindexribbon
INTEGER:: indexNyleft,indexNyright
REAL::depthsnowinit, absorb_snow, c_snow, k_snow, rho_snow
INTEGER::indexsnowinit, indexstartsnow, allocationdepth,parite
INTEGER,DIMENSION(:),ALLOCATABLE:: depthsnow, snowhelp
REAL::alpha_snow,emit_snow,theta_snow
INTEGER,DIMENSION(:),ALLOCATABLE:: indexsnowtot,indexmeltup,indexmeltdown
INTEGER,DIMENSION(:),ALLOCATABLE:: indexsnowdead,indexmeltuptot,indexmeltdowntot
REAL::c_snow_u,c_snow_f,k_snow_u,k_snow_f,c_avg_snow,k_avg_snow,avgd_T_theta_snow
REAL::c_air,k_air,rho_air,temp_melted,T_l
REAL(KIND(1.0D0)),DIMENSION(:,:),ALLOCATABLE:: rho_eff
REAL(KIND(1.0D0)),DIMENSION(:,:),ALLOCATABLE:: Taux_value, Taux_time
INTEGER,DIMENSION(:),ALLOCATABLE:: indexmeltedold
INTEGER::indexmeltedold_right, indexsnowtot_right, shift_right,shift_right_cumul
REAL(KIND(1.0D0)),DIMENSION(:,:),ALLOCATABLE:: Taux_value_right, Taux_time_right
INTEGER::cumul_lateral
INTEGER,DIMENSION(:),ALLOCATABLE::cumul_melted,cumul_vertical
REAL(KIND(1.0D0)),DIMENSION(:,:),ALLOCATABLE::Taux_new,Taux_old
REAL(KIND(1.0D0))::avgrho0
END MODULE m_global
