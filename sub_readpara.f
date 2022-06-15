SUBROUTINE readpara()
USE m_global
IMPLICIT NONE

NAMELIST /para/dt,Nz,nstep,rho,c,k,Deltat_SS,Pswr_n,Tbase_air,Tamp_air, &
               L,z_R,Delta_R,Delta_z,N_winter,sigma,eps_s,a,T_e,Tol,&
               power,c1,rho1,k1,c2,rho2,k2,z_R,L1,L2,z_e,Ly,Ny,nribbon,ribbony,&
               ribbonz,latentheat,c2unfrozen,c2frozen,k2unfrozen,k2frozen,&
               c1unfrozen,c1frozen,k1unfrozen,k1frozen,theta_0_1,theta_0_2,&
               beta,temp_ref,depthsnowinit,absorb_snow,c_snow,k_snow,rho_snow,&
               emit_snow,theta_snow,c_snow_u,c_snow_f,k_snow_u,k_snow_f,&
               c_air,k_air,rho_air,temp_melted,T_l
               

OPEN(20,FILE='para.dat')

READ(20,NML=para)
!WRITE(*,NML=para)

CLOSE(20)


END SUBROUTINE readpara
