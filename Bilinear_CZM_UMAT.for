!     Author: @Harshdeep_Sharma 
!     Affiliation: Indian Institute of Technology, Patna 
!     Find me @: Material Testing Lab, Block-3
!     Lab Incharge: Dr. Akhilendra Singh 
!     Contact me: harshsharma52@gmail.com 
!     Desc: Bilinear UMAT (2D/3D Interfacial problems) 
!     Designed for standard solver and suitable for quasi-static and static problems 
!
       SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
!
      INCLUDE 'ABA_PARAM.INC'
!
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)

!     ADDITIONAL VECTORS AND TENSORS 
      DIMENSION DSEC(NTENS,NTENS), DELTAT(NTENS), K_MAT(NTENS,NTENS) 

!      PARAMETER(MULT_FACTOR=2)

!     MATERIAL PROPERTIES FROM INPUT FILE : 
      A_KN = PROPS(1)       !NORMAL PENALTY STIFFNESS
      TAU_N = PROPS(2)      !NORMAL COHESIVE STRENGTH 
      TAU_S = PROPS(3)      !SHEAR COHESIVE STRENGTH 
      TAU_T = PROPS(4)      !TEAR COHESIVE STRENGTH 
      G_NC = PROPS(5)       !NORMAL FRACTURE TOUGHNESS
      G_SC = PROPS(6)       !SHEAR FRACTURE TOUGHNESS
      G_TC = PROPS(7)       !TEAR FRACTURE TOUGHNESS 
      ETA = PROPS(8)        !BK PARAMETER (ENERGY RELEASE RATE)

!     SHEAR PENALTY STIFFNESS 
      A_KS = (G_NC/G_SC)*((TAU_S/TAU_N)**2)*A_KN    
       
      IF(NTENS .EQ. 3) THEN 
         A_KT = (G_NC/G_TC)*((TAU_T/TAU_N)**2)*A_KN 
      ENDIF 

!     CONSTRUCT K_MAT 
      K_MAT(:,:) = 0. 
      K_MAT(1,1) = A_KN
      K_MAT(2,2) = A_KS 
      IF(NTENS .EQ. 3) THEN 
         K_MAT(3,3) = A_KS
      ENDIF  

!     RETRIEVING STATE VARIABLES 
      DMG_OLD = STATEV(1) 

!     SEPARATION @ {N+1}
      Do I = 1, NTENS  
         DELTAT(I) = STRAN(I) + DSTRAN(I)
      Enddo 

!     NORMAL AND SHEAR SEPARATION @ {N+1}
      DELTA_N = DELTAT(1)   
      DELTA_S = DELTAT(2) 
      IF(NTENS .EQ. 3) THEN !TEAR SEPARATION @ {N+1}
         DELTA_T = DELTAT(3)
      ENDIF

!     DISCRIMINATE 
!     CASE I: PURE COMPRESSIVE STATE 
      IF((NTENS .LT. 3 .AND. DELTA_N .LT. 0. .AND. ABS(DELTA_S) .LE.
     &  1.E-8 ) .OR. (NTENS .EQ. 3 .AND. DELTA_N .LT. 0. .AND.  
     &  ABS(DELTA_S) .LE. 1.E-8 .AND. ABS(DELTA_T) .LE. 1.E-8)) THEN 
        DO I=1,NTENS 
            STRESS(I) = K_MAT(I,I)*DELTAT(I)
        ENDDO 
        DDSDDE(:,:) = K_MAT(:,:)
      ELSE  
!     CASE II: DELTA_EQ WILL BE GREATER THAN 0 
!     MIXED-MODE RATIO (ENERGY)
      G_SHEAR = A_KS*DELTA_S**2.
      G_N = A_KN*MAX(DELTA_N,0.)**2
      IF(NTENS .EQ. 3) THEN 
         G_SHEAR = G_SHEAR + A_KT*DELTA_T**2 
      ENDIF 
      G_TOTAL = G_N + G_SHEAR  
      B = G_S/G_TOTAL 

!     MIXED-MODE ONSET SEPARATION VALUE 
      DELTA_NC = TAU_N/A_KN
      DELTA_SC = TAU_S/A_KS 
      A_KM = A_KN*(1.-B) + (A_KS)*B
      DELTA_MC = SQRT((A_KN*DELTA_NC**2 + (A_KS*DELTA_SC**2-
     &       A_KN*DELTA_NC**2)*(B**ETA))/A_KM) 
      IF(NTENS .EQ. 3) THEN 
         DELTA_TC = TAU_T/A_KT 
         A_KM = A_KN*(1.-B) + (A_KS+A_KT)*B
      DELTA_MC = SQRT((A_KN*DELTA_NC**2 + (A_KS*DELTA_SC**2+A_KT*
     &       DELTA_TC**2- A_KN*DELTA_NC**2)*(B**ETA))/A_KM) 
      ENDIF 
     
!     MIXED-MODE FINAL SEPARATION VALUE 
      DELTA_NF = 2.*G_NC/(A_KN*DELTA_NC)
      DELTA_SF = 2.*G_SC/(A_KS*DELTA_SC)
      DELTA_MF = ((A_KN*DELTA_NC*DELTA_NF)+((A_KS*DELTA_SC*
     & DELTA_SF)-(A_KN*DELTA_NC*DELTA_NF))*(B**ETA))/(A_KM*DELTA_MC)
      IF(NTENS .EQ. 3) THEN 
        DELTA_TF = 2.*G_TC/(A_KT*DELTA_TC)
        DELTA_MF = ((A_KN*DELTA_NC*DELTA_NF)+((A_KS*DELTA_SC*DELTA_SF+
     &       A_KT*DELTA_TC*DELTA_TF)-(A_KN*DELTA_NC*DELTA_NF))*
     &      (B**ETA))/(A_KM*DELTA_MC)
      ENDIF 

!     MIXED-MODE SEPARATION VALUE AT THE CURRENT INCREMENT 
!      DELTA_M = SQRT(DELTA_N**2. + DELTA_S**2.)
      A_NUM = A_KN*MAX(DELTA_N,0.)**2 + A_KS*DELTA_S**2 
      A_DENOM = A_KN**2*MAX(DELTA_N,0.)**2 + A_KS**2*DELTA_S**2
      IF(NTENS .EQ. 3) THEN 
        A_NUM = A_KN*MAX(DELTA_N,0.)**2 + A_KS*DELTA_S**2 + A_KT*
     &  DELTA_T**2  
        A_DENOM = A_KN**2*MAX(DELTA_N,0.)**2 + A_KS**2*DELTA_S**2+
     &  A_KT**2*DELTA_T**2  
      ENDIF 
      DELTA_M = A_NUM/SQRT(A_DENOM)
!     DAMAGE CALCULATION AT THE CURRENT INCREMENT 
      RT_OLD =(DELTA_MC*DELTA_MF)/(DELTA_MF-DMG_OLD*(DELTA_MF-DELTA_MC))
!      RT = MAX(RT_OLD,DELTA_M)
      IF(DELTA_M .GT. RT_OLD) THEN
        DMG =(DELTA_MF*(DELTA_M-DELTA_MC))/(DELTA_M*(DELTA_MF-DELTA_MC))
        IF(DMG .GE. 1.) DMG = 1.  
      ELSE
        DMG = DMG_OLD
      ENDIF

!     STORING THE STATE VARIABLES FOR THE NEXT INCREMENT 
      STATEV(1) = DMG 

!     CALCULATION OF SECANT STIFFNESS MATRIX 
      DSEC(:,:) = 0. 
      DSEC(:,:) = (1.-DMG)*K_MAT

!     CALCULATION OF TRACTION FOR THE CURRENT INCREMENT (BASED ON SECANT STIFFNESS)
      Do I = 1, NTENS  
         STRESS(I) = 0. 
         Do J = 1, NTENS  
            STRESS(I) = STRESS(I) + DSEC(I,J) * DELTAT(J) 
         Enddo 
      Enddo 

!     CALCULATION OF TANGENT STIFFNESS MATRIX 
      DDSDDE = 0.
      IF(DELTA_M .GT. RT_OLD .AND. DELTA_M .LT. DELTA_MF ) THEN
         COEFF1 = ((DELTA_MF*DELTA_MC)/((DELTA_MF-DELTA_MC)*DELTA_M**2))
         FIRST_TERM = 2.*A_KS*DELTA_S/SQRT(A_DENOM)
         SECOND_TERM=-1.*A_NUM*A_KS**2*DELTA_S/(SQRT(A_DENOM)*A_DENOM)
         DELTAD_DELTAS = COEFF1*(FIRST_TERM+SECOND_TERM)
         IF(NTENS .EQ. 3) THEN 
            FIRST_TERM = 2.*A_KT*DELTA_T/SQRT(A_DENOM)
            SECOND_TERM=-1.*A_NUM*A_KT**2*DELTA_T/(SQRT(A_DENOM)*
     &             A_DENOM)
            DELTAD_DELTAT = COEFF1*(FIRST_TERM+SECOND_TERM)
         ENDIF 
         DDSDDE(2,2) = DSEC(2,2) - A_KS*DELTA_S*DELTAD_DELTAS
         IF(NTENS .EQ. 3) THEN 
            DDSDDE(3,2) = DSEC(3,2) - A_KT*DELTA_T*DELTAD_DELTAS
            DDSDDE(2,3) = DSEC(2,3) - A_KS*DELTA_S*DELTAD_DELTAT
            DDSDDE(3,3) = DSEC(3,3) - A_KT*DELTA_T*DELTAD_DELTAT
         ENDIF 
         IF(DELTA_N.LT.0.) THEN 
            DDSDDE(:,1) = DSEC(:,1)
            DDSDDE(1,:) = DSEC(1,:)
         ELSE
           FIRST_TERM = 2.*A_KN*DELTA_N/SQRT(A_DENOM)
           SECOND_TERM=-1.*A_NUM*A_KN**2*DELTA_N/(SQRT(A_DENOM)*A_DENOM)
           DELTAD_DELTAN = COEFF1*(FIRST_TERM+SECOND_TERM)
           DDSDDE(1,1) = DSEC(1,1) - A_KN*DELTA_N*DELTAD_DELTAN
           DDSDDE(2,1) = DSEC(2,1) - A_KS*DELTA_S*DELTAD_DELTAN
           DDSDDE(1,2) = DSEC(1,2) - A_KN*DELTA_N*DELTAD_DELTAS
           IF(NTENS .EQ. 3) THEN 
             DDSDDE(3,1) = DSEC(3,1) - A_KT*DELTA_T*DELTAD_DELTAN
             DDSDDE(1,3) = DSEC(1,3) - A_KN*DELTA_N*DELTAD_DELTAT
           ENDIF 
         ENDIF 
      ELSE
         DDSDDE(:,:) = DSEC(:,:)
      ENDIF

      ENDIF 

      RETURN
      END

    