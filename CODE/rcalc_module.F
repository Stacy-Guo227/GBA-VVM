#include "definesld.com"
MODULE rcalc_module

USE kinds
USE parmsld
USE profoutld
USE const3d
USE constld
USE workcom
USE nudge_module
USE update_thermo_module
USE force_3d_module
USE advec_3d_module
USE timer
USE domain_decomposition
USE nvtx
USE bound

IMPLICIT NONE
PRIVATE

PUBLIC :: rcalc_3d

CONTAINS

!=======================================================================
   SUBROUTINE RCALC_3D ( N1, N2, ITT )
!=======================================================================

!     ALL CALCULATIONS ASSOCIATED WITH THERMODYNAMIC VARIABLES
!     EXCEPT FOR THE DIFFUSION.
!     TEMPORARY USE OF TERM1 AS AN OUTPUT OF ADVECTION.

! Argument list variables
      INTEGER (KIND=int_kind), INTENT(IN) :: &
         itt,       & ! time step count
         n1,        & ! AB forcing time index for previous timestep
         n2           ! AB forcing time index for current timestep

!      REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk3) :: term4
      REAL (KIND=dbl_kind), DIMENSION(2,mjm:mjp,nk3) :: &
      THew, QVew, QRew, QCew, QIew, NCew, NRew, NIew, QRIMew, BRIMew, QSew, uew  ! extended copy of q  ew
      REAL (KIND=dbl_kind), DIMENSION(2,mjm:mjp,nk3, ntracer) :: &
      TCew
      REAL (KIND=dbl_kind), DIMENSION(2,mjm:mjp,nk3, nchem) :: &
      CHEMew
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,2,nk3) :: &
      THns, QVns, QRns, QCns, QIns, NCns, NRns, NIns, QRIMns, BRIMns, QSns, vns  ! extended copy of q  ns
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,2,nk3, ntracer) :: &
      TCns
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,2,nk3, nchem) :: &
      CHEMns
      INTEGER (KIND=int_kind) ::  &
         i, j, k, & ! do loop indices for zonal, meridional and vertical dimensions
         l, nt      ! time level index, tracer index
      REAL (KIND=dbl_kind) :: fntddz, pi

      INTEGER (KIND=int_kind) ::  &
         tempim  ! global index 

!-----------------------------------------------------------------------
      L = N2

!-----------------------------------------------------------------------

      call nvtxStartRange('advection')
      call timer_start('advection')
      !$acc data copyin(rhou, fnu, rhoz, itypew, u3dy, w3d, u3dx, TH3D, QV3D, QR3D, QC3D, QI3D, NC3D, NR3D, NI3D, QRIM3D, BRIM3D, QS3D, TC3D(:,:,:,1:ntracer), CHEM3D(:,:,:,1:nchem)) &
      !$acc create(FTH3D(:,:,:,L), FQV3D(:,:,:,L), FQR3D(:,:,:,L), FQC3D(:,:,:,L), FQI3D(:,:,:,L), FNC3D(:,:,:,L), FNR3D(:,:,:,L), FNI3D(:,:,:,L), FQRIM3D(:,:,:,L), FBRIM3D(:,:,:,L), FQS3D(:,:,:,L), FTC3D(:,:,:,L,1:ntracer), FCHEM3D(:,:,:,L,1:nchem)) &
      !$acc create(THew, THns, QVew, QVns, QRew, QRns, QCew, QCns, QIew, QIns, NCew, NCns, NRew, NRns, NIew, NIns, QRIMew, QRIMns, BRIMew, BRIMns, QSew, QSns, TCew, TCns, CHEMew, CHEMns, uew, vns)
      
      ! pass in 3D array and call extend_ne here
      call extend_ne(nk3, TH3D, THew, THns, U3DX, uew, U3DY, vns)
      !$acc update device(THew, THns, uew, vns)
#if defined (PHYSICS)
      call extend_ne(nk3, QV3D, QVew, QVns, U3DX, uew, U3DY, vns)
      !$acc update device(QVew, QVns)
      call extend_ne(nk3, QR3D, QRew, QRns, U3DX, uew, U3DY, vns)
      !$acc update device(QRew, QRns)
      call extend_ne(nk3, QC3D, QCew, QCns, U3DX, uew, U3DY, vns)
      !$acc update device(QCew, QCns)
      call extend_ne(nk3, QI3D, QIew, QIns, U3DX, uew, U3DY, vns)
      !$acc update device(QIew, QIns)

#if defined (MICROP3)
      call extend_ne(nk3, NC3D, NCew, NCns, U3DX, uew, U3DY, vns)
      !$acc update device(NCew, NCns)
      call extend_ne(nk3, NR3D, NRew, NRns, U3DX, uew, U3DY, vns)
      !$acc update device(NRew, NRns)
      call extend_ne(nk3, NI3D, NIew, NIns, U3DX, uew, U3DY, vns)
      !$acc update device(NIew, NIns)
      call extend_ne(nk3, QRIM3D, QRIMew, QRIMns, U3DX, uew, U3DY, vns)
      !$acc update device(QRIMew, QRIMns)
      call extend_ne(nk3, BRIM3D, BRIMew, BRIMns, U3DX, uew, U3DY, vns)
      !$acc update device(BRIMew, BRIMns)
#if defined (LIQFRACP3)
      call extend_ne(nk3, QILIQ3D, QILIQew, QILIQns, U3DX, uew, U3DY, vns)
      !$acc update device(QILIQew, QILIQns)
#endif
#else

      call extend_ne(nk3, QS3D, QSew, QSns, U3DX, uew, U3DY, vns)
      !$acc update device(QSew, QSns)
#endif
      if (ntracer > 0)then
         do nt = 1,ntracer
            call extend_ne(nk3, TC3D(mim,mjm,1,nt), TCew(:,:,:,nt), TCns(:,:,:,nt), U3DX, uew, U3DY, vns)
            !$acc update device(TCew(:,:,:,nt), TCns(:,:,:,nt))
         enddo
      endif

#if defined (CHEM)
         do nt = 1,nchem
            CALL extend_ne(nk3, CHEM3D(mim,mjm,1,nt), CHEMew(:,:,:,nt), CHEMns(:,:,:,nt), U3DX, uew, U3DY, vns)
            !$acc update device(CHEMew(:,:,:,nt), CHEMns(:,:,:,nt))
         enddo
#endif
#endif

      CALL ADVEC_3D ( TH3D, FTH3D(:,:,:,L), THew, THns, uew, vns ) !$acc async(1)
#if defined (PHYSICS)

      CALL ADVEC_3D ( QV3D, FQV3D(:,:,:,L), QVew, QVns, uew, vns ) !$acc async(1)
  
      CALL ADVEC_3D ( QR3D, FQR3D(:,:,:,L), QRew, QRns, uew, vns ) !$acc async(1)
      
      CALL ADVEC_3D ( QC3D, FQC3D(:,:,:,L), QCew, QCns, uew, vns ) !$acc async(2)

      CALL ADVEC_3D ( QI3D, FQI3D(:,:,:,L), QIew, QIns, uew, vns ) !$acc async(2)

#if defined (MICROP3)
      CALL ADVEC_3D ( NC3D, FNC3D(:,:,:,L), NCew, NCns, uew, vns ) !$acc async(3)

      CALL ADVEC_3D ( NR3D, FNR3D(:,:,:,L), NRew, NRns, uew, vns ) !$acc async(3)

      CALL ADVEC_3D ( NI3D, FNI3D(:,:,:,L), NIew, NIns, uew, vns ) !$acc async(3)

      CALL ADVEC_3D ( QRIM3D, FQRIM3D(:,:,:,L), QRIMew, QRIMns, uew, vns ) !$acc async(4)

      CALL ADVEC_3D ( BRIM3D, FBRIM3D(:,:,:,L), BRIMew, BRIMns, uew, vns ) !$acc async(4)
#if defined (LIQFRACP3)

      CALL ADVEC_3D ( QILIQ3D, FQILIQ3D(:,:,:,L), QILIQew, QILIQns, uew, vns ) !$acc async(4)
#endif
#else

      CALL ADVEC_3D ( QS3D, FQS3D(:,:,:,L), QSew, QSns, uew, vns ) !$acc async(4)
!      DO 160 K = 2, NK2
!      DO 160 J = 1, MJ1
!      DO 160 I = 1, MI1
!      HYSG3D(I,J,K)=HYSG3D(I,J,K) + TERM4(I,J,K)
!  160 CONTINUE

!      CALL ADVEC_3D ( QG3D, FQG3D(1,1,1,L), term4 )
!      DO 170 K = 2, NK2
!      DO 170 J = 1, MJ1
!      DO 170 I = 1, MI1
!      HYSG3D(I,J,K)=(HYSG3D(I,J,K) + TERM4(I,J,K)) &
!                               *(HLF/CP)*86400.
!  170 CONTINUE
#endif
      !$acc wait
      if (ntracer > 0)then
         do nt = 1,ntracer
            CALL ADVEC_3D ( TC3D(mim,mjm,1,nt), FTC3D(:,:,:,L,nt), TCew(:,:,:,nt), TCns(:,:,:,nt), uew, vns )
         enddo
      endif

#if defined (CHEM)
         do nt = 1,nchem
            CALL ADVEC_3D ( CHEM3D(mim,mjm,1,nt), FCHEM3D(:,:,:,L,nt), CHEMew(:,:,:,nchem), CHEMns(:,:,:,nchem), uew, vns )
         enddo
#endif


#endif
      !$acc update self(FTH3D(:,:,:,L), FQV3D(:,:,:,L), FQR3D(:,:,:,L), FQC3D(:,:,:,L), FQI3D(:,:,:,L), FNC3D(:,:,:,L), FNR3D(:,:,:,L), FNI3D(:,:,:,L), FQRIM3D(:,:,:,L), FBRIM3D(:,:,:,L), FQS3D(:,:,:,L), FTC3D(:,:,:,L,1:ntracer), FCHEM3D(:,:,:,L,1:nchem))
      !$acc end data
      call timer_stop('advection')
      call nvtxEndRange

      call nvtxStartRange('forcing')
      call timer_start('forcing')

     IF ( START .AND. ITT .LT. ITSTOP )   &
           CALL FORCE_3D (DTHMAX,DTPERT,Z1PERT,Z2PERT,l)

!   HERE Q1LS IS DEFINED AS 'ADVECTIVE HEATING' SUCH 
!   THAT NEGATIVE VALUES COOL THE UPPER LEVELS AND
!   POSITIVE VALUES WARM THE LOWER LEVELS THEREBY
!   MAKING THE ATMOSPHERE MORE UNSTABLE, 'FORCING'
!   IT TO CONVECT. 
!   SIMILARLY, Q2LS IS DEFINED AS 'ADVECTIVE MOISTENING'
!   WHEREBY NEGATIVE VALUES INDICATE DRYING.
!
!   ADDITIONALLY, CONDENSATE AND TRACERS ARE SUBJECT TO 
!   FORCING BY THE LARGE SCALE VERTICAL WIND.
!   (Modified June 2009, TRJ)
      IF (Q1Q2) THEN
      DO 300 K = 2, NK2
      DO 309 J = 1, MJ1
      DO 309 I = 1, MI1
      FTH3D(I,J,K,L)=FTH3D(I,J,K,L) + Q1LS (K)
      FQV3D(I,J,K,L)=FQV3D(I,J,K,L) + Q2LS (K)
  309 CONTINUE
  SELECT CASE (TRIM(casename))
  CASE ('TWP-ICE','BOMEX','GATE_PHASE_III')
      fntddz = FNT(K)/DZ
      IF (WLS(K) .LT. 0.) THEN
      
      DO 301 J = 1, MJ1
      DO 301 I = 1, MI1
      ! add WLS to TH & QV      
      FTH3D(I,J,K,L)=FTH3D(I,J,K,L) - WLS(K)*         &
                (TH3D(I,J,K+1)-TH3D(I,J,K))*fntddz
      FQV3D(I,J,K,L)=FQV3D(I,J,K,L) - WLS(K)*         &
                (QV3D(I,J,K+1)-QV3D(I,J,K))*fntddz

      FQC3D(I,J,K,L)=FQC3D(I,J,K,L) - WLS(K)*         &
                (QC3D(I,J,K+1)-QC3D(I,J,K))*fntddz
      FQI3D(I,J,K,L)=FQI3D(I,J,K,L) - WLS(K)*         &
                (QI3D(I,J,K+1)-QI3D(I,J,K))*fntddz
      FQR3D(I,J,K,L)=FQR3D(I,J,K,L) - WLS(K)*         &
                (QR3D(I,J,K+1)-QR3D(I,J,K))*fntddz
#if defined (MICROP3)
      FNC3D(I,J,K,L)=FNC3D(I,J,K,L) - WLS(K)*         &
                (NC3D(I,J,K+1)-NC3D(I,J,K))*fntddz
      FNR3D(I,J,K,L)=FNR3D(I,J,K,L) - WLS(K)*         &
                (NR3D(I,J,K+1)-NR3D(I,J,K))*fntddz
      FNI3D(I,J,K,L)=FNI3D(I,J,K,L) - WLS(K)*         &
                (NI3D(I,J,K+1)-NI3D(I,J,K))*fntddz
      FQRIM3D(I,J,K,L)=FQRIM3D(I,J,K,L) - WLS(K)*         &
                (QRIM3D(I,J,K+1)-QRIM3D(I,J,K))*fntddz
      FBRIM3D(I,J,K,L)=FBRIM3D(I,J,K,L) - WLS(K)*         &
                (BRIM3D(I,J,K+1)-BRIM3D(I,J,K))*fntddz
#if defined (LIQFRACP3)
      FQILIQ3D(I,J,K,L)=FQILIQ3D(I,J,K,L) - WLS(K)*         &
                (QILIQ3D(I,J,K+1)-QILIQ3D(I,J,K))*fntddz
#endif
#else
      FQG3D(I,J,K,L)=FQG3D(I,J,K,L) - WLS(K)*         &
                (QG3D(I,J,K+1)-QG3D(I,J,K))*fntddz
      FQS3D(I,J,K,L)=FQS3D(I,J,K,L) - WLS(K)*         &
                (QS3D(I,J,K+1)-QS3D(I,J,K))*fntddz
      HYLSTD(I,J,K)=(                                  &
       (- WLS(K)*(QR3D(I,J,K+1)-QR3D(I,J,K))  &
       - WLS(K)*(QS3D(I,J,K+1)-QS3D(I,J,K))  &
       - WLS(K)*(QG3D(I,J,K+1)-QG3D(I,J,K)))*fntddz) &
       *(HLF/CP)*86400.
#endif
  301 CONTINUE

      if (ntracer > 0)then

      DO 303 nt = 1,ntracer
      DO 303 J = 1, MJ1
      DO 303 I = 1, MI1
      FTC3D(I,J,K,L,nt)=FTC3D(I,J,K,L,nt) - WLS(K)*         &
                (TC3D(I,J,K+1,nt)-TC3D(I,J,K,nt))*fntddz
  303 CONTINUE

      endif

#if defined (CHEM)
      DO  nt = 1,nchem
      DO  J = 1, MJ1
      DO  I = 1, MI1
      FCHEM3D(I,J,K,L,nt)=FCHEM3D(I,J,K,L,nt) - WLS(K)*         &
                (CHEM3D(I,J,K+1,nt)-CHEM3D(I,J,K,nt))*fntddz
      ENDDO
      ENDDO
      ENDDO
#endif

      ELSE
      DO 302 J = 1, MJ1
      DO 302 I = 1, MI1
      FTH3D(I,J,K,L)=FTH3D(I,J,K,L) - WLS(K)*         &
                (TH3D(I,J,K)-TH3D(I,J,K-1))*fntddz
      FQV3D(I,J,K,L)=FQV3D(I,J,K,L) - WLS(K)*         &
                (QV3D(I,J,K)-QV3D(I,J,K-1))*fntddz


      FQC3D(I,J,K,L)=FQC3D(I,J,K,L) - WLS(K)*         &
                (QC3D(I,J,K)-QC3D(I,J,K-1))*fntddz
      FQI3D(I,J,K,L)=FQI3D(I,J,K,L) - WLS(K)*         &
                (QI3D(I,J,K)-QI3D(I,J,K-1))*fntddz
      FQR3D(I,J,K,L)=FQR3D(I,J,K,L) - WLS(K)*         &
                (QR3D(I,J,K)-QR3D(I,J,K-1))*fntddz
#if defined (MICROP3)
      FNC3D(I,J,K,L)=FNC3D(I,J,K,L) - WLS(K)*         &
                (NC3D(I,J,K)-NC3D(I,J,K-1))*fntddz
      FNR3D(I,J,K,L)=FNR3D(I,J,K,L) - WLS(K)*         &
                (NR3D(I,J,K)-NR3D(I,J,K-1))*fntddz
      FNI3D(I,J,K,L)=FNI3D(I,J,K,L) - WLS(K)*         &
                (NI3D(I,J,K)-NI3D(I,J,K-1))*fntddz
      FQRIM3D(I,J,K,L)=FQRIM3D(I,J,K,L) - WLS(K)*         &
                (QRIM3D(I,J,K)-QRIM3D(I,J,K-1))*fntddz
      FBRIM3D(I,J,K,L)=FBRIM3D(I,J,K,L) - WLS(K)*         &
                (BRIM3D(I,J,K)-BRIM3D(I,J,K-1))*fntddz
#if defined (LIQFRACP3)
      FQILIQ3D(I,J,K,L)=FQILIQ3D(I,J,K,L) - WLS(K)*         &
                (QILIQ3D(I,J,K)-QILIQ3D(I,J,K-1))*fntddz
#endif
#else
      FQG3D(I,J,K,L)=FQG3D(I,J,K,L) - WLS(K)*         &
                (QG3D(I,J,K)-QG3D(I,J,K-1))*fntddz
      FQS3D(I,J,K,L)=FQS3D(I,J,K,L) - WLS(K)*         &
                (QS3D(I,J,K)-QS3D(I,J,K-1))*fntddz
      HYLSTD(I,J,K)=(                                  &
      ( - WLS(K)*(QR3D(I,J,K)-QR3D(I,J,K-1))  &
        - WLS(K)*(QS3D(I,J,K)-QS3D(I,J,K-1))  &
        - WLS(K)*(QG3D(I,J,K)-QG3D(I,J,K-1)))*fntddz) &
        *(HLF/CP)*86400.
#endif
  302 CONTINUE

      if (ntracer > 0)then

      DO 304 nt = 1,ntracer
      DO 304 J = 1, MJ1
      DO 304 I = 1, MI1
      FTC3D(I,J,K,L,nt)=FTC3D(I,J,K,L,nt) - WLS(K)*         &
                (TC3D(I,J,K,nt)-TC3D(I,J,K-1,nt))*fntddz
  304 CONTINUE

      endif
#if defined (CHEM)
      DO nt = 1,nchem
      DO J = 1, MJ1
      DO I = 1, MI1
      FCHEM3D(I,J,K,L,nt)=FCHEM3D(I,J,K,L,nt) - WLS(K)*         &
                (CHEM3D(I,J,K,nt)-CHEM3D(I,J,K-1,nt))*fntddz
      ENDDO
      ENDDO
      ENDDO
#endif


      END IF
  END SELECT 
  300 CONTINUE
      END IF
      call timer_stop('forcing')
      call nvtxEndRange

      
!tac-------FOLLOWING MOVED TO physics_interface------------------------
!tac#if defined (PHYSICS)
!ctac -- Calculate cloud fraction (needed for CAM radiation)
!tac      CALL CLOUD_FRAC
!tac      CALL PHYSICS_3D ( L, ITT )
!tac#endif
!tac-------------------------------------------------------------------


      call nvtxStartRange('physics')
      call timer_start('physics')
! Initialize tendency terms to zero
      DO 510 K=2,NK2
      DO 510 J=1,MJ1
      DO 510 I=1,MI1
        THAD3(I,J,K) = 0.
        QVAD3(I,J,K) = 0.
        QCAD3(I,J,K) = 0.
        QIAD3(I,J,K) = 0.
        QRAD3(I,J,K) = 0.
#if defined (MICROP3)
        NCAD3(I,J,K)   = 0.
        NRAD3(I,J,K)   = 0.
        NIAD3(I,J,K)   = 0.
        QRIMAD3(I,J,K) = 0.
        BRIMAD3(I,J,K) = 0.
#if defined (LIQFRACP3)
        QILIQAD3(I,J,K) = 0.
#endif
#else
        QSAD3(I,J,K) = 0.
        QGAD3(I,J,K) = 0.
#endif
  510 CONTINUE

! Physics calculations
      CALL physics_interface(N1,N2,ITT,DT,ZZ,ZT,PBAR,PIBAR,  &
                             RHO,PSFC,NOTHERM,NRAD,tg(1:mi1,1:mj1),      &
                             DX,DYNEW,RLAT,RLON,THBAR,THBARZ)
      call timer_stop('physics')
      call nvtxEndRange



! Update thermodynamic variables -- Replaces call to CLOUD_3D
#if defined (MICROP3)
! Update thermodynamic variables -- call before P3 microphysics
#else
      CALL update_thermodynamics (N1,N2)
#endif

! Nudge mean water vapor and mean potential temperature to mean (if so defined in nudge)
      CALL NUDGE


   END SUBROUTINE rcalc_3d

END MODULE rcalc_module
