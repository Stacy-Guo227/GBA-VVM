MODULE advec_3d_module

USE kinds
USE parmsld
USE constld
USE const3d
USE profoutld
USE bound
USE domain_decomposition
!USE nvtx


IMPLICIT NONE
PRIVATE

PUBLIC :: advec_3d

CONTAINS

      SUBROUTINE ADVEC_3D (Q,terma, termf)
!     Advection for thermodynamic variables 
!     ALADV: alpha in advection scheme

! argument list declarations
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,NK3), INTENT(in)  :: &
         Q     !  the quantity advected

      REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,NK2), INTENT(out) ::  &  
         terma      ! advective tendency
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk3), INTENT(out) ::  & 
         termf      ! vertical flux convergence

! local variables
      REAL (KIND=dbl_kind) :: &
         TEMPIZ(mim:mip,mjm:mjp,NK2),TEMPIM(mim:mip,mjm:mjp,NK2),TEMPIV(mim:mip,mjm:mjp,NK2),UPIZ(mim:mip,mjm:mjp,NK2),UPIM(mim:mip,mjm:mjp,NK2),UPIV(mim:mip,mjm:mjp,NK2),UMIZ(mim:mip,mjm:mjp,NK2),UMIM(mim:mip,mjm:mjp,NK2),UMIV(mim:mip,mjm:mjp,NK2)
      REAL (KIND=dbl_kind) :: &
         UPSRIZ(mim:mip,mjm:mjp,NK2),UPSRIM(mim:mip,mjm:mjp,NK2),UPSRIV(mim:mip,mjm:mjp,NK2),UMSRIZ(mim:mip,mjm:mjp,NK2),UMSRIM(mim:mip,mjm:mjp,NK2),UMSRIV(mim:mip,mjm:mjp,NK2)
      REAL (KIND=dbl_kind) :: &
         FLXIZ(0:MI1,0:mj1,NK2),FLXIM(0:MI1,0:mj1,NK2),FLXIV(0:MI1,0:mj1,NK2)
      REAL (KIND=dbl_kind), DIMENSION(2,mjm:mjp,nk3) :: qew,uew  ! extended copy of q  ew
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,2,nk3) :: qns,vns  ! extended copy of q  ns

      INTEGER (KIND=int_kind) :: &
         i, j, k    ! do loop indices for zonal, meridional and vertical dimensions

      call extend_ne(nk3,q,qew,qns,u3dx,uew,u3dy,vns)

!     Zonal advection
      !$acc data copyin(UEW,QEW,VNS,QNS,Q) copyout(terma) create(TEMPIZ,TEMPIM,TEMPIV,UPIZ,UPIM,UPIV,UMIZ,UMIM,UMIV,UPSRIZ,UPSRIM,UPSRIV,UMSRIZ,UMSRIM,UMSRIV,FLXIZ,FLXIM,FLXIV)

!     loop(1) 320
      !$acc parallel async(0)
      !$acc loop collapse (2)
      DO  K = 2, NK2
      DO  J = 1, MJ1
      UEW(1,J,K) = UEW(1,J,K)*RHOU(K) ! use for u(-1,J,K) =uew(1,J,K)
      UEW(1,J,K) = 0.5*(UEW(1,J,K)+ABS(UEW(1,J,K)))
      UEW(1,J,K) = SQRT(UEW(1,J,K))
      !$acc loop
      DO  I=mim,mip
      TEMPIZ(I,J,K)=U3DX(I,J,K)*RHOU(K)
      UPIZ(I,J,K)=0.5*(TEMPIZ(I,J,K)+ABS(TEMPIZ(I,J,K)))
      UMIZ(I,J,K)=0.5*(TEMPIZ(I,J,K)-ABS(TEMPIZ(I,J,K)))
      UPSRIZ(I,J,K)=SQRT(UPIZ(I,J,K))
      UMSRIZ(I,J,K)=SQRT(ABS(UMIZ(I,J,K)))
      ENDDO
      ENDDO
      ENDDO
      !$acc end parallel

!     loop(2) 330
!     loop(3) 337

!     loop(4) 340
      !$acc parallel async(0)
      !$acc loop collapse (2)
      DO  K=2,NK2
      DO  J=1,MJ1
      I = 0
      FLXIZ(I,J,K)=TEMPIZ(I,J,K)*(Q(I+1,J,K)+Q(I,J,K))     &
       -ALADV*(UPIZ(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))         &
       -UPSRIZ(I,J,K)*UEW(1,J,K)*(Q(I,J,K)-QEW(1,J,K)) &
       +UMIZ(I,J,K)*(Q(I,J,K)-Q(I+1,J,K))                &
       +UMSRIZ(I,J,K)*UMSRIZ(I+1,J,K)*(Q(I+1,J,K)-Q(I+2,J,K)))/3.
!     inner loop(1): calFLXI 341
      !$acc loop
      DO  I=1,MI1-1
      FLXIZ(I,J,K)=TEMPIZ(I,J,K)*(Q(I+1,J,K)+Q(I,J,K))      &
       -ALADV*(UPIZ(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))         &
       -UPSRIZ(I,J,K)*UPSRIZ(I-1,J,K)*(Q(I,J,K)-Q(I-1,J,K)) &
       +UMIZ(I,J,K)*(Q(I,J,K)-Q(I+1,J,K))                &
       +UMSRIZ(I,J,K)*UMSRIZ(I+1,J,K)*(Q(I+1,J,K)-Q(I+2,J,K)))/3.
      ENDDO
      I=MI1
      FLXIZ(I,J,K)=TEMPIZ(I,J,K)*(Q(I+1,J,K)+Q(I,J,K))      &
      -ALADV*(UPIZ(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))         &
      -UPSRIZ(I,J,K)*UPSRIZ(I-1,J,K)*(Q(I,J,K)-Q(I-1,J,K)) &
       +UMIZ(I,J,K)*(Q(I,J,K)-Q(I+1,J,K))                &
       +UMSRIZ(I,J,K)*UMSRIZ(I+1,J,K)*(Q(I+1,J,K)-QEW(2,J,K)))/3.
      ENDDO
      ENDDO
      !$acc end parallel

!     loop(5) 350
!      !$acc parallel
!      !$acc loop collapse (3)
!      DO  K=2,NK2
!      DO  J=1,MJ1
!      DO  I=1,MI1
!      terma(I,J,K)=-(FLXIZ(I,J,K)-FLXIZ(I-1,J,K))/(2.*DX)
!      ENDDO
!      ENDDO
!      ENDDO
!      !$acc end parallel

!      IF (.FALSE.) THEN

!     Meridional advection
!     loop(6) 420
      !$acc parallel async(1)
      !$acc loop collapse (2)
      DO  K = 2, NK2
      DO  I = 1,MI1
      VNS(I,1,K) = RHOU(K)*VNS(I,1,K)
      VNS(I,1,K) = 0.5*(VNS(I,1,K)+ABS(VNS(I,1,K)))
      VNS(I,1,K) = SQRT(VNS(I,1,K))
      !$acc loop
      DO  J = mjm,mjp
      TEMPIM(I,J,K)=U3DY(I,J,K)*RHOU(K)
      UPIM(I,J,K)=0.5*(TEMPIM(I,J,K)+ABS(TEMPIM(I,J,K)))
      UMIM(I,J,K)=0.5*(TEMPIM(I,J,K)-ABS(TEMPIM(I,J,K)))
      UPSRIM(i,J,K)=SQRT(UPIM(i,J,K))
      UMSRIM(i,J,K)=SQRT(ABS(UMIM(i,J,K)))
      ENDDO
      ENDDO
      ENDDO
      !$acc end parallel

!     loop(7) 430
!      CALL bound_ns(nk1, upi(mim,mjm,2))
!      CALL bound_ns(nk1, umi(mim,mjm,2))
!     loop(8) 437

!     loop(9) 440
      !$acc parallel async(1)
      !$acc loop collapse (2)
      DO  K=2,NK2
      DO  I=1,MI1
      J = 0 
      FLXIM(I,J,K)=TEMPIM(I,J,K)*(Q(I,J+1,K)+Q(I,J,K))       &
       -ALADV*(UPIM(I,J,K)*(Q(I,J+1,K)-Q(I,J,K))          &
       -UPSRIM(I,J,K)*VNS(I,1,K)*(Q(I,J,K)-QNS(I,1,K))  &
       +UMIM(I,J,K)*(Q(I,J,K)-Q(I,J+1,K))                 &
       +UMSRIM(I,J,K)*UMSRIM(I,J+1,K)*(Q(I,J+1,K)-Q(I,J+2,K)))/3.
!     inner loop(2): calFLXI 441
      !$acc loop
      DO  J=1,MJ1-1
      FLXIM(I,J,K)=TEMPIM(I,J,K)*(Q(I,J+1,K)+Q(I,J,K))       &
       -ALADV*(UPIM(I,J,K)*(Q(I,J+1,K)-Q(I,J,K))          &
       -UPSRIM(I,J,K)*UPSRIM(I,J-1,K)*(Q(I,J,K)-Q(I,J-1,K))  &
       +UMIM(I,J,K)*(Q(I,J,K)-Q(I,J+1,K))                 &
       +UMSRIM(I,J,K)*UMSRIM(I,J+1,K)*(Q(I,J+1,K)-Q(I,J+2,K)))/3.
      ENDDO
      J = MJ1
      FLXIM(I,J,K)=TEMPIM(I,J,K)*(Q(I,J+1,K)+Q(I,J,K))       &
       -ALADV*(UPIM(I,J,K)*(Q(I,J+1,K)-Q(I,J,K))          &
       -UPSRIM(I,J,K)*UPSRIM(I,J-1,K)*(Q(I,J,K)-Q(I,J-1,K))  &
       +UMIM(I,J,K)*(Q(I,J,K)-Q(I,J+1,K))                 &
       +UMSRIM(I,J,K)*UMSRIM(I,J+1,K)*(Q(I,J+1,K)-QNS(I,2,K)))/3.
      ENDDO
      ENDDO
      !$acc end parallel
!      CALL send_north(nk1,flxi(0,0,2))
!     loop(10) 450
!      !$acc parallel
!      !$acc loop collapse (3)
!      DO  K=2,NK2
!      DO  J=1,MJ1
!      DO  I=1,MI1
!      terma(I,J,K)=terma(I,J,K)-(FLXIM(I,J,K)-FLXIM(I,J-1,K))/(2.*DYNEW)
!      ENDDO
!      ENDDO
!      ENDDO
!      !$acc end parallel
!     Vertical advection

!     loop(11) 520
      !$acc parallel async(2)
      !$acc loop collapse (3)
      DO  K=2,NK1
      DO  J=1,MJ1
      DO  I=1,MI1
      TEMPIV(I,J,K)=W3D(I,J,K)*RHOZ(K)
      UPIV(I,J,K)=0.5*(TEMPIV(I,J,K)+ABS(TEMPIV(I,J,K)))
      UMIV(I,J,K)=0.5*(TEMPIV(I,J,K)-ABS(TEMPIV(I,J,K)))
      UPSRIV(I,J,K)=SQRT(UPIV(I,J,K))
      UMSRIV(I,J,K)=SQRT(ABS(UMIV(I,J,K)))
      ENDDO
      ENDDO
      ENDDO
      !$acc end parallel

!     loop(12) 530
!     loop(13) 537

!     loop(14) 540
      !$acc parallel async(2)
      !$acc loop collapse (3)
      DO  K=3,NK1-1
      DO  J=1,MJ1
      DO  I=1,MI1
      FLXIV(I,J,K)=TEMPIV(I,J,K)*(Q(I,J,K+1)+Q(I,J,K))       &
       -ALADV*(UPIV(I,J,K)*(Q(I,J,K+1)-Q(I,J,K))          &
       -UPSRIV(I,J,K)*UPSRIV(I,J,K-1)*(Q(I,J,K)-Q(I,J,K-1))  &
       +UMIV(I,J,K)*(Q(I,J,K)-Q(I,J,K+1))                 &
       +UMSRIV(I,J,K)*UMSRIV(I,J,K+1)*(Q(I,J,K+1)-Q(I,J,K+2)))/3.
      ENDDO
      ENDDO
      ENDDO
      !$acc end parallel

!     loop(15) 545
      !$acc parallel async(2)
      !$acc loop collapse (2)
      DO  J=1,MJ1
      DO  I=1,MI1
      IF(TEMPIV(I,J,NK1).GE.0.) THEN
         FLXIV(I,J,NK1)=TEMPIV(I,J,NK1)*(Q(I,J,NK2)+Q(I,J,NK1))   &
         -ALADV*(UPIV(I,J,NK1)*(Q(I,J,NK2)-Q(I,J,NK1))         &
         -UPSRIV(I,J,NK1)*UPSRIV(I,J,NK1-1)*(Q(I,J,NK1)-Q(I,J,NK1-1)))/3.
      ELSE
         FLXIV(I,J,NK1)=TEMPIV(I,J,NK1)*(Q(I,J,NK2)+Q(I,J,NK1))
      ENDIF

      IF(TEMPIV(I,J,2).GE.0.) THEN
         FLXIV(I,J,2)=TEMPIV(I,J,2)*(Q(I,J,3)+Q(I,J,2))
      ELSE
         FLXIV(I,J,2)=TEMPIV(I,J,2)*(Q(I,J,3)+Q(I,J,2))  &
           -ALADV*(UMIV(I,J,2)*(Q(I,J,2)-Q(I,J,3))    &
           +UMSRIV(I,J,2)*UMSRIV(I,J,3)*(Q(I,J,3)-Q(I,J,4)))/3.
      ENDIF
      FLXIV(I,J, 1)= 0.
      FLXIV(I,J,NK2)= 0.
      ENDDO
      ENDDO
      !$acc end parallel
      !$acc wait

!     loop(16) 547
!      !$acc parallel async(2)
!      !$acc loop collapse (2)
!      DO  J=1,MJ1
!      DO  I=1,MI1
!      FLXIV(I,J, 1)= 0.
!      FLXIV(I,J,NK2)= 0.
!      ENDDO
!      ENDDO
!      !$acc end parallel

!     loop(17) 550
      !$acc parallel
      !$acc loop collapse (3)
      DO  K=2,NK2
      DO  J=1,MJ1
      DO  I=1,MI1
!      zonal
      terma(I,J,K)=-(FLXIZ(I,J,K)-FLXIZ(I-1,J,K))/(2.*DX)
!      meri
      terma(I,J,K)=terma(I,J,K)-(FLXIM(I,J,K)-FLXIM(I,J-1,K))/(2.*DYNEW)
!      vertical
      terma(I,J,K)=terma(I,J,K)-(FLXIV(I,J,K)-FLXIV(I,J,K-1))*FNU(K)/(2.*DZ)
      terma(I,J,K)=terma(I,J,K)/RHOU(K)
      !termf(I,J,K)=-(FLXI(I,J,K)-FLXI(I,J,K-1))*FNU(K)/(2.*DZ)/RHOU(K)
      IF(ITYPEW(I,J,K) .NE. 1) THEN
      terma(I,J,K)=0.
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      !$acc end parallel

!     loop(18) 570

!ccwut set flux convergence on topo to zero
!     loop(19)
!      !$acc parallel
!      !$acc loop collapse (3)
!      DO K=2,maxtopo
!      DO J=1,MJ1
!      DO I=1,MI1
!      IF(ITYPEW(I,J,K) .NE. 1) THEN
!      terma(I,J,K)=0.
!      ENDIF
!      ENDDO
!      ENDDO
!      ENDDO
!      !$acc end parallel
      !$acc end data
!ccwut 
   END SUBROUTINE advec_3d
END MODULE advec_3d_module
