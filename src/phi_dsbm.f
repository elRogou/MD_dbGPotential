!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* dihedral computes the dihedral angles and the forces due to them   *
!**********************************************************************

        SUBROUTINE dbG_dihedral(dihedE_dbG)

        include 'MD.com'

      double precision lfac
      integer I3, J3, K3, L3
      real S,ftem
      DIMENSION XIJ(npa),YIJ(npa),ZIJ(npa),XKJ(npa),YKJ(npa),
     + ZKJ(npa),XKL(npa),YKL(npa),ZKL(npa),DX(npa),DY(npa),
     + DZ(npa), GX(npa),GY(npa),GZ(npa),CT(npa),CPHI(npa),
     + SPHI(npa),Z1(npa), Z2(npa),FXI(npa),FYI(npa),FZI(npa),
     + FXJ(npa),FYJ(npa),FZJ(npa), FXK(npa),FYK(npa),FZK(npa),
     + FXL(npa),FYL(npa),FZL(npa),DF(npa)
      real A1, mu1, sig1, A2, mu2, sig2, drinv, bydr12, rad,
     + elon1, elonsq1, sigsq1, twosigsq1,
     + elon2, elonsq2, sigsq2, twosigsq2,
     + A, F, G, kR,
     + ffac, gfac, rfac,
     + AF, AG, AFG, FkR, GkR, FGkR, currentEnergy,
     + dbGdihedPhase
C
      DIMENSION EPW(npa)
      DIMENSION GMUL(10)
C
      DATA GMUL/0.0d+00,2.0d+00,0.0d+00,4.0d+00,0.0d+00,6.0d+00,
     +          0.0d+00,8.0d+00,0.0d+00,10.0d+00/
      DATA TM24,TM06,tenm3/1.0d-18,1.0d-06,1.0d-03/
      data zero,one,two,four,six,twelve/0.d0,1.d0,2.d0,4.d0,6.d0,12.d0/
C
C     ---- ARRAYS GAMC = PK*COS(PHASE) AND GAMS = PK*SIN(PHASE) ----
C
      E = 0.0
      dihedE_dbG = 0.0

!	do i=1,npa

!	write(98,*) pk(i)

!	enddo 
C
C     ----- GRAND LOOP FOR THE DIHEDRAL STUFF -----
C
          DO JN = 1,nPA_dbG
            I3 = IP_dbG(JN)
            J3 = JP_dbG(JN)
            K3 = KP_dbG(JN)
            L3 = LP_dbG(JN)
            ! write(*,*) JN, nPA_dbG, I3, J3, K3, L3
C         
C           ----- CALCULATION OF ij, kj, kl VECTORS -----
C
 
	!	write(90,*) X(I3),Y(I3),Z(I3)
            XIJ(JN) = X(I3)-X(J3)
            YIJ(JN) = Y(I3)-Y(J3)
            ZIJ(JN) = Z(I3)-Z(J3)
            XKJ(JN) = X(K3)-X(J3)
            YKJ(JN) = Y(K3)-Y(J3)
            ZKJ(JN) = Z(K3)-Z(J3)
            XKL(JN) = X(K3)-X(L3)
            YKL(JN) = Y(K3)-Y(L3)
            ZKL(JN) = Z(K3)-Z(L3)                                  
          END DO
C
C         ----- GET THE NORMAL VECTOR -----
C
          ! here is the calculation of the two sets 
          ! of points - each set for each plane.
          ! why this sort of calculation?

          DO JN = 1,nPA_dbG
            DX(JN) = YIJ(JN)*ZKJ(JN)-ZIJ(JN)*YKJ(JN)
            DY(JN) = ZIJ(JN)*XKJ(JN)-XIJ(JN)*ZKJ(JN)
            DZ(JN) = XIJ(JN)*YKJ(JN)-YIJ(JN)*XKJ(JN)
            GX(JN) = ZKJ(JN)*YKL(JN)-YKJ(JN)*ZKL(JN)
            GY(JN) = XKJ(JN)*ZKL(JN)-ZKJ(JN)*XKL(JN)
            GZ(JN) = YKJ(JN)*XKL(JN)-XKJ(JN)*YKL(JN)
          END DO
C
          ! here is the following calculation:
          ! sqrt(a1^2+b1^2+c1^2) for two sets 
          ! of points - one set for each plane
          ! CT is the summ of the multiplication
          ! term of the dihedral calculation.

          DO JN = 1,nPA_dbG
            FXI(JN) = SQRT(DX(JN)*DX(JN)
     Q                    +DY(JN)*DY(JN)
     Q                    +DZ(JN)*DZ(JN)+TM24)
            FYI(JN) = SQRT(GX(JN)*GX(JN)
     Q                    +GY(JN)*GY(JN)
     Q                    +GZ(JN)*GZ(JN)+TM24)
            CT(JN) = DX(JN)*GX(JN)+DY(JN)*GY(JN)+DZ(JN)*GZ(JN)
          END DO
C
C         ----- BRANCH IF LINEAR DIHEDRAL -----
C                             
         DO JN = 1,nPA_dbG
!#ifdef CRAY_PVP
!            BIT = one/FXI(JN)
!            BIK = one/FYI(JN)
!            Z10 = CVMGT(zero,BIT,tenm3.GT.FXI(JN))
!            Z20 = CVMGT(zero,BIK,tenm3.GT.FYI(JN))
!#else
            z10 = one/FXI(jn)
            z20 = one/FYI(jn)
            if (tenm3 .gt. FXI(jn)) z10 = zero
            if (tenm3 .gt. FYI(jn)) z20 = zero
            ! what is the tenm3?
!#endif
            Z12 = Z10*Z20 ! this is the denominator
            Z1(JN) = Z10
            Z2(JN) = Z20
!#ifdef CRAY_PVP
!            FTEM = CVMGZ(zero,one,Z12)
!#else
            ftem = zero
            if (z12 .ne. zero) ftem = one
!#endif
            ! here is the numenator and denominator are combined

            FZI(JN) = FTEM
            CT0 = MIN(one,CT(JN)*Z12)
            CT1 = MAX(-one,CT0)
            S = XKJ(JN)*(DZ(JN)*GY(JN)-DY(JN)*GZ(JN))+
     Q          YKJ(JN)*(DX(JN)*GZ(JN)-DZ(JN)*GX(JN))+
     Q          ZKJ(JN)*(DY(JN)*GX(JN)-DX(JN)*GY(JN))
            AP0 = ACOS(CT1)
            AP1 = PI-SIGN(AP0,S)
            CT(JN) = AP1 ! this is the final dihedral angle
	! write(*,*) JN,AP1
            CPHI(JN) = COS(AP1)
            SPHI(JN) = SIN(AP1)
         END DO
C
C         ----- CALCULATE THE ENERGY AND THE DERIVATIVES WITH RESPECT TO
C               COSPHI -----
C
        DO JN = 1,nPA_dbG
            dbGdihedPhase = dbGPhaseDihed_array(JN)
       !      CT0 = CT(JN) + dbGdihedPhase! this is the current dihedral angle per ijkl
            CT0 = CT(JN)
            ! write(*,*) JN, CT0
            A1 = Amp1Dihed(JN) ! the epsilon value for the first well (the amplitude)
       !      mu1 = Mu1Dihed_array(JN) + dbGdihedPhase ! native angle between i-j-k-l
            mu1 = Mu1Dihed_array(JN)
            sig1 = Sig1Dihed_array(JN) ! the width of the first well
            A2 = Amp2Dihed(JN) ! the amplitude of the second well
       !      mu2 = Mu2Dihed_array(JN) + dbGdihedPhase ! native angle
            mu2 = Mu2Dihed_array(JN)
            sig2 = Sig2Dihed_array(JN) ! the width of the second well
            !rad = (RadDihed_array(JN) / 1000000) ** (12)
            ! write(*,*) CT0, A1, mu1, sig1, A2, mu2, sig2, rad
            ! write(*,*) CT0

            ! here the calculation starts
            drinv = 1.0 / CT0
            ! write(*,*) 'drinv', drinv
            bydr12 = CT0 ** (-12.0) !?
            ! write(*,*) 'bydr12', bydr12
    
            elon1 = CT0 - mu1
            elonsq1 = elon1 * elon1
            sigsq1 = sig1 * sig1
            twosigsq1 = 2.0 * sigsq1
    
            elon2 = CT0 - mu2
            elonsq2 = elon2 * elon2
            sigsq2 = sig2 * sig2
            twosigsq2 = 2.0 * sigsq2
    
            F    = exp ( - elonsq1 / twosigsq1 )
            G    = exp ( - elonsq2 / twosigsq2 )
            !kR    = rad * bydr12
    
            ! prefactors for derivative, factor (-1) missing for force
            ffac = - elon1 / sigsq1
            gfac = - elon2 / sigsq2
            rfac = 0.0000001 !-12.0 * drinv
    
            AF  = A1 * F
            AG  = A2 * G
            AFG = AF * G
            !FkR  = F * kR
            !GkR  = G * kR
            !FGkR = F * GkR
    
            currentEnergy = AFG - AF - AG! + FGkR - FkR - GkR + kR
            dihedE_dbG = E + currentEnergy
            
            DF0 = ffac * AF + gfac * AG - ( ffac + gfac ) * AFG
     Q             !+ ( ffac + rfac ) * FkR + ( gfac + rfac ) * GkR
     Q             !- ( ffac + gfac + rfac ) * FGkR - rfac * kR

             !write(*,*) DF0, CT0


! Here is the energy part
    !         E =  E + (3.0/2.0*PK(JN)-GAMC1(JN)*COSNP - 
    !  Q GAMS1(JN)*SINNP - GAMC3(JN)*COSNP3 - GAMS3(JN)*SINNP3)*FZI(JN)
    !         DF0 = -((GAMC1(JN)*SINNP-GAMS1(JN)*COSNP)
    !  Q + 3*(GAMC3(JN)*SINNP3-GAMS3(JN)*COSNP3))
! End of energy part(at least until the bottom of this routine

            DUMS = SPHI(JN)+SIGN(TM24,SPHI(JN))
            ! what is DUMS? why is the force divided by it?

!            DFLIM = GAMC(JN)*(PN(MC)-GMUL(INC)+GMUL(INC)*CPHI(JN))
! DFLIM was as is written above, but if SPhi is small ~ 0, then CPHI
! ~ 1, which means that the terms of Gmul will almost perfectly cancel.
! and all you will have left is PN, which in this case is the 1+3=4

            ! LEAVING IT FOR NOW, NEED TO UNDERSTAND IF IS NEEDED
            DFLIM = GAMC1(JN)*4.0

            df1 = (df0/CT0)
            df1 = df0/dums
            if(tm06.gt.abs(dums)) df1 = dflim
!#endif
            DF(JN) = DF1*FZI(JN)
          END DO
C                                     
C         ----- NOW DO TORSIONAL FIRST DERIVATIVES -----
C
         DO JN = 1,nPA_dbG
C
C           ----- NOW, SET UP ARRAY DC = FIRST DER. OF COSPHI W/RESPECT
C                 TO THE CARTESIAN DIFFERENCES T -----

            ! this is the tricky part to understand:
            ! how the force is calculated for ijkl in
            ! xyz directions?

            ! Why is the force calculated with respect to cosphi?
            ! should this be touched?
C
            Z11 = Z1(JN)*Z1(JN)
            Z12 = Z1(JN)*Z2(JN)
            Z22 = Z2(JN)*Z2(JN)
            DC1 = -GX(JN)*Z12-CPHI(JN)*DX(JN)*Z11
            DC2 = -GY(JN)*Z12-CPHI(JN)*DY(JN)*Z11
            DC3 = -GZ(JN)*Z12-CPHI(JN)*DZ(JN)*Z11
            DC4 =  DX(JN)*Z12+CPHI(JN)*GX(JN)*Z22
            DC5 =  DY(JN)*Z12+CPHI(JN)*GY(JN)*Z22
            DC6 =  DZ(JN)*Z12+CPHI(JN)*GZ(JN)*Z22
C
C           ----- UPDATE THE FIRST DERIVATIVE ARRAY -----
C
            DR1 = DF(JN)*(DC3*YKJ(JN) - DC2*ZKJ(JN))
            DR2 = DF(JN)*(DC1*ZKJ(JN) - DC3*XKJ(JN))
            DR3 = DF(JN)*(DC2*XKJ(JN) - DC1*YKJ(JN))
            DR4 = DF(JN)*(DC6*YKJ(JN) - DC5*ZKJ(JN))
            DR5 = DF(JN)*(DC4*ZKJ(JN) - DC6*XKJ(JN))
            DR6 = DF(JN)*(DC5*XKJ(JN) - DC4*YKJ(JN))
            DRX = DF(JN)*(-DC2*ZIJ(JN) + DC3*YIJ(JN) +
     +               DC5*ZKL(JN) - DC6*YKL(JN))
            DRY = DF(JN)*( DC1*ZIJ(JN) - DC3*XIJ(JN) -
     +               DC4*ZKL(JN) + DC6*XKL(JN))
            DRZ = DF(JN)*(-DC1*YIJ(JN) + DC2*XIJ(JN) +
     +               DC4*YKL(JN) - DC5*XKL(JN))
 
            I3 = IP_dbG(JN)
            J3 = JP_dbG(JN)
            K3 = KP_dbG(JN)
            L3 = LP_dbG(JN)


            FX(I3) = FX(I3) -  DR1
            FY(I3) = FY(I3) -  DR2
            FZ(I3) = FZ(I3) -  DR3
            FX(J3) = FX(J3) -  DRX +  DR1
            FY(J3) = FY(J3) -  DRY +  DR2
            FZ(J3) = FZ(J3) -  DRZ +  DR3
            FX(K3) = FX(K3) +  DRX +  DR4
            FY(K3) = FY(K3) +  DRY +  DR5
            FZ(K3) = FZ(K3) +  DRZ +  DR6
            FX(L3) = FX(L3) -  DR4
            FY(L3) = FY(L3) -  DR5
            FZ(L3) = FZ(L3) -  DR6

          END DO




          END


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^END of Dihedral^^^^^^^^^^^^^^^^^^^^^^^^^^^
