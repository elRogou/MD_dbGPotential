!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* CONTACTS: computes the force on all atoms due to contacts via a     *
!* 10-12 potential                                                     *
!***********************************************************************

      subroutine contacts(E, E_dbG, Conts, count, outE)
      include 'MD.com'

      integer count(NC), outE
      real E, E_dbG
      integer Conts

      ! Initialize counts and energies
      E = 0.0
      E_dbG = 0.0
      Conts = 0
      do i = 1, NC
        count(i) = 0
      end do

      ! Call the Lennard-Jones part
      call calculate_LJ(E, Conts, count, outE)

      ! Call the double Gaussian part
      if (isdbGaussian .eq. 'YES') then
        call calculate_dbG(E_dbG)
      end if

      end subroutine contacts
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      subroutine calculate_LJ(E, Conts, count, outE)
      include 'MD.com'

      integer C1, C2, count(NC), Conts, isContact, currentRangesIndex(50), currentContactInRange, outE
      real r2, rm2, rm10, currentEnergy, f_over_r
      integer i, j

      ! Initialize ranges if needed
      if (writeContactsRanges .eq. 'YES') then
        do i=1, rangesNumber
          currentRangesIndex(i) = 1
          EtotalRanges(i) = 0.0
          TwoBodyRanges(i) = 0
        end do
      endif

      do i = 1, NC
        C1 = IC(i)
        C2 = JC(i)

        dx = X(C1) - X(C2)
        dy = Y(C1) - Y(C2)
        dz = Z(C1) - Z(C2)

        r2 = dx**2 + dy**2 + dz**2

        rm2 = 1.0/r2
        rm2 = rm2 * sigma(i)
        rm10 = rm2**5

        currentEnergy = epsC(i) * rm10 * (5*rm2 - 6)
        f_over_r = epsC(i) * rm10 * (rm2 - 1) * 60 / r2

        E = E + currentEnergy

        if (r2 .le. sigma(i)*ConCut**2) then
          Conts = Conts + 1
          count(Conts) = i
          isContact = 1
        else
          isContact = 0
        endif

        if ((writeAllContacts .eq. 'YES') .and. (outE .eq. 0)) then
          write(14,*) isContact
        endif

        if (writeContactsRanges .eq. 'YES') then
          do j=1, rangesNumber
            currentContactInRange = rangeContacts(j,currentRangesIndex(j))
            if (i .eq. currentContactInRange) then
              EtotalRanges(j) = EtotalRanges(j) + currentEnergy
              TwoBodyRanges(j) = TwoBodyRanges(j) + isContact
              currentRangesIndex(j) = currentRangesIndex(j) + 1
            endif
          enddo
        endif

        ! Forces
        Fx(C1) = Fx(C1) + f_over_r * dx
        Fy(C1) = Fy(C1) + f_over_r * dy
        Fz(C1) = Fz(C1) + f_over_r * dz

        Fx(C2) = Fx(C2) - f_over_r * dx
        Fy(C2) = Fy(C2) - f_over_r * dy
        Fz(C2) = Fz(C2) - f_over_r * dz

      end do

      if ((writeContactsRanges .eq. 'YES') .and. (outE .eq. 0)) then
        do j=1, rangesNumber
          write (15,'(I5)') TwoBodyRanges(j)
          write (16,'(F9.2)') EtotalRanges(j)
        enddo
        write (15,*) ''
        write (16,*) ''
      endif

      end subroutine calculate_LJ
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      subroutine calculate_dbG(E_dbG)
      include 'MD.com'

      integer C1, C2
      real r2, dr
      real invDist_dbG, invDist12_dbG
      real deltaMu1_dbG, deltaMu1Sq_dbG, sigma1Sq_dbG, twoSigma1Sq_dbG
      real deltaMu2_dbG, deltaMu2Sq_dbG, sigma2Sq_dbG, twoSigma2Sq_dbG
      real gauss1_dbG, gauss2_dbG, repulsionTerm_dbG
      real prefactorF_dbG, prefactorG_dbG, repulsionPrefactor_dbG
      real A_times_F_dbG, A_times_G_dbG, AF_times_G_dbG
      real F_times_kR_dbG, G_times_kR_dbG, FG_times_kR_dbG
      real currentEnergy, f_over_r, rad
      real A1, A2, mu1, sig1, mu2, sig2
      integer i

      do i = 1, NC_dbG
        C1 = IC_dbG(i)
        C2 = JC_dbG(i)

        dx = X(C1) - X(C2)
        dy = Y(C1) - Y(C2)
        dz = Z(C1) - Z(C2)

        r2 = dx**2 + dy**2 + dz**2
        dr = sqrt(r2)

        A1 = Amp1(i)
        mu1 = Mu1_array(i)
        sig1 = Sig1_array(i)
        A2 = Amp2(i)
        mu2 = Mu2_array(i)
        sig2 = Sig2_array(i)
        rad = Rad_array(i)**12

        invDist_dbG = 1.0 / dr
        invDist12_dbG = dr**(-12.0)

        deltaMu1_dbG = dr - mu1
        deltaMu1Sq_dbG = deltaMu1_dbG**2
        sigma1Sq_dbG = sig1**2
        twoSigma1Sq_dbG = 2.0 * sigma1Sq_dbG

        deltaMu2_dbG = dr - mu2
        deltaMu2Sq_dbG = deltaMu2_dbG**2
        sigma2Sq_dbG = sig2**2
        twoSigma2Sq_dbG = 2.0 * sigma2Sq_dbG

        gauss1_dbG = exp(-deltaMu1Sq_dbG / twoSigma1Sq_dbG)
        gauss2_dbG = exp(-deltaMu2Sq_dbG / twoSigma2Sq_dbG)
        repulsionTerm_dbG = rad * invDist12_dbG

        prefactorF_dbG = -deltaMu1_dbG / sigma1Sq_dbG
        prefactorG_dbG = -deltaMu2_dbG / sigma2Sq_dbG
        repulsionPrefactor_dbG = -12.0 * invDist_dbG

        A_times_F_dbG = A1 * gauss1_dbG
        A_times_G_dbG = A2 * gauss2_dbG
        AF_times_G_dbG = A_times_F_dbG * gauss2_dbG
        F_times_kR_dbG = gauss1_dbG * repulsionTerm_dbG
        G_times_kR_dbG = gauss2_dbG * repulsionTerm_dbG
        FG_times_kR_dbG = gauss1_dbG * G_times_kR_dbG

        currentEnergy = AF_times_G_dbG - A_times_F_dbG - A_times_G_dbG
       &              + FG_times_kR_dbG - F_times_kR_dbG - G_times_kR_dbG
       &              + repulsionTerm_dbG
        E_dbG = E_dbG + currentEnergy

        f_over_r = prefactorF_dbG * A_times_F_dbG
       &         + prefactorG_dbG * A_times_G_dbG
       &         - (prefactorF_dbG + prefactorG_dbG)*AF_times_G_dbG
       &         + (prefactorF_dbG + repulsionPrefactor_dbG)*F_times_kR_dbG
       &         + (prefactorG_dbG + repulsionPrefactor_dbG)*G_times_kR_dbG
       &         - (prefactorF_dbG + prefactorG_dbG + repulsionPrefactor_dbG)*FG_times_kR_dbG
       &         - repulsionPrefactor_dbG * repulsionTerm_dbG
        f_over_r = f_over_r / dr

        Fx(C1) = Fx(C1) + f_over_r * dx
        Fy(C1) = Fy(C1) + f_over_r * dy
        Fz(C1) = Fz(C1) + f_over_r * dz

        Fx(C2) = Fx(C2) - f_over_r * dx
        Fy(C2) = Fy(C2) - f_over_r * dy
        Fz(C2) = Fz(C2) - f_over_r * dz

      end do

      end subroutine calculate_dbG
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<