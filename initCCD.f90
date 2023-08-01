!     This is the initialization routine.
!     It sets up several parameters, constants, and initial guess values.

Subroutine initCCD (ccd)

  include 'values.h'
  include 'global.h'
  include 'defs.h'

  type (membrane) :: ccd(0:NZ)
  double precision theta(NC),Slum(NC),Slat(NC),Sbas(NC)
  double precision Pf(NC,NC),dLA(NS,NS,NC,NC)
  double precision pos(0:NZ)
  
!---------------------------------------------------------------------72      
!---------------------------------------------------------------------72
!   Solute Indices:  1 = Na+, 2 = K+, 3 = Cl-, 4 = HCO3-, 5 = H2CO3, 6 = CO2
!	7 = HPO4(2-), 8 = H2PO4-, 9 = urea, 10 = NH3, 11 = NH4+, 12 = H+
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


!---------------------------------------------------------------------72
!	 Data for rat (AMW model, AJP Renal 2005)
!---------------------------------------------------------------------72

  ! tubular length and luminal diameter
  ccd(:)%dimL = dimLCCD
  ccd(:)%diam = DiamCCD

  SlumPinitccd = 1.2d0
  SbasPinitccd = 1.2d0
  SlatPinitccd = 6.9d0

  SlumAinitccd = 0.58d0
  SbasAinitccd = 0.58d0
  SlatAinitccd = 1.25d0
  
  SlumBinitccd = 0.17d0
  SbasBinitccd = 0.17d0
  SlatBinitccd = 1.50d0
  
  SlumEinitccd = 0.001d0
  ccd(:)%sbasEinit = 0.020d0
  
  ccd(:)%volPinit = 4.0d0 !6.0d0
  ccd(:)%volAinit = 1.8d0 !2.7d0
  ccd(:)%volBinit = 1.2d0 !1.8d0
  ccd(:)%volEinit = 0.20d0

!---------------------------------------------------------------------72
!	 Assign membrane surface area except for basal membrane of E
!---------------------------------------------------------------------72

  ccd(:)%area(1,2)=SlumPinitccd
  ccd(:)%area(1,3)=SlumAinitccd
  ccd(:)%area(1,4)=SlumBinitccd
  ccd(:)%area(1,5)=SlumEinitccd
	 
  ccd(:)%area(2,5)=SlatPinitccd
  ccd(:)%area(3,5)=SlatAinitccd
  ccd(:)%area(4,5)=SlatBinitccd
  
  ccd(:)%area(2,6)=SbasPinitccd
  ccd(:)%area(3,6)=SbasAinitccd
  ccd(:)%area(4,6)=SbasBinitccd
  
  Do K = 1,NC-1
     Do L = K+1,NC
        ccd(:)%area(L,K)=ccd(:)%area(K,L)
     End do
  End do


!---------------------------------------------------------------------72      
!---------------------------------------------------------------------72
!   Solute valence
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
					
  zval(1) = +1.0d0
  zval(2) = +1.0d0
  zval(3) = -1.0d0
  zval(4) = -1.0d0
  zval(5) = 0.0d0
  zval(6) = 0.0d0
  zval(7) = -2.0d0
  zval(8) = -1.0d0
  zval(9) = 0.0d0
  zval(10) = 0.0d0
  zval(11) = +1.0d0
  zval(12) = +1.0d0
  zval(13) = -1.0d0
  zval(14) = 0.0d0 
  zval(15) = 0.0d0
  zval(16) = +2.0d0

  ! impermeant valence
  ccd(:)%zPimp = -1.0d0
  ccd(:)%zAimp = -1.0d0
  ccd(:)%zBimp = -1.0d0


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Water permeabilities
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	The permeabilities have already been multiplied by the area (cm2/cm2 epith)

!   dLPV(K,L) = hydraulic permeability at the K-L interface 
!	dLPV(K,L) in cm/s/mmHg
!	Pf(K,L) = osmotic permeability at the K-L interface
!	Pf(K,L) in cm/s
!	Relationship between the two: LPV=Pf*(Vwbar/RT)



!  WATER PERMEABILITIES
  !!! - FEMALE ADJUSTMENT Pf x 2.4
  PfMP = 2.40*0.20d0/1.2d0
  PfPE = 2.40*0.11d0*3
  PfPS = 2.40*0.11d0*3

  PfMA = 0.22d-3
  PfAE = 5.50d-3
  PfAS = 5.50d-3

  PfMB = 0.22d-3
  PfBE = 5.50d-3
  PfBS = 5.50d-3

  PfME = 1.1d0
  PfES = 110.d0

  if (FAngAqp2) then
     PfMP = PfMP*0.46 ! AQP2-cor x 0.46-0.47
     PfPE = PfPE*0.46 ! AQP2-cor x 0.46-0.47
     PfPS = PfPS*0.46 ! AQP2-cor x 0.46-0.47
  end if

  Do K = 1, NC
     Do L = 1, NC
        Pf(K,L)=0.0d0
     End Do
  End Do

  Pf(1,2) = PfMP
  Pf(1,3) = PfMA
  Pf(1,4) = PfMB
  Pf(1,5) = PfME
  Pf(2,5) = PfPE
  Pf(3,5) = PfAE
  Pf(4,5) = PfBE
  Pf(2,6) = PfPS
  Pf(3,6) = PfAS
  Pf(4,6) = PfBS
  Pf(5,6) = PfES
  
!	  Units of dimensional water flux: cm3/s/cm2 epith
!	  Non-dimensional factor for water flux: (Pfref)*Vwbar*Cref
!	  Calculate non-dimensional dLPV = Pf*Vwbar*Cref / (Pfref*Vwbar*Cref)
!	  dLPV = Pf/Pfref 

  Do K = 1, NC
     Do L = 1, NC
        ccd(:)%dLPV(K,L) = Pf(K,L)/(Pfref)
     End Do
  End Do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     sig(I,K,L) = reflection coefficient of Ith solute between K and L
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	Initialize 
  Do I = 1,NS
     Do K = 1, NC
        Do L = 1, NC
           ccd(:)%sig(I,K,L) = 1.0d0
        End Do
     End Do
  End Do


!	Values at the basement membrane (ES)
  Do I = 1,NS	
     ccd(:)%sig(I,5,6) = 0.0d0
  End Do


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Dimensional Solute permeabilities
!     h(I,K,L) = permeability of ith solute at the K-L interface
!	h(I,K,L) in 10-5 cm/s initially
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


  PICAchlo = 10.0d0
  PICBchlo = 3.20d0
  PPCh2co3 = 130.d0
  PICh2co3 = 10.d0
  PPCco2 = 15.0d3
  PICco2 = 900.d0
  PPChpo4 = 8.00d-3
  PICAhpo4 = 4.80d-3
  PICBhpo4 = 4.80d-3
  PPCprot = 2000.d0
  PICAprot = 9.0d0
  PICBprot = 6.0d0
  PPCamon = 2000.d0
  PICamon = 900.d0
  Purea = 0.10d0
  

!	Initialize 
  Do I = 1,NS
     Do K = 1, NC
        Do L = K, NC
           ccd(:)%h(I,K,L) = 0.0d0
        End Do
     End Do
  End Do
  
  fac4 = 1.0d0/4.0d0

!	Values at the MP interface (1-2)

! ENaC activity modified by local factors, including pH. Use hENaC_CCD as basal value of expression
!!! - FEMALE ADJUSTMENT ENaC
  hENaC_CCD= FM_cENaC*35.0*fac4

  if (FAngDist) then
    hENaC_CCD = hENaC_CCD*2.90 ! cortical ENaC x 1.04 to 2.97 in grain-fed
  end if


! ROMK activity modified by local factors, including pH. Use hROMK_CCD as basal value of expression
  hROMK_CCD = 8.0d0*fac4

!  if (FAngDist) then
!    hROMK_CCD = hROMK_CCD*0.70 ! cortical ROMK x 0.52 to 0.74 in grain-fed
!  end if


  ccd(:)%h(3,1,2) = 0.80d0*fac4 
  ccd(:)%h(4,1,2) = 0.16d0*fac4 
  ccd(:)%h(5,1,2) = PPCh2co3 
  ccd(:)%h(6,1,2) = PPCco2
  ccd(:)%h(7,1,2) = 1.0d-3/1.20
  ccd(:)%h(8,1,2) = 1.0d-3/1.20
  ccd(:)%h(9,1,2) = Purea
  ccd(:)%h(10,1,2) = PPCamon  
  ccd(:)%h(11,1,2) = ccd(:)%h(2,1,2)*0.20d0
  ccd(:)%h(12,1,2) = PPCprot*fac4 
  ccd(:)%h(13,1,2) = 0.0d0
  ccd(:)%h(14,1,2) = 0.0d0 
  ccd(:)%h(15,1,2) = 0.0d0
  ccd(:)%h(16,1,2) = 0.00010d0
  
!	Values at the PE interface (2-5)
  ccd(:)%h(1,2,5) = 0.000d0
  ccd(:)%h(2,2,5) = 4.0d0*fac4
  ccd(:)%h(3,2,5) = 0.20d0*fac4
  ccd(:)%h(4,2,5) = ccd(:)%h(3,2,5)*0.20d0
  ccd(:)%h(5,2,5) = PPCh2co3
  ccd(:)%h(6,2,5) = PPCco2
  ccd(:)%h(7,2,5) = 8.00d-3
  ccd(:)%h(8,2,5) = 8.00d-3
  ccd(:)%h(9,2,5) = Purea
  ccd(:)%h(10,2,5) = PPCamon 
  ccd(:)%h(11,2,5) = ccd(:)%h(2,2,5)*0.20d0
  ccd(:)%h(12,2,5) = PPCprot*fac4 
  ccd(:)%h(13,2,5) = 1.0d-4
  ccd(:)%h(14,2,5) = 1.0d-4 
  ccd(:)%h(15,2,5) = 1.0d-4
  ccd(:)%h(16,2,5) = 0.00010d0

!	Values at the MA interface (1-3)
  ccd(:)%h(1,1,3) = 0.0d0
  ccd(:)%h(2,1,3) = 0.00d0  
  ccd(:)%h(3,1,3) = 0.0d0
  ccd(:)%h(4,1,3) = 0.0d0
  ccd(:)%h(5,1,3) = PICh2co3
  ccd(:)%h(6,1,3) = PICco2
  ccd(:)%h(7,1,3) = 0.0d0
  ccd(:)%h(8,1,3) = 0.0d0
  ccd(:)%h(9,1,3) = Purea
  ccd(:)%h(10,1,3) = PICamon
  ccd(:)%h(11,1,3) = 0.00d0 
  ccd(:)%h(12,1,3) = 0.0d0
  ccd(:)%h(13,1,3) = 0.0d0
  ccd(:)%h(14,1,3) = 0.0d0
  ccd(:)%h(15,1,3) = 0.0d0
  ccd(:)%h(16,1,3) = 0.00010d0 ! TO BE ADJUSTED
  
!	Values at the AE interface (3-5)
  ccd(:)%hCLCA=PICAchlo !!! Different in AMW model
  ccd(:)%h(1,3,5) = 0.0d0
  ccd(:)%h(2,3,5) = 0.050d0 
  ccd(:)%h(3,3,5) = 1.20d0 
  ccd(:)%h(4,3,5) = 0.18d0 
  ccd(:)%h(5,3,5) = PICh2co3
  ccd(:)%h(6,3,5) = PICco2
  ccd(:)%h(7,3,5) = 0.0120
  ccd(:)%h(8,3,5) = 0.0120
  ccd(:)%h(9,3,5) = Purea
  ccd(:)%h(10,3,5) = PICamon
  ccd(:)%h(11,3,5) = 0.030d0
  ccd(:)%h(12,3,5) = 1.50d0
  ccd(:)%h(13,3,5) = 1.0d-4
  ccd(:)%h(14,3,5) = 1.0d-4
  ccd(:)%h(15,3,5) = 1.0d-4
  ccd(:)%h(16,3,5) = 0.00010d0 ! TO BE ADJUSTED

!	Values at the MB interface (1-4)
  ccd(:)%h(1,1,4) = 0.0d0
  ccd(:)%h(2,1,4) = 0.0d0
  ccd(:)%h(3,1,4) = 0.0d0
  ccd(:)%h(4,1,4) = 0.0d0
  ccd(:)%h(5,1,4) = PICh2co3
  ccd(:)%h(6,1,4) = PICco2
  ccd(:)%h(7,1,4) = 0.0d0
  ccd(:)%h(8,1,4) = 0.0d0
  ccd(:)%h(9,1,4) = Purea
  ccd(:)%h(10,1,4) = PICamon
  ccd(:)%h(11,1,4) = 0.0d0 
  ccd(:)%h(12,1,4) = 0.0d0
  ccd(:)%h(13,1,4) = 0.0d0
  ccd(:)%h(14,1,4) = 0.0d0 
  ccd(:)%h(15,1,4) = 0.0d0
  ccd(:)%h(16,1,4) = 0.00010d0
  
!	Values at the BE interface (4-5)
  ccd(:)%hCLCB=PICBchlo
  ccd(:)%h(1,4,5) = 0.0d0
  ccd(:)%h(2,4,5) = 0.12d0*fac4 
  ccd(:)%h(3,4,5) = PICBchlo*fac4 
  ccd(:)%h(4,4,5) = 0.40d0*fac4 
  ccd(:)%h(5,4,5) = PICh2co3
  ccd(:)%h(6,4,5) = PICco2
  ccd(:)%h(7,4,5) = PICBhpo4*fac4
  ccd(:)%h(8,4,5) = PICBhpo4*fac4
  ccd(:)%h(7,4,5) = 0.120
  ccd(:)%h(8,4,5) = 0.120
  ccd(:)%h(9,4,5) = Purea
  ccd(:)%h(10,4,5) = PICamon
  ccd(:)%h(11,4,5) = 0.12d0*fac4 
  ccd(:)%h(12,4,5) = PICBprot*fac4 
  ccd(:)%h(13,4,5) = 1.0d-4
  ccd(:)%h(14,4,5) = 1.0d-4 
  ccd(:)%h(15,4,5) = 1.0d-4 
  ccd(:)%h(16,4,5) = 0.00010d0
  
!	Values at basal interfaces (2-6), (3-6), (4-6)

  Do I = 1, NS
     ccd(:)%h(I,2,6) = ccd(:)%h(I,2,5)
     ccd(:)%h(I,3,6) = ccd(:)%h(I,3,5)
     ccd(:)%h(I,4,6) = ccd(:)%h(I,4,5)
  End do


!	Values at the ME interface (1-5)
!	Account for factor 2 increase (DN model, 2008)
  ClTJperm = 1000.0d0
!!! - FEMALE ADJUSTMENT related to Claudin expression
  ccd(:)%h(1,1,5) = (1.00/FM_cENaC)*ClTJperm*1.0d0 ! Cldn8 and ENaC
  ccd(:)%h(2,1,5) = (1.00/FM_cENaC)*ClTJperm*1.2d0 ! Cldn8 and ENaC
! Cl permeability modified by local factors, including pH. Use hCltj_CNT as basal value of expression
  hCltj_CCD = 1.30*ClTJperm*1.2d0 ! Cldn7
  ccd(:)%h(4,1,5) = ClTJperm*0.30d0
  ccd(:)%h(5,1,5) = ClTJperm*2.0d0
  ccd(:)%h(6,1,5) = ClTJperm*2.0d0
  ccd(:)%h(7,1,5) = ClTJperm*0.14d0
  ccd(:)%h(8,1,5) = ClTJperm*0.14d0
  ccd(:)%h(9,1,5) = ClTJperm*0.40d0
  ccd(:)%h(10,1,5) = ClTJperm*6.0d0
  ccd(:)%h(11,1,5) = ClTJperm*1.2d0
  ccd(:)%h(12,1,5) = ClTJperm*10.d0
  ccd(:)%h(13,1,5) = ClTJperm*0.01d0
  ccd(:)%h(14,1,5) = ClTJperm*0.01d0
  ccd(:)%h(15,1,5) = ClTJperm*0.01d0
  ccd(:)%h(16,1,5) = 0.0140d0*ClTJperm ! See my 2015 AJP and Carney 1988

!  if (FAngDist) then
!    ccd(:)%h(1,1,5) = ccd(:)%h(1,1,5)/2.0 ! Cldn8 and ENaC
!    ccd(:)%h(2,1,5) = ccd(:)%h(2,1,5)/2.0 ! Cldn8 and ENaC
!  end if

!	Values at the ES interface (5-6)
  areafactor = 100.00
  ccd(:)%h(1,5,6) = 97.d0*areafactor
  ccd(:)%h(2,5,6) = 130.d0*areafactor
  ccd(:)%h(3,5,6) = 130.d0*areafactor
  ccd(:)%h(4,5,6) = 65.d0*areafactor
  ccd(:)%h(5,5,6) = 97.d0*areafactor
  ccd(:)%h(6,5,6) = 97.d0*areafactor
  ccd(:)%h(7,5,6) = 65.d0*areafactor
  ccd(:)%h(8,5,6) = 65.d0*areafactor
  ccd(:)%h(9,5,6) = 65.d0*areafactor
  ccd(:)%h(10,5,6) = 130.d0*areafactor
  ccd(:)%h(11,5,6) = 130.0d0*areafactor
  ccd(:)%h(12,5,6) = 6490.0d0*areafactor
  ccd(:)%h(13,5,6) = 4.0d0*areafactor
  ccd(:)%h(14,5,6) = 4.0d0*areafactor
  ccd(:)%h(15,5,6) = 4.0d0*areafactor
  ccd(:)%h(16,5,6) = ccd(:)%h(1,5,6)*(7.93/13.3)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Non-Dimensional Solute permeabilities
!     Divide h(I,K,L) by (href) 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  Do I = 1,NS
     Do K = 1, NC
        Do L = K, NC
           ccd(:)%h(I,K,L) = ccd(:)%h(I,K,L)*1.0d-5/(href)
        End Do
     End Do
  End Do

	
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Dimensional NET coefficients
!     dLA(I,J,K,L) = coefficient of ith and jth solute at the K-L interface
!	dLA(I,J,K,L) in mmol2/J/s initially
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Initialize 
  Do I = 1,NS
     Do J = 1, NS
        Do K = 1, NC
           Do L = K, NC
              dLA(I,J,K,L) = 0.0d0
           End Do
        End Do
     End Do
  End Do
  
!	 Na/H exchangers on all basolateral membranes	 
  dLA(1,12,2,5)=20.0d-9*fac4 
  dLA(1,12,2,6)=20.0d-9*fac4
  dLA(1,12,3,5)=24.0d-9*fac4
  dLA(1,12,3,6)=24.0d-9*fac4
  dLA(1,12,4,5)=24.0d-9*fac4
  dLA(1,12,4,6)=24.0d-9*fac4
  
!	 Na2/HPO4 co-transporters on all basolateral membranes	 
  dLA(1,7,2,5)=2.0d-9
  dLA(1,7,2,6)=2.0d-9
  dLA(1,7,3,5)=0.80d-9
  dLA(1,7,3,6)=0.80d-9
  dLA(1,7,4,5)=0.80d-9
  dLA(1,7,4,6)=0.80d-9
  
!	 Basolateral Cl/HCO3 exchanger in PC 
  dLA(3,4,2,5)=2.0d-9*fac4
  dLA(3,4,2,6)=2.0d-9*fac4
  
!	 AE1 in IC-A only
  xAE1 = 15.d-9*0.20d0

!	 Transporters in IC-B only
!!! - FEMALE ADJUSTMENT Pendrin x 1.20
  dLA(3,4,1,4)=1.20*800.d-9*fac4
  xPendrin = dLA(3,4,1,4)/1700.d0

  if (FAngDist) then
    xPendrin = xPendrin*1.47 ! Pendrin x 1.47
  end if

  xNDBCE = 0.00 !NDCBE

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  Do I = 1, NS-1
     Do J = I+1, NS
        Do K = 1, NC
           Do L = 1, NC
              dLA(J,I,K,L) = dLA(I,J,K,L)
           End Do
        End Do
     End Do
  End Do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Dimensional ATPase coefficients
!     ATPcoeff in mmol/s or mmol
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!		Na-K-ATPase (basolateral in PC and IC)
!!! - FEMALE ADJUSTMENT ENaC and NKA
 ! ATPNaKPES = FM_cENaC*3410.d-9*fac4*fscaleNaK_PC
  ATPNaKPES = FM_cENaC*3410.d-9*1.25*fac4
  ATPNaKAES = 300.d-9*fac4
  ATPNaKBES = 300.d-9*fac4 

!		H-ATPase (apical in IC-A, basolateral in IC-B)
  ATPHMA = 1000.0d-9
  ATPHBES = 400.d-9*15*fac4
!         Note: Adjusted parameter

!		H-K-ATPase (apical in PC and IC)
  ATPHKMP = 0.d0
  ATPHKMA = 60.0d-9
  ATPHKMB = 0.0d-9
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	  Non-Dimensional coefficients
!       Divide dLA(I,J,K,L) by (href)*Cref for consistency with other
!	  solute flux terms
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  Do I = 1,NS
     Do J = 1, NS
        Do K = 1, NC
           Do L = K, NC
              ccd(:)%dLA(I,J,K,L) = dLA(I,J,K,L)/(href*Cref)
           End Do
        End Do
     End Do
  End Do

  ccd(:)%xPendrin = xPendrin/(href*Cref)
  ccd(:)%xAE1 = xAE1/(href*Cref)
  ccd(:)%xNDBCE = xNDBCE/(href*Cref)
  
  ccd(:)%ATPNaK(2,5) = ATPNaKPES/(href*Cref)
  ccd(:)%ATPNaK(2,6) = ATPNaKPES/(href*Cref)
  ccd(:)%ATPNaK(3,5) = ATPNaKAES/(href*Cref)
  ccd(:)%ATPNaK(3,6) = ATPNaKAES/(href*Cref)
  ccd(:)%ATPNaK(4,5) = ATPNaKBES/(href*Cref)
  ccd(:)%ATPNaK(4,6) = ATPNaKBES/(href*Cref)
  
  ccd(:)%ATPH(1,3) = ATPHMA/(href*Cref)
  ccd(:)%ATPH(4,5) = ATPHBES/(href*Cref)
  ccd(:)%ATPH(4,6) = ATPHBES/(href*Cref)
  
  ccd(:)%ATPHK(1,2) = ATPHKMP/(href*Cref)
  ccd(:)%ATPHK(1,3) = ATPHKMA/(href*Cref)
  ccd(:)%ATPHK(1,4) = ATPHKMB/(href*Cref)
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     Other kinetic parameters
!	Units of kinetic constants are s-1
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
! KINETIC PARAMETERS
       
  dkhuncat = 0.145d0
  dkduncat = 49.6d0
  dkhP = 1.45d0
  dkdP = 496.0d0
  dkhA = 1.45d3
  dkdA = 496.0d3
  dkhB = 1.45d3
  dkdB = 496.0d3
  dkhE = 0.145d3
  dkdE = 49.6d3
  
  ccd(:)%dkd(1)=dkduncat*10.0d0
  ccd(:)%dkh(1)=dkhuncat*10.0d0
  ccd(:)%dkd(2)=dkdP
  ccd(:)%dkh(2)=dkhP
  ccd(:)%dkd(3)=dkdA
  ccd(:)%dkh(3)=dkhA
  ccd(:)%dkd(4)=dkdB
  ccd(:)%dkh(4)=dkhB
  ccd(:)%dkd(5)=dkdE
  ccd(:)%dkh(5)=dkhE

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	READ ENTERING AND PERITUBULAR CONDITIONS FROM CNT OUTLET FILE
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


  open ( unit=42, file='CNToutlet' )
  
  Do I = 1,4
     read(42,200),ccd(0)%conc(I,1),ccd(0)%conc(I,6)
  end Do
  Do I = 5,NS
     read(42,210),ccd(0)%conc(I,1),ccd(0)%conc(I,6)
  end Do
  read(42,200),ccd(0)%ph(1),ccd(0)%ph(6)
  read(42,200),ccd(0)%ep(1),ccd(0)%ep(6)
  read(42,210),ccd(0)%vol(1),ccd(0)%pres
  close ( unit=42 )
  ccd(:)%volLuminit = ccd(0)%vol(1)

  ccd(:)%ep(6)=ccd(0)%ep(6)

   Do jz=0,NZ
      pos(jz) = 1.0d0*jz/NZ
   End do

  call set_intconc ( ccd, NZ, 1, pos )

  ! reference impermeant concentrations
  ccd(:)%cPimpref = 50.0d0
  ccd(:)%cAimpref = 18.0d0
  ccd(:)%cBimpref = 18.0d0

  ! total buffer concentrations
  ccd(:)%cPbuftot = 32.0d0
  ccd(:)%cAbuftot = 40.0d0
  ccd(:)%cBbuftot = 40.0d0

  !  coalescence parameter
  ccd(:)%coalesce = 0.2d0
  
!---------------------------------------------------------------------72
!     For metabolic calculations
!---------------------------------------------------------------------72

  if (ndiabetes .eq. 0)  then
    ccd(:)%TQ = 15.0 !TNa-QO2 ratio in normal CCD
  else
    ccd(:)%TQ = 12.0 !TNa-QO2 ratio in diabetic CCD
  end if

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

 200   format (6f12.5)
 210   format (6e12.5)

  return
end Subroutine initCCD
      
