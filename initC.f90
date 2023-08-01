!     This is the initialization routine.
!     It sets up several parameters, constants, and initial guess values.

Subroutine initC (cnt)

  include 'values.h'
  include 'global.h'
  include 'defs.h'

  type (membrane) :: cnt(0:NZ)
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
  cnt(:)%dimL = dimLC
  cnt(:)%diam = DiamC


  SlumPinitcnt = 1.2d0
  SbasPinitcnt = 1.2d0
  SlatPinitcnt = 6.9d0
  
  SlumAinitcnt = 0.58d0
  SbasAinitcnt = 0.58d0
  SlatAinitcnt = 1.25d0
  
  SlumBinitcnt = 0.17d0
  SbasBinitcnt = 0.17d0
  SlatBinitcnt = 1.50d0
  
  SlumEinitcnt = 0.001d0
  cnt(:)%sbasEinit = 0.020d0
  
  cnt(:)%volPinit = 6.0d0
  cnt(:)%volAinit = 2.7d0
  cnt(:)%volBinit = 1.8d0
  cnt(:)%volEinit = 0.20d0

!---------------------------------------------------------------------72
!	 Assign membrane surface area except for basal membrane of E
!---------------------------------------------------------------------72

  cnt(:)%area(1,2)=SlumPinitcnt
  cnt(:)%area(1,3)=SlumAinitcnt
  cnt(:)%area(1,4)=SlumBinitcnt
  cnt(:)%area(1,5)=SlumEinitcnt
  
  cnt(:)%area(2,5)=SlatPinitcnt
  cnt(:)%area(3,5)=SlatAinitcnt
  cnt(:)%area(4,5)=SlatBinitcnt
  
  cnt(:)%area(2,6)=SbasPinitcnt
  cnt(:)%area(3,6)=SbasAinitcnt
  cnt(:)%area(4,6)=SbasBinitcnt
  
  Do K = 1,NC-1
     Do L = K+1,NC
        cnt(:)%area(L,K)=cnt(:)%area(K,L)
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
  cnt(:)%zPimp = -1.0d0
  cnt(:)%zAimp = -1.0d0
  cnt(:)%zBimp = -1.0d0

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Water permeabilities
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!!! - FEMALE ADJUSTMENT Pf x 2.4
  PfMP = 2.40*0.036d0
  PfPE = 2.40*0.44d0
  PfPS = 2.40*0.44d0

  PfMA = 0.22d-3
  PfAE = 5.50d-3
  PfAS = 5.50d-3

  PfMB = 0.22d-3
  PfBE = 5.50d-3
  PfBS = 5.50d-3

  PfES = 110.d0
  PfME = 1.1d0
  
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
        cnt(:)%dLPV(K,L) = Pf(K,L)/(Pfref)
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
           cnt(:)%sig(I,K,L) = 1.0d0
        End Do
     End Do
  End Do
  

!	Values at the basement membrane (ES)
  Do I = 1,NS	
     cnt(:)%sig(I,5,6) = 0.0d0
  End Do
  

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Dimensional Solute permeabilities
!     h(I,K,L) = permeability of ith solute at the K-L interface
!	h(I,K,L) in 10-5 cm/s initially
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  PICAchlo = 11.0d0
  PICBchlo = 3.20d0
  PPCh2co3 = 130.d0
  PICh2co3 = 10.d0
  PPCco2 = 15.0d3
  PICco2 = 900.d0
  PPChpo4 = 8.00d-3
  PICAhpo4 = 7.20d-3
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
           cnt(:)%h(I,K,L) = 0.0d0
        End Do
     End Do
  End Do

! Original values for ENaC expression in CNT, CCD, OMCD and IMCD multiplied by fscaleENaC
! Original values for Na/K-ATPase expression in principal cells in CNT through IMCD multiplied by fscaleNaK
fscaleENaC = 2.50
fscaleNaK = 1.25
fscaleNaK_PC = 1.25

! Female-to-male ENaC expression ratio in cortex (CNT and CCD)
FM_cENaC = 1.20

!	Values at the MP interface (1-2)

! ENaC activity modified by local factors, including pH. Use hENaC_CNT as basal value of expression
!!! - FEMALE ADJUSTMENT ENaC x 1.20
  hENaC_CNT = FM_cENaC*35.0

  if (FAngDist) then
    hENaC_CNT = hENaC_CNT*2.90 ! cortical ENaC x 1.04 to 2.97 in grain-fed
  end if


! ROMK activity modified by local factors, including pH. Use hROMK_CNT as basal value of expression
  hROMK_CNT = 8.0d0

!  if (FAngDist) then
!    hROMK_CNT = hROMK_CNT*0.70 ! cortical ROMK x 0.52 to 0.74 in grain-fed
!  end if


  cnt(:)%h(3,1,2) = 0.80d0
  cnt(:)%h(4,1,2) = 0.16d0 
  cnt(:)%h(5,1,2) = PPCh2co3 
  cnt(:)%h(6,1,2) = PPCco2
  cnt(:)%h(7,1,2) = 1.0d-3/1.20
  cnt(:)%h(8,1,2) = 1.0d-3/1.20
  cnt(:)%h(9,1,2) = Purea
  cnt(:)%h(10,1,2) = PPCamon  
  cnt(:)%h(11,1,2) = cnt(:)%h(2,1,2)*0.20d0
  cnt(:)%h(12,1,2) = PPCprot 
  cnt(:)%h(13,1,2) = 0.0d0
  cnt(:)%h(14,1,2) = 0.0d0 
  cnt(:)%h(15,1,2) = 0.0d0
  cnt(:)%h(16,1,2) = 0.0d0

!	Values at the PE interface (2-5)
  cnt(:)%h(1,2,5) = 0.000d0
  cnt(:)%h(2,2,5) = 4.d0 
  cnt(:)%h(3,2,5) = 0.2d0 
  cnt(:)%h(4,2,5) = cnt(:)%h(3,2,5)*0.20d0
  cnt(:)%h(5,2,5) = PPCh2co3
  cnt(:)%h(6,2,5) = PPCco2
  cnt(:)%h(7,2,5) = 8.0d-3
  cnt(:)%h(8,2,5) = 8.0d-3
  cnt(:)%h(9,2,5) = Purea
  cnt(:)%h(10,2,5) = PPCamon 
  cnt(:)%h(11,2,5) = cnt(:)%h(2,2,5)*0.20d0
  cnt(:)%h(12,2,5) = PPCprot 
  cnt(:)%h(13,2,5) = 1.0d-4
  cnt(:)%h(14,2,5) = 1.0d-4
  cnt(:)%h(15,2,5) = 1.0d-4
  cnt(:)%h(16,2,5) = 0.0d0

!	Values at the MA interface (1-3)
  cnt(:)%h(1,1,3) = 0.0d0
  cnt(:)%h(2,1,3) = 0.00d0
  cnt(:)%h(3,1,3) = 0.0d0
  cnt(:)%h(4,1,3) = 0.0d0
  cnt(:)%h(5,1,3) = PICh2co3
  cnt(:)%h(6,1,3) = PICco2
  cnt(:)%h(7,1,3) = 0.0d0
  cnt(:)%h(8,1,3) = 0.0d0
  cnt(:)%h(9,1,3) = Purea
  cnt(:)%h(10,1,3) = PICamon
  cnt(:)%h(11,1,3) = 0.00d0 
  cnt(:)%h(12,1,3) = 0.0d0
  cnt(:)%h(13,1,3) = 0.0d0
  cnt(:)%h(14,1,3) = 0.0d0 
  cnt(:)%h(15,1,3) = 0.0d0
  cnt(:)%h(16,1,3) = 1.0d-5
  

!	Values at the AE interface (3-5)
  cnt(:)%hCLCA=PICAchlo
  cnt(:)%h(1,3,5) = 0.0d0
  cnt(:)%h(2,3,5) = 0.448d0
  cnt(:)%h(3,3,5) = PICAchlo
  cnt(:)%h(4,3,5) = 1.50d0
  cnt(:)%h(5,3,5) = PICh2co3
  cnt(:)%h(6,3,5) = PICco2
  cnt(:)%h(7,3,5) = PICAhpo4
  cnt(:)%h(8,3,5) = PICAhpo4
  cnt(:)%h(9,3,5) = Purea
  cnt(:)%h(10,3,5) = PICamon
  cnt(:)%h(11,3,5) = 0.18d0 
  cnt(:)%h(12,3,5) = PICAprot 
  cnt(:)%h(13,3,5) = 1.0d-4
  cnt(:)%h(14,3,5) = 1.0d-4 
  cnt(:)%h(15,3,5) = 1.0d-4
  cnt(:)%h(16,3,5) = 1.0d-5
  

!	Values at the MB interface (1-4)
  cnt(:)%h(1,1,4) = 0.0d0
  cnt(:)%h(2,1,4) = 0.0d0
  cnt(:)%h(3,1,4) = 0.0d0
  cnt(:)%h(4,1,4) = 0.0d0
  cnt(:)%h(5,1,4) = PICh2co3
  cnt(:)%h(6,1,4) = PICco2
  cnt(:)%h(7,1,4) = 0.0d0
  cnt(:)%h(8,1,4) = 0.0d0
  cnt(:)%h(9,1,4) = Purea
  cnt(:)%h(10,1,4) = PICamon
  cnt(:)%h(11,1,4) = 0.0d0 
  cnt(:)%h(12,1,4) = 0.0d0
  cnt(:)%h(13,1,4) = 0.0d0
  cnt(:)%h(14,1,4) = 0.0d0 
  cnt(:)%h(15,1,4) = 0.0d0
  cnt(:)%h(16,1,4) = 1.0d-5
  
!	Values at the BE interface (4-5)
  cnt(:)%hCLCB=PICBchlo
  cnt(:)%h(1,4,5) = 0.0d0
  cnt(:)%h(2,4,5) = 0.12d0
  cnt(:)%h(3,4,5) = PICBchlo
  cnt(:)%h(4,4,5) = 0.40d0 
  cnt(:)%h(5,4,5) = PICh2co3
  cnt(:)%h(6,4,5) = PICco2
  cnt(:)%h(7,4,5) = PICBhpo4
  cnt(:)%h(8,4,5) = PICBhpo4
  cnt(:)%h(9,4,5) = Purea
  cnt(:)%h(10,4,5) = PICamon
  cnt(:)%h(11,4,5) = 0.12d0 
  cnt(:)%h(12,4,5) = PICBprot 
  cnt(:)%h(13,4,5) = 1.0d-4
  cnt(:)%h(14,4,5) = 1.0d-4
  cnt(:)%h(15,4,5) = 1.0d-4 
  cnt(:)%h(16,4,5) = 1.0d-5

!	Values at basal interfaces (2-6), (3-6), (4-6)
  
  Do I = 1, NS
     cnt(:)%h(I,2,6) = cnt(:)%h(I,2,5)
     cnt(:)%h(I,3,6) = cnt(:)%h(I,3,5)
     cnt(:)%h(I,4,6) = cnt(:)%h(I,4,5)
  End do
  

!	Values at the ME interface (1-5)

  ClTJperm = 1000.0d0 ! Assume TJ resistivity of 5 mS/cm2
!!! - FEMALE ADJUSTMENT related to Claudin expression
  cnt(:)%h(1,1,5) = (1.00/FM_cENaC)*ClTJperm*1.0d0 ! Cldn8 and ENaC
  cnt(:)%h(2,1,5) = (1.00/FM_cENaC)*ClTJperm*1.2d0 ! Cldn8 and ENaC
! Cl permeability modified by local factors, including pH. Use hCltj_CNT as basal value of expression
  hCltj_CNT = 1.30*ClTJperm*1.2d0 ! Cldn7
  cnt(:)%h(4,1,5) = ClTJperm*0.6d0
  cnt(:)%h(5,1,5) = ClTJperm*2.0d0
  cnt(:)%h(6,1,5) = ClTJperm*2.0d0
  cnt(:)%h(7,1,5) = ClTJperm*0.14d0
  cnt(:)%h(8,1,5) = ClTJperm*0.14d0
  cnt(:)%h(9,1,5) = ClTJperm*0.40d0
  cnt(:)%h(10,1,5) = ClTJperm*6.0d0
  cnt(:)%h(11,1,5) = ClTJperm*1.2d0
  cnt(:)%h(12,1,5) = ClTJperm*10.d0
  cnt(:)%h(13,1,5) = 0.01d0*ClTJperm
  cnt(:)%h(14,1,5) = 0.01d0*ClTJperm
  cnt(:)%h(15,1,5) = 0.01d0*ClTJperm
  cnt(:)%h(16,1,5) = 0.0001d0*ClTJperm

!  if (FAngDist) then
!    cnt(:)%h(1,1,5) = cnt(:)%h(1,1,5)/2.0 ! Cldn8 and ENaC
!    cnt(:)%h(2,1,5) = cnt(:)%h(2,1,5)/2.0 ! Cldn8 and ENaC
! end if

!	Values at the ES interface (5-6)
  areafactor = 50.d0
  cnt(:)%h(1,5,6) = 240.d0*areafactor
  cnt(:)%h(2,5,6) = 320.d0*areafactor
  cnt(:)%h(3,5,6) = 320.d0*areafactor
  cnt(:)%h(4,5,6) = 160.d0*areafactor
  cnt(:)%h(5,5,6) = 240.d0*areafactor
  cnt(:)%h(6,5,6) = 240.d0*areafactor
  cnt(:)%h(7,5,6) = 160.d0*areafactor
  cnt(:)%h(8,5,6) = 160.d0*areafactor
  cnt(:)%h(9,5,6) = 160.d0*areafactor
  cnt(:)%h(10,5,6) = 320.d0*areafactor
  cnt(:)%h(11,5,6) = 320.0d0*areafactor
  cnt(:)%h(12,5,6) = 16000.0d0*areafactor
  cnt(:)%h(13,5,6) = 4.0d0*areafactor
  cnt(:)%h(14,5,6) = 4.0d0*areafactor
  cnt(:)%h(15,5,6) = 4.0d0*areafactor
  cnt(:)%h(16,5,6) = cnt(:)%h(1,5,6)*(7.93/13.3)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Non-Dimensional Solute permeabilities
!     Divide h(I,K,L) by (href) 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  Do I = 1,NS
     Do K = 1, NC
        Do L = K, NC
           cnt(:)%h(I,K,L) = cnt(:)%h(I,K,L)*1.0d-5/(href)
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
  dLA(1,12,2,5)=20.0d-9 
  dLA(1,12,2,6)=20.0d-9 
  dLA(1,12,3,5)=36.0d-9
  dLA(1,12,3,6)=36.0d-9
  dLA(1,12,4,5)=24.0d-9
  dLA(1,12,4,6)=24.0d-9
  
!	 Na2/HPO4 co-transporters on all basolateral membranes	 
  dLA(1,7,2,5)=2.0d-9
  dLA(1,7,2,6)=2.0d-9
  dLA(1,7,3,5)=1.2d-9
  dLA(1,7,3,6)=1.2d-9
  dLA(1,7,4,5)=0.80d-9
  dLA(1,7,4,6)=0.80d-9

!    Apical NCC in PC !!!!!!!! Removed
  xNCC = 0

!	 Basolateral Cl/HCO3 exchanger in PC 
  dLA(3,4,2,5)=2.0d-9
  dLA(3,4,2,6)=2.0d-9
  
!	 AE1 in IC-A only
  xAE1 = 150.d-9

!	 Transporters in IC-B only
!!! - FEMALE ADJUSTMENT Pendrin x 1.20
  dLA(3,4,1,4)=1.20*800.d-9
  xPendrin = dLA(3,4,1,4)/1400.d0

  if (FAngDist) then
    xPendrin = xPendrin*1.47 ! Pendrin x 1.47
  end if

!---------------------------------------------------------------------72
! Specific calcium transporters
!---------------------------------------------------------------------72
!  Basolateral NCX exchanger
  xNCX = 25.0d-9*1.60

!  Basolateral PMCA pump
  PMCA = 2.0d-9*0.40

!  Apical TRPV5
  xTRPV5 = 8.0d6

! Apical TRPV4 channel: density of channels times single channel permeability
! See Parikh et al., Biophys J 2015
  Po_TRPV4 = 1.0/(1.0 + 780.d0/37.d0)
  PCa_TRPV4 = 4.5d-8
  xPTRPV4 = Po_TRPV4*PCa_TRPV4

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
  ATPNaKPES = FM_cENaC*3410.d-9*1.25
  ATPNaKAES = 450.d-9
  ATPNaKBES = 300.d-9

!		H-ATPase (apical in IC-A, basolateral in IC-B)
  ATPHMA = 12000.0d-9
  ATPHBES = 400.d-9*15

!		H-K-ATPase (apical in PC and IC)
  ATPHKMP = 0.d0
  ATPHKMA = 720.d-9 
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
              cnt(:)%dLA(I,J,K,L) = dLA(I,J,K,L)/(href*Cref)
           End Do
        End Do
     End Do
  End Do
  
  cnt(:)%xPendrin = xPendrin/(href*Cref)
  cnt(:)%xAE1 = xAE1/(href*Cref)

  cnt(:)%xNCX = xNCX/(href*Cref)
  cnt(:)%PMCA= PMCA/(href*Cref)
  xTRPV5_cnt = xTRPV5/(href*Cref)
  xPTRPV4_cnt = xPTRPV4/(href*Cref)

  cnt(:)%ATPNaK(2,5) = ATPNaKPES/(href*Cref)
  cnt(:)%ATPNaK(2,6) = ATPNaKPES/(href*Cref)
  cnt(:)%ATPNaK(3,5) = ATPNaKAES/(href*Cref)
  cnt(:)%ATPNaK(3,6) = ATPNaKAES/(href*Cref)
  cnt(:)%ATPNaK(4,5) = ATPNaKBES/(href*Cref)
  cnt(:)%ATPNaK(4,6) = ATPNaKBES/(href*Cref)
  
  cnt(:)%ATPH(1,3) = ATPHMA/(href*Cref)
  cnt(:)%ATPH(4,5) = ATPHBES/(href*Cref)
  cnt(:)%ATPH(4,6) = ATPHBES/(href*Cref)
  
  cnt(:)%ATPHK(1,2) = ATPHKMP/(href*Cref)
  cnt(:)%ATPHK(1,3) = ATPHKMA/(href*Cref)
  cnt(:)%ATPHK(1,4) = ATPHKMB/(href*Cref)
  
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
  
  cnt(:)%dkd(1)=dkduncat*10.0d0
  cnt(:)%dkh(1)=dkhuncat*10.0d0
  cnt(:)%dkd(2)=dkdP
  cnt(:)%dkh(2)=dkhP
  cnt(:)%dkd(3)=dkdA
  cnt(:)%dkh(3)=dkhA
  cnt(:)%dkd(4)=dkdB
  cnt(:)%dkh(4)=dkhB
  cnt(:)%dkd(5)=dkdE
  cnt(:)%dkh(5)=dkhE
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	READ ENTERING AND PERITUBULAR CONDITIONS FROM DCT OUTLET FILE
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


  open ( unit=32, file='DCToutlet' )

  Do I = 1,4
     read(32,200),cnt(0)%conc(I,1),cnt(0)%conc(I,6)
  end Do
  Do I = 5,NS
     read(32,210),cnt(0)%conc(I,1),cnt(0)%conc(I,6)
  end Do
  read(32,200),cnt(0)%ph(1),cnt(0)%ph(6)
  read(32,200),cnt(0)%ep(1),cnt(0)%ep(6)
  read(32,210),cnt(0)%vol(1),cnt(0)%pres
  close ( unit=32 )
  cnt(:)%volLuminit = cnt(0)%vol(1)

  cnt(:)%ep(6)=cnt(0)%ep(6)
  
   pos(:) = 0.0  ! stays near cortex surface

  call set_intconc ( cnt, NZ, 1, pos )

  ! reference impermeant concentrations
  cnt(:)%cPimpref = 50.0d0
  cnt(:)%cAimpref = 18.0d0
  cnt(:)%cBimpref = 18.0d0

  ! total buffer concentrations
  cnt(:)%cPbuftot = 32.0d0
  cnt(:)%cAbuftot = 40.0d0
  cnt(:)%cBbuftot = 40.0d0
  
  !  initialize coalescence parameter
  do jz = 0, NZ
     cnt(jz)%coalesce =  2.00d0**(-2.32d0*(jz)/NZ)
  end do
  
!---------------------------------------------------------------------72
!     For metabolic calculations
!---------------------------------------------------------------------72

  if (ndiabetes .eq. 0)  then
    cnt(:)%TQ = 15.0 !TNa-QO2 ratio in normal CNT
  else
    cnt(:)%TQ = 12.0 !TNa-QO2 ratio in diabetic CNT
  end if


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

 200   format (6f12.5)
 210   format (6e12.5)
  
  return
end Subroutine initC
      
