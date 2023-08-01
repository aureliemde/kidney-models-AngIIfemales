!     This is the initialization routine.
!     It sets up several parameters, constants, and initial guess values.

Subroutine initOMC (omcd)

  include 'values.h'
  include 'global.h'
  include 'defs.h'
  
  type (membrane) :: omcd(0:NZ)
  double precision theta(NC),Slum(NC),Slat(NC),Sbas(NC)
  double precision Pf(NC,NC),dLA(NS,NS,NC,NC)
  double precision pos(0:NZ)
  

!---------------------------------------------------------------------72
!	 Data for rat (AMW model, AJP Renal 2001)
!---------------------------------------------------------------------72

  SlumPinitomc = 1.2d0
  SbasPinitomc = 1.2d0
  SlatPinitomc = 6.9d0

  SlumAinitomc = 2.0d0 
  SbasAinitomc = 2.0d0 
  SlatAinitomc = 6.0d0 
  
  SlumEinitomc = 0.001d0
  omcd(:)%sbasEinit = 0.020d0
  
  omcd(:)%volPinit = 4.0d0 
  omcd(:)%volAinit = 3.0d0
  omcd(:)%volEinit = 0.80d0

  !	 Assign areas and volume to non-existent B cells
  SlumBinitomc = SlumAinitomc
  SbasBinitomc = SbasAinitomc
  SlatBinitomc = SlatAinitomc
  omcd(:)%volBinit = omcd(0)%volAinit

!---------------------------------------------------------------------72
!	 Assign membrane surface area except for basal membrane of E
!---------------------------------------------------------------------72

  omcd(:)%area(1,2)=SlumPinitomc
  omcd(:)%area(1,3)=SlumAinitomc
  omcd(:)%area(1,4)=SlumBinitomc
  omcd(:)%area(1,5)=SlumEinitomc
	 
  omcd(:)%area(2,5)=SlatPinitomc
  omcd(:)%area(3,5)=SlatAinitomc
  omcd(:)%area(4,5)=SlatBinitomc
  
  omcd(:)%area(2,6)=SbasPinitomc
  omcd(:)%area(3,6)=SbasAinitomc
  omcd(:)%area(4,6)=SbasBinitomc
  
  Do K = 1,NC-1
     Do L = K+1,NC
        omcd(:)%area(L,K)=omcd(:)%area(K,L)
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
  omcd(:)%zPimp = -1.0d0
  omcd(:)%zAimp = -1.0d0
  omcd(:)%zBimp = -1.0d0

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Water permeabilities
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!!! - FEMALE ADJUSTMENT Pf x 1.0
  PfMP = 0.20d0/1.2d0
  PfPE = 0.11d0*2
  PfPS = 0.11d0*2
  
  PfMA = 0.22d-3
  PfAE = 5.50d-3
  PfAS = 5.50d-3
  
  PfMB = 0.22d-3
  PfBE = 5.50d-3
  PfBS = 5.50d-3
  
  PfME = 28.d0
  PfES = 41.d0
  
  if (FAngAqp2) then
     PfMP = PfMP*0.59 ! AQP2-med x 0.53-0.66
     PfPE = PfPE*0.59 ! AQP2-med x 0.53-0.66
     PfPS = PfPS*0.59 ! AQP2-med x 0.53-0.66
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
        omcd(:)%dLPV(K,L) = Pf(K,L)/(Pfref)
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
           omcd(:)%sig(I,K,L) = 1.0d0
        End Do
     End Do
  End Do


!	Values at the basement membrane (ES)
  Do I = 1,NS	
     omcd(:)%sig(I,5,6) = 0.0d0
  End Do
  

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Dimensional Solute permeabilities
!     h(I,K,L) = permeability of ith solute at the K-L interface
!	h(I,K,L) in 10-5 cm/s initially
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  PICAchlo = 1.2d0
  PPCh2co3 = 130.d0
  PICh2co3 = 10.d0 
  PPCco2 = 2.0d3
  PICco2 = 900.d0 
  PPChpo4 = 8.0d-3
  PPCprot = 2000.d0
  PICAprot = 1.5d0
  PPCamon = 2000.d0
  PICamon = 900.d0 
  PPCurea = 0.10d0
  PICurea = 1.0d0

!	Initialize 
  Do I = 1,NS
     Do K = 1, NC
        Do L = K, NC
           omcd(:)%h(I,K,L) = 0.0d0
        End Do
     End Do
  End Do
  
  fac4 = 1.0d0/4.0d0

!	Values at the MP interface (1-2)

! Female-to-male ENaC expression ratio in medulla (OMCD and IMCD)
FM_mENaC = 1.20

! ENaC activity modified by local factors, including pH. Use hENaC_OMC as basal value of expression
!!! - FEMALE ADJUSTMENT ENaC
  hENaC_OMC = FM_mENaC*14.0*2.50*fac4

  if (FAngDist) then
    hENaC_OMC = hENaC_OMC*1.50 !ENaC-gamma CL x 0.84 - 1.51 = x 1.175 in grain-fed
!    hENaC_OMC = hENaC_OMC*1.35! 1.70 !ENaC-gamma CL x 0.9 - 2.1 = x 1.7 in casein-fed
  end if

! ROMK activity modified by local factors, including pH. Use hROMK_OMC as basal value of expression
  hROMK_OMC = 8.0d0*fac4

  omcd(:)%h(3,1,2) = 0.80d0*fac4
  omcd(:)%h(4,1,2) = 0.16d0*fac4
  omcd(:)%h(5,1,2) = PPCh2co3 
  omcd(:)%h(6,1,2) = PPCco2
  omcd(:)%h(7,1,2) = 1.0d-3/1.20
  omcd(:)%h(8,1,2) = 1.0d-3/1.20
  omcd(:)%h(9,1,2) = PPCurea
  omcd(:)%h(10,1,2) = PPCamon  
  omcd(:)%h(11,1,2) = omcd(:)%h(2,1,2)*0.20d0
  omcd(:)%h(12,1,2) = PPCprot*fac4 
  omcd(:)%h(13,1,2) = 0.0d0
  omcd(:)%h(14,1,2) = 0.0d0
  omcd(:)%h(15,1,2) = 1.0d-10
  omcd(:)%h(16,1,2) = 0.00010d0

!	Values at the PE interface (2-5)
  omcd(:)%h(1,2,5) = 0.000d0
  omcd(:)%h(2,2,5) = 4.0d0*fac4
  omcd(:)%h(3,2,5) = 0.20d0*fac4
  omcd(:)%h(4,2,5) = omcd(:)%h(3,2,5)*0.20d0
  omcd(:)%h(5,2,5) = PPCh2co3
  omcd(:)%h(6,2,5) = PPCco2
  omcd(:)%h(7,2,5) = 8.0d-3
  omcd(:)%h(8,2,5) = 8.0d-3
  omcd(:)%h(9,2,5) = PPCurea
  omcd(:)%h(10,2,5) = PPCamon 
  omcd(:)%h(11,2,5) = omcd(:)%h(2,2,5)*0.20d0
  omcd(:)%h(12,2,5) = PPCprot*fac4 
  omcd(:)%h(13,2,5) = 1.0d-4
  omcd(:)%h(14,2,5) = 1.0d-4
  omcd(:)%h(14,2,5) = 1.0d-4
  omcd(:)%h(16,2,5) = 0.00010d0

!	Values at the MA interface (1-3)
  omcd(:)%h(1,1,3) = 0.0d0
  omcd(:)%h(2,1,3) = 0.00d0  
  omcd(:)%h(3,1,3) = 0.0d0
  omcd(:)%h(4,1,3) = 0.0d0
  omcd(:)%h(5,1,3) = PICh2co3
  omcd(:)%h(6,1,3) = PICco2
  omcd(:)%h(7,1,3) = 0.0d0
  omcd(:)%h(8,1,3) = 0.0d0
  omcd(:)%h(9,1,3) = PICurea
  omcd(:)%h(10,1,3) = PICamon
  omcd(:)%h(11,1,3) = 0.10d-4 
  omcd(:)%h(12,1,3) = 0.0d0
  omcd(:)%h(13,1,3) = 0.0d0
  omcd(:)%h(14,1,3) = 0.0d0
  omcd(:)%h(15,1,3) = 0.0d0
  omcd(:)%h(16,1,3) = 0.00010d0

!	Values at the AE interface (3-5)
  omcd(:)%hCLCA=PICAchlo
  omcd(:)%h(1,3,5) = 0.0d0
  omcd(:)%h(2,3,5) = 0.12d0
  omcd(:)%h(3,3,5) = PICAchlo
  omcd(:)%h(4,3,5) = 0.15d0
  omcd(:)%h(5,3,5) = PICh2co3
  omcd(:)%h(6,3,5) = PICco2
  omcd(:)%h(7,3,5) = 1.20d-3
  omcd(:)%h(8,3,5) = 1.20d-3
  omcd(:)%h(9,3,5) = PICurea
  omcd(:)%h(10,3,5) = PICamon
  omcd(:)%h(11,3,5) = 0.030d0
  omcd(:)%h(12,3,5) = 1.50d0
  omcd(:)%h(13,3,5) = 1.0d-4
  omcd(:)%h(14,3,5) = 1.0d-4
  omcd(:)%h(15,3,5) = 1.0d-4
  omcd(:)%h(16,3,5) = 0.00010d0

!	Values at other cellular interfaces (2-6), (3-6), (4-1), (4-5) and (4-6)
  
  Do I = 1, NS
     omcd(:)%h(I,2,6) = omcd(:)%h(I,2,5)
     omcd(:)%h(I,3,6) = omcd(:)%h(I,3,5)
     omcd(:)%h(I,1,4) = 0.0d0
     omcd(:)%h(I,4,5) = 0.0d0
     omcd(:)%h(I,4,6) = 0.0d0
  End do
  

!	Values at the ME interface (1-5)

  ClTJperm = 1000.0d0 ! Assume TJ resistivity of 5 mS/cm2
!!! - FEMALE ADJUSTMENT related to Claudin expression
  omcd(:)%h(1,1,5) = (1.00/FM_mENaC)*ClTJperm*0.80d0 ! Cldn8 and ENaC
  omcd(:)%h(2,1,5) = (1.00/FM_mENaC)*ClTJperm*1.20d0 ! Cldn8 and ENaC
! Cl permeability modified by local factors, including pH. Use hCltj_OMC as basal value of expression
  hCltj_OMC = 1.30*ClTJperm*1.00d0  ! Cldn7
  omcd(:)%h(4,1,5) = ClTJperm*0.30d0
  omcd(:)%h(5,1,5) = ClTJperm*1.20d0
  omcd(:)%h(6,1,5) = ClTJperm*1.20d0
  omcd(:)%h(7,1,5) = ClTJperm*0.10d0
  omcd(:)%h(8,1,5) = ClTJperm*0.10d0
  omcd(:)%h(9,1,5) = ClTJperm*2.0d0 !Per AMW email on 03/08/13
  omcd(:)%h(10,1,5) = ClTJperm*1.0d0
  omcd(:)%h(11,1,5) = ClTJperm*1.5d0
  omcd(:)%h(12,1,5) = ClTJperm*6.0d0
  omcd(:)%h(13,1,5) = ClTJperm*0.01d0
  omcd(:)%h(14,1,5) = ClTJperm*0.01d0
  omcd(:)%h(15,1,5) = ClTJperm*0.01d0
  omcd(:)%h(16,1,5) = ClTJperm*0.0140d0 ! See my 2015 AJP and Carney 1988

!	Values at the ES interface (5-6)
  areafactor = 50.d0
  omcd(:)%h(1,5,6) = 89.d0*areafactor
  omcd(:)%h(2,5,6) = 118.d0*areafactor
  omcd(:)%h(3,5,6) = 118.d0*areafactor
  omcd(:)%h(4,5,6) = 59.d0*areafactor
  omcd(:)%h(5,5,6) = 89.d0*areafactor
  omcd(:)%h(6,5,6) = 89.d0*areafactor
  omcd(:)%h(7,5,6) = 59.d0*areafactor
  omcd(:)%h(8,5,6) = 59.d0*areafactor
  omcd(:)%h(9,5,6) = 59.d0*areafactor
  omcd(:)%h(10,5,6) = 59.d0*areafactor
  omcd(:)%h(11,5,6) = 118.0d0*areafactor
  omcd(:)%h(12,5,6) = 590.0d0*areafactor
  omcd(:)%h(13,5,6) = 1.0d0*areafactor
  omcd(:)%h(14,5,6) = 1.0d0*areafactor
  omcd(:)%h(15,5,6) = 1.0d0*areafactor
  omcd(:)%h(16,5,6) = omcd(:)%h(1,5,6)*(7.93/13.3)


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Non-Dimensional Solute permeabilities
!     Divide h(I,K,L) by (href) 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  Do I = 1,NS
     Do K = 1, NC
        Do L = K, NC
           omcd(:)%h(I,K,L) = omcd(:)%h(I,K,L)*1.0d-5/(href)
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
  dLA(1,12,2,5)=2.50d-9
  dLA(1,12,2,6)=2.50d-9
  dLA(1,12,3,5)=6.0d-9
  dLA(1,12,3,6)=6.0d-9
  
  xNHE1P = 6.95d-10*680/695*fac4

!	 Na2/HPO4 co-transporters on all basolateral membranes	 
  dLA(1,7,2,5)=2.0d-9
  dLA(1,7,2,6)=2.0d-9
  dLA(1,7,3,5)=0.20d-9
  dLA(1,7,3,6)=0.20d-9
  
!	 Basolateral Cl/HCO3 exchanger in PC 
  dLA(3,4,2,5)=2.0d-9*fac4
  dLA(3,4,2,6)=2.0d-9*fac4

!	 AE1 in IC-A only
  xAE1 = 12.0d-9

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
!	Dimensional ATPa!e coefficients
!     ATPcoeff in mmol/s or mmol
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


!		Na-K-ATPase (basolateral in PC and IC)
  !!! - FEMALE ADJUSTMENT ENaC and NKA
  ATPNaKPES = FM_mENaC*3410.d-9*1.25*fac4
  ATPNaKAES = 75.d-9

!		H-ATPase (apical in IC-A, basolateral in IC-B)
  ATPHMA = 750.0d-9

!		H-K-ATPase (apical in PC and IC)
  ATPHKMP = 0.d0
  ATPHKMA = 150.d-9
  

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
              omcd(:)%dLA(I,J,K,L) = dLA(I,J,K,L)/(href*Cref)
           End Do
        End Do
     End Do
  End Do
  
  omcd(:)%xNHE1(2) = xNHE1P/(href*Cref)  
  omcd(:)%xAE1 = xAE1/(href*Cref)
		
  omcd(:)%ATPNaK(2,5) = ATPNaKPES/(href*Cref)
  omcd(:)%ATPNaK(2,6) = ATPNaKPES/(href*Cref)
  omcd(:)%ATPNaK(3,5) = ATPNaKAES/(href*Cref)
  omcd(:)%ATPNaK(3,6) = ATPNaKAES/(href*Cref)
  
  omcd(:)%ATPH(1,3) = ATPHMA/(href*Cref)
  
  omcd(:)%ATPHK(1,2) = ATPHKMP/(href*Cref)
  omcd(:)%ATPHK(1,3) = ATPHKMA/(href*Cref)
  
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
  dkhE = 0.145d0
  dkdE = 49.6d0

  omcd(:)%dkd(1)=dkduncat 
  omcd(:)%dkh(1)=dkhuncat
  omcd(:)%dkd(2)=dkdP
  omcd(:)%dkh(2)=dkhP
  omcd(:)%dkd(3)=dkdA
  omcd(:)%dkh(3)=dkhA
  omcd(:)%dkd(4)=dkdB
  omcd(:)%dkh(4)=dkhB
  omcd(:)%dkd(5)=dkdE
  omcd(:)%dkh(5)=dkhE
  

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     BOUNDARY CONDITIONS	IN PERITUBULAR SOLUTION
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  omcd(0)%ep(6) = -0.001d-3/EPref 
  
  Do jz=0,NZ
     
     omcd(jz)%ep(6)=omcd(0)%ep(6)

     pos(jz) = 1.0d0*jz/NZ

  End do

  call set_intconc (omcd, NZ, 2, pos)
  
  !  coalesence parameter
  omcd(:)%coalesce = 0.2d0
  

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	READ ENTERING CONDITIONS FROM CCD OUTLET FILE
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  open ( unit=52, file='CCDoutlet' )
  
  Do I = 1,4
     read(52,200),omcd(0)%conc(I,1),Comc6
  end Do
  Do I = 5,NS
     read(52,210),omcd(0)%conc(I,1),Comc6
  end Do
  read(52,200),omcd(0)%ph(1),phomc6
  read(52,200),omcd(0)%ep(1),EPomc6
  read(52,210),omcd(0)%vol(1),omcd(0)%pres
  close ( unit=52 )

  omcd(:)%volLuminit = omcd(0)%vol(1)
  
!---------------------------------------------------------------------72
!     For metabolic calculations
!---------------------------------------------------------------------72

  if (ndiabetes .eq. 0)  then
    omcd(:)%TQ = 15.0 !TNa-QO2 ratio in normal OMCD
  else
    omcd(:)%TQ = 12.0 !TNa-QO2 ratio in diabetic OMCD
  end if

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

200 format (6f12.5)
210 format (6e12.5)

  return
end Subroutine initOMC
      
