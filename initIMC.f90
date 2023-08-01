!     This is the initialization routine.
!     It sets up several parameters, constants, and initial guess values.

Subroutine initIMC (imcd)

  include 'values.h'
  include 'global.h'
  include 'defs.h'

  ! passed variables
  type (membrane) :: imcd(0:NZ)
  
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
!---------------------------------------------------------------------72
!	 Calculate initial surfaces and volumes
!	 Volumes, surfaces, and angles implicitly include the number of
!	 each type of cell (but the PtoIratio doesn't).
!	 Data for rat (AMW model, AJP Renal 2001)
!---------------------------------------------------------------------72

  SlumPinitimc = 0.6d0
  SbasPinitimc = 0.6d0
  SlatPinitimc = 5.4d0
  
  SlumEinitimc = 0.001d0
  imcd(:)%sbasEinit = 0.020d0
  
  imcd(:)%volPinit = 8.0d0 
  imcd(:)%volEinit = 0.80d0 
  
  imcd(:)%volAinit = 8.0d0
  imcd(:)%volBinit = 8.0d0

  !	 Assign areas and volume to non-existent A and B cells
  SlumAinitimc = SlumPinitimc
  SbasAinitimc = SbasPinitimc
  SlatAinitimc = SlatPinitimc

  SlumBinitimc = SlumPinitimc
  SbasBinitimc = SbasPinitimc
  SlatBinitimc = SlatPinitimc

!---------------------------------------------------------------------72
!	 Assign membrane surface area except for basal membrane of E
!---------------------------------------------------------------------72

  imcd(:)%area(1,2)=SlumPinitimc
  imcd(:)%area(1,3)=SlumAinitimc
  imcd(:)%area(1,4)=SlumBinitimc
  imcd(:)%area(1,5)=SlumEinitimc
  
  imcd(:)%area(2,5)=SlatPinitimc
  imcd(:)%area(3,5)=SlatAinitimc
  imcd(:)%area(4,5)=SlatBinitimc
  
  imcd(:)%area(2,6)=SbasPinitimc
  imcd(:)%area(3,6)=SbasAinitimc
  imcd(:)%area(4,6)=SbasBinitimc
  
  Do K = 1,NC-1
     Do L = K+1,NC
        imcd(:)%area(L,K)=imcd(:)%area(K,L)
     End do
  End do

!---------------------------------------------------------------------72      
!---------------------------------------------------------------------72
!   Solute valence
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!				
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
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Water permeabilities
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  !!! - FEMALE ADJUSTMENT Pf x 1.0
  PfMP = 0.750d-3/Vwbar
  PfPE = 0.150d-3/Vwbar 
  PfPS = 0.150d-3/Vwbar 
  
  PfME = 36.0d-3/Vwbar
  PfES = 600.d-3/Vwbar 
   
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
  Pf(1,3) = 0.d0
  Pf(1,4) = 0.d0
  Pf(1,5) = PfME
  Pf(2,5) = PfPE
  Pf(3,5) = 0.d0
  Pf(4,5) = 0.d0
  Pf(2,6) = PfPS
  Pf(3,6) = 0.d0
  Pf(4,6) = 0.d0
  Pf(5,6) = PfES

!	  Units of dimensional water flux: cm3/s/cm2 epith
!	  Non-dimensional factor for water flux: (Pfref)*Vwbar*Cref
!	  Calculate non-dimensional dLPV = Lp*RT*Cref / (Pfref*Vwbar*Cref)
!	  dLPV = Pf/Pfref 

  Do K = 1, NC
     Do L = 1, NC
        imcd(:)%dLPV(K,L) = Pf(K,L)/(Pfref)
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
           imcd(:)%sig(I,K,L) = 1.0d0
        End Do
     End Do
  End Do

!	Values at the basement membrane (ES)
  Do I = 1,NS	
     imcd(:)%sig(I,5,6) = 0.0d0
  End Do
    
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Dimensional Solute permeabilities
!     h(I,K,L) = permeability of ith solute at the K-L interface
!	h(I,K,L) in 10-5 cm/s initially
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


!	Initialize 
  Do I = 1,NS
     Do K = 1, NC
        Do L = K, NC
           imcd(:)%h(I,K,L) = 0.0d0
        End Do
     End Do
  End Do

!	Values at the MP interface (1-2)

! ENaC activity modified by local factors, including pH. Use hENaC_IMC as basal value of expression
! imcd(:)%hNaMP = 1.0*fscaleENaC ! Not used, replaced with hENaC
!!! - FEMALE ADJUSTMENT ENaC
  hENaC_IMC = FM_mENaC*2.50

  if (FAngDist) then
    hENaC_IMC = hENaC_IMC*1.50 ! ENaC-gamma CL x 0.84 - 1.51 = x 1.175 in grain-fed
!     hENaC_IMC = hENaC_IMC*1.35 !1.70 !ENaC-gamma CL x 0.9 - 2.1 = x 1.7 in casein-fed
  end if

! ROMK activity modified by local factors, including pH. Use hROMK_IMC as basal value of expression
  hROMK_IMC = 1.0d0


  imcd(:)%h(2,1,2) = 1.0d0
  imcd(:)%h(3,1,2) = 2.0d-4
  imcd(:)%h(4,1,2) = 2.0d-4 
  imcd(:)%h(5,1,2) = 130.d0
  imcd(:)%h(6,1,2) = 15.0d3/6.5
  imcd(:)%h(7,1,2) = 0.20d-3
  imcd(:)%h(8,1,2) = 0.20d-3
  imcd(:)%h(9,1,2) = 300.d0
  imcd(:)%h(10,1,2) = 400.d0 
  imcd(:)%h(11,1,2) = 0.20d0
  imcd(:)%h(12,1,2) = 2.0d-4
  imcd(:)%h(13,1,2) = 0.d0
  imcd(:)%h(14,1,2) = 0.d0
  imcd(:)%h(15,1,2) = 0.d0
  imcd(:)%h(16,1,2) = 0.00010d0


!	Values at the PE interface (2-5)
  imcd(:)%h(1,2,5) = 2.0d-4
  imcd(:)%h(2,2,5) = 1.50d0 
  imcd(:)%h(3,2,5) = 0.050d0 
  imcd(:)%h(4,2,5) = 0.025d0
  imcd(:)%h(5,2,5) = 130.0d0
  imcd(:)%h(6,2,5) = 15.0d3/6.5 !See DN model (AJP 2008)
  imcd(:)%h(7,2,5) = 8.0d-3
  imcd(:)%h(8,2,5) = 8.0d-3
  imcd(:)%h(9,2,5) = 15.0d0
  imcd(:)%h(10,2,5) = 400.d0 
  imcd(:)%h(11,2,5) = 0.30d0
  imcd(:)%h(12,2,5) = 2.0d-4 
  imcd(:)%h(13,2,5) = 1.0d-4
  imcd(:)%h(14,2,5) = 1.0d-4 
  imcd(:)%h(15,2,5) = 1.0d-4 
  imcd(:)%h(16,2,5) = 0.00010d0

!	Values at other cellular interfaces (2-6)
!	Dummy values for non-existent intercalated cells (3-1),(3-5),(3-6),(4-1),(4-5),(4-6)

  Do I = 1, NS
     imcd(:)%h(I,2,6) = imcd(:)%h(I,2,5)
     imcd(:)%h(I,1,3) = 0.0d0
     imcd(:)%h(I,3,5) = 0.0d0
     imcd(:)%h(I,3,6) = 0.0d0
     imcd(:)%h(I,1,4) = 0.0d0
     imcd(:)%h(I,4,5) = 0.0d0
     imcd(:)%h(I,4,6) = 0.0d0
  End do
  
!	Values at the ME interface (1-5)
  ClTJperm = 1000.0d0
!!! - FEMALE ADJUSTMENT related to Claudin expression
  imcd(:)%h(1,1,5) = (1.00/FM_mENaC)*ClTJperm*0.90d0 ! Cldn8 and ENaC
  imcd(:)%h(2,1,5) = (1.00/FM_mENaC)*ClTJperm*0.70d0 ! Cldn8 and ENaC
! Cl permeability modified by local factors, including pH. Use hCltj_IMC as basal value of expression
  hCltj_IMC = 1.30*ClTJperm*1.60d0  ! Cldn7
  imcd(:)%h(4,1,5) = ClTJperm*0.40d0
  imcd(:)%h(5,1,5) = ClTJperm*1.20d0
  imcd(:)%h(6,1,5) = ClTJperm*1.20d0
  imcd(:)%h(7,1,5) = ClTJperm*0.20d0
  imcd(:)%h(8,1,5) = ClTJperm*0.20d0
  imcd(:)%h(9,1,5) = ClTJperm*0.80d0 
  imcd(:)%h(10,1,5) = ClTJperm*1.0d0
  imcd(:)%h(11,1,5) = ClTJperm*0.7d0
  imcd(:)%h(12,1,5) = ClTJperm*6.0d0
  imcd(:)%h(13,1,5) = ClTJperm*0.01d0
  imcd(:)%h(14,1,5) = ClTJperm*0.01d0
  imcd(:)%h(15,1,5) = ClTJperm*0.01d0
  imcd(:)%h(16,1,5) = ClTJperm*0.0140d0

  if (FAngAqp2) then
      imcd(:)%h(9,1,5) = imcd(:)%h(9,1,5)*0.50 ! AVP effects on urea permeability
  end if

  !	Values at the ES interface (5-6)
  areafactor = 50.d0
  imcd(:)%h(1,5,6) = 60.d0*areafactor
  imcd(:)%h(2,5,6) = 80.d0*areafactor
  imcd(:)%h(3,5,6) = 80.d0*areafactor
  imcd(:)%h(4,5,6) = 40.d0*areafactor
  imcd(:)%h(5,5,6) = 60.d0*areafactor
  imcd(:)%h(6,5,6) = 60.d0*areafactor
  imcd(:)%h(7,5,6) = 40.d0*areafactor
  imcd(:)%h(8,5,6) = 40.d0*areafactor
  imcd(:)%h(9,5,6) = 40.d0*areafactor
  imcd(:)%h(10,5,6) = 50.d0*areafactor
  imcd(:)%h(11,5,6) = 80.0d0*areafactor
  imcd(:)%h(12,5,6) = 400.0d0*areafactor
  imcd(:)%h(13,5,6) = 1.0*areafactor
  imcd(:)%h(14,5,6) = 1.0*areafactor
  imcd(:)%h(15,5,6) = 1.0*areafactor
  imcd(:)%h(16,5,6) = imcd(:)%h(1,5,6)*(7.93/13.3)


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Non-Dimensional Solute permeabilities
!     Divide h(I,K,L) by (href) 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  Do I = 1,NS
     Do K = 1, NC
        Do L = K, NC
           imcd(:)%h(I,K,L) = imcd(:)%h(I,K,L)*1.0d-5/(href)
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
  
  !	 Na/H exchangers on basolateral membrane	 
  dLA(1,12,2,5)=6.0d-9
  dLA(1,12,2,6)=6.0d-9
  
  !	 Na2/HPO4 co-transporters on basolateral membrane	 
  dLA(1,7,2,5)=2.0d-9
  dLA(1,7,2,6)=2.0d-9
  
  !	 NaKCl2 cotransporter on basolateral membrane
  dLA(1,2,2,5)=4.0d-9
  dLA(1,2,2,6)=4.0d-9
  
  !	 KCl cotransporter on basolateral membrane
  dLA(2,3,2,5)=250.0d-9
  dLA(2,3,2,6)=250.0d-9
  
  !	 Cl/HCO3 exchanger on basolateral membrane
  dLA(3,4,2,5)=20.0d-9
  dLA(3,4,2,6)=20.0d-9
  
  !	 Na/Cl exchanger on apical membrane
  dLA(1,3,1,2)=1200.0d-9
  
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

  !		Na-K-ATPase (basolateral)
  !!! - FEMALE ADJUSTMENT ENaC and NKA
  ATPNaKPES = FM_mENaC*2000.d-9*1.25

  !		H-K-ATPase (apical)
  ATPHKMP = 2000.d-9*0.10*0.50

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
              imcd(:)%dLA(I,J,K,L) = dLA(I,J,K,L)/(href*Cref)
           End Do
        End Do
     End Do
  End Do
  
  imcd(:)%ATPNaK(2,5) = ATPNaKPES/(href*Cref)
  imcd(:)%ATPNaK(2,6) = ATPNaKPES/(href*Cref)
  imcd(:)%ATPHK(1,2) = ATPHKMP/(href*Cref)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     Other kinetic parameters
!	Units of kinetic constants are s-1
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
! KINETIC PARAMETERS
       
  dkhuncat = 0.145d0
  dkduncat = 49.6d0
  dkhP = dkhuncat*100.d0
  dkdP = dkduncat*100.d0
  dkhE = dkhuncat
  dkdE = dkduncat
  
  imcd(:)%dkd(1)=dkduncat
  imcd(:)%dkh(1)=dkhuncat
  imcd(:)%dkd(2)=dkdP
  imcd(:)%dkh(2)=dkhP
  imcd(:)%dkd(3)=dkdP
  imcd(:)%dkh(3)=dkhP
  imcd(:)%dkd(4)=dkdP
  imcd(:)%dkh(4)=dkhP
  imcd(:)%dkd(5)=dkdE
  imcd(:)%dkh(5)=dkhE
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     BOUNDARY CONDITIONS IN PERITUBULAR SOLUTION
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  imcd(0)%ep(6) = -0.001d-3/EPref 

  Do jz=0,NZIMC
     imcd(jz)%ep(6)=imcd(0)%ep(6)

! coalescence parameter
     imcd(jz)%coalesce = 0.20*(2.00d0**(-6.00d0*jz/NZIMC))
     
     pos(jz) = 1.0d0*jz/NZIMC
  End do

  call set_intconc ( imcd, NZIMC, 3, pos )
     
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	READ ENTERING CONDITIONS FROM OMCD OUTLET FILE
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  open ( unit=62, file='OMCoutlet' )

  Do I = 1,4
     read(62,200),imcd(0)%conc(I,1),Cimc6
  end Do
  Do I = 5,NS
     read(62,210),imcd(0)%conc(I,1),Cimc6
  end Do
  read(62,200),imcd(0)%ph(1),phimc6
  read(62,200),imcd(0)%ep(1),EPimc6
  read(62,210),imcd(0)%vol(1),imcd(0)%pres
  close ( unit=62 )
  imcd(:)%volLuminit = imcd(0)%vol(1)
  
!---------------------------------------------------------------------72
!     For metabolic calculations
!---------------------------------------------------------------------72

  if (ndiabetes .eq. 0)  then
    imcd(:)%TQ = 15.0 !TNa-QO2 ratio in normal OMCD
  else
    imcd(:)%TQ = 12.0 !TNa-QO2 ratio in diabetic OMCD
  end if

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

 200   format (6f12.5)
 210   format (6e12.5)

  return
end Subroutine initIMC
      
