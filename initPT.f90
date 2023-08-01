!     This is the initialization routine for the PT
!     It sets up several parameters, constants, and initial guess values.

Subroutine initPT (pt)

  include 'values.h'
  include 'global.h'
  include 'defs.h'

  ! passed variables
  type(membrane) :: pt(0:NZ)
  
!  local variables
  double precision Pf(NC,NC)
  logical :: bdiabetes  ! true if simulating a diabetic kidney
  double precision pos(0:NZ)
  integer ind

  bdiabetes = .false.
  if (.not. bdiabetes)  then
     ndiabetes = 0
  else
     ndiabetes = 1
  end if


!---------------------------------------------------------------------72      
!---------------------------------------------------------------------72
!   Solute Indices:  1 = Na+, 2 = K+, 3 = Cl-, 4 = HCO3-, 5 = H2CO3, 6 = CO2
!	7 = HPO4(2-), 8 = H2PO4-, 9 = urea, 10 = NH3, 11 = NH4+, 12 = H+
!   13 = HCO2-, 14 = H2CO2, 15 = glucose
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72      
!---------------------------------------------------------------------72
!	 Calculate initial surfaces and volumes
!	 Volumes, surfaces, and angles implicitly include the number of
!	 each type of cell (but the PtoIratio doesn't).
!---------------------------------------------------------------------72      
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	 Data for rat (AMW model, AJP Renal 2007)
!---------------------------------------------------------------------72

  SlumPinitpt = 36.0d0
  SbasPinitpt = 1.0d0
  SlatPinitpt = 36.0d0

  SlumEinitpt = 0.001d0
  pt(:)%sbasEinit = 0.020d0

  pt(:)%volPinit = 10.0d0
  pt(:)%volEinit = 0.7d0

  pt(:)%volAinit = 10.0d0 ! Needed for water flux subroutine (AURELIE)
  pt(:)%volBinit = 10.0d0 ! Needed for water flux subroutine (AURELIE)

!---------------------------------------------------------------------72
!	 Assign membrane surface area except for basal membrane of E
!---------------------------------------------------------------------72

  if (.not. bdiabetes) then
    pt(:)%diam = 0.00250d0   ! normal kidney
  else
    pt(:)%diam = 0.0030d0    ! diabetic kidney, hypertrophied PT
  end if

  pt(:)%area(1,2)=SlumPinitpt 
  pt(:)%area(1,3)=SlumPinitpt
  pt(:)%area(1,4)=SlumPinitpt
  pt(:)%area(1,5)=SlumEinitpt
  
  pt(:)%area(2,5)=SlatPinitpt
  pt(:)%area(3,5)=SlatPinitpt
  pt(:)%area(4,5)=SlatPinitpt
  
  pt(:)%area(2,6)=SbasPinitpt
  pt(:)%area(3,6)=SbasPinitpt
  pt(:)%area(4,6)=SbasPinitpt

  do J = 0, NZ
     if (J >= xS3*NZ) then
        pt(J)%area(1,2) = 0.5*pt(J)%area(1,2) !*1.50/1.75
        pt(J)%area(2,5) = 0.5*pt(J)%area(2,5) !*1.50/1.75
        pt(J)%area(2,6) = 0.5*pt(J)%area(2,6) !*1.50/1.75
     end if
  end do
 
  Do K = 1,NC-1
     Do L = K+1,NC
        pt(:)%area(L,K)=pt(:)%area(K,L)
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

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Water permeabilities
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!   The permeabilities have already been multiplied by the area (cm2/cm2 epith)
!   dLPV(K,L) = hydraulic permeability at the K-L interface, in cm/s/mmHg
!   Pf(K,L) = osmotic permeability at the K-L interface, in cm/s
!   Relationship between the two: LPV=Pf*(Vwbar/RT)

! WATER PERMEABILITIES
!!! FEMALE ADJUSTMENT - Pf x 0.64
  PfMP = 0.64*0.40/36.0
  PfMA = 0.0d0
  PfMB = 0.0d0
  PfME = 0.22/0.0010
  PfPE = 0.64*0.40/36.0
  PfAE = 0.0d0
  PfBE = 0.0d0
  PfPS = PfPE
  PfAS = 0.0d0
  PfBS = 0.0d0
  PfES = 6.60/0.020
  

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
        pt(:)%dLPV(K,L) = Pf(K,L)/(Pfref)
     End Do
  End Do
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     sig(I,K,L) = reflection coefficient of Ith solute between K and L
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	Initialize 
  Do I = 1,NSPT
     Do K = 1, NC
        Do L = 1, NC
           pt(:)%sig(I,K,L) = 1.0d0
        End Do
     End Do
  End Do

!    Reflection coefficient = 0 at the basement membrane (ES)
  Do I = 1,NSPT
     pt(:)%sig(I,5,6) = 0.0d0
  End Do
  
!	Reflection coefficient varies at the tight junction (ME)
  pt(:)%sig(1,1,5) = 0.750
  pt(:)%sig(2,1,5) = 0.600
  pt(:)%sig(3,1,5) = 0.300
  pt(:)%sig(4,1,5) = 0.900
  pt(:)%sig(5,1,5) = 0.900
  pt(:)%sig(6,1,5) = 0.900
  pt(:)%sig(7,1,5) = 0.900
  pt(:)%sig(8,1,5) = 0.900
  pt(:)%sig(9,1,5) = 0.700
  pt(:)%sig(10,1,5) = 0.300
  pt(:)%sig(11,1,5) = 0.600
  pt(:)%sig(12,1,5) = 0.200
  pt(:)%sig(13,1,5) = 0.300
  pt(:)%sig(14,1,5) = 0.700
  pt(:)%sig(15,1,5) = 1.000
  pt(:)%sig(16,1,5) = 0.890

  Do I = 1,NSPT	
     pt(:)%sig(I,5,1) = pt(:)%sig(I,1,5)
     pt(:)%sig(I,6,5) = pt(:)%sig(I,5,6)
  End Do
  

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Dimensional Solute permeabilities
!	Alan's permeability values are multiplied by the area (cm2/cm2 epith)
!   h(I,K,L) = permeability of ith solute at the K-L interface
!	h(I,K,L) in 10-5 cm/s initially
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	Initialize 
  Do I = 1,NSPT
     Do K = 1, NC
        Do L = K, NC
           pt(:)%h(I,K,L) = 0.0d0
        End Do
     End Do
  End Do
  
!	Values at the lumen-cell interface (MP, 1-2)
  pt(:)%h(1,1,2) = 0.0
  pt(:)%h(2,1,2) = 0.90/36.0
  pt(:)%h(3,1,2) = 0.0
  pt(:)%h(4,1,2) = 0.036/36.0
  pt(:)%h(5,1,2) = 2.34d3/36.0
  pt(:)%h(6,1,2) = 4.32d4/36.0
  pt(:)%h(7,1,2) =  0.0
  pt(:)%h(8,1,2) =  0.000/36.0
  pt(:)%h(9,1,2) =  3.78/36.0
  pt(:)%h(10,1,2) = 3.06d3/36.0
  pt(:)%h(11,1,2) = 0.774/36.0
  pt(:)%h(12,1,2) = 3.06d4/36.0
  pt(:)%h(13,1,2) = 0.0
  pt(:)%h(14,1,2) = 1.80d5/36.0
  pt(:)%h(15,1,2) = 0.0
  pt(:)%h(16,1,2) = 0.005d0

!	Values at the cell-LIS interface (PE, 2-5)
  pt(:)%h(1,2,5) = 0.0
  pt(:)%h(2,2,5) = 0.20
  pt(:)%h(3,2,5) = 0.01
  pt(:)%h(4,2,5) = 0.0
  pt(:)%h(5,2,5) = 65.0
  pt(:)%h(6,2,5) = 1.20d3
  pt(:)%h(7,2,5) = 0.00225
  pt(:)%h(8,2,5) = 0.0330
  pt(:)%h(9,2,5) = 0.100
  pt(:)%h(10,2,5) = 100.0 
  pt(:)%h(11,2,5) = 0.060
  pt(:)%h(12,2,5) = 850.0 
  pt(:)%h(13,2,5) = 0.0190
  pt(:)%h(14,2,5) = 6.0d3
  pt(:)%h(15,2,5) = 0.0
  pt(:)%h(16,2,5) = 0.0d0

!	Values at the cell-interstitium interface (PS, 2-6)
  Do I = 1, NSPT
     pt(:)%h(I,2,6) = pt(:)%h(I,2,5)
  End do

!	Values at the MA interface (1-3)
!	Values at the AE interface (3-5)
!	Values at the AS interface (3-6)
!	Values at the MB interface (1-4)
!	Values at the BE interface (4-5)
!	Values at the BS interface (4-6)

  Do I = 1, NSPT
     pt(:)%h(I,1,3) = 0.d0
     pt(:)%h(I,3,5) = 0.d0
     pt(:)%h(I,3,6) = 0.d0
     pt(:)%h(I,1,4) = 0.d0
     pt(:)%h(I,4,5) = 0.d0
     pt(:)%h(I,4,6) = 0.d0
  End do


!	Values at the lumen-LIS interface (ME, 1-5)
  ClTJperm = 1.d0/0.0010d0 ! Factor for area = 1/0.001
!!! FEMALE ADJUSTMENT - PNa and PK (not PCl) x 0.38
  pt(:)%h(1,1,5) = 0.38*26.0*ClTJperm
  pt(:)%h(2,1,5) = 0.38*29.0*ClTJperm
  pt(:)%h(3,1,5) = 20.0*ClTJperm
  pt(:)%h(4,1,5) =  8.0*ClTJperm
  pt(:)%h(5,1,5) =  8.0*ClTJperm
  pt(:)%h(6,1,5) =  8.0*ClTJperm
  pt(:)%h(7,1,5) =  4.0*ClTJperm
  pt(:)%h(8,1,5) =  4.0*ClTJperm
  pt(:)%h(9,1,5) =  8.0*ClTJperm
  pt(:)%h(10,1,5) = 50.0*ClTJperm
  pt(:)%h(11,1,5) = 50.0*ClTJperm
  pt(:)%h(12,1,5) = 600.0*ClTJperm
  pt(:)%h(13,1,5) = 14.0*ClTJperm
  pt(:)%h(14,1,5) = 28.0*ClTJperm
  pt(:)%h(15,1,5) = 0.31*ClTJperm
  pt(:)%h(16,1,5) = 20.0d0*ClTJperm


!	Values at the LIS-interstitium interface (ES, 5-6)
  areafactor = 1.0d0/0.020d0
  pt(:)%h(1,5,6) = 100.0*areafactor
  pt(:)%h(2,5,6) = 140.0*areafactor
  pt(:)%h(3,5,6) = 120.0*areafactor
  pt(:)%h(4,5,6) = 100.0*areafactor
  pt(:)%h(5,5,6) = 100.0*areafactor
  pt(:)%h(6,5,6) = 100.0*areafactor
  pt(:)%h(7,5,6) =  80.0*areafactor
  pt(:)%h(8,5,6) =  80.0*areafactor
  pt(:)%h(9,5,6) =  160.0*areafactor
  pt(:)%h(10,5,6) = 400.0*areafactor
  pt(:)%h(11,5,6) = 400.0*areafactor
  pt(:)%h(12,5,6) = 6.0d4*areafactor
  pt(:)%h(13,5,6) = 100.0*areafactor
  pt(:)%h(14,5,6) = 180.0*areafactor
  pt(:)%h(15,5,6) =  60.0*areafactor
  pt(:)%h(16,5,6) = pt(:)%h(1,5,6)*(7.93/13.3)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Non-Dimensional Solute permeabilities
!   Divide h(I,K,L) by (href) 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  Do I = 1,NSPT
     Do K = 1, NC
        Do L = K, NC
           pt(:)%h(I,K,L) = pt(:)%h(I,K,L)*1.0d-5/(href)
        End Do
     End Do
  End Do
	
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Dimensional NET coefficients
!   dLA(I,J,K,L) = coefficient of ith and jth solute at the K-L interface
!	dLA(I,J,K,L) in mmol2/J/s initially
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Initialize 
  Do I = 1,NSPT
     Do J = 1, NSPT
        Do K = 1, NC
           Do L = K, NC
              pt(:)%dLA(I,J,K,L) = 0.0d0
           End Do
        End Do
     End Do
  End Do

!   Apical NHE3
!!! FEMALE ADJUSTMENT - NHE3 x 0.83 or 0.78
  xNHE3 = 0.780*1000.0d-9/36.0

  if (FAngProx) then
     xNHE3 = xNHE3*0.70! 0.90d0 ! NH3p x 1.09 in grain-fed
  end if

!	SGLT1 and SGLT2 on apical membrane
!!! FEMALE ADJUSTMENT - SGLT2 x 1.16
  CTsglt1 = 0.30d-5/36.0 ! SGLT1
  pt(:)%dLA(1,15,1,2) = 1.16*1.300d-9 ! SGLT2 (original AMW value=270/36=7.50)

  if (FAngProx) then
     CTsglt1 = Ctsglt1*2.05 ! SLGT1-med x 2.05
     pt(:)%dLA(1,15,1,2) =  pt(:)%dLA(1,15,1,2)*0.65 ! SGLT2 x 0.65
  end if

!	 glucose transporter on basolateral membrane 
  CTglut1 = 1.2*0.2500d-5
  CTglut2 = 1.3*0.1250d-5

!	 Na-H2PO4 co-transporter on apical membrane
!!! FEMALE ADJUSTMENT - NaPi x 0.75
  pt(:)%dLA(1,8,1,2)=1.25d-9
  CTnapiIIa = 0.75*0.30d-9
  CTnapiIIb = 0.75*0.0d0
  CTnapiIIc = 0.75*0.250d-9
  CTpit2 = 0.10d-9

!	 Cl/HCO3 exchanger on apical membrane
  pt(:)%dLA(3,4,1,2)=0.0
  
!	 Cl/HCO2 exchanger on apical membrane
  pt(:)%dLA(3,13,1,2)=5.0d-9*2.0

!	 K-Cl co-transporter on basolateral membrane	 
  pt(:)%dLA(2,3,2,5)=5.0d-9

!  if (FAngProx) then
!     pt(:)%dLA(2,3,2,5) = pt(:)%dLA(2,3,2,5)*1.28 ! KCC3/4 x 1.28-1.22 in grain-fed
!  end if

  pt(:)%dLA(2,3,2,6)=pt(:)%dLA(1,12,2,5)

!	 Na/3HCO3 exchanger on basolateral membrane
  pt(:)%dLA(1,4,2,5)=5.0d-9
  pt(:)%dLA(1,4,2,6)=pt(:)%dLA(1,4,2,5)

!   NDCBE (Na-2HCO3/Cl cotransporter) on basolateral membrane
  pt(:)%dLA(1,3,2,5)=35.0d-9
  pt(:)%dLA(1,3,2,6)=pt(:)%dLA(1,3,2,5)
! Note that 3 solutes are transported, even though there are only 2 solute indices

!    Basolateral Na,K-APTase
!!! FEMALE ADJUSTMENT - NHE3 and NKA x 0.78
  ATPNaKPES = 0.780*300.0d-9

!    Apical H-ATPase
  ATPHMP = 50.d-9

!    Basolateral PMCA
  PMCA = 0.50d-9

!    Rate of ammoniagenesis
  QNH4 = 0.250d-6

!    Apical H+/HCO2- cotransporter - ADDITION TO AMW MODEL - REMOVED
  pt(:)%dLA(12,13,1,2) = 0.0d-9

!    Apical SGLT2/SGLT1 and GLUT2/GLUT1 expression levels vary in S1-S2 vs S3
  do J = 0, NZ
     if (J < xS3*NZ) then
        pt(J)%xSGLT2 = 1.0d0
        pt(J)%xSGLT1 = 0.0d0
        pt(J)%xGLUT2 = 1.0d0
        pt(J)%xGLUT1 = 0.0d0
     else
        pt(J)%xSGLT2 = 0.0d0
        pt(J)%xSGLT1 = 1.0d0
        pt(J)%xGLUT2 = 0.0d0
        pt(J)%xGLUT1 = 1.0d0
     end if
  end do

!---------------------------------------------------------------------72
! set up parameter to scale torque effect  
  do J = 0, NZ
     if (J < xS3*NZ) then
        pt(J)%scaleT = 1.0
     else
        pt(J)%scaleT = 0.5
     end if
  end do
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	  Non-Dimensional coefficients
!       Divide dLA(I,J,K,L) by (href)*Cref for consistency with other
!	  solute flux terms
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


  Do I = 1,NSPT
     Do J = 1, NSPT
        Do K = 1, NC
           Do L = K, NC
              pt(:)%dLA(I,J,K,L) = pt(:)%dLA(I,J,K,L)/(href*Cref)
           End Do
        End Do
     End Do
  End Do

  pt(:)%xNHE3 = xNHE3/(href*Cref)

  pt(:)%ATPNaK(2,5) = ATPNaKPES/(href*Cref)
  pt(:)%ATPNaK(2,6) = ATPNaKPES/(href*Cref)
  
  pt(:)%ATPH(1,2) = ATPHMP/(href*Cref)

  pt(:)%qnh4 = QNH4/(href*Cref)
  
  pt(:)%CTsglt1 = CTsglt1/(href*Cref)

  pt(:)%CTglut1 = CTglut1/(href*Cref)
  pt(:)%CTglut2 = CTglut2/(href*Cref)

  pt(:)%PMCA = PMCA/(href*Cref)

  xNaPiIIaPT = CTnapiIIa/(href*Cref)
  xNaPiIIbPT = CTnapiIIb/(href*Cref)
  xNaPiIIcPT = CTnapiIIc/(href*Cref)
  xPit2PT = CTpit2/(href*Cref)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     Other kinetic parameters
!	Units of kinetic constants are s-1
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  dkhcat = 1.450d3
  dkdcat = 496.d3

  do J = 0, NZ
     if (J < xS3*NZ) then
        pt(J)%dkd = dkdcat
        pt(J)%dkh = dkhcat
     else  ! reduced reaction rates along the S3
        pt(J)%dkd = dkdcat/100.0
        pt(J)%dkh = dkhcat/100.0
     end if
  end do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     BOUNDARY CONDITIONS	IN PERITUBULAR SOLUTION
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!   Solute concentration in compartment S (mmole/liter)
!	Since reference concentration is 1 mmole/liter = 1d-3 mmol/cm3,
!	there is no need to convert to non-dimensional values


  BathImperm = 2.0d0 !Impermeant in bath (not present in interspace)
    
  xIS = 0.60/2.00d0

  pt(:)%ep(6) = -0.000d-3/EPref  ! qnewton2PT assumes this is 0

!---------------------------------------------------------------------72
!     Account for gradients in peritubular space
!---------------------------------------------------------------------72
! pos is a vector with the the non-dimensional distance from the upper boundary

  do jz = 0, int(xS3*NZ)
     pos(jz) = jz/(xS3*NZ)
  end do
  call set_intconc ( pt, int(xS3*NZ), 1, pos )

  do jz = int(xS3*NZ+1), NZ
     pos(jz) = xIS*(jz-(xS3*NZ))/((1-xS3)*NZ)
  end do
  ind = int(xS3*NZ)+1
  indrev = int(NZ-xS3*NZ)
  call set_intconc ( pt(ind:NZ), indrev, 2, pos(ind:NZ))
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     INITIAL CONDITIONS AT LUMEN ENTRANCE
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  pt(0)%vol(1)=sngfr/Vref ! Normalize flows and volumes by Vref
  pt(:)%volLuminit = pt(0)%vol(1)

  pt(0)%pres = PMinitPT ! initialize inflow pressure

  LumImperm = 0.0d0 ! No impermeant in lumen under normal conditions
  
  do I = 1,NSPT
     pt(0)%conc(I,1) = pt(0)%conc(I,6)
  end do
  pt(0)%ph(1) = pt(0)%ph(6)
  
  pt(0)%ep(1) = -0.16d-3/EPref ! Initial guess for lumen potential

!---------------------------------------------------------------------72
!     Check electroneutrality
!---------------------------------------------------------------------72
  
  elecM = 0.0d0
  elecS = 0.0d0
  elecS_out = 0.d0
  Do I = 1,NSPT
     elecM = elecM + zval(I)*pt(0)%conc(I,1)
     elecS = elecS + zval(I)*pt(0)%conc(I,6)
	 elecS_out = elecS_out +  zval(I)*pt(NZ)%conc(I,6)
  end do 

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!
!     FOR TORQUE-MODULATED EFFECTS ON TRANSCELLULAR TRANSPORTER DENSITY
!
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!     Baseline Torque (Eq. 37 in 2007 PT model)
!	  Radref = DiamPT/2.0d0
  Radref = 0.00200d0/2.0d0
  flowref = 24.0d-6/60.0d0 ! 24 nl/min converted to cm3/s
  factor1 = 8.0*visc*flowref*torqL/(Radref**2)
  factor2 = 1.0 + (torqL+torqd)/Radref + 0.50*((torqL/Radref)**2)
  pt(:)%TM0= factor1*factor2

  ncompl = 1 !Set to zero for non-compliant tubule
  ntorq = 1 !Set to zero in the absence of torque effects

!---------------------------------------------------------------------72
!     For metabolic calculations
!---------------------------------------------------------------------72

  if (.not. bdiabetes)  then
    pt(:)%TQ = 15.0 !TNa-QO2 ratio in normal PT
  else
    pt(:)%TQ = 12.0 !TNa-QO2 ratio in diabetic PT
  end if

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     INITIAL GUESSES FOR CONCENTRATIONS, VOLUMES, AND POTENTIALS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

     open ( unit=11, file='PTresults')
     Do I = 1,4
        read(11,200),Ci1,pt(0)%conc(I,2),pt(0)%conc(I,5)
     end Do
     Do I = 5,NSPT
        read(11,210),Ci1,pt(0)%conc(I,2),pt(0)%conc(I,5)
     end Do
     read(11,200),phi1,pt(0)%ph(2),pt(0)%ph(5)
     read(11,200),Vol1,pt(0)%vol(2),pt(0)%vol(5)
     read(11,200),pt(0)%ep(1),pt(0)%ep(2),pt(0)%ep(5)
     close ( unit=11 )


!---------------------------------------------------------------------72
!	Assign arbitrary values to compartments A (3) and B (4)
!---------------------------------------------------------------------72

  Do I = 1,NSPT
     pt(0)%conc(I,3) = pt(0)%conc(I,2)
     pt(0)%conc(I,4) = pt(0)%conc(I,2)
  End Do
  
  pt(0)%ph(3) = pt(0)%ph(2)
  pt(0)%ph(4) = pt(0)%ph(2)
  
  pt(0)%vol(3) = 0 ! VolAinitmtal
  pt(0)%vol(4) = 0 ! VolBinitmtal

  pt(0)%ep(3) = pt(0)%ep(2)
  pt(0)%ep(4) = pt(0)%ep(2)

!---------------------------------------------------------------------72
!	Define initial PT volume to determine ONC(1) in water flux calculations
!---------------------------------------------------------------------72

  PTinitVol = pt(0)%vol(1)


200 format (6f12.5)
210 format (6e12.5)

  return

end Subroutine initPT
