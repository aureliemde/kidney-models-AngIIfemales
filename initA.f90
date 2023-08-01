!     This is the initialization routine for the mTAL
!     It sets up several parameters, constants, and initial guess values.

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!   Solute Indices:  1 = Na+, 2 = K+, 3 = Cl-, 4 = HCO3-, 5 = H2CO3, 6 = CO2
!   7 = HPO4(2-), 8 = H2PO4-, 9 = urea, 10 = NH3, 11 = NH4+, 12 = H+
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

Subroutine initA (mtal)

  include 'values.h'
  include 'global.h'
  include 'defs.h'

  ! passed variable
  type(membrane) :: mtal(0:NZ), sdl(0:NZ)
  !  local variables
  double precision Pf(NC,NC),dLA(NS,NS,NC,NC)
  double precision pos(0:NZ)


!---------------------------------------------------------------------72      
!---------------------------------------------------------------------72
!	 Calculate initial surfaces and volumes
!	 Volumes, surfaces, and angles implicitly include the number of
!	 each type of cell (but the PtoIratio doesn't).
!---------------------------------------------------------------------72      
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	 Data for rat (AMW model, AJP Renal 2005)
!---------------------------------------------------------------------72

  SlumPinitmtal = 2.0d0
  SbasPinitmtal = 2.0d0
  SlatPinitmtal = 10.0d0
  
  SlumEinitmtal = 0.001d0
  mtal(:)%sbasEinit = 0.020d0
  
  mtal(:)%volPinit = 5.0d0
  mtal(:)%volEinit = 0.40d0
  
  mtal(:)%volAinit = 5.0d0 ! Needed for water flux subroutine
  mtal(:)%volBinit = 5.0d0 ! Needed for water flux subroutine

!---------------------------------------------------------------------72
!	 Assign membrane surface area except for basal membrane of E
!---------------------------------------------------------------------72

  mtal(:)%area(1,2)=SlumPinitmtal
  mtal(:)%area(1,3)=SlumPinitmtal
  mtal(:)%area(1,4)=SlumPinitmtal
  mtal(:)%area(1,5)=SlumEinitmtal
  
  mtal(:)%area(2,5)=SlatPinitmtal
  mtal(:)%area(3,5)=SlatPinitmtal
  mtal(:)%area(4,5)=SlatPinitmtal
  
  mtal(:)%area(2,6)=SbasPinitmtal
  mtal(:)%area(3,6)=SbasPinitmtal
  mtal(:)%area(4,6)=SbasPinitmtal
  
  Do K = 1,NC-1
     Do L = K+1,NC
        mtal(:)%area(L,K)=mtal(:)%area(K,L)
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
!	The permeabilities have already been multiplied by the area (cm2/cm2 epith)
!   dLPV(K,L) = hydraulic permeability at the K-L interface, in cm/s/mmHg
!	Pf(K,L) = osmotic permeability at the K-L interface, in cm/s
!	Relationship between the two: LPV=Pf*(Vwbar/RT)

! WATER PERMEABILITIES

  PfMP = 33.0d-4/2.0d0*1.0d-3
  PfMA = 0.0d0
  PfMB = 0.0d0
  PfME = 0.70d-4/0.0010d0
  PfPE = 170.d-4/10.0d0
  PfAE = 0.0d0
  PfBE = 0.0d0
  PfPS = 33.0d-4/2.0d0
  PfAS = 0.0d0
  PfBS = 0.0d0
  PfES = 8000.0d-4/0.020d0
  
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
        mtal(:)%dLPV(K,L) = Pf(K,L)/(Pfref)
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
           mtal(:)%sig(I,K,L) = 1.0d0
        End Do
     End Do
  End Do
  
!	Values at the basement membrane (ES)
  Do I = 1,NS	
     mtal(:)%sig(I,5,6) = 0.0d0
     mtal(:)%sig(I,6,5) = 0.0d0
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
           mtal(:)%h(I,K,L) = 0.0d0
        End Do
     End Do
  End Do
  
!	Values at the MP interface (1-2)
  mtal(:)%h(1,1,2) = 0.000d0
  mtal(:)%h(2,1,2) = 20.0d0
  mtal(:)%h(3,1,2) = 0.0d0
  mtal(:)%h(4,1,2) = 0.0d0 
  mtal(:)%h(5,1,2) = 130.d0
  mtal(:)%h(6,1,2) = 1.50d4
  mtal(:)%h(7,1,2) = 0.0d0
  mtal(:)%h(8,1,2) = 0.0d0
  mtal(:)%h(9,1,2) = 0.06d0
  mtal(:)%h(10,1,2) = 1.5d3  
  mtal(:)%h(11,1,2) = 4.0d0
  mtal(:)%h(12,1,2) = 2.0d3 
  mtal(:)%h(13,1,2) = 0.d0
  mtal(:)%h(14,1,2) = 0.d0
  mtal(:)%h(15,1,2) = 0.001d0
  mtal(:)%h(16,1,2) = 0.0001d0

!  if (FAngProx) then
!    mtal(:)%h(2,1,2) = mtal(:)%h(2,1,2)*0.70! Inhibition of medullary ROMK
!  end if

!	Values at the PE interface (2-5)
  mtal(:)%h(1,2,5) = 0.000d0
  mtal(:)%h(2,2,5) = 2.00d0 
  mtal(:)%h(3,2,5) = 0.50d0 
  mtal(:)%h(4,2,5) = 0.10d0
  mtal(:)%h(5,2,5) = 130.d0
  mtal(:)%h(6,2,5) = 1.50d4
  mtal(:)%h(7,2,5) = 0.008d0
  mtal(:)%h(8,2,5) = 0.008d0
  mtal(:)%h(9,2,5) = 0.06d0
  mtal(:)%h(10,2,5) = 1.0d3 
  mtal(:)%h(11,2,5) = 0.40d0
  mtal(:)%h(12,2,5) = 2.0d3 
  mtal(:)%h(13,2,5) = 1.0d-4
  mtal(:)%h(14,2,5) = 1.0d-4
  mtal(:)%h(15,2,5) = 0.001d0
  mtal(:)%h(16,2,5) = 0.0001d0

!	Values at the PS interface (2-6)
  Do I = 1, NS
     mtal(:)%h(I,2,6) = mtal(:)%h(I,2,5)
  End do

!	Values at the MA interface (1-3)
!	Values at the AE interface (3-5)
!	Values at the AS interface (3-6)
!	Values at the MB interface (1-4)
!	Values at the BE interface (4-5)
!	Values at the BS interface (4-6)

  Do I = 1, NS
     mtal(:)%h(I,1,3) = 0.d0
     mtal(:)%h(I,3,5) = 0.d0
     mtal(:)%h(I,3,6) = 0.d0
     mtal(:)%h(I,1,4) = 0.d0
     mtal(:)%h(I,4,5) = 0.d0
     mtal(:)%h(I,4,6) = 0.d0
  End do


!	Values at the ME interface (1-5)
!   Values in Weinstein models are multiplied by the area, 0.001. Divide by 0.001.
  ClTJperm = 1000.0d0
  mtal(:)%h(1,1,5) = 2.80d0*ClTJperm
  mtal(:)%h(2,1,5) = 3.00d0*ClTJperm
  mtal(:)%h(3,1,5) = 1.40d0*ClTJperm
  mtal(:)%h(4,1,5) = 0.50d0*ClTJperm
  mtal(:)%h(5,1,5) = 2.00d0*ClTJperm
  mtal(:)%h(6,1,5) = 2.00d0*ClTJperm
  mtal(:)%h(7,1,5) = 0.20d0*ClTJperm
  mtal(:)%h(8,1,5) = 0.20d0*ClTJperm
  mtal(:)%h(9,1,5) = 0.40d0*ClTJperm
  mtal(:)%h(10,1,5) = 6.00d0*ClTJperm
  mtal(:)%h(11,1,5) = 3.00d0*ClTJperm
  mtal(:)%h(12,1,5) = 6.00d0*ClTJperm
  mtal(:)%h(13,1,5) = 0.01d0*ClTJperm
  mtal(:)%h(14,1,5) = 0.01d0*ClTJperm
  mtal(:)%h(15,1,5) = 0.01d0*ClTJperm
  mtal(:)%h(16,1,5) = 4.20d0*ClTJperm

! if (FAngProx) then
!       mtal(:)%h(1,1,5) = mtal(:)%h(1,1,5)*0.54 ! mClaudin 10 x 0.54 in grain-fed
!       mtal(:)%h(3,1,5) = mtal(:)%h(3,1,5)*0.54 ! mClaudin 10 x 0.54 in grain-fed
! end if

!	Values at the ES interface (5-6)
  areafactor = 1.d0/0.02d0
  mtal(:)%h(1,5,6) = 72.14d0*areafactor
  mtal(:)%h(2,5,6) = 96.18d0*areafactor
  mtal(:)%h(3,5,6) = 96.18d0*areafactor
  mtal(:)%h(4,5,6) = 48.09d0*areafactor
  mtal(:)%h(5,5,6) = 72.14d0*areafactor
  mtal(:)%h(6,5,6) = 72.14d0*areafactor
  mtal(:)%h(7,5,6) = 48.09d0*areafactor
  mtal(:)%h(8,5,6) = 48.09d0*areafactor
  mtal(:)%h(9,5,6) = 48.09d0*areafactor
  mtal(:)%h(10,5,6) = 60.11d0*areafactor
  mtal(:)%h(11,5,6) = 96.18d0*areafactor
  mtal(:)%h(12,5,6) = 480.90d0*areafactor
  mtal(:)%h(13,5,6) = 1.0d0*areafactor
  mtal(:)%h(14,5,6) = 1.0d0*areafactor
  mtal(:)%h(15,5,6) = 1.0d0*areafactor
  mtal(:)%h(16,5,6) = mtal(:)%h(1,5,6)*(7.93/13.3)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Non-Dimensional Solute permeabilities
!     Divide h(I,K,L) by (href) 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  Do I = 1,NS
     Do K = 1, NC
        Do L = K, NC
           mtal(:)%h(I,K,L) = mtal(:)%h(I,K,L)*1.0d-5/(href)
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

!	 Cl/HCO3 exchanger on basolateral membrane
  dLA(3,4,2,5)=3.0d-9
  dLA(3,4,2,6)=dLA(3,4,2,5)

!	 Na/H exchangers on basolateral membrane	 
  dLA(1,12,2,5)=10.0d-9
  dLA(1,12,2,6)=dLA(1,12,2,5)

!	 Na2/HPO4 co-transporters on basolateral membrane	 
  dLA(1,7,2,5)=0.50d-9
  dLA(1,7,2,6)=dLA(1,7,2,5)

!	 Na/3HCO3 exchanger on basolateral membrane
  dLA(1,4,2,5)=0.50d-9
  dLA(1,4,2,6)=dLA(1,4,5,6)

!	 Apical NKCC2
!!! FEMALE ADJUSTMENT - mNKCC2 x 1.30
  xNKCC2A = 1.30*8.0d-9
  xNKCC2F = 1.30*50.0d-9

  if (FAngProx) then
     xNKCC2A = xNKCC2A*0.50 ! NKCC2 x 0.43 - NKCC2p-med x 1.13 in grain-fed
     xNKCC2F = xNKCC2F*0.50 ! NKCC2 x 0.43 - NKCC2p-med x 1.13 in grain-fed
  end if

!	 Apical NHE3 
  xNHE3 = 6.0d-9 

 if (FAngProx) then
     xNHE3 = xNHE3*0.62 ! NHE3p-med x 0.62 in grain fed
  end if

!	 Basolateral KCC4
  xKCC4 = 0.70d-9

!	 Basolateral Na-K-ATPase 
!!! FEMALE ADJUSTMENT - mNKCC2 x 1.3 and NKA x 1.3
  ATPNaKPES = 1.30*1300.d-9

  if (FAngProx) then
     ATPNaKPES = ATPNaKPES*0.77 ! NKCC2 x 0.43 and mNKA x 0.77 in grain-fed
  end if

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	  Non-Dimensional coefficients
!       Divide dLA(I,J,K,L) by (href)*Cref for consistency with other
!	  solute flux terms
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
  
  Do I = 1,NS
     Do J = 1, NS
        Do K = 1, NC
           Do L = K, NC
              mtal(:)%dLA(I,J,K,L) = dLA(I,J,K,L)/(href*Cref)
           End Do
        End Do
     End Do
  End Do

  mtal(:)%xNKCC2F = xNKCC2F/(href*Cref)
  mtal(:)%xNKCC2A = xNKCC2A/(href*Cref)
  
  mtal(:)%xNHE3 = xNHE3/(href*Cref)
  mtal(:)%xKCC4 = xKCC4/(href*Cref)

  mtal(:)%ATPNaK(2,5) = ATPNaKPES/(href*Cref)
  mtal(:)%ATPNaK(2,6) = ATPNaKPES/(href*Cref)
  
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     Other kinetic parameters
!	Units of kinetic constants are s-1
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  dkhcat = 1.450d3
  dkdcat = 496.d3

  mtal(:)%dkd(1)=dkdcat
  mtal(:)%dkh(1)=dkhcat
  mtal(:)%dkd(2)=dkdcat
  mtal(:)%dkh(2)=dkhcat
  mtal(:)%dkd(5)=dkdcat
  mtal(:)%dkh(5)=dkhcat

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     BOUNDARY CONDITIONS IN PERITUBULAR SOLUTION
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!   Solute concentration in compartment S (mmole/liter)
!	Since reference concentration is 1 mmole/liter = 1d-3 mmol/cm3,
!	there is no need to convert to non-dimensional values

  mtal(0)%EP(6) = -0.001d-3/EPref 

  Do jz=0,NZ
     mtal(jz)%ep(6)=mtal(0)%EP(6)
     pos(jz) = 1.0d0*(NZ-jz)/NZ ! Non-dimensional position along mTAL
  End do

  call set_intconc ( mtal, NZ, 2, pos )
  

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     INITIAL CONDITIONS AT LUMEN ENTRANCE
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

50   open ( unit=12, file='SDLoutlet' )

  Do I = 1,4
     read(12,200),mtal(0)%conc(I,1),Ci6
  end Do
  Do I = 5,NS
     read(12,210),mtal(0)%conc(I,1),Ci6
  end Do
  read(12,200),mtal(0)%ph(1),phi6
  read(12,200),mtal(0)%EP(1),EPi6
  read(12,210),mtal(0)%vol(1),mtal(0)%pres
  close ( unit=12 )

  mtal(:)%volLuminit = mtal(0)%vol(1)

!---------------------------------------------------------------------72
!     For metabolic calculations
!---------------------------------------------------------------------72

  if (ndiabetes .eq. 0)  then
    mtal(:)%TQ = 15.0 !TNa-QO2 ratio in normal mTAL
  else
    mtal(:)%TQ = 12.0 !TNa-QO2 ratio in diabetic cTAL
  end if

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

 200   format (6f12.5)
 210   format (6e12.5)

  return
end Subroutine initA
      
