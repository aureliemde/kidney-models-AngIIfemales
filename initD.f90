!     This is the initialization routine for the DCT.
!     It sets up several parameters, constants, and initial guess values.

Subroutine initD (dct)

  include 'values.h'
  include 'global.h'
  include 'defs.h'

  ! passed variables
  type (membrane) :: dct(0:NZ)
  
!  local variables
  double precision Pf(NC,NC),dLA(NS,NS,NC,NC)
  double precision pos(0:NZ)

!---------------------------------------------------------------------72      
!---------------------------------------------------------------------72
!	 Calculate initial surfaces and volumes
!	 Volumes, surfaces, and angles implicitly include the number of
!	 each type of cell (but the PtoIratio doesn't).
!
!	 Data for rat (AMW model, AJP Renal 2005)
!---------------------------------------------------------------------72

  SlumPinitdct = 4.7d0
  SbasPinitdct = 4.7d0
  SlatPinitdct = 64.3d0

  SlumEinitdct = 0.001d0
  dct(:)%sbasEinit = 0.020d0
  
  VolPinitdct = 7.5d0
  VolEinitdct = 0.80d0

  dct(:)%volPinit = 7.5d0
  dct(:)%volEinit = 0.80d0

  dct(:)%volAinit = 7.5d0
  dct(:)%volBinit = 7.5d0

  
  SlumAinitdct = 1.0d0
  SbasAinitdct = 1.0d0
  SlatAinitdct = 1.0d0
  SlumBinitdct = 1.0d0
  SbasBinitdct = 1.0d0
  SlatBinitdct = 1.0d0
  VolAinitdct = 1.0d0
  VolBinitdct = 1.0d0
  
!---------------------------------------------------------------------72
!	 Assign membrane surface area except for basal membrane of E
!---------------------------------------------------------------------72

  dct(:)%area(1,2)=SlumPinitdct
  dct(:)%area(1,3)=SlumAinitdct
  dct(:)%area(1,4)=SlumBinitdct
  dct(:)%area(1,5)=SlumEinitdct
  
  dct(:)%area(2,5)=SlatPinitdct
  dct(:)%area(3,5)=SlatAinitdct
  dct(:)%area(4,5)=SlatBinitdct
  
  dct(:)%area(2,6)=SbasPinitdct
  dct(:)%area(3,6)=SbasAinitdct
  dct(:)%area(4,6)=SbasBinitdct
  
  Do K = 1,NC-1
     Do L = K+1,NC
        dct(:)%area(L,K)=dct(:)%area(K,L)
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

!!! - FEMALE ADJUSTMENT Pf x 2.4
  PfMP = 2.4*0.00117d0
  PfMA = 0.0d0
  PfMB = 0.0d0
  PfME = 2.0d0
  PfPE = 2.4*0.00835d0
  PfAE = 0.0d0
  PfBE = 0.0d0
  PfPS = 2.4*0.00835d0
  PfAS = 0.0d0
  PfBS = 0.0d0
  PfES = 35.5d0
  
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
        dct(:)%dLPV(K,L) = Pf(K,L)/(Pfref)
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
           dct(:)%sig(I,K,L) = 1.0d0
        End Do
     End Do
  End Do
  

!	Values at the basement membrane (ES)
  Do I = 1,NS	
     dct(:)%sig(I,5,6) = 0.0d0
     dct(:)%sig(I,6,5) = 0.0d0
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
           dct(:)%h(I,K,L) = 0.0d0
        End Do
     End Do
  End Do
  
  !	Values at the MP interface (1-2)
  dct(:)%hNaMP = 14.00d0*0.25d0 ! maximum ENaC perm

  dct(:)%h(1,1,2) = 0.072d0
  dct(:)%h(2,1,2) = 0.60d0
  dct(:)%h(3,1,2) = 0.0d0
  dct(:)%h(4,1,2) = 0.0d0 
  dct(:)%h(5,1,2) = 130.d0
  dct(:)%h(6,1,2) = 1.50d4
  dct(:)%h(7,1,2) = 0.0d0
  dct(:)%h(8,1,2) = 0.0d0
  dct(:)%h(9,1,2) = 0.20d0
  dct(:)%h(10,1,2) = 200.d0  
  dct(:)%h(11,1,2) = 0.12d0
  dct(:)%h(12,1,2) = 0.20d0 
  dct(:)%h(13,1,2) = 0.0d0
  dct(:)%h(14,1,2) = 0.0d0
  dct(:)%h(15,1,2) = 0.0d0
  dct(:)%h(16,1,2) = 0.0d0


!	Values at the PE interface (2-5)
  dct(:)%h(1,2,5) = 0.000d0
  dct(:)%h(2,2,5) = 0.12d0 
  dct(:)%h(3,2,5) = 0.04d0 
  dct(:)%h(4,2,5) = 0.02d0
  dct(:)%h(5,2,5) = 130.d0
  dct(:)%h(6,2,5) = 1.50d4
  dct(:)%h(7,2,5) = 0.002d0
  dct(:)%h(8,2,5) = 0.002d0
  dct(:)%h(9,2,5) = 0.20d0
  dct(:)%h(10,2,5) = 200.d0 
  dct(:)%h(11,2,5) = 0.0234d0
  dct(:)%h(12,2,5) = 0.20d0 
  dct(:)%h(13,2,5) = 1.0d-4
  dct(:)%h(14,2,5) = 1.0d-4
  dct(:)%h(15,2,5) = 1.0d-4
  dct(:)%h(16,2,5) = 0.0d0

!	Values at the PS interface (2-6)
  Do I = 1, NS
     dct(:)%h(I,2,6) = dct(:)%h(I,2,5)
  End do
  
!	Values at the MA interface (1-3)
!	Values at the AE interface (3-5)
!	Values at the AS interface (3-6)
!	Values at the MB interface (1-4)
!	Values at the BE interface (4-5)
!	Values at the BS interface (4-6)

  Do I = 1, NS
     dct(:)%h(I,1,3) = 0.d0
     dct(:)%h(I,3,5) = 0.d0
     dct(:)%h(I,3,6) = 0.d0
     dct(:)%h(I,1,4) = 0.d0
     dct(:)%h(I,4,5) = 0.d0
     dct(:)%h(I,4,6) = 0.d0
  End do


!	Values at the ME interface (1-5)
  ClTJperm = 1000.0d0 ! Factor for area = 1/0.001
  !!! - FEMALE ADJUSTMENT related to Claudin expression
  dct(:)%h(1,1,5) = 0.80d0*ClTJperm
  dct(:)%h(2,1,5) = 0.80d0*ClTJperm
  dct(:)%h(3,1,5) = 1.30*0.50d0*ClTJperm ! Cldn7
  dct(:)%h(4,1,5) = 0.50d0*ClTJperm
  dct(:)%h(5,1,5) = 0.50d0*ClTJperm
  dct(:)%h(6,1,5) = 0.50d0*ClTJperm
  dct(:)%h(7,1,5) = 0.10d0*ClTJperm
  dct(:)%h(8,1,5) = 0.10d0*ClTJperm
  dct(:)%h(9,1,5) = 0.20d0*ClTJperm
  dct(:)%h(10,1,5) = 0.80d0*ClTJperm
  dct(:)%h(11,1,5) = 0.80d0*ClTJperm
  dct(:)%h(12,1,5) = 0.80d0*ClTJperm
  dct(:)%h(13,1,5) = 0.01d0*ClTJperm
  dct(:)%h(14,1,5) = 0.01d0*ClTJperm
  dct(:)%h(15,1,5) = 0.01d0*ClTJperm
  dct(:)%h(16,1,5) = 0.0001d0*ClTJperm ! TO BE ADJUSTED
  

!	Values at the ES interface (5-6)
  areafactor = 1.d0/0.02d0
  dct(:)%h(1,5,6) = 63.d0*areafactor
  dct(:)%h(2,5,6) = 84.d0*areafactor
  dct(:)%h(3,5,6) = 84.d0*areafactor
  dct(:)%h(4,5,6) = 42.d0*areafactor
  dct(:)%h(5,5,6) = 63.d0*areafactor
  dct(:)%h(6,5,6) = 63.d0*areafactor
  dct(:)%h(7,5,6) = 42.d0*areafactor
  dct(:)%h(8,5,6) = 42.d0*areafactor
  dct(:)%h(9,5,6) = 42.d0*areafactor
  dct(:)%h(10,5,6) = 52.d0*areafactor
  dct(:)%h(11,5,6) = 84.0d0*areafactor
  dct(:)%h(12,5,6) = 419.0d0*areafactor
  dct(:)%h(13,5,6) = 1.0d0*areafactor
  dct(:)%h(14,5,6) = 1.0d0*areafactor
  dct(:)%h(15,5,6) = 1.0d0*areafactor
  dct(:)%h(16,5,6) = dct(:)%h(1,5,6)*(7.93/13.3)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Non-Dimensional Solute permeabilities
!     Divide h(I,K,L) by (href) 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  Do I = 1,NS
     Do K = 1, NC
        Do L = K, NC
           dct(:)%h(I,K,L) = dct(:)%h(I,K,L)*1.0d-5/(href)
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
  
!	 Apical NHE2/3 in PC
  xNHE3 = 200.d-9

!	 Apical NCC in PC
!!! - FEMALE ADJUSTMENT NCC x 1.80
  xNCC = 1.80*15.0d-9

  if (FAngDist) then
     xNCC = xNCC*2.47 ! NCCp-cor x 2.47 in grain-fed
  end if

!	 Apical and  basolateral K-Cl in PC
  dLA(2,3,1,2)=4.0d-9
  dLA(2,3,2,5)=20.0d-9
  dLA(2,3,2,6)=20.0d-9
  
!	 Basolateral Na/H exchanger
  dLA(1,12,2,5)=4.0d-9 
  dLA(1,12,2,6)=4.0d-9 
  
  xNHE1 = 6.95d-10

!	 Basolateral Cl/HCO3 exchanger  
  dLA(3,4,2,5)=15.0d-9
  dLA(3,4,2,6)=15.0d-9
  
!	 Basolateral Na/HPO4 co-transporter	 
  dLA(1,7,2,5)=0.2d-9
  dLA(1,7,2,6)=0.2d-9
  
!	 Na-K-ATPase (basolateral in PC and IC)
!!! - FEMALE ADJUSTMENT NCC and NKA x 1.80
  ATPNaKPES = 1.80*400.d-9

!---------------------------------------------------------------------72
! Specific calcium transporters
!---------------------------------------------------------------------72
 !  Basolateral NCX exchanger
  xNCX = 50.0d-9*1.40

 ! Basolateral PMCA pump
  PMCA = 0.70d-9*0.60

 ! Apical TRPV5
  xTRPV5 = 13.5d6

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
              dct(:)%dLA(I,J,K,L) = dLA(I,J,K,L)/(href*Cref)
           End Do
        End Do
     End Do
  End Do
  
  dct(:)%xNCC = xNCC/(href*Cref)
  dct(:)%xNHE3 = xNHE3/(href*Cref)
  
  dct(:)%xNCX = xNCX/(href*Cref)
  dct(:)%PMCA= PMCA/(href*Cref)
  xTRPV5_dct = xTRPV5/(href*Cref)

  dct(:)%ATPNaK(2,5) = ATPNaKPES/(href*Cref)
  dct(:)%ATPNaK(2,6) = ATPNaKPES/(href*Cref)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     Other kinetic parameters
!	Units of kinetic constants are s-1
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  dKhuncat = 0.145d0
  dKduncat = 49.6d0
  dkhP = 145.d0
  dkdP = 49600.d0
  dkhE = 145.d0
  dkdE = 49600.d0
  
  dct(:)%dkd(1)=dkduncat*10.0d0
  dct(:)%dkh(1)=dkhuncat*10.0d0
  dct(:)%dkd(2)=dkdP
  dct(:)%dkh(2)=dkhP
  dct(:)%dkd(5)=dkdE
  dct(:)%dkh(5)=dkhE

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	READ ENTERING AND PERITUBULAR CONDITIONS FROM TAL OUTLET FILE
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  open ( unit=22, file='cTALoutlet' )

  Do I = 1,4
     read(22,200),dct(0)%conc(I,1),dct(0)%conc(I,6)
  end Do
  Do I = 5,NS
     read(22,210),dct(0)%conc(I,1),dct(0)%conc(I,6)
  end Do
  read(22,200),dct(0)%ph(1),dct(0)%ph(6)
  read(22,200),dct(0)%ep(1),dct(0)%ep(6)
  read(22,210),dct(0)%vol(1),dct(0)%pres
  close ( unit=22 )
  dct(:)%volLuminit = dct(0)%vol(1)
  dct(:)%ep(6)=dct(0)%ep(6)

  pos(:) = 0.0d0 ! twirls around the cortical surface

  call set_intconc ( dct, NZ, 1, pos )
  
!---------------------------------------------------------------------72
!     For metabolic calculations
!---------------------------------------------------------------------72

  if (ndiabetes .eq. 0)  then
    dct(:)%TQ = 15.0 !TNa-QO2 ratio in normal DCT
  else
    dct(:)%TQ = 12.0 !TNa-QO2 ratio in diabetic DCT
  end if

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
 
200 format (6f12.5)
210 format (6e12.5)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  return
end Subroutine initD
      
