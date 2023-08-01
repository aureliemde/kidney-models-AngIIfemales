!     This is the initialization routine for the mTAL
!     It sets up several parameters, constants, and initial guess values.

Subroutine initSDL (sdl)

  include 'values.h'
  include 'global.h'
  include 'defs.h'

  ! passed variables
  type(membrane) :: sdl(0:NZ)
  
!  local variables
  double precision Pf(NC,NC)
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
!---------------------------------------------------------------------72      
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	 Data for rat (AMW model, AJP Renal 2005)
!---------------------------------------------------------------------72

  SlumPinitsdl = 2.0d0
  SbasPinitsdl = 2.0d0
  SlatPinitsdl = 10.0d0
  
  SlumEinitsdl = 0.001d0
  SbasEinitsdl = 0.020d0
  
!---------------------------------------------------------------------72
!	 Assign membrane surface area except for basal membrane of E
!---------------------------------------------------------------------72

  sdl(:)%area(1,2)=SlumPinitsdl
  sdl(:)%area(1,3)=SlumPinitsdl !!! Replaced (SlumAinitsld not defined)
  sdl(:)%area(1,4)=SlumPinitsdl !!! Replaced (SlumBinitsld not defined)
  sdl(:)%area(1,5)=SlumEinitsdl
  
  sdl(:)%area(2,5)=SlatPinitsdl
  sdl(:)%area(3,5)=SlatPinitsdl !!! (idem)
  sdl(:)%area(4,5)=SlatPinitsdl !!! (idem)
  
  sdl(:)%area(2,6)=SbasPinitsdl
  sdl(:)%area(3,6)=SbasPinitsdl !!! (idem)
  sdl(:)%area(4,6)=SbasPinitsdl !!! (idem)
  
  Do K = 1,NC-1
     Do L = K+1,NC
        sdl(:)%area(L,K)=sdl(:)%area(K,L)
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


  PfMP = 33.0d-4/2.0d0
  PfMA = 0.0d0
  PfMB = 0.0d0
  PfME = 0.70d-4/0.0010d0
  PfPE = 170.d-4/10.0d0
  PfAE = 0.0d0
  PfBE = 0.0d0
  PfPS = 0.1 * 33.0d-4/2.0d0
  PfAS = 0.0d0
  PfBS = 0.0d0
  PfES = 8000.0d-4/0.020d0

  PfMP = 0.40/36.0
  PfMA = 0.0d0
  PfMB = 0.0d0
  PfME = 0.22/0.0010
  PfPE = 0.40/36.0
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
        sdl(:)%dLPV(K,L) = Pf(K,L)/(Pfref)
     End Do
  End Do

  do jz = 1, NZ	  	   
     if (jz>=NZ*0.46) then   ! terminal water impermeable segment
        sdl(jz)%dLPV(1,2) = 0.001*sdl(jz)%dLPV(1,2)
        sdl(jz)%dLPV(1,5) = 0.001*sdl(jz)%dLPV(1,5)
        sdl(jz)%dLPV(2,5) = 0.001*sdl(jz)%dLPV(2,5)
        sdl(jz)%dLPV(2,6) = 0.001*sdl(jz)%dLPV(2,6)        
     end if
  end do
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     sig(I,K,L) = reflection coefficient of Ith solute between K and L
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	Initialize 
  Do I = 1,NS
     Do K = 1, NC
        Do L = 1, NC
           sdl(:)%sig(I,K,L) = 1.0d0
        End Do
     End Do
  End Do
  
!	Values at the basement membrane (ES)
  Do I = 1,NS	
     sdl(:)%sig(I,5,6) = 0.0d0
     sdl(:)%sig(I,6,5) = 0.0d0
  End Do
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     BOUNDARY CONDITIONS	IN PERITUBULAR SOLUTION
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  sdl(0)%ep(6) = -0.001d-3/EPref 

  xIS = 0.6/2.00                     ! outer-stripe 0.6 mm
  do jz = 0, NZ
     sdl(jz)%ep(6)=sdl(0)%ep(6)
     pos(jz) = xIS + (1-xIS)*jz/NZ
 End do

 call set_intconc ( sdl, NZ, 2, pos )

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     INITIAL CONDITIONS AT LUMEN ENTRANCE
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  open ( unit=12, file='PToutlet' )

  Do I = 1,4
     read(12,200),sdl(0)%conc(I,1),Ci6
  end Do
  Do I = 5,NS
     read(12,210),sdl(0)%conc(I,1),Ci6
  end Do
  read(12,200),sdl(0)%ph(1),phi6
  read(12,200),sdl(0)%ep(1),EPi6
  read(12,210),sdl(0)%vol(1),sdl(0)%pres
  close ( unit=12 )

!---------------------------------------------------------------------72
!     For metabolic calculations
!---------------------------------------------------------------------72

  if (ndiabetes .eq. 0)  then
    sdl(:)%TQ = 15.0 !TNa-QO2 ratio in normal SDL
  else
    sdl(:)%TQ = 12.0 !TNa-QO2 ratio in diabetic SDL
  end if

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


 200   format (6e12.5)
 210   format (6e12.5)

  return
end Subroutine initSDL
      
