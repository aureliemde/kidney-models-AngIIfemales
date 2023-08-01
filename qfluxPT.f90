!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!      This is a subroutine to compute PT fluxes	

Subroutine qfluxPT(x,Jvol,Jsol,Cext,EPext,pt,vol0,dcompl,dtorq)

! INPUT ARGUMENTS:
!	 x: vector of 4*NSPT concentrations, 4 volumes, and 4 potentials
!
! OUTPUT ARGUMENTS:
!     Jvol:    volume fluxes
!     Jsol:    solute fluxes

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  include 'values.h'
  include 'global.h'
  include 'defs.h'

  type (membrane) :: pt
  double precision x(NDPT),Cext(NSPT),EPext
  double precision Jvol(NC,NC), Jsol(NSPT,NC,NC)
  double precision :: dcompl  ! compliant PT (>0)
  double precision :: dtorq   ! torque effects on transport (>0)
  double precision :: vol0    ! inflow volume

  double precision ONC(NC), PRES(NC)
  double precision C(NSPT,NC),dmu(NSPT,NC),ph(NC),Vol(NC),EP(NC)
  double precision delmu(NSPT,NC,NC)
  double precision hkconc(4),Amat(Natp,Natp)
  double precision nagluparam(7),fluxsglt1(2),fluxsglt2(2)

!	for HKATPase matrix inversion
  integer info, lwork
  parameter (lwork=10000)
  integer ipiv(Natp)
  double precision work(lwork)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	ASSIGN CONCENTRATIONS, VOLUMES, AND POTENTIALS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  Do J = 1, NSPT
     C(J,6) = Cext(J) 
  End Do
  EP(6) = EPext 

  Do I = 1,NSPT
     C(I,1)=x(1+3*(I-1))
     C(I,2)=x(2+3*(I-1))
     C(I,5)=x(3+3*(I-1))
  End do

  Vol(1) = x(1+3*NSPT) 
  Vol(2) = x(2+3*NSPT) 
  Vol(5) = x(3+3*NSPT) 
  EP(1) = x(4+3*NSPT) 
  EP(2) = x(5+3*NSPT) 
  EP(5) = x(6+3*NSPT) 
  PM = x(7+3*NSPT)


!	 Assign dummy values to compartments 3 and 4
  Do I = 1,NSPT
     C(I,3) = C(I,2)
     C(I,4) = C(I,2)
  End do

  Do K = 1,NC
     ph(K)=-dlog(C(12,K)/1.0d3)/dlog(10.0d0)
  End do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!		DETERMINE THE NEW SURFACE AREAS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  pt%area(5,6)=pt%sbasEinit*max(Vol(5)/pt%volEinit,1.0d0)
  pt%area(6,5)=pt%area(5,6)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	INITIALIZE ALL FLUXES 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  Do K = 1, NC
     Do L = 1, NC
        Jvol(K,L) = 0.0d0
        Do I = 1, NSPT
           Jsol(I,K,L) = 0.0d0
        End Do
     End Do
  End Do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE WATER FLUXES
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  call compute_water_fluxes (C,PM,PbloodPT,Vol,pt%volLuminit,pt%volEinit, &
       pt%volPinit,CPimprefPT,pt%volAinit,CAimprefPT,pt%volBinit,&
	   CBimprefPT,pt%area,pt%sig,pt%dLPV,complPT,Jvol)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     COMPUTE SOLUTE FLUXES
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!     FOR TORQUE-MODULATED EFFECTS
!---------------------------------------------------------------------72
!  Compliant tubule radius (Eq. 38 in 2007 PT model)
  if (dcompl < 0) then
     RMtorq = pt%diam/2.0d0 !non-compliant tubule
  else
     RMtorq = torqR*(1.0d0+torqvm*(PM - PbloodPT))
  endif

! Torque (Eq. 37 in 2007 PT model)
  factor1 = 8.0*visc*(Vol(1)*Vref)*torqL/(RMtorq**2)
  factor2 = 1.0 + (torqL+torqd)/RMtorq + 0.50*((torqL/RMtorq)**2)
  Torque = factor1*factor2

! Torque scaling parameter (Eq. 39 in 2007 PT model)
  if (dtorq < 0) then
     Scaletorq = 1.0 !No torque effect on transporter density
  else
     Scaletorq = 1.0 + TS*pt%scaleT*(Torque/pt%TM0-1.0)
     PTtorque(LzPT+1) = Scaletorq
  endif


!	 Non-dimensional solute fluxes are first converted to dimensional fluxes
!      (in mmol/s/cm2 epith) by multiplying by href*Cref
!	 Then they are converted to units of pmol/min/mm tubule
  convert=href*Cref*PI*pt%diam*60/10*1.0d9*Scaletorq !Conversion factor

!---------------------------------------------------------------------72
!     ELECTRO-CONVECTIVE-DIFFUSIVE FLUXES
!---------------------------------------------------------------------72

  call compute_ecd_fluxes (C,EP,pt%area,pt%sig,pt%h,Jvol,Jsol,delmu)


!---------------------------------------------------------------------72	
!---------------------------------------------------------------------72
!	 FLUXES ACROSS EXCHANGERS AND COTRANSPORTERS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	 Luminal Na/glucose cotransporters 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!   Compute SGLT1 fluxes
!---------------------------------------------------------------------72

!   Use kinetic description of SGLT1
  nagluparam(1) = C(1,1)
  nagluparam(2) = C(1,2)
  nagluparam(3) = C(15,1)
  nagluparam(4) = C(15,2)
  nagluparam(5) = EP(1)
  nagluparam(6) = EP(2)
  n = 1
  nagluparam(7) = pt%CTsglt1
  call sglt(n,nagluparam,pt%area(1,2),fluxsglt1)
  Jsol(1,1,2)=Jsol(1,1,2)+fluxsglt1(1)*pt%xSGLT1  
  Jsol(15,1,2)=Jsol(15,1,2)+fluxsglt1(2)*pt%xSGLT1   
  fluxnasglt1=fluxsglt1(1)*pt%xSGLT1*convert
  fluxglusglt1=fluxsglt1(2)*pt%xSGLT1*convert


!---------------------------------------------------------------------72
!   Compute SGLT2 fluxes
!---------------------------------------------------------------------72
!   NET formulation from Alan's 2007 PT model
!---------------------------------------------------------------------72
!	dJNaglu1=SAreapt2(1,2)*dLApt(1,15,1,2)*(delmu(1,1,2)+delmu(15,1,2)) 

  zeta = F*EPref*(EP(1)-EP(2))/RT
  affglu=4.90d0
  affna=25.0d0
  snal = C(1,1)/affna
  snac = C(1,2)/affna
  glul = C(15,1)/affglu
  gluc = C(15,2)/affglu
  denom = (1.0 + snal + glul + snal*glul)*(1.0 + snac*gluc) + &
       (1.0 + snac + gluc + snac*gluc)*(1.0 + snal*glul*exp(zeta))
  dJNaglu=pt%area(1,2)*pt%dLA(1,15,1,2)*(C(1,1)*C(15,1)*exp(zeta) - &
       C(1,2)*C(15,2))/denom

  Jsol(1,1,2)=Jsol(1,1,2)+dJNaglu*pt%xSGLT2  
  Jsol(15,1,2)=Jsol(15,1,2)+dJNaglu*pt%xSGLT2 
  fluxnasglt2=dJNaglu*pt%xSGLT2*convert
  fluxglusglt2=dJNaglu*pt%xSGLT2*convert


!---------------------------------------------------------------------72
!	 Basolateral glucose transporters 
!---------------------------------------------------------------------72

  Gi = C(15,2)
  Go5 = C(15,5)
  Go6 = C(15,6)

!   Compute GLUT1 fluxes
!---------------------------------------------------------------------72
  
  affglut1 = 2.0
  Ro5 = affglut1*(Gi-Go5)/(affglut1+Gi)/(affglut1+Go5)
  fluxglut1PE = pt%CTglut1*pt%area(2,5)*Ro5
  Ro6 = affglut1*(Gi-Go6)/(affglut1+Gi)/(affglut1+Go6)
  fluxglut1PS = pt%CTglut1*pt%area(2,6)*Ro6
  Jsol(15,2,5) = Jsol(15,2,5) + fluxglut1PE*pt%xGLUT1
  Jsol(15,2,6) = Jsol(15,2,6) + fluxglut1PS*pt%xGLUT1
  fluxglut1 = (fluxglut1PE+fluxglut1PS)*pt%xGLUT1*convert

!   Compute GLUT2 fluxes
!---------------------------------------------------------------------72

  affglut2 = 17.0
  Ro5 = affglut2*(Gi-Go5)/(affglut2+Gi)/(affglut2+Go5)
  fluxglut2PE = pt%CTglut2*pt%area(2,5)*Ro5
  Ro6 = affglut2*(Gi-Go6)/(affglut2+Gi)/(affglut2+Go6)
  fluxglut2PS = pt%CTglut2*pt%area(2,6)*Ro6
  Jsol(15,2,5) = Jsol(15,2,5) + fluxglut2PE*pt%xGLUT2
  Jsol(15,2,6) = Jsol(15,2,6) + fluxglut2PS*pt%xGLUT2
  fluxglut2 = (fluxglut2PE+fluxglut2PS)*pt%xGLUT2*convert

!---------------------------------------------------------------------72
!	 Luminal NaH2PO4 cotransporter 
!---------------------------------------------------------------------72

  dJNaP=pt%area(1,2)*pt%dLA(1,8,1,2)*(delmu(1,1,2)+delmu(8,1,2)) 
!  Jsol(1,1,2)=Jsol(1,1,2)+dJNaP
!  Jsol(8,1,2)=Jsol(8,1,2)+dJNaP
  fluxNaP=dJNaP*convert

!   This Na-Pi cotransporter represents NaPiIIa: 3 Na+ with 1 HPO4^{2-}
    dJNaPiIIa=pt%area(1,2)*xNaPiIIaPT*(3*delmu(1,1,2)+delmu(7,1,2))
    Jsol(1,1,2)=Jsol(1,1,2)+3.0d0*dJNaPiIIa
    Jsol(7,1,2)=Jsol(7,1,2)+dJNaPiIIa
    fluxNaPiIIa=dJNaPiIIa*convert

!   This Na-Pi cotransporter represents NaPiIIc: 2 Na+ with 1 HPO4^{2-}
!   xNaPiIIcPT is used owing to the non-homogeneous distribution of NaPiII-c
    dJNaPiIIc=pt%area(1,2)*xNaPiIIcPT*(2*delmu(1,1,2)+delmu(7,1,2))
    dJNaPiIIc=dJNaPiIIc*pt%xSGLT2
    Jsol(1,1,2)=Jsol(1,1,2)+2.0d0*dJNaPiIIc
    Jsol(7,1,2)=Jsol(7,1,2)+dJNaPiIIc
    fluxNaPiIIc=dJNaPiIIc*convert

!   This Na-Pi cotransporter represents PiT-2: 2 Na+ with 1 H2PO4-
!   xPit2PT is used owing to the non-homoegeneous distribution of Pit-2
    dJPit2=pt%area(1,2)*xPit2PT*(2*delmu(1,1,2)+delmu(8,1,2))
    dJPit2=dJPit2*pt%xSGLT2
    Jsol(1,1,2)=Jsol(1,1,2)+2.0d0*dJPit2
    Jsol(8,1,2)=Jsol(8,1,2)+dJPit2
    fluxPit2=dJPit2*convert

!---------------------------------------------------------------------72
!	 Luminal Cl/HCO3 exchanger 
!---------------------------------------------------------------------72

  dJClBic=pt%area(1,2)*pt%dLA(3,4,1,2)*(delmu(3,1,2)-delmu(4,1,2)) 
  Jsol(3,1,2)=Jsol(3,1,2)+dJClBic
  Jsol(4,1,2)=Jsol(4,1,2)-dJClBic
  fluxClBic=dJClBic*convert

!---------------------------------------------------------------------72
!	 Luminal Cl/HCO2 exchanger 
!---------------------------------------------------------------------72

  dJClHco2=pt%area(1,2)*pt%dLA(3,13,1,2)*(delmu(3,1,2)-delmu(13,1,2)) 
  Jsol(3,1,2)=Jsol(3,1,2)+dJClHco2
  Jsol(13,1,2)=Jsol(13,1,2)-dJClHco2
  fluxClHco2=dJClHco2*convert

!---------------------------------------------------------------------72
!	 Luminal H-HCO2 cotransporter - ADDITION TO AMW MODEL
!---------------------------------------------------------------------72

  dJH_Hco2=pt%area(1,2)*pt%dLA(12,13,1,2)*(delmu(12,1,2)+delmu(13,1,2)) 
  Jsol(12,1,2)=Jsol(12,1,2)+dJH_Hco2
  Jsol(13,1,2)=Jsol(13,1,2)+dJH_Hco2
  fluxH_Hco2=dJH_Hco2*convert

!---------------------------------------------------------------------72
!	 Basolateral KCl cotransporter 
!---------------------------------------------------------------------72

  sumJES=0.0d0
  Do L = 5,6
     dJKCl=pt%area(2,L)*pt%dLA(2,3,2,L)*(delmu(2,2,L)+delmu(3,2,L)) 
     Jsol(2,2,L)=Jsol(2,2,L)+dJKCl
     Jsol(3,2,L)=Jsol(3,2,L)+dJKCl
     sumJES=sumJES+dJKCl
  End do
  fluxKCl=sumJES*convert

!---------------------------------------------------------------------72
!	 Basolateral Na(1)-HCO3(3) cotransporter 
!---------------------------------------------------------------------72

  sumJES = 0.d0
  Do L = 5,6
     dJNaBic=pt%area(2,L)*pt%dLA(1,4,2,L)*(delmu(1,2,L)+3*delmu(4,2,L))
     Jsol(1,2,L)=Jsol(1,2,L)+dJNaBic
     Jsol(4,2,L)=Jsol(4,2,L)+3*dJNaBic
     sumJES=sumJES+dJNaBic
  End do
  fluxNaBic = sumJES*convert

!---------------------------------------------------------------------72
!	 Basolateral Na(1)-HCO3(2)/Cl(1) cotransporter (NDCBE)
!---------------------------------------------------------------------72

  sumJES = 0.d0
  Do L = 5,6
     dJNDCBE=pt%area(2,L)*pt%dLA(1,3,2,L)*(delmu(1,2,L) - &
     delmu(3,2,L) + 2*delmu(4,2,L))
     Jsol(1,2,L)=Jsol(1,2,L)+dJNDCBE
     Jsol(3,2,L)=Jsol(3,2,L)-dJNDCBE
     Jsol(4,2,L)=Jsol(4,2,L)+2*dJNDCBE
     sumJES=sumJES+dJNDCBE
  End do
  fluxNDCBE = sumJES*convert

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	 FLUXES ACROSS EXCHANGERS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!		NHE3 EXCHANGER AT LUMINAL MEMBRANE OF CELL
!		n is for Na, c is for Cl
!		p (or prime) is for the luminal compartment (M)
!	    pp (or double prime) is for the cytosolic compartment (I)
!---------------------------------------------------------------------72

  call compute_nhe3_fluxes (C,pt%area(1,2),pt%xNHE3,dJNHEsod,dJNHEprot,dJNHEamm)
  
!	  The rate constants are different for the PT. They need to be divided by 8000./792.
  dJNHEsod = dJNHEsod*792.d0/8000.d0
  dJNHEprot = dJNHEprot*792.d0/8000.d0
  dJNHEamm = dJNHEamm*792.d0/8000.d0
  Jsol(1,1,2)=Jsol(1,1,2)+dJNHEsod
  Jsol(12,1,2)=Jsol(12,1,2)+dJNHEprot
  Jsol(11,1,2)=Jsol(11,1,2)+dJNHEamm
  fluxNHEsod=dJNHEsod*convert
  fluxNHEprot=dJNHEprot*convert
  fluxNHEamm=dJNHEamm*convert
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	ACTIVE TRANSPORT VIA ATPases
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	Na-K-ATPase
!---------------------------------------------------------------------72

  AffNa = 0.2d0*(1.0d0+C(2,2)/8.33d0)
  actNa = C(1,2)/(C(1,2)+AffNa)

  AffK5 = 0.1d0*(1.0d0+C(1,5)/18.5d0)
  AffNH5 =  AffK5
  actK5 = (C(2,5)+C(11,5))/(C(2,5)+C(11,5)+AffK5)
 
  AffK6 = 0.1d0*(1.0d0+C(1,6)/18.5d0)
  AffNH6 =  AffK6
  actK6 = (C(2,6)+C(11,6))/(C(2,6)+C(11,6)+AffK6)

  ro5 = (C(11,5)/AffNH5)/(C(2,5)/AffK5)
  ro6 = (C(11,6)/AffNH6)/(C(2,6)/AffK6)

  dJactNa5 = pt%area(2,5)*pt%ATPNaK(2,5)*(actNa**3.0d0)*(actK5**2.0d0)
  dJactNa6 = pt%area(2,6)*pt%ATPNaK(2,6)*(actNa**3.0d0)*(actK6**2.0d0)

  dJactK5 = -2.0d0/3.0d0*dJactNa5/(1.d0+ro5)
  dJactK6 = -2.0d0/3.0d0*dJactNa6/(1.d0+ro6)

  Jsol(1,2,5) = Jsol(1,2,5)+dJactNa5
  Jsol(1,2,6) = Jsol(1,2,6)+dJactNa6

  Jsol(2,2,5) = Jsol(2,2,5)+dJactK5
  Jsol(2,2,6) = Jsol(2,2,6)+dJactK6

  Jsol(11,2,5) = Jsol(11,2,5)+dJactK5*ro5
  Jsol(11,2,6) = Jsol(11,2,6)+dJactK6*ro6

  fluxNaKsod = (dJactNa5+dJactNa6)*convert
  fluxNaKpot = (dJactK5+dJactK6)*convert
  fluxNaKamm = (dJactK5*ro5+dJactK6*ro6)*convert

!---------------------------------------------------------------------72
!	 H-ATPase
!	 See Strieter & Weinstein paper for signs (AJP 263, 1992)
!---------------------------------------------------------------------72

  DactH=1.0d0+dexp(steepATPH*(delmu(12,1,2)-dmuATPHPT))
  dJactH=-pt%area(1,2)*pt%ATPH(1,2)/DactH

  Jsol(12,1,2) = Jsol(12,1,2)+dJactH

  fluxHATPase = dJactH*convert

!---------------------------------------------------------------------72
!   Putative basolateral PMCA pump
!---------------------------------------------------------------------72

  AffCa = 75.6d-6 !Tsukamoto et al., Biochim Biophys Acta 1992
  dJPMCA5 = pt%area(2,5)*pt%PMCA*C(16,2)/(C(16,2)+AffCa)
  dJPMCA6 = pt%area(2,6)*pt%PMCA*C(16,2)/(C(16,2)+AffCa)

  Jsol(16,2,5) = Jsol(16,2,5)+dJPMCA5
  Jsol(16,2,6) = Jsol(16,2,6)+dJPMCA6

  fluxPMCA = (dJPMCA5+dJPMCA6)*convert

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!
!     ADD TORQUE-MODULATED EFFECTS ON TRANSCELLULAR TRANSPORTER DENSITY
!
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!		Torque scaling parameter (Eq. 39 in 2007 PT model)
!       Scale transcellular water and solute fluxes

  Do I = 1, NSPT
     Jsol(I,1,2) = Scaletorq*Jsol(I,1,2)
     Jsol(I,2,5) = Scaletorq*Jsol(I,2,5)
     Jsol(I,2,6) = Scaletorq*Jsol(I,2,6)
  End do

  Jvol(1,2) = Scaletorq*Jvol(1,2)
  Jvol(2,5) = Scaletorq*Jvol(2,5)
  Jvol(2,6) = Scaletorq*Jvol(2,6)

!---------------------------------------------------------------------72
!
!     STORE LOCAL FLUX VALUES FOR LATER INTEGRATION
!
!---------------------------------------------------------------------72

  cv=href*Cref*PI*pt%diam*60/10*1.0d9 !for dimensional solute fluxes
  cvw=Pfref*Cref*Vwbar*PI*pt%diam*60/10*1.0d6 !for dimensional water flux

  pt%FNatrans = Jsol(1,1,2)*cv  ! FNatrans
  pt%FNapara = Jsol(1,1,5)*cv  ! FNapara
  pt%FNaK = fluxNaKsod  ! FNaK, already converted to pmol/min/mm
  pt%FHase = fluxHATPase ! FHase
  pt%FGluPara = Jsol(15,1,5)*cv
  pt%FGluSGLT1 = fluxglusglt1
  pt%FGluSGLT2 = fluxglusglt2

  pt%FKtrans = Jsol(2,1,2)*cv
  pt%FKpara = Jsol(2,1,5)*cv

  FCaTransPT(LzPT+1) = Jsol(16,1,2)*cv
  FCaParaPT(LzPT+1) = Jsol(16,1,5)*cv

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

500 nlocal = 1

  return
end Subroutine qfluxPT

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
