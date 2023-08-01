!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
 
!  *************************** FLUXES ***************************
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!      This subroutine computes the fluxes in the IMCD

Subroutine qflux2IMC(x,Jvol,Jsol,Cext,EPext,imcd)

! INPUT ARGUMENTS:
!	 x: vector of 4*NS concentrations, 4 volumes, and 4 potentials
! 
! OUTPUT ARGUMENTS:
!     Jvol:    volume fluxes
!     Jsol:    solute fluxes

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  include 'values.h'
  include 'global.h'
  include 'defs.h'

  type (membrane) :: imcd
  double precision x(7+3*NS2)
  double precision Jvol(NC,NC), Jsol(NS,NC,NC)
  double precision :: Cext(NS), EPext
  double precision ONC(NC), PRES(NC)
  double precision C(NS,NC),dmu(NS,NC),ph(NC),Vol(NC),EP(NC)
  double precision delmu(NS,NC,NC)
  double precision hkconc(4),Amat(Natp,Natp)
  double precision theta(NC),Slum(NC),Slat(NC),Sbas(NC)

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

  Do J = 1, NS
     C(J,6) = Cext(J)
  End Do
  EP(6) = EPext 
	
  Do I = 1,NS2
     C(I,1)=x(1+3*(I-1))
     C(I,2)=x(2+3*(I-1))
     C(I,5)=x(3+3*(I-1))
  End do

  Vol(1) = x(1+3*NS2) 
  Vol(2) = x(2+3*NS2) 
  Vol(5) = x(3+3*NS2) 
  EP(1) = x(4+3*NS2) 
  EP(2) = x(5+3*NS2) 
  EP(5) = x(6+3*NS2) 
  PM = x(7+3*NS2)

!	 Assign dummy values to compartments 3 and 4
  Do I = 1,NS
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

  imcd%area(5,6)=imcd%sbasEinit*max(Vol(5)/imcd%volEinit,1.0d0)
  imcd%area(6,5)=imcd%area(5,6)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	INITIALIZE ALL FLUXES 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  Do K = 1, NC
     Do L = 1, NC
        Jvol(K,L) = 0.0d0
        Do I = 1, NS
           Jsol(I,K,L) = 0.0d0
        End Do
     End Do
  End Do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE WATER FLUXES
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  call compute_water_fluxes (C,PM,0.0d0,Vol,imcd%volLuminit,imcd%volEinit,  &
       imcd%volPinit,CPimprefIMC,imcd%volAinit,CAimprefIMC,imcd%volBinit,   &
	   CBimprefIMC,imcd%area,imcd%sig,imcd%dLPV,complIMC,Jvol)


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     COMPUTE SOLUTE FLUXES
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	 Non-dimensional solute fluxes are first converted to dimensional fluxes 
!      (in mmol/s/cm2 epith)by multiplying by href*Cref
!	 Then they are converted to units of pmol/min/mm tubule 
  convert=href*Cref*PI*DiamIMC*60/10*1.0d9 !Conversion factor


!---------------------------------------------------------------------72
!     ELECTRO-CONVECTIVE-DIFFUSIVE FLUXES
!---------------------------------------------------------------------72

!    dependence of ENaC, ROMK, and apical paracellular Cl permeability on pH

   facphMP = 1.0*(0.1 + 2.0d0/ (1+dexp(-6.0d0*(pH(2)-7.50d0))))
   facphTJ = 2.0/(1.0 + dexp(10.0*(ph(5)-7.32d0)))

   facNaMP=(30.d0/(30.d0+C(1,1)))*(50.d0/(50.d0+C(1,2)))
!   imcd%h(1,1,2)=imcd%hNaMP*facNaMP*facphMP
   imcd%h(1,1,2)=hENaC_IMC*facNaMP*facphMP

   imcd%h(2,1,2)=hROMK_IMC*facphMP

   imcd%h(3,1,5)=hCltj_IMC*facphTJ


  call compute_ecd_fluxes (C,EP,imcd%area,imcd%sig,imcd%h,Jvol,Jsol,delmu)

!	    dimensional fluxes through channels, for print out

  fluxNachMP=Jsol(1,1,2)*convert
  fluxNachPES=(Jsol(1,2,5)+Jsol(1,2,6))*convert
  
  fluxKchMP=Jsol(2,1,2)*convert
  fluxKchPES=(Jsol(2,2,5)+Jsol(2,2,6))*convert
  
  fluxClchMP=Jsol(3,1,2)*convert
  fluxClchPES=(Jsol(3,2,5)+Jsol(3,2,6))*convert
  fluxBichPES=(Jsol(4,2,5)+Jsol(4,2,6))*convert
  
  fluxH2CO3MP=Jsol(5,1,2)*convert
  fluxCO2MP=Jsol(6,1,2)*convert
  
  fluxH2CO3PES=(Jsol(5,2,5)+Jsol(5,2,6))*convert
  fluxCO2PES=(Jsol(6,2,5)+Jsol(6,2,6))*convert
  
  fluxHP2mchPES=(Jsol(7,2,5)+Jsol(7,2,6))*convert
  fluxHPmPES=(Jsol(8,2,5)+Jsol(8,2,6))*convert
  
  fluxNH3MP=Jsol(10,1,2)*convert
  fluxNH3PES=(Jsol(10,2,5)+Jsol(10,2,6))*convert

!---------------------------------------------------------------------72	
!---------------------------------------------------------------------72
!	 FLUXES ACROSS COTRANSPORTERS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	 Apical NaCl cotransporter 
!---------------------------------------------------------------------72

  dJNaCl=imcd%area(1,2)*imcd%dLA(1,3,1,2)*(delmu(1,1,2)+delmu(3,1,2)) 
  Jsol(1,1,2)=Jsol(1,1,2)+dJNaCl
  Jsol(3,1,2)=Jsol(3,1,2)+dJNaCl
  fluxNaClMP=dJNaCl*convert
  
!  fNCC(LzIMC) = dJNaCl

!---------------------------------------------------------------------72
!	 Basolateral Na2HPO4 cotransporter 
!---------------------------------------------------------------------72

!		Factor 1/2 to be consistent with AMW model results
  sumJES=0.0d0
  Do L = 5,6
     dJNaP=imcd%area(2,L)*imcd%dLA(1,7,2,L)/2.d0*(2*delmu(1,2,L)+delmu(7,2,L)) 
     Jsol(1,2,L)=Jsol(1,2,L)+2*dJNaP
     Jsol(7,2,L)=Jsol(7,2,L)+dJNaP
     sumJES=sumJES+2*dJNaP
  End do
  fluxNaPatPES=sumJES*convert
  
!---------------------------------------------------------------------72
!	 Basolateral KCl cotransporter
!---------------------------------------------------------------------72

  sumJES = 0.d0
  Do L = 5,6
     dJKCl=imcd%area(2,L)*imcd%dLA(2,3,2,L)  &
          *(delmu(2,2,L)+delmu(3,2,L))
     Jsol(2,2,L)=Jsol(2,2,L)+dJKCl
     Jsol(3,2,L)=Jsol(3,2,L)+dJKCl
     sumJES=sumJES+dJKCl
  End do
  fluxKClPES = sumJES*convert

!---------------------------------------------------------------------72
!	 Basolateral NKCl2 cotransporter
!---------------------------------------------------------------------72

  sumJES = 0.d0
  Do L = 5,6
     dJNKCl2=imcd%area(2,L)*imcd%dLA(1,2,2,L)  &
          *(delmu(1,2,L)+delmu(2,2,L)+2*delmu(3,2,L))
     Jsol(1,2,L)=Jsol(1,2,L)+dJNKCl2
     Jsol(2,2,L)=Jsol(2,2,L)+dJNKCl2
     Jsol(3,2,L)=Jsol(3,2,L)+2*dJNKCl2
     sumJES=sumJES+dJNKCl2
  End do
  fluxNKCl2PES = sumJES*convert

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	 FLUXES ACROSS EXCHANGERS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	 Basolateral NaH exchangers 
!---------------------------------------------------------------------72

  sumJES=0.0d0
  Do L = 5,6
     dJNaH=imcd%area(2,L)*imcd%dLA(1,12,2,L)  &
          *(delmu(1,2,L)-delmu(12,2,L))
     Jsol(1,2,L) = Jsol(1,2,L)+dJNaH
     Jsol(12,2,L)=Jsol(12,2,L)-dJNaH
     sumJES=sumJES+dJNaH
  End do
  fluxNaHPES = sumJES*convert

!---------------------------------------------------------------------72
!	 Cl/HCO3 exchanger at PE,PS interfaces
!---------------------------------------------------------------------72

  sumJES=0.0d0
  Do L = 5,6
     dJClHCO3=imcd%area(2,L)*imcd%dLA(3,4,2,L)*  &
          (delmu(3,2,L)-delmu(4,2,L))
     Jsol(3,2,L) = Jsol(3,2,L)+dJClHCO3
     Jsol(4,2,L) = Jsol(4,2,L)-dJClHCO3
     sumJES=sumJES+dJClHCO3
  End do
  fluxClHCO3exPES = sumJES*convert
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	ACTIVE TRANSPORT VIA ATPases
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	Basolateral Na-K-ATPase
!---------------------------------------------------------------------72

  AffNa = 0.2d0*(1.0d0+C(2,2)/8.33d0)
  actNa = C(1,2)/(C(1,2)+AffNa)
  AffK = 0.1d0*(1.0d0+C(1,6)/18.5d0)
  AffNH4 =  5.0d0*AffK
  actK5 = C(2,5)/(C(2,5)+AffK)
  actK6 = C(2,6)/(C(2,6)+AffK)
  
  dJact5 = imcd%area(2,5)*imcd%ATPNaK(2,5)*(actNa**3.0d0) *(actK5**2.0d0)
  dJact6 = imcd%area(2,6)*imcd%ATPNaK(2,6)*(actNa**3.0d0) *(actK6**2.0d0)
  
  ro5 = (C(11,5)/AffNH4)/(C(2,5)/AffK)
  ro6 = (C(11,6)/AffNH4)/(C(2,6)/AffK)
  
  Jsol(1,2,5) = Jsol(1,2,5)+dJact5
  Jsol(1,2,6) = Jsol(1,2,6)+dJact6
  
  Jsol(2,2,5) = Jsol(2,2,5)-2.0d0/3.0d0*dJact5/(1.d0+ro5)
  Jsol(2,2,6) = Jsol(2,2,6)-2.0d0/3.0d0*dJact6/(1.d0+ro6)
  
  Jsol(11,2,5) = Jsol(11,2,5)-2.0d0/3.0d0*dJact5*ro5/(1.d0+ro5)
  Jsol(11,2,6) = Jsol(11,2,6)-2.0d0/3.0d0*dJact6*ro6/(1.d0+ro6)
  
  fluxNaKPESsod = (dJact5+dJact6)*convert
  
  fluxNaKPESpot = -2.d0/3.0d0*(dJact5/(1.d0+ro5)+dJact6/(1.d0+ro6))*convert
  fluxNaKPESamm = -2.d0/3.0d0*(dJact5*ro5/(1.d0+ro5)+dJact6*ro6/(1.d0+ro6))*convert

!---------------------------------------------------------------------72
!	 Luminal H-K-ATPase	
!---------------------------------------------------------------------72

  hkconc(1) = C(2,2) ![K]in = K in P
  hkconc(2) = C(2,1) ![K]out = K in M
  hkconc(3) = C(12,2) ![H]in = H in P
  hkconc(4) = C(12,1) ![H]out = H in M
  
  call fatpase(Natp,hkconc,Amat)

!     invert Jacobian matrix
  call dgetrf( Natp, Natp, Amat, Natp, ipiv, info )
  call dgetri( Natp, Amat, Natp, ipiv, work, lwork, info )

!     determine needed variables for calculating H efflux = K influx
!	  Need c(7) = H2-E1-P and c(8) = H2-E2-P

  c7 = Amat(7,1)
  c8 = Amat(8,1)
  dkf5 = 4.0d1
  dkb5 = 2.0d2
  hefflux = imcd%area(1,2)*imcd%ATPHK(1,2)*(dkf5*c7-dkb5*c8)
  Jsol(2,1,2) = Jsol(2,1,2)+hefflux
  Jsol(12,1,2) = Jsol(12,1,2)-hefflux
  fluxHKATPaseMP = hefflux*convert

!---------------------------------------------------------------------72
!     STORE LOCAL FLUX VALUES FOR LATER INTEGRATION
!---------------------------------------------------------------------72

  cv=href*Cref*PI*DiamIMC*60/10*1.0d9 !for dimensional solute fluxes
  cvw=Pfref*Cref*Vwbar*PI*DiamIMC*60/10*1.0d6 !for dimensional water flux

  imcd%FNatrans = Jsol(1,1,2)*cv*imcd%coalesce
  imcd%FNapara = Jsol(1,1,5)*cv*imcd%coalesce
  imcd%FNaK = fluxNaKPESsod*cv*imcd%coalesce
  imcd%FHase = 0

  imcd%FKtrans = Jsol(2,1,2)*cv*imcd%coalesce
  imcd%FKpara = Jsol(2,1,5)*cv*imcd%coalesce


  return
end Subroutine qflux2IMC

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
