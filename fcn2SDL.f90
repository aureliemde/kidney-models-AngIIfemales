!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	 This subroutine yields the non-linear equations to be solved at any
!	 point below the SDL inlet, in order to determine the *luminal*
!        concentrations, volumes, and EP
!    It is used to compute values below the entry to the SDL.
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
subroutine fcn2SDL(n,x,fvec,iflag,numpar,pars,sdl,id)
  
  include 'values.h'
  include 'global.h'
  include 'defs.h'

  type (membrane) :: sdl
  integer n,iflag,numpar,id
  double precision x(n),fvec(n),S(NDA),pars(numpar)
  double precision fsola(n),fsolb(n)
  double precision Ca(NS,NC),Vola(NC) !,EPa(NC),pha(NC)
  double precision Cb(NS,NC),Volb(NC),EPb(NC),phb(NC)
  double precision PMa,PMb,VolEinit,VolPinit,Vol0
  double precision Jva(NC,NC), Jsa(NS,NC,NC),Jcasra
  double precision Jvb(NC,NC), Jsb(NS,NC,NC),Jcasrb
  double precision CaBT(NC), dkh(NC), dkd(NC), Diam
  double precision dimL, CPimpref

!	 Ca(K) denotes the known concentrations at Lz (same with pH, Vol, EP)
!	 y is the associated vector, of known variables

!	 Cb(K) denotes the concentrations at Lz+1	(same with pH, Vol, EP)
!	 x is the associated vector, of unknown variables

!---------------------------------------------------------------------72
!	 Assign concentrations, volumes, and potentials 
!	 in P, A, B, E, and the lumen M at Lz
!---------------------------------------------------------------------72

  Do I = 1,NS
     Ca(I,1)=pars(I)
     Ca(I,2)=pars(I+NS)
     Ca(I,5)=pars(I+2*NS)
  End do
  
  Vola(1)=pars(1+3*NS)
  Vola(2)=pars(2+3*NS)
  Vola(5)=pars(3+3*NS)
  
  PMa=pars(4+3*NS)

!---------------------------------------------------------------------72
!	 Assign concentrations and potentials in peritubular solution
!---------------------------------------------------------------------72
  Do J = 1,NS
     Ca(J,6) = pars(J+4+3*NS)
  End Do
  Cb(:,6) = sdl%conc(:,6)
  EPb(6) = sdl%ep(6)
  
!---------------------------------------------------------------------72
!	 Extract other parameters
!---------------------------------------------------------------------72
!  VolEinit = pars(5+4*NS)
!  VolPinit = pars(6+4*NS)
  dimL = pars(5+4*NS)
  CPimpref = pars(6+4*NS)
  Diam   = pars(7+4*NS)
  Vol0   = pars(8+4*NS)

!---------------------------------------------------------------------72
!	 Assign concentrations, volumes, and potentials 
!	 in the lumen at Lz+1; don't care about other compartments
!---------------------------------------------------------------------72

  Do I = 1,NS2
     Cb(I,1)=x(I)
  End do
  Volb(1) = x(1+NS2) 
  PMb = x(2+NS2)
  phb(1)=-dlog(Cb(12,1)/1.0d3)/dlog(10.0d0)

  !  SDL is a little water permeable, but completely solute impermeable
  Call compute_sdl_water_fluxes (Cb,PMb,Volb,Vol0,sdl%area,sdl%sig,sdl%dLPV,Jvb)
  
!---------------------------------------------------------------------72
!	Initialize  source terms
!---------------------------------------------------------------------72

  Do K = 1, n
     S(K) = 0.0d0
     fvec(K) = 0.0d0
  End Do	 

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE SOURCE TERMS FOR FLOW OF SOLUTE I IN THE LUMEN
!	The non-dimensional volume flux needs to be multiplied by 
!	 (Vref/href)
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  Bm = PI*Diam
  Am = PI*(Diam**2)/4.0d0
  
  Do I = 1, NS2
     sumJsb=0!  assume solute impermeable 
     fsola(I)=Vola(1)*Ca(I,1)*Vref/href
     fsolb(I)=Volb(1)*Cb(I,1)*Vref/href
     S(I)=fsolb(I)-fsola(I)+Bm*dimL*sumJsb/NZ 
  End do

 ! For non-reacting species
     fvec(1) = S(1)
	 fvec(2) = S(2)
	 fvec(3) = S(3)
	 fvec(9) = S(9)
	 fvec(15) = S(15)
	 fvec(16) = S(16)

!---------------------------------------------------------------------72
!	 For CO2/HCO3/H2CO3
!---------------------------------------------------------------------72
  dkhuncat = 0.145d0
  dkduncat = 49.60d0

  fvec(4) = S(4)+S(5)+S(6)
  fvec(5) = phb(1)-pKHCO3-dlog(Cb(4,1)/Cb(5,1))/dlog(10.0d0)	 
  fklum=(dkhuncat*Cb(6,1)-dkduncat*Cb(5,1))
  fvec(6) = S(6) + Am*dimL*fklum/NZ/href


!---------------------------------------------------------------------72
!	 For HPO4(2-)/H2PO4(-)
!---------------------------------------------------------------------72

  fvec(7) = S(7)+S(8)
  fvec(8) = phb(1)-pKHPO4-dlog(Cb(7,1)/Cb(8,1))/dlog(10.0d0)

!---------------------------------------------------------------------72
!	 For NH3/NH4
!---------------------------------------------------------------------72
  
   fvec(10) = S(10)+S(11)
   fvec(11) = phb(1)-pKNH3-dlog(Cb(10,1)/Cb(11,1))/dlog(10.0d0)

!---------------------------------------------------------------------72
!	 For HCO2-/H2CO2
!---------------------------------------------------------------------72

  fvec(13) = S(13)+S(14)
  fvec(14) = phb(1)-pKHCO2-dlog(abs(Cb(13,1)/Cb(14,1)))/dlog(10.0d0)

!---------------------------------------------------------------------72
!	 For pH
!---------------------------------------------------------------------72

  fvec(12) = S(12)+S(11)-S(4)-S(7)-S(13)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE SOURCE TERMS FOR VOLUME FLOW IN THE LUMEN
!	The non-dimensional volume flux needs to be multiplied by 
!	Vref/(Pfref.Vwbar.Cref)
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  fvmult = (Pfref*Vwbar*Cref)
  sumJvb=Jvb(1,6)*fvmult
  fvola=Vola(1)*Vref
  fvolb=Volb(1)*Vref
  fvec(1+NS2)=fvolb-fvola+Bm*dimL*sumJvb/NZ

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!		Equation for pressure
!		Flow must be multiplied by Vref
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  ratio=8.0d0*visc/(PI*(0.5d0*Diam)**4)
  fvec(2+NS2)=PMb-PMa+ratio*Volb(1)*Vref*dimL/NZ


  return
end subroutine fcn2SDL



!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE SDL WATER FLUXES
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	 The hydraulic and oncotic pressures are made non-dimensional 
!	 by dividing by RT*Cref

subroutine compute_sdl_water_fluxes (C,PM,Vol,Vol0,area,sig,dLPV,Jvol)

  include 'values.h'
  include 'global.h'

  ! input variables
  double precision C(NS2,NC),area(NC,NC),sig(NS,NC,NC),dLPV(NC,NC)
  double precision Vol(NC),Vol0,PM

  ! outputs
  double precision Jvol(NC,NC)

  ! local variables
  double precision PRES(NC),ONC(NC)
  integer J
  
  !	 The hydraulic and oncotic pressures are made non-dimensional 
  !	 by dividing by RT*Cref

  PRES(1)=PM/(RTosm*Cref)
  PRES(6)=PbloodPT/(RTosm*Cref)
  ONC(1)=LumImperm*PTinitVol/Vol(1)
  ONC(6)=BathImperm

  OSM=0.0d0
  Do J = 1, NS2
     OSM=OSM+sig(J,1,6)*(C(J,1)-C(J,6))
  End Do
  Jvol(1,6)=(area(1,2)*dLPV(1,2)+area(1,5)*dLPV(1,5))*(PRES(1)-PRES(6)-ONC(1)+ONC(6)-OSM)
  Jvol(6,1)=-Jvol(1,6)



  return
end subroutine compute_sdl_water_fluxes


