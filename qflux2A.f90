!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!      This subroutine computes the fluxes in the mTAL

Subroutine qflux2A(x,Jvol,Jsol,Cext,EPext,mtal)

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

  type (membrane) :: mtal
  double precision x(NDA)
  double precision Jvol(NC,NC), Jsol(NS,NC,NC)
  double precision Cext(NS),EPext
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

  mtal%area(5,6)=mtal%sbasEinit*max(Vol(5)/mtal%volEinit,1.0d0)
  mtal%area(6,5)=mtal%area(5,6)

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

  call compute_water_fluxes (C,PM,0.0d0,Vol,mtal%volLuminit,mtal%volEinit,  &
       mtal%volPinit,CPimprefA,mtal%volAinit,CAimprefA,mtal%volBinit, & 
	   CBimprefA,mtal%area,mtal%sig,mtal%dLPV,complA,Jvol)


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     COMPUTE SOLUTE FLUXES
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	 Non-dimensional solute fluxes are first converted to dimensional fluxes 
!      (in mmol/s/cm2 epith) by multiplying by href*Cref
!	 Then they are converted to units of pmol/min/mm tubule 
  convert=href*Cref*PI*DiamA*60/10*1.0d9 !Conversion factor
  eps = 1.0d-4 !For exponential (Peclet) terms

!---------------------------------------------------------------------72
!     ELECTRO-CONVECTIVE-DIFFUSIVE FLUXES
!---------------------------------------------------------------------72

  call compute_ecd_fluxes (C,EP,mtal%area,mtal%sig,mtal%h,Jvol,Jsol,delmu)

 fluxkap = Jsol(2,1,2)*convert
 fluxkbas = (Jsol(2,2,5)+Jsol(2,2,6))*convert
!---------------------------------------------------------------------72	
!---------------------------------------------------------------------72
!	 FLUXES ACROSS COTRANSPORTERS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	 Basolateral Na2HPO4 cotransporter 
!---------------------------------------------------------------------72

  sumJES=0.0d0
  Do L = 5,6
     dJNaP=mtal%area(2,L)*mtal%dLA(1,7,2,L)*(2*delmu(1,2,L)+delmu(7,2,L)) 
     Jsol(1,2,L)=Jsol(1,2,L)+2*dJNaP
     Jsol(7,2,L)=Jsol(7,2,L)+dJNaP
     sumJES=sumJES+2*dJNaP
  End do
  fluxNaPatPES=sumJES*convert

!---------------------------------------------------------------------72
!	 Basolateral Na(1)-HCO3(3) cotransporter 
!---------------------------------------------------------------------72

  sumJES = 0.d0
  Do L = 5,6
     dJNaBic=mtal%area(2,L)*mtal%dLA(1,4,2,L)*(delmu(1,2,L)+3*delmu(4,2,L))
     Jsol(1,2,L)=Jsol(1,2,L)+dJNaBic
     Jsol(4,2,L)=Jsol(4,2,L)+3*dJNaBic
     sumJES=sumJES+dJNaBic
     
  End do
  fluxNaBicPES = sumJES*convert
 
!---------------------------------------------------------------------72
!	 NKCC2 cotransporter on MP interface
!---------------------------------------------------------------------72

!	    F isoform
!---------------------------------------------------------------------72

  call compute_nkcc2_flux ( C,mtal%area(1,2),mtal%xNKCC2F,bn2F,bk2F,bc2F,bm2F,  &
       popnkccF,pnkccpF,pnmccpF,poppnkccF,pnkccppF,pnmccppF, &
       dJnNKCC2F,dJkNKCC2F,dJcNKCC2F,dJmNKCC2F )
  
  Jsol(1,1,2) = Jsol(1,1,2)+dJnNKCC2F
  Jsol(2,1,2) = Jsol(2,1,2)+dJkNKCC2F
  Jsol(3,1,2) = Jsol(3,1,2)+dJcNKCC2F
  Jsol(11,1,2) = Jsol(11,1,2)+dJmNKCC2F

  fluxnNKCC2F = dJnNKCC2F*convert
  fluxkNKCC2F = dJkNKCC2F*convert
  fluxmNKCC2F = dJmNKCC2F*convert
  fluxcNKCC2F = dJcNKCC2F*convert


!	    A isoform
!---------------------------------------------------------------------72

  call compute_nkcc2_flux ( C,mtal%area(1,2),mtal%xNKCC2A,bn2A,bk2A,bc2A,bm2A, &
       popnkccA,pnkccpA,pnmccpA,poppnkccA,pnkccppA,pnmccppA, &
       dJnNKCC2A,dJkNKCC2A,dJcNKCC2A,dJmNKCC2A )
  
  Jsol(1,1,2) = Jsol(1,1,2)+dJnNKCC2A
  Jsol(2,1,2) = Jsol(2,1,2)+dJkNKCC2A
  Jsol(3,1,2) = Jsol(3,1,2)+dJcNKCC2A
  Jsol(11,1,2) = Jsol(11,1,2)+dJmNKCC2A
  
  fluxnNKCC2A = dJnNKCC2A*convert
  fluxkNKCC2A = dJkNKCC2A*convert
  fluxmNKCC2A = dJmNKCC2A*convert
  fluxcNKCC2A = dJcNKCC2A*convert

  fluxnNKCC2 = fluxnNKCC2A+fluxnNKCC2F
  fluxkNKCC2 = fluxkNKCC2A+fluxkNKCC2F
  fluxmNKCC2 = fluxmNKCC2A+fluxmNKCC2F
  fluxcNKCC2 = fluxcNKCC2A+fluxcNKCC2F

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

  call compute_nhe3_fluxes (C,mtal%area(1,2),mtal%xNHE3,dJNHEsod,dJNHEprot,dJNHEamm)

  Jsol(1,1,2)=Jsol(1,1,2)+dJNHEsod
  Jsol(12,1,2)=Jsol(12,1,2)+dJNHEprot
  Jsol(11,1,2)=Jsol(11,1,2)+dJNHEamm
  fluxNHEsod=dJNHEsod*convert
  fluxNHEprot=dJNHEprot*convert
  fluxNHEamm=dJNHEamm*convert


!---------------------------------------------------------------------72
!	 Basolateral NaH exchangers 
!---------------------------------------------------------------------72


  Do K = 2,2
     sumJES=0.0d0
     Do L = 5,6
        dJNaH=mtal%area(K,L)*mtal%dLA(1,12,K,L)*(delmu(1,K,L)-delmu(12,K,L))
        Jsol(1,K,L) = Jsol(1,K,L)+dJNaH
        Jsol(12,K,L)=Jsol(12,K,L)-dJNaH
        sumJES=sumJES+dJNaH
     End do
     fluxNaHPES = sumJES*convert
  End Do


!---------------------------------------------------------------------72
!	 Cl/HCO3 exchanger at PE,PS interfaces
!---------------------------------------------------------------------72

  sumJES=0.0d0
  Do L = 5,6
     dJClHCO3=mtal%area(2,L)*mtal%dLA(3,4,2,L)*(delmu(3,2,L)-delmu(4,2,L))
     Jsol(3,2,L) = Jsol(3,2,L)+dJClHCO3
     Jsol(4,2,L) = Jsol(4,2,L)-dJClHCO3
     sumJES=sumJES+dJClHCO3
  End do
  fluxClHCO3exPES = sumJES*convert
  
  
!---------------------------------------------------------------------72
!	 KCC cotransporter on PE & PS interfaces
!---------------------------------------------------------------------72

  call compute_kcc_fluxes (C,mtal%area(2,5),mtal%area(2,6), &
       mtal%xKCC4,dJk5,dJc5,dJm5, dJk6,dJc6,dJm6)

  Jsol(2,2,5) = Jsol(2,2,5)+dJk5
  Jsol(3,2,5) = Jsol(3,2,5)+dJc5
  Jsol(11,2,5) = Jsol(11,2,5)+dJm5

  Jsol(2,2,6) = Jsol(2,2,6)+dJk6
  Jsol(3,2,6) = Jsol(3,2,6)+dJc6
  Jsol(11,2,6) = Jsol(11,2,6)+dJm6
  
  fluxkKCC = (dJk5 + dJk6)*convert
  fluxcKCC = (dJc5 + dJc6)*convert
  fluxmKCC = (dJm5 + dJm6)*convert


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
  AffK = 0.1d0*(1.0d0+C(1,6)/18.5d0)
  AffNH4 =  AffK
  actK5 = C(2,5)/(C(2,5)+AffK)
  actK6 = C(2,6)/(C(2,6)+AffK)
  
  dJact5 = mtal%area(2,5)*mtal%ATPNaK(2,5)*(actNa**3.0d0)*(actK5**2.0d0)
  dJact6 = mtal%area(2,6)*mtal%ATPNaK(2,6)*(actNa**3.0d0)*(actK6**2.0d0)
  
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
!     STORE LOCAL FLUX VALUES FOR LATER INTEGRATION
!---------------------------------------------------------------------72

  cv=href*Cref*PI*DiamA*60/10*1.0d9 !for dimensional solute fluxes
  cvw=Pfref*Cref*Vwbar*PI*DiamA*60/10*1.0d6 !for dimensional water flux

  mtal%FNatrans = Jsol(1,1,2)*cv
  mtal%FNapara  = Jsol(1,1,5)*cv
  mtal%FNaK     = fluxNaKPESsod
  mtal%FHase    = 0.0

  mtal%FKtrans = Jsol(2,1,2)*cv
  mtal%FKpara  = Jsol(2,1,5)*cv

  return
end Subroutine qflux2A

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
