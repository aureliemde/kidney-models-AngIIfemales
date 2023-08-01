!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
 
!  *************************** FLUXES ***************************
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!      This subroutine computes the fluxes in the OMCD

Subroutine qflux2OMC(x,Jvol,Jsol,omcd)

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
  
  type (membrane) :: omcd
  double precision x(NDC)
  double precision Jvol(NC,NC), Jsol(NS,NC,NC)
  double precision ONC(NC), PRES(NC)
  double precision C(NS,NC),dmu(NS,NC),ph(NC),Vol(NC),EP(NC)
  double precision delmu(NS,NC,NC), delmu2(NS,NC,NC)
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
     C(J,6) = omcd%conc(J,6)
  End Do
  EP(6) = omcd%ep(6)
	
  Do I = 1,NS2
     C(I,1)=x(1+5*(I-1))
     C(I,2)=x(2+5*(I-1))
     C(I,3)=x(3+5*(I-1))
     C(I,4)=x(4+5*(I-1))
     C(I,5)=x(5+5*(I-1))
  End do
   
  Do K = 1,NC-1
     ph(K)=-dlog(C(12,K)/1.0d3)/dlog(10.0d0)
     Vol(K)=x(K+5*NS2)
     EP(K)=x(K+5+5*NS2)
  End do

  PM=x(11+5*NS2) !Hydraulic pressure in lumen

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!		DETERMINE THE NEW SURFACE AREAS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  omcd%area(5,6)=omcd%sbasEinit*max(Vol(5)/omcd%volEinit,1.0d0)
  omcd%area(6,5)=omcd%area(5,6)
	
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


  call compute_water_fluxes (C,PM,0.0d0,Vol,omcd%volLuminit,omcd%volEinit,  &
       omcd%volPinit,CPimprefOMC,omcd%volAinit,CAimprefOMC,omcd%volBinit,   &
	   CBimprefOMC,omcd%area,omcd%sig,omcd%dLPV,complOMC,Jvol)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     COMPUTE SOLUTE FLUXES
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	 Non-dimensional solute fluxes are first converted to dimensional fluxes 
!      (in mmol/s/cm2 epith)by multiplying by href*Cref
!	 Then they are converted to units of pmol/min/mm tubule 
  convert=href*Cref*PI*DiamOMC*60/10*1.0d9 !Conversion factor

!---------------------------------------------------------------------72
!     ELECTRO-CONVECTIVE-DIFFUSIVE FLUXES
!---------------------------------------------------------------------72

!    dependence of ENaC, ROMK, and apical paracellular Cl permeability on pH

   facphMP = 1.0*(0.1 + 2.0d0/ (1+dexp(-6.0d0*(pH(2)-7.50d0))))
   facphTJ = 2.0/(1.0 + dexp(10.0*(ph(5)-7.32d0)))

   facNaMP=(30.d0/(30.d0+C(1,1)))*(50.d0/(50.d0+C(1,2)))
!   omcd%h(1,1,2)=omcd%hNaMP*facNaMP*facphMP
   omcd%h(1,1,2)=hENaC_OMC*facNaMP*facphMP

   omcd%h(2,1,2)=hROMK_OMC*facphMP

   omcd%h(3,1,5)=hCltj_OMC*facphTJ
  
  !	 define electrochemical potential of each solute in each compartment

  call compute_ecd_fluxes (C,EP,omcd%area,omcd%sig,omcd%h,Jvol,Jsol,delmu)

  XI = zval(1)*F*EPref/RT*(EP(1)-EP(2))
  dint = dexp(-XI)
  if (dabs(1.0d0-dint).lt. 1.d-6) then
            JNa=omcd%area(1,2)*omcd%h(1,1,2)*(C(1,1)-C(1,2))
  else
            JNa=omcd%area(1,2)*omcd%h(1,1,2)*XI &
                    *(C(1,1)-C(1,2)*dint) /(1.0d0-dint)
  end if

  !	    dimensional fluxes through channels, for print out 
  
  fluxENaC=Jsol(1,1,2)*convert
  fluxROMK=Jsol(2,1,2)*convert

  fluxKchPES=(Jsol(2,2,5)+Jsol(2,2,6))*convert
  fluxKchAES=(Jsol(2,3,5)+Jsol(2,3,6))*convert
  fluxClchPES=(Jsol(3,2,5)+Jsol(3,2,6))*convert
  fluxClchAES=(Jsol(3,3,5)+Jsol(3,3,6))*convert
  fluxBichPES=(Jsol(4,2,5)+Jsol(4,2,6))*convert
  fluxBichAES=(Jsol(4,3,5)+Jsol(4,3,6))*convert
  
  fluxH2CO3MP=Jsol(5,1,2)*convert
  fluxH2CO3MA=Jsol(5,1,3)*convert
  fluxCO2MP=Jsol(6,1,2)*convert
  fluxCO2MA=Jsol(6,1,3)*convert
  
  fluxH2CO3PES=(Jsol(5,2,5)+Jsol(5,2,6))*convert
  fluxH2CO3AES=(Jsol(5,3,5)+Jsol(5,3,6))*convert
  fluxCO2PES=(Jsol(6,2,5)+Jsol(6,2,6))*convert
  fluxCO2AES=(Jsol(6,3,5)+Jsol(6,3,6))*convert
  
  fluxHP2mchPES=(Jsol(7,2,5)+Jsol(7,2,6))*convert
  fluxHP2mchAES=(Jsol(7,3,5)+Jsol(7,3,6))*convert
  fluxHPmPES=(Jsol(8,2,5)+Jsol(8,2,6))*convert
  fluxHPmAES=(Jsol(8,3,5)+Jsol(8,3,6))*convert
  
!---------------------------------------------------------------------72	
!---------------------------------------------------------------------72
!	 FLUXES ACROSS COTRANSPORTERS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


!	 Na2HPO4 cotransporter at PE,PS,AE,AS,BE,BS interfaces
  Do K = 2,3
     sumJES=0.0d0
     Do L = 5,6
        dJNaP=omcd%area(K,L)*omcd%dLA(1,7,K,L)*(2*delmu(1,K,L)+delmu(7,K,L)) 
        Jsol(1,K,L)=Jsol(1,K,L)+2*dJNaP
        Jsol(7,K,L)=Jsol(7,K,L)+dJNaP
        sumJES=sumJES+2*dJNaP
     End do
     if (K.eq.2) fluxNaPatPES=sumJES*convert
     if (K.eq.3) fluxNaPatAES=sumJES*convert
  End Do


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	 FLUXES ACROSS EXCHANGERS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	 NaH exchanger at PE,PS,AE,AS interfaces
!	 NHE1 model for principal cell ! Model of Fuster et al (J Gen Physiol 2009) - NOT USED

  Do K = 2,2
     sumJES=0.0d0
     Do L = 5,6
        
        affnao = 34.0d0
        affnai = 102.0d0
        affho = 0.0183d-3
        affhi = 0.054d-3
        fnai = C(1,K)
        hi = C(12,K)
        fnao = C(1,L)
        ho = C(12,L)
        Fno = (fnao/affnao)/(1.0d0+fnao/affnao+ho/affho)
        Fni = (fnai/affnai)/(1.0d0+fnai/affnai+hi/affhi)
        Fho = (ho/affho)/(1.0d0+fnao/affnao+ho/affho)
        Fhi = (hi/affhi)/(1.0d0+fnai/affnai+hi/affhi)
        E2mod1 = (Fni+Fhi)/(Fni+Fhi+Fno+Fho)
        E1mod1 = 1.0d0-E2mod1
        E2mod2 = (Fni**2+Fhi**2)/(Fni**2+Fhi**2+Fno**2+Fho**2)
        E1mod2 = 1.0d0-E2mod2
        Fmod1 = (hi**2)/(hi**2+(0.3d-3)**2)
        Rnhe = 1.0d3*(1.0d0*(1.0d0-Fmod1)*(E2mod2*(Fno**2)-E1mod2*(Fni**2)) &
             +Fmod1*(E2mod1*Fno-E1mod1*Fni))
        dJNaH = - omcd%area(K,L)*omcd%xNHE1(K)*Rnhe
!        Jsol(1,K,L) = Jsol(1,K,L)+dJNaH
!        Jsol(12,K,L)=Jsol(12,K,L)-dJNaH
        sumJES=sumJES+dJNaH
        
     End do
     fluxNaHPES = sumJES*convert
  End Do
  
!  NaH exchanger at PE and PS interfaces
  sumJES = 0.0d0
  Do L = 5,6
     dJNaHPES=omcd%area(2,L)*omcd%dLA(1,12,2,L)*(delmu(1,2,L)-delmu(12,2,L))
     Jsol(1,2,L) = Jsol(1,2,L)+dJNaHPES
     Jsol(12,2,L)=Jsol(12,2,L)-dJNaHPES
     sumJES=sumJES+dJNaHPES
  End do
  fluxNaHPES = sumJES*convert

!  NaH exchanger at AE and AS interfaces
  sumJES = 0.0d0
  Do L = 5,6
     dJNaHAES=omcd%area(3,L)*omcd%dLA(1,12,3,L)*(delmu(1,3,L)-delmu(12,3,L))
     Jsol(1,3,L) = Jsol(1,3,L)+dJNaHAES
     Jsol(12,3,L)=Jsol(12,3,L)-dJNaHAES
     sumJES=sumJES+dJNaHAES
  End do
  fluxNaHAES = sumJES*convert
  
  
  !	 Cl/HCO3 exchanger at PE,PS interfaces
  sumJES=0.0d0
  Do L = 5,6
     dJClHCO3=omcd%area(2,L)*omcd%dLA(3,4,2,L) *(delmu(3,2,L)-delmu(4,2,L))
     Jsol(3,2,L) = Jsol(3,2,L)+dJClHCO3
     Jsol(4,2,L) = Jsol(4,2,L)-dJClHCO3
     sumJES=sumJES+dJClHCO3
  End do
  fluxClHCO3exPES = sumJES*convert
  
  
!---------------------------------------------------------------------72
!		AE1 EXCHANGER AT PERITUBULAR MEMBRANE OF ALPHA CELL
!		b is for HCO3, c is for Cl
!		p (or prime) is for the internal compartment (A, or 3)
!	    pp (or double prime) is for the external compartment (E/S)
!---------------------------------------------------------------------72

  bpp = C(4,3)
  cpp = C(3,3)
  betapp = bpp/dKbpp
  gampp = cpp/dKcpp
  sumJES=0.0d0
  Do K = 5,6
     bp = C(4,K)
     cp = C(3,K)
     betap = bp/dKbp
     gamp = cp/dKcp
     xT=omcd%xAE1/(1+bpp/172.0d0)
     sum=(1+betap+gamp)*(Pbpp*betapp+Pcpp*gampp)
     sum=sum+(1+betapp+gampp)*(Pbp*betap+Pcp*gamp)
     befflux=omcd%area(3,K)*xT/sum*(Pbpp*betapp*Pcp*gamp-Pbp*betap*Pcpp*gampp)
     Jsol(4,3,K)=Jsol(4,3,K)+befflux
     Jsol(3,3,K)=Jsol(3,3,K)-befflux
     sumJES=sumJES-befflux
  End Do
  fluxAE1 = sumJES*convert

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	ACTIVE TRANSPORT VIA ATPases
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	Na-K-ATPase
!---------------------------------------------------------------------72


  Do M = 2, 3 ! For compartments P, A, B
     
     AffNa = 0.2d0*(1.0d0+C(2,M)/8.33d0)
     !	    AffNa = 16.3d0
     actNa = (C(1,M)/(C(1,M)+AffNa))
     AffK = 0.1d0*(1.0d0+C(1,6)/18.5d0)
     AffNH4 =  AffK/0.20d0 !Per Alan's email
     actK5 = (C(2,5)/(C(2,5)+AffK))
     actK6 = (C(2,6)/(C(2,6)+AffK))
     
     dJact5 = omcd%area(M,5)*omcd%ATPNaK(M,5)*(actNa**3.0d0) *(actK5**2.0d0)
     dJact6 = omcd%area(M,6)*omcd%ATPNaK(M,6)*(actNa**3.0d0) *(actK6**2.0d0)
     ro5 = (C(11,5)/AffNH4)/(C(2,5)/AffK)
     ro6 = (C(11,6)/AffNH4)/(C(2,6)/AffK)
     
     Jsol(1,M,5) = Jsol(1,M,5)+dJact5
     Jsol(1,M,6) = Jsol(1,M,6)+dJact6
     
     Jsol(2,M,5) = Jsol(2,M,5)-2.0d0/3.0d0*dJact5/(1+ro5)
     Jsol(2,M,6) = Jsol(2,M,6)-2.0d0/3.0d0*dJact6/(1+ro6)
     
     Jsol(11,M,5) = Jsol(11,M,5)-2.0d0/3.0d0*dJact5*ro5/(1+ro5)
     Jsol(11,M,6) = Jsol(11,M,6)-2.0d0/3.0d0*dJact6*ro6/(1+ro6)
     
     if (M.eq.2) fluxNaKasePES = (dJact5+dJact6)*convert
     if (M.eq.3) fluxNaKaseAES = (dJact5+dJact6)*convert
     if (M.eq.4) fluxNaKaseBES = (dJact5+dJact6)*convert
     
  End Do
  
!---------------------------------------------------------------------72
!	 H-ATPase
!	 See Strieter & Weinstein paper for signs (AJP 263, 1992)
!---------------------------------------------------------------------72

  denom13=1.0d0+dexp(steepA*(delmu(12,1,3)-dmuATPH))
  dJact13=-omcd%area(1,3)*omcd%ATPH(1,3)/denom13
  Jsol(12,1,3) = Jsol(12,1,3)+dJact13
  fluxHATPaseMA = (dJact13)*convert
  
!---------------------------------------------------------------------72
!	 H-K-ATPase	at MP and MA interfaces
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	 H-K-ATPase	at MP interface
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
  hefflux = omcd%area(1,2)*omcd%ATPHK(1,2)*(dkf5*c7-dkb5*c8)
  Jsol(2,1,2) = Jsol(2,1,2)+2*hefflux
  Jsol(12,1,2) = Jsol(12,1,2)-2*hefflux
  fluxHKATPaseMP = 2*hefflux*convert
  
!---------------------------------------------------------------------72
!	 H-K-ATPase	at MA interface
!---------------------------------------------------------------------72

  hkconc(1) = C(2,3) ![K]in = K in A
  hkconc(2) = C(2,1) ![K]out = K in M
  hkconc(3) = C(12,3) ![H]in = H in A
  hkconc(4) = C(12,1) ![H]out = H in M
  
  call fatpase(Natp,hkconc,Amat)
  
  call dgetrf( Natp, Natp, Amat, Natp, ipiv, info )
  call dgetri( Natp, Amat, Natp, ipiv, work, lwork, info )
  
  c7 = Amat(7,1)
  c8 = Amat(8,1)
  dkf5 = 4.0d1
  dkb5 = 2.0d2
  hefflux = omcd%area(1,3)*omcd%ATPHK(1,3)*(dkf5*c7-dkb5*c8)
  Jsol(2,1,3) = Jsol(2,1,3)+2*hefflux
  Jsol(12,1,3) = Jsol(12,1,3)-2*hefflux
  fluxHKATPaseMA = 2*hefflux*convert


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	 Cancel all IC-B fluxes
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  Jvol(1,4) = 0.0d0
  Jvol(4,5) = 0.0d0
  Jvol(4,6) = 0.0d0
  Do I = 1, NS
     Jsol(I,1,4) = 0.0d0	
     Jsol(I,4,5) = 0.0d0
     Jsol(I,4,6) = 0.0d0
  End Do

!---------------------------------------------------------------------72
!     STORE LOCAL FLUX VALUES FOR LATER INTEGRATION
!---------------------------------------------------------------------72

  cv=href*Cref*PI*DiamOMC*60/10*1.0d9 !for dimensional solute fluxes
  cvw=Pfref*Cref*Vwbar*PI*DiamOMC*60/10*1.0d6 !for dimensional water flux

  omcd%FNatrans = (Jsol(1,1,2)+Jsol(1,1,3))*cv*omcd%coalesce
  omcd%FNapara = Jsol(1,1,5)*cv*omcd%coalesce
  omcd%FNaK = (fluxNaKasePES+fluxNaKaseAES+fluxNaKaseBES)*omcd%coalesce
  
  omcd%FKtrans = (Jsol(2,1,2)+Jsol(2,1,3))*cv*omcd%coalesce
  omcd%FKpara = Jsol(2,1,5)*cv*omcd%coalesce

  return
end Subroutine qflux2OMC


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
