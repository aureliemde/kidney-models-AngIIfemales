!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
 
!  *************************** FLUXES ***************************
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!      This subroutine computes the fluxes in the CNT

Subroutine qflux2C(x,Jvol,Jsol,cnt)

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

  ! passed variables
  type (membrane) :: cnt
  double precision x(NDC)
  double precision Jvol(NC,NC), Jsol(NS,NC,NC)
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

!   for NCX fluxes
  double precision var_ncx(16)
  double precision dJNCXca5,dJNCXca6

!   for TRPV5 current
  double precision k1v5,k2v5,k3v5,k4v5
  double precision pcv5,pfv5,psv5,gfv5,gsv5,gv5
  double precision ECa

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	ASSIGN CONCENTRATIONS, VOLUMES, AND POTENTIALS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  Do J = 1, NS
     C(J,6) = cnt%conc(J,6)  
  End Do
  EP(6) = cnt%ep(6)  
	
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

  cnt%area(5,6)=cnt%sbasEinit*max(Vol(5)/cnt%volEinit,1.0d0)
  cnt%area(6,5)=cnt%area(5,6)
	
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

  call compute_water_fluxes (C,PM,0.0d0,Vol,cnt%volLuminit,cnt%volEinit,  &
       cnt%volPinit,CPimprefC,cnt%volAinit,CAimprefC,cnt%volBinit,   &
	   CBimprefC,cnt%area,cnt%sig,cnt%dLPV,complC,Jvol)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     COMPUTE SOLUTE FLUXES
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	 Non-dimensional solute fluxes are first converted to dimensional fluxes 
!      (in mmol/s/cm2 epith) by multiplying by href*Cref
!	 Then they are converted to units of pmol/min/mm tubule 
  convert=href*Cref*PI*DiamC*60/10*1.0d9 !Conversion factor

!---------------------------------------------------------------------72
!     ELECTRO-CONVECTIVE-DIFFUSIVE FLUXES
!---------------------------------------------------------------------72

!    dependence of ENaC, ROMK, and apical paracellular Cl permeability on pH

   facphMP = 1.0*(0.1 + 2.0d0/ (1+dexp(-6.0d0*(pH(2)-7.50d0))))
   facphTJ = 2.0/(1.0 + dexp(10.0*(ph(5)-7.32d0)))

   facNaMP=(30.d0/(30.d0+C(1,1)))*(50.d0/(50.d0+C(1,2)))
 !  cnt%h(1,1,2)=cnt%hNaMP*facNaMP*facphMP
   cnt%h(1,1,2)=hENaC_CNT*facNaMP*facphMP

   cnt%h(2,1,2)=hROMK_CNT*facphMP

   cnt%h(3,1,5)=hCltj_CNT*facphTJ

!	 define electrochemical potential of each solute in each compartment

  call compute_ecd_fluxes (C,EP,cnt%area,cnt%sig,cnt%h,Jvol,Jsol,delmu)

!	    dimensional fluxes through channels, for print out


  fluxENaC=Jsol(1,1,2)*convert
  fluxROMK=Jsol(2,1,2)*convert
  
  fluxKchPES=(Jsol(2,2,5)+Jsol(2,2,6))*convert
  fluxKchAES=(Jsol(2,3,5)+Jsol(2,3,6))*convert
  fluxKchBES=(Jsol(2,4,5)+Jsol(2,4,6))*convert
  fluxClchPES=(Jsol(3,2,5)+Jsol(3,2,6))*convert
  fluxClchAES=(Jsol(3,3,5)+Jsol(3,3,6))*convert
  fluxClchBES=(Jsol(3,4,5)+Jsol(3,4,6))*convert
  fluxBichPES=(Jsol(4,2,5)+Jsol(4,2,6))*convert
  fluxBichAES=(Jsol(4,3,5)+Jsol(4,3,6))*convert
  fluxBichBES=(Jsol(4,4,5)+Jsol(4,4,6))*convert
  
  fluxH2CO3MP=Jsol(5,1,2)*convert
  fluxH2CO3MA=Jsol(5,1,3)*convert
  fluxH2CO3MB=Jsol(5,1,4)*convert
  fluxCO2MP=Jsol(6,1,2)*convert
  fluxCO2MA=Jsol(6,1,3)*convert
  fluxCO2MB=Jsol(6,1,4)*convert
  
  fluxH2CO3PES=(Jsol(5,2,5)+Jsol(5,2,6))*convert
  fluxH2CO3AES=(Jsol(5,3,5)+Jsol(5,3,6))*convert
  fluxH2CO3BES=(Jsol(5,4,5)+Jsol(5,4,6))*convert
  fluxCO2PES=(Jsol(6,2,5)+Jsol(6,2,6))*convert
  fluxCO2AES=(Jsol(6,3,5)+Jsol(6,3,6))*convert
  fluxCO2BES=(Jsol(6,4,5)+Jsol(6,4,6))*convert
  
  fluxHP2mchPES=(Jsol(7,2,5)+Jsol(7,2,6))*convert
  fluxHP2mchAES=(Jsol(7,3,5)+Jsol(7,3,6))*convert
  fluxHP2mchBES=(Jsol(7,4,5)+Jsol(7,4,6))*convert
  fluxHPmPES=(Jsol(8,2,5)+Jsol(8,2,6))*convert
  fluxHPmAES=(Jsol(8,3,5)+Jsol(8,3,6))*convert
  fluxHPmBES=(Jsol(8,4,5)+Jsol(8,4,6))*convert


!---------------------------------------------------------------------72	
!---------------------------------------------------------------------72
!	 FLUXES ACROSS COTRANSPORTERS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	 Na2HPO4 cotransporter at PE,PS,AE,AS,BE,BS interfaces
  Do K = 2,4
     sumJES=0.0d0
     Do L = 5,6
        dJNaP=cnt%area(K,L)*cnt%dLA(1,7,K,L)  &
             *(2*delmu(1,K,L)+delmu(7,K,L)) 
        Jsol(1,K,L)=Jsol(1,K,L)+2*dJNaP
        Jsol(7,K,L)=Jsol(7,K,L)+dJNaP
        sumJES=sumJES+2*dJNaP
     End do
     if (K.eq.2) fluxNaPatPES=sumJES*convert
     if (K.eq.3) fluxNaPatAES=sumJES*convert
     if (K.eq.4) fluxNaPatBES=sumJES*convert
  End Do


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	 FLUXES ACROSS EXCHANGERS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


!	 NaH exchanger at PE,PS,AE,AS,BE,BS interfaces

  Do K = 2,4		
     sumJES=0.0d0
     Do L = 5,6
        dJNaH=cnt%area(K,L)*cnt%dLA(1,12,K,L)  &
             *(delmu(1,K,L)-delmu(12,K,L))
        
        Jsol(1,K,L) = Jsol(1,K,L)+dJNaH
        Jsol(12,K,L)=Jsol(12,K,L)-dJNaH
        sumJES=sumJES+dJNaH
     End do
     if (K.eq.2) fluxNaHPES = sumJES*convert
     if (K.eq.3) fluxNaHAES = sumJES*convert
     if (K.eq.4) fluxNAHBES = sumJES*convert
  End Do


!	 Cl/HCO3 exchanger at PE,PS interfaces
  sumJES=0.0d0
  Do L = 5,6
     dJClHCO3=cnt%area(2,L)*cnt%dLA(3,4,2,L)  &
          *(delmu(3,2,L)-delmu(4,2,L))
     Jsol(3,2,L) = Jsol(3,2,L)+dJClHCO3
     Jsol(4,2,L) = Jsol(4,2,L)-dJClHCO3
     sumJES=sumJES+dJClHCO3
  End do
  fluxClHCO3exPES = sumJES*convert

!	 Cl/HCO3 exchanger at BE,BS interfaces - REMOVED
  sumJES=0.0d0
  Do L = 5,6
     dJClHCO3=cnt%area(4,L)*cnt%dLA(3,4,4,L)  &
          *(delmu(3,4,L)-delmu(4,4,L))
     Jsol(3,4,L) = Jsol(3,4,L)+dJClHCO3
     Jsol(4,4,L) = Jsol(4,4,L)-dJClHCO3
     sumJES=sumJES+dJClHCO3
  End do
  fluxClHCO3exBES = sumJES*convert


!---------------------------------------------------------------------72
!	 Pendrin exchanger at the apical membrane of the beta cell
!---------------------------------------------------------------------72

  Pbiiclepd = Pbieclepd*Pcliclepd !assuming internal = external affinity
  Pohiclepd = Poheclepd*Pcliclepd !assuming internal = external affinity

  ohe = 1.0d-11/C(12,1) !H+ concentrations are in mM, OH- are in M
  ohi = 1.0d-11/C(12,4) !H+ concentrations in mM, OH- in M 
  alpe = ohe/dKohpd
  alpi = ohi/dKohpd
  game = C(3,1)/dKclpd
  gami = C(3,4)/dKclpd
  bete = C(4,1)/dKbipd
  beti = C(4,4)/dKbipd
  
  dele = 1.0d0+game+bete+alpe
  deli = 1.0d0+gami+beti+alpi
  etae = (1.0d0*game+Pbieclepd*bete+Poheclepd*alpe)
  etai = (Pcliclepd*gami+Pbiiclepd*beti+Pohiclepd*alpi)
  sigma = dele*etai+deli*etae
  
!		Stimulation by low pHi
  if (ph(4).lt.7.9d0) then
     adj=dexp((1.d0/6.4189)*dlog((7.9d0-ph(4))/0.0115))
  else 
     adj = 1.0d0
  end if
  
  factpd = cnt%area(1,4)*cnt%xPendrin*Pclepd*adj
  
  dJclpd=factpd*(1.0d0*game*etai-Pcliclepd*gami*etae)/sigma
  dJbipd=factpd*(Pbieclepd*bete*etai-Pbiiclepd*beti*etae)/sigma
  dJohpd=factpd*(Poheclepd*alpe*etai-Pohiclepd*alpi*etae)/sigma
  
  Jsol(3,1,4) = Jsol(3,1,4)+dJclpd
  Jsol(4,1,4) = Jsol(4,1,4)+dJbipd
  fluxPdMB = dJclpd*convert

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
     xT=cnt%xAE1/(1+bpp/172.0d0)
     sum=(1+betap+gamp)*(Pbpp*betapp+Pcpp*gampp)
     sum=sum+(1+betapp+gampp)*(Pbp*betap+Pcp*gamp)
     befflux=cnt%area(3,K)*xT/sum  &
          *(Pbpp*betapp*Pcp*gamp-Pbp*betap*Pcpp*gampp)
     Jsol(4,3,K)=Jsol(4,3,K)+befflux
     Jsol(3,3,K)=Jsol(3,3,K)-befflux
     sumJES=sumJES-befflux
  End Do
  fluxAE1 = sumJES*convert

!---------------------------------------------------------------------72
!       NCX EXCHANGER AT BASOLATERAL MEMBRANE OF CELL
!---------------------------------------------------------------------72

  var_ncx(1) = C(1,2)
  var_ncx(2) = C(1,5)
  var_ncx(3) = C(1,6)
  var_ncx(4) = C(16,2)
  var_ncx(5) = C(16,5)
  var_ncx(6) = C(16,6)
  var_ncx(7) = EP(2)
  var_ncx(8) = EP(5)
  var_ncx(9) = EP(6)
  var_ncx(10) = cnt%area(2,5)
  var_ncx(11) = cnt%area(2,6)
  var_ncx(12) = cnt%xNCX

  call compute_nCX_fluxes(var_ncx,dJNCXca5,dJNCXca6)

  Jsol(1,2,5) = Jsol(1,2,5) - 3.0d0*dJNCXca5
  Jsol(1,2,6) = Jsol(1,2,6) - 3.0d0*dJNCXca6
  Jsol(16,2,5) = Jsol(16,2,5) + dJNCXca5
  Jsol(16,2,6) = Jsol(16,2,6) + dJNCXca6
  fluxNCXca = (dJNCXca5+dJNCXca6)*convert


!---------------------------------------------------------------------72
!    Ca2+ flux across TRPV5
!---------------------------------------------------------------------72

    k1v5 = 42.7d0
    k3v5 = 0.1684*dexp(0.6035*ph(2))
    k4v5 = 58.7d0
    if (ph(2) .lt. 7.4d0) then
        k2v5 = 55.9d0 + (173.3-55.9)/(7.0-7.4)*(ph(2)-7.4d0)
     else
        k2v5 = 55.9d0 + (30.4-55.9)/(8.4-7.4)*(ph(2)-7.4d0)
    end if
    psv5 = 1.0d0/(1.0d0+k3v5/k4v5+k2v5/k1v5)
    pcv5 = k2v5/k1v5*psv5
    pfv5 = k3v5/k4v5*psv5

    gfv5 = 59.d0 + 59.d0/(59.d0+29.0d0)*(91.0-58.0)/(7.4-5.4)*(ph(1)-7.4)
    gsv5 = 29.d0 + 29.d0/(59.d0+29.0d0)*(91.0-58.0)/(7.4-5.4)*(ph(1)-7.4)
    gv5 = (gfv5*pfv5+gsv5*psv5)*1.0d-12

    ECa = RT/(2*F)*dlog(C(16,2)/C(16,1)) ! in volts (RT/F in J/C)
    dfv5 = (EP(1)-EP(2))*EPref - ECa
    finhib_v5 = 1.0d0/(1.0d0+C(16,2)/Cinhib_v5)

    dJTRPV5 = cnt%area(1,2)*xTRPV5_cnt*finhib_v5*gv5*dfv5/(2*F)

    Jsol(16,1,2) = Jsol(16,1,2)+dJTRPV5

    fluxTRPV5 = dJTRPV5*convert

!---------------------------------------------------------------------72
!    Ca2+ flux across TRPV4
!---------------------------------------------------------------------72

    XICa = zval(16)*F*EPref/RT*(EP(1)-EP(2))
    dint = dexp(-XICa)
    if (dabs(1.0d0-dint) .lt. 1.0d-6) then
       Df_TRPV4 = xPTRPV4_cnt * (C(16,1)-C(16,2))
       else
       Df_TRPV4 = xPTRPV4_cnt * XICa *(C(16,1)-C(16,2)*dint) /(1.0d0-dint)
    end if

    dJTRPV4 = cnt%area(1,2)*Df_TRPV4

    Jsol(16,1,2) = Jsol(16,1,2)+dJTRPV4

    fluxTRPV4 = dJTRPV4*convert


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	ACTIVE TRANSPORT VIA ATPases
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	Na-K-ATPase
!---------------------------------------------------------------------72


  Do M = 2, 4 ! For compartments P, A, B

     AffNa = 0.2d0*(1.0d0+C(2,M)/8.33d0)
     actNa = (C(1,M)/(C(1,M)+AffNa))
     AffK = 0.1d0*(1.0d0+C(1,6)/18.5d0)
     AffNH4 =  AffK/0.20d0 !Per Alan's email
     actK5 = (C(2,5)/(C(2,5)+AffK))
     actK6 = (C(2,6)/(C(2,6)+AffK))
     
     dJact5 = cnt%area(M,5)*cnt%ATPNaK(M,5)*(actNa**3.0d0)  &
          *(actK5**2.0d0)
     dJact6 = cnt%area(M,6)*cnt%ATPNaK(M,6)*(actNa**3.0d0)  &
          *(actK6**2.0d0)
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
  denom45=1.0d0+dexp(-steepB*(delmu(12,4,5)-dmuATPH)) 
  denom46=1.0d0+dexp(-steepB*(delmu(12,4,6)-dmuATPH))
  dJact13=-cnt%area(1,3)*cnt%ATPH(1,3)/denom13
  dJact45=cnt%area(4,5)*cnt%ATPH(4,5)/denom45
  dJact46=cnt%area(4,6)*cnt%ATPH(4,6)/denom46
  
  Jsol(12,1,3) = Jsol(12,1,3)+dJact13
  Jsol(12,4,5) = Jsol(12,4,5)+dJact45
  Jsol(12,4,6) = Jsol(12,4,6)+dJact46
  
  fluxHATPaseMA = (dJact13)*convert
  fluxHATPaseBES = (dJact45+dJact46)*convert
  
!---------------------------------------------------------------------72
!	 H-K-ATPase	at MP, MA and MB interfaces
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
  hefflux = cnt%area(1,2)*cnt%ATPHK(1,2)*(dkf5*c7-dkb5*c8)
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
  hefflux = cnt%area(1,3)*cnt%ATPHK(1,3)*(dkf5*c7-dkb5*c8)
  Jsol(2,1,3) = Jsol(2,1,3)+2*hefflux
  Jsol(12,1,3) = Jsol(12,1,3)-2*hefflux
  fluxHKATPaseMA = 2*hefflux*convert
  
!---------------------------------------------------------------------72
!	 H-K-ATPase	at MB interface
!---------------------------------------------------------------------72

  hkconc(1) = C(2,4) ![K]in = K in B
  hkconc(2) = C(2,1) ![K]out = K in M
  hkconc(3) = C(12,4) ![H]in = H in B
  hkconc(4) = C(12,1) ![H]out = H in M
  
  call fatpase(Natp,hkconc,Amat)
  
  call dgetrf( Natp, Natp, Amat, Natp, ipiv, info )
  call dgetri( Natp, Amat, Natp, ipiv, work, lwork, info )
  
  c7 = Amat(7,1)
  c8 = Amat(8,1)
  dkf5 = 4.0d1
  dkb5 = 2.0d2
  hefflux = cnt%area(1,4)*cnt%ATPHK(1,4)*(dkf5*c7-dkb5*c8)
  Jsol(2,1,4) = Jsol(2,1,4)+2*hefflux
  Jsol(12,1,4) = Jsol(12,1,4)-2*hefflux
  fluxHKATPaseMB = 2*hefflux*convert

!---------------------------------------------------------------------72
!   PMCA Ca2+ on basolateral membrane
!---------------------------------------------------------------------72

   ratio = C(16,2)/(C(16,2)+dKmPMCA)
   dJPMCA5 = cnt%PMCA*cnt%area(2,5)*ratio
   dJPMCA6 = cnt%PMCA*cnt%area(2,6)*ratio

   Jsol(16,2,5) = Jsol(16,2,5) + dJPMCA5
   Jsol(16,2,6) = Jsol(16,2,6) + dJPMCA6

   fluxPMCA = (dJPMCA5+dJPMCA6)*convert

!---------------------------------------------------------------------72
!     STORE LOCAL FLUX VALUES FOR LATER INTEGRATION
!---------------------------------------------------------------------72


  cv=href*Cref*PI*DiamC*60/10*1.0d9 !for dimensional solute fluxes
  cvw=Pfref*Cref*Vwbar*PI*DiamC*60/10*1.0d6 !for dimensional water flux

  cnt%FNatrans = (Jsol(1,1,2)+Jsol(1,1,3)+Jsol(1,1,4))*cv*cnt%coalesce
  cnt%FNapara = Jsol(1,1,5)*cv*cnt%coalesce
  cnt%FNaK = (fluxNaKasePES+fluxNaKaseAES+fluxNaKaseBES)*cnt%coalesce

  cnt%FKtrans = (Jsol(2,1,2)+Jsol(2,1,3)+Jsol(2,1,4))*cv*cnt%coalesce
  cnt%FKpara = Jsol(2,1,5)*cv*cnt%coalesce


  fTRPV5_cnt(LzC+1) = fluxTRPV5*cnt%coalesce
  fNCX_cnt(LzC+1) = fluxNCXca*cnt%coalesce
  fPMCA_cnt(LzC+1) = fluxPMCA*cnt%coalesce
  fTRPV4_cnt(LzC+1) = fluxTRPV4*cnt%coalesce

  return
end Subroutine qflux2C


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
