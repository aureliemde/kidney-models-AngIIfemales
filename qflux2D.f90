!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
 
!  *************************** FLUXES ***************************
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!      This is a subroutine to compute the fluxes in the DCT

Subroutine qflux2D(x,Jvol,Jsol,Cext,EPext,dct)

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

  type (membrane) :: dct
  double precision Cext(NS),EPext
  double precision x(NDA)
  double precision Jvol(NC,NC), Jsol(NS,NC,NC)
  double precision ONC(NC), PRES(NC)
  double precision C(NS,NC),dmu(NS,NC),ph(NC),Vol(NC),EP(NC)
  double precision delmu(NS,NC,NC)
  double precision hkconc(4),Amat(Natp,Natp)
  double precision theta(NC),Slum(NC),Slat(NC),Sbas(NC)
  double precision xdct2,xl,xn,ENaCexp,NCCexp,CaTexp

!	for HKATPase matrix inversion
  integer info, lwork
  parameter (lwork=10000)
  integer ipiv(Natp)
  double precision work(lwork)

!   for NCX fluxes
  double precision var_ncx(16)
  double precision dJNCXca5,dJNCXca6


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

  dct%area(5,6)=dct%sbasEinit*max(Vol(5)/dct%volEinit,1.0d0)
  dct%area(6,5)=dct%area(5,6)

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

  call compute_water_fluxes (C,PM,0.0d0,Vol,dct%volLuminit,dct%volEinit,  &
       dct%volPinit,CPimprefD,dct%volAinit,CAimprefD,dct%volBinit,   &
	   CBimprefD,dct%area,dct%sig,dct%dLPV,complD,Jvol)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     COMPUTE SOLUTE FLUXES
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	 Non-dimensional solute fluxes are first converted to dimensional fluxes 
!      (in mmol/s/cm2 epith)by multiplying by href*Cref
!	 Then they are converted to units of pmol/min/mm tubule 
  convert=href*Cref*PI*DiamD*60/10*1.0d9 !Conversion factor
  eps = 1.0d-6 !For exponential (Peclet) terms

!---------------------------------------------------------------------72
!     PARTIAL OVERLAP OF NCC AND ENAC IN THE DCT2
!     NO Calcium transporters in DCT1
!---------------------------------------------------------------------72

  xl = LzD+1.0d0 !convert to real number
  xn = NZ  !convert to real number
  xdct2 = 2.0d0/3.d0 ! DCT1: first 2/3rds of length, DCT2 last 1/3rd
  if (xl/xn .lt. xdct2) then
     NCCexp = 1.0d0
     ENaCexp = 0.0d0
     CaTexp = 0.001d0
  else
     NCCexp = (1.0d0 - (xl/xn - xdct2)/(1.0d0-xdct2))
     ENaCexp = (xl/xn - xdct2)/(1.0d0-xdct2) 
     CaTexp = 1.0d0
  endif

!---------------------------------------------------------------------72
!     ELECTRO-CONVECTIVE-DIFFUSIVE FLUXES
!---------------------------------------------------------------------72

!	 define electrochemical potential of each solute in each compartment
	
  Do I = 1, NS
     Do K = 1, NC
        dmu(I,K) = RT*dlog(abs(C(I,K)))+zval(I)*F*EPref*EP(K)
     End Do
  End Do

!	 begin flux calculation
  
  Do I = 1, NS2
     Do K = 1, NC-1
        Do L = K+1, NC

!		electrodiffusive component

           XI = zval(I)*F*EPref/RT*(EP(K)-EP(L))
           dint = dexp(-XI)
           if (dabs(1.0d0-dint).lt.eps) then
              Jsol(I,K,L)=dct%area(K,L)*dct%h(I,K,L)*(C(I,K)-C(I,L))
           else 
              Jsol(I,K,L)=dct%area(K,L)*dct%h(I,K,L)*XI  &
                   *(C(I,K)-C(I,L)*dint)/(1.0d0-dint)
           end if


!		convective component

		 concdiff=C(I,K)-C(I,L)
		 if (dabs(concdiff) > eps)	then
			concmean=(C(I,K)-C(I,L))/dlog(abs(C(I,K)/C(I,L)))
			dimless=(Pfref*Vwbar*Cref)/href
			convect=(1.0d0-dct%sig(I,K,L))*concmean*Jvol(K,L)*dimless
			Jsol(I,K,L) = Jsol(I,K,L)+convect
	     end if


!		define driving force for coupled fluxes below

		 delmu(I,K,L)=dmu(I,K)-dmu(I,L) !for coupled flux calculations
	  
		End Do
	   End Do
	 End Do

!	    dimensional fluxes through channels, for print out

		fluxNachMP = Jsol(1,1,2)*convert

		fluxKchMP=Jsol(2,1,2)*convert
		fluxKchPES=(Jsol(2,2,5)+Jsol(2,2,6))*convert

		fluxClchPES=(Jsol(3,2,5)+Jsol(3,2,6))*convert
		fluxBichPES=(Jsol(4,2,5)+Jsol(4,2,6))*convert

		fluxH2CO3MP=Jsol(5,1,2)*convert
		fluxCO2MP=Jsol(6,1,2)*convert

		fluxH2CO3PES=(Jsol(5,2,5)+Jsol(5,2,6))*convert
		fluxCO2PES=(Jsol(6,2,5)+Jsol(6,2,6))*convert

		fluxHP2mchPES=(Jsol(7,2,5)+Jsol(7,2,6))*convert
		fluxHPmPES=(Jsol(8,2,5)+Jsol(8,2,6))*convert

        Jsol(16,1,2) = CaTexp*Jsol(16,1,2) ! Ca2+ transport in DCT2 only
        Jsol(16,2,5) = CaTexp*Jsol(16,2,5) ! Ca2+ transport in DCT2 only
        Jsol(16,2,5) = CaTexp*Jsol(16,2,6) ! Ca2+ transport in DCT2 only


!---------------------------------------------------------------------72
!	 ENaC in DCT2 only
!---------------------------------------------------------------------72

!	 dependence of ENaC permeability to Na 

	 facNaMP=(30.d0/(30.d0+C(1,1)))*(50.d0/(50.d0+C(1,2)))
     facphMP = 1.0*(0.1 + 2.0d0/ (1+dexp(-6.0d0*(pH(2)-7.50d0))))
	 hENaC=ENaCexp*dct%hNaMP*facNaMP*facphMP !ENaCexp is zero in early DCT

	 XI = zval(1)*F*EPref/RT*(EP(1)-EP(2))
	 dint = dexp(-XI)
	 if (dabs(1.0d0-dint).lt.eps) then
		 dJENaC=dct%area(1,2)*hENaC*(C(1,1)-C(1,2))
	  else 
		 dJENaC=dct%area(1,2)*hENaC*XI*(C(1,1)-C(1,2)*dint)/(1.0d0-dint)
	 end if

	 Jsol(1,1,2) = Jsol(1,1,2) + dJENaC

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
			dJNaP=dct%area(2,L)*dct%dLA(1,7,2,L)* &
      			  (2*delmu(1,2,L)+delmu(7,2,L)) 
			Jsol(1,2,L)=Jsol(1,2,L)+2*dJNaP
			Jsol(7,2,L)=Jsol(7,2,L)+dJNaP
			sumJES=sumJES+2*dJNaP
		End do
		fluxNaPatPES=sumJES*convert

!---------------------------------------------------------------------72
!	 Apical and basolateral KCl cotransporter 
!---------------------------------------------------------------------72

		dJKClMP=dct%area(1,2)*dct%dLA(2,3,1,2)* (delmu(2,1,2)+delmu(3,1,2))
		Jsol(2,1,2)=Jsol(2,1,2)+dJKClMP
		Jsol(3,1,2)=Jsol(3,1,2)+dJKClMP
		fluxKClMP=dJKClMP*convert

		sumJES = 0.d0
		Do L = 5,6
			dJKCl=dct%area(2,L)*dct%dLA(2,3,2,L)* (delmu(2,2,L)+delmu(3,2,L))
			Jsol(2,2,L)=Jsol(2,2,L)+dJKCl
			Jsol(3,2,L)=Jsol(3,2,L)+dJKCl
			sumJES = sumJES + dJKCl
		End do
		fluxKClPES=sumJES*convert

!---------------------------------------------------------------------72
!	 NCC cotransporter at the apical membrane of principal cells
!---------------------------------------------------------------------72

		alp = C(1,1)/dKnncc
		alpp = C(1,2)/dKnncc
		betp = C(3,1)/dKcncc
	    betpp = C(3,2)/dKcncc

	    gamp = C(1,1)*C(3,1)/(dKnncc*dKncncc)
	    gampp = C(1,2)*C(3,2)/(dKnncc*dKncncc)

		rhop = 1.d0 + alp + betp + gamp
	    rhopp = 1.d0 + alpp + betpp + gampp
		sigma = rhop*(poppncc+gampp*pnppncc) +rhopp*(popncc+gamp*pnpncc)
	    dJNCC = NCCexp*dct%area(1,2)*dct%xNCC*pnpncc*poppncc*(gamp-gampp)/sigma

		Jsol(1,1,2) = Jsol(1,1,2) + dJNCC
	    Jsol(3,1,2) = Jsol(3,1,2) + dJNCC
	    fluxNCC = dJNCC*convert

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	 FLUXES ACROSS EXCHANGERS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


!---------------------------------------------------------------------72
!		NHE EXCHANGER AT LUMINAL MEMBRANE OF CELL
!		n is for Na, c is for Cl
!		p (or prime) is for the luminal compartment (M)
!	    pp (or double prime) is for the cytosolic compartment (I)
!---------------------------------------------------------------------72

  call compute_nhe3_fluxes (C,dct%area(1,2),dct%xNHE3,dJNHEsod,dJNHEprot,dJNHEamm)

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

        			dJNaH=dct%area(K,L)*dct%dLA(1,12,K,L) &
      				  *(delmu(1,K,L)-delmu(12,K,L))
			goto 10	! Do not use Fuster model

!		NHE1 model for principal cell
!		Model of Fuster et al (J Gen Physiol 2009)
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
				Rnhe = 1.0d3*(1.0d0*(1.0d0-Fmod1)* &
      				(E2mod2*(Fno**2)-E1mod2*(Fni**2)) &
      				+Fmod1*(E2mod1*Fno-E1mod1*Fni))

10				Jsol(1,K,L) = Jsol(1,K,L)+dJNaH
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
			dJClHCO3=dct%area(2,L)*dct%dLA(3,4,2,L)* (delmu(3,2,L)-delmu(4,2,L))
			Jsol(3,2,L) = Jsol(3,2,L)+dJClHCO3
			Jsol(4,2,L) = Jsol(4,2,L)-dJClHCO3
			sumJES=sumJES+dJClHCO3
		End do
		fluxClHCO3exPES = sumJES*convert

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
  var_ncx(10) = dct%area(2,5)
  var_ncx(11) = dct%area(2,6)
  var_ncx(12) = dct%xNCX

  call compute_nCX_fluxes(var_ncx,dJNCXca5,dJNCXca6)

  dJNCXca5 = CaTexp*dJNCXca5
  dJNCXca6 = CaTexp*dJNCXca6

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

    dJTRPV5 = CaTexp*dct%area(1,2)*xTRPV5_dct*finhib_v5*gv5*dfv5/(2*F)

    Jsol(16,1,2) = Jsol(16,1,2)+dJTRPV5

    fluxTRPV5 = dJTRPV5*convert

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

		dJact5 = dct%area(2,5)*dct%ATPNaK(2,5)*(actNa**3.0d0) *(actK5**2.0d0)
		dJact6 = dct%area(2,6)*dct%ATPNaK(2,6)*(actNa**3.0d0) *(actK6**2.0d0)

		ro5 = (C(11,5)/AffNH4)/(C(2,5)/AffK)
		ro6 = (C(11,6)/AffNH4)/(C(2,6)/AffK)

		Jsol(1,2,5) = Jsol(1,2,5)+dJact5
                Jsol(1,2,6) = Jsol(1,2,6)+dJact6

		Jsol(2,2,5) = Jsol(2,2,5)-2.0d0/3.0d0*dJact5/(1.d0+ro5)
		Jsol(2,2,6) = Jsol(2,2,6)-2.0d0/3.0d0*dJact6/(1.d0+ro6)

		Jsol(11,2,5) = Jsol(11,2,5)-2.0d0/3.0d0*dJact5*ro5/(1.d0+ro5)
		Jsol(11,2,6) = Jsol(11,2,6)-2.0d0/3.0d0*dJact6*ro6/(1.d0+ro6)

		fluxNaKPESsod = (dJact5+dJact6)*convert

		fluxNaKPESpot = -2.d0/3.0d0*(dJact5/(1.d0+ro5)+ &
      					  dJact6/(1.d0+ro6))*convert
		fluxNaKPESamm = -2.d0/3.0d0*(dJact5*ro5/(1.d0+ro5)+  &
      					  dJact6*ro6/(1.d0+ro6))*convert


!---------------------------------------------------------------------72
!   PMCA Ca2+ PUMP AT BASOLATERAL MEMBRANE OF CELL
!---------------------------------------------------------------------72

       ratio = C(16,2)/(C(16,2)+dKmPMCA)
       dJPMCA5 = CaTexp*dct%PMCA*dct%area(2,5)*ratio
       dJPMCA6 = CaTexp*dct%PMCA*dct%area(2,6)*ratio

       Jsol(16,2,5) = Jsol(16,2,5) + dJPMCA5
       Jsol(16,2,6) = Jsol(16,2,6) + dJPMCA6

       fluxPMCA = (dJPMCA5+dJPMCA6)*convert

!---------------------------------------------------------------------72
!     STORE LOCAL FLUX VALUES FOR LATER INTEGRATION
!---------------------------------------------------------------------72

    cv=href*Cref*PI*DiamD*60/10*1.0d9 !for dimensional solute fluxes
    cvw=Pfref*Cref*Vwbar*PI*DiamD*60/10*1.0d6 !for dimensional water flux

	dct%FNatrans = Jsol(1,1,2)*cv
	dct%FNapara = Jsol(1,1,5)*cv
    dct%FNaK = fluxNaKPESsod
    dct%FHase = 0.0

	dct%FKtrans = Jsol(2,1,2)*cv
	dct%FKpara = Jsol(2,1,5)*cv

    fTRPV5_dct(LzD+1) = fluxTRPV5
    fNCX_dct(LzD+1) = fluxNCXca
    fPMCA_dct(LzD+1) = fluxPMCA

    return
  end Subroutine qflux2D


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
