!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	 This subroutine yields the non-linear equations to be solved at any
!	 point below the inlet, in order to determine the concentrations, volumes, 
!	 and EP in P, A, B, E, and the lumen M.  
!    It is used to compute values below the entry of segments that only have
!    principal cells (except for the PT): mTAL, cTAL, DCT and IMCD.

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
subroutine fcn2b(n,x,fvec,iflag,numpar,pars,tube,id)
  
  include 'values.h'
  include 'global.h'
  include 'defs.h'

  integer n,iflag,numpar,id  
  double precision x(n),fvec(n),S(NDA),pars(numpar)
  type (membrane) :: tube
  double precision Ca(NS,NC),Vola(NC) !,EPa(NC),pha(NC)
  double precision Cb(NS,NC),Volb(NC),EPb(NC),phb(NC)
  double precision PMa,PMb
  double precision Jva(NC,NC), Jsa(NS,NC,NC),Jcasra
  double precision Jvb(NC,NC), Jsb(NS,NC,NC),Jcasrb
  double precision CaBT(NC), Diam
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
     Cb(J,6) = pars(J+4+4*NS)
  End Do
  EPb(6) = pars(5+5*NS)

!---------------------------------------------------------------------72
!	 Extract other parameters
!---------------------------------------------------------------------72
  dimL = pars(6+5*NS)
  CPimpref = pars(7+5*NS)
  Diam   = pars(8+5*NS)

  if (id==3) then
		CPbuffer = CPbuftotA 
   elseif (id==4) then
		CPbuffer = CPbuftotT
   elseif (id==5) then
		CPbuffer= CPbuftotD
   elseif (id==9) then
		CPbuffer= CPbuftotIMC
   else
     print *, "wrong id", id
	 pause
  end if

!---------------------------------------------------------------------72
!	 Assign concentrations, volumes, and potentials 
!	 in P, A, B, E, and the lumen M at Lz+1
!---------------------------------------------------------------------72

  Do I = 1,NS2
     Cb(I,1)=x(1+3*(I-1))
     Cb(I,2)=x(2+3*(I-1))
     Cb(I,5)=x(3+3*(I-1))
  End do

  phb(1)=-dlog(Cb(12,1)/1.0d3)/dlog(10.0d0)
  phb(2)=-dlog(Cb(12,2)/1.0d3)/dlog(10.0d0)
  phb(5)=-dlog(Cb(12,5)/1.0d3)/dlog(10.0d0)

  Volb(1) = x(1+3*NS2) 
  Volb(2) = x(2+3*NS2) 
  Volb(5) = x(3+3*NS2) 
  EPb(1) = x(4+3*NS2) 
  EPb(2) = x(5+3*NS2) 
  EPb(5) = x(6+3*NS2) 
  PMb = x(7+3*NS2)

  if (id .eq. 3) then !mTAL
     Call qflux2A(x,Jvb,Jsb,tube%conc(:,6),tube%ep(6),tube)
  elseif (id .eq. 4) then ! cTAL
     Call qflux2T(x,Jvb,Jsb,tube%conc(:,6),tube%ep(6),tube)
  elseif (id .eq. 5) then ! DCT
     Call qflux2D(x,Jvb,Jsb,tube%conc(:,6),tube%ep(6),tube)
  elseif (id .eq. 9) then ! IMCD
     Call qflux2IMC(x,Jvb,Jsb,tube%conc(:,6),tube%ep(6),tube)
  end if

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

  if (id==9) then
     Bm = Bm*tube%coalesce
     Am = Am*tube%coalesce
  end if
  
  Do I = 1, NS2
     sumJsb=(Jsb(I,1,2)+Jsb(I,1,5))
     fsola=Vola(1)*Ca(I,1)*Vref/href
     fsolb=Volb(1)*Cb(I,1)*Vref/href
     S(1+3*(I-1))=fsolb-fsola+Bm*dimL*sumJsb/NZ
  End do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE SOURCE TERMS FOR SOLUTE I IN EACH CELLULAR COMPARTMENT
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  Do I = 1, NS2
     S(2+3*(I-1)) = Jsb(I,2,5)+Jsb(I,2,6)-Jsb(I,1,2) !solute I in P
     S(3+3*(I-1)) = Jsb(I,5,6)-Jsb(I,1,5)-Jsb(I,2,5) !solute I in E
  End Do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE SOURCE TERMS FOR VOLUME FLOW IN THE LUMEN
!	The non-dimensional volume flux needs to be multiplied by 
!	Vref/(Pfref.Vwbar.Cref)
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  fvmult = (Pfref*Vwbar*Cref)
  sumJvb=(Jvb(1,2)+Jvb(1,5))*fvmult
  fvola=Vola(1)*Vref
  fvolb=Volb(1)*Vref
  S(1+3*NS2)=fvolb-fvola+Bm*dimL*sumJvb/NZ

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE SOURCE TERMS FOR VOLUME IN EACH CELLULAR COMPARTMENT
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  S(2+3*NS2) = Jvb(2,5)+Jvb(2,6)-Jvb(1,2) !volume in P
  S(3+3*NS2) = Jvb(5,6)-Jvb(1,5)-Jvb(2,5) !volume in E

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Convert to equations of the form F(X) = 0
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	 For I = 1-3 & 9 & 13 & 16 (non-reacting solutes)
!---------------------------------------------------------------------72
		
  Do M=1,9 ! Na+, K+, Cl-
     fvec(M)=S(M)
  End Do
  
  fvec(25)=S(25) ! Urea - solute 9; index = 1+3*(9-1)
  fvec(26)=S(26)
  fvec(27)=S(27)

  fvec(43)=S(43) ! Glucose - solute 15; index = 1+3*(15-1)
  fvec(44)=S(44)
  fvec(45)=S(45)

  fvec(46)=S(46) ! Calcium - solute 16; index = 1+3*(16-1)
  fvec(47)=S(47)
  fvec(48)=S(48)


!---------------------------------------------------------------------72
!	 For CO2/HCO3/H2CO3
!	 The dimensional factor for the kinetic term in the lumen is 
!	 Am/Bm = Pi R^2 / 2 PI R = R / 2 = D / 4
!---------------------------------------------------------------------72

  fvec(10) = S(10)+S(13)+S(16)
  fvec(11) = S(11)+S(14)+S(17)
  fvec(12) = S(12)+S(15)+S(18)
  fvec(13) = phb(1)-pKHCO3-dlog(Cb(4,1)/Cb(5,1))/dlog(10.0d0)
  fvec(14) = phb(2)-pKHCO3-dlog(Cb(4,2)/Cb(5,2))/dlog(10.0d0)
  fvec(15) = phb(5)-pKHCO3-dlog(Cb(4,5)/Cb(5,5))/dlog(10.0d0)
	 
  fkin1=(tube%dkh(1)*Ca(6,1)-tube%dkd(1)*Ca(5,1))
  fkin2=(tube%dkh(1)*Cb(6,1)-tube%dkd(1)*Cb(5,1))
  fvec(16) = S(16) + Am*dimL*fkin2/NZ/href
  
  facnd=Vref/href
  fvec(17) = S(17)+Volb(2)*(tube%dkh(2)*Cb(6,2) -tube%dkd(2)*Cb(5,2))*facnd
  fvec(18) = S(18)+max(Volb(5),tube%volEinit)*(tube%dkh(5)*Cb(6,5) &
       -tube%dkd(5)*Cb(5,5))*facnd

!---------------------------------------------------------------------72
!	 For HPO4(2-)/H2PO4(-)
!---------------------------------------------------------------------72

  fvec(19) = S(19)+S(22)
  fvec(20) = S(20)+S(23)
  fvec(21) = S(21)+S(24)
  fvec(22) = phb(1)-pKHPO4-dlog(Cb(7,1)/Cb(8,1))/dlog(10.0d0)
  fvec(23) = phb(2)-pKHPO4-dlog(Cb(7,2)/Cb(8,2))/dlog(10.0d0)
  fvec(24) = phb(5)-pKHPO4-dlog(Cb(7,5)/Cb(8,5))/dlog(10.0d0)

!---------------------------------------------------------------------72
!	 For NH3/NH4
!---------------------------------------------------------------------72
  
  fvec(28) = S(28)+S(31)
  fvec(29) = S(29)+S(32)
  fvec(30) = S(30)+S(33)
  fvec(31) = phb(1)-pKNH3-dlog(Cb(10,1)/Cb(11,1))/dlog(10.0d0)
  fvec(32) = phb(2)-pKNH3-dlog(Cb(10,2)/Cb(11,2))/dlog(10.0d0)
  fvec(33) = phb(5)-pKNH3-dlog(Cb(10,5)/Cb(11,5))/dlog(10.0d0)


!---------------------------------------------------------------------72
!	 For HCO2-/H2CO2
!---------------------------------------------------------------------72

  fvec(37) = S(37)+S(40)
  fvec(38) = S(38)+S(41)
  fvec(39) = S(39)+S(42)
  fvec(40) = phb(1)-pKHCO2-dlog(abs(Cb(13,1)/Cb(14,1)))/dlog(10.0d0)
  fvec(41) = phb(2)-pKHCO2-dlog(abs(Cb(13,2)/Cb(14,2)))/dlog(10.0d0)
  fvec(42) = phb(5)-pKHCO2-dlog(abs(Cb(13,5)/Cb(14,5)))/dlog(10.0d0)


!---------------------------------------------------------------------72
!	 For pH
!---------------------------------------------------------------------72

  fvec(34) = S(34)+S(31)-S(19)-S(10)-S(37)	 !pH in M
  fvec(35) = S(35)+S(32)-S(20)-S(11)-S(38)	 !pH in P
  fvec(36) = S(36)+S(33)-S(21)-S(12)-S(39)	 !pH in E


!---------------------------------------------------------------------72
!	 For volume
!---------------------------------------------------------------------72

  fvec(1+3*NS2)=S(1+3*NS2)
  fvec(2+3*NS2)=S(2+3*NS2)
  fvec(3+3*NS2)=S(3+3*NS2)

!---------------------------------------------------------------------72
!	 For EP, need to satisfy electroneutrality in epithelial compartments
!	 and zero-current condition in lumen
!---------------------------------------------------------------------72
	
  currM=0.0d0
  do I=1,NS
     currM=currM+zval(I)*(Jsb(I,1,2)+Jsb(I,1,5))
  end do
  fvec(4+3*NS2)=currM
  
  volPrat=tube%volPinit/Volb(2)
  CimpP=CPimpref*volPrat
  facP=dexp(dlog(10.0d0)*(phb(2)-pKbuf))
  CbufP=CPbuffer*volPrat*facP/(facP+1)

  elecP=zPimp*CimpP-CbufP
  elecE=0.0d0
  do I=1,NS2
     elecP=elecP+zval(I)*Cb(I,2)
     elecE=elecE+zval(I)*Cb(I,5)
  end do
  fvec(5+3*NS2)=elecP
  fvec(6+3*NS2)=elecE

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!		Equation for pressure - Poiseuille's law
!		Flow must be multiplied by Vref
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!    Am is the luminal cross-sectional area, given by PI*(DiamPT**2)/4.0d0
  ratio=8.0d0*PI*visc/(Am**2)
  if (id==9) ratio = ratio*tube%coalesce*2  ! correct extra coalescence factor in Am, add extra resistance from merging

  fvec(7+3*NS2)=PMb-PMa+ratio*Volb(1)*Vref*dimL/NZ


  return
end subroutine fcn2b
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
