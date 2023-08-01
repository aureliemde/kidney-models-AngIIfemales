!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	 This subroutine yields the non-linear equations to be solved at any
!	 point below the inlet, in order to determine the concentrations, volumes, 
!	 and EP in P, A, B, E, and the lumen M.  
!    It is used to compute values below the entry to the PT.
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
subroutine fcn2PT(n,x,fvec,iflag,numpar,pars,pt,id)


  include 'values.h'
  include 'global.h'
  include 'defs.h'

  type (membrane) :: pt
  integer n,iflag,numpar,id
  double precision x(n),fvec(n),S(NDPT),pars(numpar)
  double precision Ca(NSPT,NC),Vola(NC) !,EPa(NC),pha(NC)
  double precision Cb(NSPT,NC),Volb(NC),EPb(NC),phb(NC)
  double precision PMa,PMb,VolEinit,VolPinit,vol0
  double precision Jva(NC,NC), Jsa(NSPT,NC,NC)
  double precision Jvb(NC,NC), Jsb(NSPT,NC,NC)
  double precision dimL, CPimpref
  double precision dkd(NC),dkh(NC)
  double precision :: dcompl, dtorq  ! compliance and torque effects
  
!	 Ca(K) denotes the known concentrations at Lz (same with pH, Vol, EP)
!	 y is the associated vector, of known variables

!	 Cb(K) denotes the concentrations at Lz+1	(same with pH, Vol, EP)
!	 x is the associated vector, of unknown variables

!---------------------------------------------------------------------72
!	 Assign concentrations, volumes, and potentials 
!	 in P, A, B, E, and the lumen M at Lz
!---------------------------------------------------------------------72

  Ca(:,1)=pars(1:NS)
  Ca(:,2)=pars(1+NS:2*NS)
  Ca(:,5)=pars(1+2*NS:3*NS)
  Ca(:,6)=pars(1+3*NS:4*NS)
  
  Vola(1)=pars(1+4*NSPT)
  Vola(2)=pars(2+4*NSPT)
  Vola(5)=pars(3+4*NSPT)
  
  PMa=pars(4+4*NSPT)

!---------------------------------------------------------------------72
!	 Extract other parameters
!---------------------------------------------------------------------72

  VolEinit = pars(5+4*NSPT)
  VolPinit = pars(6+4*NSPT)
  dimL = pars(7+4*NSPT)
  CPimpref = pars(8+4*NSPT)
  CPbuftot1 = pars(9+4*NSPT)
  dkd(1) = pars(10+4*NSPT) 
  dkd(2) = pars(11+4*NSPT) 
  dkd(5) = pars(12+4*NSPT) 
  dkh(1) = pars(13+4*NSPT) 
  dkh(2) = pars(14+4*NSPT)  
  dkh(5) = pars(15+4*NSPT)  

  dcompl = pars(16+4*NSPT)  ! compliant PT?
  dtorq  = pars(17+4*NSPT)  ! compliant PT?
  vol0   = pars(18+4*NSPT)  ! inflow volume
  
!---------------------------------------------------------------------72
!	 Assign concentrations, volumes, and potentials 
!	 in P, A, B, E, and the lumen M at Lz+1
!---------------------------------------------------------------------72

  Do I = 1,NSPT
     Cb(I,1)=x(1+3*(I-1))
     Cb(I,2)=x(2+3*(I-1))
     Cb(I,5)=x(3+3*(I-1))
  End do
  Cb(:,6) = pt%conc(:,6)
  EPb(6) = pt%ep(6)

  phb(1)=-dlog(Cb(12,1)/1.0d3)/dlog(10.0d0)
  phb(2)=-dlog(Cb(12,2)/1.0d3)/dlog(10.0d0)
  phb(5)=-dlog(Cb(12,5)/1.0d3)/dlog(10.0d0)
  
  Volb(1) = x(1+3*NSPT) 
  Volb(2) = x(2+3*NSPT) 
  Volb(5) = x(3+3*NSPT) 
  EPb(1) = x(4+3*NSPT) 
  EPb(2) = x(5+3*NSPT) 
  EPb(5) = x(6+3*NSPT) 
  PMb = x(7+3*NSPT)

  Call qfluxPT(x,Jvb,Jsb,Cb(:,6),EPb(6),pt,vol0,dcompl,dtorq)
  
!---------------------------------------------------------------------72
!	Initialize  source terms
!---------------------------------------------------------------------72

  Do K = 1, NDPT
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

  Bm = PI*pt%diam
  Am = PI*(pt%diam**2)/4.0d0
  
  RMcompl=torqR*(1.0d0+torqvm*(PMb - PbloodPT))
  Bmcompl=2*PI*RMcompl ! Do not assume a change in perimeter for absorption
  ! Already accounted for via torque effects

  Do I = 1, NSPT
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

  Do I = 1, NSPT
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
  S(1+3*NSPT)=fvolb-fvola+Bm*dimL*sumJvb/NZ

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE SOURCE TERMS FOR VOLUME IN EACH CELLULAR COMPARTMENT
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  S(2+3*NSPT) = Jvb(2,5)+Jvb(2,6)-Jvb(1,2) !volume in P
  S(3+3*NSPT) = Jvb(5,6)-Jvb(1,5)-Jvb(2,5) !volume in E

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Convert to equations of the form F(X) = 0
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	 For I = 1-3, and 9 (non-reacting solutes)
!---------------------------------------------------------------------72
		
  Do M=1,9
     fvec(M)=S(M)
  End Do

  fvec(25)=S(25)
  fvec(26)=S(26)
  fvec(27)=S(27)
  
!---------------------------------------------------------------------72
!	 For glucose, with consumption term linked to Na,K-ATPase activity
!---------------------------------------------------------------------72

  ! flux_rec(3) = FNaK
  Pumpactivity = pt%FNaK/(href*Cref*PI*pt%diam*60/10*1.0d9)
  Gluconsumption = Pumpactivity/(pt%TQ*6)
  fvec(43)=S(43)
  fvec(44)=S(44)+Gluconsumption*0
  fvec(45)=S(45)

!---------------------------------------------------------------------72
!    For Ca2+
!---------------------------------------------------------------------72

  fvec(46) = S(46)
  fvec(47) = S(47)
  fvec(48) = S(48)

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

	 
  fkin1=(dkh(1)*Ca(6,1)-dkd(1)*Ca(5,1))
  fkin2=(dkh(1)*Cb(6,1)-dkd(1)*Cb(5,1))
  fvec(16) = S(16) + Am*dimL*fkin2/NZ/href

  facnd=Vref/href
  fvec(17) = S(17)+Volb(2)*(dkh(2)*Cb(6,2) -dkd(2)*Cb(5,2))*facnd
  fvec(18) = S(18)+max(Volb(5),VolEinit)*(dkh(5)*Cb(6,5) &
       -dkd(5)*Cb(5,5))*facnd

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
!	 For NH3/NH4 with ammoniagenesis within lumen
!---------------------------------------------------------------------72

  if (dcompl < 0) then
     RMtorq = pt%diam/2.0d0 !non-compliant tubule
  else
     RMtorq = torqR*(1.0d0+torqvm*(PMb - PbloodPT))
!     RMtorq = torqR*(1.0d0+0.50*dtanh(torqvm*(PMb - PbloodPT)))
  endif
  factor1 = 8.0*visc*(Volb(1)*Vref)*torqL/(RMtorq**2)
  factor2 = 1.0 + (torqL+torqd)/RMtorq + 0.50*((torqL/RMtorq)**2)
  Torque = factor1*factor2

  if (dtorq < 0) then
     Scaletorq = 1.0 !No torque effect on transporter density
  else
     Scaletorq = 1.0 + TS*pt%scaleT*(Torque/pt%TM0-1.0)
  endif
  Scaletorq = max(Scaletorq,0.001d0)

!    The rate of ammoniagenesis also scales with the torque
  if (id == 0) then
     Qnh4 = pt%qnh4*Scaletorq
  else
     Qnh4 = 0.d0
  endif
  
  fvec(28) = S(28)+S(31)
  fvec(29) = S(29)+S(32)-Qnh4
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

  fvec(1+3*NSPT)=S(1+3*NSPT)
  fvec(2+3*NSPT)=S(2+3*NSPT)
  fvec(3+3*NSPT)=S(3+3*NSPT)

!---------------------------------------------------------------------72
!	 For EP, need to satisfy electroneutrality in epithelial compartments
!	 and zero-current condition in lumen
!---------------------------------------------------------------------72

	
  currM=0.0d0
  do I=1,NSPT
     currM=currM+zval(I)*(Jsb(I,1,2)+Jsb(I,1,5))
  end do
  fvec(4+3*NSPT)=currM

  volPrat=VolPinit/Volb(2)
  CimpP=CPimpref*volPrat
  facP=dexp(dlog(10.0d0)*(phb(2)-pKbuf))
  CbufP=CPbuftot1*volPrat*facP/(facP+1)
  
  elecP=zPimpPT*CimpP-CbufP
  elecE=0.0d0
  do I=1,NSPT
     elecP=elecP+zval(I)*Cb(I,2)
     elecE=elecE+zval(I)*Cb(I,5)
  end do
  fvec(5+3*NSPT)=elecP
  fvec(6+3*NSPT)=elecE

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!		Equation for pressure - Poiseuille law
!		Flow must be multiplied by Vref
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!    For a compliant tubule, the radius and area vary with pressure

  RMcompl=torqR*(1.0d0+torqvm*(PMb - PbloodPT))
  Amcompl=PI*(RMcompl**2)

  if (dcompl < 0) then   ! non-compliant PT
     factor1=8.0d0*PI*visc/(Am**2)
  else
     factor1=8.0d0*PI*visc/(Amcompl**2)
  endif
  
  fvec(7+3*NSPT)=PMb-PMa+factor1*Volb(1)*Vref*dimL/NZ
  

  return
end subroutine fcn2PT
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
