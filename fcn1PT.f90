!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	 This subroutine yields the non-linear equations to be solved at the inlet
!	 in order to determine the concentrations, volumes, and EP 
!	 in P, A, B, and E	   
!    It is used to compute values at the entry of the PT.
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
	   
subroutine fcn1PT(n,x,fvec,iflag,numpar,pars,pt,id)


  include 'values.h'
  include 'global.h'
  include 'defs.h'

  type (membrane) :: pt
  integer n,iflag,numpar,id
  double precision x(NDPT),fvec(n), pars(numpar)
!  double precision dLA(NSPT,NSPT,NC,NC)
  double precision S(NUPT),y(NDPT)
  double precision C(NSPT,NC),Vol(NC),EP(NC),ph(NC),Vol0
  double precision Jvol(NC,NC), Jsol(NSPT,NC,NC)
  double precision VolEinit,VolPinit
  double precision dkd(NC),dkh(NC)

  double precision :: dcompl, dtorq  ! compliance and torque effects
  
!	 Assign concentrations, volumes, and potentials in lumen and blood
  C(:,1) = pt%conc(:,1)
  C(:,6) = pt%conc(:,6)
  ph(1) = -dlog(C(12,1)/1.0d3)/dlog(10.0d0)
  ph(6) = -dlog(C(12,6)/1.0d3)/dlog(10.0d0)

  EP(6) = pt%ep(6)

  VolEinit = pars(1)
  VolPinit = pars(2)
  CPimpref = pars(3)
  CPbuftot1 = pars(4)
  dkd(2) = pars(5)
  dkd(5) = pars(6)
  dkh(2) = pars(7)
  dkh(5) = pars(8)
  dcompl = pars(9)
  dtorq  = pars(10)
  
!	 Assign concentrations, volumes, and potentials in P, A, B, and E

  Do I = 1,NSPT
     C(I,2)=x(1+2*(I-1))
     C(I,5)=x(2+2*(I-1))
  End do
	 
  ph(2) = -dlog(C(12,2)/1.0d3)/dlog(10.0d0)
  ph(5) = -dlog(C(12,5)/1.0d3)/dlog(10.0d0)

  Vol(2)=x(1+2*NSPT)
  Vol(5)=x(2+2*NSPT)
  EP(1)=x(3+2*NSPT)
  EP(2)=x(4+2*NSPT)
  EP(5)=x(5+2*NSPT)
  
!---------------------------------------------------------------------72
!	Initialize  source terms
!---------------------------------------------------------------------72

  Do K = 1, NUPT
     S(K) = 0.0d0
     fvec(K) = 0.0d0
  End Do
	 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE SOURCE TERMS FOR SOLUTE I IN EACH CELLULAR COMPARTMENT
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  Do I = 1,NSPT
     y(1+3*(I-1))=pt%conc(I,1)
     y(2+3*(I-1))=C(I,2)
     y(3+3*(I-1))=C(I,5)
  End do
  Vol0 = pt%vol(1)
  y(1+3*NSPT) = Vol0
  y(2+3*NSPT) = Vol(2)
  y(3+3*NSPT) = Vol(5)
  y(4+3*NSPT) = EP(1)
  y(5+3*NSPT) = EP(2)
  y(6+3*NSPT) = EP(5)
  y(7+3*NSPT) = PMinitPT

  Call qfluxPT(y,Jvol,Jsol,C(:,6),EP(6),pt,Vol0,dcompl,dtorq)

  Do I = 1, NSPT
     S(1+2*(I-1)) = -Jsol(I,1,2)+Jsol(I,2,5)+Jsol(I,2,6) 
     S(2+2*(I-1)) = -Jsol(I,1,5)-Jsol(I,2,5)+Jsol(I,5,6)
  End Do
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE SOURCE TERMS FOR VOLUME IN EACH CELLULAR COMPARTMENT
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  S(1+2*NSPT) = -Jvol(1,2)+Jvol(2,5)+Jvol(2,6)
  S(2+2*NSPT) = -Jvol(1,5)-Jvol(2,5)+Jvol(5,6) 
  
!---------------------------------------------------------------------72
!	Convert to equations of the form F(X) = 0
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	 For I = 1-3, and 9 (non-reacting solutes)
!---------------------------------------------------------------------72
		
  Do M=1,6
     fvec(M)=S(M)
  End Do

  fvec(17)=S(17)
  fvec(18)=S(18)
  
!---------------------------------------------------------------------72
!	 For glucose, with consumption term linked to Na,K-ATPase activity
!---------------------------------------------------------------------72

  Pumpactivity = pt%FNaK/(href*Cref*PI*pt%diam*60/10*1.0d9)
  Gluconsumption = Pumpactivity/(pt%TQ*6)
  fvec(29)=S(29)+Gluconsumption*0 ! Multiply by 0 to not account for consumption
  fvec(30)=S(30)

!---------------------------------------------------------------------72
!    For Ca2+
!---------------------------------------------------------------------72

  fvec(31) = S(31)
  fvec(32) = S(32)

!---------------------------------------------------------------------72
!	 For CO2/HCO3/H2CO3
!---------------------------------------------------------------------72

  facnd=Vref/href

  fvec(7) = S(7)+S(9)+S(11)
  fvec(8) = S(8)+S(10)+S(12)
  fvec(9) = ph(2)-pKHCO3-dlog(abs(C(4,2)/C(5,2)))/dlog(10.0d0)
  fvec(10) = ph(5)-pKHCO3-dlog(abs(C(4,5)/C(5,5)))/dlog(10.0d0)
  fvec(11) = S(11)+Vol(2)*(dkh(2)*C(6,2)-dkd(2)*C(5,2))*facnd
  fvec(12) = S(12)+max(Vol(5),VolEinit)*(dkh(5)*C(6,5)-dkd(5)*C(5,5))*facnd

!---------------------------------------------------------------------72
!	 For HPO4(2-)/H2PO4(-)
!---------------------------------------------------------------------72

  fvec(13) = S(13)+S(15)
  fvec(14) = S(14)+S(16)
  fvec(15) = ph(2)-pKHPO4-dlog(abs(C(7,2)/C(8,2)))/dlog(10.0d0)
  fvec(16) = ph(5)-pKHPO4-dlog(abs(C(7,5)/C(8,5)))/dlog(10.0d0)
  
!---------------------------------------------------------------------72
!	 For NH3/NH4 with ammoniagenesis within lumen
!---------------------------------------------------------------------72

  if (dcompl < 0) then
     RMtorq = pt%diam/2.0d0 !non-compliant tubule
  else
    RMtorq = torqR*(1.0d0+torqvm*(PMinitPT - PbloodPT))
!    RMtorq = torqR*(1.0d0+0.50*dtanh(torqvm*(PMinitPT - PbloodPT)))
  endif
  factor1 = 8.0*visc*(Vol0*Vref)*torqL/(RMtorq**2)
  factor2 = 1.0 + (torqL+torqd)/RMtorq + 0.50*((torqL/RMtorq)**2)
  Torque = factor1*factor2

  if (dtorq < 0) then
     Scaletorq = 1.0 !No torque effect on transporter density
  else
     Scaletorq = 1.0 + TS*pt%scaleT*(Torque/pt%TM0-1.0)
  endif
  Scaletorq = max(Scaletorq,0.001d0) ! To avoid issues with negative values
  
  !    The rate of ammoniagenesis also scales with the torque
  if (id == 0) then
     Qnh4 = pt%qnh4*Scaletorq
  else
     Qnh4 = 0.d0
  endif
  
  fvec(19) = S(19)+S(21)-Qnh4
  fvec(20) = S(20)+S(22)
  fvec(21) = ph(2)-pKNH3-dlog(abs(C(10,2)/C(11,2)))/dlog(10.0d0)
  fvec(22) = ph(5)-pKNH3-dlog(abs(C(10,5)/C(11,5)))/dlog(10.0d0)

!---------------------------------------------------------------------72
!	 For HCO2-/H2CO2
!---------------------------------------------------------------------72

  fvec(25) = S(25)+S(27)
  fvec(26) = S(26)+S(28)
  fvec(27) = ph(2)-pKHCO2-dlog(abs(C(13,2)/C(14,2)))/dlog(10.0d0)
  fvec(28) = ph(5)-pKHCO2-dlog(abs(C(13,5)/C(14,5)))/dlog(10.0d0)


!---------------------------------------------------------------------72
!	 For pH
!---------------------------------------------------------------------72
  
  fvec(23) = S(23)+S(21)-S(13)-S(7)-S(25)	 
  fvec(24) = S(24)+S(22)-S(14)-S(8)-S(26)	 
	
!---------------------------------------------------------------------72
!	 For volume
!---------------------------------------------------------------------72

  fvec(1+2*NSPT)=S(1+2*NSPT)
  fvec(2+2*NSPT)=S(2+2*NSPT)

!---------------------------------------------------------------------72
!	 For EP, need to satisfy electroneutrality
!---------------------------------------------------------------------72

  volPrat=VolPinit/Vol(2)
  CimpP=CPimpref*volPrat
  facP=dexp(dlog(10.0d0)*(ph(2)-pKbuf))
  CbufP=CPbuftot1*volPrat*facP/(facP+1)

  elecP=zPimpPT*CimpP-CbufP
  elecE=0.0d0
  do I=1,NSPT
     elecP=elecP+zval(I)*C(I,2)
     elecE=elecE+zval(I)*C(I,5)
  end do
  fvec(3+2*NSPT)=elecP
  fvec(4+2*NSPT)=elecE
  
!	 Electroneutrality/Zero-current in lumen
	
  currM=0.0d0
  do I=1,NSPT
     currM=currM+zval(I)*(Jsol(I,1,2)+Jsol(I,1,5))
  end do
  fvec(5+2*NSPT)=currM


  return
end subroutine fcn1PT


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
