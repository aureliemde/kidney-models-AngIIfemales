!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	 This subroutine yields the non-linear equations to be solved at any
!	 point below the inlet, in order to determine the concentrations, volumes, 
!	 and EP in P, A, B, E, and the lumen M.  
!    It is used to compute values below the entry to the CNT and the CCD.
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
subroutine fcn2C(n,x,fvec,iflag,numpar,pars,tube,id)

  include 'values.h'
  include 'global.h'
  include 'defs.h'

  integer n,iflag,nummpar
  type(membrane) :: tube
  double precision x(n),fvec(n),S(n),pars(numpar)
  double precision y(n) 
  double precision Ca(NS,NC),Vola(NC),EPa(NC),pha(NC)
  double precision Cb(NS,NC),Volb(NC),EPb(NC),phb(NC)
  double precision PMa,PMb
  double precision Jva(NC,NC), Jsa(NS,NC,NC),Jcasra
  double precision Jvb(NC,NC), Jsb(NS,NC,NC),Jcasrb
  double precision CaBT(NC)
  
!	 Ca(K) denotes the known concentrations at Lz (same with pH, Vol, EP)
!	 y is the associated vector, of known variables

!	 Cb(K) denotes the concentrations at Lz+1	(same with pH, Vol, EP)
!	 x is the associated vector, of unknown variables

!---------------------------------------------------------------------72
!	 Assign concentrations, volumes, and potentials in peritubular solution
!---------------------------------------------------------------------72
  Do J = 1,NS
     Ca(J,6) = tube%conc(J,6)
     Cb(J,6) = tube%conc(J,6)
  End Do
  EPa(6) = tube%ep(6)
  EPb(6) = tube%ep(6)

!---------------------------------------------------------------------72
!	 Assign concentrations, volumes, and potentials 
!	 in P, A, B, E, and the lumen M at Lz
!---------------------------------------------------------------------72

  Do K = 1,NC-1
     Do I = 1,NS2
        Ca(I,K)=pars(I+(K-1)*NS) 
        y(K+5*(I-1))=Ca(I,K)     
     End do
     pha(K)=-dlog(Ca(12,K)/1.0d3)/dlog(10.0d0)
     Vola(K)=pars(K+(NC-1)*NS) !where Vol(1) is the non-dim volume flux
     y(K+5*NS2)=Vola(K)        
     EPa(K)=pars(K+(NC-1)*(NS+1)) 
     y(K+5+5*NS2)=EPa(K)          
  End do
  PMa=pars(1+(NC-1)*(NS+2))   
  y(11+5*NS2)=PMa             
  
  if (id==6) then
     Call qflux2C(y,Jva,Jsa,tube)
  elseif (id==7) then
     Call qflux2CCD(y,Jva,Jsa,tube)
  else
     print *,"wrong id",id
     stop
  end if

!---------------------------------------------------------------------72
!	 Assign concentrations, volumes, and potentials 
!	 in P, A, B, E, and the lumen M at Lz+1
!---------------------------------------------------------------------72

  Do I = 1,NS2
     Cb(I,1)=x(1+5*(I-1))
     Cb(I,2)=x(2+5*(I-1))
     Cb(I,3)=x(3+5*(I-1))
     Cb(I,4)=x(4+5*(I-1))
     Cb(I,5)=x(5+5*(I-1))
  End do
  Do K = 1,NC
     phb(K)=-dlog(Cb(12,K)/1.0d3)/dlog(10.0d0)
     Volb(K)=x(K+5*NS2)
     EPb(K)=x(K+5+5*NS2)
  End do
  PMb=x(11+5*NS2) !Hydraulic pressure in lumen
  
  if (id==6) then
     Call qflux2C(x,Jvb,Jsb,tube)
  elseif (id==7) then
     Call qflux2CCD(x,Jvb,Jsb,tube)
  else
     print *,"wrong id",id
     stop
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

  coalesce = tube%coalesce

  Bm = PI*tube%diam*coalesce
  Am = PI*(tube%diam**2)/4.0d0*coalesce

  Do I = 1, NS2
     sumJsa=(Jsa(I,1,2)+Jsa(I,1,3)+Jsa(I,1,4)+Jsa(I,1,5))
     sumJsb=(Jsb(I,1,2)+Jsb(I,1,3)+Jsb(I,1,4)+Jsb(I,1,5))
     fsola=Vola(1)*Ca(I,1)*Vref/href
     fsolb=Volb(1)*Cb(I,1)*Vref/href
     S(1+5*(I-1))=fsolb-fsola+Bm*tube%dimL*sumJsb/NZ
  End do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE SOURCE TERMS FOR SOLUTE I IN EACH CELLULAR COMPARTMENT
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  Do I = 1, NS2
     S(2+5*(I-1)) = Jsb(I,2,5)+Jsb(I,2,6)-Jsb(I,1,2) !solute I in P
     S(3+5*(I-1)) = Jsb(I,3,5)+Jsb(I,3,6)-Jsb(I,1,3) !solute I in A
     S(4+5*(I-1)) = Jsb(I,4,5)+Jsb(I,4,6)-Jsb(I,1,4) !solute I in B
     S(5+5*(I-1)) = -Jsb(I,1,5)-Jsb(I,2,5)-Jsb(I,3,5)-Jsb(I,4,5)+Jsb(I,5,6) !solute I in E
  End Do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE SOURCE TERMS FOR VOLUME FLOW IN THE LUMEN
!	The non-dimensional volume flux needs to be multiplied by 
!	Vref/(Pfref.Vwbar.Cref)
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  fvmult = (Pfref*Vwbar*Cref)
  sumJva=(Jva(1,2)+Jva(1,3)+Jva(1,4)+Jva(1,5))*fvmult
  sumJvb=(Jvb(1,2)+Jvb(1,3)+Jvb(1,4)+Jvb(1,5))*fvmult
  fvola=Vola(1)*Vref
  fvolb=Volb(1)*Vref
  S(1+5*NS2)=fvolb-fvola+Bm*tube%dimL*sumJvb/NZ

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE SOURCE TERMS FOR VOLUME IN EACH CELLULAR COMPARTMENT
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  S(2+5*NS2) = Jvb(2,5)+Jvb(2,6)-Jvb(1,2) !volume in P
  S(3+5*NS2) = Jvb(3,5)+Jvb(3,6)-Jvb(1,3) !in A
  S(4+5*NS2) = Jvb(4,5)+Jvb(4,6)-Jvb(1,4) !in B
  S(5+5*NS2) = -Jvb(1,5)-Jvb(2,5)-Jvb(3,5)-Jvb(4,5)+Jvb(5,6) !in E
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Convert to equations of the form F(X) = 0
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	 For I = 1-3 and 9 (non-reacting solutes)
!---------------------------------------------------------------------72
		
  Do M=1,15
     fvec(M)=S(M)
  End Do
  
  fvec(41)=S(41) ! Urea - solute 9; index = 1+5*(9-1)
  fvec(42)=S(42)
  fvec(43)=S(43)
  fvec(44)=S(44)
  fvec(45)=S(45)

  fvec(71)=S(71) ! Glucose - solute 15; index = 1+5*(15-1)
  fvec(72)=S(72)
  fvec(73)=S(73)
  fvec(74)=S(74)
  fvec(75)=S(75)

  fvec(76)=S(76) ! Calcium - solute 16; index = 1+5*(16-1)
  fvec(77)=S(77)
  fvec(78)=S(78)
  fvec(79)=S(79)
  fvec(80)=S(80)

! HCO2-/H2CO2 equations modified below
  Do M=13,16
     fvec(1+5*(M-1))=S(1+5*(M-1))
     fvec(2+5*(M-1))=S(2+5*(M-1))
     fvec(3+5*(M-1))=S(3+5*(M-1))
     fvec(4+5*(M-1))=S(4+5*(M-1))
     fvec(5+5*(M-1))=S(5+5*(M-1))
  End Do

!---------------------------------------------------------------------72
!	 For CO2/HCO3/H2CO3
!	 The dimensional factor for the kinetic term in the lumen is 
!	 Am/Bm = Pi R^2 / 2 PI R = R / 2 = D / 4
!---------------------------------------------------------------------72

  fvec(16) = S(16)+S(21)+S(26)
  fvec(17) = S(17)+S(22)+S(27)
  fvec(18) = S(18)+S(23)+S(28)
  fvec(19) = S(19)+S(24)+S(29)
  fvec(20) = S(20)+S(25)+S(30)
  fvec(21) = phb(1)-pKHCO3-dlog(Cb(4,1)/Cb(5,1))/dlog(10.0d0)
  fvec(22) = phb(2)-pKHCO3-dlog(Cb(4,2)/Cb(5,2))/dlog(10.0d0)
  fvec(23) = phb(3)-pKHCO3-dlog(Cb(4,3)/Cb(5,3))/dlog(10.0d0)
  fvec(24) = phb(4)-pKHCO3-dlog(Cb(4,4)/Cb(5,4))/dlog(10.0d0)
  fvec(25) = phb(5)-pKHCO3-dlog(Cb(4,5)/Cb(5,5))/dlog(10.0d0)
  
  fkin1=(tube%dkh(1)*Ca(6,1)-tube%dkd(1)*Ca(5,1))
  fkin2=(tube%dkh(1)*Cb(6,1)-tube%dkd(1)*Cb(5,1))
  fvec(26) = S(26) + Am*tube%dimL*fkin2/NZ/href
  
  facnd=Vref/href
  fvec(27) = S(27)+Volb(2)*(tube%dkh(2)*Cb(6,2)-tube%dkd(2)*Cb(5,2))*facnd
  fvec(28) = S(28)+Volb(3)*(tube%dkh(3)*Cb(6,3)-tube%dkd(3)*Cb(5,3))*facnd
  fvec(29) = S(29)+Volb(4)*(tube%dkh(4)*Cb(6,4)-tube%dkd(4)*Cb(5,4))*facnd
  fvec(30) = S(30)+max(Volb(5),tube%volEinit)*(tube%dkh(5)*Cb(6,5)-tube%dkd(5)*Cb(5,5))*facnd

!---------------------------------------------------------------------72
!	 For HPO4(2-)/H2PO4(-)
!---------------------------------------------------------------------72

  fvec(31) = S(31)+S(36)
  fvec(32) = S(32)+S(37)
  fvec(33) = S(33)+S(38)
  fvec(34) = S(34)+S(39)
  fvec(35) = S(35)+S(40)
  fvec(36) = phb(1)-pKHPO4-dlog(Cb(7,1)/Cb(8,1))/dlog(10.0d0)
  fvec(37) = phb(2)-pKHPO4-dlog(Cb(7,2)/Cb(8,2))/dlog(10.0d0)
  fvec(38) = phb(3)-pKHPO4-dlog(Cb(7,3)/Cb(8,3))/dlog(10.0d0)
  fvec(39) = phb(4)-pKHPO4-dlog(Cb(7,4)/Cb(8,4))/dlog(10.0d0)
  fvec(40) = phb(5)-pKHPO4-dlog(Cb(7,5)/Cb(8,5))/dlog(10.0d0)

!---------------------------------------------------------------------72
!	 For NH3/NH4
!---------------------------------------------------------------------72

  fvec(46) = S(46)+S(51)
  fvec(47) = S(47)+S(52)
  fvec(48) = S(48)+S(53)
  fvec(49) = S(49)+S(54)
  fvec(50) = S(50)+S(55)
  fvec(51) = phb(1)-pKNH3-dlog(Cb(10,1)/Cb(11,1))/dlog(10.0d0)
  fvec(52) = phb(2)-pKNH3-dlog(Cb(10,2)/Cb(11,2))/dlog(10.0d0)
  fvec(53) = phb(3)-pKNH3-dlog(Cb(10,3)/Cb(11,3))/dlog(10.0d0)
  fvec(54) = phb(4)-pKNH3-dlog(Cb(10,4)/Cb(11,4))/dlog(10.0d0)
  fvec(55) = phb(5)-pKNH3-dlog(Cb(10,5)/Cb(11,5))/dlog(10.0d0)

!---------------------------------------------------------------------72
!	 For HCO2-/H2CO2
!---------------------------------------------------------------------72

  fvec(61) = S(61)+S(66)
  fvec(62) = S(62)+S(67)
  fvec(63) = S(63)+S(68)
  fvec(64) = S(64)+S(69)
  fvec(65) = S(65)+S(70)

  fvec(66) = phb(1)-pKHCO2-dlog(abs(Cb(13,1)/Cb(14,1)))/dlog(10.0d0)
  fvec(67) = phb(2)-pKHCO2-dlog(abs(Cb(13,2)/Cb(14,2)))/dlog(10.0d0)
  fvec(68) = phb(3)-pKHCO2-dlog(abs(Cb(13,3)/Cb(14,3)))/dlog(10.0d0)
  fvec(69) = phb(4)-pKHCO2-dlog(abs(Cb(13,4)/Cb(14,4)))/dlog(10.0d0)
  fvec(70) = phb(5)-pKHCO2-dlog(abs(Cb(13,5)/Cb(14,5)))/dlog(10.0d0)


!---------------------------------------------------------------------72
!	 For pH
!---------------------------------------------------------------------72

  fvec(56) = S(56)+S(51)-S(16)-S(31)-S(61)	 !pH in M
  fvec(57) = S(57)+S(52)-S(17)-S(32)-S(62)	 !pH in P
  fvec(58) = S(58)+S(53)-S(18)-S(33)-S(63)	 !pH in A
  fvec(59) = S(59)+S(54)-S(19)-S(34)-S(64)	 !pH in B
  fvec(60) = S(60)+S(55)-S(20)-S(35)-S(65)	 !pH in E

!---------------------------------------------------------------------72
!	 For volume
!---------------------------------------------------------------------72

  fvec(1+5*NS2)=S(1+5*NS2)
  fvec(2+5*NS2)=S(2+5*NS2)
  fvec(3+5*NS2)=S(3+5*NS2)
  fvec(4+5*NS2)=S(4+5*NS2)
  fvec(5+5*NS2)=S(5+5*NS2)
  
!---------------------------------------------------------------------72
!	 For EP, need to satisfy electroneutrality in epithelial compartments
!	 and zero-current condition in lumen
!---------------------------------------------------------------------72
	
  currM=0.0d0
  do I=1,NS2
     do K=2,5
        currM=currM+zval(I)*Jsb(I,1,K)
     end do
  end do
  fvec(6+5*NS2)=currM
  
  volPrat=tube%volPinit/Volb(2)
  volArat=tube%volAinit/Volb(3)
  volBrat=tube%volBinit/Volb(4)
  
  CimpP=tube%cPimpref*volPrat
  CimpA=tube%cAimpref*VolArat
  CimpB=tube%cBimpref*VolBrat
  
  facP=dexp(dlog(10.0d0)*(phb(2)-pKbuf))
  facA=dexp(dlog(10.0d0)*(phb(3)-pKbuf))
  facB=dexp(dlog(10.0d0)*(phb(4)-pKbuf))
  CbufP=tube%cPbuftot*volPrat*facP/(facP+1)
  CbufA=tube%cAbuftot*volArat*facA/(facA+1)
  CbufB=tube%cBbuftot*volBrat*facB/(facB+1)
  
  elecP=tube%zPimp*CimpP-CbufP
  elecA=tube%zAimp*CimpA-CbufA
  elecB=tube%zBimp*CimpB-CbufB
  elecE=0.0d0
  do I=1,NS2
     elecP=elecP+zval(I)*Cb(I,2)
     elecA=elecA+zval(I)*Cb(I,3)
     elecB=elecB+zval(I)*Cb(I,4)
     elecE=elecE+zval(I)*Cb(I,5)
  end do
  fvec(7+5*NS2)=elecP
  fvec(8+5*NS2)=elecA
  fvec(9+5*NS2)=elecB
  fvec(10+5*NS2)=elecE

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!		Equation for pressure
!		Flow must be multiplied by Vref
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  ratio=8.0d0*PI*visc/(Am**2)
  if (id==6) then
     ratio=ratio*coalesce*2 ! correct extra coalescence factor in Am, add extra resistance from merging
  elseif (id==7) then
     ratio=ratio*coalesce ! correct extra coalescence factor in Am
  endif
  ! ratio=ratio*coalesce*2
  fvec(11+5*NS2)=PMb-PMa+ratio*Volb(1)*Vref*tube%dimL/NZ


  
  return
end subroutine fcn2C
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
