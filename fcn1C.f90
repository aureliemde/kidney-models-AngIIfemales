!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	 This subroutine yields the non-linear equations to be solved at the inlet
!	 in order to determine the concentrations, volumes, and EP 
!	 in P, A, B, and E	   
!    It is used to compute values at the entry of the CNT.
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
	   
subroutine fcn1C(n,x,fvec,iflag,numpar,pars,tube,id)

  include 'values.h'
  include 'global.h'
  include 'defs.h'

  integer n,iflag,numpar,id
  type(membrane) :: tube
  double precision x(n),fvec(n),S(n),pars(numpar)
  double precision C(NS,NC),Vol(NC),EP(NC),ph(NC)
  double precision Jvol(NC,NC), Jsol(NS,NC,NC),Jcasr
  double precision CaBT(NC)
  double precision :: y(NDC)  ! variable needed to use qflux2C

  !	 Assign concentrations, volumes, and potentials in lumen and blood
  Do J = 1, NS
     C(J,1) = tube%conc(J,1)
     C(J,6) = tube%conc(J,6)
  End Do
  EP(6) = tube%ep(6)

  !	 Assign concentrations, volumes, and potentials in P, A, B, and E
  
  Do I = 1,NS2
     C(I,2)=x(1+4*(I-1))
     C(I,3)=x(2+4*(I-1))
     C(I,4)=x(3+4*(I-1))
     C(I,5)=x(4+4*(I-1))
  End do
  
  Do K = 2,5
     ph(K) = -dlog(C(12,K)/1.0d3)/dlog(10.0d0)
     Vol(K)=x(K-1+4*NS2)
     EP(K)=x(K+3+4*NS2)
  End do
  
  EP(1)=x(9+4*NS2)

  ! move x to y, add more varaibles, so we can use qflux2C
  do I = 1, NS
     y(1+5*(I-1)) = C(I,1)
     y(2+5*(I-1)) = C(I,2)
     y(3+5*(I-1)) = C(I,3)
     y(4+5*(I-1)) = C(I,4)
     y(5+5*(I-1)) = C(I,5)
  end do
  y(1+5*NS) = tube%vol(1)
  y(2+5*NS) = Vol(2)
  y(3+5*NS) = Vol(3)
  y(4+5*NS) = Vol(4)
  y(5+5*NS) = Vol(5)
  do K = 1, NC-1
     y(K+5+5*NS) = EP(K)
  end do
  y(11+5*NS) = tube%pres
  
!---------------------------------------------------------------------72
!	Initialize  source terms
!---------------------------------------------------------------------72
  
  Do K = 1, NUC
     S(K) = 0.0d0
     fvec(K) = 0.0d0
  End Do
	 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE SOURCE TERMS FOR SOLUTE I IN EACH CELLULAR COMPARTMENT
!---------------------------------------------------------------------72
  !---------------------------------------------------------------------72
  if (id==5) then
     Call qflux2C(y,Jvol,Jsol,tube)
  elseif (id==6) then
     Call qflux2CCD(y,Jvol,Jsol,tube)
  else
     print*,"wrong id!",id
     stop
  end if
  
  Do I = 1, NS2
     S(1+4*(I-1)) = Jsol(I,2,5)+Jsol(I,2,6)-Jsol(I,1,2) !solute I in P
     S(2+4*(I-1)) = Jsol(I,3,5)+Jsol(I,3,6)-Jsol(I,1,3) !solute I in A
     S(3+4*(I-1)) = Jsol(I,4,5)+Jsol(I,4,6)-Jsol(I,1,4) !solute I in B
     S(4+4*(I-1)) = -Jsol(I,1,5)-Jsol(I,2,5)-Jsol(I,3,5)-Jsol(I,4,5)+Jsol(I,5,6) !solute I in E
  End Do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE SOURCE TERMS FOR VOLUME IN EACH CELLULAR COMPARTMENT
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  S(1+4*NS2) = Jvol(2,5)+Jvol(2,6)-Jvol(1,2) !volume in P
  S(2+4*NS2) = Jvol(3,5)+Jvol(3,6)-Jvol(1,3) !in A
  S(3+4*NS2) = Jvol(4,5)+Jvol(4,6)-Jvol(1,4) !in B
  S(4+4*NS2) = -Jvol(1,5)-Jvol(2,5)-Jvol(3,5)-Jvol(4,5)+Jvol(5,6) !in E
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	Convert to equations of the form F(X) = 0
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	 For I = 1-3 and 9 (non-reacting solutes)
!---------------------------------------------------------------------72
		
  Do M=1,12
     fvec(M)=S(M)
  End Do

  fvec(33)=S(33)
  fvec(34)=S(34)
  fvec(35)=S(35)
  fvec(36)=S(36)

! HCO2-/H2CO2 equations modified below
  Do M = 13,15
     fvec(1+4*(M-1))=S(1+4*(M-1))
     fvec(2+4*(M-1))=S(2+4*(M-1))
     fvec(3+4*(M-1))=S(3+4*(M-1))
     fvec(4+4*(M-1))=S(4+4*(M-1))
  End Do

!---------------------------------------------------------------------72
!	 For CO2/HCO3/H2CO3
!---------------------------------------------------------------------72

  facnd=Vref/href

  fvec(13) = S(13)+S(17)+S(21)
  fvec(14) = S(14)+S(18)+S(22)
  fvec(15) = S(15)+S(19)+S(23)
  fvec(16) = S(16)+S(20)+S(24)
  fvec(17) = ph(2)-pKHCO3-dlog(abs(C(4,2)/C(5,2)))/dlog(10.0d0)
  fvec(18) = ph(3)-pKHCO3-dlog(abs(C(4,3)/C(5,3)))/dlog(10.0d0)
  fvec(19) = ph(4)-pKHCO3-dlog(abs(C(4,4)/C(5,4)))/dlog(10.0d0)
  fvec(20) = ph(5)-pKHCO3-dlog(abs(C(4,5)/C(5,5)))/dlog(10.0d0)
  fvec(21) = S(21)+Vol(2)*(tube%dkh(2)*C(6,2)-tube%dkd(2)*C(5,2))*facnd
  fvec(22) = S(22)+Vol(3)*(tube%dkh(3)*C(6,3)-tube%dkd(3)*C(5,3))*facnd
  fvec(23) = S(23)+Vol(4)*(tube%dkh(4)*C(6,4)-tube%dkd(4)*C(5,4))*facnd
  fvec(24) = S(24)+max(Vol(5),tube%volEinit)*(tube%dkh(5)*C(6,5)-tube%dkd(5)*C(5,5))*facnd

!---------------------------------------------------------------------72
!	 For HPO4(2-)/H2PO4(-)
!---------------------------------------------------------------------72
  
  fvec(25) = S(25)+S(29)
  fvec(26) = S(26)+S(30)
  fvec(27) = S(27)+S(31)
  fvec(28) = S(28)+S(32)
  fvec(29) = ph(2)-pKHPO4-dlog(abs(C(7,2)/C(8,2)))/dlog(10.0d0)
  fvec(30) = ph(3)-pKHPO4-dlog(abs(C(7,3)/C(8,3)))/dlog(10.0d0)
  fvec(31) = ph(4)-pKHPO4-dlog(abs(C(7,4)/C(8,4)))/dlog(10.0d0)
  fvec(32) = ph(5)-pKHPO4-dlog(abs(C(7,5)/C(8,5)))/dlog(10.0d0)

!---------------------------------------------------------------------72
!	 For NH3/NH4
!---------------------------------------------------------------------72

  fvec(37) = S(37)+S(41)
  fvec(38) = S(38)+S(42)
  fvec(39) = S(39)+S(43)
  fvec(40) = S(40)+S(44)
  fvec(41) = ph(2)-pKNH3-dlog(abs(C(10,2)/C(11,2)))/dlog(10.0d0)
  fvec(42) = ph(3)-pKNH3-dlog(abs(C(10,3)/C(11,3)))/dlog(10.0d0)
  fvec(43) = ph(4)-pKNH3-dlog(abs(C(10,4)/C(11,4)))/dlog(10.0d0)
  fvec(44) = ph(5)-pKNH3-dlog(abs(C(10,5)/C(11,5)))/dlog(10.0d0)

!---------------------------------------------------------------------72
!	 For HCO2-/H2CO2
!---------------------------------------------------------------------72

  fvec(49) = S(49)+S(53)
  fvec(50) = S(50)+S(54)
  fvec(51) = S(51)+S(55)
  fvec(52) = S(52)+S(56)
  fvec(53) = ph(2)-pKHCO2-dlog(abs(C(13,2)/C(14,2)))/dlog(10.0d0)
  fvec(54) = ph(3)-pKHCO2-dlog(abs(C(13,3)/C(14,3)))/dlog(10.0d0)
  fvec(55) = ph(4)-pKHCO2-dlog(abs(C(13,4)/C(14,4)))/dlog(10.0d0)
  fvec(56) = ph(5)-pKHCO2-dlog(abs(C(13,5)/C(14,5)))/dlog(10.0d0)

!---------------------------------------------------------------------72
!	 For pH
!---------------------------------------------------------------------72

  
  fvec(45) = S(45)+S(41)-S(13)-S(25)-S(49)	 !pH in P
  fvec(46) = S(46)+S(42)-S(14)-S(26)-S(50)	 !pH in A
  fvec(47) = S(47)+S(43)-S(15)-S(27)-S(51)	 !pH in B
  fvec(48) = S(48)+S(44)-S(16)-S(28)-S(52)	 !pH in E


!---------------------------------------------------------------------72
!	 For volume
!---------------------------------------------------------------------72

  fvec(1+4*NS2)=S(1+4*NS2)
  fvec(2+4*NS2)=S(2+4*NS2)
  fvec(3+4*NS2)=S(3+4*NS2)
  fvec(4+4*NS2)=S(4+4*NS2)

!---------------------------------------------------------------------72
!	 For EP, need to satisfy electroneutrality
!---------------------------------------------------------------------72

  volPrat=tube%volPinit/Vol(2)
  volArat=tube%volAinit/Vol(3)
  volBrat=tube%volBinit/Vol(4)

  CimpP=tube%cPimpref*volPrat
  CimpA=tube%cAimpref*VolArat
  CimpB=tube%cBimpref*VolBrat
  
  facP=dexp(dlog(10.0d0)*(ph(2)-pKbuf))
  facA=dexp(dlog(10.0d0)*(ph(3)-pKbuf))
  facB=dexp(dlog(10.0d0)*(ph(4)-pKbuf))
  CbufP=tube%cPbuftot*volPrat*facP/(facP+1)
  CbufA=tube%cAbuftot*volArat*facA/(facA+1)
  CbufB=tube%cBbuftot*volBrat*facB/(facB+1)
  
  elecP=tube%zPimp*CimpP-CbufP
  elecA=tube%zAimp*CimpA-CbufA
  elecB=tube%zBimp*CimpB-CbufB
  elecE=0.0d0
  do I=1,NS2
     elecP=elecP+zval(I)*C(I,2)
     elecA=elecA+zval(I)*C(I,3)
     elecB=elecB+zval(I)*C(I,4)
     elecE=elecE+zval(I)*C(I,5)
  end do
  fvec(5+4*NS2)=elecP
  fvec(6+4*NS2)=elecA
  fvec(7+4*NS2)=elecB
  fvec(8+4*NS2)=elecE

!	 Electroneutrality/Zero-current in lumen
  
  currM=0.0d0
  do I=1,NS2
     do K=2,5
        currM=currM+zval(I)*Jsol(I,1,K)
     end do
  end do
  fvec(9+4*NS2)=currM
  

  return
end subroutine fcn1C


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
