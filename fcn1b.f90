!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!	 This subroutine yields the non-linear equations to be solved at the inlet
!	 in order to determine the concentrations, volumes, and EP 
!	 in P, A, B, and E.
!    It is used to compute values at the entry of segments that only have
!    principal cells (except for the PT): mTAL, cTAL, DCT and IMCD.

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
	   
subroutine fcn1b(n,x,fvec,iflag,numpar,pars,tube,id)

  include 'values.h'
  include 'global.h'
  include 'defs.h'

  integer n,iflag,numpar,id
  double precision x(NDA),fvec(n), pars(numpar)
  type(membrane) :: tube
  double precision S(NUA),y(NDA)
  double precision C(NS,NC),Vol(NC),EP(NC),ph(NC)
  double precision Jvol(NC,NC), Jsol(NS,NC,NC),Jcasr
  double precision CaBT(NC)


!	 Assign concentrations, volumes, and potentials in lumen and blood
  Do J = 1, NS
     C(J,1) = pars(J) 
     C(J,6) = pars(J+NS) 
  End Do
  ph(1) = -dlog(C(12,1)/1.0d3)/dlog(10.0d0)
  ph(6) = -dlog(C(12,6)/1.0d3)/dlog(10.0d0)
  EP(6) = pars(2*NS+1)
  CPimpref = pars(2*NS+2)

!!!!!!!!!!!!!!!!   ADDED TO ACCOUNT FOR DIFFERENT BUFFER CONCENTRATIONS (AURELIE)
  if (id==1) then
		CPbuffer = CPbuftotA 
   elseif (id==2) then
		CPbuffer = CPbuftotT
   elseif (id==3) then
		CPbuffer= CPbuftotD
   elseif (id==8) then
		CPbuffer= CPbuftotIMC
   else
     print *, "wrong id", id
	 pause
  end if

!	 Assign concentrations, volumes, and potentials in P, A, B, and E
  
  Do I = 1,NS2
     C(I,2)=x(1+2*(I-1))
     C(I,5)=x(2+2*(I-1))
  End do

  ph(2) = -dlog(C(12,2)/1.0d3)/dlog(10.0d0)
  ph(5) = -dlog(C(12,5)/1.0d3)/dlog(10.0d0)
  
  Vol(2)=x(1+2*NS2)
  Vol(5)=x(2+2*NS2)
  EP(1)=x(3+2*NS2)
  EP(2)=x(4+2*NS2)
  EP(5)=x(5+2*NS2)
  

!---------------------------------------------------------------------72
!	Initialize  source terms
!---------------------------------------------------------------------72

  Do K = 1, n
     S(K) = 0.0d0
     fvec(K) = 0.0d0
  End Do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE SOURCE TERMS FOR SOLUTE I IN EACH CELLULAR COMPARTMENT
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  Do I = 1,NS2
     if (id==1 .or. id==2 .or. id==3 .or. id==8) then
        y(1+3*(I-1))=tube%conc(I,1)
     else
        print *, "wrong id", id
     end if
     y(2+3*(I-1))=C(I,2)
     y(3+3*(I-1))=C(I,5)
  End do

  if (id==1 .or. id==2 .or. id==3 .or. id==8) then
     y(1+3*NS2) = tube%vol(1)
  else
     print *, "wrong id", id
  end if
  y(2+3*NS2) = Vol(2)
  y(3+3*NS2) = Vol(5)
  y(4+3*NS2) = EP(1)
  y(5+3*NS2) = EP(2)
  y(6+3*NS2) = EP(5)
  if (id==1 .or. id==2 .or. id==3 .or. id==8) then
     y(7+3*NS2) = tube%pres
  elseif (id==-1) then
     y(7+3*NS2) = PMinit
  else
     print *,"wrong id",id
  end if


  if (id == 1) then ! mTAL
     Call qflux2A(y,Jvol,Jsol,tube%conc(:,6),tube%ep(6),tube)
  elseif (id == 2) then ! cTAL
     Call qflux2T(y,Jvol,Jsol,tube%conc(:,6),tube%ep(6),tube)
  elseif (id == 3) then ! DCT
     Call qflux2D(y,Jvol,Jsol,tube%conc(:,6),tube%ep(6),tube)
  elseif (id == 5) then ! DCT
     Call qflux2C(y,Jvol,Jsol,tube)
  elseif (id == 8) then ! IMCD
     Call qflux2IMC(y,Jvol,Jsol,tube%conc(:,6),tube%ep(6),tube)
  else
     print*, "wrong id", id
  end if

  Do I = 1, NS2
     S(1+2*(I-1)) = -Jsol(I,1,2)+Jsol(I,2,5)+Jsol(I,2,6) 
     S(2+2*(I-1)) = -Jsol(I,1,5)-Jsol(I,2,5)+Jsol(I,5,6)
  End Do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE SOURCE TERMS FOR VOLUME IN EACH CELLULAR COMPARTMENT
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  
  S(1+2*NS2) = -Jvol(1,2)+Jvol(2,5)+Jvol(2,6)
  S(2+2*NS2) = -Jvol(1,5)-Jvol(2,5)+Jvol(5,6) 

!---------------------------------------------------------------------72
!	Convert to equations of the form F(X) = 0
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!	 For I = 1-3 and 9 (non-reacting solutes)
!---------------------------------------------------------------------72
		
  Do M=1,6
     fvec(M)=S(M)
  End Do
  
  fvec(17)=S(17)
  fvec(18)=S(18)

! HCO2-/H2CO2 equations modified below
  Do M=13,15
    fvec(2*M-1)=S(2*M-1)
    fvec(2*M)=S(2*M)
  End Do

!---------------------------------------------------------------------72
!	 For CO2/HCO3/H2CO3
!---------------------------------------------------------------------72
  
  facnd=Vref/href
  
  fvec(7) = S(7)+S(9)+S(11)
  fvec(8) = S(8)+S(10)+S(12)
  fvec(9) = ph(2)-pKHCO3-dlog(abs(C(4,2)/C(5,2)))/dlog(10.0d0)
  fvec(10) = ph(5)-pKHCO3-dlog(abs(C(4,5)/C(5,5)))/dlog(10.0d0)
  fvec(11) = S(11)+Vol(2)*(tube%dkh(2)*C(6,2)-tube%dkd(2)*C(5,2))*facnd
  fvec(12) = S(12)+max(Vol(5),tube%volEinit)*(tube%dkh(5)*C(6,5)-tube%dkd(5)*C(5,5))*facnd
  
!---------------------------------------------------------------------72
!	 For HPO4(2-)/H2PO4(-)
!---------------------------------------------------------------------72

  fvec(13) = S(13)+S(15)
  fvec(14) = S(14)+S(16)
  fvec(15) = ph(2)-pKHPO4-dlog(abs(C(7,2)/C(8,2)))/dlog(10.0d0)
  fvec(16) = ph(5)-pKHPO4-dlog(abs(C(7,5)/C(8,5)))/dlog(10.0d0)

!---------------------------------------------------------------------72
!	 For NH3/NH4
!---------------------------------------------------------------------72

  fvec(19) = S(19)+S(21)
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

  fvec(1+2*NS2)=S(1+2*NS2)
  fvec(2+2*NS2)=S(2+2*NS2)

!---------------------------------------------------------------------72
!	 For EP, need to satisfy electroneutrality
!---------------------------------------------------------------------72

  volPrat=tube%volPinit/Vol(2)
  CimpP=CPimpref*volPrat
  facP=dexp(dlog(10.0d0)*(ph(2)-pKbuf))
  CbufP=CPbuffer*volPrat*facP/(facP+1)
  
  elecP=zPimp*CimpP-CbufP
  elecE=0.0d0
  do I=1,NS2
     elecP=elecP+zval(I)*C(I,2)
     elecE=elecE+zval(I)*C(I,5)
  end do
  fvec(3+2*NS2)=elecP
  fvec(4+2*NS2)=elecE
  
!	 Electroneutrality/Zero-current in lumen
  
  currM=0.0d0
  do I=1,NS2
     currM=currM+zval(I)*(Jsol(I,1,2)+Jsol(I,1,5))
  end do
  fvec(5+2*NS2)=currM
  
  return
end subroutine fcn1b


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
