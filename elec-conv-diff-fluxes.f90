! Compute electro-convective-diffusive fluxes

subroutine compute_ecd_fluxes (C,EP,area,sig,h,Jvol,Jsol,delmu)

  include 'values.h'
  include 'global.h'

! passed variables
  double precision C(NS,NC), EP(NC)
  double precision area(NC,NC), h(NS,NC,NC), sig(NS,NC,NC)
  
! outputs
  double precision Jsol(NS,NC,NC), Jvol(NC,NC)
  double precision delmu(NS,NC,NC)

! local variables
  double precision dmu(NS,NC)
  double precision dint, XI
  double precision concdiff, concmean, convect, dimless
  integer I, K, L

  double precision eps

  eps = 1.0d-6 !For exponential (Peclet) terms

!	 define electrochemical potential of each solute in each compartment
	
  Do I = 1, NS2
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
              Jsol(I,K,L)=area(K,L)*h(I,K,L)*(C(I,K)-C(I,L))
           else 
              Jsol(I,K,L)=area(K,L)*h(I,K,L)*XI *(C(I,K)-C(I,L)*dint) /(1.0d0-dint)
           end if

!		convective component

           concdiff=C(I,K)-C(I,L)
           if (dabs(concdiff) >  eps)	then
              concmean=(C(I,K)-C(I,L))/dlog(abs(C(I,K)/C(I,L)))
              dimless=(Pfref*Vwbar*Cref)/href ! Account for normalizing factors
              convect=(1.0d0-sig(I,K,L))*concmean*Jvol(K,L)*dimless
              Jsol(I,K,L) = Jsol(I,K,L)+convect
           end if

!		define driving force for coupled fluxes below

           delmu(I,K,L)=dmu(I,K)-dmu(I,L) !for coupled flux calculations
           
        End Do
     End Do
  End Do


  return
end subroutine compute_ecd_fluxes
