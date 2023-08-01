!---------------------------------------------------------------------72
! This subroutine computes K, Cl, and NH4 fluxes across basolateral
! KCC4 cotransporters
!---------------------------------------------------------------------72

subroutine compute_kcc_fluxes (C,area5,area6,xKCC4,dJk5,dJc5,dJm5,dJk6,dJc6,dJm6)

include 'values.h'

! passed parameters and values
  double precision C(NS,NC) ! concentrations
  double precision area5, area6, xKCC4

! output (K,Cl,NH4 fluxes)
  double precision dJk5,dJc5,dJm5,dJk6,dJc6,dJm6

		bet2 = C(2,2)/bkkcc
		bet5 = C(2,5)/bkkcc
		bet6 = C(2,6)/bkkcc
	    gam2 = C(3,2)/bckcc
		gam5 = C(3,5)/bckcc
		gam6 = C(3,6)/bckcc
		dnu2 = C(11,2)/bmkcc
		dnu5 = C(11,5)/bmkcc
		dnu6 = C(11,6)/bmkcc

		sig2 = 1.d0+gam2*(1.d0+bet2+dnu2)
		sig5= 1.d0+bet5*(1.d0+gam5)+dnu5*(1.d0+gam5)
		sig6= 1.d0+bet6*(1.d0+gam6)+dnu6*(1.d0+gam6)

		rho2 = poppkcc + pkccpp*bet2*gam2 + pmccpp*dnu2*gam2 !rhopp
		rho5 = popkcc + pkccp*bet5*gam5 + pmccp*dnu5*gam5 !rhop for 5
		rho6 = popkcc + pkccp*bet6*gam6 + pmccp*dnu6*gam6 !rhop for 6

		bigsum5 = sig5*rho2 + sig2*rho5
		bigsum6 = sig6*rho2 + sig2*rho6

!		Flux across lateral membrane

		t1for5 = poppkcc*pkccp*bet5*gam5-popkcc*pkccpp*bet2*gam2
		t2for5 = pmccpp*dnu2*gam2*pkccp*bet5*gam5 &
      			-pmccp*dnu5*gam5*pkccpp*bet2*gam2
		t3for5 = poppkcc*pmccp*dnu5*gam5-popkcc*pmccpp*dnu2*gam2
	    t4for5 = pmccp*dnu5*gam5*pkccpp*bet2*gam2 &
      			-pmccpp*dnu2*gam2*pkccp*bet5*gam5
		
		dJk5 = -xKCC4*area5*(t1for5+t2for5)/bigsum5
		dJm5 = -xKCC4*area5*(t3for5+t4for5)/bigsum5
		dJc5 = dJk5 + dJm5

!		Flux across basal membrane

		t1for6 = poppkcc*pkccp*bet6*gam6-popkcc*pkccpp*bet2*gam2
		t2for6 = pmccpp*dnu2*gam2*pkccp*bet6*gam6 &
      			-pmccp*dnu6*gam6*pkccpp*bet2*gam2
		t3for6 = poppkcc*pmccp*dnu6*gam6-popkcc*pmccpp*dnu2*gam2
	    t4for6 = pmccp*dnu6*gam6*pkccpp*bet2*gam2 &
      			-pmccpp*dnu2*gam2*pkccp*bet6*gam6

		dJk6 = -xKCC4*area6*(t1for6+t2for6)/bigsum6
		dJm6 = -xKCC4*area6*(t3for6+t4for6)/bigsum6
		dJc6 = dJk6 + dJm6

return
end
