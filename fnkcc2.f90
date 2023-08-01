!---------------------------------------------------------------------72
 
!  *************************** NKCC2 kinetics **************************

!---------------------------------------------------------------------72

!    This is a subroutine to compute the NKCC2 fluxes
!	 The units for the kinetic constants are M-1 or M-1.s-1
!	 The kinetic equations and parameters come from
!		1. Mercano et al. AJP 296: F369, 2009
!		2. Nieves-Gonzales et al., AJP , 2012	

	   subroutine fnkcc2(n,aconc,Akcc2,nisoform)

	   include 'values.h'

       integer n
	   double precision nain,naout,kin,kout,clin,clout
       double precision aconc(8),Akcc2(n,n)

!---------------------------------------------------------------------72
!	 Initialize matrix
!---------------------------------------------------------------------72

	 Do j=1,n
			Do i=1,n
				Akcc2(i,j)=0.0d0
			end do
	 end do

!---------------------------------------------------------------------72
!	 Fix reaction rates
!---------------------------------------------------------------------72

	 dkon = 1.0d8

	 if (nisoform .eq. 1) then

!	 A isoform
		dK1 = 430.d0
		dK2 = 137.d0
		dK3 = 0.709d0
		dK4 = dK2

		dkff = 100000.d0
		dkbf = 1000.d0
		dkfe = 1200.70d0
!		dkbe = 120070.d0
		dkbe = dkff*dkfe/dkbf

		dK3n = 2.69d0
!	    dkffn = 1406.9d0
	    dkbfn = 140.69d0
	    dkffn = dkbe*dkbfn/dkfe

			dkoff1 =    2.322988d5
	        dkoff2 =    7.299803d5
	        dkoff3 = 1410.297996d5
	        dkoff4 =    7.299803d5
	        dkoff3n = 371.747212d5

!			dkoff1 =    3.02701686d5
!	        dkoff2 =    7.57281937d5
!	        dkoff3 = 1410.29908357d5
!	        dkoff4 =    7.71627760d5
!	        dkoff3n = 371.74721190d5

		else if (nisoform .eq. 2) then

!	 B isoform
			dK1 = 728.d0
			dK2 = 1000d0
			dK3 = 0.135d0
			dK4 = dK2

			dkff = 100000.d0
			dkbf = 1049.6d0
			dkfe = 7626.6d0
!			dkbe = 727000.d0
			dkbe = dkff*dkfe/dkbf

			dK3n = 1.75d0
!			dkffn = 3730.d0
			dkbfn = 39.1d0
	        dkffn = dkbe*dkbfn/dkfe

			dkoff1 =    1.373683d5
	        dkoff2 =    1.0000000d5
	        dkoff3 = 7407.4074074d5
	        dkoff4 =    1.0000000d5
	        dkoff3n = 571.8206770d5


		else

!	 F isoform
			dK1 = 138.d0
			dK2 = 751.d0
			dK3 = 0.10d0
			dK4 = dK2

			dkff = 100000.d0
			dkbf = 1000.d0
			dkfe = 4872.8d0
!			dkbe = 487000.d0
			dkbe = dkff*dkfe/dkbf

			dK3n = 1.33d0
!			dkffn = 4180.d0
			dkbfn = 41.8d0
	        dkffn = dkbe*dkbfn/dkfe
  
			dkoff1 = 0.72727272d6
			dkoff2 = 0.13315757d6
			dkoff3 = 1000.00000d6
			dkoff4 = 0.13315757d6
			dkoff3n = 75.39772299d6

		end if


!	 determine off constants
!	 dkoff1 = dkon/dK1
!	 dkoff2 = dkon/dK2
!	 dkoff3 = dkon/dK3
!	 dkoff4 = dkon/dK4
!	 dkoff3n = dkon/dK3n

!---------------------------------------------------------------------72
!	 Change dimensions of concentrations (convert from mM to M)	
!---------------------------------------------------------------------72

	 naout = aconc(1)*1.d-3
	 nain  = aconc(2)*1.d-3		
	 kout  = aconc(3)*1.d-3	 
	 kin   = aconc(4)*1.d-3
	 clout = aconc(5)*1.d-3
	 clin  = aconc(6)*1.d-3
	 ammout = aconc(7)*1.d-3
	 ammin  = aconc(8)*1.d-3

!---------------------------------------------------------------------72
!	  Determine matrix to inverse
!	  Indices 1-10 represent E0,E1,..E9. 
!	  Indices 11-15 represent E3N, E4N, E5N, E6N and E7N
!---------------------------------------------------------------------72
!
!	  The first equation expressses conservation of species
!	  E0 + E1 + ... + E15 = 1. In other words, Akcc2(1,i = 1-15) = 1

	 Do j=1,n
			Akcc2(1,j) = 1.0d0
	 end do

		Akcc2(2,1) = + dkon*naout 
		Akcc2(2,2) = - dkon*clout - dkoff1
		Akcc2(2,3) = +dkoff2

		Akcc2(3,2) = + dkon*clout
	    Akcc2(3,3) = - dkon*kout -dkon*ammout - dkoff2
	    Akcc2(3,4) = + dkoff3
		Akcc2(3,11) = + dkoff3n

		Akcc2(4,3) = + dkon*kout 
		Akcc2(4,4) = - dkon*clout - dkoff3
		Akcc2(4,5) = + dkoff4
		
		Akcc2(5,4) = + dkon*clout
	    Akcc2(5,5) = - dkff - dkoff4
	    Akcc2(5,6) = + dkbf

		Akcc2(6,5) = + dkff
	    Akcc2(6,6) = - dkoff1 - dkbf
	    Akcc2(6,7) = + dkon*nain

		Akcc2(7,6) = + dkoff1
	    Akcc2(7,7) = - dkon*nain - dkoff2
	    Akcc2(7,8) = + dkon*clin

		Akcc2(8,7) = + dkoff2
		Akcc2(8,8) = - dkon*clin - dkoff3
	    Akcc2(8,9) = + dkon*kin 

		Akcc2(9,8) = + dkoff3
		Akcc2(9,9) = - dkon*kin - dkon*ammin - dkoff4
		Akcc2(9,10) = + dkon*clin	
		Akcc2(9,15) = + dkoff3n

		Akcc2(10,9) = + dkoff4
	    Akcc2(10,10) = - dkon*clin - dkfe
	    Akcc2(10,1) = + dkbe

!		Additional equations with NH4+ 
		 
	 	Akcc2(11,3) = + dkon*ammout
		Akcc2(11,11) = - dkon*clout - dkoff3n
		Akcc2(11,12) = + dkoff4

		Akcc2(12,11) = + dkon*clout
	    Akcc2(12,12) = - dkffn - dkoff4
	    Akcc2(12,13) = + dkbfn

		Akcc2(13,12) = + dkffn
	    Akcc2(13,13) = - dkoff1 - dkbfn
	    Akcc2(13,14) = + dkon*nain

		Akcc2(14,13) = + dkoff1
	    Akcc2(14,14) = - dkon*nain - dkoff2
	    Akcc2(14,15) = + dkon*clin

		Akcc2(15,14) = + dkoff2
		Akcc2(15,15) = - dkon*clin - dkoff3n
	    Akcc2(15,9) = + dkon*ammin 

		 return
       end


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
subroutine compute_nkcc2_flux ( C,area,xNKCC2,bn2,bk2,bc2,bm2, popnkcc,pnkccp,&
        pnmccp,poppnkcc,pnkccpp,pnmccpp, dJnNKCC2,dJkNKCC2,dJcNKCC2,dJmNKCC2 )

	   include 'values.h'
	   include 'global.h'

! passed variables
  double precision C(NS,NC),area
  double precision bn2,bk2,bc2,bm2
  double precision popnkcc,pnkccp,pnmccp,poppnkcc,pnkccpp,pnmccpp

! output variables
  ! fluxes (Na, K, Cl, NH4)
  double precision dJnNKCC2,dJkNKCC2,dJcNKCC2,dJmNKCC2 

         integer n
	    alp1 = C(1,1)/bn2
		alp2 = C(1,2)/bn2
		bet1 = C(2,1)/bk2
		bet2 = C(2,2)/bk2
	    gam1 = C(3,1)/bc2
		gam2 = C(3,2)/bc2
		dnu1 = C(11,1)/bm2
		dnu2 = C(11,2)/bm2

		sig1 = 1.d0 +alp1 + alp1*gam1*(1.d0 + bet1 + bet1*gam1&
      			+ dnu1 +dnu1*gam1)
		sig2 = 1.d0 + gam2*(1.d0 + bet2 + bet2*gam2 + bet2*gam2*alp2&
      			+ dnu2 + dnu2*gam2 + dnu2*gam2*alp2)

		rho1 = popnkcc + pnkccp*alp1*bet1*(gam1**2) + &
      			pnmccp*alp1*dnu1*(gam1**2)
		rho2 = poppnkcc + pnkccpp*alp2*bet2*(gam2**2) + &
      			pnmccpp*alp2*dnu2*(gam2**2)

		bigsum = sig1*rho2 + sig2*rho1

!		Flux across lateral membrane

		t1 = poppnkcc*pnkccp*alp1*bet1*(gam1**2) - &
      			popnkcc*pnkccpp*alp2*bet2*(gam2**2)
		t2 = poppnkcc*pnmccp*alp1*dnu1*(gam1**2) - &
      			popnkcc*pnmccpp*alp2*dnu2*(gam2**2)
		t3 = pnmccpp*alp2*dnu2*(gam2**2)*pnkccp*alp1*bet1*(gam1**2)
	    t4 = pnmccp*alp1*dnu1*(gam1**2)*pnkccpp*alp2*bet2*(gam2**2)
		
		dJnNKCC2 = xNKCC2*area*(t1+t2)/bigsum
		dJkNKCC2 = xNKCC2*area*(t1+t3-t4)/bigsum
		dJmNKCC2 = xNKCC2*area*(t2+t4-t3)/bigsum
		dJcNKCC2 = 2*dJnNKCC2


return
end
