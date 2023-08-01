! NCX exchanger module

!---------------------------------------------------------------------72
!		NCX EXCHANGER AT BASOLATERAL MEMBRANE OF CELL
!		n is for Na, c is for Ca
!		p (or prime) is for the luminal compartment (M)
!	    pp (or double prime) is for the cytosolic compartment (I)
!---------------------------------------------------------------------72

subroutine compute_ncx_fluxes (var_ncx,dJNCXca5,dJNCXca6)

include 'values.h'

! input variables
double precision var_ncx(16)

! output
double precision dJNCXca5,dJNCXca6

! local variables
double precision nai,nao5,nao6,cai,cao5,cao6,epi,epo5,epo6,area5,area6,xNCX
double precision phir5, phir6, phif5, phif6, g5, g6, fmod

		nai = var_ncx(1)  !Na+
		nao5 = var_ncx(2)
		nao6 = var_ncx(3)
		cai = var_ncx(4) !Ca2+
		cao5 = var_ncx(5)
		cao6 = var_ncx(6)
		epi = var_ncx(7)
		epo5 = var_ncx(8)
		epo6 = var_ncx(9)
        area5 = var_ncx(10)
        area6 = var_ncx(11)
        xNCX = var_ncx(12)

		fmod = (cai/dKm_ncx)**2/(1.0d0 + (cai/dKm_ncx)**2)

        phir5 = dexp((gamma_ncx - 1.0d0)*F*EPref/RT*(epi-epo5))
        phir6 = dexp((gamma_ncx - 1.0d0)*F*EPref/RT*(epi-epo6))
        phif5 = dexp(gamma_ncx*F*EPref/RT*(epi-epo5))
        phif6 = dexp(gamma_ncx*F*EPref/RT*(epi-epo6))

!        g5 = (cna5**3)*cai + (cnai**3)*cao5 + (dKmNao**3)*cai + (nai**3)*dKmCao &
!            + (dKmNai**3)*cao5*(1.0d0 + cai/dKmCai) &
!            + (cnao5**3)*dKmCai*(1.d0 + (cnai/dKmNai)**3)
!        g6 = (cna6**3)*cai + (cnai**3)*cao6 + (dKmNao**3)*cai + (nai**3)*dKmCao &
!            + (dKmNai**3)*cao6*(1.0d0 + cai/dKmCai) &
!            + (cnao6**3)*dKmCai*(1.d0 + (cnai/dKmNai)**3)

        g5 = (nao5**3)*cai + (nai**3)*cao5 + (dKmNao**3)*cai + (nai**3)*dKmCao &
            + (dKmNai**3)*cao5*(1.0d0 + cai/dKmCai) &
            + (nao5**3)*dKmCai*(1.d0 + (nai/dKmNai)**3)
        g6 = (nao6**3)*cai + (nai**3)*cao6 + (dKmNao**3)*cai + (nai**3)*dKmCao &
            + (dKmNai**3)*cao6*(1.0d0 + cai/dKmCai) &
            + (nao6**3)*dKmCai*(1.d0 + (nai/dKmNai)**3)

		fac5 = (phir5*(nao5**3)*cai - phif5*(nai**3)*cao5)/g5/(1+dKsat_ncx*phir5)
        fac6 = (phir6*(nao6**3)*cai - phif6*(nai**3)*cao6)/g6/(1+dKsat_ncx*phir6)

		dJNCXca5 = area5*xNCX*fmod*fac5
		dJNCXca6 = area6*xNCX*fmod*fac6


return
end


