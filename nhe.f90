! NHE3 exchanger modules

!---------------------------------------------------------------------72
!		NHE3 EXCHANGER AT LUMINAL MEMBRANE OF CELL
!		n is for Na, c is for Cl
!		p (or prime) is for the luminal compartment (M)
!	    pp (or double prime) is for the cytosolic compartment (I)
!---------------------------------------------------------------------72

subroutine compute_nhe3_fluxes (C,area,xNHE3,dJNHEsod,dJNHEprot,dJNHEamm)

include 'values.h'

! input variables
double precision C(NS,NC)
double precision area, xNHE3

! output
double precision dJNHEsod,dJNHEprot,dJNHEamm

		ap = C(1,1)	 !Na+
		bp = C(12,1) !H+
		cp = C(11,1) !NH4+
	    app = C(1,2)
		bpp = C(12,2)
		cpp = C(11,2)
		alp=ap/dKaNH
		alpp=app/dKaNH
		betap=bp/dKbNH
		betapp=bpp/dKbNH
		gamp=cp/dKcNH
		gampp=cpp/dKcNH

		fmod=fMNH*C(12,2)/(C(12,2)+dKINH)
		sum1=(1+alp+betap+gamp)*(PaNH*alpp+PbNH*betapp+PcNH*gampp)
		sum2=(1+alpp+betapp+gampp)*(PaNH*alp+PbNH*betap+PcNH*gamp)
		sum=sum1+sum2
		termNaH = fmod*PaNH*PbNH*(alp*betapp-alpp*betap)
		termNaNH4 = fmod*PaNH*PcNH*(alp*gampp-alpp*gamp)
		termHNH4 = fmod*PbNH*PcNH*(betap*gampp-betapp*gamp)
		dJNHEsod = area*xNHE3*(termNaH+termNaNH4)/sum
		dJNHEprot = area*xNHE3*(-termNaH+termHNH4)/sum
		dJNHEamm = area*xNHE3*(-termNaNH4-termHNH4)/sum

return
end


