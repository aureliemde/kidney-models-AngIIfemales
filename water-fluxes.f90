!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	COMPUTE WATER FLUXES
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	 The hydraulic and oncotic pressures are made non-dimensional 
!	 by dividing by RT*Cref

subroutine compute_water_fluxes (C,PM,PB,Vol,Vol0,VolEinit,VolPinit,CPimpref, &
				VolAinit,CAimpref,VolBinit,CBimpref,area,sig,dLPV,compl,Jvol)

include 'values.h'
include 'global.h'

! input variables
double precision C(NS,NC),area(NC,NC),sig(NS,NC,NC),dLPV(NC,NC)
double precision Vol(NC),Vol0,VolEinit,VolPinit,VolAinit,VolBinit
double precision CPimpref,CAimpref,CBimpref
double precision PM,PB,compl

! outputs
double precision Jvol(NC,NC)

! local variables
double precision PRES(NC),ONC(NC)
integer K, L, J

!	 The hydraulic and oncotic pressures are made non-dimensional 
!	 by dividing by RT*Cref

	 PRES(1)=PM/(RTosm*Cref)
	 PRES(2)=PRES(1)
	 PRES(3)=PRES(1)
	 PRES(4)=PRES(1)
	 PRES(6)=PbloodPT/(RTosm*Cref)

!	 PRES(5)=PRES(2)+(Vol(5)/VolEinit-1.0d0)/(compl)/(RTosm*Cref)
	 PRES(5)=PRES(6)+(Vol(5)/VolEinit-1.0d0)/(compl)/(RTosm*Cref)

	 ONC(1)=LumImperm*PTinitVol/Vol(1)
	 ONC(2)=CPimpref*VolPinit/Vol(2)
	 ONC(3)=CAimpref*VolAinit/Vol(3)
	 ONC(4)=CBimpref*VolBinit/Vol(4)
	 ONC(5)=0.0d0
     ONC(6)=BathImperm


	 Do K = 1, NC-1
		Do L = K+1, NC
			OSM=0.0d0
			Do J = 1, NS2
			    if (J .lt. 13 .or. J .gt. 14) then
				    OSM=OSM+sig(J,K,L)*(C(J,K)-C(J,L))
				end if
			End Do
			Jvol(K,L)=area(K,L)*dLPV(K,L)*(PRES(K)-PRES(L)-ONC(K)+ONC(L)-OSM)
		End Do
	 end Do

	 Do K = 1, NC-1
		Do L = K+1, NC
			Jvol(L,K)=-Jvol(K,L)
		End Do
	 End Do


return
end


