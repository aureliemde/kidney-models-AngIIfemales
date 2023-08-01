!---------------------------------------------------------------------72
!  *************************** SGLT fluxes ***************************
!---------------------------------------------------------------------72
!    This subroutine computes the fluxes across SGLT1 (Na-glucose)
!	 cotransporters, based on the model of Parent et al. (J Memb Biol 1992)
!	 with the modifications described by Eskandari et al. (J Memb Biol 2005)
	
subroutine sglt(n,nagluparam,area,fluxsglt)
  include 'values.h'
  integer n
  double precision nao,nai,gluo,glui,epo,epi,area
  double precision nagluparam(7),fluxsglt(2)

!    Assign concentrations (in M, not mM) and potentials
  nao = nagluparam(1)*1.0d-3 
  nai = nagluparam(2)*1.0d-3
  gluo = nagluparam(3)*1.0d-3
  glui = nagluparam(4)*1.0d-3
  epo = nagluparam(5)
  epi = nagluparam(6)
  CT = nagluparam(7)
  pot = F*(epi-epo)*EPref/RT
!	 The reaction rates for SGLT1 come from Wright et al. Physiol Rev 2011

  if (n .eq. 1) then
  delta = 0.70
  alphap = 0.30
  alphapp = 0.00
  dk12 = 140000.0*dexp(-pot*alphap)*(nao**2)
  dk21 = 300.0*dexp(+pot*alphap)
  dk23 = 45000.0*gluo
  dk32 = 20.0
  dk34 = 50.0
  dk43 = 50.0
  dk45 = 800.0
  dk54 = 190000.0*glui 
  dk56 = 5.0*dexp(-pot*alphapp)
  dk65 = 2250.0*dexp(+pot*alphapp)*(nai**2) 
  dk61 = 25.0*dexp(-pot*delta) 
  dk16 = 600.0*dexp(+pot*delta) 
  dk25 = 0.01
  dk52 = 0.0005
  nstoich = 2
  end if
!	 The reaction rates for SGLT2 come from Mackenzie (J Biol Chem 1996)
  if (n .eq. 2) then
  delta = 0.70
  alphap = 0.30
  alphapp = 0.00
  dk12 = 20000.0*dexp(-pot*alphap/2.0)*nao
  dk21 = 400.0*dexp(+pot*alphap/2.0)
  dk23 = 10000.0*gluo
  dk32 = 20.0
  dk34 = 50.0
  dk43 = 50.0
  dk45 = 800.0
  dk54 = 6700000.0*glui
  dk56 = 48.0*dexp(-pot*alphapp/2.0)
  dk65 = 50.0*dexp(+pot*alphapp/2.0)*nai
  dk61 = 35.0*dexp(-pot*delta/2.0)
  dk16 = 100.0*dexp(+pot*delta/2.0)
  dk25 = 0
  dk52 = 0
  nstoich = 1
  end if

!	 Compute the KAT terms and the concentrations of each enzyme state

  delC1 = dk21*dk32*dk43*dk52*dk61 + dk21*dk32*dk45*dk52*dk61 + dk21*dk34*dk45*dk52*dk61 + & 
		  dk21*dk32*dk43*dk54*dk61 + dk21*dk32*dk43*dk56*dk61 + dk25*dk32*dk43*dk56*dk61 + &
		  dk21*dk32*dk45*dk56*dk61 + dk25*dk32*dk45*dk56*dk61 + dk21*dk34*dk45*dk56*dk61 + &
		  dk23*dk34*dk45*dk56*dk61 + dk25*dk34*dk45*dk56*dk61 + dk21*dk32*dk43*dk52*dk65 + &
		  dk21*dk32*dk45*dk52*dk65 + dk21*dk34*dk45*dk52*dk65 + dk21*dk32*dk43*dk54*dk65
  delC2 = dk12*dk32*dk43*dk52*dk61 + dk12*dk32*dk45*dk52*dk61 + dk12*dk34*dk45*dk52*dk61 + & 
		  dk12*dk32*dk43*dk54*dk61 + dk12*dk32*dk43*dk56*dk61 + dk12*dk32*dk45*dk56*dk61 + &
		  dk12*dk34*dk45*dk56*dk61 + dk12*dk32*dk43*dk52*dk65 + dk16*dk32*dk43*dk52*dk65 + &
		  dk12*dk32*dk45*dk52*dk65 + dk16*dk32*dk45*dk52*dk65 + dk12*dk34*dk45*dk52*dk65 + &
		  dk16*dk34*dk45*dk52*dk65 + dk12*dk32*dk43*dk54*dk65 + dk16*dk32*dk43*dk54*dk65
  delC3 = dk12*dk23*dk43*dk52*dk61 + dk12*dk23*dk45*dk52*dk61 + dk12*dk23*dk43*dk54*dk61 + & 
		  dk12*dk25*dk43*dk54*dk61 + dk12*dk23*dk43*dk56*dk61 + dk12*dk23*dk45*dk56*dk61 + &
		  dk12*dk23*dk43*dk52*dk65 + dk16*dk23*dk43*dk52*dk65 + dk12*dk23*dk45*dk52*dk65 + &
		  dk16*dk23*dk45*dk52*dk65 + dk16*dk21*dk43*dk54*dk65 + dk12*dk23*dk43*dk54*dk65 + & 
		  dk16*dk23*dk43*dk54*dk65 + dk12*dk25*dk43*dk54*dk65 + dk16*dk25*dk43*dk54*dk65

  delC4 = dk12*dk23*dk34*dk52*dk61 + dk12*dk25*dk32*dk54*dk61 + dk12*dk23*dk34*dk54*dk61 + & 
		  dk12*dk25*dk34*dk54*dk61 + dk12*dk23*dk34*dk56*dk61 + dk12*dk23*dk34*dk52*dk65 + &
		  dk16*dk23*dk34*dk52*dk65 + dk16*dk21*dk32*dk54*dk65 + dk12*dk25*dk32*dk54*dk65 + &
		  dk16*dk25*dk32*dk54*dk65 + dk16*dk21*dk34*dk54*dk65 + dk12*dk23*dk34*dk54*dk65 + &
		  dk16*dk23*dk34*dk54*dk65 + dk12*dk25*dk34*dk54*dk65 + dk16*dk25*dk34*dk54*dk65
  delC5 = dk12*dk25*dk32*dk43*dk61 + dk12*dk25*dk32*dk45*dk61 + dk12*dk23*dk34*dk45*dk61 + & 
		  dk12*dk25*dk34*dk45*dk61 + dk16*dk21*dk32*dk43*dk65 + dk12*dk25*dk32*dk43*dk65 + &
		  dk16*dk25*dk32*dk43*dk65 + dk16*dk21*dk32*dk45*dk65 + dk12*dk25*dk32*dk45*dk65 + &
		  dk16*dk25*dk32*dk45*dk65 + dk16*dk21*dk34*dk45*dk65 + dk12*dk23*dk34*dk45*dk65 + &
		  dk16*dk23*dk34*dk45*dk65 + dk12*dk25*dk34*dk45*dk65 + dk16*dk25*dk34*dk45*dk65
  delC6 = dk16*dk21*dk32*dk43*dk52 + dk16*dk21*dk32*dk45*dk52 + dk16*dk21*dk34*dk45*dk52 + & 
		  dk16*dk21*dk32*dk43*dk54 + dk16*dk21*dk32*dk43*dk56 + dk12*dk25*dk32*dk43*dk56 + &
		  dk16*dk25*dk32*dk43*dk56 + dk16*dk21*dk32*dk45*dk56 + dk12*dk25*dk32*dk45*dk56 + &
		  dk16*dk25*dk32*dk45*dk56 + dk16*dk21*dk34*dk45*dk56 + dk12*dk23*dk34*dk45*dk56 + &
		  dk16*dk23*dk34*dk45*dk56 + dk12*dk25*dk34*dk45*dk56 + dk16*dk25*dk34*dk45*dk56
  SumdelC = delC1 + delC2 + delC3 + delC4 + delC5 + delC6
  C1 = delC1/SumdelC
  C2 = delC2/SumdelC
  C3 = delC3/SumdelC
  C4 = delC4/SumdelC
  C5 = delC5/SumdelC
  C6 = delC6/SumdelC
!	 Compute the sodium and glucose fluxes
  fluxsglt(1) = area*CT*nstoich*(dk34*C3 - dk43*C4 + dk25*C2 - dk52*C5)
  fluxsglt(2) = area*CT*1.0*(dk34*C3 - dk43*C4)
  return
end subroutine sglt


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
