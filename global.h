!---------------------------------------------------------------------72

!	Dimension the concentration and valence arrays:
!   zval(n) is the valence of solute n

	Dimension zval(NS)
    Dimension fENaC(0:NZ),fNCC(0:NZ)
    Dimension fNCX_dct(0:NZ),fPMCA_dct(0:NZ),fTRPV5_dct(0:NZ)
    Dimension fNCX_cnt(0:NZ),fPMCA_cnt(0:NZ),fTRPV5_cnt(0:NZ)
    Dimension fTRPV4_cnt(0:NZ)
    Dimension fCaTransPT(0:NZ), fCaParaPT(0:NZ)
	Dimension PTtorque(0:NZ)
	Double precision LumImperm, BathImperm
	logical :: FAng, FAngAqp2
	logical :: FAngProx, FAngDist

!---------------------------------------------------------------------72
!     Specify common blocks   
!---------------------------------------------------------------------72

       Common /length/ LzPT,LzA,LzT,LzD,LzC
	   Common /length/ LzCCD,LzOMC,LzIMC

       Common zval

	   Common /AngIIgen/ FAng, FAngAqp2
	   Common /AngIIexp/ FAngProx,FAngDist
	   Common /options/ ndiabetes
	   Common/ FMratios/ FM_cENaC,FM_mENaC

	   Common /scaling/ fscaleENaC,fscaleNaK,fscaleNaK_PC
	   Common /pheffects/ hENaC_CNT,hENaC_CCD,hENaC_OMC,hENaC_IMC
	   Common /pheffects/ hROMK_CNT,hROMK_CCD,hROMK_OMC,hROMK_IMC
	   Common /pheffects/ hCltj_CNT,hCltj_CCD,hCltj_OMC,hCltj_IMC

	   Common /cond1/ sngfr, dctout_ref,dctout_flow
	   Common /cond2/ BathImperm,LumImperm
	   Common /cond3/ PTinitVol,PTtorque

	   Common /NaPO4/ xNaPiIIaPT,xNaPiIIbPT,xNaPiIIcPT,xPit2PT
	   Common /NaIMCD/ fENaC,fNCC
	   Common /TRPV5/ xTRPV5_dct,xTRPV5_cnt
	   Common /Calcium/ fNCX_dct,fPMCA_dct,fTRPV5_dct
	   Common /Calcium/ fNCX_cnt,fPMCA_cnt,fTRPV5_cnt
	   Common /Calcium/ xPTRPV4_cnt,fTRPV4_cnt
	   Common /Calcium/ fCaTransPT,fCaParaPT

	   Common /gradOM/ TotSodOI,TotSodCM,TotPotOI,TotPotCM
	   Common /gradOM/ TotCloOI,TotCloCM,TotBicOI,TotBicCM
	   Common /gradOM/ TotHcoOI,TotHcoCM,TotCo2OI,TotCo2CM
	   Common /gradOM/ TotPhoOI,TotPhoCM,TotureaOI,TotureaCM
	   Common /gradOM/ TotAmmOI,TotAmmCM,TotHco3OI,TotHco3CM
	   Common /gradOM/ TotGluOI,TotGluCM,TotHco2OI,TotHco2CM
	   Common /gradCJ/ TotCloCJ,TotCloCT,TotAmmCJ,TotAmmCT
	   Common /gradIM/ TotSodPap,TotPotPap,TotCloPap
	   Common /gradIM/ TotBicPap,TotHcoPap,TotCo2Pap
	   Common /gradIM/ TotPhoPap,TotAmmPap,TotureaPap
	   Common /gradIM/ TotHco2Pap,TotGluPap,GluInt,Hco2Int

	  
