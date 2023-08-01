  !   This subroutine sets interstitial concentrations for a given nephron segment

subroutine set_intconc ( tube, iN, ireg, pos )

  include 'values.h'
  include 'global.h'
  include 'defs.h'

  ! passed variables
  integer :: iN  ! number of spatial intervals
  type (membrane) :: tube(0:iN)
  integer :: ireg  ! specify region, cortex (1), OM (2), IM (3)
  double precision :: pos(0:iN)  ! nondimensional distance from region upper boundary

  ! local variables
  integer :: j
  logical :: bdiabetes  ! true if simulating a diabetic kidney
  logical :: furo ! true in the presence of furosemide
  double precision :: facbic, facpho, facamm, fachco2

  bdiabetes = .false.

  ! Total ammonia at the uppermost cortex
  TotAmmCT = 0.203d0

  ! concentrations in cortex and at cortical medullary boundary
  TotSodCM = 144.00d0
  TotPotCM = 4.0d0
  TotCloCM = 121.8d0  !!! MODIFIED BELOW FOR ELECTRONEUTRALITY
  TotBicCM = 25.0d0 
  TotHcoCM = 4.41d-3
  TotCo2CM = 1.50d0
  TotPhoCM = 2.60d0
  TotureaCM = 5.0d0
  TotAmmCM = 1.0d0 
  TotHco2CM = 1.0d0
  TotCaCM = 1.25d0
  if (.not. bdiabetes) then
    TotgluCM = 5.0d0   ! normal kidney
  else
    TotgluCM = 25.0d0  ! diabetic kidney
  end if

  ! concentrations at OM-IM junction 
  TotSodOI = 284.0d0
  TotPotOI = 10.0d0
  TotCloOI = 264.94d0  !!! MODIFIED BELOW FOR ELECTRONEUTRALITY
  TotBicOI = 25.0d0 
  TotHcoOI = 4.41d-3
  TotCo2OI = 1.50d0
  TotPhoOI = 3.90d0
  TotureaOI = 20.0d0 
  TotAmmOI = 3.90d0
  TotHco2OI = TotHco2CM*1.d0 !!!!! Interstitial gradient of HCO2-/H2CO2 species
  TotCaOI = 2.50d0
  xIS = 0.60/2.00d0
  TotGluOI = (5.0 + 1.0/xIS)*(TotGluCM/5.0d0) ! TotGluOI is set so that TotGluIS equals 6.0

  ! concentrations at papillary tip
  TotSodPap = 284.0d0 + 15.0d0 !!!! ADDED 15 mM AS AT OI JUNCTION
  TotPotPap = 20.d0
  TotCloPap = 279.92d0 + 15.0d0 !!!! MODIFIED BELOW FOR ELECTRONEUTRALITY
  TotBicPap = 25.0d0
  TotHcoPap = 4.41d-3
  TotCo2Pap = 1.50d0
  TotPhoPap = 3.90d0
  TotAmmPap = 8.95d0
  TotureaPap = 500.d0
  TotHco2Pap = TotHco2CM*1.0d0
  TotCaPap = 4.0d0
  TotGluPap = 8.50d0*(TotGluCM/5.0d0) !!!!! Based on Hervy and Thomas, AJP Renal 2003

  if (FAngAqp2) then
    TotureaPap = 200.d0
  end if


  if (ireg == 1) then  ! cortex

     tube(0)%ph(6) = 7.323d0

     facbic=dexp(dlog(10.0d0)*(tube(0)%ph(6)-pKHCO3))
     facpot=dexp(dlog(10.0d0)*(tube(0)%ph(6)-pKHPO4))
     facamm=dexp(dlog(10.0d0)*(tube(0)%ph(6)-pKNH3))
     fachco2=dexp(dlog(10.0d0)*(tube(0)%ph(6)-pKHCO2))

     do j = 0, iN
        tube(j)%pH(6) = tube(0)%pH(6)

        tube(j)%conc(1,6) = TotSodCM
        tube(j)%conc(2,6) = TotPotCM
        tube(j)%conc(3,6) = TotCloCM 
        tube(j)%conc(4,6) = TotBicCM
        tube(j)%conc(5,6) = TotHcoCM
        tube(j)%conc(6,6) = TotCo2CM

        tube(j)%conc(7,6) = TotPhoCM*facpot/(1.d0+facpot)
        tube(j)%conc(8,6) = TotPhoCM/(1.d0+facpot)

        tube(j)%conc(9,6) = TotureaCM

        Ammtotz=(TotAmmCT+(TotAmmCM-TotAmmCT)*pos(j))
        tube(j)%conc(3,6)=TotCloCM + (Ammtotz-TotAmmCM)
        tube(j)%conc(10,6)=Ammtotz*facamm/(1.0d0+facamm)
        tube(j)%conc(11,6)=Ammtotz/(1.0d0+facamm)

        tube(j)%conc(12,6) = dexp(-dlog(10.0d0)*tube(j)%ph(6))*1d3 !Convert to mM
        
        tube(j)%conc(13,6) = TotHco2CM*fachco2/(1.d0+fachco2)
        tube(j)%conc(14,6) = TotHco2CM/(1.d0+fachco2)

        tube(j)%conc(15,6) = TotGluCM

        tube(j)%conc(16,6) = TotCaCM

!       Set chloride concentration so as to satisfy electroneutrality
        ElecS = 0.0d0
        Do I = 1, NS
            ElecS = ElecS + zval(I)*tube(j)%conc(I,6)
        end do
        tube(j)%conc(3,6) = tube(j)%conc(3,6) + ElecS

     end do

  elseif (ireg == 2) then  ! OM

     tube(0)%ph(6) = 7.323d0
     
     facbic=dexp(dlog(10.0d0)*(tube(0)%ph(6)-pKHCO3))
     facpho=dexp(dlog(10.0d0)*(tube(0)%ph(6)-pKHPO4))
     facamm=dexp(dlog(10.0d0)*(tube(0)%ph(6)-pKNH3))
     fachco2=dexp(dlog(10.0d0)*(tube(0)%ph(6)-pKHCO2))

     do j = 0, iN
        tube(j)%ph(6) = tube(0)%ph(6)

        tube(j)%conc(1,6) = TotSodCM + (TotSodOI-TotSodCM)*pos(j)
        tube(j)%conc(2,6) = TotPotCM + (TotPotOI-TotPotCM)*pos(j)
        tube(j)%conc(3,6) = TotCloCM + (TotCloOI-TotCloCM)*pos(j)
        tube(j)%conc(4,6) = TotBicCM + (TotBicOI-TotBicCM)*pos(j)
        tube(j)%conc(5,6) = TotHcoCM + (TotHcoOI-TotHcoCM)*pos(j)
        tube(j)%conc(6,6) = TotCo2CM + (TotCo2OI-TotCo2CM)*pos(j)
        
        Phototz=(TotPhoCM+(TotPhoOI-TotPhoCM)*pos(j))
        tube(j)%conc(7,6) = Phototz*facpho/(1.0d0+facpho)
        tube(j)%conc(8,6) = Phototz/(1.0d0+facpho)        

        tube(j)%conc(9,6) = TotureaCM + (TotureaOI-TotureaCM)*pos(j)
        
        Ammtotz=TotAmmCM+(TotAmmOI-TotAmmCM)*pos(j)
        tube(j)%conc(10,6)=Ammtotz*facamm/(1.0d0+facamm)
        tube(j)%conc(11,6)=Ammtotz/(1.0d0+facamm)

        tube(j)%conc(12,6) = dexp(-dlog(10.0d0)*tube(j)%ph(6))*1d3 !Convert to mM

        Hco2totz = TotHco2CM + (TotHco2OI-TotHco2CM)*pos(j)
        tube(j)%conc(13,6) = Hco2totz*fachco2/(1.d0+fachco2)
        tube(j)%conc(14,6) = Hco2totz/(1.d0+fachco2)

        tube(j)%conc(15,6) = TotGluCM + (TotGluOI-TotGluCM)*pos(j)

        tube(j)%conc(16,6) = TotCaCM + (TotCaOI-TotCaCM)*pos(j)

!       Set chloride concentration so as to satisfy electroneutrality
        ElecS = 0.0d0
        Do I = 1, NS
            ElecS = ElecS + zval(I)*tube(j)%conc(I,6)
        end do
        tube(j)%conc(3,6) = tube(j)%conc(3,6) + ElecS
     end do

  elseif (ireg == 3) then  ! IM

     tube(0)%ph(6) = 7.323d0

     facbic=dexp(dlog(10.0d0)*(tube(0)%ph(6)-pKHCO3))
     facpho=dexp(dlog(10.0d0)*(tube(0)%ph(6)-pKHPO4))
     facamm=dexp(dlog(10.0d0)*(tube(0)%ph(6)-pKNH3))
     fachco2=dexp(dlog(10.0d0)*(tube(0)%ph(6)-pKHCO2))

     do j = 0, iN
        tube(j)%ph(6) = tube(0)%ph(6)
        
        tube(j)%conc(1,6) = TotSodOI + (TotSodPap-TotSodOI)*pos(j)
        tube(j)%conc(2,6) = TotPotOI + (TotPotPap-TotPotOI)*pos(j)
        tube(j)%conc(3,6) = TotCloOI + (TotCloPap-TotCloOI)*pos(j)
        tube(j)%conc(4,6) = TotBicOI + (TotBicPap-TotBicOI)*pos(j)
        tube(j)%conc(5,6) = TotHcoOI + (TotHcoPap-TotHcoOI)*pos(j)
        tube(j)%conc(6,6) = TotCo2OI + (TotCo2Pap-TotCo2OI)*pos(j)
        
        Phototz=(TotPhoOI+(TotPhoPap-TotPhoOI)*pos(j))
        tube(j)%conc(7,6) = Phototz*facpho/(1.0d0+facpho)
        tube(j)%conc(8,6) = Phototz/(1.0d0+facpho)

        tube(j)%conc(9,6) = TotureaOI+(TotureaPap-TotureaOI)*pos(j)
        
        Ammtotz=TotAmmOI+(TotAmmPap-TotAmmOI)*pos(j)
        tube(j)%conc(10,6)=Ammtotz*facamm/(1.0d0+facamm)
        tube(j)%conc(11,6)=Ammtotz/(1.0d0+facamm)

        tube(j)%conc(12,6) = dexp(-dlog(10.0d0)*tube(0)%ph(6))*1d3

        Hco2totz = TotHco2OI + (TotHco2pap-TotHco2OI)*pos(j)
        tube(j)%conc(13,6) = Hco2totz*fachco2/(1.d0+fachco2)
        tube(j)%conc(14,6) = Hco2totz/(1.d0+fachco2)

        tube(j)%conc(15,6) = TotGluOI + (TotGluPap-TotGluOI)*pos(j)

        tube(j)%conc(16,6) = TotCaOI + (TotCaPap-TotCaOI)*pos(j)

!       Set chloride concentration so as to satisfy electroneutrality
        ElecS = 0.0d0
        Do I = 1, NS
            ElecS = ElecS + zval(I)*tube(j)%conc(I,6)
        end do
        tube(j)%conc(3,6) = tube(j)%conc(3,6) + ElecS

     end do

  else
     print *, "wrong region ID", id
  end if

end subroutine set_intconc
