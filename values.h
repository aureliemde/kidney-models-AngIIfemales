!---------------------------------------------------------------------72

      Implicit Double Precision (a-h,o-y)
      Implicit Integer (i-n)

!     NS = Number of solutes
      Parameter (NS = 16)
      Parameter (NS2 = 16)
      Parameter (NSPT = 16)

!     NC = Number of compartments
      Parameter (NC = 6)

!     NZ = Number of intervals in space grid
      Parameter (NZ = 200)
      Parameter (NZIMC = 200)

!     NUD = Number of unknowns in 2 epithelial compartments + potential in lumen
	  Parameter (NUA= 5+2*NS2)
	  Parameter (NUPT= NUA)

!     NUC = Number of unknowns in all 4 epithelial compartments + potential in lumen
      Parameter (NUC = 9+4*NS2) 

!     NDD = number of unknowns in lumen and 2 epithelial compartments
      Parameter (NDA = 7+3*NS2) 
      Parameter (NDPT = NDA)
 
!     NDC = number of unknowns in lumen and 4 epithelial compartments
      Parameter (NDC = 11+5*NS2)
      Parameter (NDIMC = 7+3*NS2)

!     Nneph = number of nephrons in the rat
      Parameter (Nneph = 36000)

!     Natp = Number of variables for H-K-ATPase
      Parameter (Natp = 14)
  
!     Nkcc2 = Number of variables for NKCC2
      Parameter (Nkcc2 = 15)

!     types of fluxes to record (NaK, Hase, Na-trans, Na-para)
      Parameter (NFR = 4)
 
!  PHYSICAL PARAMETERS
      Parameter (PI = 3.14159265359d0)
      Parameter (visc = 6.4d-6)

      Parameter (DiamPT = 0.00200d0)
      Parameter (dimLPT = 0.880d0)
      Parameter (complPT = 0.1)
      Parameter (PMinitPT = 15.00d0)
      Parameter (PbloodPT = 9.0d0)
      Parameter (xS3 = 0.88d0)

      Parameter (DiamD = 0.001275d0)
      Parameter (dimLD = 0.085d0)
      Parameter (SRvol = 0.185d0)
      Parameter (complD = 0.1)
      Parameter (PMinit = 10.0d0)

      Parameter (DiamC = 0.001530d0)
      Parameter (dimLC = 0.170d0)
      Parameter (complC = 0.3)

      Parameter (DiamCCD = 0.002125d0)
      Parameter (dimLCCD = 0.170d0)
      Parameter (complCCD = 0.3)

      Parameter (DiamOMC = 0.002125d0)
      Parameter (dimLOMC = 0.170d0)
      Parameter (complOMC = 0.1)

      Parameter (DiamIMC = 0.00238d0)
      Parameter (dimLIMC = 0.4250d0)
      Parameter (complIMC = 0.1)

      Parameter (DiamT = 0.00170d0)
      Parameter (dimLT = 0.170d0)
      Parameter (complT = 0.1)
      Parameter (DiamMD= 0.00170d0)
      Parameter (dimMD = 0.170d0)
      Parameter (complMD = 0.1)

      Parameter (DiamA = 0.00170d0)
      Parameter (dimLA = 0.170d0)
      Parameter (complA = 0.1)
  
      Parameter (DiamSDL = 0.00170d0)
      Parameter (dimLSDL = 0.1190d0)

! BUFFER AND IMPERMEANT PROPERTIES	

      Parameter (zPimp = -1.0d0)

      Parameter (CPimprefPT = 60.0d0)
      Parameter (CAimprefPT = 60.0d0)
      Parameter (CBimprefPT = 60.0d0)
      Parameter (zPimpPT = -1.0d0)
      Parameter (CPbuftotPT = 60.0d0)

      Parameter (CPimprefSDL = 200.0d0)
      Parameter (CPbuftot = 55.0d0)

      Parameter (CPimprefA = 200.0d0)
	  Parameter (CAimprefA = 200.0d0)
      Parameter (CBimprefA = 200.0d0)      
	  Parameter (zPimpA = -1.0d0)
      Parameter (CPbuftotA = 55.0d0)

      Parameter (CPimprefT = 100.0d0)
	  Parameter (CAimprefT = 100.0d0)
      Parameter (CBimprefT = 100.0d0)
      Parameter (zPimpT = -1.0d0)
      Parameter (CPbuftotT = 55.0d0)

      Parameter (CPimprefMD = 100.0d0)
      Parameter (zPimpMD = -1.0d0)
      Parameter (CPbuftotMD = 55.0d0)

      Parameter (CPimprefD = 70.0d0)
	  Parameter (CAimprefD = 70.0d0)
      Parameter (CBimprefD = 70.0d0)
      Parameter (zPimpD = -1.0d0)
      Parameter (CPbuftotD = 40.0d0)
		
      Parameter (CPimprefC = 50.0d0)
      Parameter (CAimprefC = 18.0d0)
      Parameter (CBimprefC = 18.0d0)
      Parameter (zPimpC = -1.0d0)
      Parameter (zAimpC = -1.0d0)
      Parameter (zBimpC = -1.0d0)
      Parameter (CPbuftotC = 32.0d0)
      Parameter (CAbuftotC = 40.0d0)
      Parameter (CBbuftotC = 40.0d0)

      Parameter (CPimprefCCD = 50.0d0)
      Parameter (CAimprefCCD = 18.0d0)
      Parameter (CBimprefCCD = 18.0d0)
      Parameter (zPimpCCD = -1.0d0)
      Parameter (zAimpCCD = -1.0d0)
      Parameter (zBimpCCD = -1.0d0)
      Parameter (CPbuftotCCD = 32.0d0)
      Parameter (CAbuftotCCD = 40.0d0)
      Parameter (CBbuftotCCD = 40.0d0)
  
      Parameter (CPimprefOMC = 50.0d0)
      Parameter (CAimprefOMC = 60.0d0)
	  Parameter (CBimprefOMC = 60.0d0)
      Parameter (zPimpOMC = -1.0d0)
      Parameter (zAimpOMC = -1.0d0)
      Parameter (zBimpOMC = -1.0d0)
      Parameter (CPbuftotOMC = 32.0d0)
      Parameter (CAbuftotOMC = 40.0d0)
      Parameter (CBbuftotOMC = 40.0d0)
  
      Parameter (CPimprefIMC = 50.0d0)
      Parameter (CAimprefIMC = 50.0d0)
      Parameter (CBimprefIMC = 50.0d0)
      Parameter (zPimpIMC = -1.0d0)
      Parameter (zAimpIMC = -1.0d0)
      Parameter (zBimpIMC = -1.0d0)
      Parameter (CPbuftotIMC = 32.0d0)
      Parameter (CAbuftotIMC = 32.0d0)
      Parameter (CBbuftotIMC = 32.0d0)


! REFERENCE concentrations, volume, permeabilities, membrane potentials
! for non-dimensionalization
! Since Cref = 1 mM, concentrations in mM are already non-dimensional

  parameter (Cref = 1.d-3) !Equals 1 mM = 0.001 mmol/cm3
  parameter (Vref = 1.d-4) !In cm3/cm2 epith
  parameter (Pfref = 1.d0) !In cm/s
  parameter (href = 1.d-5) !In cm/s
  parameter (EPref = 1.0d-3) !In Volts

!---------------------------------------------------------------------72
!     Vwbar = Molar volume of water  (cm3/mmole)
!     Formula weight of H20 = 2*1.00797 + 15.9994 = 18.0153
!     There is 1 mole H20/18.0153 gm of H20.
!     Density of water at 37 C = 0.99335 gm/cm3 (CRC Handbook)
!     Vwbar = (gm/millimole)/(gm/cm3)=(cm3/millimole)

  parameter (Vwbar = (18.0153d0/1000.0d0)/0.99335d0)

! OTHER PARAMETERS
      Parameter (RT = 2.57d0)			
      Parameter (RTosm = 1.93d4)
      Parameter (F = 96.5d0)
      Parameter (pKHCO3 = 3.57d0)
      Parameter (pKHPO4 = 6.80d0)
      Parameter (pKNH3 = 9.15d0)
      Parameter (pKbuf = 7.5d0)
      Parameter (pKHCO2 = 3.76d0)

!  TORQUE PARAMETERS for proximal tubule
      Parameter (torqvm = 0.030d0)
      Parameter (torqL = 2.50d-4)
      Parameter (torqd = 1.50d-5)
!      Parameter (torqR = 0.00106d0)
!	  Parameter (torqR = 0.000910d0)
	  Parameter (torqR = 0.0008475d0)
      Parameter (TS = 1.40d0)

! TRANSPORT PARAMETERS FOR NCC
      Parameter (popncc = 4.295d6)
      Parameter (poppncc = 0.100d6)
      Parameter (pnpncc = 7692.d0)
      Parameter (pnppncc = 179.d0)
      Parameter (dKnncc = 0.293d0)
      Parameter (dKcncc = 112.7d0)
      Parameter (dKncncc = 0.565d0)
      Parameter (dKcnncc = 1.47d-3)

! TRANSPORTER PARAMETERS
      !Parameter (dmuATPH = 2.1d0)
	  Parameter (dmuATPH = 1.450d0)
      Parameter (steepA = 0.40d0)
      Parameter (steepB = 0.40d0)
      Parameter (nNDBCEa = 0)
      Parameter (nNDBCEb = 1)
      Parameter (nAE4a = 0)
      Parameter (nAE4b = 1)
      Parameter (nstoch = 3)
  
! TRANSPORTER PARAMETERS FOR AE1
      Parameter (Pbp = 1247.0d0)
      Parameter (Pbpp = 135.0d0)
      Parameter (Pcp = 562.0d0)
      Parameter (Pcpp = 61.0d0)
      Parameter (dKbp = 198.0d0)
      Parameter (dKbpp = 198.0d0)
      Parameter (dKcp = 50.0d0)
      Parameter (dKcpp = 50.0d0)
  
! TRANSPORTER PARAMETERS FOR PENDRIN
      Parameter (Pclepd = 10000.0d0)
      Parameter (Pcliclepd = 1.239d0)
      Parameter (Pbieclepd = 10.76d0)
      Parameter (Poheclepd = 0.262d0)
      Parameter (dKclpd = 3.01d0)
      Parameter (dKbipd = 5.94d0)
      Parameter (dKohpd = 1.38d-6)

! TRANSPORTER PARAMETERS FOR NDBCE
      Parameter (pBCE = 10000.0d0)
      Parameter (dKnaiBCE = 20.d0)
      Parameter (dKnaeBCE = 20.0d0)
      Parameter (dKcliBCE = 20.d0)
      Parameter (dKcleBCE = 20.0d0)
      Parameter (dKbiiBCE = 15.d0)
      Parameter (dKbieBCE = 15.0d0)

! TRANSPORT PARAMETERS FOR NHE1 
      Parameter (dknhe1na = 15.0d0)
      Parameter (dknhe1h = 1.7d-5)
      Parameter (dknhe1l = 3.6d-3)
      Parameter (dlo = 5952.d0)

! TRANSPORTER PARAMETERS FOR KCC (KCC4 isoform)
      Parameter (poppkcc = 39280.d0)
      Parameter (popkcc = 357700.d0)
      Parameter (pkccp = 10000.0d0)
      Parameter (pkccpp = 1098.d0)
      Parameter (pmccp = 2000.d0)
      Parameter (pmccpp = 219.6d0)
      Parameter (bckcc = 21.08d0)
      Parameter (bkkcc = 1.45d0)
      Parameter (bmkcc = 1.45d0)

! TRANSPORTER PARAMETERS FOR NKCC2 (F isoform)
      Parameter (poppnkccF = 39280.d0)
      Parameter (popnkccF = 357800.d0)
      Parameter (pnkccpF = 10000.0d0)
      Parameter (pnkccppF = 1098.d0)
      Parameter (pnmccpF = 2000.d0)
      Parameter (pnmccppF = 219.6d0)
      Parameter (bn2F = 58.93d0)
      Parameter (bc2F = 13.12d0)
      Parameter (bk2F = 9.149d0)
      Parameter (bm2F = 9.149d0)

! TRANSPORTER PARAMETERS FOR NKCC2 (A isoform)
      Parameter (poppnkccA = 75350.d0)
      Parameter (popnkccA = 259400.d0)
      Parameter (pnkccpA = 10000.0d0)
      Parameter (pnkccppA = 2904.d0)
      Parameter (pnmccpA = 2000.d0)
      Parameter (pnmccppA = 580.8d0)
      Parameter (bn2A = 118.8d0)
      Parameter (bc2A = 0.08834d0)
      Parameter (bk2A = 18710.d0)
      Parameter (bm2A = 18710.d0)

! TRANSPORTER PARAMETERS FOR NKCC2 (B isoform)
      Parameter (poppnkccB = 251700.d0)
      Parameter (popnkccB = 259600.d0)
      Parameter (pnkccpB = 10000.0d0)
      Parameter (pnkccppB = 9695.d0)
      Parameter (pnmccpB = 2000.d0)
      Parameter (pnmccppB = 1939.d0)
      Parameter (bn2B = 275.0d0)
      Parameter (bc2B = 0.08157d0)
      Parameter (bk2B = 5577.d0)
      Parameter (bm2B = 5577.d0)

! TRANSPORTER PARAMETERS FOR NHE
      Parameter (PaNH = 8000.d0) 
      Parameter (PbNH = 8000.d0*0.48d0/1.60d0)
      Parameter (PcNH = 8000.d0)
      Parameter (dKaNH = 30.d0)
      Parameter (dKbNH = 72.0d-6)
      Parameter (dKcNH = 27.d0)
      Parameter (fMNH = 2.0d0)
      Parameter (dKINH = 1.0d-3)
  
! TRANSPORT PARAMETERS FOR H-ATPase
      Parameter (dmuATPHPT = 1.45d0)  ! for proximal tubule
      Parameter (steepATPH = 0.40d0)
      Parameter (dmuHATP = 2.1d0)
      Parameter (steepHATP = 0.40d0)

! TRANSPORTER PARAMETERS FOR NCX
	  Parameter (dKmNao = 87.50d0)
	  Parameter (dKmNai = 12.29d0)
	  Parameter (dKmCao = 1.30d0)
	  Parameter (dKmCai = 3.59d-3)
	  Parameter (gamma_ncx = 0.35d0)
	  Parameter (dKsat_ncx = 0.27d0)
	  Parameter (dKm_ncx = 0.125d-3)

! TRANSPORTER PARAMETERS FOR PMCA
	  Parameter (dKmPMCA = 42.6d-6)

! TRANSPORTER PARAMETERS FOR TRPV5
	  Parameter (Cinhib_v5 = 74.0d-6)

! PARAMETERS FOR CASR
	  Parameter (EC50 = 1.25d0)


