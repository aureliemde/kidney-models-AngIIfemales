!  Define structure for membrane
TYPE membrane

  !!! key predicted variables
  ! solute concentrations
  double precision :: conc(NSPT,NC)
  ! pH
  double precision :: ph(NC)
  ! volume
  double precision :: vol(NC)
  ! membrane potential
  double precision :: ep(NC)
  ! fluid pressure
  double precision :: pres

  ! initial volumes
  double precision :: volEinit, volPinit, volAinit, volBinit, volLuminit
  ! initial area
  double precision :: sbasEinit

  !!! tubular characteristics
  ! tubule length
  double precision :: dimL
  ! luminal diameter
  double precision :: diam
  ! membrane surface area  
  double precision :: area(NC,NC)
  ! water permeability
  double precision :: dLPV(NC,NC)
  ! solute reflection coefficient
  double precision :: sig(NS,NC,NC)
  ! solute permeability
  double precision :: h(NS,NC,NC)
  ! NET coefficient
  double precision :: dLA(NS,NS,NC,NC)
  ! maximum ATPase flux
  double precision :: ATPNaK(NC,NC), ATPH(NC,NC), ATPHK(NC,NC)
  ! NHE3 expression
  double precision :: xNHE3, xNHE1(NC)
  ! SGLT2/SGLT1 expression parameters
  double precision :: xSGLT2, xSGLT1, CTsglt1, CTsglt2
  ! GLUT2/GLUT1 expressions
  double precision :: xGLUT2, xGLUT1, CTglut1, CTglut2
  ! NKCC2 expression
  double precision :: xNKCC2A, xNKCC2B, xNKCC2F
  ! NCC expression
  double precision :: xNCC
  ! KCC4 expression
  double precision :: xKCC4, xKCC4A, xKCC4B
  ! CO2/HCO3/H2CO3 reaction rates
  double precision :: dkd(NC),dkh(NC)
  ! rate of ammoniagenesis
  double precision :: qnh4
  ! AE1 exchanger at peritubular membrane of alpha cell
  double precision :: xAE1
  ! Pendrin exchanger at apical membrane of beta cell
  double precision :: xPendrin
  ! NDBCE exchanger at MB interface
  double precision :: xNDBCE
  ! NCX exchanger on basolateral membranes
  double precision :: xNCX
  ! maximum flux across PMCA pump on basolateral membranes
  double precision :: PMCA
  ! various maximum permeabilities 
  double precision :: hNaMP, hCLCA, hCLCB

  ! total buffer concentrations
  double precision :: cPbuftot, cAbuftot, cBbuftot
  ! reference impermeant celullar concentrations 
  double precision :: cPimpref, cAimpref, cBimpref
  ! impermeant valence
  double precision :: zPimp, zAimp, zBimp
		    
  ! parameter to scale torque effect
  double precision :: scaleT
  double precision :: TM0  ! reference torque

  ! parameter for coalescing tubules
  double precision :: coalesce

  !!! record for keeping tracking of fluxes
  double precision :: FNatrans, FNapara, FNaK, FHase
  double precision :: FGluPara,FGluSGLT1,FGluSGLT2
  double precision :: FKtrans, FKpara

  !!! calculations for metabolism
  double precision :: TQ
  double precision :: nephronAct, nephronTNa, nephronQO2
  
END TYPE membrane
