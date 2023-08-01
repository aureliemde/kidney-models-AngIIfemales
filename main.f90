!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!
!  File:        main.f
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!  Description: This is the main module for the entire nephron.
!  This model yields values of concentrations, volumes, and EP
!   in the lumen and epithelial compartments at equilibrium.

!---------------------------------------------------------------------72
!
!  The units are those of the cgs system: cm, mmol, mmol/cm3 = M, s.
!  The morphological parameters are those of the rat.
!
!---------------------------------------------------------------------72
!
!   Solute Indices:
!    1 = Na+, 2 = K+, 3 = Cl-, 4 = HCO3-, 5 = H2CO3, 6 = CO2
!    7 = HPO4(2-), 8 = H2PO4-, 9 = urea, 10 = NH3, 11 = NH4+, 12 = H+
!   13 = HCO2(-), 14 = H2CO2, 15 = glucose, 16 = Ca2+
!---------------------------------------------------------------------72


Program main

  include 'values.h'
  include 'global.h'
  include 'defs.h'


  type (membrane) :: pt(0:NZ), sdl(0:NZ), mtal(0:NZ), ctal(0:NZ), dct(0:NZ)
  type (membrane) :: cnt(0:NZ), ccd(0:NZ), omcd(0:NZ), imcd(0:NZIMC)

  integer :: ncompl  ! compliant PT or not
  integer :: ntorq   ! accounts for torque effects on transport or not
  integer :: niter   ! iteration number for TGF

! variables used in printing output
  double precision inlet(NS),outlet(NS),fluxs(NS)
  double precision deliv(NS),fracdel(NS)
  double precision :: totalAct, totalTNa, o2consum
  double precision :: nephronAct, nephronTNa, nephronQO2

  ! indices
  integer I,jz,K

!  open (unit=131, file='Flowdata_12F_gen')
!  open (unit=132, file='Flowdata_12F_gen.dat')
!  open (unit=133, file='o2stats_12F_gen')
!  open (unit=134, file='o2stats_12F_gen.dat')
!  open (unit=135, file='Na&Kfluxes_12F_gen')
!  open (unit=136, file='Na&Kfluxes_12F_gen.dat')

  open (unit=131, file='Flowdata_AngII_o')
  open (unit=132, file='Flowdata_AngII_o.dat')
  open (unit=133, file='o2stats_AngII_o')
  open (unit=134, file='o2stats_AngII_o.dat')
  open (unit=135, file='Na&Kfluxes_AngII_o')
  open (unit=136, file='Na&Kfluxes_AngII_o.dat')

!  Scaling factors for volume flows, which are normalized by Vref
!  Assume that there are 36,000 nephrons that coalesce to 7,200 CDs
cw = Vref*60.d0*Nneph ! to convert cm3/s to ml/min with Nneph nephrons
!  cw = Vref*60.d0*1d6 ! If basis is one nephron

  ! initalize metabolic counters to zero
  nephronAct = 0.0 ! Quantifies active sodium reabsorption
  nephronTNa = 0.0 ! Quantifies total sodium reabsorption
  nephronQO2 = 0.0 ! Quantifies total O2 consumption
  nephronTK = 0.0  ! Quantifies total potassium reabsorption

  ! Specify AngII status
  FAng = .true.
  FAngProx = .true.
  FAngDist = .true.
  FAngAqp2 =  .true.

  ! Specify single nephron filtration rate (SNGFR)
  sngfr0 = 1.00*24.0d-6/60.0d0 ! 24 nl/min converted to cm3/s
  sngfr = sngfr0

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!   DETERMINE PROFILES ALONG THE PROXIMAL TUBULE
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  print*,"****************************************************"
  print*,"                    PT RESULTS"
  print*,"****************************************************"
  print*

  ! initialize flags
  ncompl = 1  ! compliant PT (>0)
  ntorq = 1   ! accounts for torque effects on transport (>0)

  ! initialize tubular transport parameters and inflow conditions
  call initPT (pt)

  ! first solve for the first cell: given BC, need only solve for values in cells
  LzPT = -1
  call qnewton1PT (pt(0),CPimprefPT,CPbuftotPT,0,ncompl,ntorq)


  ! then the rest of the PT: solve for values in lumen and cells
  LzPT = 0
  do J = 1, NZ
     call qnewton2PT (pt(J-1),pt(J),J-1,0,pt(0)%vol,CPimprefPT,CPbuftotPT,NDPT,ncompl,ntorq)
     LzPT = LzPt+1
  End Do

  ! output results
  deliv(:) = pt(0)%conc(:,1)*pt(0)%vol(1)*cw
  call out_data_PT (pt,nephronAct,nephronTNa,nephronQO2)
  volreab = (pt(0)%vol(1)-pt(NZ)%vol(1))*cw
  Do I = 1, NS
     outlet(I) = pt(NZ)%vol(1)*pt(NZ)%conc(I,1)*cw
     inlet(I) = pt(0)%vol(1)*pt(0)%conc(I,1)*cw
     fluxs(I) = (inlet(I) - outlet(I))
  !     write(*,'(a,4g12.4)')," Solute",I,inlet(I),outlet(I),fluxs(I)
  end do

  write(*,'(a,4g12.4)'),"H2O  in/reab/out",pt(0)%vol(1)*cw*1.0d3, &
    -volreab*1.0d3,pt(NZ)%vol(1)*cw*1.0d3
  write(*,'(a,4g12.4)'),"Na+  in/reab/out",inlet(1),-fluxs(1),outlet(1)
  write(*,'(a,4g12.4)'),"K+   in/reab/out",inlet(2),-fluxs(2),outlet(2)
  write(*,'(a,4g12.4)'),"Cl-  in/reab/out",inlet(3),-fluxs(3),outlet(3)

  write(131,'(a,4g12.4)'),"PT",pt(0)%vol(1)*cw*1.0d3,inlet(1),inlet(2)
  write(131,'(a,4g12.4)'),"PT",volreab*1.0d3,fluxs(1),fluxs(2)
  write(132,'(4g12.4)'),pt(0)%vol(1)*cw*1.0d3,inlet(1),inlet(2)


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!   DETERMINE PROFILES ALONG THE SDL
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  print*
  print*,"****************************************************"
  print*,"                    SDL RESULTS"
  print*,"****************************************************"
  print*

!   Initialization routine

  call initSDL (sdl)

! If starting from PT, replace luminal concentrations with PT outflow values
   sdl(0)%conc(:,1) = pt(NZ)%conc(:,1)
   sdl(0)%ph(1) = pt(NZ)%ph(1)
   sdl(0)%ep(1) = pt(NZ)%ep(1)
   sdl(0)%vol(1) = pt(NZ)%vol(1)
   sdl(0)%pres = pt(NZ)%pres

!     Determine luminal flow and concentrations along the SDL
  do jz = 1, NZ
     call qnewton2SDL (sdl(jz-1),sdl(jz),jz-1,2,sdl(0)%vol,NDA )
  End Do

!---------------------------------------------------------------------72
!     Output solute flux
!---------------------------------------------------------------------72


  print*,"****************************************************"
  print*
  volreab = (sdl(0)%vol(1)-sdl(NZ)%vol(1))*cw
  Do I = 1, NS
     outlet(I) = sdl(NZ)%vol(1)*sdl(NZ)%conc(I,1)*cw
     inlet(I) = sdl(0)%vol(1)*sdl(0)%conc(I,1)*cw
     fluxs(I) = inlet(I) - outlet(I)
  !     write(*,'(a,4g12.4)')," Solute",I,inlet(I),outlet(I),fluxs(I)
  end do
  print*



!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!   DETERMINE PROFILES ALONG THE mTAL
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  print*
  print*,"****************************************************"
  print*,"                    mTAL RESULTS"
  print*,"****************************************************"
  print*

!   Initialization routine

  call initA (mtal)

! Luminal concentrations are set to SDL outflow values in initA
! Initial guesses for cellular and LIS concentrations are in mTALresults

  open ( unit=31, file='mTALresults')
  Do I = 1,4
     read(31,200),mtal(0)%conc(I,2),mtal(0)%conc(I,5)
  end Do
  Do I = 5,NS
     read(31,210),mtal(0)%conc(I,2),mtal(0)%conc(I,5)
  end Do
  read(31,210),mtal(0)%ph(2),mtal(0)%ph(5)
  read(31,210),mtal(0)%vol(2),mtal(0)%vol(5)
  read(31,210),mtal(0)%ep(2),mtal(0)%ep(5)
  close ( unit=31 )

!     Determine luminal flow and concentrations  along the mTAL

  LzA = 0  ! current (as opposed to next) position
  do jz = 1, NZ
     call qnewton2b (mtal(jz-1),mtal(jz),mtal(0)%vol,jz-1,3,NDA)
     LzA = LzA + 1
  End Do

  ! output first cell solution
  open ( unit=31, file='mTALresults')
  Do I = 1,4
     write(31,200),mtal(1)%conc(I,2),mtal(1)%conc(I,5)
  end Do
  Do I = 5,NS
     write(31,210),mtal(1)%conc(I,2),mtal(1)%conc(I,5)
  end Do
  write(31,210),mtal(1)%ph(2),mtal(1)%ph(5)
  write(31,210),mtal(1)%vol(2),mtal(1)%vol(5)
  write(31,210),mtal(1)%ep(2),mtal(1)%ep(5)
  close ( unit=31 )

!---------------------------------------------------------------------72
!     Output solute flux
!---------------------------------------------------------------------72

  print*,"****************************************************"
  print*
  volreab = (mtal(0)%vol(1)-mtal(NZ)%vol(1))*cw
  Do I = 1, NS
     outlet(I) = mtal(NZ)%vol(1)*mtal(NZ)%conc(I,1)*cw
     inlet(I) = mtal(0)%vol(1)*mtal(0)%conc(I,1)*cw
     fluxs(I) = inlet(I) - outlet(I)
  !     write(*,'(a,4g12.4)')," Solute",I,inlet(I),outlet(I),fluxs(I)
  end do
  print*

  call compute_o2_consumption (mtal,'mTAL    ',dimLA,1,NZ,nephronAct,nephronTNa,nephronQO2)
  write(*,'(a,4g12.4)'),"H2O  in/reab/out",mtal(0)%vol(1)*cw*1.0d3, &
    -volreab*1.0d3,mtal(NZ)%vol(1)*cw*1.0d3
  write(*,'(a,4g12.4)'),"Na+  in/reab/out",inlet(1),-fluxs(1),outlet(1)
  write(*,'(a,4g12.4)'),"K+   in/reab/out",inlet(2),-fluxs(2),outlet(2)
  write(*,'(a,4g12.4)'),"Cl-  in/reab/out",inlet(3),-fluxs(3),outlet(3)

  write(131,'(a,4g12.4)'),"mTAL",mtal(0)%vol(1)*cw*1.0d3,inlet(1),inlet(2)
  write(131,'(a,4g12.4)'),"mTAL",volreab*1.0d3,fluxs(1),fluxs(2)
  write(132,'(4g12.4)'),mtal(0)%vol(1)*cw*1.0d3,inlet(1),inlet(2)


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!   DETERMINE PROFILES ALONG THE cTAL
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  print*
  print*,"****************************************************"
  print*,"                    cTAL RESULTS"
  print*,"****************************************************"
  print*

!   Initialization routine

  call initT (ctal)

! Luminal concentrations are set to mTAL outflow values in initT
! Initial guesses for cellular and LIS concentrations are in cTALresults


  open ( unit=41, file='cTALresults')
  Do I = 1,4
     read(41,200),ctal(0)%conc(I,2),ctal(0)%conc(I,5)
  end Do
  Do I = 5,NS
     read(41,210),ctal(0)%conc(I,2),ctal(0)%conc(I,5)
  end Do
  read(41,210),ctal(0)%ph(2),ctal(0)%ph(5)
  read(41,210),ctal(0)%vol(2),ctal(0)%vol(5)
  read(41,210),ctal(0)%ep(2),ctal(0)%ep(5)
  close ( unit=41 )

!     Determine luminal flow and concentrations along the cTAL
  LzT = 0
  do jz = 1, NZ
     call qnewton2b (ctal(jz-1),ctal(jz),ctal(0)%vol,LzT,4,NDA )
     LzT = LzT + 1
  End Do

  ! output first cell solution
  open ( unit=41, file='cTALresults')
  Do I = 1,4
     write(41,200),ctal(1)%conc(I,2),ctal(1)%conc(I,5)
  end Do
  Do I = 5,NS
     write(41,210),ctal(1)%conc(I,2),ctal(1)%conc(I,5)
  end Do
  write(41,210),ctal(1)%ph(2),ctal(1)%ph(5)
  write(41,210),ctal(1)%vol(2),ctal(1)%vol(5)
  write(41,210),ctal(1)%ep(2),ctal(1)%ep(5)
  close ( unit=41 )

!---------------------------------------------------------------------72
!     Output solute flux
!---------------------------------------------------------------------72


  print*,"****************************************************"
  print*
  volreab = (ctal(0)%vol(1)-ctal(NZ)%vol(1))*cw
!  print*,"Solute in/out and reabsorbed (umol/min)"
  Do I = 1, NS
     outlet(I) = ctal(NZ)%vol(1)*ctal(NZ)%conc(I,1)*cw
     inlet(I) = ctal(0)%vol(1)*ctal(0)%conc(I,1)*cw
     fluxs(I) = inlet(I) - outlet(I)
!     write(*,'(a,4g12.4)')," Solute",I,inlet(I),outlet(I),fluxs(I)
  end do
  print*

  call compute_o2_consumption (ctal,'cTAL    ',dimLT,1,NZ,nephronAct,nephronTNa,nephronQO2)
  write(*,'(a,4g12.4)'),"H2O  in/reab/out",ctal(0)%vol(1)*cw*1.0d3, &
    -volreab*1.0d3,ctal(NZ)%vol(1)*cw*1.0d3
  write(*,'(a,4g12.4)'),"Na+  in/reab/out",inlet(1),-fluxs(1),outlet(1)
  write(*,'(a,4g12.4)'),"K+   in/reab/out",inlet(2),-fluxs(2),outlet(2)
  write(*,'(a,4g12.4)'),"Cl-  in/reab/out",inlet(3),-fluxs(3),outlet(3)

  write(131,'(a,4g12.4)'),"cTAL",ctal(0)%vol(1)*cw*1.0d3,inlet(1),inlet(2)
  write(131,'(a,4g12.4)'),"cTAL",volreab*1.0d3,fluxs(1),fluxs(2)
  write(132,'(4g12.4)'),ctal(0)%vol(1)*cw*1.0d3,inlet(1),inlet(2)


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!   DETERMINE PROFILES ALONG THE DCT
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

   print*
   print*,"****************************************************"
   print*,"                    DCT RESULTS"
   print*,"****************************************************"
   print*

   call initD (dct)

! Luminal concentrations are set to cTAL outflow values in initD
! Initial guesses for cellular and LIS concentrations are in DCTresults

   open ( unit=51, file='DCTresults')

   Do I = 1,4
     read(51,200),dct(0)%conc(I,2),dct(0)%conc(I,5)
   end Do
   Do I = 5,NS
     read(51,210),dct(0)%conc(I,2),dct(0)%conc(I,5)
   end Do
   read(51,210),dct(0)%ph(2),dct(0)%ph(5)
   read(51,210),dct(0)%vol(2),dct(0)%vol(5)
   read(51,210),dct(0)%ep(2),dct(0)%ep(5)
   close ( unit=51 )

!     Determine luminal flow and concentrations along the DCT
   LzD = 0
   do jz = 1, NZ
      call qnewton2b (dct(jz-1),dct(jz),dct(0)%vol,LzD,5,NDA )
      LzD = LzD + 1
   End Do

   ! output first cell solution
  open ( unit=51, file='DCTresults')
  Do I = 1,4
     write(51,200),dct(1)%conc(I,2),dct(1)%conc(I,5)
  end Do
  Do I = 5,NS
     write(51,210),dct(1)%conc(I,2),dct(1)%conc(I,5)
  end Do
  write(51,210),dct(1)%ph(2),dct(1)%ph(5)
  write(51,210),dct(1)%vol(2),dct(1)%vol(5)
  write(51,210),dct(1)%ep(2),dct(1)%ep(5)
  close ( unit=51 )

!---------------------------------------------------------------------72
!     Output solute flux
!---------------------------------------------------------------------72

   print*,"****************************************************"
   print*
   volreab = (dct(0)%vol(1)-dct(NZ)%vol(1))*cw
   Do I = 1, NS
      outlet(I) = dct(NZ)%vol(1)*dct(NZ)%conc(I,1)*cw
      inlet(I) = dct(0)%vol(1)*dct(0)%conc(I,1)*cw
      fluxs(I) = inlet(I) - outlet(I)
 !     write(*,'(a,4g12.4)')," Solute",I,inlet(I),outlet(I),fluxs(I)
   end do
   print*

   call compute_o2_consumption (dct,'DCT     ',dimLD,1,NZ,nephronAct,nephronTNa,nephronQO2)
   write(*,'(a,4g12.4)'),"H2O  in/reab/out",dct(0)%vol(1)*cw*1.0d3, &
    -volreab*1.0d3,dct(NZ)%vol(1)*cw*1.0d3
  write(*,'(a,4g12.4)'),"Na+  in/reab/out",inlet(1),-fluxs(1),outlet(1)
  write(*,'(a,4g12.4)'),"K+   in/reab/out",inlet(2),-fluxs(2),outlet(2)
  write(*,'(a,4g12.4)'),"Cl-  in/reab/out",inlet(3),-fluxs(3),outlet(3)

  write(131,'(a,4g12.4)'),"DCT",dct(0)%vol(1)*cw*1.0d3,inlet(1),inlet(2)
  write(131,'(a,4g12.4)'),"DCT",volreab*1.0d3,fluxs(1),fluxs(2)
  write(132,'(4g12.4)'),dct(0)%vol(1)*cw*1.0d3,inlet(1),inlet(2)


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!   DETERMINE PROFILES ALONG THE CNT
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

   print*
   print*,"****************************************************"
   print*,"                   CNT RESULTS"
   print*,"****************************************************"
   print*

   call initC (cnt)

! Luminal concentrations are set to DCT outflow values in initC
! Initial guesses for cellular and LIS concentrations are in CNTresults

   open ( unit=61, file='CNTresults')
   Do I = 1,4
      read(61,200),cnt(0)%conc(I,2),cnt(0)%conc(I,3),cnt(0)%conc(I,4),cnt(0)%conc(I,5)
   end Do
   Do I = 5,NS
      read(61,210),cnt(0)%conc(I,2),cnt(0)%conc(I,3),cnt(0)%conc(I,4),cnt(0)%conc(I,5)
   end Do
   read(61,210),cnt(0)%ph(2),cnt(0)%ph(3),cnt(0)%ph(4),cnt(0)%ph(5)
   read(61,210),cnt(0)%vol(2),cnt(0)%vol(3),cnt(0)%vol(4),cnt(0)%vol(5)
   read(61,210),cnt(0)%ep(2),cnt(0)%ep(3),cnt(0)%ep(4),cnt(0)%ep(5)
   close ( unit=61 )

!     Determine luminal flow and concentrations along the CNT
   LzC = 0
   do jz = 1, NZ
      call qnewton2icb (cnt(jz-1),cnt(jz),jz-1, &
           6,cnt(0)%vol, cnt(jz)%volEinit,cnt(jz)%volPinit, &
           cnt(jz)%volAinit,cnt(jz)%volBinit,NDC )
      LzC = LzC+1
   End Do

   ! output first cell solution
   open ( unit=61, file='CNTresults')
   Do I = 1,4
      write(61,200),cnt(1)%conc(I,2),cnt(1)%conc(I,3),cnt(1)%conc(I,4),cnt(1)%conc(I,5)
   end Do
   Do I = 5,NS
      write(61,210),cnt(1)%conc(I,2),cnt(1)%conc(I,3),cnt(1)%conc(I,4),cnt(1)%conc(I,5)
   end Do
   write(61,210),cnt(1)%ph(2),cnt(1)%ph(3),cnt(1)%ph(4),cnt(1)%ph(5)
   write(61,210),cnt(1)%vol(2),cnt(1)%vol(3),cnt(1)%vol(4),cnt(1)%vol(5)
   write(61,210),cnt(1)%ep(2),cnt(1)%ep(3),cnt(1)%ep(4),cnt(1)%ep(5)
   close ( unit=61 )


!---------------------------------------------------------------------72
!     Output solute flux
!---------------------------------------------------------------------72

   print*,"****************************************************"
   print*
   volreab = (cnt(0)%vol(1)-cnt(NZ)%vol(1))*cw
   Do I = 1, NS
      outlet(I) = cnt(NZ)%vol(1)*cnt(NZ)%conc(I,1)*cw
      inlet(I) = cnt(0)%vol(1)*cnt(0)%conc(I,1)*cw
      fluxs(I) = inlet(I) - outlet(I)
 !     write(*,'(a,4g12.4)')," Solute",I,inlet(I),outlet(I),fluxs(I)
   end do
   print*

   call compute_o2_consumption (cnt,'CNT     ',dimLC,1,NZ,nephronAct,nephronTNa,nephronQO2)
   write(*,'(a,4g12.4)'),"H2O  in/reab/out",cnt(0)%vol(1)*cw*1.0d3, &
    -volreab*1.0d3,cnt(NZ)%vol(1)*cw*1.0d3
  write(*,'(a,4g12.4)'),"Na+  in/reab/out",inlet(1),-fluxs(1),outlet(1)
  write(*,'(a,4g12.4)'),"K+   in/reab/out",inlet(2),-fluxs(2),outlet(2)
  write(*,'(a,4g12.4)'),"Cl-  in/reab/out",inlet(3),-fluxs(3),outlet(3)

  write(131,'(a,4g12.4)'),"CNT",cnt(0)%vol(1)*cw*1.0d3,inlet(1),inlet(2)
  write(131,'(a,4g12.4)'),"CNT",volreab*1.0d3,fluxs(1),fluxs(2)
  write(132,'(4g12.4)'),cnt(0)%vol(1)*cw*1.0d3,inlet(1),inlet(2)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!   DETERMINE PROFILES ALONG THE CCD
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

   print*
   print*,"****************************************************"
   print*,"                   CCD RESULTS"
   print*,"****************************************************"
   print*

   call initCCD (ccd)

! Luminal concentrations are set to CNT outflow values in initD
! Initial guesses for cellular and LIS concentrations are in CCDresults

   open ( unit=71, file='CCDresults')
   Do I = 1,4
      read(71,200),ccd(0)%conc(I,2),ccd(0)%conc(I,3),ccd(0)%conc(I,4),ccd(0)%conc(I,5)
   end Do
   Do I = 5,NS
      read(71,210),ccd(0)%conc(I,2),ccd(0)%conc(I,3),ccd(0)%conc(I,4),ccd(0)%conc(I,5)
   end Do
   read(71,210),ccd(0)%ph(2),ccd(0)%ph(3),ccd(0)%ph(4),ccd(0)%ph(5)
   read(71,210),ccd(0)%vol(2),ccd(0)%vol(3),ccd(0)%vol(4),ccd(0)%vol(5)
   read(71,210),ccd(0)%ep(2),ccd(0)%ep(3),ccd(0)%ep(4),ccd(0)%ep(5)
   close ( unit=71 )

!     Determine luminal flow and concentrations along the CCD
   LzCCD = 0
   do jz = 1, NZ
      call qnewton2icb (ccd(jz-1),ccd(jz),jz-1, &
           7 ,ccd(0)%vol, ccd(jz)%volEinit,ccd(jz)%volPinit, &
           ccd(jz)%volAinit,ccd(jz)%volBinit,NDC )
      LzCCD = LzCCD + 1
   End Do

   ! output first cell solution
   open ( unit=71, file='CCDresults')
   Do I = 1,4
      write(71,200),ccd(1)%conc(I,2),ccd(1)%conc(I,3),ccd(1)%conc(I,4),ccd(1)%conc(I,5)
   end Do
   Do I = 5,NS
      write(71,210),ccd(1)%conc(I,2),ccd(1)%conc(I,3),ccd(1)%conc(I,4),ccd(1)%conc(I,5)
   end Do
   write(71,210),ccd(1)%ph(2),ccd(1)%ph(3),ccd(1)%ph(4),ccd(1)%ph(5)
   write(71,210),ccd(1)%vol(2),ccd(1)%vol(3),ccd(1)%vol(4),ccd(1)%vol(5)
   write(71,210),ccd(1)%ep(2),ccd(1)%ep(3),ccd(1)%ep(4),ccd(1)%ep(5)
   close ( unit=71 )

!---------------------------------------------------------------------72
!     Output solute flux
!---------------------------------------------------------------------72

  print*,"****************************************************"
  print*
  volreab = (ccd(0)%vol(1)-ccd(NZ)%vol(1))*cw
  Do I = 1, NS
        outlet(I) = ccd(NZ)%vol(1)*ccd(NZ)%conc(I,1)*cw
        inlet(I) = ccd(0)%vol(1)*ccd(0)%conc(I,1)*cw
        fluxs(I) = inlet(I) - outlet(I)
  !     write(*,'(a,4g12.4)')," Solute",I,inlet(I),outlet(I),fluxs(I)
  end do
  print*

  call compute_o2_consumption (ccd,'CCD     ',dimLCCD,1,NZ,nephronAct,nephronTNa,nephronQO2)
   write(*,'(a,4g12.4)'),"H2O  in/reab/out",ccd(0)%vol(1)*cw*1.0d3, &
    -volreab*1.0d3,ccd(NZ)%vol(1)*cw*1.0d3
  write(*,'(a,4g12.4)'),"Na+  in/reab/out",inlet(1),-fluxs(1),outlet(1)
  write(*,'(a,4g12.4)'),"K+   in/reab/out",inlet(2),-fluxs(2),outlet(2)
  write(*,'(a,4g12.4)'),"Cl-  in/reab/out",inlet(3),-fluxs(3),outlet(3)

  write(131,'(a,4g12.4)'),"CCD",ccd(0)%vol(1)*cw*1.0d3,inlet(1),inlet(2)
  write(131,'(a,4g12.4)'),"CCD",volreab*1.0d3,fluxs(1),fluxs(2)
  write(132,'(4g12.4)'),ccd(0)%vol(1)*cw*1.0d3,inlet(1),inlet(2)


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!   DETERMINE PROFILES ALONG THE OMCD
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

   print*
   print*,"****************************************************"
   print*,"                   OMCD RESULTS"
   print*,"****************************************************"
   print*

   call initOMC (omcd)

! Luminal concentrations are set to CCD outflow values in initOMC
! Initial guesses for cellular and LIS concentrations are in OMCresults

   open ( unit=81, file='OMCresults')
   Do I = 1,4
      read(81,200),omcd(0)%conc(I,2),omcd(0)%conc(I,3),omcd(0)%conc(I,4),omcd(0)%conc(I,5)
   end Do
   Do I = 5,NS
      read(81,210),omcd(0)%conc(I,2),omcd(0)%conc(I,3),omcd(0)%conc(I,4),omcd(0)%conc(I,5)
   end Do
   read(81,210),omcd(0)%ph(2),omcd(0)%ph(3),omcd(0)%ph(4),omcd(0)%ph(5)
   read(81,210),omcd(0)%vol(2),omcd(0)%vol(3),omcd(0)%vol(4),omcd(0)%vol(5)
   read(81,210),omcd(0)%ep(2),omcd(0)%ep(3),omcd(0)%ep(4),omcd(0)%ep(5)
   close ( unit=81 )

!     Determine luminal flow and concentrations along the OMCD
  LzOMC = 0
  do jz = 1, NZ
     call qnewton2icb (omcd(jz-1),omcd(jz),jz-1, &
          8 ,omcd(0)%vol, omcd(jz)%volEinit,omcd(jz)%volPinit, &
          omcd(jz)%volAinit,omcd(jz)%volBinit,NDC )
     LzOMC = LzOMC + 1

  End Do

  ! output first cell
  open ( unit=81, file='OMCresults')
  Do I = 1,4
     write(81,200),omcd(1)%conc(I,2),omcd(1)%conc(I,3),omcd(1)%conc(I,4),omcd(1)%conc(I,5)
  end Do
  Do I = 5,NS
     write(81,210),omcd(1)%conc(I,2),omcd(1)%conc(I,3),omcd(1)%conc(I,4),omcd(1)%conc(I,5)
  end Do
  write(81,210),omcd(1)%ph(2),omcd(1)%ph(3),omcd(1)%ph(4),omcd(1)%ph(5)
  write(81,210),omcd(1)%vol(2),omcd(1)%vol(3),omcd(1)%vol(4),omcd(1)%vol(5)
  write(81,210),omcd(1)%ep(2),omcd(1)%ep(3),omcd(1)%ep(4),omcd(1)%ep(5)
  close ( unit=81 )

!---------------------------------------------------------------------72
!     Output solute flux
!---------------------------------------------------------------------72

  print*,"****************************************************"
  print*
  volreab = (omcd(0)%vol(1)-omcd(NZ)%vol(1))*cw

  Do I = 1, NS
        outlet(I) = omcd(NZ)%vol(1)*omcd(NZ)%conc(I,1)*cw
        inlet(I) = omcd(0)%vol(1)*omcd(0)%conc(I,1)*cw
        fluxs(I) = inlet(I) - outlet(I)
  !     write(*,'(a,4g12.4)')," Solute",I,inlet(I),outlet(I),fluxs(I)
  end do

  call compute_o2_consumption (omcd,'OMCD    ',dimLOMC,1,NZ,nephronAct,nephronTNa,nephronQO2)
   write(*,'(a,4g12.4)'),"H2O  in/reab/out",omcd(0)%vol(1)*cw*1.0d3, &
    -volreab*1.0d3,omcd(NZ)%vol(1)*cw*1.0d3
  write(*,'(a,4g12.4)'),"Na+  in/reab/out",inlet(1),-fluxs(1),outlet(1)
  write(*,'(a,4g12.4)'),"K+   in/reab/out",inlet(2),-fluxs(2),outlet(2)
  write(*,'(a,4g12.4)'),"Cl-  in/reab/out",inlet(3),-fluxs(3),outlet(3)

  write(131,'(a,4g12.4)'),"OMCD",omcd(0)%vol(1)*cw*1.0d3,inlet(1),inlet(2)
  write(131,'(a,4g12.4)'),"OMCD",volreab*1.0d3,fluxs(1),fluxs(2)
  write(132,'(4g12.4)'),omcd(0)%vol(1)*cw*1.0d3,inlet(1),inlet(2)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!   DETERMINE PROFILES ALONG THE IMCD
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

      print*
      print*,"****************************************************"
      print*,"                   IMCD RESULTS"
      print*,"****************************************************"
      print*

   call initIMC (imcd)

! Luminal concentrations are set to CCD outflow values in initIMC
! Initial guesses for cellular and LIS concentrations are in IMCresults

    open ( unit=91, file='IMCresults')
    Do I = 1,4
       read(91,200),imcd(0)%conc(I,2),imcd(0)%conc(I,5)
    end Do
    Do I = 5,NS
       read(91,210),imcd(0)%conc(I,2),imcd(0)%conc(I,5)
    end Do
    read(91,210),imcd(0)%ph(2),imcd(0)%ph(5)
    read(91,210),imcd(0)%vol(2),imcd(0)%vol(5)
    read(91,210),imcd(0)%ep(2),imcd(0)%ep(5)
    close ( unit=91 )

!     Determine luminal flow and concentrations along the IMC
    nw1imc = 1
    LzIMC = 0
    do jz = 1, NZIMC
       call qnewton2b (imcd(jz-1),imcd(jz),imcd(0)%vol,LzIMC,9,NDIMC )
       LzIMC = LzIMC + 1
    End Do


  ! output first cell
  open ( unit=91, file='IMCresults')
  Do I = 1,4
     write(91,200),imcd(1)%conc(I,2),imcd(1)%conc(I,5)
  end Do
  Do I = 5,NS
     write(91,210),imcd(1)%conc(I,2),imcd(1)%conc(I,5)
  end Do
  write(91,210),imcd(1)%ph(2),imcd(1)%ph(5)
  write(91,210),imcd(1)%vol(2),imcd(1)%vol(5)
  write(91,210),imcd(1)%ep(2),imcd(1)%ep(5)
  close ( unit=91 )

!---------------------------------------------------------------------72
!     Output solute flux
!---------------------------------------------------------------------72
  TotalENaC = 0.0d0
  TotalNCC = 0.0d0
  FluxK_trans = imcd(0)%FKtrans/NZIMC
  FluxK_para = imcd(0)%FKpara/NZIMC
  Do jz = 1, NZIMC
    TotalENaC = TotalENaC + fENaC(jz)
    TotalNCC = TotalNCC + fNCC(jz)
    FluxK_trans = FluxK_trans + imcd(jz)%FKtrans/NZ
    FluxK_para = FluxK_para + imcd(jz)%FKpara/NZ
  end do

  print*,"****************************************************"
  print*
  volreab = (imcd(0)%vol(1)-imcd(NZ)%vol(1))*cw
  Urine_osm = 0
  Do I = 1, NS
        outlet(I) = imcd(NZ)%vol(1)*imcd(NZ)%conc(I,1)*cw
        inlet(I) = imcd(0)%vol(1)*imcd(0)%conc(I,1)*cw
        fluxs(I) = inlet(I) - outlet(I)
!        write(*,'(a,4g12.4)')," Solute in/out/reab",I,inlet(I),outlet(I),fluxs(I)
        Urine_osm = Urine_osm + imcd(NZ)%conc(I,1)
  end do
  print*

  call compute_o2_consumption (imcd,'IMCD    ',dimLIMC,1,NZ,nephronAct,nephronTNa,nephronQO2)
  write(133,'(a,4g12.5)'), "Total", nephronAct, nephronTNa, nephronQO2
  write(134,'(4g12.5)'), nephronAct, nephronTNa, nephronQO2

  write(*,'(a,4g12.4)'),"H2O  in/reab/out",imcd(0)%vol(1)*cw*1.0d3, &
    -volreab*1.0d3,imcd(NZIMC)%vol(1)*cw*1.0d3
  write(*,'(a,4g12.4)'),"Na+  in/reab/out",inlet(1),-fluxs(1),outlet(1)
  write(*,'(a,4g12.4)'),"K+   in/reab/out",inlet(2),-fluxs(2),outlet(2)
  write(*,'(a,4g12.4)'),"Cl-  in/reab/out",inlet(3),-fluxs(3),outlet(3)

  write(131,'(a,4g12.4)'),"IMCD",imcd(0)%vol(1)*cw*1.0d3,inlet(1),inlet(2)
  write(131,'(a,4g12.4)'),"IMCD",volreab*1.0d3,fluxs(1),fluxs(2)
  write(132,'(4g12.4)'),imcd(0)%vol(1)*cw*1.0d3,inlet(1),inlet(2)
  write(131,'(a,4g12.4)'),"Urine",imcd(NZIMC)%vol(1)*cw*1.0d3,outlet(1),outlet(2)
  write(132,'(4g12.4)'),imcd(NZIMC)%vol(1)*cw*1.0d3,outlet(1),outlet(2)

!  call output_solute(pt,sdl,mtal,ctal,dct,cnt,ccd,omcd,imcd)


     close (unit=131)
     close (unit=132)
     close (unit=133)
     close (unit=134)
     close (unit=135)
     close (unit=136)

200  format (6f12.5)
210  format (6e12.5)

123  STOP
END Program main

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!     Sodium Transport and Oxygen Consumption
!---------------------------------------------------------------------72

subroutine compute_o2_consumption (tube,tubename,dimL,ind1,ind2,nephronAct,nephronTNa,nephronQO2)

  include 'values.h'
  include 'defs.h'

  ! passed variables
  type(membrane) :: tube(0:NZ)
  character :: tubename(8)
  double precision :: dimL  ! tubule length
  integer :: ind1, ind2  ! beginning and end indices of segment of interest
  double precision :: nephronAct, nephronTNa, nephronTK, nephronQO2

  ! local variables
  double precision :: totalAct, totalTNa, totalHpump, o2consum, totalTK
  integer nHATpase

  ! nHATPase = 1 if we account for ATP consumption by H-ATPase pumps in the PT
  nHATPase = 0

  totalAct = (0.5*(tube(ind1)%FNaK+tube(ind2)%FNaK)+sum(tube(ind1+1:ind2-1)%FNaK))*dimL/(ind2-ind1)
  totalTNa = (0.5*(tube(ind1)%FNatrans+tube(ind2)%FNatrans)+sum(tube(ind1+1:ind2-1)%FNatrans))*dimL/(ind2-ind1)
  totalTNa = totalTNa + (0.5*(tube(ind1)%FNapara+tube(ind2)%FNapara)+sum(tube(ind1+1:ind2-1)%FNapara))*dimL/(ind2-ind1)
  totalHpump = (0.5*(tube(ind1)%FHase+tube(ind2)%FHase)+sum(tube(ind1+1:ind2-1)%FHase))*dimL/(ind2-ind1)

  o2consum = totalAct/tube(0)%TQ + abs(totalHpump)/10.0d0*nHATPase

  totalTK = (0.5*(tube(ind1)%FKtrans+tube(ind2)%FKtrans)+sum(tube(ind1+1:ind2-1)%FKtrans))*dimL/(ind2-ind1)
  totalTK = totalTK+ (0.5*(tube(ind1)%FKpara+tube(ind2)%FKpara)+sum(tube(ind1+1:ind2-1)%FKpara))*dimL/(ind2-ind1)

 write(133,'(a,4g12.5)'),tubename, totalAct*10.0d0, totalTNa*10.0d0, o2consum*10.0d0
!  write(133,*),tubename, totalAct*10.0d0, totalTNa*10.0d0, o2consum*10.0d0
  write(134,'(4g12.5)'),totalAct*10.0d0, totalTNa*10.0d0, o2consum*10.0d0

  write(135,*),tubename, totalTNa*10.0d0, totalTK*10.0d0, o2consum*10.0d0
  write(136,'(4g12.5)'),totalTNa*10.0d0, totalTK*10.0d0, o2consum*10.0d0

  nephronAct = nephronAct + totalAct*10.0d0
  nephronTNa = nephronTNa + totalTNa*10.0d0
  nephronQO2 = nephronQO2 + o2consum*10.0d0
  nephronTK = nephronTK + totalTK*10.0d0

  return
end subroutine compute_o2_consumption





