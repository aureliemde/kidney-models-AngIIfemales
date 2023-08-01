subroutine output_tal_profiles (mtal,ctal)

  include 'values.h'
  include 'global.h'
  include 'defs.h'
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	OUTPUT CONCENTRATION PROFILES ALONG ENTIRE TAL 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  !passed variables
  type (membrane) :: mtal(0:NZ), ctal(0:NZ)
  
      double precision osmol
      double precision cw
      cw = Vref*60.d6 !to convert to nl/min

!     peritubular profiles
      open ( unit=21, file='talNa_p' )
      open ( unit=22, file='talK_p' )
      open ( unit=23, file='talCl_p' )
      open ( unit=24, file='talHCO3_p' )
      open ( unit=25, file='talH2CO3_p' )
      open ( unit=26, file='talCO2_p' )
      open ( unit=27, file='talHPO4_p' )
      open ( unit=28, file='talH2PO4_p' )
      open ( unit=29, file='talurea_p' )
      open ( unit=30, file='talNH3_p' )
      open ( unit=31, file='talNH4_p' )
      open ( unit=32, file='talH_p' )
      open ( unit=33, file='talHCO2_p' )
      open ( unit=34, file='talH2CO2_p' )
      open ( unit=35, file='talGlu_p' )

!     luminal profiles
      open ( unit=41, file='talNa_l' )
      open ( unit=42, file='talK_l' )
      open ( unit=43, file='talCl_l' )
      open ( unit=44, file='talHCO3_l' )
      open ( unit=45, file='talH2CO3_l' )
      open ( unit=46, file='talCO2_l' )
      open ( unit=47, file='talHPO4_l' )
      open ( unit=48, file='talH2PO4_l' )
      open ( unit=49, file='talurea_l' )
      open ( unit=50, file='talNH3_l' )
      open ( unit=51, file='talNH4_l' )
      open ( unit=52, file='talH_l' )
      open ( unit=53, file='talHCO2_l' )
      open ( unit=54, file='talH2CO2_l' )
      open ( unit=55, file='talGlu_l' )

!     cellular profiles
      open ( unit=61, file='talNa_c' )
      open ( unit=62, file='talK_c' )
      open ( unit=63, file='talCl_c' )
      open ( unit=64, file='talHCO3_c' )
      open ( unit=65, file='talH2CO3_c' )
      open ( unit=66, file='talCO2_c' )
      open ( unit=67, file='talHPO4_c' )
      open ( unit=68, file='talH2PO4_c' )
      open ( unit=69, file='talurea_c' )
      open ( unit=70, file='talNH3_c' )
      open ( unit=71, file='talNH4_c' )
      open ( unit=72, file='talH_c' )
      open ( unit=73, file='talHCO2_c' )
      open ( unit=74, file='talH2CO2_c' )
      open ( unit=75, file='talGlu_c' )

      open ( unit=80, file='talVol' )
      open ( unit=81, file='talKflow' )
      open ( unit=82, file='talNaflow' )
      open ( unit=83, file='talClflow' )
      open ( unit=84, file='talNH4flow' )
      open ( unit=85, file='talos_c')
      open ( unit=86, file='talpH_l')
      open ( unit=87, file='talpH_c')
      open ( unit=88, file='talflow')
      open ( unit=89, file='talHCO3flow' )
      open ( unit=90, file='talH2CO3flow' )
      open ( unit=91, file='talHPO4flow' )
      open ( unit=92, file='talH2PO4flow' )
      open ( unit=93, file='talureaflow' )
      open ( unit=94, file='talNH3flow' )
      open ( unit=95, file='talpres' )

!     flux values
!     open ( unit=74, file='fluxNKCC2' )
!     open ( unit=75, file='fluxNHE' )

      dx = dimLA/NZ
      do j = 0,NZ
        pos = dimLPT + dimLSDL + j*dx
        osmol = mtal(j)%conc(NS,2)
        do i = 1, NS
          write(20+i,200), pos,mtal(j)%conc(i,6)
          write(40+i,200), pos,mtal(j)%conc(i,1)
          write(60+i,200), pos,mtal(j)%conc(i,2)
          osmol = osmol + mtal(j)%conc(i,2) 
        end do 
        write(80,200) pos,mtal(j)%vol(2)*cw 
        write(81,200) pos,mtal(j)%vol(1)*mtal(j)%conc(2,1)*cw
        write(82,200) pos,mtal(j)%vol(1)*mtal(j)%conc(1,1)*cw
        write(83,200) pos,mtal(j)%vol(1)*mtal(j)%conc(3,1)*cw
        write(84,200) pos,mtal(j)%vol(1)*mtal(j)%conc(11,1)*cw
        write(85,200) pos,osmol
        write(86,200) pos,-log10(mtal(j)%conc(12,1))+3
        write(87,200) pos,-log10(mtal(j)%conc(12,2))+3
        write(88,200) pos,mtal(j)%vol(1) *cw
        write(89,200) pos,mtal(j)%vol(1)*mtal(j)%conc(4,1)*cw
        write(90,200) pos,mtal(j)%vol(1)*mtal(j)%conc(5,1)*cw
        write(91,200) pos,mtal(j)%vol(1)*mtal(j)%conc(7,1)*cw
        write(92,200) pos,mtal(j)%vol(1)*mtal(j)%conc(8,1)*cw
        write(93,200) pos,mtal(j)%vol(1)*mtal(j)%conc(9,1)*cw
        write(94,200) pos,mtal(j)%vol(1)*mtal(j)%conc(10,1)*cw
        write(95,200) pos,mtal(j)%pres

      end do
      dx = dimLT/NZ
      do j = 1,NZ
        osmol = ctal(j)%conc(NS,2)
        pos = j*dx+dimLPT+dimLSDL+dimLA
        do i = 1, NS
          write(20+i,200), pos,ctal(j)%conc(i,6)
          write(40+i,200), pos,ctal(j)%conc(i,1)
          write(60+i,200), pos,ctal(j)%conc(i,2)
          osmol = osmol + ctal(j)%conc(i,2) 
        end do
        write(80,200) pos,ctal(j)%vol(2)*cw 
        write(81,200) pos,ctal(j)%vol(1)*ctal(j)%conc(2,1)*cw
        write(82,200) pos,ctal(j)%vol(1)*ctal(j)%conc(1,1)*cw
        write(83,200) pos,ctal(j)%vol(1)*ctal(j)%conc(3,1)*cw
        write(84,200) pos,ctal(j)%vol(1)*ctal(j)%conc(11,1)*cw
        write(85,200) pos,osmol
        write(86,200) pos,-log10(ctal(j)%conc(12,1))+3
        write(87,200) pos,-log10(ctal(j)%conc(12,2))+3
        write(88,200) pos,ctal(j)%vol(1) *cw
        write(89,200) pos,ctal(j)%vol(1)*ctal(j)%conc(4,1)*cw
        write(90,200) pos,ctal(j)%vol(1)*ctal(j)%conc(5,1)*cw
        write(91,200) pos,ctal(j)%vol(1)*ctal(j)%conc(7,1)*cw
        write(92,200) pos,ctal(j)%vol(1)*ctal(j)%conc(8,1)*cw
        write(93,200) pos,ctal(j)%vol(1)*ctal(j)%conc(9,1)*cw
        write(94,200) pos,ctal(j)%vol(1)*ctal(j)%conc(10,1)*cw
        write(95,200) pos,ctal(j)%pres

      end do
      do i = 1, NS
        close(unit=20+i) 
        close(unit=40+i) 
        close(unit=60+i)
      end do
      close(74)
      close(75)
      close(76)
      close(77)
      close(78)
      close(79)
      close(80)
      close(81)
      close(82)
      close(83)
      close(84)
      close(85)
      close(86)
      close(87)
      close(88)
      close(90)
 200   format (6f18.11)

      return
      end

!-----------------------------------------------------------------
      subroutine output_dnt_profiles (dct,cnt)

       include 'values.h'
       include 'global.h'
       include 'defs.h'

       ! passed variables
       type (membrane) :: dct(0:NZ), cnt(0:NZ)
       
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	OUTPUT CONCENTRATION PROFILES ALONG DISTAL AND CONNECTING TUBULES
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

      double precision osmol
      double precision cw
      cw = Vref*60.d6 !to convert to nl/min

!     peritubular profiles
      open ( unit=21, file='dntNa_p' )
      open ( unit=22, file='dntK_p' )
      open ( unit=23, file='dntCl_p' )
      open ( unit=24, file='dntHCO3_p' )
      open ( unit=25, file='dntH2CO3_p' )
      open ( unit=26, file='dntCO2_p' )
      open ( unit=27, file='dntHPO4_p' )
      open ( unit=28, file='dntH2PO4_p' )
      open ( unit=29, file='dnturea_p' )
      open ( unit=30, file='dntNH3_p' )
      open ( unit=31, file='dntNH4_p' )
      open ( unit=32, file='dntH_p' )
      open ( unit=33, file='dntHCO2_p' )
      open ( unit=34, file='dntH2CO2_p' )
	  open ( unit=35, file='dntGlu_p' )

!     luminal profiles
      open ( unit=41, file='dntNa_l' )
      open ( unit=42, file='dntK_l' )
      open ( unit=43, file='dntCl_l' )
      open ( unit=44, file='dntHCO3_l' )
      open ( unit=45, file='dntH2CO3_l' )
      open ( unit=46, file='dntCO2_l' )
      open ( unit=47, file='dntHPO4_l' )
      open ( unit=48, file='dntH2PO4_l' )
      open ( unit=49, file='dnturea_l' )
      open ( unit=50, file='dntNH3_l' )
      open ( unit=51, file='dntNH4_l' )
      open ( unit=52, file='dntH_l' )
      open ( unit=53, file='dntHCO2_l' )
      open ( unit=54, file='dntH2CO2_l' )
	  open ( unit=55, file='dntGlu_l' )

!     cellular profiles
      open ( unit=61, file='dntNa_c' )
      open ( unit=62, file='dntK_c' )
      open ( unit=63, file='dntCl_c' )
      open ( unit=64, file='dntHCO3_c' )
      open ( unit=65, file='dntH2CO3_c' )
      open ( unit=66, file='dntCO2_c' )
      open ( unit=67, file='dntHPO4_c' )
      open ( unit=68, file='dntH2PO4_c' )
      open ( unit=69, file='dnturea_c' )
      open ( unit=70, file='dntNH3_c' )
	  open ( unit=71, file='dntNH4_c' )
      open ( unit=72, file='dntH_c' )
      open ( unit=73, file='dntHCO2_c' )
      open ( unit=74, file='dntH2CO2_c' )
	  open ( unit=75, file='dntGlu_c' )

      open ( unit=80, file='dntVol' )
      open ( unit=81, file='dntKflow' )
      open ( unit=82, file='dntNaflow' )
      open ( unit=83, file='dntClflow' )
      open ( unit=84, file='dntNH4flow' )
      open ( unit=85, file='dntos_c')
      open ( unit=86, file='dntpH_l')
      open ( unit=87, file='dntpH_c')
      open ( unit=88, file='dntflow')
      open ( unit=89, file='dntHCO3flow' )
      open ( unit=90, file='dntH2CO3flow' )
      open ( unit=91, file='dntHPO4flow' )
      open ( unit=92, file='dntH2PO4flow' )
      open ( unit=93, file='dntureaflow' )
      open ( unit=94, file='dntNH3flow' )
      open ( unit=95, file='dntpres' )

      dx = dimLD/NZ
      do j = 1,NZ
        pos = dimLPT + dimLSDL + dimLA + dimLT + j*dx
        osmol = dct(j)%conc(NS,2)
        do i = 1, NS
          write(20+i,200), pos,dct(j)%conc(i,6)
          write(40+i,200), pos,dct(j)%conc(i,1)
          write(60+i,200), pos,dct(j)%conc(i,2)
          osmol = osmol + dct(j)%conc(i,2)
        end do
        write(80,200) pos,dct(j)%vol(2)*cw
        write(81,200) pos,dct(j)%vol(1)*dct(j)%conc(2,1)*cw
        write(82,200) pos,dct(j)%vol(1)*dct(j)%conc(1,1)*cw
        write(83,200) pos,dct(j)%vol(1)*dct(j)%conc(3,1)*cw
        write(84,200) pos,dct(j)%vol(1)*dct(j)%conc(11,1)*cw
        write(85,200) pos,osmol
        write(86,200) pos,-log10(dct(j)%conc(12,1))+3
        write(87,200) pos,-log10(dct(j)%conc(12,2))+3
        write(88,200) pos,dct(j)%vol(1)*cw
        write(89,200) pos,dct(j)%vol(1)*dct(j)%conc(4,1)*cw
        write(90,200) pos,dct(j)%vol(1)*dct(j)%conc(5,1)*cw
        write(91,200) pos,dct(j)%vol(1)*dct(j)%conc(7,1)*cw
        write(92,200) pos,dct(j)%vol(1)*dct(j)%conc(8,1)*cw
        write(93,200) pos,dct(j)%vol(1)*dct(j)%conc(9,1)*cw
        write(94,200) pos,dct(j)%vol(1)*dct(j)%conc(10,1)*cw
        write(95,200) pos,dct(j)%pres
      end do
      dx = dimLC/NZ
      do j = 1,NZ
        osmol = cnt(j)%conc(NS,2)
        pos = dimLPT + dimLSDL + dimLA + dimLT + dimLD + j*dx
        do i = 1, NS
          write(20+i,200), pos,cnt(j)%conc(i,6)
          write(40+i,200), pos,cnt(j)%conc(i,1)
          write(60+i,200), pos,cnt(j)%conc(i,2)
          osmol = osmol + cnt(j)%conc(i,2)
        end do
        write(80,200) pos,cnt(j)%vol(2)*cw
        write(81,200) pos,cnt(j)%vol(1)*cnt(j)%conc(2,1)*cw
        write(82,200) pos,cnt(j)%vol(1)*cnt(j)%conc(1,1)*cw
        write(83,200) pos,cnt(j)%vol(1)*cnt(j)%conc(3,1)*cw
        write(84,200) pos,cnt(j)%vol(1)*cnt(j)%conc(11,1)*cw
        write(85,200) pos,osmol
        write(86,200) pos,-log10(cnt(j)%conc(12,1))+3
        write(87,200) pos,-log10(cnt(j)%conc(12,2))+3
        write(88,200) pos,cnt(j)%vol(1)*cw
        write(89,200) pos,cnt(j)%vol(1)*cnt(j)%conc(4,1)*cw
        write(90,200) pos,cnt(j)%vol(1)*cnt(j)%conc(5,1)*cw
        write(91,200) pos,cnt(j)%vol(1)*cnt(j)%conc(7,1)*cw
        write(92,200) pos,cnt(j)%vol(1)*cnt(j)%conc(8,1)*cw
        write(93,200) pos,cnt(j)%vol(1)*cnt(j)%conc(9,1)*cw
        write(94,200) pos,cnt(j)%vol(1)*cnt(j)%conc(10,1)*cw
        write(95,200) pos,cnt(j)%pres

      if(j==1) then
!        print*,"CNT inlet"
!        print*,"Na",cnt(j)%vol(1)*cnt(j)%conc(1,1)*cw
!        print*,"K",cnt(j)%vol(1)*cnt(j)%conc(2,1)*cw
!        print*,"Cl",cnt(j)%vol(1)*cnt(j)%conc(3,1)*cw
!        print*,"HCO3-",cnt(j)%vol(1)*cnt(j)%conc(4,1)*cw
!        print*,"NH4+",cnt(j)%vol(1)*cnt(j)%conc(11,1)*cw
      end if

      end do
      do i = 1, NS
        close(unit=20+i)
        close(unit=40+i)
        close(unit=60+i)
      end do
      close(80)
      close(81)
      close(82)
      close(83)
      close(84)
      close(85)
      close(86)
      close(87)
      close(88)
      close(90)
	  close(91)
	  close(92)
	  close(93)
	  close(94)
	  close(95)
 200   format (6f18.11)

      return
      end


!-----------------------------------------------------------------
      subroutine output_cd_profiles (ccd,omcd,imcd)

       include 'values.h'
       include 'global.h'
       include 'defs.h'

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!       OUTPUT CONCENTRATION PROFILES ALONG DISTAL AND CONNECTING TUBULES
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

       ! passed variables
       type (membrane) :: ccd(0:NZ), omcd(0:NZ), imcd(0:NZ)
       
      double precision osmol
      double precision cw
      cw = Vref*60.d6 !to convert to nl/min

!     peritubular profiles
      open ( unit=21, file='cdNa_p' )
      open ( unit=22, file='cdK_p' )
      open ( unit=23, file='cdCl_p' )
      open ( unit=24, file='cdHCO3_p' )
      open ( unit=25, file='cdH2CO3_p' )
      open ( unit=26, file='cdCO2_p' )
      open ( unit=27, file='cdHPO4_p' )
      open ( unit=28, file='cdH2PO4_p' )
      open ( unit=29, file='cdurea_p' )
      open ( unit=30, file='cdNH3_p' )
      open ( unit=31, file='cdNH4_p' )
      open ( unit=32, file='cdH_p' )
      open ( unit=33, file='cdHCO2_p' )
      open ( unit=34, file='cdH2CO2_p' )
      open ( unit=35, file='cdGlu_p' )

!     luminal profiles
      open ( unit=41, file='cdNa_l' )
      open ( unit=42, file='cdK_l' )
      open ( unit=43, file='cdCl_l' )
      open ( unit=44, file='cdHCO3_l' )
      open ( unit=45, file='cdH2CO3_l' )
      open ( unit=46, file='cdCO2_l' )
      open ( unit=47, file='cdHPO4_l' )
      open ( unit=48, file='cdH2PO4_l' )
      open ( unit=49, file='cdurea_l' )
      open ( unit=50, file='cdNH3_l' )
      open ( unit=51, file='cdNH4_l' )
      open ( unit=52, file='cdH_l' )
      open ( unit=53, file='cdHCO2_l' )
      open ( unit=54, file='cdH2CO2_l' )
      open ( unit=55, file='cdGlu_l' )

!     cellular profiles
      open ( unit=61, file='cdNa_c' )
      open ( unit=62, file='cdK_c' )
      open ( unit=63, file='cdCl_c' )
      open ( unit=64, file='cdHCO3_c' )
      open ( unit=65, file='cdH2CO3_c' )
      open ( unit=66, file='cdCO2_c' )
      open ( unit=67, file='cdHPO4_c' )
      open ( unit=68, file='cdH2PO4_c' )
      open ( unit=69, file='cdurea_c' )
      open ( unit=70, file='cdNH3_c' )
      open ( unit=71, file='cdNH4_c' )
      open ( unit=72, file='cdH_c' )
      open ( unit=73, file='cdHCO2_c' )
      open ( unit=74, file='cdH2CO2_c' )
      open ( unit=55, file='cdGlu_c' )

      open ( unit=80, file='cdVol' )
      open ( unit=81, file='cdKflow' )
      open ( unit=82, file='cdNaflow' )
      open ( unit=83, file='cdClflow' )
      open ( unit=84, file='cdNH4flow' )
      open ( unit=85, file='cdos_c')
      open ( unit=86, file='cdpH_l')
      open ( unit=87, file='cdpH_c')
      open ( unit=88, file='cdflow')
      open ( unit=89, file='cdHCO3flow' )
      open ( unit=90, file='cdH2CO3flow' )
      open ( unit=91, file='cdHPO4flow' )
      open ( unit=92, file='cdH2PO4flow' )
      open ( unit=93, file='cdureaflow' )
      open ( unit=94, file='cdNH3flow' )
      open ( unit=95, file='cdpres' )

      dx = dimLCCD/NZ
      do j = 1,NZ
        pos = dimLPT + dimLSDL + dimLA + dimLT + dimLD + dimLC + j*dx
        osmol = ccd(j)%conc(NS,2)
        do i = 1, NS
          write(20+i,200), pos,ccd(j)%conc(i,6)
          write(40+i,200), pos,ccd(j)%conc(i,1)
          write(60+i,200), pos,ccd(j)%conc(i,2)
          osmol = osmol + ccd(j)%conc(i,2)
        end do
        write(80,200) pos,ccd(j)%vol(2)*cw
        write(81,200) pos,ccd(j)%vol(1)*ccd(j)%conc(2,1)*cw
        write(82,200) pos,ccd(j)%vol(1)*ccd(j)%conc(1,1)*cw
        write(83,200) pos,ccd(j)%vol(1)*ccd(j)%conc(3,1)*cw
        write(84,200) pos,ccd(j)%vol(1)*ccd(j)%conc(11,1)*cw
        write(85,200) pos,osmol
        write(86,200) pos,-log10(ccd(j)%conc(12,1))+3
        write(87,200) pos,-log10(ccd(j)%conc(12,2))+3
        write(88,200) pos,ccd(j)%vol(1)*cw
        write(89,200) pos,ccd(j)%vol(1)*ccd(j)%conc(4,1)*cw
        write(90,200) pos,ccd(j)%vol(1)*ccd(j)%conc(5,1)*cw
        write(91,200) pos,ccd(j)%vol(1)*ccd(j)%conc(7,1)*cw
        write(92,200) pos,ccd(j)%vol(1)*ccd(j)%conc(8,1)*cw
        write(93,200) pos,ccd(j)%vol(1)*ccd(j)%conc(9,1)*cw
        write(94,200) pos,ccd(j)%vol(1)*ccd(j)%conc(10,1)*cw
        write(95,200) pos,ccd(j)%pres

      if(j==1) then
!       print*,"CCD inlet"
!       print*,"Na",ccd(j)%vol(1)*ccd(j)%conc(1,1)*cw
!       print*,"K",ccd(j)%vol(1)*ccd(j)%conc(2,1)*cw
!       print*,"Cl",ccd(j)%vol(1)*ccd(j)%conc(3,1)*cw
!       print*,"HCO3-",ccd(j)%vol(1)*ccd(j)%conc(4,1)*cw
!       print*,"NH4+",ccd(j)%vol(1)*ccd(j)%conc(11,1)*cw
      end if

      end do
      dx = dimLOMC/NZ
      do j = 1,NZ
        osmol = omcd(j)%conc(NS,2)
        pos = dimLPT + dimLSDL + dimLA + dimLT + dimLD + dimLC + dimLCCD + j*dx
        do i = 1, NS
          write(20+i,200), pos,omcd(j)%conc(i,6)
          write(40+i,200), pos,omcd(j)%conc(i,1)
          write(60+i,200), pos,omcd(j)%conc(i,2)
          osmol = osmol + omcd(j)%conc(i,2)
        end do
        write(80,200) pos,omcd(j)%vol(2)*cw
        write(81,200) pos,omcd(j)%vol(1)*omcd(j)%conc(2,1)*cw
        write(82,200) pos,omcd(j)%vol(1)*omcd(j)%conc(1,1)*cw
        write(83,200) pos,omcd(j)%vol(1)*omcd(j)%conc(3,1)*cw
        write(84,200) pos,omcd(j)%vol(1)*omcd(j)%conc(11,1)*cw
        write(85,200) pos,osmol
        write(86,200) pos,-log10(omcd(j)%conc(12,1))+3
        write(87,200) pos,-log10(omcd(j)%conc(12,2))+3
        write(88,200) pos,omcd(j)%vol(1)*cw
        write(89,200) pos,omcd(j)%vol(1)*omcd(j)%conc(4,1)*cw
        write(90,200) pos,omcd(j)%vol(1)*omcd(j)%conc(5,1)*cw
        write(91,200) pos,omcd(j)%vol(1)*omcd(j)%conc(7,1)*cw
        write(92,200) pos,omcd(j)%vol(1)*omcd(j)%conc(8,1)*cw
        write(93,200) pos,omcd(j)%vol(1)*omcd(j)%conc(9,1)*cw
        write(94,200) pos,omcd(j)%vol(1)*omcd(j)%conc(10,1)*cw
        write(95,200) pos,omcd(j)%pres


      if(j==1) then
 !       print*,"OMCD inlet"
 !       print*,"Na",omcd(j)%vol(1)*omcd(j)%conc(1,1)*cw
 !       print*,"K",omcd(j)%vol(1)*omcd(j)%conc(2,1)*cw
 !       print*,"Cl",omcd(j)%vol(1)*omcd(j)%conc(3,1)*cw
 !       print*,"HCO3-",omcd(j)%vol(1)*omcd(j)%conc(4,1)*cw
 !       print*,"NH4+",omcd(j)%vol(1)*omcd(j)%conc(11,1)*cw
      end if

!       write(74,200) j*dx+dimLA,fNKCC2omc(j)
!       write(75,200) j*dx+dimLA,fNHEomc(j)
      end do
      dx = dimLIMC/NZIMC
      do j = 1,NZIMC
        osmol = imcd(j)%conc(NS,2)
        pos = dimLPT + dimLSDL + dimLA + dimLT + dimLD + dimLC + dimLCCD + dimLOMC + j*dx
        do i = 1, NS
          write(20+i,200), pos,imcd(j)%conc(i,6)
          write(40+i,200), pos,imcd(j)%conc(i,1)
          write(60+i,200), pos,imcd(j)%conc(i,2)
          osmol = osmol + imcd(j)%conc(i,2)
        end do
        write(80,200) pos,imcd(j)%vol(2)*cw
        write(81,200) pos,imcd(j)%vol(1)*imcd(j)%conc(2,1)*cw
        write(82,200) pos,imcd(j)%vol(1)*imcd(j)%conc(1,1)*cw
        write(83,200) pos,imcd(j)%vol(1)*imcd(j)%conc(3,1)*cw
        write(84,200) pos,imcd(j)%vol(1)*imcd(j)%conc(11,1)*cw
        write(85,200) pos,osmol
        write(86,200) pos,-log10(imcd(j)%conc(12,1))+3
        write(87,200) pos,-log10(imcd(j)%conc(12,2))+3
        write(88,200) pos,imcd(j)%vol(1)*cw
        write(89,200) pos,imcd(j)%vol(1)*imcd(j)%conc(4,1)*cw
        write(90,200) pos,imcd(j)%vol(1)*imcd(j)%conc(5,1)*cw
        write(91,200) pos,imcd(j)%vol(1)*imcd(j)%conc(7,1)*cw
        write(92,200) pos,imcd(j)%vol(1)*imcd(j)%conc(8,1)*cw
        write(93,200) pos,imcd(j)%vol(1)*imcd(j)%conc(9,1)*cw
        write(94,200) pos,imcd(j)%vol(1)*imcd(j)%conc(10,1)*cw
        write(95,200) pos,imcd(j)%pres


      if(j==1) then
 !       print*,"IMCD inlet"
 !       print*,"Na",imcd(j)%vol(1)*imcd(j)%conc(1,1)*cw
 !       print*,"K",imcd(j)%vol(1)*imcd(j)%conc(2,1)*cw
 !       print*,"Cl",imcd(j)%vol(1)*imcd(j)%conc(3,1)*cw
 !       print*,"HCO3-",imcd(j)%vol(1)*imcd(j)%conc(4,1)*cw
 !       print*,"NH4+",imcd(j)%vol(1)*imcd(j)%conc(11,1)*cw
      end if
      end do
      do i = 1, NS
        close(unit=20+i)
        close(unit=40+i)
        close(unit=60+i)
      end do
      close(74)
      close(75)
      close(76)
      close(77)
      close(78)
      close(79)
      close(80)
      close(81)
      close(82)
      close(83)
      close(84)
      close(85)
      close(86)
      close(87)
      close(88)
      close(90)
 200   format (6f18.11)

      return
      end

!----------------------------------------------------------------------
      subroutine output_sdl_profiles (sdl)
        
        include 'values.h'
        include 'global.h'
        include 'defs.h'

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	OUTPUT CONCENTRATION PROFILES ALONG SDL 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  type (membrane) :: sdl(0:NZ)
      double precision osmol
      double precision cw
      cw = Vref*60.d6 !to convert to nl/min

!     peritubular profiles
      open ( unit=21, file='sdlNa_p' )
      open ( unit=22, file='sdlK_p' )
      open ( unit=23, file='sdlCl_p' )
      open ( unit=24, file='sdlHCO3_p' )
      open ( unit=25, file='sdlH2CO3_p' )
      open ( unit=26, file='sdlCO2_p' )
      open ( unit=27, file='sdlHPO4_p' )
      open ( unit=28, file='sdlH2PO4_p' )
      open ( unit=29, file='sdlurea_p' )
      open ( unit=30, file='sdlNH3_p' )
      open ( unit=31, file='sdlNH4_p' )
      open ( unit=32, file='sdlH_p' )
      open ( unit=33, file='sdlHCO2_p' )

!     luminal profiles
      open ( unit=41, file='sdlNa_l' )
      open ( unit=42, file='sdlK_l' )
      open ( unit=43, file='sdlCl_l' )
      open ( unit=44, file='sdlHCO3_l' )
      open ( unit=45, file='sdlH2CO3_l' )
      open ( unit=46, file='sdlCO2_l' )
      open ( unit=47, file='sdlHPO4_l' )
      open ( unit=48, file='sdlH2PO4_l' )
      open ( unit=49, file='sdlurea_l' )
      open ( unit=50, file='sdlNH3_l' )
      open ( unit=51, file='sdlNH4_l' )
      open ( unit=52, file='sdlH_l' )
      open ( unit=53, file='sdlHCO2_l' )

!     cellular profiles
      open ( unit=61, file='sdlNa_c' )
      open ( unit=62, file='sdlK_c' )
      open ( unit=63, file='sdlCl_c' )
      open ( unit=64, file='sdlHCO3_c' )
      open ( unit=65, file='sdlH2CO3_c' )
      open ( unit=66, file='sdlCO2_c' )
      open ( unit=67, file='sdlHPO4_c' )
      open ( unit=68, file='sdlH2PO4_c' )
      open ( unit=69, file='sdlurea_c' )
      open ( unit=70, file='sdlNH3_c' )
      open ( unit=71, file='sdlNH4_c' )
      open ( unit=72, file='sdlH_c' )
      open ( unit=73, file='sdlHCO2_c' )

      open ( unit=74, file='sdlVol' )
      open ( unit=75, file='sdlKflow' )
      open ( unit=76, file='sdlNaflow' )
      open ( unit=77, file='sdlClflow' )
      open ( unit=78, file='sdlNH4flow' )
      open ( unit=79, file='sdlos_c')
      open ( unit=80, file='sdlpH_l')
      open ( unit=81, file='sdlpH_c')
      open ( unit=82, file='sdlflow')
      open ( unit=83, file='sdlHCO3flow' )
      open ( unit=84, file='sdlH2CO3flow' )
      open ( unit=85, file='sdlHPO4flow' )
      open ( unit=86, file='sdlH2PO4flow' )
      open ( unit=87, file='sdlureaflow' )
      open ( unit=88, file='sdlNH3flow' )
      open ( unit=90, file='sdlpres' )

!     flux values
!     open ( unit=74, file='fluxNKCC2' )
!     open ( unit=75, file='fluxNHE' )

      dx = dimLSDL/NZ
      do j = 1,NZ
        osmol = sdl(j)%conc(NS,2)
        pos = dimLPT + j*dx
        do i = 1, NS-2
          write(20+i,200), pos,sdl(j)%conc(i,6)
          write(40+i,200), pos,sdl(j)%conc(i,1)
          write(60+i,200), pos,sdl(j)%conc(i,2)
          osmol = osmol + sdl(j)%conc(i,2) 
        end do 
        write(74,200) pos,sdl(j)%vol(2) *cw
        write(75,200) pos,sdl(j)%vol(1)*sdl(j)%conc(2,1)*cw
        write(76,200) pos,sdl(j)%vol(1)*sdl(j)%conc(1,1)*cw
        write(77,200) pos,sdl(j)%vol(1)*sdl(j)%conc(3,1)*cw
        write(78,200) pos,sdl(j)%vol(1)*sdl(j)%conc(11,1)*cw
        write(79,200) pos,osmol
 !       write(80,200) pos,-log10(sdl(j)%conc(12,1))+3
 !       write(81,200) pos,-log10(sdl(j)%conc(12,2))+3
        write(82,200) pos,sdl(j)%vol(1) *cw
        write(83,200) pos,sdl(j)%vol(1)*sdl(j)%conc(4,1)*cw
        write(84,200) pos,sdl(j)%vol(1)*sdl(j)%conc(5,1)*cw
        write(85,200) pos,sdl(j)%vol(1)*sdl(j)%conc(7,1)*cw
        write(86,200) pos,sdl(j)%vol(1)*sdl(j)%conc(8,1)*cw
        write(87,200) pos,sdl(j)%vol(1)*sdl(j)%conc(9,1)*cw
        write(88,200) pos,sdl(j)%vol(1)*sdl(j)%conc(10,1)*cw
!       write(74,200) j*dx,fNKCC2sdl(j)
!       write(75,200) j*dx,fNHEsdl(j)
        write(90,200) pos,sdl(j)%pres
      end do
      do i = 1, NS-2
        close(unit=20+i) 
        close(unit=40+i) 
        close(unit=60+i)
      end do
      close(74)
      close(75)
      close(76)
      close(77)
      close(78)
      close(79)
      close(80)
      close(81)
      close(82)
      close(83)
      close(84)
      close(85)
      close(86)
      close(87)
      close(88)
      close(90)
 200   format (6f18.11)

      return
      end

!-----------------------------------------------------------------

!----------------------------------------------------------------------
      subroutine output_pt_profiles (pt)

  include 'values.h'
  include 'global.h'
  include 'defs.h'
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	OUTPUT CONCENTRATION PROFILES ALONG SDL 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  ! passed variables
  type(membrane) :: pt(0:NZ)
  
      double precision osmol
      double precision cw
      cw = Vref*60.d6 !to convert to nl/min

!     peritubular profiles
      open ( unit=21, file='ptNa_p' )
      open ( unit=22, file='ptK_p' )
      open ( unit=23, file='ptCl_p' )
      open ( unit=24, file='ptHCO3_p' )
      open ( unit=25, file='ptH2CO3_p' )
      open ( unit=26, file='ptCO2_p' )
      open ( unit=27, file='ptHPO4_p' )
      open ( unit=28, file='ptH2PO4_p' )
      open ( unit=29, file='pturea_p' )
      open ( unit=30, file='ptNH3_p' )
      open ( unit=31, file='ptNH4_p' )
      open ( unit=32, file='ptH_p' )
      open ( unit=33, file='ptHCO2_p' )
      open ( unit=34, file='ptG_p' )

!     luminal profiles
      open ( unit=41, file='ptNa_l' )
      open ( unit=42, file='ptK_l' )
      open ( unit=43, file='ptCl_l' )
      open ( unit=44, file='ptHCO3_l' )
      open ( unit=45, file='ptH2CO3_l' )
      open ( unit=46, file='ptCO2_l' )
      open ( unit=47, file='ptHPO4_l' )
      open ( unit=48, file='ptH2PO4_l' )
      open ( unit=49, file='pturea_l' )
      open ( unit=50, file='ptNH3_l' )
      open ( unit=51, file='ptNH4_l' )
      open ( unit=52, file='ptH_l' )
      open ( unit=53, file='ptHCO2_l' )
      open ( unit=54, file='ptG_l' )

!     cellular profiles
      open ( unit=61, file='ptNa_c' )
      open ( unit=62, file='ptK_c' )
      open ( unit=63, file='ptCl_c' )
      open ( unit=64, file='ptHCO3_c' )
      open ( unit=65, file='ptH2CO3_c' )
      open ( unit=66, file='ptCO2_c' )
      open ( unit=67, file='ptHPO4_c' )
      open ( unit=68, file='ptH2PO4_c' )
      open ( unit=69, file='pturea_c' )
      open ( unit=70, file='ptNH3_c' )
      open ( unit=71, file='ptNH4_c' )
      open ( unit=72, file='ptH_c' )
      open ( unit=73, file='ptHCO2_c' )

      open ( unit=74, file='ptVol' )
      open ( unit=75, file='ptKflow' )
      open ( unit=76, file='ptNaflow' )
      open ( unit=77, file='ptClflow' )
      open ( unit=78, file='ptNH4flow' )
      open ( unit=79, file='ptos_c')
      open ( unit=80, file='ptpH_l')
      open ( unit=81, file='ptpH_c')
      open ( unit=82, file='ptflow')
      open ( unit=83, file='ptHCO3flow' )
      open ( unit=84, file='ptH2CO3flow' )
      open ( unit=85, file='ptHPO4flow' )
      open ( unit=86, file='ptH2PO4flow' )
      open ( unit=87, file='ptureaflow' )
      open ( unit=88, file='ptNH3flow' )
      open ( unit=89, file='ptGflow' )
      open ( unit=90, file='ptpres' )

!     flux values
!     open ( unit=74, file='fluxNKCC2' )
!     open ( unit=75, file='fluxNHE' )

      dx = dimLPT/NZ
      do j = 1,NZ
        osmol = pt(j)%conc(NS,2)
        do i = 1, NSPT-2
          write(20+i,200), j*dx,pt(j)%conc(i,6)
          write(40+i,200), j*dx,pt(j)%conc(i,1)
          write(60+i,200), j*dx,pt(j)%conc(i,2)
          osmol = osmol + pt(j)%conc(i,2) 
        end do 
        write(54,200), j*dx,pt(j)%conc(15,1)
        write(74,200) j*dx,pt(j)%vol(2) *cw
        write(75,200) j*dx,pt(j)%vol(1)*pt(j)%conc(2,1)*cw
        write(76,200) j*dx,pt(j)%vol(1)*pt(j)%conc(1,1)*cw
        write(77,200) j*dx,pt(j)%vol(1)*pt(j)%conc(3,1)*cw
        write(78,200) j*dx,pt(j)%vol(1)*pt(j)%conc(11,1)*cw
        write(79,200) j*dx,osmol
        write(80,200) j*dx,-log10(pt(j)%conc(12,1))+3
        write(81,200) j*dx,-log10(pt(j)%conc(12,2))+3
        write(82,200) j*dx,pt(j)%vol(1) *cw
        write(83,200) j*dx,pt(j)%vol(1)*pt(j)%conc(4,1)*cw
        write(84,200) j*dx,pt(j)%vol(1)*pt(j)%conc(5,1)*cw
        write(85,200) j*dx,pt(j)%vol(1)*pt(j)%conc(7,1)*cw
        write(86,200) j*dx,pt(j)%vol(1)*pt(j)%conc(8,1)*cw
        write(87,200) j*dx,pt(j)%vol(1)*pt(j)%conc(9,1)*cw
        write(88,200) j*dx,pt(j)%vol(1)*pt(j)%conc(10,1)*cw
        write(89,200) j*dx,pt(j)%vol(1)*pt(j)%conc(15,1)*cw
!       write(74,200) j*dx,fNKCC2sdl(j)
!       write(75,200) j*dx,fNHEsdl(j)
        write(90,200) j*dx,pt(j)%pres
      end do
      do i = 1, NSPT-1
        close(unit=20+i) 
        close(unit=40+i) 
        close(unit=60+i)
      end do
      close(54)
      close(74)
      close(75)
      close(76)
      close(77)
      close(78)
      close(79)
      close(80)
      close(81)
      close(82)
      close(83)
      close(84)
      close(85)
      close(86)
      close(87)
      close(88)
      close(89)
      close(90)
 200   format (6f18.11)

      return
      end

!-----------------------------------------------------------------

!  outputs for Na transport efficiency

      subroutine output_TNa (n,jNaKATPase,jNacell,jNapara,dx,pos0,filetail)

      implicit none

!  passed variables
      integer :: n
      double precision :: dx, pos0
      double precision :: jNaKATPase(0:n),jNacell(0:n),jNapara(0:n)
      !double precision :: jNaKATPase(0:N), jNa(0:N), dx, pos0
      character (len=4) :: filetail

! local variables
      integer j
      character (len=11) :: name1, name2, name3, name4

      name1 = 'jNaactt'  // filetail
      name2 = 'jNacell'  // filetail
      name3 = 'jNapara'  // filetail
      name4 = 'jNaallt'  // filetail

      open (unit=11, file=name1)
      open (unit=12, file=name2)
      open (unit=13, file=name3)
      open (unit=14, file=name4)

      do j = 1, n
         write(11,200) j*dx+pos0, jNaKATPase(j)
         write(12,200) j*dx+pos0, jNacell(j)
         write(13,200) j*dx+pos0, jNapara(j)
         write(14,200) j*dx+pos0, jNacell(j)+jNapara(j)
      end do

      close (unit=11)
      close (unit=12)
      close (unit=13)
      close (unit=14)

 200   format (6f18.11)

      return
      end

