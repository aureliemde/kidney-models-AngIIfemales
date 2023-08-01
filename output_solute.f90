subroutine output_solute (pt,sdl,mtal,ctal,dct,cnt,ccd,omcd,imcd)

  include 'values.h'
  include 'global.h'
  include 'defs.h'
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	OUTPUT CONCENTRATION PROFILES ALONG ENTIRE TAL 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  !passed variables
  type (membrane) :: pt(0:NZ)
  type (membrane) :: sdl(0:NZ)
  type (membrane) :: mtal(0:NZ), ctal(0:NZ)
  type (membrane) :: dct(0:NZ), cnt(0:NZ)
  type (membrane) :: ccd(0:NZ), omcd(0:NZ), imcd(0:NZ)
  
      double precision osmol
      double precision cw
      cw = Vref*60.d6 !to convert to nl/min
!     cw = Vref*60.d0*Nneph !to convert to ml/min


!     peritubular profiles
      open ( unit=21, file='ConcNa_p' )
      open ( unit=22, file='ConcK_p' )
      open ( unit=23, file='ConcCl_p' )
      open ( unit=24, file='ConcS_p' )

!     luminal profiles
      open ( unit=41, file='ConcNa_l' )
      open ( unit=42, file='ConcK_l' )
      open ( unit=43, file='ConcCl_l' )
      open ( unit=44, file='ConcS_l' )


!     cellular profiles
      open ( unit=61, file='ConcNa_c' )
      open ( unit=62, file='ConcK_c' )
      open ( unit=63, file='ConcCl_c' )
      open ( unit=64, file='ConcS_c' )

      open ( unit=80, file='FlowH2O' )
      open ( unit=81, file='FlowNa' )
      open ( unit=82, file='FlowK' )
      open ( unit=83, file='FlowCl' )
      open ( unit=84, file='FlowS' )


      dx = dimLPT/NZ
      do j = 1,NZ
        do i = 1, 3
          write(20+i,200), j*dx,pt(j)%conc(i,6)
          write(40+i,200), j*dx,pt(j)%conc(i,1)
          write(60+i,200), j*dx,pt(j)%conc(i,2)
        end do
        write(24,200), j*dx,pt(j)%conc(7,6)
        write(44,200), j*dx,pt(j)%conc(7,1)
        write(64,200), j*dx,pt(j)%conc(7,2)
        write(80,200), j*dx,pt(j)%vol(1)*cw
        write(81,200), j*dx,pt(j)%vol(1)*pt(j)%conc(1,1)*cw
        write(82,200), j*dx,pt(j)%vol(1)*pt(j)%conc(2,1)*cw
        write(83,200), j*dx,pt(j)%vol(1)*pt(j)%conc(3,1)*cw
        write(84,200), j*dx,pt(j)%vol(1)*pt(j)%conc(7,1)*cw
      end do


      dx = dimLSDL/NZ
      do j = 1,NZ
        pos = dimLPT + j*dx
        do i = 1, 3
          write(20+i,200), pos,sdl(j)%conc(i,6)
          write(40+i,200), pos,sdl(j)%conc(i,1)
          write(60+i,200), pos,sdl(j)%conc(i,2)
        end do
        write(24,200), pos,sdl(j)%conc(7,6)
        write(44,200), pos,sdl(j)%conc(7,1)
        write(64,200), pos,sdl(j)%conc(7,2)
        write(80,200), pos,sdl(j)%vol(1)*cw
        write(81,200), pos,sdl(j)%vol(1)*sdl(j)%conc(1,1)*cw
        write(82,200), pos,sdl(j)%vol(1)*sdl(j)%conc(2,1)*cw
        write(83,200), pos,sdl(j)%vol(1)*sdl(j)%conc(3,1)*cw
        write(84,200), pos,sdl(j)%vol(1)*sdl(j)%conc(7,1)*cw
      end do


      dx = dimLA/NZ
      do j = 0,NZ
        pos = dimLPT + dimLSDL + j*dx
        do i = 1, 3
          write(20+i,200), pos,mtal(j)%conc(i,6)
          write(40+i,200), pos,mtal(j)%conc(i,1)
          write(60+i,200), pos,mtal(j)%conc(i,2)
        end do 
        write(24,200), pos,mtal(j)%conc(7,6)
        write(44,200), pos,mtal(j)%conc(7,1)
        write(64,200), pos,mtal(j)%conc(7,2)
        write(80,200), pos,mtal(j)%vol(1)*cw
        write(81,200), pos,mtal(j)%vol(1)*mtal(j)%conc(1,1)*cw
        write(82,200), pos,mtal(j)%vol(1)*mtal(j)%conc(2,1)*cw
        write(83,200), pos,mtal(j)%vol(1)*mtal(j)%conc(3,1)*cw
        write(84,200), pos,mtal(j)%vol(1)*mtal(j)%conc(7,1)*cw
      end do


      dx = dimLT/NZ
      do j = 1,NZ
        pos = j*dx+dimLPT+dimLSDL+dimLA
        do i = 1, 3
          write(20+i,200), pos,ctal(j)%conc(i,6)
          write(40+i,200), pos,ctal(j)%conc(i,1)
          write(60+i,200), pos,ctal(j)%conc(i,2)
        end do
        write(24,200), pos,ctal(j)%conc(7,6)
        write(44,200), pos,ctal(j)%conc(7,1)
        write(64,200), pos,ctal(j)%conc(7,2)
        write(80,200), pos,ctal(j)%vol(1)*cw
        write(81,200), pos,ctal(j)%vol(1)*ctal(j)%conc(1,1)*cw
        write(82,200), pos,ctal(j)%vol(1)*ctal(j)%conc(2,1)*cw
        write(83,200), pos,ctal(j)%vol(1)*ctal(j)%conc(3,1)*cw
        write(84,200), pos,ctal(j)%vol(1)*ctal(j)%conc(7,1)*cw
      end do


      dx = dimLD/NZ
      do j = 1,NZ
        pos = dimLPT + dimLSDL + dimLA + dimLT + j*dx
        do i = 1, 3
          write(20+i,200), pos,dct(j)%conc(i,6)
          write(40+i,200), pos,dct(j)%conc(i,1)
          write(60+i,200), pos,dct(j)%conc(i,2)
        end do
        write(24,200), pos,dct(j)%conc(7,6)
        write(44,200), pos,dct(j)%conc(7,1)
        write(64,200), pos,dct(j)%conc(7,2)
        write(80,200), pos,dct(j)%vol(1)*cw
        write(81,200), pos,dct(j)%vol(1)*dct(j)%conc(1,1)*cw
        write(82,200), pos,dct(j)%vol(1)*dct(j)%conc(2,1)*cw
        write(83,200), pos,dct(j)%vol(1)*dct(j)%conc(3,1)*cw
        write(84,200), pos,dct(j)%vol(1)*dct(j)%conc(7,1)*cw
      end do


      dx = dimLC/NZ
      do j = 1,NZ
        pos = dimLPT + dimLSDL + dimLA + dimLT + dimLD + j*dx
        do i = 1, 3
          write(20+i,200), pos,cnt(j)%conc(i,6)
          write(40+i,200), pos,cnt(j)%conc(i,1)
          write(60+i,200), pos,cnt(j)%conc(i,2)
        end do
        write(24,200), pos,cnt(j)%conc(7,6)
        write(44,200), pos,cnt(j)%conc(7,1)
        write(64,200), pos,cnt(j)%conc(7,2)
        write(80,200), pos,cnt(j)%vol(1)*cw
        write(81,200), pos,cnt(j)%vol(1)*cnt(j)%conc(1,1)*cw
        write(82,200), pos,cnt(j)%vol(1)*cnt(j)%conc(2,1)*cw
        write(83,200), pos,cnt(j)%vol(1)*cnt(j)%conc(3,1)*cw
        write(84,200), pos,cnt(j)%vol(1)*cnt(j)%conc(7,1)*cw
     end do


      dx = dimLCCD/NZ
      do j = 1,NZ
        pos = dimLPT + dimLSDL + dimLA + dimLT + dimLD + dimLC + j*dx
        do i = 1, 3
          write(20+i,200), pos,ccd(j)%conc(i,6)
          write(40+i,200), pos,ccd(j)%conc(i,1)
          write(60+i,200), pos,ccd(j)%conc(i,2)
        end do
        write(24,200), pos,ccd(j)%conc(7,6)
        write(44,200), pos,ccd(j)%conc(7,1)
        write(64,200), pos,ccd(j)%conc(7,2)
        write(80,200), pos,ccd(j)%vol(1)*cw
        write(81,200), pos,ccd(j)%vol(1)*ccd(j)%conc(1,1)*cw
        write(82,200), pos,ccd(j)%vol(1)*ccd(j)%conc(2,1)*cw
        write(83,200), pos,ccd(j)%vol(1)*ccd(j)%conc(3,1)*cw
        write(84,200), pos,ccd(j)%vol(1)*ccd(j)%conc(7,1)*cw
      end do


      dx = dimLOMC/NZ
      do j = 1,NZ
        pos = dimLPT + dimLSDL + dimLA + dimLT + dimLD + dimLC + dimLCCD + j*dx
        do i = 1, 3
          write(20+i,200), pos,omcd(j)%conc(i,6)
          write(40+i,200), pos,omcd(j)%conc(i,1)
          write(60+i,200), pos,omcd(j)%conc(i,2)
        end do
        write(24,200), pos,omcd(j)%conc(7,6)
        write(44,200), pos,omcd(j)%conc(7,1)
        write(64,200), pos,omcd(j)%conc(7,2)
        write(80,200), pos,omcd(j)%vol(1)*cw
        write(81,200), pos,omcd(j)%vol(1)*omcd(j)%conc(1,1)*cw
        write(82,200), pos,omcd(j)%vol(1)*omcd(j)%conc(2,1)*cw
        write(83,200), pos,omcd(j)%vol(1)*omcd(j)%conc(3,1)*cw
        write(84,200), pos,omcd(j)%vol(1)*omcd(j)%conc(7,1)*cw
      end do


      dx = dimLIMC/NZIMC
      do j = 1,NZIMC
        pos = dimLPT + dimLSDL + dimLA + dimLT + dimLD + dimLC + &
              dimLCCD + dimLOMC + j*dx
        do i = 1, 3
          write(20+i,200), pos,imcd(j)%conc(i,6)
          write(40+i,200), pos,imcd(j)%conc(i,1)
          write(60+i,200), pos,imcd(j)%conc(i,2)
        end do
        write(24,200), pos,imcd(j)%conc(7,6)
        write(44,200), pos,imcd(j)%conc(7,1)
        write(64,200), pos,imcd(j)%conc(7,2)
        write(80,200), pos,imcd(j)%vol(1)*cw
        write(81,200), pos,imcd(j)%vol(1)*imcd(j)%conc(1,1)*cw
        write(82,200), pos,imcd(j)%vol(1)*imcd(j)%conc(2,1)*cw
        write(83,200), pos,imcd(j)%vol(1)*imcd(j)%conc(3,1)*cw
        write(84,200), pos,imcd(j)%vol(1)*imcd(j)%conc(7,1)*cw

     end do

      do i = 1, 3
        close(unit=20+i) 
        close(unit=40+i) 
        close(unit=60+i)
      end do
      close(24)
      close(44)
      close(64)
      close(80)
      close(81)
      close(82)
      close(83)
      close(84)

 200   format (6f18.11)

      return
      end



