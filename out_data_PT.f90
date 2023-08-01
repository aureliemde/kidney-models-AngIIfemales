!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!
!  File:        out_data_PT.f
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

Subroutine out_data_PT (pt,nephronAct,nephronTNa,nephronQO2)

  include 'values.h'
  include 'global.h'
  include 'defs.h'

  ! passed variables
  type (membrane) :: pt(0:NZ)
  double precision :: nephronAct, nephronTNa, nephronQO2
  
  ! variables used in printing output
  double precision inlet(NSPT),outlet(NSPT),fluxs(NSPT)
  double precision deliv(NSPT),fracdel(NSPT)
  double precision :: totalAct, totalCell, totalPass
  double precision :: nephronActPCT, nephronTNaPCT, nephronQO2PCT
  double precision :: nephronActS3, nephronTNaS3, nephronQO2S3
  ! indices
  integer I,J


!---------------------------------------------------------------------72
!	  Output solute flux 
!---------------------------------------------------------------------72

!	  convert from cm3/s to nl/min/mm tubule
  cwplPT = Vref*60.d0*1.0d6/(dimLPT*10.d0) !to convert to nl/min/mm

!	 If the basis is one nephron
 cw = Vref*60.d0*1d6	!to convert to nl/min

  do I = 1, NSPT
     deliv(I) = pt(0)%conc(I,1)*pt(0)%vol(1)*cw
  end do


  volreab = (pt(0)%vol(1)-pt(NZ)%vol(1))*cw
!  print*
!  print*,"Solute in, reabsorbed, out (pmol/min)"
  print*
  Do I = 1, NSPT
     outlet(I) = pt(NZ)%vol(1)*pt(NZ)%conc(I,1)*cw
     inlet(I) = pt(0)%vol(1)*pt(0)%conc(I,1)*cw
     fluxs(I) = inlet(I) - outlet(I)
     fracdel(I) = fluxs(I)/deliv(I)
!     write(*,'(a,4g12.4)')," Solute",I,inlet(I),-fluxs(I),outlet(I)
  end do
  print*

!---------------------------------------------------------------------72
!	  Results summary
!         reabsorption along PT, and DL delivery
!---------------------------------------------------------------------72

  ReabNapt = (pt(0)%vol(1)*pt(0)%conc(1,1)-pt(NZ)%vol(1)*pt(NZ)%conc(1,1))*cw
  DelivNa = pt(NZ)%vol(1)*pt(NZ)%conc(1,1)*cw

  ReabClpt = (pt(0)%vol(1)*pt(0)%conc(3,1)-pt(NZ)%vol(1)*pt(NZ)%conc(3,1))*cw
  DelivCl = pt(NZ)%vol(1)*pt(NZ)%conc(3,1)*cw

  ReabKpt = (pt(0)%vol(1)*pt(0)%conc(2,1)-pt(NZ)%vol(1)*pt(NZ)%conc(2,1))*cw
  DelivK = pt(NZ)%vol(1)*pt(NZ)%conc(2,1)*cw

  print *,"PT reabs: absolute vs. percentage"
  write(*,'(a,4g12.4)')," Na",ReabNapt, &
			ReabNapt/(pt(0)%vol(1)*pt(0)%conc(1,1)*cw)
  write(*,'(a,4g12.4)')," Cl",ReabClpt, &
			ReabClpt/(pt(0)%vol(1)*pt(0)%conc(3,1)*cw)
  write(*,'(a,4g12.4)')," K ",ReabKpt, &
			ReabKpt/(pt(0)%vol(1)*pt(0)%conc(2,1)*cw)
  write(*,'(a,4g12.4)')," pH",pt(0)%ph(1),pt(NZ)%ph(1)

  ! generate output files
!  call output_pt_profiles (pt)
  
!---------------------------------------------------------------------72
!	  Results summary
!         integrated fluxes along PT
!---------------------------------------------------------------------72

  call compute_o2_consumption (pt,'PT      ',dimLPT,0,NZ,nephronAct,nephronTNa,nephronQO2)

  print*
  NZS3 = 88*NZ/100-1

  nephronActPCT = 0
  nephronTNaPCT = 0
  nephronQO2PCT = 0
  call compute_o2_consumption (pt,'PCT     ',dimLPT*NZS3/NZ,0,NZS3,nephronActPCT,nephronTNaPCT,nephronQO2PCT)

  nephronActS3= 0
  nephronTNaS3 = 0
  nephronQO2S3 = 0
  call compute_o2_consumption (pt,'S3      ',dimLPT-dimLPT*NZS3/NZ,NZS3,NZ,nephronActS3,nephronTNaS3,nephronQO2S3)

  return
END subroutine out_data_PT

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72








