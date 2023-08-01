!---------------------------------------------------------------------72
 
!  *************************** H-K-ATPase ***************************

!---------------------------------------------------------------------72

!    This subroutine computes the H-K-ATPase concentration variables
!	 The units for the kinetic constants are M (moles) - per Alan's email	

subroutine fatpase(n,hkconc,Amat)

  include 'values.h'

  integer n
  double precision kin,kout
  double precision hkconc(4),Amat(n,n)

  Do i=1,n
     Amat(1,i) = 1.0d0
     Do j=2,n
        Amat(j,i)=0.0d0
     end do
  end do
  
!	 Fix reaction rates

  dkf1 = 1.30d7
  dkb1 = 6.50d0
  dkf2a = 8.90d3
  dkb2a = 7.30d4
  dkf2b = 8.90d3
  dkb2b = 7.30d4
  dkf3a = 5.30d9
  dkb3a = 6.60d2
  dkf3b = 5.30d9
  dkb3b = 6.60d2
  dkf4 = 5.0d1
  dkb4 = 2.50d6
  dkf5 = 4.0d1
  dkb5 = 2.0d2
  dkf6a = 5.0d7
  dkb6a = 8.0d12
  dkf6b = 5.0d7
  dkb6b = 8.0d12
  dkf7a = 2.60d10
  dkb7a = 1.50d8
  dkf7b = 2.60d10
  dkb7b = 1.50d8
  dkf8 = 5.40d1
  dkb8 = 3.20d1
  dkf9 = 1.75d0
  dkb9 = 3.50d1
  dkf10 = 5.0d4
  dkb10 = 5.0d1
  dkf11 = 5.0d2
  dkb11 = 5.0d0
  
!	 Fix cellular concentrations of ATP, ADP, Pi

  catp = 2.0d0*1e-3 ! convert from mM to M
  cadp = 0.04d0*1e-3
  cpi = 5.0d0*1e-3

!	 Assign 14 H-K-ATPase variables
!	 species 1 = K2-E1, 2 = K2-E1-ATP, 3 = K-E1-ATP, 4 = E1-ATP, 5 = H-E1-AT!
!	 species 6 = H2-E1-ATP, 7 = H2-E1-P, 8=H2-E2-P, 9 = H-E2-P, 10 = E2-P
!	 species 11 = K-E2-P, 12 = K2-E2-P, 13 = K2-E2, 14 = K2-E2-ATP	
				
  kin = hkconc(1)*1e-3	! convert from mM to M
  kout = hkconc(2)*1e-3
  hin = hkconc(3)*1e-3
  hout = hkconc(4)*1e-3
  
!	  Determine matrix to inverse
	
  Amat(2,2) = +dkf2a
  Amat(2,3) = -(dkb2a*kin+dkf2b)
  Amat(2,4) = +dkb2b*kin
  Amat(3,3) = +dkf2b
  Amat(3,4) = -(dkb2b*kin+dkf3a*hin)
  Amat(3,5) = +dkb3a
  Amat(4,4) = +dkf3a*hin
  Amat(4,5) = -(dkb3a+dkf3b*hin)
  Amat(4,6) = +dkb3b
  Amat(5,5) = +dkf3b*hin
  Amat(5,6) = -(dkb3b+dkf4)
  Amat(5,7) = +dkb4*cadp
  Amat(6,6) = +dkf4
  Amat(6,7) = -(dkb4*cadp+dkf5)
  Amat(6,8) = +dkb5
  Amat(7,7) = +dkf5
  Amat(7,8) = -(dkb5+dkf6a)
  Amat(7,9) = dkb6a*hout
  Amat(8,8) = +dkf6a 
  Amat(8,9) = -(dkb6a*hout+dkf6b) 
  Amat(8,10) = +dkb6b*hout
  Amat(9,9) = +dkf6b 
  Amat(9,10) = -(dkb6b*hout+dkf7a*kout) 
  Amat(9,11) = +dkb7a
  Amat(10,10) = +dkf7a*kout 
  Amat(10,11) = -(dkb7a+dkf7b*kout) 
  Amat(10,12) = +dkb7b
  Amat(11,11) = +dkf7b*kout 
  Amat(11,12) = -(dkb7b+dkf8) 
  Amat(11,13) = +dkb8*cpi
  Amat(12,1) =  +dkb9 
  Amat(12,12) = +dkf8 
  Amat(12,13) = -(dkb8*cpi+dkf9+dkf10*catp)
  Amat(12,14) = +dkb10
  Amat(13,2) = +dkb11 
  Amat(13,13) = +dkf10*catp 
  Amat(13,14) = -(dkb10+dkf11) 
  Amat(14,1) = +dkf1*catp
  Amat(14,2) = -(dkb1+dkf2a+dkb11)
  Amat(14,3) = +dkb2a*kin
  Amat(14,14) = dkf11
  
  return
end subroutine fatpase


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
