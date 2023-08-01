!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
 
!  *************************** NEWTON SOLVER *************************
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

!  Description: This is the subroutine for the proximal tubule segment
!  This model yields values of concentrations, volumes, and EP 
!	in the lumen and epithelial compartments at equilibrium. 

!---------------------------------------------------------------------72
!									
!  The units are those of the cgs system: cm, mmol, mmol/cm3 = M, s.  
!  The morphological parameters are those of the rat. Parameters are
!  taken from the PT model of A.M. Weinstein (AJP 2007).
!
!---------------------------------------------------------------------72

!   Solute Indices:  1 = Na+, 2 = K+, 3 = Cl-, 4 = HCO3-, 5 = H2CO3, 6 = CO2
!	 7 = HPO4(2-), 8 = H2PO4-, 9 = urea, 10 = NH3, 11 = NH4+, 12 = H+
!	 13 = HCO2-, 14 = H2CO2, 15 = glucose, 16 = Ca2+

!---------------------------------------------------------------------72
!
!  The June 2 version differs from the AMW 2007 model in that:
!   - the pressure equation has been corrected. TS is taken as 1.6 instead of 2.2.
!	- it doesn't include an apical Cl/HCO3 exchnager. See Planelles (2004)
!		and Aronson (1997, 2002). No evidence for its presence in the rat PT.
!   - it includes a small basolateral Cl conductance. See Aronson and Giebisch
!		(1997).
!   - it includes a glucose consumption term, which depends on the Na-K-ATPase
!		activity.
!   - it includes the kinetic description of the Na,glucose cotransporters
!		SGLT1 and SGLT2
!   - it also includes a kinetic description of the GLUT1 & GLUT2 transporters
!   - the S3 segment, where SGLT1 and GLUT1 are located, is taken into account.
!
! 
!---------------------------------------------------------------------72

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     This is the Newton solver used to determine the epithelial 
!	concentrations, volumes, and EP at the inlet (jz = 0)

subroutine qnewton1PT (pt,CPimpref,CPbuftot1,id,ncompl2,ntorq2)

  include 'values.h'
  include 'global.h'
  include 'defs.h'
  
  external fcn1PT

  ! passed variables
  type(membrane) :: pt
  double precision CPimpref,CPbuftot1
  integer id
  integer :: ncompl2  ! compliant PT or not
  integer :: ntorq2   ! torque effects on transport or not
  
  double precision AVF(NUPT),A(NUPT,NUPT)
  double precision Cprev(NSPT,NC),Volprev(NC),EPprev(NC),phprev(NC)
  double precision osmol(NC)


  ! for fjac subroutine
  parameter (num=37)
  parameter (numpar=10) 
  integer iflag,ml,mu
  double precision epsfcn, pars(numpar)
  double precision x(NUPT),fvec(num),fjac(num,num),wa1(num),wa2(num)
  double precision fjacold(num,num),fjacinv(num,num)

  ! for matrix inversion
  integer info, lwork
  parameter (lwork=10000)
  integer ipiv(NUPT)
  double precision work(lwork)

!---------------------------------------------------------------------72
!     private variables
!---------------------------------------------------------------------72

  TOL = 1.0d-5

  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     store initial concentrations
!	CPREV's are used in computing iteration residue
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
       
  Do K = 1, NC
     Do I = 1, NSPT
        Cprev(I,K) = pt%conc(I,K)
     End Do
     phprev(K)=pt%ph(K)
     Volprev(K)=pt%vol(K)
     EPprev(K)=pt%ep(K)
  End Do
   
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     START ITERATION LOOP TO FIND SOLUTION
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  res = 1.0d0
  iter = 0
         
  Do While ( res .gt. TOL ) 
		    
     res = 0.0d0   
     iter = iter + 1
     
!---------------------------------------------------------------------72
!	Create vector x of variables to solve for
!---------------------------------------------------------------------72

     Do I = 1,NSPT
        x(1+2*(I-1))=pt%conc(I,2)
        x(2+2*(I-1))=pt%conc(I,5)
     End do
     x(1+2*NSPT)=pt%vol(2)
     x(2+2*NSPT)=pt%vol(5)
     x(3+2*NSPT)=pt%ep(1)
     x(4+2*NSPT)=pt%ep(2)
     x(5+2*NSPT)=pt%ep(5)
     
!---------------------------------------------------------------------72
!	determine vector fvec of non-linear equations to be solved
!	determine Jacobian matrix fjac corresponding to fvec
!	fjac is a (NUxNU) matrix
!---------------------------------------------------------------------72

     ldfjac = num
     iflag = 1
     ml = num
     mu = num
     epsfcn = 1.d-5

     pars(1) = pt%volEinit
     pars(2) = pt%volPinit
     pars(3) = CPimpref
     pars(4) = CPbuftot1
     pars(5) = pt%dkd(2)
     pars(6) = pt%dkd(5)
     pars(7) = pt%dkh(2)
     pars(8) = pt%dkh(5)
     pars(9) = ncompl2*1.0
     pars(10) = ntorq2*1.0

     call fcn1PT(num,x,fvec,iflag,numpar,pars,pt,id)

     call jacobi2(fcn1PT,num,x,fvec,fjac,ldfjac,iflag,epsfcn,numpar,pars,pt,id)


!---------------------------------------------------------------------72
!     invert Jacobian matrix
!---------------------------------------------------------------------72

     call dgetrf( num, num, fjac, num, ipiv, info )
     
     call dgetri( num, fjac, num, ipiv, work, lwork, info )

!---------------------------------------------------------------------72
!     Calculate AVF: product of fjac  (inverted Jacobian) and fvec
!---------------------------------------------------------------------72

     Do L = 1, num
        AVF(L)=0.d0
        Do M = 1, num
           AVF(L) = AVF(L) + 1.0*fjac(L,M)*fvec(M) ! original = 0.500
        End Do
     End do
     
     
!---------------------------------------------------------------------72
!     Calculate updated concentrations, volumes, and potentials
!	Then update "old" values
!---------------------------------------------------------------------72

     Do I = 1,NSPT
        pt%conc(I,2) = Cprev(I,2)-AVF(1+2*(I-1))
        pt%conc(I,5) = Cprev(I,5)-AVF(2+2*(I-1))
        res = dmax1(res,dabs(pt%conc(I,2)/Cprev(I,2)-1.d0))
        res = dmax1(res,dabs(pt%conc(I,5)/Cprev(I,5)-1.d0))
     End Do
     
     pt%vol(2) = Volprev(2)-AVF(1+2*NSPT)
     res = dmax1(res,dabs(pt%vol(2)-Volprev(2)))
     pt%vol(5) = Volprev(5)-AVF(2+2*NSPT)
     res = dmax1(res,dabs(pt%vol(5)-Volprev(5)))
          
     pt%ep(1) = EPprev(1)-AVF(3+2*NSPT)
     res = dmax1(res,dabs(pt%ep(1)-EPprev(1)))
     pt%ep(2) = EPprev(2)-AVF(4+2*NSPT)
     res = dmax1(res,dabs(pt%ep(2)-EPprev(2)))
     pt%ep(5) = EPprev(5)-AVF(5+2*NSPT)
     res = dmax1(res,dabs(pt%ep(5)-EPprev(5)))
     
!	update "previous" concentrations, volumes, and potentials

     Do K = 2,5
        Do I = 1, NSPT
           Cprev(I,K) = pt%conc(I,K)
        End Do
        pt%ph(K) = -dlog(pt%conc(12,K)/1.d3)/dlog(10.0d0)
        phprev(K) = pt%ph(K)
        Volprev(K) = pt%vol(K)
        EPprev(K) = pt%ep(K)
     End Do
     EPprev(1)=pt%ep(1)
	 
     if (iter >= 11) then
        if (iter < 20) then
           TOL = 1.0d-3
!           write(*,'(a,4g12.4)'),"Qnewton1PT - iter > 10",res
        else
           TOL = res*1.010
!           write(*,'(a,4g12.4)'),"Qnewton1PT - iter > 20",res,TOL
        endif
     end if

  End Do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	                           FINAL RESULTS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

! This variable is needed to determine impermeant concentrations along
! the entire nephron

  VolLumInit = pt%vol(1)

!---------------------------------------------------------------------72
!	 print to output file
!---------------------------------------------------------------------72

  open ( unit=11, file='PTresults' )
  
  Do I = 1,4
     write(11,200),pt%conc(I,1),pt%conc(I,2),pt%conc(I,5)
  end Do
  Do I = 5,NSPT
     write(11,210),pt%conc(I,1),pt%conc(I,2),pt%conc(I,5)
  end Do
  write(11,200),pt%ph(1),pt%ph(2),pt%ph(5)
  write(11,200),pt%vol(1),pt%vol(2),pt%vol(5)
  write(11,200),pt%ep(1)-pt%ep(6),pt%ep(2)-pt%ep(6),pt%ep(5)-pt%ep(6)
  
  close ( unit=11 )
 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

200 format (6f12.5)
210 format (6e12.5)

  return
end subroutine qnewton1PT


