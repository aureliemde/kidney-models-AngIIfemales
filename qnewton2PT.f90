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
!   - the S3 segment, where SGLT1 and GLUT1, is accounted for.
!
! 
!---------------------------------------------------------------------72
 
!  *************************** NEWTON SOLVER *************************
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     This is the Newton solver used to determine the luminal and epithelial 
!	concentrations, volumes, and EP below the inlet 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

subroutine qnewton2PT (ptprev,pt,Lz,id,Vol0,CPimpref,CPbuftot1,ND,ncompl2,ntorq2)


  include 'values.h'
  include 'global.h'
  include 'defs.h'
  
  external fcn2PT

  ! passed variables
  type(membrane) :: ptprev, pt
  integer ND
  double precision EPzprev(NC),EPz(NC)
  double precision Vol0(NC),CPimpref,CPbuftot1
  integer Lz,id
  integer :: ncompl2  ! compliant PT?
  integer :: ntorq2   ! torque effects on transport?

  double precision Cb(NSPT,NC),Volb(NC),EPb(NC),phb(NC)
  double precision AVF(ND),A(ND,ND)
  double precision Cprev(NSPT,NC),Volprev(NC),EPprev(NC),phprev(NC)
  double precision osmol(NC)

!      for fjac subroutine
  parameter (num=55)
  parameter (numpar=18+4*NSPT)
  integer iflag,ml,mu
  double precision epsfcn
  double precision x(num),fvec(num),fjac(num,num),wa1(num),wa2(num)
  double precision pars(numpar)

!	for matrix inversion
  integer info, lwork
  parameter (lwork=10000)
  integer ipiv(ND)
  double precision work(lwork)

  TOL = 1.0d-5

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     INITIAL GUESS: ASSUME SAME VALUES AS AT PREVIOUS STEP
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     store initial concentrations
!	CPREV's are used in computing iteration residue
!---------------------------------------------------------------------72

  Do K = 1, NC-1
     
     Do I = 1, NSPT        
        Cb(I,K)=ptprev%conc(I,K)
        Cprev(I,K) = Cb(I,K)
     End Do

     phb(K)=ptprev%ph(K)
     phprev(K)=phb(K)

     Volb(K)=ptprev%vol(K)
     Volprev(K)=Volb(K)

     EPb(K)=ptprev%ep(K)
     EPprev(K)=EPb(K)

  End Do
         
  PMb=ptprev%pres
  PMprev=PMb

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     START ITERATION LOOP TO DETERMINE SOLUTION
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  
  res = 1.0d0
  iter = 0
         
  Do While ( res > TOL ) 
     
     res = 0.0d0   
     iter = iter + 1
!---------------------------------------------------------------------72
!	Create vector x of variables to solve for
!---------------------------------------------------------------------72

     Do I = 1,NSPT
        x(1+3*(I-1))=Cb(I,1)
        x(2+3*(I-1))=Cb(I,2)
        x(3+3*(I-1))=Cb(I,5)
     End do
     x(1+3*NSPT) = Volb(1)
     x(2+3*NSPT) = Volb(2)
     x(3+3*NSPT) = Volb(5)
     x(4+3*NSPT) = EPb(1)
     x(5+3*NSPT) = EPb(2)
     x(6+3*NSPT) = EPb(5)
     x(7+3*NSPT) = PMb
     
!---------------------------------------------------------------------72
!	determine vector fvec of non-linear equations to be solved
!	determine Jacobian matrix fjac corresponding to fvec
!	fjac is a (NUxNU) matrix
!---------------------------------------------------------------------72

     ldfjac = num
     iflag = 1
     ml = num
     mu = num
     epsfcn = 1.d-6
     
     pars(1:NSPT) = ptprev%conc(:,1)
     pars(1+NSPT:2*NSPT) = ptprev%conc(:,2)
     pars(1+2*NSPT:3*NSPT) = ptprev%conc(:,5)
     pars(1+3*NSPT:4*NSPT) = ptprev%conc(:,6)
     pars(1+4*NSPT) = ptprev%vol(1)
     pars(2+4*NSPT) = ptprev%vol(2)
     pars(3+4*NSPT) = ptprev%vol(5)
     pars(4+4*NSPT) = ptprev%pres

     pars(5+4*NSPT) = pt%volEinit
     pars(6+4*NSPT) = pt%volPinit
     if (id == 0) then
        pars(7+4*NSPT) = dimLPT
     endif
     pars(8+4*NSPT) = CPimpref
     pars(9+4*NSPT) = CPbuftot1
     pars(10+4*NSPT) = pt%dkd(1)
     pars(11+4*NSPT) = pt%dkd(2)
     pars(12+4*NSPT) = pt%dkd(5)
     pars(13+4*NSPT) = pt%dkh(1)
     pars(14+4*NSPT) = pt%dkh(2)
     pars(15+4*NSPT) = pt%dkh(5)
     pars(16+4*NSPT) = ncompl2*1.0
     pars(17+4*NSPT) = ntorq2*1.0
     pars(18+4*NSPT) = Vol0(1)

     call fcn2PT(num,x,fvec,iflag,numpar,pars,pt,id)

     call jacobi2(fcn2PT,num,x,fvec,fjac,ldfjac,iflag,epsfcn,numpar,pars,pt,id)

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
        Cb(I,1) = Cprev(I,1)-AVF(1+3*(I-1))
        Cb(I,2) = Cprev(I,2)-AVF(2+3*(I-1))
        Cb(I,5) = Cprev(I,5)-AVF(3+3*(I-1))
        res = dmax1(res,dabs(Cb(I,1)/Cprev(I,1)-1.d0))
        res = dmax1(res,dabs(Cb(I,2)/Cprev(I,2)-1.d0))
        res = dmax1(res,dabs(Cb(I,5)/Cprev(I,5)-1.d0))        
!		print*,"solute",I, res
!		print*,Cb(I,1),Cprev(I,1)
!		print*,Cb(I,2),Cprev(I,2)
!		print*,Cb(I,5),Cprev(I,5)
     End Do
  	 
     phb(1) = -dlog(Cb(12,1)/1.0d3)/dlog(10.0d0)
     phb(2) = -dlog(Cb(12,2)/1.0d3)/dlog(10.0d0)
     phb(5) = -dlog(Cb(12,5)/1.0d3)/dlog(10.0d0)
     
     volb(1) = Volprev(1)-AVF(1+3*NSPT)
     res = dmax1(res,dabs(Volb(1)-Volprev(1)))
     volb(2) = Volprev(2)-AVF(2+3*NSPT)
     res = dmax1(res,dabs(Volb(2)-Volprev(2)))
     volb(5) = Volprev(5)-AVF(3+3*NSPT)
     res = dmax1(res,dabs(Volb(5)-Volprev(5)))

     EPb(1) = EPprev(1)-AVF(4+3*NSPT)
     res = dmax1(res,1.d-3*dabs(EPb(1)-EPprev(1)))
     EPb(2) = EPprev(2)-AVF(5+3*NSPT)
     res = dmax1(res,1.d-3*dabs(EPb(2)-EPprev(2)))
     EPb(5) = EPprev(5)-AVF(6+3*NSPT)
     res = dmax1(res,1.d-3*dabs(EPb(5)-EPprev(5)))

     PMb = PMprev-AVF(7+3*NSPT)
     res = dmax1(res,dabs(PMb-PMprev))

!	update "previous" concentrations, volumes, and potentials

     Do I = 1, NSPT
        Cprev(I,1) = Cb(I,1)
        Cprev(I,2) = Cb(I,2)
        Cprev(I,5) = Cb(I,5)
     End Do

     phprev(1) = phb(1)
     Volprev(1) = Volb(1)
     EPprev(1) = EPb(1)
   
     phprev(2) = phb(2)
     Volprev(2) = Volb(2)
     EPprev(2) = EPb(2)

     phprev(5) = phb(5)
     Volprev(5) = Volb(5)
     EPprev(5) = EPb(5)

     PMprev=PMb

     if (iter >= 11) then
        if (iter < 20) then
        	TOL = 1.0d-3
!        	write(*,'(a,4g12.4)'),"Qnewton2PT - iter > 10",Lz+1,res
        else
			TOL = res*1.010
!			write(*,'(a,4g12.4)'),"Qnewton2PT - iter > 20",Lz+1,res,TOL
        endif
     end if

  End Do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	 Update concentrations at Lz + 1
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  Do I = 1, NSPT
     pt%conc(I,1) = Cb(I,1)
     pt%conc(I,2) = Cb(I,2)
     pt%conc(I,5) = Cb(I,5)
  End Do

  pt%ph(1) = phb(1)
  pt%vol(1) = Volb(1)
  pt%ep(1) = EPb(1)

  pt%ph(2) = phb(2)
  pt%vol(2) = Volb(2)
  pt%ep(2) = EPb(2)

  pt%ph(5) = phb(5)
  pt%vol(5) = Volb(5)
  pt%ep(5) = EPb(5)

  pt%pres = PMb
  

!---------------------------------------------------------------------72
!	 print to output file
!---------------------------------------------------------------------72

  if (id==0) then  
     open ( unit=12, file='PToutlet' )
  end if
  
  Do I = 1,4
     write(12,200),pt%conc(I,1),pt%conc(I,6)
  end Do
  Do I = 5,NSPT
     write(12,210),pt%conc(I,1),pt%conc(I,6)
  end Do
  write(12,200),pt%ph(1),pt%ph(6)
  write(12,200),ptprev%ep(1),ptprev%ep(6)
  write(12,210),pt%vol(1),pt%pres
  
  close ( unit=12 )

!---------------------------------------------------------------------72
!	 print to output file
!---------------------------------------------------------------------72

  if (id==0) then
     open ( unit=14, file='PToutlet_all' )
  end if
  
  Do I = 1,4
     write(14,200),pt%conc(I,2),pt%conc(I,5)
  end Do
  Do I = 5,NSPT
     write(14,210),pt%conc(I,2),pt%conc(I,5)
  end Do
  write(14,200),pt%ph(2),pt%ph(5)
  write(14,200),ptprev%ep(2),ptprev%ep(5)
  write(14,210),pt%vol(2),pt%vol(5)
  
  close ( unit=14 )
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

 200   format (6f12.5)
 210   format (6e12.5)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  return
end subroutine qnewton2PT


