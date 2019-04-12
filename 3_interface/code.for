!************************************************************************
! User element (UEL) for small-deformation phase field-elasticity in 2D 
!  -- plane-strain
!************************************************************************
! Element details:
!************************************************************************
! 
! This subroutine is for a two-dimensional 4-node isoparametric
!  quadrilateral element as shown below with 4pt (full) integration, 
!  as well as plane-strain or axisymmetric settings.
!
! Solution variables (or nodal variables) are the displacements (DOFs 1-2) and DOF 11 for phi (phase field variable).
!
! Material behavior is isotropic linear elasticity.
!
!
! Mechanical body forces and traction- and pressure-type boundary conditions 
!  may be applied to the dummy mesh using the Abaqus built-in commands *Dload 
!  or *Dsload.
!     
!
!              A eta (=xi_2)
!  4-node      |
!   quad       | Face 3
!        4-----------3
!        |     |     |
!        |     |     |
! Face 4 |     ------|---> xi (=xi_1)
!        |           |
!        |           |  Face 2
!        1-----------2
!            Face 1
!
! Kaushik Vijaykumar, November 2014
!
!************************************************************************
! Usage:
!************************************************************************
!
! User element statement in the input file:
!  *User Element,Nodes=4,Type=U1,Iproperties=1,Properties=6,Coordinates=2,Variables=40,Unsymm
!  1,2,11
!
!     State Variables
!     --------------------------------------------------------------
!     Global SDV's (for use in UVARM to make contour plots)
!       1) 11-component of stress                   (sig_11)
!       2) 22-component of stress                   (sig_22)
!       3) 33-component of stress                   (sig_33)
!       4) 23-component of stress                   (sig_23)
!       5) 13-component of stress                   (sig_13)
!       6) 12-component of stress                   (sig_12)
!       7) Trace of strain                          (eps_kk)
!       8) Phase field variable                     (phi)
!       9) Pressure = (-1/3)(sig_11+sig_22+sig_33)  (P)
!      10) Mises Stress                             (sig_mises)
!    
!     In the input file, set '*User output variables' = 10
!    
!     No local state variables are used in this element. However, commented sections
!      demonstrate how to use local state variables for constitutive models that
!      require them (think of plasticity theory). The above parameter 
!      'Variables' should be set to (nInt*nlSdv).
!
!     Material Properties Vector
!     --------------------------------------------------------------
!     Gshear  = props(1)  ! Shear modulus
!     Kbulk   = props(2)  ! Bulk modulus
!     Gc      = props(3)  ! Energy release rate
!     l_0     = props(4)  ! Length scale
!     Kvar    = props(5)  ! Kvar for element stiffness
!     Eta     = prosp(6)  ! Penalty parameter
!      n      = jprops(1) ! positive integer
!
!************************************************************************

      module global

      ! This module is used to transfer SDV's from the UEL
      !  to the UVARM so that SDV's can be visualized on a
      !  dummy mesh
      !
      !  globalSdv(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      !
      !  numElem
      !   Total number of elements in the real mesh, the dummy
      !   mesh needs to have the same number of elements, and 
      !   the dummy mesh needs to have the same number of integ
      !   points. You must set that parameter value here.
      !
      !  offset
      !   Offset in the real mesh numbers and the dummy mesh 
      !   numbers. You must set that parameter value here.
      !
      !  elemPtf
      !   Element pointer
      !

      integer numElem,offset,elemPt,elCount,err
      parameter(numElem=38487,offset=38487)

      real*8, allocatable :: globalSdv(:,:,:)

      
      end module global

!***********************************************************************

      subroutine UVARM(uvar,direct,t,time,dtime,cmname,orname,
     + nuvarm,noel,npt,layer,kspt,kstep,kinc,ndi,nshr,coord,
     + jmac,jmatyp,matlayo,laccfla)

      ! This subroutine is used to transfer SDV's from the UEL
      !  onto the dummy mesh for viewing.
     
      use global
     
      include 'ABA_PARAM.INC'

      character*80 cmname,orname
      character*3 flgray(15)
      dimension uvar(nuvarm),direct(3,3),t(3,3),time(2)
      dimension array(15),jarray(15),jmac(*),jmatyp(*),coord(*)

      ! The dimensions of the variables flgray, array and jarray
      !  must be set equal to or greater than 15.
      
      elCount = noel - offset

      uvar(1) = globalSdv(elCount,npt,1)
      uvar(2) = globalSdv(elCount,npt,2)
      uvar(3) = globalSdv(elCount,npt,3)
      uvar(4) = globalSdv(elCount,npt,4)
      uvar(5) = globalSdv(elCount,npt,5)
      uvar(6) = globalSdv(elCount,npt,6)
      uvar(7) = globalSdv(elCount,npt,7)
      uvar(8) = globalSdv(elCount,npt,8)
      uvar(9) = globalSdv(elCount,npt,9)
      uvar(10) = globalSdv(elCount,npt,10)


      return
      end subroutine uvarm

!****************************************************************************

      subroutine UEL(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     +     props,nprops,coords,mcrd,nnode,uall,duall,vel,accn,jtype,
     +     time,dtime,kstep,kinc,jelem,params,ndload,jdltyp,adlmag,
     +     predef,npredf,lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,
     +     njprop,period)

      use global
      !
      implicit none
      !
      ! variables defined in uel, passed back to Abaqus
      !
      real*8 rhs(mlvarx,*),amatrx(ndofel,ndofel),svars(*),energy(8),
     +  pnewdt
      !
      ! variables passed into UEL
      !
      integer ndofel,nrhs,nsvars,nprops,mcrd,nnode,jtype,kstep,kinc,
     +  jelem,ndload,jdltyp(mdload,*),npredf,lflags(*),mlvarx,mdload,
     +  jprops(*),njprop
      !
      real*8 props(*),coords(mcrd,nnode),uall(ndofel),duall(mlvarx,*),
     +  vel(ndofel),accn(ndofel),time(2),dtime,params(*),
     +  adlmag(mdload,*),predef(2,npredf,nnode),ddlmag(mdload,*),period
      !
      ! variables defined and used in the UEL
      !
      integer i,j,k,jj,nInt,nIntPt,intpt,nDim,ngSdv,nlSdv,stat,pe,nDof,
     +        A11,A12,B11,B12,n,istat,Hs(3),Hs_tr 
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      parameter(nDim=2)  ! number of coordinate dimensions, do not change
      parameter(nDof=3)  ! number of degrees of freedom
      parameter(ngSdv=10) ! number of global state variables per integ pt, do not change
      parameter(nlSdv=4) ! number of local state variables per integ pt, do not change
      parameter(nInt=4)  ! number of integration points, do not change
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      real*8 u(nNode,2),Ru(2*nNode,1),Kuu(2*nNode,2*nNode),
     +  xi(nInt,2),w(nInt),sh0(nNode),sh(nNode),dsh0(nNode,2),
     +  dsh(nNode,2),dshxi(nNode,2),detMapJ0,detMapJ,Rc,urc,Hc_tau(3,3),
     +  epsc_tau(3,3),tr_epsc,R,ur,H_tau(3,3),eps_tau(3,3),tr_eps,
     +  sig_tau(3,3),Ctan(3,3,3,3),enden,Bmat(3,2*nNode),Smat(3,1),
     +  B0mat(3,2*nNode),Amat(3,3),Qmat(3,3),BmatAx(4,2*nNode),
     +  SmatAx(4,1),B0MatAx(4,2*nNode),AmatAx(4,4),QmatAx(4,4),
     +  Rphi(nNode,1),Kphiphi(nNode,nNode),Kuphi(2*nNode,nNode),
     +  phi_grad_tau(nDim,1),phi_tau_0,phi(nNode,1),Kbulk,
     +  Nmat(1,nNode),Pressure,dphi_tau_i,Const1,Const2,
     +  sig_dev_tau(3,3),Iden(3,3),mises,phi_tau_i,kvar,
     +  Bmatphi(2,4),Kphiu(nNode,2*nNode),Gc,l_0,eta,dphi(nNode,1),
     +  Temp1(4,1),Temp2(4,1),Temp3(4,1),sgn,Var1(4,4),Var2(4,4)
     
      real*8 lambda,eigvec(3,3),omega(3),omegamatp(3,3),omegamatn(3,3),
     + ga(3),dab(3,3),Gshear,Smatp(3,1),epsp(3,3),epsn(3,3),
     + sig_taup(3,3),omegamat_trp,omega_tr,endenp,tr_eps0,eps_tau_0(3,3)
     +,sig_taun(3,3),endenn,omegamat(3,3),eps1(3,3),sig_taup1(3,1),Y_t
      
      
      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,
     +     Pi=3.141592653d0,three=3.d0,third=1.d0/3.d0)


      ! Check the procedure type; this should be a 
      !*Coupled temperature-displacement step, which is either 72 or 73
      ! so that we have access to degree of freedom 11
      if((lflags(1).eq.72).or.(lflags(1).eq.73)) then
         !
         ! correct procedure specified
         !
      else
         write(*,*) 'Abaqus does not have the right procedure'
         write(*,*) 'go back and check the procedure type'
         write(*,*) 'lflags(1)=',lflags(1)
         call xit
      endif


      ! Make sure Abaqus knows you are doing a small
      !  deformation problem
      !
      if(lflags(2).eq.1) then
         !
         ! lflags(2)=0 -> small disp.
         ! lflags(2)=1 -> large disp.
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a large displacement analysis'
         write(*,*) 'go in and set nlgeom=no'
         call xit
      endif


      ! Check to see if you are doing a general
      !  step or a linear perturbation step
      !
      if(lflags(4).eq.1) then
         !
         ! lflags(4)=0 -> general step
         ! lflags(4)=1 -> linear perturbation step
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a linear perturbation step'
         call xit         
      endif


      ! Do nothing if a ``dummy'' step
      !
      if(dtime.eq.zero) return


      ! Allocate memory for the globalSdv's used for viewing
      !  results on the dummy mesh
      !
      if(.not.allocated(globalSdv)) then
         !
         ! allocate memory for the globalSdv's
         !
         ! numElem needs to be set in the MODULE
         ! nInt needs to be set in the UEL
         !
         stat=0
         allocate(globalSdv(numElem,nInt,ngSdv),stat=err)
         if(stat.ne.0) then
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) 'error when allocating globalSdv'
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) '   stat=',stat
            write(*,*) '  ngSdv=',ngSdv
            write(*,*) '   nInt=',nInt
            write(*,*) 'numElem=',numElem
            write(*,*) '  nNode=',nNode
            write(*,*) 'lbound(globalSdv)',lbound(globalSdv)
            write(*,*) 'ubound(globalSdv)',ubound(globalSdv)
            write(*,*) '//////////////////////////////////////////////'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) 'error when allocating globalSdv'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) '   stat=',stat
            write(80,*) '  ngSdv=',ngSdv
            write(80,*) '   nInt=',nInt
            write(80,*) 'numElem=',numElem
            write(80,*) '  nNode=',nNode
            write(80,*) 'lbound(globalSdv)=',lbound(globalSdv)
            write(80,*) 'ubound(globalSdv)=',ubound(globalSdv)
            write(80,*) '//////////////////////////////////////////////'
            call xit
         endif
         write(*,*) '-------------------------------------------------'
         write(*,*) '----------- globalSDV ALLOCATED -----------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF ELEMENTS -----------'
         write(*,*) '---------- numElem=',numElem
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF POINTS -------------'
         write(*,*) '---------- nInt=',nInt
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF SDVs ---------------'
         write(*,*) '---------- ngSdv=',ngSdv
         write(*,*) '-------------------------------------------------'
      endif


      ! Counting for the SDV output
      !
      elemPt = jelem


      ! Read input properties 
      Gshear = props(1)
      Kbulk = props(2)
      Gc    = props(3) 
      l_0   = props(4)  
      kvar  = props(5)
      eta   = props(6)
      n     = jprops(1)
      lambda = Kbulk - (two/three)*Gshear !lame parameter lambda
      
      

      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero
      Rphi = zero
      Kuu = zero
      Kphiphi = zero
      Kuphi = zero
      Kphiu = zero
      Energy = zero


      ! Obtain nodal displacements and temperatures.
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
         enddo
         k = k + 1
         phi(i,1)  = Uall(k)
         dphi(i,1) = DUall(k,1)
      enddo
      
      
      !----------------------------------------------------------------
      ! 
      ! Perform calculations at the element centroid in order to 
      !  implement selectively reduced integration on the volumetric
      !  response.
      !
      ! Obtain shape functions and their local gradients at the element
      !  centroid, that means xi_1=xi_2=0.0, and nIntPt=1
      !
      if(nNode.eq.4) then
         call calcShape2DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif


      ! Map shape functions from local to global coordinate system
      !
      if(mcrd.eq.2) then
         call mapShape2Da(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      elseif(mcrd.eq.3) then
         call mapShape2D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
      else
         ! big problem
         write(*,*) 'Unexpected error, mcrd=',mcrd
         call xit
      endif
      
      !  Calculate the displacement gradient at the element centroid 
      !  at the end of the increment.  The subscript tau denotes 
      !  the time at the end of the increment.
      !
      Hc_tau = zero
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Hc_tau(i,j) = Hc_tau(i,j) + dsh0(k,j)*u(k,i)
            enddo
         enddo
      enddo
      eps_tau_0 = half*(Hc_tau + transpose(Hc_tau))
      tr_eps0 = eps_tau_0(1,1) + eps_tau_0(2,2) + eps_tau_0(3,3)
      
      ! Calculating the phase-field variable at the centroidal gauss point
      phi_tau_0 = zero
         
      do k = 1,nNode
         phi_tau_0 = phi_tau_0 + sh0(k)*phi(k,1)
      end do
      !
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      ! Begin the loop over integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.4) then
         !
         ! gauss integration for a rectangular element
         !
         if(nInt.eq.4) then
            call xint2D4pt(xi,w,nIntPt) ! 4-pt integration, nInt=4 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.4'
         call xit
      endif
      


      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt
        

         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.4) then
            call calcShape2DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.4'
            call xit
         endif
        

         ! Map shape functions from local to global coordinate system
         !
         if(mcrd.eq.2) then
            call mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         elseif(mcrd.eq.3) then
            call mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
            if(stat.eq.0) then
               pnewdt = 0.5
               return
            endif
         else
            ! big problem
            write(*,*) 'Unexpected error, mcrd=',mcrd
            call xit 
         endif


         ! Calculate the displacement gradient at this integration point.
         !
         H_tau = zero
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  H_tau(i,j) = H_tau(i,j) + dsh(k,j)*u(k,i)
               enddo
            enddo
         enddo
         
         
         
         ! Calculating the phase-field variable at the integration point
         phi_tau_i = zero
         dphi_tau_i = zero
         
         do k = 1,nNode
            phi_tau_i = phi_tau_i + sh(k)*phi(k,1)
            dphi_tau_i = dphi_tau_i + sh(k)*dphi(k,1)
         end do
         
         ! Obtain the gradient of the phase-field variable at the integration point
         phi_grad_tau = zero

         do i = 1,nDim
            do k = 1,nNode
            phi_grad_tau(i,1) = phi_grad_tau(i,1) + dsh(k,i)*phi(k,1)
            end do
         end do
         

         ! Calculate the strain tensor and volumetric strain at this 
         ! integration point.
         !
         eps_tau = half*(H_tau + transpose(H_tau))
         tr_eps = eps_tau(1,1) + eps_tau(2,2) + eps_tau(3,3)
         
         ! Degrading factor
            Const1 = (one - phi_tau_i)**two + kvar
         ! Calculate the eigenvalues and eigenvectors of Strain eps
            
          call spectral(eps_tau,omega,eigvec,istat)
         
          
          ! Calculate the positive eigenvectors and eigenvalues
           omegamatp = zero
           omegamat_trp = zero
           ga = zero
           dab = zero
           Hs = 0
           omegamat = zero
           
           omega_tr = omega(1) + omega(2) + omega(3)
           
           if(omega_tr.gt.zero) then
              omegamat_trp = omega_tr
              Hs_tr = 1
           else
              omegamat_trp = zero
              Hs_tr = 0
           end if
            
           do i =1,3
           omegamat(i,i) = omega(i)
            if(omega(i).gt.zero) then
               omegamatp(i,i) = omega(i)
               Hs(i) = 1
            else
               omegamatn(i,i) = omega(i)
            end if
           end do  
            
            
          do i =1,3
              ga(i) = const1*(lambda*omegamat_trp*Hs_tr 
     +              + two*Gshear*omegamatp(i,i)*Hs(i)) 
     +              + lambda*omega_tr*(1-Hs_tr)
     +              + two*Gshear*omega(i)*(1-Hs(i))        
          end do

           call onem(Iden)
           do i = 1,3
              do j =1,3
                 dab(i,j) = const1*(lambda*Hs_tr 
     +                    + two*Gshear*Hs(i)*Hs(j)*Iden(i,j))
     +                    + lambda*(1-Hs_tr) 
     +                    + two*Gshear*(1-Hs(i))*(1-Hs(j))*Iden(i,j)
              end do     
           end do
           ! Postive and negative strain tensors
           epsp = matmul(matmul(eigvec,omegamatp),Transpose(eigvec))
           epsn = eps_tau - epsp

         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive update at this integ. point
         !
         call elastic(props,nprops,epsp,epsn,sig_tau,Ctan,enden,omega
     +               ,eigvec,dab,ga,kinc,sig_taup,const1,endenp,eps_tau
     +               ,sig_taun,endenn,Hs)
         !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         
         Pressure = third*(sig_tau(1,1) + sig_tau(2,2) + sig_tau(3,3)) 
         sig_dev_tau = sig_tau - Pressure*Iden
         Mises = ((three/two)*sum(sig_dev_tau*sig_dev_tau))**(half)

         
         ! Compute/update the displacement residual vector.
         !
            !
            ! This is plane strain.
            !
            ! First, compile the B-matrix for the displacment, the matrix of shape function 
            !  derivatives, at this integration point.
            !
            Bmat = zero
            do k=1,nNode
               Bmat(1,1+nDim*(k-1)) = dsh(k,1)
               Bmat(2,2+nDim*(k-1)) = dsh(k,2)
               Bmat(3,1+nDim*(k-1)) = dsh(k,2)
               Bmat(3,2+nDim*(k-1)) = dsh(k,1)
            enddo
           
            !
            ! Next, compile the matrix of stress components,
            !  called the S-matrix.
            !
            Smat(1,1) = sig_tau(1,1)
            Smat(2,1) = sig_tau(2,2)
            Smat(3,1) = sig_tau(1,2)  
            
            if((kinc.le.2)) then
            Y_t = 0.d0
            sig_taup1(1,1) = 0.d0
            sig_taup1(2,2) = 0.d0
            sig_taup1(1,2) = 0.d0
            else
            Y_t = svars(1+jj)
            end if
            
            if(Y_t.gt.endenp) then
            svars(1+jj) = Y_t
            svars(2+jj) = 0.d0
            svars(3+jj) = 0.d0
            svars(4+jj) = 0.d0
            else
            svars(1+jj) = endenp
            svars(2+jj) = sig_taup(1,1)
            svars(3+jj) = sig_taup(2,2)
            svars(4+jj) = sig_taup(1,2)
            endif
            Smatp(1,1) = svars(2+jj)
            Smatp(2,1) = svars(3+jj)
            Smatp(3,1) = svars(4+jj)
            !
            !
            ! Finally, add the contribution of this integration point
            !  to the residual vector.
            !

            Ru = Ru 
     +           - matmul(transpose(Bmat),Smat)*detmapJ*w(intpt)
     
            ! Compute the B-matrix for the phase-field variable phi
            
            Bmatphi = zero
            do k=1,nNode
               Bmatphi(1,k) = dsh(k,1)
               Bmatphi(2,k) = dsh(k,2)
            enddo
            
            ! Shape functions
            do i=1,nNode
               Nmat(1,i) = sh(i)
            end do
            
            if(dphi_tau_i.lt.zero) then
            sgn = -one
            else
            sgn = one
            end if
            
            
            
            ! Use the centroidal value for phi
          Temp1 = Gc*l_0*matmul(Transpose(Bmatphi),phi_grad_tau)
          Temp2 = ((Gc/l_0) + two*svars(1+jj))*Transpose(Nmat)*phi_tau_i
          Temp3 = Transpose(Nmat)*(-two*svars(1+jj) + (eta/(n*dtime))*
     +                                                sgn*dphi_tau_i**n)
            
          
               
            !
            ! Compute the residual vector for the phase-field variable
          
            Rphi = Rphi - (Temp1 + Temp2 + Temp3)*detmapJ*w(intpt)
                     


         ! Compute/update the displacement tangent matrix.
         !
         
            ! This is plane strain
            !
            ! We will need the B-matrix at the element centroid.
            !
            B0mat = zero
            do k=1,nNode
               B0mat(1,1+nDim*(k-1)) = dsh0(k,1)
               B0mat(2,2+nDim*(k-1)) = dsh0(k,2)
               B0mat(3,1+nDim*(k-1)) = dsh0(k,2)
               B0mat(3,2+nDim*(k-1)) = dsh0(k,1)
            enddo
            !
            !  Next, compile the matrix of tangent components -- call 
            !  it the A-matrix. Although the A-matrix is symmetric for linear
            !  elasticity, we won't utilize this property here, so that this
            !  code may be used with more complex constitutive models.
            !
            Amat = zero
            Amat(1,1) = Ctan(1,1,1,1)
            Amat(1,2) = Ctan(1,1,2,2)
            Amat(1,3) = Ctan(1,1,1,2)
            Amat(2,1) = Ctan(2,2,1,1)
            Amat(2,2) = Ctan(2,2,2,2)
            Amat(2,3) = Ctan(2,2,1,2)
            Amat(3,1) = Ctan(1,2,1,1)
            Amat(3,2) = Ctan(1,2,2,2)
            Amat(3,3) = Ctan(1,2,1,2)

            
            ! Finally, add the contribution of this integration point
            ! to the tangent matrix.
            !
            Kuu = Kuu + detmapJ*w(intPt)*
     +            (
     +            matmul(matmul(transpose(Bmat),Amat),Bmat)
     +            )
     
     
         ! Compute the displacement/phi tangent matrix

            Const2 = -two*(one - phi_tau_i)
           
            Kuphi = Kuphi + matmul(transpose(Bmat),matmul(Smatp,
     +                          Nmat))*detmapJ*w(intPt)*Const2 
     
          
          
     
         ! Compute the phi/displacement tangent matrix
         !   
            Kphiu = Transpose(Kuphi)


        Var1 = matmul(transpose(Bmatphi),Bmatphi)*Gc*l_0
        Var2 = ((Gc/l_0) + two*svars(1+jj) + (eta/dtime)*dphi_tau_i*sgn)
     +           *matmul(Transpose(Nmat),Nmat) 
                 
 
         ! Evaluate Stiffness contribution from phase-field equation
           Kphiphi = Kphiphi + detmapJ*w(intPt)*(Var1 + Var2)
           
           ! Save the state variables at this integ point in the
         !  global array used for plotting field output.
         !
         globalSdv(elemPt,intPt,1) = sig_tau(1,1)
         globalSdv(elemPt,intPt,2) = sig_tau(2,2)
         globalSdv(elemPt,intPt,3) = sig_tau(3,3)
         globalSdv(elemPt,intPt,4) = sig_tau(2,3)
         globalSdv(elemPt,intPt,5) = sig_tau(1,3)
         globalSdv(elemPt,intPt,6) = sig_tau(1,2)
         globalSdv(elemPt,intPt,7) = Pressure
         globalSdv(elemPt,intPt,8) = phi_tau_i
         globalSdv(elemPt,intPt,9) = Mises
         globalSdv(elemPt,intPt,10) = enden
         Energy(2) = Energy(2) + enden*detmapJ*w(intpt)

                            
      jj = jj + nlSdv ! setup for the next intPt
      enddo
      
      
      
       
      
      !
      ! End the loop over integration points
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      ! Return to Abaqus the RHS vector and the stiffness matrix.  This
      !  is essentially giving Abaqus the residual and the tangent matrix.
      !
      ! Return to Abaqus the right hand side vector
      !
      do i=1,nNode
           A11 = (nDim+1)*(i-1)+1
           A12 = nDim*(i-1)+1
         !
         ! displacement
         !
         rhs(A11,1) = Ru(A12,1)
         rhs(A11+1,1) = Ru(A12+1,1)
         !
         ! Temperature
         !
         rhs(A11+2,1) = Rphi(i,1)
      enddo
        
      !
      ! Return to Abaqus the tangent matrix and add contribution from displacement and Temperature   
      !
      amatrx = zero
      do i=1,nNode
         do j=1,nNode
            A11 = (nDim+1)*(i-1)+1
            A12 = nDim*(i-1)+1
            B11 = (nDim+1)*(j-1)+1
            B12 = nDim*(j-1)+1
            !
            ! displacement
            !
            amatrx(A11,B11) = Kuu(A12,B12)
            amatrx(A11,B11+1) = Kuu(A12,B12+1)
            amatrx(A11+1,B11) = Kuu(A12+1,B12)
            amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
            !
            ! phase-field 
            !
            amatrx(A11+2,B11+2) = Kphiphi(i,j)
            ! displacement/phase-field
            !
            amatrx(A11,B11+2) = Kuphi(A12,j)
            amatrx(A11+1,B11+2) = Kuphi(A12+1,j)
            ! phase-field/displacement
            !
            amatrx(A11+2,B11) = Kphiu(i,B12)
            amatrx(A11+2,B11+1) = Kphiu(i,B12+1)
         enddo
      enddo
      
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------


      return
      end subroutine uel

!************************************************************************
!     Material subroutines
!************************************************************************

      subroutine elastic(props,nprops,epsp,epsn,sig_tau,Ctan,enden,omega
     +               ,eigvec,dab,ga,kinc,sig_taup,const1,endenp,eps_tau,
     +               sig_taun,endenn,hs)

      implicit none
      !
      integer i,j,k,l,nprops,a,b,c,kinc,counter,jelem
      !
      real*8 props(nprops),sig_tau(3,3),Ctan(3,3,3,3),
     +  enden,Iden(3,3),Gshear,Kbulk,omega(3),const1,eps_tau(3,3),
     +  Eigvec(3,3),lambda,ga(3),dab(3,3),endenp,endenn,epsp(3,3),
     +  epsn(3,3),tr_eps,sig_taup(3,3),sig_taun(3,3),tr_epsp,term1,
     +  hs(3),tol1
      !
      real*8 zero,one,two,three,fourth,third,half,four
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,fourth=1.d0/4.d0,
     +     third=1.d0/3.d0,half=1.d0/2.d0,four = 4.d0,tol1 = 1.d-10)
 

      ! Identity matrix
      !
      call onem(Iden)
 

      ! Obtain material properties
      !
      Gshear = props(1)    ! Shear modulus (MPa)
      Kbulk  = props(2)    ! Bulk modulus (MPa)
      lambda = Kbulk - (two/three)*Gshear !lame parameter lambda (MPa)

      tr_eps = eps_tau(1,1) + eps_tau(2,2) + eps_tau(3,3)
      if(tr_eps.gt.zero) then
      tr_epsp = tr_eps
      else
      tr_epsp = zero
      end if
     
      
      
      ! Compute the stress
      !
      sig_taup = two*Gshear*epsp + lambda*tr_epsp*Iden
      sig_taun = two*Gshear*epsn + lambda*(tr_eps-tr_epsp)*Iden
      sig_tau = const1*sig_taup + sig_taun
   
      

      ! Calculate the tangent modulus
      !
      Ctan = zero
      !
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                
                do a = 1,3
                  do b = 1,3
                  Ctan(i,j,k,l) = Ctan(i,j,k,l)+
     +              dab(a,b)*eigvec(i,a)
     +              *eigvec(j,a)*eigvec(k,b)*eigvec(l,b)
                  
                  end do
                 end do
                 
                 do a = 1,3
                   do c = 1,2
                     if(a.eq.1) then
                      b = c+1
                     end if
                     if(a.eq.2.and.c.eq.1) then
                      b = a-1
                     end if
                     if(a.eq.2.and.c.eq.2) then
                      b = a+1
                     end if 
                     if(a.eq.3) then
                      b = c
                     end if
                     if(abs(omega(a)-omega(b)).le.tol1) then
                        term1 =(const1*two*Gshear*hs(a)
     +                        + two*Gshear*(1-hs(a)))
                     else    
                        term1 =  (ga(a)-ga(b))/(omega(a)-omega(b)) 
                     end if
                    Ctan(i,j,k,l) = Ctan(i,j,k,l)
     +              + term1*
     +              eigvec(i,a)*eigvec(j,b)*
     +              (
     +               eigvec(k,a)*eigvec(l,b) +  
     +               eigvec(k,b)*eigvec(l,a)
     +              )
                    
                    end do
                   end do                       
               enddo
            enddo
         enddo
      enddo
     
      if(kinc.lt.10) then
      Ctan = zero
      
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  Ctan(i,j,k,l) = Ctan(i,j,k,l)
     +                    + Kbulk*Iden(i,j)*Iden(k,l)
     +                    + Gshear*
     +                       (Iden(i,k)*Iden(j,l) + 
     +                        Iden(i,l)*Iden(j,k) - 
     +                        (two/three)*Iden(i,j)*Iden(k,l))
               enddo
            enddo
         enddo
      enddo
     
      end if

      ! Calculate the strain energy density
      !
      endenp = Gshear*sum(epsp*epsp) + half*lambda*(tr_epsp**two) 
      endenn = Gshear*sum(epsn*epsn) 
     +         + half*lambda*((tr_eps-tr_epsp)**two) 
      enden  = endenp*const1 + endenn
      

      return
      end subroutine elastic
      
!****************************************************************************
!     Element subroutines
!****************************************************************************

      subroutine xint2D4pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 4 gauss points for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(4,2), w(4)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 4


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint2D4pt
      
!************************************************************************
      subroutine xint2D1pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 2D elements
      !  using 1 gauss point for integration
      !
      !  xi(nIntPt,2): xi,eta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(1,2), w(1)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w = 4.d0
      

      ! Gauss pt location in master element
      !
      xi(1,1) = -1.d0
      xi(1,2) = -1.d0


      return
      end subroutine xint2D1pt
      
!************************************************************************

      subroutine calcShape2DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer intpt,nDim,nIntPt
      !
      real*8 xi_int(nIntPt,2),sh(4),dshxi(4,2),xi,eta
      !
      real*8 zero,one,fourth
      parameter(zero=0.d0,one=1.d0,fourth=1.d0/4.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      
      ! The shape functions
      !
      sh(1) = fourth*(one - xi)*(one - eta)
      sh(2) = fourth*(one + xi)*(one - eta)
      sh(3) = fourth*(one + xi)*(one + eta)
      sh(4) = fourth*(one - xi)*(one + eta)
      
      
      ! The first derivatives
      !
      dshxi(1,1) = -fourth*(one - eta)
      dshxi(1,2) = -fourth*(one - xi)
      dshxi(2,1) = fourth*(one - eta)
      dshxi(2,2) = -fourth*(one + xi)
      dshxi(3,1) = fourth*(one + eta)
      dshxi(3,2) = fourth*(one + xi)
      dshxi(4,1) = -fourth*(one + eta)
      dshxi(4,2) = fourth*(one - xi)
      

      return
      end subroutine calcShape2DLinear

!************************************************************************

      subroutine mapShape2D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(3,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)
      


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2D

!*************************************************************************

      subroutine mapShape2Da(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      ! This subroutine is exactly the same as the regular mapShape2D
      !  with the exception that coords(2,nNode) here and coords(3,nNode)
      !  in the regular.  I have noticed that a "heat transfer" and 
      !  "static" step uses MCRD=2, but for "coupled-temperature-displacement"
      !  you will get MCRD=3, even for a plane analysis.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,2),dsh(nNode,2),coords(2,nNode),mapJ(2,2),
     +  mapJ_inv(2,2),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)


      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,2
        do j=1,2
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv2D(mapJ,mapJ_inv,detMapJ,stat)

      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))
      

      return
      end subroutine mapShape2Da

!****************************************************************************
!     Utility subroutines
!****************************************************************************

!****************************************************************************

      subroutine matInv2D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse, and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(2,2),A_inv(2,2),det_A,det_A_inv

      
      istat = 1
      
      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv2D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
            
      det_A_inv = 1.d0/det_A
          
      A_inv(1,1) =  det_A_inv*A(2,2)
      A_inv(1,2) = -det_A_inv*A(1,2)
      A_inv(2,1) = -det_A_inv*A(2,1)
      A_inv(2,2) =  det_A_inv*A(1,1)


      return
      end subroutine matInv2D
	
!****************************************************************************

      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)


      do i=1,3
         do J=1,3
	    if (i .eq. j) then
              A(i,j) = 1.0
            else
              A(i,j) = 0.0
            end if
         end do
      end do


      return
      end subroutine onem

!****************************************************************************
!****************************************************************************
!     The following subroutines calculate the spectral
!      decomposition of a symmetric 3 by 3 matrix
!****************************************************************************

      subroutine spectral(A,D,V,istat)
      !
      ! This subroutine calculates the eigenvalues and eigenvectors of
      !  a symmetric 3 by 3 matrix A.
      !
      ! The output consists of a vector D containing the three
      !  eigenvalues in ascending order, and a matrix V whose
      !  columns contain the corresponding eigenvectors.
      !
      implicit none
      !
      integer np,nrot,i,j,istat
      parameter(np=3)
      !
      real*8 D(3),V(3,3),A(3,3),E(3,3)


      E = A
      !
      call jacobi(E,3,np,D,V,nrot,istat)
      !call eigsrt(D,V,3,np)
	

      return
      end subroutine spectral
	
!****************************************************************************

      subroutine jacobi(A,n,np,D,V,nrot,istat)
      !
      ! Computes all eigenvalues and eigenvectors of a real symmetric
      !  matrix A, which is of size n by n, stored in a physical
      !  np by np array.  On output, elements of A above the diagonal
      !  are destroyed, but the diagonal and sub-diagonal are unchanged
      !  and give full information about the original symmetric matrix.
      !  Vector D returns the eigenvalues of A in its first n elements.
      !  V is a matrix with the same logical and physical dimensions as
      !  A whose columns contain, upon output, the normalized
      !  eigenvectors of A.  nrot returns the number of Jacobi rotation
      !  which were required.
      !
      ! This subroutine is taken from 'Numerical Recipes.'
      !
      implicit none
      !
      integer ip,iq,n,nmax,np,nrot,i,j,istat
      parameter (nmax=100)
      !
      real*8 A(np,np),D(np),V(np,np),B(nmax),Z(nmax),
     +  sm,tresh,G,T,H,theta,S,C,tau
     
      
      ! Initialize V to the identity matrix
      !
      call onem(V)
      
      
      ! Initialize B and D to the diagonal of A, and Z to zero.
      !  The vector Z will accumulate terms of the form T*A_PQ as
      !  in equation (11.1.14)
      !
      do ip = 1,n
	B(ip) = A(ip,ip)
	D(ip) = B(ip)
	Z(ip) = 0.d0
      end do
      
      
      ! Begin iteration
      !
      nrot = 0
      do i=1,50
          !
          ! Sum off-diagonal elements
          !
          sm = 0.d0
          do ip=1,n-1
            do iq=ip+1,n
	      sm = sm + dabs(A(ip,iq))
            end do
          end do
          !
          ! If sm = 0., then return.  This is the normal return,
          !  which relies on quadratic convergence to machine
          !  underflow.
          !
          if (sm.eq.0.d0) return
          !
          ! In the first three sweeps carry out the PQ rotation only if
          !  |A_PQ| > tresh, where tresh is some threshold value,
          !  see equation (11.1.25).  Thereafter tresh = 0.
          !
          if (i.lt.4) then
            tresh = 0.2d0*sm/n**2
          else
            tresh = 0.d0
          end if
          !
          do ip=1,n-1
            do iq=ip+1,n
              G = 100.d0*dabs(A(ip,iq))
              !
              ! After four sweeps, skip the rotation if the 
              !  off-diagonal element is small.
              !
	      if ((i.gt.4).and.(dabs(D(ip))+G.eq.dabs(D(ip)))
     +            .and.(dabs(D(iq))+G.eq.dabs(D(iq)))) then
                A(ip,iq) = 0.d0
              else if (dabs(A(ip,iq)).gt.tresh) then
                H = D(iq) - D(ip)
                if (dabs(H)+G.eq.dabs(H)) then
                  !
                  ! T = 1./(2.*theta), equation (11.1.10)
                  !
	          T =A(ip,iq)/H
	        else
	          theta = 0.5d0*H/A(ip,iq)
	          T =1.d0/(dabs(theta)+dsqrt(1.d0+theta**2.d0))
	          if (theta.lt.0.d0) T = -T
	        end if
	        C = 1.d0/dsqrt(1.d0 + T**2.d0)
	        S = T*C
	        tau = S/(1.d0 + C)
	        H = T*A(ip,iq)
	        Z(ip) = Z(ip) - H
	        Z(iq) = Z(iq) + H
	        D(ip) = D(ip) - H
	        D(iq) = D(iq) + H
	        A(ip,iq) = 0.d0
                !
                ! Case of rotations 1 <= J < P
		!		
	        do j=1,ip-1
	          G = A(j,ip)
	          H = A(j,iq)
	          A(j,ip) = G - S*(H + G*tau)
	          A(j,iq) = H + S*(G - H*tau)
	        end do
                !
                ! Case of rotations P < J < Q
                !
	        do j=ip+1,iq-1
	          G = A(ip,j)
	          H = A(j,iq)
	          A(ip,j) = G - S*(H + G*tau)
	          A(j,iq) = H + S*(G - H*tau)
	        end do
                !
                ! Case of rotations Q < J <= N
                !
	        do j=iq+1,n
                  G = A(ip,j)
	          H = A(iq,j)
	          A(ip,j) = G - S*(H + G*tau)
	          A(iq,j) = H + S*(G - H*tau)
	        end do
	        do j = 1,n
	          G = V(j,ip)
	          H = V(j,iq)
	          V(j,ip) = G - S*(H + G*tau)
	          V(j,iq) = H + S*(G - H*tau)
	        end do
	        nrot = nrot + 1
              end if
	    end do
	  end do
          !
          ! Update D with the sum of T*A_PQ, and reinitialize Z
          !
	  do ip=1,n
	    B(ip) = B(ip) + Z(ip)
	    D(ip) = B(ip)
	    Z(ip) = 0.d0
	  end do
	end do


      ! If the algorithm has reached this stage, then there
      !  are too many sweeps.  Print a diagnostic and cut the 
      !  time increment.
      !
      write (*,'(/1X,A/)') '50 iterations in jacobi should never happen'
      istat = 0
      

      return
      end subroutine jacobi
	
!****************************************************************************

      subroutine eigsrt(D,V,n,np)
      !
      ! Given the eigenvalues D and eigenvectors V as output from
      !  jacobi, this subroutine sorts the eigenvales into ascending
      !  order and rearranges the colmns of V accordingly.
      !
      ! The subroutine is taken from 'Numerical Recipes.'
      !
      implicit none
      !
      integer n,np,i,j,k
      !
      real*8 D(np),V(np,np),P
      

      do i=1,n-1
	k = i
	P = D(i)
	do j=i+1,n
	  if (D(j).ge.P) then
	    k = j
	    P = D(j)
	  end if
	end do
	if (k.ne.i) then
	  D(k) = D(i)
	  D(i) = P
	  do j=1,n
	    P = V(j,i)
	    V(j,i) = V(j,k)
	    V(j,k) = P
	  end do
  	end if
      end do
      

      return
      end subroutine eigsrt

!****************************************************************************
