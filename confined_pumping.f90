!**********************************************************************
! program confined_pumping:
!     Simulates the dynamics in a hydrodynamically interacting confined active suspension
!     using a continuum model.

! Generates the numerical simulation data in 
! chapter X of my thesis
! Please see (insert link) for details
! of the underlying model and confined fluid solver.

! Last used : August 10, 2015 Barath Ezhilan

!******************************************************************	     
! Global variables 
!******************************************************************

	     module glob

! Parameter definitions	     
	     double precision, parameter :: Pi1 = 1.D0/6.D0, Pes = 0.5D0 !Pi1 is \Lambda in this version of the paper
	     double precision, parameter :: beta = 1.D0, alpha = -40.D0   ! \beta shape parameter and \alpha stresslet strength
	     double precision, parameter :: H = 1.D0, Ly = 5.D0		  ! H channel half-height Ly period in the y direction
	     double precision, parameter :: pi = 4.D0*atan(1.D0), dt = 0.0000075D0 ! dt time step
	     integer, parameter :: ny = 128, nz = 256 ! ny, nz number of points in the y and z directions; FV scheme converges exponentially in a periodic domain
	     integer, parameter :: ntime = 39000000, nsave = 4000 ! nsave : the time period for saving stuff

! (:,ny,nz) arrays
	     
	     double precision, dimension(0:ny+1,0:nz+1) :: c ! concentration field
	     double precision, dimension(2,0:ny+1,0:nz+1) :: m, u ! polarization and velocity field that satisfies the BC
	     double precision, dimension(2,2,0:ny+1,0:nz+1) :: D ! Nematic Order Parameter tensor
	     double precision, dimension(2,2,1:ny,1:nz) :: E, W ! Rate of strain and vorticity tensor
	     double precision, dimension(1:ny,1:nz+1) :: upy, upz ! upy and upz are the y and z components of the periodic Hasimoto solution
	     double precision, dimension(1:ny,1:nz+1) :: ucy, ucz, pre  !ucy and ucz are the y and z components of the velocity correction; pre is the pressure
	     
! 1D (ny) or (nz) arrays

		 double precision, dimension(1:nz+1) :: zzU ! z discretization centered at the velocity grid points.
		 double precision, dimension(0:nz+1) :: zzM ! z discretization centered at the moment grid points.
		 double precision, dimension(1:ny)   :: yy  ! y discretization

		 double precision, dimension(1:ny,1:nz)     :: Fc, Fco !Concentration flux at current and previous time step
	         double precision, dimension(2,1:ny,1:nz)   :: Fm, Fmo !Polarization flux at current and previous time step
	   	 double precision, dimension(2,2,1:ny,1:nz) :: FD, FDo !Nematicity flux at current and previous time step

		 double precision, dimension(2,1:nz)	    :: Favg
		 double precision, dimension(nz)	    :: Dyzavg
		 double precision, dimension(nz+1) 	    :: Umean

		 double precision FTcsum,FBcsum, FTmsum(2),FBmsum(2),FTDsum(2,2), FBDsum(2,2)
		 double precision conFluxSum, mFluxSum(2), DFluxSum(2,2),VelZWallFourierSum
		 double precision Stokescorrect, dF, dvor

		 double precision, dimension(6,6) :: delta 	       ! identity tensor
! 0D

		 double precision dy, dz, time, cmax, cmean, sumdiv, sumdiv2
		 double precision stokes1, stokes2, stokes1sum, stokes2sum, uysum, uzmag
		 integer j, k, a, aa, ab, itime, ifile  
		 character*20 filename1
		 end module glob		 
	     
!**********************************************************************
! 			MAIN PROGRAM
!**********************************************************************
	     
	     program turbulence
	     
	     use Glob

             use ieee_arithmetic
	     
	     implicit none

! Perturbed initial condition

		 call IC
	 
		 call fileoutput
	         
		 call metrics

		 do itime = 1, ntime
		 
		 	time = dble(itime - 1.D0)*dt
! Set ghost cell values 

	     		call moment_ghost
        
! Now we have all the ghost cell values, now calculate velocity field

           		call velfield

! Calculate fluxes
       
         		call flux

!  Time march to find the values at the next time step.

        if (itime.eq.1) then
! Euler for the first time step
           do k = 1, nz
           do j = 1, ny
           
               c(j,k) = c(j,k) + dt*Fc(j,k)
               Fco(j,k) = Fc(j,k)
               
               m(:,j,k) = m(:,j,k) + dt*Fm(:,j,k)
               Fmo(:,j,k) = Fm(:,j,k)

               D(:,:,j,k) = D(:,:,j,k) + dt*FD(:,:,j,k)
               FDo(:,:,j,k) = FD(:,:,j,k)        
               
            enddo
            enddo
        
        else

! adam's Bashforth

           do k = 1, nz
           do j = 1, ny
        
            c(j,k) = c(j,k) + 0.5D0*dt*(3.d0*Fc(j,k)-Fco(j,k))
            Fco(j,k) = Fc(j,k)

            m(:,j,k) = m(:,j,k) + 0.5d0*dt*(3.d0*Fm(:,j,k)-Fmo(:,j,k))
            Fmo(:,j,k) = Fm(:,j,k)

            D(:,:,j,k) = D(:,:,j,k) + 0.5d0*dt*(3.d0*FD(:,:,j,k)-FDo(:,:,j,k))
            FDo(:,:,j,k) = FD(:,:,j,k)
            
            enddo
            enddo    

        endif
        
! File output; c,m,D,u        


		if (mod(itime,nsave) == 0) then
		
			ifile = itime/nsave
		
			call fileoutput

			call metrics
			
			if(ieee_is_nan(cmax)) then
				print*,'Stop - solution blowing up at step',itime
				exit
			endif 	
			
		
		endif
		
		enddo

		end

!***********************************************************************
! Subroutine IC 
! Define initial conditions and variables
!***********************************************************************

	subroutine metrics
	
		Use Glob
		
		implicit none
		
		cmax = abs(c(1,1))
		
	    cmean = 0.D0 
            sumdiv2 = 0.d0
	    uzmag   = 0.d0
	    uysum    = 0.D0
	    stokes1sum = 0.d0
	    stokes2sum = 0.d0
	    
      do k = 1,nz
         do j = 1,ny
         
         if(cmax.lt.abs(c(j,k))) then
            cmax = abs(c(j,k))
         endif
      
		cmean = cmean + c(j,k)

         enddo
      enddo 
	  	
	  	
	 do k = 2,nz
         do j = 2,ny-1
            
	    sumdiv2 = sumdiv2 + abs((u(1,j+1,k) - u(1,j-1,k))/(2.d0*dy) + (u(2,j,k+1)-u(2,j,k-1))/(2.d0*dz))
		stokes1 = -(pre(j+1,k)-pre(j-1,k))/(2.d0*dy) + (ucy(j+1,k) - 2*ucy(j,k) + ucy(j-1,k))/(dy**(2.d0)) + (ucy(j,k+1) - 2*ucy(j,k) + ucy(j,k-1))/(dz**(2.d0)) !+ FA(1,j,k)
		stokes2 = -(pre(j,k+1)-pre(j,k-1))/(2.d0*dz) + (ucz(j+1,k) - 2*ucz(j,k) + ucz(j-1,k))/(dy**(2.d0)) + (ucz(j,k+1) - 2*ucz(j,k) + ucz(j,k-1))/(dz**(2.d0)) !+ FA(2,j,k)

		stokes1sum = stokes1sum + abs(stokes1)
		stokes2sum = stokes2sum + abs(stokes2)	    

		 enddo
	 enddo

	    do k = 1,nz+1
               do j = 1,ny

            	uysum = uysum + u(1,j,k)
				uzmag = uzmag   + abs(u(2,j,k)) 
				
               enddo
            enddo

		sumdiv2 = sumdiv2/dble(ny*nz)
		uzmag   = uzmag/dble(ny*nz)
		stokes1sum = stokes1sum/dble(ny*nz)
		stokes2sum = stokes2sum/dble(ny*nz)
        	uysum      = uysum/dble(ny*(nz+1))
       		uzmag      = uzmag/dble(ny*nz)
	  	cmean = cmean/dble(NY*NZ)

! Most of this stuff is for debugging for the fluid velocity solver

		print *, '==============================================================='	

		print *, 'time =', itime*dt, 'file =',ifile
		print *, 'cmax =', cmax, 'cmean =', cmean

		print *, '=========stokes equation check for the correction part=========' 
		print *, 'average abs divergence = ', sumdiv2
		print *, 'average abs y-momentum of the velocity correction = ', stokes1sum
		print *, 'average abs z-momentum of the velocity correction = ', stokes2sum

		print *, '============velocity characteristics==========================='
		print *, 'average abs z-velocity = ', uzmag		
		print *, 'average y-velocity (net fluid pumping) = ', uysum
		print *, 'VelZWallFourierSum =',VelZWallFourierSum
		print *, 'stokes equation curl sum =', Stokescorrect

		print *, '======Upper and Lower wall Boundary satisfaction==============='
		print *, 'c fluxes',  FTcsum, FBcsum
		print *, 'my fluxes', FTmsum(1), FBmsum(1)
		print *, 'mz fluxes', FTmsum(2), FTmsum(2)
		print *, 'Dyy fluxes', FTDsum(1,1), FBDsum(1,1)
		print *, 'Dyz fluxes', FTDsum(1,2), FTDsum(1,2)
		print *, 'Dzz fluxes', FTDsum(2,2), FTDsum(2,2)

		print *, '======Upper and Lower wall Boundary satisfaction==============='
		print *, 'concentration',  conFluxSum
		print *, 'my',  mFluxSum(1)
		print *, 'mz',  mFluxSum(2)
		print *, 'Dyy', DFluxSum(1,1)
		print *, 'Dyz', DFluxSum(1,2)
		print *, 'Dzz', DFluxSum(2,2)
		print *, '==============================================================='

	return
	end subroutine metrics

!***********************************************************************
! Subroutine IC 
! Define initial conditions and variables
!***********************************************************************

	subroutine IC
	
		Use Glob
		
		implicit none
		
		double precision eps, eps1(16,0:nz+1), eps2(16,0:nz+1),eps3(16,0:nz+1), eps4(16,0:nz+1),Bsq !actually Bsq = B (As in BD_channel_flow.pdf)
		integer nk, idum, ikky   !nk is number of perturbation modes in y direction
		
		nk = 8
		idum = -854881
		eps = 0.015d0

		dy = dble(Ly)/ny
		dz = dble(2.D0*H)/nz

		Bsq = sqrt((6.D0*Pi1+1)/(12.D0*(Pi1**2.D0)*(Pes**2.D0)))
		
		delta(:,:) = 0.D0
        
	        do a = 1,6
        	    delta(a,a) = 1.D0
        	enddo
		
		do k = 1, nz+1
	
		zzU(k) = -1.D0*H + (k-1)*dz 	
		
		enddo
		
		do k = 0, nz+1
	
		zzM(k) = -1.D0*H + (k-1)*dz + dz/2.D0
		
		enddo
		
		do j = 1,ny
		
		yy(j) = (j-1)*dy
		
		enddo
		
! Uniform isotropic	

		c(:,:) = 1.D0
		m(:,:,:) = 0.D0
		D(:,:,:,:) = 0.D0

		do k = 0, nz+1
 	    	do ikky = 1,2*nk

           eps1(ikky,k) = eps*2.D0*(ran(idum)-0.50)
		   eps2(ikky,k) = eps*2.D0*(ran(idum)-0.50)
		   eps3(ikky,k) = eps*2.D0*(ran(idum)-0.50)
		   eps4(ikky,k) = eps*2.D0*(ran(idum)-0.50)
            	enddo
		enddo

		do k = 0, nz+1

! This is the analytical concentration and wall normal polarization from Ezhilan and Saintillan (2015)
			c(:,k)	= (6.D0*Pi1*cosh(Bsq) + cosh(Bsq*zzM(k)))/(6.D0*Pi1*cosh(Bsq) + sinh(Bsq)/Bsq)
			m(2,:,k) = 2.D0*Pi1*Pes*Bsq*sinh(Bsq*zzM(k))/(6.D0*Pi1*cosh(Bsq) + sinh(Bsq)/Bsq)

           do j = 1, ny
		   do ikky = 1,2*nk
		
			c(j,k) = c(j,k) !+ eps1(ikky,k)*cos(2.D0*pi*yy(j)/Ly) + eps2(ikky,k)*sin(2.D0*pi*yy(j)/Ly)
			m(1,j,k) = m(1,j,k) + eps1(ikky,k)*cos(2.D0*pi*yy(j)/Ly) + eps2(ikky,k)*sin(2.D0*pi*yy(j)/Ly)
			D(1,2,j,k) = D(1,2,j,k) + eps3(ikky,k)*cos(2.D0*pi*yy(j)/Ly) + eps4(ikky,k)*sin(2.D0*pi*yy(j)/Ly)

		   enddo
		   enddo

			   D(2,1,j,k) = D(1,2,j,k)

		enddo

		call moment_ghost

		E(:,:,:,:) = 0.D0
		W(:,:,:,:) = 0.D0
		
		upy(:,:) = 0.D0
		upz(:,:) = 0.D0
		
		Fc(:,:) = 0.D0
		Fco(:,:) = 0.D0
		
		Fm(:,:,:) = 0.D0
		Fmo(:,:,:) = 0.D0
		
		FD(:,:,:,:) = 0.D0
		FDo(:,:,:,:) = 0.D0
	
		end subroutine IC

!**********************************************************************
! Subroutine velfield:
!     Obtain velocity field spectrally (Hasimoto solution)
!**********************************************************************

      subroutine velfield

	  use Glob

      implicit none
      
      integer nn(2)

      double precision, dimension(2,1:ny,1:nz)  :: fcompy,fcompz
      double precision, dimension(2,1:ny,1:nz):: FA
      double precision, dimension(2,1:ny,1:nz)  :: ukz, uky
      double precision, dimension(nz)		  :: S13ii	
      double precision kky, kkz, ksq, cst
      double precision fk(2,2), dot(2), dU(2,2)
      double precision S13it
! Copy forces to complex arrays

		VelZWallFourierSum = 0.d0
		stokescorrect      = 0.d0
		Favg		   = 0.d0        
		Dyzavg             = 0.d0

         do k = 1, nz+1 
	   do j = 1, ny

            
      FA(:,j,k) = alpha*(((D(1,:,j+1,k) + D(1,:,j+1,k-1)) - (D(1,:,j-1,k) + D(1,:,j-1,k-1)))/(4.D0*dy) + (D(2,:,j,k) - D(2,:,j,k-1))/(dz))
	  Favg(:,k) = Favg(:,k) + FA(:,j,k)
	  Dyzavg(k) = Dyzavg(k) + D(2,1,j,k)  
            
           enddo
        enddo

	 Favg(:,:) = Favg(:,:)/dble(ny)
         Dyzavg(:) = Dyzavg(:)/dble(ny)


	S13it    = 0.d0 ! \int_{-1}^{1} \Sigma_{13} dz
	S13ii(:) = 0.d0 !\int_{-1}^{z} \Sigma_{13} dz

	do k = 1,nz

	   S13it = S13it + alpha*Dyzavg(k)*dz 
 	   S13ii(k) = S13it !- alpha*D(1,2,1,k)*dz/2.d0 ! To ensure symmetric, the endpoint of the integral have a factor of half

	enddo

	Umean(1) = 0.d0

	do k = 1,nz

	   Umean(k+1) = -1.d0*S13ii(k)+ 0.5d0*(zzU(k)+1.d0)*S13it ! see roman stocker new notes      

	enddo

         do k = 1,nz
         do j = 1,ny

            fcompy(1,j,k) = FA(1,j,k) - Favg(1,k)
            fcompy(2,j,k) = 0.D0
            fcompz(1,j,k) = FA(2,j,k) - Favg(2,k)
            fcompz(2,j,k) = 0.D0

         enddo
      enddo

      nn(1) = ny
      nn(2) = nz

! Take the Fourier transform of force distribution 
      call fourn(fcompy,nn,2,1)
      call fourn(fcompz,nn,2,1)

! Get Fourier coefficients of the velocity field (Hasimoto)
      do k = 1,nz
         kkz = k-1
         if (k.gt.(nz/2)) kkz = k-1-nz
         kkz = -kkz/dble(2.D0*H) ! H is the half-length

         do j = 1,ny
            kky = j-1
            if (j.gt.(ny/2)) kky = j-1-ny
            kky = -kky/dble(Ly)
            
            if ((j.eq.1)) then ! k = 0 is taken care of in the analytical solution
                  
               uky(1,j,k) = 0.d0
               uky(2,j,k) = 0.d0
               ukz(1,j,k) = 0.d0
               ukz(2,j,k) = 0.d0
                
            else
               
               ksq = kky**(2.d0) + kkz**(2.d0)

               fk(1,1) = fcompy(1,j,k)
               fk(1,2) = fcompz(1,j,k)
               fk(2,1) = fcompy(2,j,k)
               fk(2,2) = fcompz(2,j,k)
                  
               dot(1) = fk(1,1)*kky+fk(1,2)*kkz
               dot(2) = fk(2,1)*kky+fk(2,2)*kkz
                  
               uky(1,j,k) = (fk(1,1)-dot(1)*kky/ksq)/ksq/(4*pi**2.d0)
               uky(2,j,k) = (fk(2,1)-dot(2)*kky/ksq)/ksq/(4*pi**2.d0)
               
               ukz(1,j,k) = (fk(1,2)-dot(1)*kkz/ksq)/ksq/(4*pi**2.d0)
               ukz(2,j,k) = (fk(2,2)-dot(2)*kkz/ksq)/ksq/(4*pi**2.d0) 
               
            endif
            
         enddo
      enddo


! Call the inverse Fourier transform
      call fourn(uky,nn,2,-1)
      call fourn(ukz,nn,2,-1)

! Copy arrays
      cst = dble(ny*nz)
      
      do k = 1,nz
         do j = 1,ny

            upy(j,k) = uky(1,j,k)/cst
            upz(j,k) = ukz(1,j,k)/cst 
         
         enddo
      enddo

	 upy(:,nz+1) = upy(:,1)
	 upz(:,nz+1) = upz(:,1)

	  call velcorrection

      do k = 1,nz
         do j = 1,ny
            
            dU(1,1) = ((u(1,j+1,k+1) + u(1,j+1,k)) - (u(1,j-1,k+1) + u(1,j-1,k)))/(4.d0*dy)
     
            dU(2,1) = ((u(2,j+1,k+1) + u(2,j+1,k)) - (u(2,j-1,k+1) + u(2,j-1,k)))/(4.d0*dy)
     
            dU(1,2) = (u(1,j,k+1)-u(1,j,k))/(dz)   

            dU(2,2) = (u(2,j,k+1)-u(2,j,k))/(dz)

            do a = 1,2
               do aa = 1,2
               
! E and W are needed only in the moment grid points      
     
               E(a,aa,j,k) = (dU(a,aa)+dU(aa,a))/2.D0
               W(a,aa,j,k) = (dU(a,aa)-dU(aa,a))/2.D0
      	
               enddo
            enddo

! check the stuff without any flow first
            
         enddo
      enddo


	do k = 5,nz-5
           do j = 5,ny-5
                 
            dF = ((FA(2,j+1,k+1) + FA(2,j+1,k)) - (FA(2,j-1,k+1) + FA(2,j-1,k)))/(4.d0*dy) -  (FA(1,j,k+1)-FA(1,j,k))/(dz)   

!            dF = ( FA(2,j+1,k) - FA(2,j-1,k))/(2.d0*dy) -  (FA(1,j,k+1)-FA(1,j,k-1))/(2*dz)   
     
            dvor = 2.d0*((W(1,2,j+1,k) + W(1,2,j-1,k) - 2.d0*W(1,2,j,k))/dy**(2.d0) + (W(1,2,j,k+1) + W(1,2,j,k-1) - 2.d0*W(1,2,j,k))/dz**(2.d0))
	    
	    Stokescorrect = Stokescorrect + abs(dF-dvor)

	   enddo
	enddo

		VelZWallFourierSum = VelZWallFourierSum/dble(ny*(nz+1))
        Stokescorrect = Stokescorrect/dble((ny-2)*(nz-2))

      return
      end subroutine velfield


!*************************************************************
!Subroutine to calculate the correction in velocity
!*************************************************************
! See the notes and Stokes.nb for more details.

	   SUBROUTINE velcorrection
	   
		Use Glob
		
		implicit none

		 integer nny(1),jmax

		 double precision, dimension(2,1:ny) :: uykT,uykB,uykM,uyksum,uykdiff
		 double precision, dimension(2,1:ny) :: uzkT,uzkB,uzkM,uzksum,uzkdiff
		 double precision, dimension(2,1:ny) :: pres
		 double precision, dimension(2,1:ny) :: aaa,bbb

	     double precision, dimension(2,1:ny,1:nz+1) :: uc1y, uc1z, uc2y, uc2z
	     double precision, dimension(2) :: dwdz,div
	     double precision CCC,SSS,GM1,GM2
	     double precision CCCD,SSSD,GM1D,GM2D
	     double precision cst, kky, kkz, q
		
! First take fourier transform of the y-dependence as a function of z 

        do j = 1,ny
! This has to be modified in all the codes, I had top and bottom reversed before.

            uykT(1,j) = upy(j,NZ+1)
            uykT(2,j) = 0.D0

            uykB(1,j) = upy(j,1)
            uykB(2,j) = 0.D0

            uzkT(1,j) = upz(j,NZ+1)
            uzkT(2,j) = 0.D0
            
            uzkB(1,j) = upz(j,1)
            uzkB(2,j) = 0.D0

	    VelZWallFourierSum = VelZWallFourierSum + abs(upz(j,NZ+1)-upz(j,1))

        enddo

      nny(1) = ny

! Take the Fourier transform of force distribution 
      call fourn(uykT,nny,1,1)
      call fourn(uykB,nny,1,1)
      call fourn(uzkT,nny,1,1)
      call fourn(uzkB,nny,1,1)
      
      cst = dble(ny)
      
      uykT = uykT/cst
      uykB = uykB/cst
      uzkT = uzkT/cst
      uzkB = uzkB/cst
      
      uyksum(:,:)  = (uykT(:,:) + uykB(:,:)) ! Do we actually need /2 here?
      uzksum(:,:)  = (uzkT(:,:) + uzkB(:,:))
      uykdiff(:,:) = (uykT(:,:) - uykB(:,:))
      uzkdiff(:,:) = (uzkT(:,:) - uzkB(:,:))

        
        uc1y(:,:,:) = 0.D0
        uc1z(:,:,:) = 0.D0
        uc2y(:,:,:) = 0.D0
        uc2z(:,:,:) = 0.D0
        
        aaa(:,:)    = 0.D0
	    bbb(:,:)    = 0.D0
	    j = 1
! for kky = 0 is taken care of in the analytical solution
	  
      do j = 2,ny

            kky = j-1
            if (j.gt.(ny/2)) kky = j-1-ny
            kky = -kky/dble(Ly)
            
            q = 2*pi*kky
        
        aaa(1,j) = q*(-uyksum(2,j)*sinh(q*H) + uzkdiff(1,j)*cosh(q*H))/(q*H - sinh(q*H)*cosh(q*H))
     
        aaa(2,j) = q*(uyksum(1,j)*sinh(q*H) + uzkdiff(2,j)*cosh(q*H))/(q*H - sinh(q*H)*cosh(q*H))

        bbb(1,j) = q*(uykdiff(2,j)*cosh(q*H) - uzksum(1,j)*sinh(q*H))/(q*H + sinh(q*H)*cosh(q*H))
     
        bbb(2,j) = -1.D0*q*(uykdiff(1,j)*cosh(q*H) + uzksum(2,j)*sinh(q*H))/(q*H + sinh(q*H)*cosh(q*H))

		enddo

! q = 0 mode, set to zero
     
           do k = 1, nz+1
              do j = 2,ny

            kky = j-1
            if (j.gt.(ny/2)) kky = j-1-ny
            kky = -kky/dble(Ly)
            
            q = 2*pi*kky

            CCC = cosh(q*zzU(k))/cosh(q*H)
            SSS = sinh(q*zzU(k))/sinh(q*H)
        
            GM1 = sinh(q*H)*(zzU(k)*SSS - H*CCC)
            GM2 = cosh(q*H)*(zzU(k)*CCC - H*SSS)
        
        uc1y(:,j,k) = (CCC*uyksum(:,j)/2.D0 + SSS*uykdiff(:,j)/2.D0)
        uc1z(:,j,k) = (CCC*uzksum(:,j)/2.D0 + SSS*uzkdiff(:,j)/2.D0)

        uc2y(1,j,k) = -1.D0/2.D0*(aaa(2,j)*GM1 + bbb(2,j)*GM2)
        uc2y(2,j,k) = 1.D0/2.D0*(aaa(1,j)*GM1  + bbb(1,j)*GM2)

        uc2z(1,j,k) = 1.D0/2.D0*(aaa(1,j)*GM2 + bbb(1,j)*GM1)
        uc2z(2,j,k) = 1.D0/2.D0*(aaa(2,j)*GM2 + bbb(2,j)*GM1)
        
           enddo
        enddo

      do k = 1,nz+1
         do j = 1,ny
         
			kky = j-1
            if (j.gt.(ny/2)) kky = j-1-ny
            kky = -kky/dble(Ly)
            
            q = 2*pi*kky
        
            uykM(1,j) = uc1y(1,j,k) + uc2y(1,j,k) 
            uykM(2,j) = uc1y(2,j,k) + uc2y(2,j,k)

            uzkM(1,j) = uc1z(1,j,k) + uc2z(1,j,k)
            uzkM(2,j) = uc1z(2,j,k) + uc2z(2,j,k)
            
            pres(:,j) = aaa(:,j)*cosh(q*zzU(k)) + bbb(:,j)*sinh(q*zzU(k))

	 enddo
         	pres(:,1) = 0.D0	
            nny(1) = ny

! Take the Fourier transform of velocity
         call fourn(uykM,nny,1,-1)
         call fourn(uzkM,nny,1,-1)
         call fourn(pres,nny,1,-1)
         
         do j = 1,ny
         
            ucy(j,k) = uykM(1,j)
            ucz(j,k) = uzkM(1,j)
            pre(j,k) = pres(1,j)
          
         enddo
      enddo

      do k = 1, nz + 1
         do j = 1,ny
         
            u(1,j,k) = upy(j,k) - ucy(j,k) + Umean(k)
            u(2,j,k) = upz(j,k) - ucz(j,k)
         
         enddo
      enddo
      
! Set velocity boundary conditions (assign values at y-grid points)

		u(:,0,:)  = u(:,NY,:)
    	u(:,NY+1,:) = u(:,1,:)

		return
		end subroutine velcorrection
	   
!*******************************************************************
! Subroutine to set ghost cell points
!*******************************************************************

		subroutine moment_ghost

		Use Glob
		
		implicit none

		double precision  Pi2, check1
		double precision B(6,6), BMP(6,6), BMM(6,6), BMPI(6,6), BMMI(6,6)
        double precision s0(6), s1(6), sN(6), sNp1(6)
		double precision mT, mB, dcdzB, dcdzT, FT, FB, DTT, cT, dmdzT
		integer bmi,bmj
	    
! Boundary matrix definition

        Pi2 = dz/(4.D0*Pi1*Pes)

! B matrix

        B(:,:) = 0.D0

        B(1,3) = 1.D0
        B(2,5) = 1.D0
        B(3,1) = 1.D0/3.D0
        B(3,6) = 1.D0
        B(4,3) = -2.D0/15.D0
        B(5,2) = 1.D0/5.D0
        B(6,3) = 4.D0/15.D0
        
        BMP(:,:) = delta(:,:) + Pi2*B(:,:)
        BMM(:,:) = delta(:,:) - Pi2*B(:,:)       
        
        call inverse(BMP,BMPI,6)
        call inverse(BMM,BMMI,6)
		
! Set ghost cell points 

        do j = 1,ny
        
        s1(1) = c(j,1)
        s1(2) = m(1,j,1)
        s1(3) = m(2,j,1)
        s1(4) = D(1,1,j,1)
        s1(5) = D(1,2,j,1)
        s1(6) = D(2,2,j,1)    
        
        sN(1) = c(j,NZ)
        sN(2) = m(1,j,NZ)
        sN(3) = m(2,j,NZ)
        sN(4) = D(1,1,j,NZ)
        sN(5) = D(1,2,j,NZ)
        sN(6) = D(2,2,j,NZ)    
        
        s0(:) = 0.D0
        sNp1(:) = 0.D0

           do a = 1,6
              do aa = 1,6
                 do ab = 1,6

! Figure out the most optimized way to do this. Shouldn't matter because this is just a 6*6 matrix.
            
                 s0(ab) = s0(ab) + BMPI(ab,aa)*BMM(aa,a)*s1(a)
				 sNp1(ab) = sNp1(ab) + BMMI(ab,aa)*BMP(aa,a)*sN(a)
				 
                 enddo
              enddo
           enddo
           
! z (wall-normal ghost cells)
           
        c(j,0)    = s0(1)
        m(1,j,0)  = s0(2)
        m(2,j,0)  = s0(3)
        D(1,1,j,0) = s0(4)
        D(1,2,j,0) = s0(5)
        D(2,1,j,0) = s0(5)
        D(2,2,j,0) = s0(6)   
        
        c(j,NZ+1)    = sNp1(1)
        m(1,j,NZ+1)  = sNp1(2)
        m(2,j,NZ+1)  = sNp1(3)
        D(1,1,j,NZ+1) = sNp1(4)
        D(1,2,j,NZ+1) = sNp1(5)
        D(2,1,j,NZ+1) = sNp1(5)
        D(2,2,j,NZ+1) = sNp1(6) 

        enddo

! y ghost cells

        do k = 0, nz+1

        c(0,k)    = c(NY,k)
        m(:,0,k)  = m(:,NY,k)
        D(:,:,0,k)  = D(:,:,NY,k)

        c(NY+1,k)    = c(1,k)
        m(:,NY+1,k)  = m(:,1,k)
        D(:,:,NY+1,k)  = D(:,:,1,k)
        
        enddo

		 return
		 end subroutine moment_ghost

!**********************************************************************
! Subroutine flux
!***********************************************************************

		subroutine flux
		
		 Use Glob
		 
		 implicit none
		
		 double precision cB, cT, cR, cL
		 double precision dcdzB, dcdzT, dcdyR, dcdyL
		 double precision mB(2), mT(2), mR(2), mL(2)
		 double precision dmdzB(2), dmdzT(2), dmdyR(2), dmdyL(2)
		 double precision DB(2,2), DTT(2,2), DL(2,2), DR(2,2)
		 double precision dDdzB(2,2), dDdzT(2,2), dDdyR(2,2), dDdyL(2,2)
		 
		 double precision uB(2), uT(2), uR(2), uL(2)
		 
		 double precision FR, FL, FT, FB
		 
		 integer b

		FTcsum	      = 0.d0
		FTmsum(:)     = 0.d0
		FTDsum(:,:)   = 0.d0	 
   
		conFluxSum    = 0.d0
		mFluxSum(:)   = 0.d0
		DFluxSum(:,:) = 0.d0

		FBcsum        = 0.d0
		FBmsum(:)     = 0.d0
		FBDsum(:,:)   = 0.d0


        do j = 1, NY
            do k = 1, NZ

        uB(:) = u(:,j,k)
        uT(:) = u(:,j,k+1)
        uL(:) = (u(:,j,k+1) + u(:,j-1,k+1) + u(:,j,k) + u(:,j-1,k))/4.D0
        uR(:) = (u(:,j,k+1) + u(:,j+1,k+1) + u(:,j,k) + u(:,j+1,k))/4.D0
     
        cB = (c(j,k) + c(j,k-1))/2.D0
        cT = (c(j,k) + c(j,k+1))/2.D0
        cL = (c(j,k) + c(j-1,k))/2.D0
        cR = (c(j,k) + c(j+1,k))/2.D0
        
        dcdzB = (c(j,k) - c(j,k-1))/dz
        dcdzT = (c(j,k+1) - c(j,k))/dz
        dcdyL = (c(j,k) - c(j-1,k))/dy
        dcdyR = (c(j+1,k) - c(j,k))/dy

        mB(:) = (m(:,j,k) + m(:,j,k-1))/2.D0
        mT(:) = (m(:,j,k) + m(:,j,k+1))/2.D0
        mL(:) = (m(:,j,k) + m(:,j-1,k))/2.D0
        mR(:) = (m(:,j,k) + m(:,j+1,k))/2.D0

        dmdzB(:) = (m(:,j,k) - m(:,j,k-1))/dz
        dmdzT(:) = (m(:,j,k+1) - m(:,j,k))/dz
        dmdyL(:) = (m(:,j,k) - m(:,j-1,k))/dy
        dmdyR(:) = (m(:,j+1,k) - m(:,j,k))/dy

        DB(:,:) = (D(:,:,j,k) + D(:,:,j,k-1))/2.D0
        DTT(:,:) = (D(:,:,j,k) + D(:,:,j,k+1))/2.D0
        DL(:,:) = (D(:,:,j,k) + D(:,:,j-1,k))/2.D0
        DR(:,:) = (D(:,:,j,k) + D(:,:,j+1,k))/2.D0
        
        dDdzB(:,:) = (D(:,:,j,k) - D(:,:,j,k-1))/dz
        dDdzT(:,:) = (D(:,:,j,k+1) - D(:,:,j,k))/dz
        dDdyL(:,:) = (D(:,:,j,k) - D(:,:,j-1,k))/dy
        dDdyR(:,:) = (D(:,:,j+1,k) - D(:,:,j,k))/dy

! concentration equation flux terms

        FR = uR(1)*cR + 2.D0*Pes*mR(1) - 4*pi1*(Pes)**(2.D0)*dcdyR
        FL = uL(1)*cL + 2.D0*Pes*mL(1) - 4*pi1*(Pes)**(2.D0)*dcdyL
        
        FT = uT(2)*cT + 2.D0*Pes*mT(2) - 4*pi1*(Pes)**(2.D0)*dcdzT
        FB = uB(2)*cB + 2.D0*Pes*mB(2) - 4*pi1*(Pes)**(2.D0)*dcdzB

!================================== Set boundary fluxes to zero and calculate them to verify the ghost point calculations==========================

	if (k.eq.1) then
	FBcsum = FBcsum + abs(2.D0*Pes*mB(2) - 4*pi1*(Pes)**(2.D0)*dcdzB)
	FB = 0.D0
	endif

	if (k.eq.NZ) then
	FTcsum = FTcsum + abs(2.D0*Pes*mT(2) - 4*pi1*(Pes)**(2.D0)*dcdzT)
	FT = 0.D0
	endif

!====================================================================================================================================================

        Fc(j,k) = -((FR - FL)/dy + (FT - FB)/dz)


	conFluxSum = conFluxSum - ((FR - FL)/dy + (FT - FB)/dz)

! m equation flux terms

        do a = 1,2

        FR = uR(1)*mR(a) + 2.D0*Pes*(DR(a,1) + cR*delta(a,1)/3.D0) - 4*Pi1*(Pes)**(2.D0)*dmdyR(a)

        FL = uL(1)*mL(a) + 2.D0*Pes*(DL(a,1) + cL*delta(a,1)/3.D0) - 4*Pi1*(Pes)**(2.D0)*dmdyL(a)
        
        FT = uT(2)*mT(a) + 2.D0*Pes*(DTT(a,2) + cT*delta(a,2)/3.D0) - 4*Pi1*(Pes)**(2.D0)*dmdzT(a)   
 
        FB = uB(2)*mB(a) + 2.D0*Pes*(DB(a,2) + cB*delta(a,2)/3.D0) - 4*Pi1*(Pes)**(2.D0)*dmdzB(a)   

!================================== Set boundary fluxes to zero and calculate them to verify the ghost point calculations==========================

	if (k.eq.1) then
	FBmsum(a) = FBmsum(a) + abs(FB - uB(2)*mB(a))
	FB = 0.D0
	endif

	if (k.eq.NZ) then
	FTmsum(a) = FTmsum(a) + abs(FT - uT(2)*mT(a))
	FT = 0.D0
	endif
!====================================================================================================================================================

	mFluxSum(a) = mFluxSum(a) - ((FR - FL)/dy + (FT - FB)/dz)

        Fm(a,j,k) = -((FR - FL)/dy + (FT - FB)/dz) - 2.D0*m(a,j,k)
        


            do aa = 1,2
        
                Fm(a,j,k) = Fm(a,j,k) + (3.D0/5.D0*beta*E(a,aa,j,k) + W(a,aa,j,k))*m(aa,j,k) 
          
            enddo
        enddo

! D equation flux terms
        
        do a = 1,2 !a is i 
           do b = 1,2 !b is j
        
        FR = uR(1)*DR(a,b) + 2.D0*Pes*(1.D0/5.D0*(mR(a)*delta(b,1)  + mR(b)*delta(a,1)) - 2.D0/15.D0*mR(1)*delta(a,b)) 
        FR = FR - 4.D0*Pi1*(Pes)**(2.D0)*dDdyR(a,b) 

        FL = uL(1)*DL(a,b) + 2.D0*Pes*(1.D0/5.D0*(mL(a)*delta(b,1)  + mL(b)*delta(a,1)) - 2.D0/15.D0*mL(1)*delta(a,b)) 
        FL = FL - 4.D0*pi1*(Pes)**(2.D0)*dDdyL(a,b) 
        
        FT = uT(2)*DTT(a,b) + 2.D0*Pes*(1.D0/5.D0*(mT(a)*delta(b,2)  + mT(b)*delta(a,2)) - 2.D0/15.D0*mT(2)*delta(a,b)) 
        FT = FT - 4.D0*pi1*(Pes)**(2.D0)*dDdzT(a,b) 
 
        FB = uB(2)*DB(a,b) + 2.D0*Pes*(1.D0/5.D0*(mB(a)*delta(b,2)  + mB(b)*delta(a,2)) - 2.D0/15.D0*mB(2)*delta(a,b)) 
        FB = FB - 4.D0*pi1*(Pes)**(2.D0)*dDdzB(a,b) 

!================================== Set boundary fluxes to zero and calculate them to verify the ghost point calculations==========================

	if (k.eq.1) then
	FBDsum(a,b) = FBDsum(a,b) + abs(FB - uB(2)*DB(a,b))
	FB = 0.D0
	endif

	if (k.eq.NZ) then
	FTDsum(a,b) = FTDsum(a,b) + abs(FT - uT(2)*DTT(a,b))
	FT = 0.D0
	endif
!====================================================================================================================================================

	DFluxSum(a,b) = DFluxSum(a,b) -((FR - FL)/dy + (FT - FB)/dz)
     
        FD(a,b,j,k) = -((FR - FL)/dy + (FT - FB)/dz) - 6.D0*D(a,b,j,k) + beta*(2.D0/5.D0*c(j,k)*E(a,b,j,k))
     
        do aa = 1,2 
        
        FD(a,b,j,k) = FD(a,b,j,k) + 3.D0*beta/7.D0*(E(a,aa,j,k)*D(aa,b,j,k) +  D(a,aa,j,k)*E(aa,b,j,k)) 
	FD(a,b,j,k) = FD(a,b,j,k) + (W(a,aa,j,k)*D(aa,b,j,k) - D(a,aa,j,k)*W(aa,b,j,k))
		
        enddo
      
           enddo
        enddo
        
        
        enddo
        enddo

		FTcsum	   = FTcsum/dble(ny)
		FTmsum     = FTmsum/dble(ny)
		FTDsum     = FTDsum/dble(ny)

		FBcsum     = FBcsum/dble(ny)
		FBmsum     = FBmsum/dble(ny)
		FBDsum     = FBDsum/dble(ny)
        
        return
        end subroutine flux

!*******************************************************************
! File output
!*******************************************************************

	  subroutine fileoutput
		
		Use Glob

		implicit none
		
		double precision sum

! replace by ifile		
    	write(filename1,'(A3,I5.5,A4)') 'con',ifile,'.txt'
        open(unit=5,file=filename1,status='unknown')
         
         do k = 1,nz
            do j = 1,ny

            write(5,*)  yy(j), zzM(k), c(j,k)         

               enddo
            enddo
        close(5)

		write(filename1,'(A3,I5.5,A4)') 'mnp',ifile,'.txt'
        open(unit=5,file=filename1,status='unknown')
         
         do k = 1,nz
            do j = 1,ny

            write(5,*)  yy(j), zzM(k), m(1,j,k), m(2,j,k)         

               enddo
            enddo
        close(5)
        
        write(filename1,'(A3,I5.5,A4)') 'nem',ifile,'.txt'
        open(unit=5,file=filename1,status='unknown')
         
         do k = 1,nz
            do j = 1,ny

            write(5,*)  yy(j), zzM(k), D(1,1,j,k), D(1,2,j,k), D(2,1,j,k), D(2,2,j,k)

               enddo
            enddo
        close(5)
        
        write(filename1,'(A3,I5.5,A4)') 'vel',ifile,'.txt'
        open(unit=5,file=filename1,status='unknown')
         
         do k = 1,nz+1
            do j = 1,ny

            write(5,*)  yy(j), zzU(k), u(1,j,k), u(2,j,k)

               enddo
            enddo
        close(5)

		return
		end subroutine fileoutput
		
! Don't use glob in this subroutine				
!*******************************************************************
! Function fourn:
!     Computes the fast Fourier transform of a multidimensional
!     array of complex numbers.
!     From Numerical Recipes, Press et al., 2001.
!*******************************************************************

      SUBROUTINE fourn(data,nn,ndim,isign) 

      INTEGER isign,ndim,nn(ndim) 
      double precision data(*) 
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2 
      INTEGER ip3,k1,k2,n,nprev,nrem,ntot 
      double precision tempi,tempr 
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp 

      ntot=1 
      do idim=1,ndim 
         ntot=ntot*nn(idim) 
      enddo 
!      print *, 'f1'
      nprev=1 
      do idim=1,ndim
!		 print *,'f2',idim
         n=nn(idim) 
         nrem=ntot/(n*nprev) 
         ip1=2*nprev 
         ip2=ip1*n 
         ip3=ip2*nrem 
         i2rev=1 
         do i2=1,ip2,ip1 
!         print *, 'f3',i2
            if(i2.lt.i2rev) then 
               do i1=i2,i2+ip1-2,2 
!               print *, 'f4',i1
                  do i3=i1,ip3,ip2 
!                  print *, 'f5',i3
                     i3rev=i2rev+i3-i2 
                     tempr=data(i3) 
                     tempi=data(i3+1) 
                     data(i3)=data(i3rev) 
                     data(i3+1)=data(i3rev+1) 
                     data(i3rev)=tempr 
                     data(i3rev+1)=tempi 
                  enddo 
               enddo 
            endif 
            ibit=ip2/2 
 1          if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then 
               i2rev=i2rev-ibit 
               ibit=ibit/2 
               goto 1 
            endif 
            i2rev=i2rev+ibit 
         enddo 
         ifp1=ip1  
 2       if(ifp1.lt.ip2) then 
            ifp2=2*ifp1 
            theta=isign*6.28318530717959d0/(ifp2/ip1) 
            wpr=-2.d0*sin(0.5d0*theta)**2 
            wpi=sin(theta) 
            wr=1.d0 
            wi=0.d0 
            do i3=1,ifp1,ip1 
               do i1=i3,i3+ip1-2,2 
                  do i2=i1,ip3,ifp2 
                     k1=i2 
                     k2=k1+ifp1 
                     tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1) 
                     tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2) 
                     data(k2)=data(k1)-tempr 
                     data(k2+1)=data(k1+1)-tempi
                     data(k1)=data(k1)+tempr 
                     data(k1+1)=data(k1+1)+tempi 
                  enddo 
               enddo 
               wtemp=wr  
               wr=wr*wpr-wi*wpi+wr 
               wi=wi*wpr+wtemp*wpi+wi 
            enddo 
            ifp1=ifp2 
            goto 2 
         endif 
         nprev=n*nprev 
      enddo 
      return 
      END SUBROUTINE FOURN


!*******************************************************************
! Function ran:   
!     Random number generator.
!*******************************************************************

      function ran(idum)

      implicit none
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      double precision ran,am,eps,rnmx
      parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,ia1=40014)
      parameter (ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791)
      parameter (ntab=32,ndiv=1+imm1/ntab,eps=1.2E-7,rnmx=1.-eps)
      integer idum2,j,k,iv(ntab),iy
      save iv,iy,idum2
      data idum2/123456789/, iv/ntab*0/, iy/0/
      
      if (idum.le.0) then
         idum = max(-idum,1)
         idum2 = idum
         do j = ntab+8,1,-1
            k = idum/iq1
            idum = ia1*(idum-k*iq1)-k*ir1
            if (idum.lt.0) idum = idum+im1
            if (j.le.ntab) iv(j) = idum
         enddo
         iy = iv(1)
      endif
      k = idum/iq1
      idum = ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum = idum+im1
      k = idum2/iq2
      idum2 = ia2*(idum2-k*iq2)-k*ir2
      if ( idum2.lt.0) idum2 = idum2+im2
      j = 1+iy/ndiv
      iy = iv(j)-idum2
      iv(j) = idum
      if (iy.lt.1) iy = iy+imm1
      ran = min(am*iy,rnmx)

      return
      end

! Don't use glob in this!!
	   subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
		implicit none 
		integer n
		double precision a(n,n), c(n,n)
		double precision L(n,n), U(n,n), b(n), d(n), x(n)
		double precision coeff
		integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
		L=0.0
		U=0.0
		b=0.0

! step 1: forward elimination
		do k=1, n-1
		do i=k+1,n
			coeff=a(i,k)/a(k,k)
			L(i,k) = coeff
		do j=k+1,n
			a(i,j) = a(i,j)-coeff*a(k,j)
        end do
		end do
		end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
		do i=1,n
			L(i,i) = 1.0
		end do
! U matrix is the upper triangular part of A

		do j=1,n
		   do i=1,j
				U(i,j) = a(i,j)
		   end do
		end do

! Step 3: compute columns of the inverse matrix C
		do k=1,n
			b(k)=1.0
			d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
			do i=2,n
			d(i)=b(i)
				do j=1,i-1
					d(i) = d(i) - L(i,j)*d(j)
				end do
			end do
! Step 3b: Solve Ux=d using the back substitution
			x(n)=d(n)/U(n,n)
			do i = n-1,1,-1
				x(i) = d(i)
				do j=n,i+1,-1
					x(i)=x(i)-U(i,j)*x(j)
				end do
				x(i) = x(i)/u(i,i)
			end do
! Step 3c: fill the solutions x(n) into column k of C
			do i=1,n
				c(i,k) = x(i)
			end do
				b(k)=0.0
	    end do
		end subroutine inverse
