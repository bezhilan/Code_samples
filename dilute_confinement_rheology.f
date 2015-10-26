C**********************************************************************
C Program dilute_confinement_rheology:
C     A finite volume simulation of a dilute active suspension
C     between two parallel plates. Extension of the dilute_confinement
C     to include rheology calculations

C This code uses the wall boundary condition, ie., it ignores steric exclusion
C For computational details see appendix C and C1 of
C "Transport of a dilute active suspension in pressure-driven channel flow"
C  B. Ezhilan, D. Saintillan, Journal of Fluid Mechanics, 777 482-522 (2015)

C Last used : August 10, 2015 Barath Ezhilan

C Serial code, use ifort for faster simulation runs 
C*********************************************************************
C***************************MODULE GLOB*******************************
C*********************************************************************

      Module Glob

      double precision,dimension(:,:,:),allocatable ::c,F,Fo,PPS,PTS,CPP
      double precision,dimension(:),allocatable :: ci,Myz,Myyzz
      double precision,dimension(:),allocatable :: vx,dvx,ux,dux,eta,exs
      double precision,dimension(:,:,:),allocatable :: pp,XST
      double precision,dimension(:,:),allocatable :: delta,ni
      double precision,dimension(:),allocatable :: theta,thetam,thetap
      double precision,dimension(:),allocatable :: phi,zz,Fi,Ti
      double precision,dimension(:),allocatable :: dph2,dphm,dphp,AR
      
      double precision :: H,time,alpha,beta,dfx,dfr,Gamma,tau,V0
      double precision :: pi,dz,dph,dth,dt,mean,meanF,cmax
      double precision :: Pe_s,Pe_f

      integer :: P,T,N
      integer*8 :: itime
 
      end module

C************************************************************************
C   MAIN PROGRAM
C************************************************************************

      program rheology

	use Glob

	implicit none

	integer k

C Define simulation parameters
       call param

C Define initial configuration field
       call initial

C Perform Euler step
       call euler

C Time marching
      do itime = 1,5000000
         print *,'itime =',itime
       call output

C Update concentration field
       call update
       call intphi(c,ci)

      cmax = abs(ci(1))

      do k = 2,N
        
         if(cmax.lt.abs(ci(k))) then
            cmax = abs(ci(k))
         endif

      enddo

      if (cmax.gt.50) then
      print*,'Stop - solution blowing up at step',itime
      exit
      endif

	print *,'cmax =',cmax

      enddo

      end

C**********************************************************************
C Subroutine param:
C     Define simulation parameters
C**********************************************************************

      subroutine param

      use Glob

      implicit none

      integer i,j,k,l,m

      H = 1.d0

C Number of points in x, y and phi
      N  = 100
      T  = 16
      P  = 32

      allocate(c(P,T,N),F(P,T,N),Fo(P,T,N),CPP(P,T,N))
      allocate(PPS(P,T,N),PTS(P,T,N))
      allocate(ci(N),Myz(N),Myyzz(N))
      allocate(Fi(N),Ti(N),vx(N),dvx(N),ux(N),dux(N),eta(N),exs(N))
      allocate(pp(3,P,T),delta(3,3),XST(3,3,N),ni(3,N))
      allocate(theta(T),thetam(T),thetap(T),phi(P),zz(N))
      allocate(dphp(T),dph2(T),dphm(T),AR(T))

      pi  = 4.d0*atan(1.d0)

      dth = dble(pi/T)
      dph = dble(2*pi/P)
      dz  = dble(2*H/N) ! channel boundaries are -H and +H

	print *,dth,dph,dz

      do k = 1,N
	zz(k) = (k-0.5d0)*dz - 1.d0*H
      enddo

      do j = 1,T

	theta(j) = dble((j-0.5d0)*dth)
	thetam(j) = theta(j)-dth/2.d0
	thetap(j) = theta(j)+dth/2.d0
	dph2(j) = dph*sin(theta(j))
	dphm(j)= dph*sin(thetam(j))
	dphp(j)= dph*sin(thetap(j))
	AR(j) = dph*(cos(thetam(j))-cos(thetap(j)))

      enddo

	print *, thetam(1), thetap(T)

      do i = 1,P
	phi(i) = dble((i-1.d0)*dph)
      enddo

C Only variable
      Pe_s = 0.5d0
      Pe_f = 0.25d0
C Center of mass and orientation diffusivities

      dfx = (Pe_s)**(2.d0)/3.d0 ! translational diffusion
      dfr = 0.5d0 !rotational diffusion
      tau = 10000.d0 ! tumbling rate

C Time marching

      dt   = 0.0001d0
      time = 0.d0

      V0    = Pe_s
      Gamma = Pe_f/2.d0

      do k = 1,N

!	vx(k)  = Gamma*(1-zz(k)**2.d0)
        dvx(k) = -Gamma*zz(k)

      enddo

	do i = 1,P
           do j = 1,T

	pp(1,i,j) = sin(theta(j))*cos(phi(i))
	pp(2,i,j) = sin(theta(j))*sin(phi(i))
	pp(3,i,j) = cos(theta(j))

	   enddo
	enddo

	do m = 1,3
         do l = 1,3

           delta(l,m) = 0.d0
         enddo

          delta(m,m) = 1.d0
      enddo

      return
      end

C**********************************************************************
C Subroutine output:
C     Define simulation parameters
C**********************************************************************

      subroutine output

      Use Glob

      implicit none

      character*32 filename1,filename2,filename3
      integer ifile,i,j,k

      if (100000.0*(itime/100000).eq.1.0*itime.OR.itime.eq.1) then 
      ifile = itime/100000

      call visc

      write(filename1,'(A3,I4.4,A4)') 'vis',ifile,'.txt'
      open(unit=5,file=filename1,status='unknown')
      write(5,*) 'Variables = "z","exs","eta","ux" '

      do k = 1,N
      write(5,*) zz(k),exs(k),eta(k),ux(k)
      enddo

      close(5)


      write(filename1,'(A3,I4.4,A4)') 'con',ifile,'.txt'
      open(unit=5,file=filename1,status='unknown')
      write(5,*) 'Variables = "z","con" '
      
      do k = 1,N
      write(5,*) zz(k),ci(k)
      enddo

      close(5)

      do i = 1,P
	 do j = 1,T
	    do k = 1,N

      CPP(i,j,k) = pp(1,i,j)*c(i,j,k)

	     enddo
	 enddo
      enddo

      call intphi(CPP,ni(1,:))

      do i = 1,P
         do j = 1,T
            do k = 1,N

      CPP(i,j,k) = pp(2,i,j)*c(i,j,k)

             enddo
         enddo
      enddo

      call intphi(CPP,ni(2,:))

      do i = 1,P
         do j = 1,T
            do k = 1,N

      CPP(i,j,k) = pp(3,i,j)*c(i,j,k)

             enddo
         enddo
      enddo

      call intphi(CPP,ni(3,:))

      write(filename3,'(A3,I4.4,A4)') 'mnp',ifile,'.txt'
      open(unit=5,file=filename3,status='unknown')
      write(5,*) 'Variables = "k","n1", "n2","n3"'

      do k = 1,N
      write(5,*) zz(k),ni(1,k),ni(2,k),ni(3,k)
      enddo

      close(5)

      call pdot_components

      write(filename3,'(A3,I4.4,A4)') 'str',ifile,'.txt'
      open(unit=5,file=filename3,status='unknown')
      write(5,*) 'Variables = "k","S11", "S22","S33",
     $"S12","S13","S23"'

      do k = 1,N

      write(5,*) zz(k),XST(1,1,k),XST(2,2,k),XST(3,3,k),
     $ XST(1,2,k),XST(1,3,k),XST(2,3,k)

      enddo

      close(5)


      write(filename3,'(A3,I4.4,A4)') 'ori',ifile,'.txt'
      open(unit=5,file=filename3,status='unknown')
      write(5,*) 'Variables = "phi","theta", "z=z1","z=z2",
     $"z=z3","z=z4","z=z5","z=z6"'

      do i = 1,P
	 do j = 1,T

      write(5,*) phi(i),theta(j),c(i,j,1),c(i,j,11),
     $	         c(i,j,21),c(i,j,31),c(i,j,41),c(i,j,51)

	 enddo
      enddo

      close(5)

	endif

          mean = 0.d0
	  meanF = 0.d0

        do k = 1,N

         mean = mean + ci(k)
	 meanF = meanF + Fi(k)

        enddo

          mean = mean/dble(N)
	  meanF = meanF/dble(N)

          print *, 'mean =',mean,'meanF=',meanF


      end

C**********************************************************************
C Subroutine initial:
C     Define intial configuration field
C**********************************************************************

      subroutine initial

      use Glob

      implicit none

      integer i,j,k,l,m

        do i = 1,P
           do j = 1,T
              do k = 1,N
           
	     c(i,j,k) = 1.d0

              enddo
           enddo
        enddo
	  call intphi(c,ci)

	  mean = 0.d0

	do k = 1,N
	 mean = mean + ci(k)
	enddo

	  mean = mean/dble(N)

	  print *, 'mean =',mean

          c = c/mean

	call intphi(c,ci)
 
	ux(:)  = 0.d0
	dux(:) = 0.d0

	return
	end

C**********************************************************************
C Subroutine pdot_components
C     Find the velocity from the stresses
C*********************************************************************

	subroutine pdot_components

	use Glob

	implicit none

	integer i,j,k,l,m
	double precision XS(P,T,N),XSi(N),S13it,S13ii(N),pdot(3)

	XS = 0.d0

	do m = 1,3
	do l = 1,3

	do i = 1,P
	do j = 1,T
	do k = 1,N

C	Extra stress
	XS(i,j,k) = c(i,j,k)*(pp(l,i,j)*pp(m,i,j)-delta(l,m)/3.d0)

	enddo
	enddo
	enddo

	call intphi(XS,XSi)
C     XST extra stress tensor
	XST(l,m,:) = XSi(:)

	enddo
	enddo

	return
	end
    
C**********************************************************************
C Subroutine intphi:
C     Calculates the integral over phi using trapezoidal rule
C**********************************************************************

      subroutine intphi(cll,clli)

      Use Glob

      implicit none

      integer i,j,k
      double precision cll(P,T,N),clli(N)

      clli(:) = 0.d0

      do i = 1,P
         do j = 1,T
            do k = 1,N
               clli(k) = clli(k) + cll(i,j,k)*AR(j) !area function
            enddo
         enddo
      enddo

      return
      end

C**********************************************************************
C Subroutine euler:
C     Perform initial step using Euler algorithm
C**********************************************************************

      subroutine euler

      use Glob

      implicit none

      integer i,j,k,l
      double precision pdot(3)

C Get fluxes
	PPS = 0.d0
	PTS = 0.d0

          do i = 1,P
          do j = 1,T
          do k = 1,N

          do l = 1,3

          pdot(l) = (delta(l,1)-pp(l,i,j)*pp(1,i,j))*
     $               (dvx(k))*pp(3,i,j)

          enddo

C phi component of the pdot vector
       PPS(i,j,k)=(-sin(phi(i))*pdot(1)+cos(phi(i))*pdot(2))*
     $            c(i,j,k)
C theta component of the pdot vector
       PTS(i,j,k)= (cos(theta(j))*cos(phi(i))*pdot(1) +
     $         cos(theta(j))*sin(phi(i))*pdot(2) -
     $         sin(theta(j))*pdot(3))*c(i,j,k)

          enddo
          enddo
          enddo


	call flux

C Update configuration field using explicit Euler scheme
C and copy flux

         do i = 1,P
            do j = 1,T
	       do k = 1,N

               c(i,j,k) = c(i,j,k) - dt*F(i,j,k)
               Fo(i,j,k)= F(i,j,k)

            enddo
         enddo
      enddo

C Update time

	call intphi(c,ci)

      time = time + dt

      return
      end

C**********************************************************************
C Subroutine update:
C     Advance configuration field using 2nd order Adams-Bashforth
C**********************************************************************

      subroutine update

      use Glob

      implicit none

      integer i,j,k,l
      double precision pdot(3)

C Get fluxes

      PPS = 0.d0
      PTS = 0.d0

          do i = 1,P
          do j = 1,T
          do k = 1,N

          do l = 1,3

          pdot(l) = (delta(l,1)-pp(l,i,j)*pp(1,i,j))*
     $               (dvx(k))*pp(3,i,j)

          enddo

C phi component of the pdot vector
       PPS(i,j,k)=(-sin(phi(i))*pdot(1)+cos(phi(i))*pdot(2))*
     $            c(i,j,k)
C theta component of the pdot vector
       PTS(i,j,k)= (cos(theta(j))*cos(phi(i))*pdot(1) +
     $         cos(theta(j))*sin(phi(i))*pdot(2) -
     $         sin(theta(j))*pdot(3))*c(i,j,k)

          enddo
          enddo
          enddo

 
      call flux

C Update configuration field using explicit Euler scheme
C and copy flux

         do i = 1,P
            do j = 1,T
		do k = 1,N

                 c(i,j,k) = c(i,j,k) - 
     $         0.5d0*dt*(3.d0*F(i,j,k)-Fo(i,j,k))

               Fo(i,j,k) = F(i,j,k)

            enddo
         enddo
      enddo

C Update time

	call intphi(c,ci)
	call intphi(F,Fi)
      time = time + dt

      return
      end

C********************************************************************

      subroutine flux

      use Glob

      implicit none

      integer i,ip,im,j,jp,jm,k
      double precision INS

C Initialize flux

               F = 0.d0

	call intphi(c,ci)

	do i = 1,P
	
	   ip = i+1
	   if (ip.gt.P) ip = 1
	   im = i-1
	   if (im.lt.1) im = P

           do j = 1,T
	      jp = j+1
	      jm = j-1
!	      if (jm.lt.1) jp = T+1     !jp at j=T,jm at j=1 elements are anyway multiplied by sin(0) and sin(pi). 

	   !k = 1 

       F(i,j,1) = F(i,j,1) +
     $	          V0*cos(theta(j))*(c(i,j,1)+c(i,j,2))/(2.d0*dz) -
     $		  dfx*(c(i,j,2)-c(i,j,1))/(dz)**(2.d0)

	   !k = N

	F(i,j,N) = F(i,j,N) - 
     $		    V0*cos(theta(j))*(c(i,j,N)+c(i,j,N-1))/(2.d0*dz) +
     $		    dfx*(c(i,j,N)-c(i,j,N-1))/(dz)**(2.d0)

	      do k = 2,N-1

		   F(i,j,k) = F(i,j,k) + 
     $	   V0*cos(theta(j))*(c(i,j,k+1)-c(i,j,k-1))/(2.d0*dz) -
     $	   dfx*(c(i,j,k+1)-2.d0*c(i,j,k)+c(i,j,k-1))/(dz)**(2.d0)

	      enddo
	    enddo
	  enddo

	do i = 1,P

           ip = i+1
           if (ip.gt.P) ip = 1
           im = i-1
           if (im.lt.1) im = P

	   j = 1
	   jp = j+1

             do k = 1,N

        F(i,j,k)   = F(i,j,k)
     $       + (PPS(ip,j,k)-PPS(im,j,k))*dth/(2.d0*AR(j))       
     $       + (PTS(i,jp,k)*dphp(j)
     $       +  PTS(i,j,k)*dphp(j))/(2.d0*AR(j))
     $  - dfr/AR(j)*((c(ip,j,k)-2.d0*c(i,j,k)+c(im,j,k))*dth/dph2(j)
     $              +(c(i,jp,k)*dphp(j)-c(i,j,k)*dphp(j)
     $               )/dth)
!     $              +(c(i,j,k)-ci(k)/(4.d0*pi))/tau
 

              enddo


	   do j = 2,T-1
	     
 	      jp = j+1
              jm = j-1

	      do k = 1,N

	F(i,j,k)   = F(i,j,k) 
     $	     + (PPS(ip,j,k)-PPS(im,j,k))*dth/(2.d0*AR(j))	
     $	     + (PTS(i,jp,k)*dphp(j)-PTS(i,jm,k)*dphm(j)
     $	     +  PTS(i,j,k)*(dphp(j)-dphm(j)))/(2.d0*AR(j)) 
     $	- dfr/AR(j)*((c(ip,j,k)-2.d0*c(i,j,k)+c(im,j,k))*dth/dph2(j)
     $              +(c(i,jp,k)*dphp(j)-c(i,j,k)*(dphm(j)+dphp(j))
     $               +c(i,jm,k)*dphm(j))/dth)
!     $		     +(c(i,j,k)-ci(k)/(4.d0*pi))/tau
	

	      enddo
	   enddo

	   j = T
	   jm = j-1

	   do k = 1,N

        F(i,j,k)   = F(i,j,k)
     $       + (PPS(ip,j,k)-PPS(im,j,k))*dth/(2.d0*AR(j))       
     $       + (-PTS(i,jm,k)*dphm(j)
     $       +  PTS(i,j,k)*(-dphm(j)))/(2.d0*AR(j))
     $  - dfr/AR(j)*((c(ip,j,k)-2.d0*c(i,j,k)+c(im,j,k))*dth/dph2(j)
     $              +(-c(i,j,k)*(dphm(j))
     $              +c(i,jm,k)*dphm(j))/dth)
!     $		    +(c(i,j,k)-ci(k)/(4.d0*pi))/tau
	      enddo
	enddo 

! tumble first checking if it sums up to zero

	call intphi(F,Fi)

      return
      end

C**********************************************************************
C Subroutine visc:
C     Perform initial step using Euler algorithm
C**********************************************************************

      subroutine visc

      use Glob

      implicit none

      integer i,j,k,l
      double precision a1,b1,S13it,S13ii(N)

      a1 = 0.54d0
      b1 = 0.12d0

	eta(:) = 0.d0
	exs(:) = 0.d0

! calculate Myz

	do i = 1,P
         do j = 1,T
            do k = 1,N

             CPP(i,j,k) = pp(1,i,j)*pp(3,i,j)*c(i,j,k)

             enddo
         enddo
      enddo

      call intphi(CPP,Myz)

! calculate Myyzz

        do i = 1,P
         do j = 1,T
            do k = 1,N

             CPP(i,j,k) = ((pp(1,i,j)*pp(3,i,j))**(2.0))*c(i,j,k)

             enddo
         enddo
      enddo

      call intphi(CPP,Myyzz)

	do k = 1,N

	exs(k) = a1*Myz(k) + b1*Pe_f*zz(k)*Myyzz(k)

	if(Pe_f.gt.0) then
		eta(k) = exs(k)/(Pe_f*zz(k))
	endif

	enddo

        S13it = 0.d0 ! \int_{-1}^{1} \Sigma_{13} dz
        S13ii(:) = 0.d0 !\int_{-1}^{z} \Sigma_{13} dz

        do k = 1,N

           S13it = S13it + exs(k)*dz
           S13ii(k) = S13it - exs(k)*dz/2.d0 ! To ensure symmetric, the endpoint of the integral have a factor of half

        enddo

        do k = 1,N

           ux(k) = -1.d0*(S13ii(k)- 0.5d0*(zz(k)+1.d0)*S13it) ! see roman stocker new notes

        enddo


	return
	end
