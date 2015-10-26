C**********************************************************************
C Program dilute_confinement_hypersurface:
C     A finite volume simulation of a dilute active suspension
C     between two parallel plates.

C Generates the numerical simulation data in Figure 14
C "Transport of a dilute active suspension in pressure-driven channel flow"
C  B. Ezhilan, D. Saintillan, Journal of Fluid Mechanics, 777 482-522 (2015)

C This code uses the hypersurface boundary condition, ie., it includes steric exclusion
C For simulation details see appendix C and C2

C Please see http://stokeslet.ucsd.edu/publications_files/poiseuille.pdf for details
C of the underlying model 

C Last used : August 10, 2015 Barath Ezhilan

C Serial code, use ifort for faster simulation runs 

C*********************************************************************
C***************************MODULE GLOB*******************************
C*********************************************************************

      Module Glob

      double precision,dimension(:,:,:),allocatable ::c,F,Fo,PPS,PTS,CPP
      double precision,dimension(:),allocatable :: ci,vx,dvx,ux,dux
      double precision,dimension(:,:,:),allocatable :: pp,XST
      double precision,dimension(:,:),allocatable :: delta,ni,FAC
      double precision,dimension(:),allocatable :: theta,thetam,thetap
      double precision,dimension(:),allocatable :: phi,zz,Fi,Ti
      double precision,dimension(:),allocatable :: dlph,dlphm,dlphp,dlth
      double precision,dimension(:),allocatable :: rr,dARt 
      double precision :: H,time,alpha,beta,dfx,dfr,Gamma,tau,V0
      double precision :: pi,dz,dph,dth,dt,mean,meanF,lbac,dr,dAR
      integer,dimension(:),allocatable :: J1,J2
      integer :: P,N,T,itime
 
      end module

C************************************************************************
C   MAIN PROGRAM
C************************************************************************

      program rheology

	use Glob

	implicit none

C Define simulation parameters call param

	call param
C Define initial configuration field
       call initial
C Perform Euler step
       call euler
C Time marching
      do itime = 1,50000000
         print *,itime
       call output

C Update concentration field
       call update
       call intphi(c,ci)

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
      T  = 16
      P  = 16
      N  = 100*T

      allocate(c(P,N,T),F(P,N,T),Fo(P,N,T),CPP(P,N,T))
      allocate(PPS(P,N,T),PTS(P,N,T))
      allocate(ci(N),Fi(N),Ti(N),vx(N),dvx(N),ux(N),dux(N))
      allocate(pp(3,P,T),delta(3,3),XST(3,3,N),ni(3,N))
      allocate(theta(T),thetam(T),thetap(T),phi(P),zz(N))
      allocate(dlphp(T),dlph(T),dlphm(T),dARt(T),dlth(T),rr(T))
      allocate(J1(N),J2(N),FAC(N,T))
      pi  = 4.d0*atan(1.d0)

      lbac = dble(1.d0/100)
      dr = dble(2.d0/T)
      dph = dble(2*pi/P)
      dz  = dble(2.d0/N) ! channel boundaries are -H and +H

	print *,dth,dph,dz

C Boundary points 

C k=n                *******CC*******
C		     ******C**C******
C                    *****C****C*****
C                    ****C******C****
C                    ***C********C***
C                    **C**********C**
C                    *C************C*
C k=n-t/2+1          C**************C

C Interior points    ****************  totally N-T interior points j=1,j=T still need special treatment
C 		     ****************
C 		     ****************
C 		     ****************
C 		     ****************

C k = T/2	     C**************C
C                    *C************C*
C                    **C**********C**
C                    ***C********C***
C                    ****C******C****
C                    *****C****C*****
C                    ******C**C******
C k=1                *******CC*******

	FAC(:,:) = 1.d0

      do k = 1,N
        	zz(k) = dble((2.d0*k-1)/N) - 1.d0
	if (zz(k).le.(-1+lbac)) then
		J1(k) = T/2 + 1 - k
		J2(k) = T/2 + k
		FAC(k,J1(k)) = 0.5d0
		FAC(k,J2(k)) = 0.5d0
	else if	 (zz(k).ge.(1-lbac)) then
	J1(k) = T/2 - n + k
	J2(k) = T/2 + n + 1 - k
	FAC(k,J1(k)) = 0.5d0
	FAC(k,J2(k)) = 0.5d0
	else 
	J1(k) = 1
	J2(k) = T
	endif
      enddo

	do k = 1,N
	   do j = 1,T

		print *, k,j,FAC(k,j)

	   enddo
	enddo

      do j = 1,T

	rr(j)     = dble((2.d0*j-1)/T)-1.d0 
	theta(j)  = acos(rr(j))
	thetam(j) = acos(rr(j)-dr/2.d0) 
	thetap(j) = acos(rr(j)+dr/2.d0)
	dlth(j)   = thetam(j)-thetap(j) 
	dlph(j)   = dph*sin(theta(j))
	dlphm(j)   = dph*sin(thetam(j))
	dlphp(j)   = dph*sin(thetap(j))
	dARt(j)   = dz*dlth(j)/2.d0
     $	+ (rr(j)*dlth(j) - 2*sin(dlth(j)/2.d0)*
     $	cos((thetap(j)+thetam(j))/2.d0))
      enddo

	dAR = dble(4*pi/P/T)

c	print *, rr(1)-dr/2.d0,thetam(1), thetap(T)

      do i = 1,P
	phi(i) = dble((i-1.d0)*dph)
      enddo

C Dimensionless stresslet 

      beta  = 0.d0

C Center of mass and orientation diffusivities

      dfx = 0.33d0
      dfr = 0.5d0
      tau = 10000.d0

C Time marching

      dt   = 0.0000025d0
      time = 0.d0

      V0    = 1.d0
      Gamma = 0.d0

      do k = 1,N

	vx(k)  = Gamma*(1-zz(k)**2.d0)
        dvx(k) = -2.d0*Gamma*zz(k)

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

      if (40000.0*(itime/40000).eq.1.0*itime.OR.itime.eq.1) then 
      ifile = itime/40000

      write(filename1,'(A3,I4.4,A4)') 'con',ifile,'.txt'
      open(unit=5,file=filename1,status='unknown')
      write(5,*) 'Variables = "z","con" '
      
      do k = 1,N
      write(5,*) zz(k),ci(k)
      enddo

      close(5)

      write(filename2,'(A3,I4.4,A4)') 'str',ifile,'.txt'
      open(unit=5,file=filename2,status='unknown')
      write(5,*) 'Variables = "z","S11","S12","S13","S22","S23","S33"'

      do k = 1,N

      write(5,*) zz(k),XST(1,1,k),XST(1,2,k),
     $      XST(1,3,k),XST(2,2,k),XST(2,3,k),XST(3,3,k)

      enddo

      close(5)

      do i = 1,P
	 do k = 1,N
	    do j = J1(k),J2(k)

      CPP(i,k,j) = pp(1,i,j)*c(i,k,j)

	     enddo
	 enddo
      enddo

      call intphi(CPP,ni(1,:))

      do i = 1,P
         do k = 1,N
	    do j = J1(k),J2(k)

      CPP(i,k,j) = pp(2,i,j)*c(i,k,j)

             enddo
         enddo
      enddo

      call intphi(CPP,ni(2,:))

      do i = 1,P
         do k = 1,N
	    do j = J1(k),J2(k)

      CPP(i,k,j) = pp(3,i,j)*c(i,k,j)

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

      write(filename3,'(A3,I4.4,A4)') 'ori',ifile,'.txt'
      open(unit=5,file=filename3,status='unknown')
      write(5,*) 'Variables = "phi","theta", "z=-0.96875","z=-0.46875",
     $"z=0.03125"'
      
      do i = 1,P
	 do j = 1,T

         write(5,*) phi(i),theta(j),c(i,1,j),c(i,9,j),c(i,17,j)

         enddo
      enddo

      close(5)

      write(filename3,'(A3,I4.4,A4)') 'vel',ifile,'.txt'
      open(unit=5,file=filename3,status='unknown')
      write(5,*) 'Variables = "k","HI", "EXT","TOT"'

      do k = 1,N
      write(5,*) zz(k),ux(k),vx(k),ux(k)+vx(k)
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

          print *, 'mean =',mean,meanF


      end

C**********************************************************************
C Subroutine initial:
C     Define intial configuration field
C**********************************************************************

      subroutine initial

      use Glob

      implicit none

      integer i,j,k,l,m

	c = 0.d0

        do i = 1,P
           do k = 1,N
c	      print *, J1(k),J2(k)
	      do j = J1(k),J2(k)

	     c(i,k,j) = 1.d0

              enddo
           enddo
        enddo

	  call intphi(c,ci)
	  do k = 1,N
		c(:,k,:) = c(:,k,:)/ci(k)
	  enddo
	  call intphi(c,ci)

	  mean = 0.d0

	do k = 1,N
	 mean = mean + ci(k)
	enddo

	  print *, mean
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

	integer i,k,j,l,m
	double precision XS(P,N,T),XSi(N),S13it,S13ii(N),pdot(3)

	XS = 0.d0

	do m = 1,3
	do l = 1,3

	do i = 1,P
	do k = 1,N
	do j = J1(k),J2(k)

C	Extra stress
	XS(i,k,j) = c(i,k,j)*pp(l,i,j)*pp(m,i,j)

	enddo
	enddo
	enddo

	call intphi(XS,XSi)
C     XST extra stress tensor
	XST(l,m,:) = XSi(:)

	enddo
	enddo

	S13it = 0.d0 ! \int_{-1}^{1} \Sigma_{13} dz
	S13ii(:) = 0.d0 !\int_{-1}^{z} \Sigma_{13} dz

	do k = 1,N

	   S13it = S13it + XST(1,3,k)*dz 
 	   S13ii(k) = S13it

	enddo

	do k = 1,N

	   ux(k) = beta*(S13ii(k)- 0.5d0*(zz(k)+1.d0)*S13it)
	   dux(k)= beta*(XST(1,3,k) - 0.5d0*S13it)

	enddo

	  ux(1) = ux(1)/2.d0

	do k = 2,N
	   ux(k) = (ux(k-1)+ux(k))/2.d0
	enddo

	   do i = 1,P
	   do j = 1,T
	   do k = 1,N

	   do l = 1,3

	   pdot(l) = (delta(l,1)-pp(l,i,j)*pp(1,i,j))*
     $		      (dux(k)+dvx(k))*pp(3,i,j)

	   enddo

C phi component of the pdot vector
	PPS(i,k,j)=(-sin(phi(i))*pdot(1)+cos(phi(i))*pdot(2))*
     $	           c(i,k,j)
C theta component of the pdot vector
	PTS(i,k,j)= (cos(theta(j))*cos(phi(i))*pdot(1) +
     $         cos(theta(j))*sin(phi(i))*pdot(2) -
     $         sin(theta(j))*pdot(3))*c(i,k,j)

	   enddo
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

      integer i,k,j
      double precision cll(P,N,T),clli(N)

      clli(:) = 0.d0

      do i = 1,P
         do k = 1,N
            do j = J1(k),J2(k)
               clli(k) = clli(k) + cll(i,k,j)*dAR*FAC(k,j)
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

      integer i,k,j

C Get fluxes
	call pdot_components
	call flux
C Update configuration field using explicit Euler scheme
C and copy flux

         do i = 1,P
            do k = 1,N
	       do j = J1(k),J2(k)

               c(i,k,j) = c(i,k,j) - dt*F(i,k,j)
               Fo(i,k,j)= F(i,k,j)

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

      integer i,k,j

C Get fluxes

      call pdot_components
      call flux

C Update configuration field using explicit Euler scheme
C and copy flux

         do i = 1,P
            do k = 1,N

		do j = J1(k),J2(k)

                 c(i,k,j) = c(i,k,j) - 
     $         0.5d0*dt*(3.d0*F(i,k,j)-Fo(i,k,j))

               Fo(i,k,j) = F(i,k,j)

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
	
	   do k = 1,N
	  if (zz(k).le.(-1+lbac)) then 
	      call fluxBL(i,k)
	      call fluxBR(i,k)
	  else if (zz(k).ge.(1-lbac)) then
	      call fluxTL(i,k)
              call fluxTR(i,k)

	  else 
	      call flux1(i,k)
	      call fluxT(i,k)
          endif
		do j = J1(k)+1,J2(k)-1
            	      call fluxINTERIOR(i,k,j)
		enddo
	      enddo	
	enddo

	return
	end

C********************************************************************
C Discretization for the top left boundary point
C********************************************************************

      subroutine fluxBL(i,k)

      Use Glob

      implicit none

      integer i,ip,im,j,jp,jm,k
      double precision INS

           ip = i+1
           if (ip.gt.P) ip = 1
           im = i-1
           if (im.lt.1) im = P

           j  = J1(k)
           jp = j+1
           jm = j-1

           F(i,k,j) = F(i,k,j) + 2.d0*(
     $     V0*cos(theta(j))*(c(i,k+1,j)+c(i,k,j))/(2.d0*dz) -
     $     dfx*(c(i,k+1,j)-c(i,k,j))/(dz)**(2.d0)
     $       + (PPS(ip,k,j)-PPS(im,k,j))*dARt(j)/(2.d0*dAR*dz)
     $       - ((PTS(i,k,jp)+PTS(i,k,j))*dlphp(j))/(2.d0*dAR)
     $  - dfr/dAR*((c(ip,k,j)-2.d0*c(i,k,j)+c(im,k,j))*dARt(j)/(dph*dz)
     $            +(c(i,k,jp)-c(i,k,j))*((dlphp(j))**2.d0)/(dAR)))

	return
	end

C********************************************************************
C Discretization for the bottom right boundary point (j=J2(k))
C********************************************************************
      subroutine fluxBR(i,k)

      Use Glob

      implicit none

      integer i,ip,im,j,jp,jm,k
      double precision INS

           ip = i+1
           if (ip.gt.P) ip = 1
           im = i-1
           if (im.lt.1) im = P

	   j  = J2(k)
	   jp = j+1
           jm = j-1

           F(i,k,j) = F(i,k,j) + 2.d0*(
     $     V0*cos(theta(j))*(c(i,k+1,j)+c(i,k,j))/(2.d0*dz) -
     $     dfx*(c(i,k+1,j)-c(i,k,j))/(dz)**(2.d0)
     $       + (PPS(ip,k,j)-PPS(im,k,j))*dARt(j)/(2.d0*dAR*dz)
     $       + ((PTS(i,k,j)+PTS(i,k,jm))*dlphm(j))/(2.d0*dAR)
     $  - dfr/dAR*((c(ip,k,j)-2.d0*c(i,k,j)+c(im,k,j))*dARt(j)/(dph*dz)
     $            -(c(i,k,j)-c(i,k,jm))*((dlphm(j))**2.d0)/(dAR)))

	return
        end

C********************************************************************
C Discretization for the Top Left boundary point (j=J1(k))
C********************************************************************
      
      subroutine fluxTL(i,k)

      Use Glob

      implicit none

      integer i,ip,im,j,jp,jm,k
      double precision INS

           ip = i+1
           if (ip.gt.P) ip = 1
           im = i-1
           if (im.lt.1) im = P

	   j  = J1(k)
	   jp = j+1
           jm = j-1

           F(i,k,j) = F(i,k,j) + 2.d0*(-
     $     V0*cos(theta(j))*(c(i,k-1,j)+c(i,k,j))/(2.d0*dz) +
     $     dfx*(c(i,k,j)-c(i,k-1,j))/(dz)**(2.d0)
     $       + (PPS(ip,k,j)-PPS(im,k,j))*dARt(j)/(2.d0*dAR*dz)
     $       - ((PTS(i,k,jp)+PTS(i,k,j))*dlphp(j))/(2.d0*dAR)
     $  - dfr/dAR*((c(ip,k,j)-2.d0*c(i,k,j)+c(im,k,j))*dARt(j)/(dph*dz)
     $            +(c(i,k,jp)-c(i,k,j))*((dlphp(j))**2.d0)/(dAR)))

	return
        end

C********************************************************************
C Discretization for Top right boundary point (j=J2(k))
C********************************************************************
      
      subroutine fluxTR(i,k)

      Use Glob

      implicit none

      integer i,ip,im,j,jp,jm,k
      double precision INS

           ip = i+1
           if (ip.gt.P) ip = 1
           im = i-1
           if (im.lt.1) im = P

	      j  = J2(k)
	      jp = j+1
              jm = j-1

           F(i,k,j) = F(i,k,j) + 2.d0*(-
     $     V0*cos(theta(j))*(c(i,k-1,j)+c(i,k,j))/(2.d0*dz) +
     $     dfx*(c(i,k,j)-c(i,k-1,j))/(dz)**(2.d0)
     $       + (PPS(ip,k,j)-PPS(im,k,j))*dARt(j)/(2.d0*dAR*dz)
     $       + ((PTS(i,k,j)+PTS(i,k,jm))*dlphm(j))/(2.d0*dAR)
     $  - dfr/dAR*((c(ip,k,j)-2.d0*c(i,k,j)+c(im,k,j))*dARt(j)/(dph*dz)
     $            -(c(i,k,j)-c(i,k,jm))*((dlphm(j))**2.d0)/(dAR)))

	return
        end

C********************************************************************
C Discretization when the |Z| < 1-l and for \theta = dth/2 (j = 1) boundary point
C********************************************************************
      subroutine flux1(i,k)

      Use Glob

      implicit none

      integer i,ip,im,j,jp,jm,k
      double precision INS

           ip = i+1
           if (ip.gt.P) ip = 1
           im = i-1
           if (im.lt.1) im = P

	      j  = 1
              jp = j+1

           F(i,k,j) = F(i,k,j) +
     $     V0*cos(theta(j))*(c(i,k+1,j)-c(i,k-1,j))/(2.d0*dz) -
     $     dfx*(c(i,k+1,j)-2.d0*c(i,k,j)+c(i,k-1,j))/(dz)**(2.d0)
     $       + (PPS(ip,k,j)-PPS(im,k,j))*dlth(j)/(2.d0*dAR)
     $       - (PTS(i,k,jp)+PTS(i,k,j))*dlphp(j)/(2.d0*dAR)
     $  - dfr/dAR*((c(ip,k,j)-2.d0*c(i,k,j)+c(im,k,j))*dlth(j)/dph
     $            +(c(i,k,jp)-c(i,k,j))*(dlphp(j))**2.d0/(dAR))

	return
        end

C********************************************************************
C Discretization when the |Z| < 1-l and for \theta = pi - dth/2 (j = T) boundary point
C********************************************************************

      subroutine fluxT(i,k) 

      Use Glob

      implicit none

      integer i,ip,im,j,jp,jm,k
      double precision INS

           ip = i+1
           if (ip.gt.P) ip = 1
           im = i-1
           if (im.lt.1) im = P

              j  = T
              jm = j-1


           F(i,k,j) = F(i,k,j) +
     $     V0*cos(theta(j))*(c(i,k+1,j)-c(i,k-1,j))/(2.d0*dz) -
     $     dfx*(c(i,k+1,j)-2.d0*c(i,k,j)+c(i,k-1,j))/(dz)**(2.d0)
     $       + (PPS(ip,k,j)-PPS(im,k,j))*dlth(j)/(2.d0*dAR)
     $       + (PTS(i,k,jm)+PTS(i,k,j))*dlphm(j)/(2.d0*dAR)
     $  - dfr/dAR*((c(ip,k,j)-2.d0*c(i,k,j)+c(im,k,j))*dlth(j)/dph
     $           -(c(i,k,j)-c(i,k,jm))*(dlphm(j))**2.d0/(dAR))

	return
        end

C********************************************************************
C Discretization for the interior boundary points
C********************************************************************
      subroutine fluxINTERIOR(i,k,j)

      Use Glob

      implicit none

      integer i,ip,im,j,jp,jm,k
      double precision INS

           ip = i+1
           if (ip.gt.P) ip = 1
           im = i-1
           if (im.lt.1) im = P

              jp = j+1
              jm = j-1

           F(i,k,j) = F(i,k,j) +
     $     V0*cos(theta(j))*(c(i,k+1,j)-c(i,k-1,j))/(2.d0*dz) -
     $     dfx*(c(i,k+1,j)-2.d0*c(i,k,j)+c(i,k-1,j))/(dz)**(2.d0)
     $       + (PPS(ip,k,j)-PPS(im,k,j))*dlth(j)/(2.d0*dAR)       
     $       - (PTS(i,k,jp)*dlphp(j)-PTS(i,k,jm)*dlphm(j)
     $          +PTS(i,k,j)*(dlphp(j)-dlphm(j)))/(2.d0*dAR)
     $  - dfr/dAR*((c(ip,k,j)-2.d0*c(i,k,j)+c(im,k,j))*dlth(j)/dph
     $            +(c(i,k,jp)*(dlphp(j))**2.d0 
     $		    -c(i,k,j)*((dlphm(j))**2.d0 + (dlphp(j))**2.d0)
     $		    +c(i,k,jm)*(dlphm(j))**2.d0)/(dAR))

	return
        end

