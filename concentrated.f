C**********************************************************************
C Program concentrated:
C     Simulates the dynamics in a 3D suspension of concentrated
C     active suspension in periodic domain using a continuum model.

C Generates the numerical simulation data in 
C "Instabilities and nonlinear dynamics of concentrated active suspensions"
C  B. Ezhilan, M. J. Shelley, D. Saintillan, Physics of Fluids, 25 070607 (2013)

C Please see http://stokeslet.ucsd.edu/publications_files/concentrated.pdf for details
C of the underlying model 

C Last used : August 10, 2015 Barath Ezhilan

C MPI Parallel code usually run on 64/128 processors (spatial direction z is parallelized)
C Download the appropriate MPIF and lfftw libraries before compiling 
C To compile use "mpif90 -O3 filename.f -lfftw_mpi -lfftw -lm"
C To execute use "mpirun -np 8 a.out"

C I am not sure if this is the final version of the code used in the paper
C**********************************************************************



C*********************************************************************
C 				MODULE GLOB
C********************************************************************

      Module Glob

      real,dimension(:,:,:,:,:),allocatable :: c,bEW,M,MM,ing,ft 
      ! c = \psi; bEW = \beta E + W; M = second moment \int \psi p p dp; MM = M\cdot M (see eq.(22) of Ezhilan et al. POF 2013)
      ! ing is a dummy variable for receiving (3,3,nx,ny,nz) data from other processors.
      
      real,dimension(:,:,:,:,:),allocatable :: E,f1
      real,dimension(:,:,:,:,:),allocatable :: F,Fo

      real,dimension(:,:,:),allocatable :: cremesh,ppE,ppM


      real,dimension(:,:,:),allocatable     :: f1tot,f2tot,f3tot

      real,dimension(:,:),allocatable       :: cl,clt

      real,dimension(:,:),allocatable       :: uxwait,uywait,uzwait


      real,dimension(:,:,:,:),  allocatable :: fcompx,fcompy,cwait
      real,dimension(:,:,:,:),  allocatable :: EMwait
      real,dimension(:,:,:,:),  allocatable :: fcompz,ukx,uky,ukz,divMM

      real,dimension(:,:,:,:),  allocatable :: csend
 
      real,dimension(:,:,:),    allocatable :: ci,ux,uy,uz,S11i,S12i
      real,dimension(:,:,:),    allocatable :: fx,fy,fz 

      real,dimension(:),    allocatable     :: theta,phi 
      
      real :: Lx,Ly,Lz,time,alpha,beta,U0,nu,dfx,dfr,gamma,dt,dy2
      real :: pi,dx,dy,dz,dph,dth

      integer :: nx,ny,nz,nph,nth,itime
      integer :: totnx,totny,totnz,req,req1,req2,req3

      real,dimension(:,:,:),allocatable :: ciphi

      integer :: Scoords(1),Dcoords(1),rootR,rankR,rankL
      integer :: rankF,rankO,coords(1),rank,ierr,reorder,tag,cartcomm
      integer :: periods(1),dims(1),numtasks,rc
      integer :: coordsR(1),coordsL(1) 
      integer :: Imin,Imax,Jmin,Jmax,Kmin,Kmax

      real,dimension(:,:,:),allocatable :: fxT,fyT,fzT
      real,dimension(:,:,:,:),allocatable :: ftotcompx,ftotcompy
      real,dimension(:,:,:,:),allocatable :: ftotcompz,fTcompx,fTcompy
      real,dimension(:,:,:,:),allocatable :: fTcompz,ukxT,ukyT,ukzT
      real,dimension(:,:,:,:),allocatable :: ukxtot,ukytot,ukztot
      real,dimension(:,:,:),allocatable :: ciT,LET,pwd
      real,dimension(:,:,:),allocatable :: pxiT,pxi,pyi,pzi
      real,dimension(:,:,:),allocatable :: pyiT,pziT,pow,LEI
      real,dimension(:,:,:),allocatable :: uxT,uyT,uzT

      integer,allocatable,dimension(:) :: stat

      integer*8 :: plan,iplan

      integer :: lxs,lnx,lnyt,lyst,lsize

      integer :: fftw_fwd, fftw_bkwd,fftw_est,fftw_norm_order

      double complex,dimension(:),allocatable :: fxft,fyft,fzft,work
      double complex,dimension(:),allocatable :: uxft,uyft,uzft


C      real :: MPI_REAL

      end module


C*********************************************************************
C 				MAIN PROGRAM
C********************************************************************


      program continuum
    
      use Glob

      implicit none

      include 'mpif.h'

C      integer stat(MPI_STATUS_SIZE)

      allocate (stat(MPI_STATUS_SIZE))

      call MPI_INIT(ierr)
      if (ierr .ne. MPI_SUCCESS) then
         print *,'Error starting MPI program. Terminating.'
         call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
      end if

      call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

      periods(1) = 1

      dims(1)    = 64

      reorder    = 0

      call MPI_CART_CREATE(MPI_COMM_WORLD,1,dims,periods,reorder,
     $                     cartcomm,ierr)

      call MPI_COMM_RANK(cartcomm, rank, ierr)

      call MPI_CART_COORDS(cartcomm,rank,1,coords,ierr)

      print*,'rank = ',rank


      nx = 64 
      ny = 64
      nz = 1

      nz = nz + 2

         
C Define simulation parameters

      call param

C Start MPI

C      Call Init_MPI

C Define initial configuration field

      call initial


C Perform Euler step
      call euler
C Time marching
      do itime = 1,1000000


C Update concentration field
         call update


         call stressmagnitude

      enddo

      Call MPI_Barrier(cartcomm,ierr)

      call MPI_Finalize(ierr)

      end

C**********************************************************************
C Subroutine param:
C     Define simulation parameters
C**********************************************************************

      subroutine param
      
      use Glob

      implicit none


C Box dimensions in x, y and z
      Lx = 50.0d0
      Ly = 50.0d0
      Lz = 50.0d0

C Number of points in x, y, z and phi, theta
      totnx = 64
      totny = 64
      totnz = 64
      nph = 16
      nth = 16

C Dimensionless stresslet 
      alpha = -1.d0 !alpha = active stress; -ve for pushers and positive for pullers
      beta  = 1.75d0 !beta = flow stress; For more details see equations (22) and (23) of Ezhilan, Saintillan and Shelley (2013)

C     U0 = 67.3*dfr

      U0 = 0.1346d0 ! Steric coefficient
      nu = 0.2d0 ! Volume fraction
      
C Center of mass and orientation diffusivities

      dfx = 2.0d0   ! translational diffusion
      dfr = 0.002d0 ! rotational diffusion

      gamma = 1.d0  ! shape ratio = 1 for slender objects

C Time marching
      dt = 0.005d0
      time = 0.0d0

      pi = 4.d0*atan(1.d0)

      dx = Lx/(totnx)
      dy = Ly/(totny)
      dz = Lz/(totnz)
      dph = 2.d0*pi/real(nph)
      dth = pi/real(nth)

      dy2 = dy

      return
      end


C**********************************************************************
C Subroutine initial:
C     Define initial configuration field
C**********************************************************************

      subroutine initial
 
      use Glob

      implicit none
      include 'mpif.h'


      integer ix,iy,iz,iph,ith,i,j,k,is
      integer nk,ik,idum,ikx,iky,ikz,ind
      character*32 filename1

      real fftheta(nth),ffph(nph)
      real kx,ky,kz,eps,kys
      real ran
      real eps1(40,40,40),eps2(40,40,40)
      real a(8),mean


      allocate(uxwait(nx,ny),uywait(nx,ny),uzwait(nx,ny))
      allocate(EMwait(3,3,nx,ny))


      if (rank.eq.0) then

      print*,'rank 0 allocation started'

      allocate(ft(3,3,totnx,totny,totnz))
      allocate(f3tot(totnx,totny,totnz))
      allocate(f1tot(totnx,totny,totnz),f2tot(totnx,totny,totnz))

      allocate(fxT(nx,ny,nz),fyT(nx,ny,nz))
      allocate(fzT(nx,ny,nz))
      allocate(ftotcompx(2,totnx,totny,totnz))
      allocate(ftotcompy(2,totnx,totny,totnz))
      allocate(ftotcompz(2,totnx,totny,totnz))
      allocate(ukxtot(2,totnx,totny,totnz),ukytot(2,totnx,totny,totnz))
      allocate(ukztot(2,totnx,totny,totnz))
      allocate(ciT(nx,ny,nz))
      allocate(ukxT(2,nx,ny,nz),ukyT(2,nx,ny,nz))
      allocate(ukzT(2,nx,ny,nz))
      allocate(fTcompx(2,nx,ny,nz))
      allocate(fTcompy(2,nx,ny,nz))
      allocate(fTcompz(2,nx,ny,nz))
      allocate(pxiT(nx,ny,nz))
      allocate(pyiT(nx,ny,nz))
      allocate(pziT(nx,ny,nz))
      allocate(uxT(nx,ny,nz))
      allocate(uyT(nx,ny,nz))
      allocate(uzT(nx,ny,nz))

      allocate(LET(nx,ny,nz))
      allocate(pwd(nx,ny,nz))

      print*,'rank 0 allocation done'

      end if

      call MPI_Barrier(cartcomm,ierr)


      allocate(cremesh(nx,ny,nz))
      allocate(cl(nph,nth),clt(nph,nth))


      allocate(LEI(nx,ny,nz))
      allocate(pow(nx,ny,nz))!,sft(nx,ny,nz),sfd(3,nx,ny,nz))
      allocate(cwait(nx,ny,nph,nth))
      allocate(csend(nx,ny,nph,nth))

      allocate(E(3,3,nx,ny,nz),divMM(3,nx,ny,nz))

      allocate(c(nx,ny,nz,nph,nth),ci(nx,ny,nz))

      allocate(fx(nx,ny,nz),fy(nx,ny,nz),fz(nx,ny,nz))
      allocate(ux(nx,ny,nz),uy(nx,ny,nz),uz(nx,ny,nz))
      allocate(phi(nph),theta(nth),ppM(nx,ny,nz))

      allocate(F(nx,ny,nz,nph,nth),Fo(nx,ny,nz,nph,nth))

      allocate(bEW(3,3,nx,ny,nz),M(3,3,nx,ny,nz),ing(3,3,nx,ny,nz))
      allocate(ciphi(nx,ny,nz),MM(3,3,nx,ny,nz),ppE(nx,ny,nz))
      allocate(pxi(nx,ny,nz),pyi(nx,ny,nz),pzi(nx,ny,nz))
      allocate(fcompx(2,nx,ny,nz),fcompy(2,nx,ny,nz),fcompz(2,nx,ny,nz))
      allocate(ukx(2,nx,ny,nz),uky(2,nx,ny,nz),ukz(2,nx,ny,nz))

      allocate(f1(nx,ny,nz,nph,nth))

      call MPI_Barrier(cartcomm,ierr)

      print*,'all ranks allocation done'


      do iph = 1,nph
         phi(iph) = (iph-1.0)/nph*2*pi
      end do
 
      do ith = 1,nth
         theta(ith) = (ith-0.5d0)/nth*pi
      end do

C  Consider nearly uniform/isotropic suspension
      eps = 0.15


C Build random initial distribution
      c = 1.d0       
 
C Number of spatial modes in x and y directions
      nk = 8
C Get amplitudes of the various modes
      idum = -854881
      do ikz = 1,2*nk
         do iky = 1,2*nk
            do ikx = 1,2*nk
               eps1(ikx,iky,ikz) = eps*2.0*(ran(idum)-0.50)
               eps2(ikx,iky,ikz) = eps*2.0*(ran(idum)-0.50)
            enddo
         enddo
      enddo
      do ind = 1,8
         a(ind) = 2.0*(ran(idum)-0.5)
      enddo

!3rd Order Polynomial for fi

      do iph = 1,nph
         ffph(iph) = (a(1)+a(2)*cos(phi(iph))+
     $   a(3)*cos(phi(iph))**2.0+a(4)*cos(phi(iph))**3.0)*
     $   (a(5)+a(6)*sin(phi(iph))
     $   +a(7)*sin(phi(iph))**2.0+a(8)*sin(phi(iph))**3.0)
      enddo

!3rd Order Polynomial for Theta

      do ith = 1,nth
         fftheta(ith) = (a(1)+a(2)*cos(theta(ith))
     $   +a(3)*cos(theta(ith))**2.0
     $   +a(4)*cos(theta(ith))**3.0)*
     $   (a(5)+a(6)*sin(theta(ith))
     $   +a(7)*sin(theta(ith))**2.0
     $   +a(8)*sin(theta(ith))**3.0)
      enddo
 
C For each mode, find angle dependence

      call MPI_CART_COORDS(cartcomm,rank,1,Scoords,ierr)

      Kmin = Scoords(1)*(nz-2) + 1
      Kmax = Scoords(1)*(nz-2) + nz-2



      do ikz = 1,2*nk
         if(rank.eq.0) then
            print *,'wave number = ',ikz
         end if
         do iky = 1,2*nk
            do ikx = 1,2*nk
          
            do ith = 1,nth
               do iph = 1,nph
                  do iz = Kmin,Kmax
                     do iy = 1,ny
                        do ix = 1,nx




              k = iz + 1 - Scoords(1)*(nz-2)
             
!      if (sqrt((ikx-nx)**2.00+(iky-ny)**2.00+
!     $   (ikz-nz)**2.00).gt.7.00) then
      c(ix,iy,k,iph,ith) = c(ix,iy,k,iph,ith)+
     $(eps1(ikx,iky,ikz)*cos(2.0*pi*(ikx-nk)*(ix-1.0)
     $/real(totnx)+2.0*pi*(iky-nk)*(iy-1.0)/real(totny)
     $+2.0*pi*(ikz-nk)*(iz-1.0)/real(totnz))
     $+eps2(ikx,iky,ikz)*sin(2.0*pi*(ikx-nk)*(ix-1.0)
     $/real(totnx)+2.0*pi*(iky-nk)*(iy-1.0)/real(totny)
     $+2.0*pi*(ikz-nk)*(iz-1.0)/real(totnz)))
     $*ffph(iph)*fftheta(ith)
!      endif

c              end if

            enddo
            enddo
            enddo
            enddo
            enddo 
       
           enddo
         enddo
      enddo


C Make sure mean concentration is 1
 
      call intphi(c,ci)

      mean = 0.d0
      do iz = 2,nz-1
         do iy = 1,ny
            do ix = 1,nx
               mean = mean + ci(ix,iy,iz)
            enddo
         enddo
      enddo


      mean = mean/real(nx*ny*(nz-2))

      print*,'mean = ',mean

      
      c = c / mean  

      call intphi(c,ci)

      call MPI_Barrier(cartcomm,ierr)
 

C COPY TO GHOST CELLS FROM NEIGHBORS


C******************************************************************
C	        COORDINATES AND RANKS OF NEIGHBORS
 
         call MPI_CART_COORDS(cartcomm,rank,1,Scoords,ierr)

         coordsR(:) = Scoords(:)
         coordsR(1) = Scoords(1) + 1

         call MPI_CART_RANK(cartcomm, coordsR, rankR, ierr)


         coordsL(:) = Scoords(:)
         coordsL(1) = Scoords(1) - 1

         call MPI_CART_RANK(cartcomm, coordsL, rankL, ierr)


C********************************************************************

         tag = 0
         Call MPI_IRECV(cwait(:,:,:,:),nx*ny*nph*nth,MPI_REAL,rankL,
     $                          tag,cartcomm,req,ierr)
 
         csend(:,:,:,:) = c(:,:,nz-1,:,:) 

         tag = 0
         Call MPI_SEND(csend(:,:,:,:),nx*ny*nph*nth,MPI_REAL,rankR,
     $                          tag,cartcomm,ierr)


         call MPI_WAIT(req,stat,ierr) 

         c(:,:,1,:,:) = cwait(:,:,:,:)

C     	RECEIVE DATA FROM THE LEFT NEIGHBOR 

C*********************************************************************

C********************************************************************

         tag = 1
         Call MPI_IRECV(cwait(:,:,:,:),nx*ny*nph*nth,MPI_REAL,rankR,
     $                          tag,cartcomm,req,ierr)
 
         csend(:,:,:,:) = c(:,:,2,:,:) 


         tag = 1
         Call MPI_SEND(csend(:,:,:,:),nx*ny*nph*nth,MPI_REAL,rankL,
     $                          tag,cartcomm,ierr)

         call MPI_WAIT(req,stat,ierr)

         c(:,:,nz,:,:) = cwait(:,:,:,:)


      call intphi(c,ci)



      fftw_fwd  = -1
      fftw_bkwd =  1
      fftw_est  =  0
      fftw_norm_order = 0


       dy2 = dy
c      dy2 = dy*sqrt(1.d0+(Shear*Stime)**2.d0)

C*********************************************************************************
C*********************************************************************************
C*********************************************************************************




       Call fftw3d_f77_mpi_create_plan(plan,cartcomm,totnx,totny,totnz,
     $                                 fftw_fwd,fftw_est)

       Call fftw3d_f77_mpi_create_plan(iplan,cartcomm,totnx,totny,totnz,
     $                                 fftw_bkwd,fftw_est)

 
       Call fftwnd_f77_mpi_local_sizes(plan,lnx,lxs,lnyt,lyst,lsize)


       do is = 1,numtasks

          if (rank.eq.is-1) then
             print*,'rank = ',rank
             print*,'lxs  = ',lxs
             print*,'lnx  = ',lnx
             print*,'lsize =',lsize
          end if
          call MPI_Barrier(cartcomm,ierr)
       end do

      allocate(fxft(lsize),fyft(lsize),fzft(lsize))
      allocate(uxft(lsize),uyft(lsize),uzft(lsize))
      allocate(work(lsize))

      
      return
      end


C**********************************************************************
C Subroutine intphi:
C     Calculates the integral over the sphere of orientation using trapezoidal rule
C**********************************************************************

      subroutine intphi(cll,clli)

      use Glob

      implicit none

      integer ix,iy,iz,iph,ith

      real cll(nx,ny,nz,nph,nth),clli(nx,ny,nz)
   
      clli = 0.0

      do ith = 1,nth
         ciphi = 0.0
         do iph = 1,nph
            do iz = 2,nz-1
               do iy = 1,ny
                  do ix = 1,nx
                     ciphi(ix,iy,iz) = ciphi(ix,iy,iz)
     $               + cll(ix,iy,iz,iph,ith)*dph
                  enddo
               enddo
            enddo
         enddo
         clli = clli + ciphi*sin(theta(ith))*dth
      enddo
   

      return
      end

C**********************************************************************
C Subroutine stressmagnitude:
C     Get shear stress magnitude in suspension
C**********************************************************************

      subroutine stressmagnitude
    
      use Glob

      implicit none
      include 'mpif.h'


      real TotEntropy,TEntropy,TEntropyT,Entropy,EntropyZ,EntropyY
      real PowerZ,PowerY,Power,TPower
      integer ix,iy,iz,iph,ith,ifile,iorth,ipp,ipt
      integer i,j,k,is !stat(MPI_STATUS_SIZE),is
      character*32 filename1


      real csum


      if (200.0*(itime/200).eq.1.0*itime.OR.itime.eq.1) then
         ifile = itime/200
c1
 
          if (rank.eq.0) then
            print *,'file',ifile
         end if


         call intphi(c,ci)
                  call meanp

C Total Entropy

        f1 = 4.d0*pi*c*(log(4.d0*pi*abs(c)))

        call intphi(f1,LEI)

        Do ith = 1,nth
           Do iph = 1,nph
              csum = 0.d0
              Do iz = 2,nz-1
                 Do iy = 1,ny
                    Do ix = 1,nx

                       csum = csum + c(ix,iy,iz,iph,ith)

                    End Do
                 End Do
              End Do
              cl(iph,ith) = csum/real(nx*ny*(nz-2))
           End Do
        End Do

        Call MPI_REDUCE (cl,clt,nph*nth,MPI_REAL,MPI_SUM,
     $                   rootR,cartcomm,ierr)


C   Total power

        do iz= 1,nz
          do iy = 1,ny
            do ix = 1,nx

	   pow(ix,iy,iz)=0.d0

 	     do i = 1,3
	      do j =1,3

! POW is the term within the integral of equation (75)
	 pow(ix,iy,iz)=pow(ix,iy,iz)+E(i,j,ix,iy,iz)*M(i,j,ix,iy,iz)

	     enddo
	    enddo
	   enddo
	  enddo
	enddo

C*****************************************************************

C     SEND DATA FROM ALL PROCESSORS TO RANK ZERO 

      rootR = 0


      Do is = 1,numtasks-1

         if (rank.eq.is.OR.rank.eq.0) then
c2

         if (rank.eq.rootR) then

            tag   = 0
            Call MPI_IRECV(ciT(:,:,:),nx*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req1,ierr)


            tag   = 9
            Call MPI_IRECV(LET(:,:,:),nx*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req2,ierr)

	    tag   = 10
            Call MPI_IRECV(pwd(:,:,:),nx*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req3,ierr)

         end if


 
         if (rank.eq.is) then

            tag   = 0
            Call MPI_SEND(ci(:,:,:),nx*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)

            tag   = 9
            Call MPI_SEND(LEI(:,:,:),nx*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)
c
            tag   = 10
            Call MPI_SEND(pow(:,:,:),nx*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)
        
	 end if

         if (rank.eq.rootR) then

            call MPI_WAIT(req1,stat,ierr)
            call MPI_WAIT(req2,stat,ierr)
            call MPI_WAIT(req3,stat,ierr)

         end if


C	RANK ZERO RECEIVES DATA FROM ALL PROCESSORS ONE BY ONE


         if (rank.eq.0) then

            call MPI_CART_COORDS(cartcomm,is,1,Scoords,ierr)

C	RANK ZERO COPIES THE RECEIVED DATA INTO A GLOBAL SIZE ARRAY

            Do iz = 2,nz-1
               Do iy = 1,ny
                  Do ix = 1,nx

                     k = Scoords(1)*(nz-2) + iz - 1

                     f1tot(ix,iy,k)    = ciT(ix,iy,iz) 

                     f2tot(ix,iy,k)    = LET(ix,iy,iz)

		     f3tot(ix,iy,k)    = pwd(ix,iy,iz)

                 End Do
               End Do
            End Do

         end if

         end if
c1
         Call MPI_Barrier(cartcomm,ierr)

      end do

      Call MPI_Barrier(cartcomm,ierr)


C	INCLUDING THE DATA FROM RANK ZERO ITSELF

      if (rank.eq.rootR) then
c2
         call MPI_CART_COORDS(cartcomm,rootR,1,Scoords,ierr)

         Do iz = 2,nz-1
            Do iy = 1,ny
               Do ix = 1,nx

                  k = Scoords(1)*(nz-2) + iz - 1

                  f1tot(ix,iy,k)  = ci(ix,iy,iz)

                  f2tot(ix,iy,k)  = LEI(ix,iy,iz)

		  f3tot(ix,iy,k)  = pow(ix,iy,iz)

              End Do
            End Do
         End Do



C       ENTROPY

         TEntropy = 0.d0
	 TPower = 0.d0
         do ix = 1,totnx
            EntropyY = 0.d0
	    PowerY = 0.d0
            do iy = 1,totny
               EntropyZ = 0.d0
	       PowerZ = 0.d0
               do iz = 1,totnz
                  Entropy = f2tot(ix,iy,iz)
		  Power = f3tot(ix,iy,iz)
                  EntropyZ = EntropyZ + Entropy*dz
		  PowerZ = PowerZ + Power*dz/Lz
               end do
               EntropyY = EntropyY + EntropyZ*dy
	       PowerY = PowerY + PowerZ*dy/Ly
            end do
            TEntropy = TEntropy + EntropyY*dx
	    TPower = TPower + PowerY*dx/Lx
         end do

! Spatially and orientationally averaged Entropy and Power
         write(56,*) time, TEntropy,TPower

! Output concentration profile
         write(filename1,'(A3,I4.4,A4)') 'con',ifile,'.txt'

         open(unit=5,file=filename1,status='unknown')
         write(5,*) 'Variables = "x", "y", "z" , "citot"'
         write(5,*) 'Zone I= ',totnx/2, ' J= ', totny/2, ' K= ', totnz/2

         do iz = 1,totnz,2
            do iy = 1,totny,2
               do ix = 1,totnx,2

             write(5,*) (ix-1)*(Lx/totnx),
     $          (Ly/totny)*(iy-1),(Lz/totnz)*(iz-1),f1tot(ix,iy,iz)
               enddo

            enddo
         enddo
         close(5)

! OUTPUT power 
         write(filename1,'(A3,I4.4,A4)') 'POW',ifile,'.txt'

         open(unit=5,file=filename1,status='unknown')
         write(5,*) 'Variables = "x", "y", "z" , "pow"'
         write(5,*) 'Zone I= ',totnx/2, ' J= ', totny/2, ' K= ', totnz/2

         do iz = 1,totnz,2
            do iy = 1,totny,2
               do ix = 1,totnx,2

            write(5,*)  (ix-1)*(Lx/totnx),
     $         (Ly/totny)*(iy-1),(Lz/totnz)*(iz-1), f3tot(ix,iy,iz)
               enddo

            enddo
         enddo

         close(5)

      end if

      call MPI_Barrier(cartcomm,ierr)

      Do is = 1,numtasks-1

         if (rank.eq.is.OR.rank.eq.0) then

         if (rank.eq.rootR) then

            tag   = 6
            Call MPI_IRECV(uxT(:,:,:),nx*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req1,ierr)

            tag   = 7
            Call MPI_IRECV(uyT(:,:,:),nx*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req2,ierr)

            tag   = 8
            Call MPI_IRECV(uzT(:,:,:),nx*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req3,ierr)

         end if


         if (rank.eq.is) then


            tag   = 6
            Call MPI_SEND(ux(:,:,:),nx*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)

            tag   = 7
            Call MPI_SEND(uy(:,:,:),nx*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)

            tag   = 8
            Call MPI_SEND(uz(:,:,:),nx*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)

c
         end if

         if (rank.eq.rootR) then

            call MPI_WAIT(req1,stat,ierr)
            call MPI_WAIT(req2,stat,ierr)
            call MPI_WAIT(req3,stat,ierr)

         end if

C	RANK ZERO RECEIVES DATA FROM ALL PROCESSORS ONE BY ONE


         if (rank.eq.0) then

            call MPI_CART_COORDS(cartcomm,is,1,Scoords,ierr)

C	RANK ZERO COPIES THE RECEIVED DATA INTO A GLOBAL SIZE ARRAY

            Do iz = 2,nz-1
               Do iy = 1,ny
                  Do ix = 1,nx

                     k = Scoords(1)*(nz-2) + iz - 1

                     f1tot(ix,iy,k)   = uxT(ix,iy,iz)
                     f2tot(ix,iy,k)   = uyT(ix,iy,iz)
                     f3tot(ix,iy,k)   = uzT(ix,iy,iz)

                 End Do
               End Do
            End Do

         end if

         end if
c1
         Call MPI_Barrier(cartcomm,ierr)

      end do

         Call MPI_Barrier(cartcomm,ierr)


C	INCLUDING THE DATA FROM RANK ZERO ITSELF

      if (rank.eq.rootR) then

         call MPI_CART_COORDS(cartcomm,rootR,1,Scoords,ierr)

         Do iz = 2,nz-1
            Do iy = 1,ny
               Do ix = 1,nx

                  k = Scoords(1)*(nz-2) + iz - 1

                  f1tot(ix,iy,k)   = ux(ix,iy,iz)
                  f2tot(ix,iy,k)   = uy(ix,iy,iz)
                  f3tot(ix,iy,k)   = uz(ix,iy,iz)

 
              End Do
            End Do
         End Do

! Output velocity
         write(filename1,'(A3,I4.4,A4)') 'vel',ifile,'.txt'

         open(unit=5,file=filename1,status='unknown')
          write(5,*) 'Variables = "x", "y", "z" , "ux", "uy", "uz"'
         write(5,*) 'Zone I= ',totnx/2, ' J= ', totny/2, ' K= ', totnz/2

         do iz = 1,totnz,2
            do iy = 1,totny,2
               do ix = 1,totnx,2

            write(5,*)  (ix-1)*(Lx/totnx),
     $          (Ly/totny)*(iy-1),(Lz/totnz)*(iz-1),f1tot(ix,iy,iz)
     $          ,f2tot(ix,iy,iz),f3tot(ix,iy,iz)
               enddo

            enddo
         enddo
         close(5)

      end if

      call MPI_Barrier(cartcomm,ierr)

      Do is = 1,numtasks-1

         if (rank.eq.is.OR.rank.eq.0) then
 
         if (rank.eq.rootR) then

            tag   = 3
            Call MPI_IRECV(pxiT(:,:,:),nx*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req1,ierr)

            tag   = 4
            Call MPI_IRECV(pyiT(:,:,:),nx*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req2,ierr)

            tag   = 5
            Call MPI_IRECV(pziT(:,:,:),nx*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req3,ierr)

         end if


         if (rank.eq.is) then

            tag   = 3
            Call MPI_SEND(pxi(:,:,:),nx*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)

            tag   = 4
            Call MPI_SEND(pyi(:,:,:),nx*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)

            tag   = 5
            Call MPI_SEND(pzi(:,:,:),nx*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)

         end if

         if (rank.eq.rootR) then

            call MPI_WAIT(req1,stat,ierr)
            call MPI_WAIT(req2,stat,ierr)
            call MPI_WAIT(req3,stat,ierr)
 
         end if

C	RANK ZERO RECEIVES DATA FROM ALL PROCESSORS ONE BY ONE


         if (rank.eq.0) then

            call MPI_CART_COORDS(cartcomm,is,1,Scoords,ierr)

C	RANK ZERO COPIES THE RECEIVED DATA INTO A GLOBAL SIZE ARRAY

            Do iz = 2,nz-1
               Do iy = 1,ny
                  Do ix = 1,nx

                     k = Scoords(1)*(nz-2) + iz - 1

                     f1tot(ix,iy,k)  = pxiT(ix,iy,iz)
                     f2tot(ix,iy,k)  = pyiT(ix,iy,iz)
                     f3tot(ix,iy,k)  = pziT(ix,iy,iz)

                 End Do
               End Do
            End Do

         end if

         end if
c1
         Call MPI_Barrier(cartcomm,ierr)

      end do

      Call MPI_Barrier(cartcomm,ierr)


C	INCLUDING THE DATA FROM RANK ZERO ITSELF

      if (rank.eq.rootR) then
c2
         call MPI_CART_COORDS(cartcomm,rootR,1,Scoords,ierr)

        Do iz = 2,iz-1
            Do iy = 1,ny
               Do ix = 1,nx

                  k = Scoords(1)*(nz-2) + iz - 1

                  f1tot(ix,iy,k)  = pxi(ix,iy,iz)
                  f2tot(ix,iy,k)  = pyi(ix,iy,iz)
                  f3tot(ix,iy,k)  = pzi(ix,iy,iz)
 
              End Do
            End Do
         End Do

! MNP = first orientational moment / polarization

         write(filename1,'(A3,I4.4,A4)') 'mnp',ifile,'.txt'

         open(unit=5,file=filename1,status='unknown')
         write(5,*) 'Variables = "x", "y", "z" , "px", "py", "pz"'
         write(5,*) 'Zone I= ',totnx/2, ' J= ', totny/2, ' K= ', totnz/2

         do iz = 1,totnz,2
            do iy = 1,totny,2
               do ix = 1,totnx,2

            write(5,*)  (ix-1)*(Lx/totnx),
     $          (Ly/totny)*(iy-1),(Lz/totnz)*(iz-1),f1tot(ix,iy,iz),
     $          f2tot(ix,iy,iz),f3tot(ix,iy,iz)
               enddo

            enddo
         enddo
         close(5)

! LCON = spatially averaged orientational distribution

         write(filename1,'(A4,I4.4,A4)') 'lcon',ifile,'.txt'

         open(unit=5,file=filename1,status='unknown')
         write(5,*) 'Variables = "phi", "theta", "lcon" '
         write(5,*) 'Zone I = ',nph,' J= ', nth

         do ith = 1,nth
            do iph = 1,nph

               write(5,*) phi(iph),theta(ith),clt(iph,ith)/numtasks

            enddo
         enddo
         
         close(5)

        end if
c1
        call MPI_Barrier(cartcomm,ierr)

        rootR = 0

 	Do is = 1,numtasks-1

         if (rank.eq.is.OR.rank.eq.0) then
c2

         if (rank.eq.rootR) then

            tag   = 0
            Call MPI_IRECV(ing(:,:,:,:,:),3*3*nx*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req1,ierr)

         end if

         if (rank.eq.is) then

            tag   = 0
            Call MPI_SEND(M(:,:,:,:,:),3*3*nx*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)
c	    print *,Scoords(1),M(:,:,32,32,5)
        endif

         if (rank.eq.rootR) then

            call MPI_WAIT(req1,stat,ierr)

        endif
	
	         if (rank.eq.0) then

            call MPI_CART_COORDS(cartcomm,is,1,Scoords,ierr)
c	    print *,Scoords(1),ing(:,:,32,32,5)
C       RANK ZERO COPIES THE RECEIVED DATA INTO A GLOBAL SIZE ARRAY

            Do iz = 2,nz-1
               Do iy = 1,ny
                  Do ix = 1,nx

                     k = Scoords(1)*(nz-2) + iz - 1

                     ft(:,:,ix,iy,k)    = ing(:,:,ix,iy,iz)

                 End Do
               End Do
            End Do

         end if

         end if
c1
         Call MPI_Barrier(cartcomm,ierr)

      end do

      Call MPI_Barrier(cartcomm,ierr)

      if (rank.eq.rootR) then
c2
         call MPI_CART_COORDS(cartcomm,rootR,1,Scoords,ierr)

         Do iz = 2,nz-1
            Do iy = 1,ny
               Do ix = 1,nx

                  k = Scoords(1)*(nz-2) + iz - 1

                  ft(:,:,ix,iy,k)  = M(:,:,ix,iy,iz)

              End Do
            End Do
         End Do

         write(filename1,'(A3,I4.4,A4)') 'str',ifile,'.txt'

         open(unit=5,file=filename1,status='unknown')
         write(5,*) 'Variables = "x", "y", "z" , "11", "22", "33", 
     $    "21", "31", "32"'
         write(5,*) 'Zone I= ',totnx/2, ' J= ', totny/2, ' K= ', totnz/2

         do iz = 1,totnz,2
            do iy = 1,totny,2
               do ix = 1,totnx,2

             write(5,*) (ix-1)*(Lx/totnx),
     $          (Ly/totny)*(iy-1),(Lz/totnz)*(iz-1),ft(1,1,ix,iy,iz),
     $          ft(2,2,ix,iy,iz),ft(3,3,ix,iy,iz),ft(1,2,ix,iy,iz),
     $          ft(1,3,ix,iy,iz),ft(2,3,ix,iy,iz)
               enddo
	     enddo
	 enddo

	close(5)

	endif

       call MPI_Barrier(cartcomm,ierr)

      end if

      return
      end

C**********************************************************************
C Subroutine meanp:
C     Get polarization field
C**********************************************************************

      subroutine meanp

      use Glob

      implicit none

      integer ix,iy,iz,iph,ith
      real p(3)

C Get orientation vector

      pxi = 0.d0
      pyi = 0.d0
      pzi = 0.d0

      do ith = 1,nth
         do iph = 1,nph

            p(1) = sin(theta(ith))*cos(phi(iph))
            p(2) = sin(theta(ith))*sin(phi(iph))
            p(3) = cos(theta(ith))

            do iz = 2,nz-1
               do iy = 1,ny
                  do ix = 1,nx

                     pxi(ix,iy,iz) = pxi(ix,iy,iz)
     $               + c(ix,iy,iz,iph,ith)
     $                  *p(1)
     $                  *sin(theta(ith))*dph*dth

                     pyi(ix,iy,iz) = pyi(ix,iy,iz)
     $               + c(ix,iy,iz,iph,ith)
     $                  *p(2)
     $                  *sin(theta(ith))*dph*dth

                     pzi(ix,iy,iz) = pzi(ix,iy,iz)
     $               + c(ix,iy,iz,iph,ith)
     $                  *p(3)
     $                  *sin(theta(ith))*dph*dth

                  enddo
               enddo
            enddo
         enddo
      enddo

C Integrate over phi to obtain stress field
      call intphi(c,ci)


      pxi = pxi / ci
      pyi = pyi / ci
      pzi = pzi / ci


      return
      end


C**********************************************************************
C       Calculate the second moment
C       M = \int \psi p p dp, 
C       Note that this is a crucial term that appears both the extra stress term (22) 
C       and steric interaction term in equation (20) in Ezhilan et al. (2013)
C**********************************************************************
      
      subroutine secondmoment

      use Glob

      implicit none

      integer ix,iy,iz,iph,ith
      real    p(3)
      include 'mpif.h'


C Get orientation vector
      do ith = 1,nth
         do iph = 1,nph

            p(1) = sin(theta(ith))*cos(phi(iph))
            p(2) = sin(theta(ith))*sin(phi(iph))
            p(3) = cos(theta(ith))

            do iz = 1,nz
               do iy = 1,ny
                  do ix = 1,nx


            f1(ix,iy,iz,iph,ith) =
     $      c(ix,iy,iz,iph,ith)*(p(1)*p(1))



                  enddo
               enddo
            enddo
         enddo
      enddo


      call intphi(f1,M(1,1,:,:,:))
      
      do ith = 1,nth
         do iph = 1,nph

            p(1) = sin(theta(ith))*cos(phi(iph))
            p(2) = sin(theta(ith))*sin(phi(iph))
            p(3) = cos(theta(ith))

            do iz = 1,nz
               do iy = 1,ny
                  do ix = 1,nx


            f1(ix,iy,iz,iph,ith) =
     $      c(ix,iy,iz,iph,ith)*(p(1)*p(2))


                  enddo
               enddo
            enddo
         enddo
      enddo


      call intphi(f1,M(1,2,:,:,:))
      M(2,1,:,:,:)=M(1,2,:,:,:)

      do ith = 1,nth
         do iph = 1,nph

            p(1) = sin(theta(ith))*cos(phi(iph))
            p(2) = sin(theta(ith))*sin(phi(iph))
            p(3) = cos(theta(ith))

            do iz = 1,nz
               do iy = 1,ny
                  do ix = 1,nx


            f1(ix,iy,iz,iph,ith) =
     $      c(ix,iy,iz,iph,ith)*(p(1)*p(3))


                  enddo
               enddo
            enddo
         enddo
      enddo


      call intphi(f1,M(1,3,:,:,:))
      M(3,1,:,:,:)=M(1,3,:,:,:)

      do ith = 1,nth
         do iph = 1,nph

            p(1) = sin(theta(ith))*cos(phi(iph))
            p(2) = sin(theta(ith))*sin(phi(iph))
            p(3) = cos(theta(ith))

            do iz = 1,nz
               do iy = 1,ny
                  do ix = 1,nx


            f1(ix,iy,iz,iph,ith) =
     $      c(ix,iy,iz,iph,ith)*(p(2)*p(2))

                  enddo
               enddo
            enddo
         enddo
      enddo


C Integrate over phi to obtain stress field
      call intphi(f1,M(2,2,:,:,:))


      do ith = 1,nth
         do iph = 1,nph

            p(1) = sin(theta(ith))*cos(phi(iph))
            p(2) = sin(theta(ith))*sin(phi(iph))
            p(3) = cos(theta(ith))

            do iz = 1,nz
               do iy = 1,ny
                  do ix = 1,nx



            f1(ix,iy,iz,iph,ith) =
     $      c(ix,iy,iz,iph,ith)*(p(2)*p(3))



                  enddo
               enddo
            enddo
         enddo
      enddo


C Integrate over phi to obtain stress field
      call intphi(f1,M(2,3,:,:,:))
      M(3,2,:,:,:) = M(2,3,:,:,:)

      do ith = 1,nth
         do iph = 1,nph

            p(1) = sin(theta(ith))*cos(phi(iph))
            p(2) = sin(theta(ith))*sin(phi(iph))
            p(3) = cos(theta(ith))

            do iz = 1,nz
               do iy = 1,ny
                  do ix = 1,nx


            f1(ix,iy,iz,iph,ith) =
     $      c(ix,iy,iz,iph,ith)*(p(3)*p(3))

                 enddo
               enddo
            enddo
         enddo
      enddo


C Integrate over phi to obtain stress field
      call intphi(f1,M(3,3,:,:,:))

C COPY TO GHOST CELLS FROM NEIGHBORS

!         call MPI_CART_COORDS(cartcomm,rank,1,Scoords,ierr)

C	SEND DATA TO THE RIGHT NEIGHBOR 
!         coordsR(:) = Scoords(:)
!         coordsR(1) = Scoords(1) + 1

!         call MPI_CART_RANK(cartcomm, coordsR, rankR, ierr)

!         coordsL(:) = Scoords(:)
!         coordsL(1) = Scoords(1) - 1

!         call MPI_CART_RANK(cartcomm, coordsL, rankL, ierr)

C     	RECEIVE DATA FROM THE LEFT NEIGHBOR 

!         tag = 0
!         Call MPI_IRECV(EMwait(:,:,:,:),3*3*nx*ny,MPI_REAL,rankL,
!     $                          tag,cartcomm,req1,ierr)

!         tag = 0
!         Call MPI_SEND(M(:,:,:,:,nz-1),3*3*nx*ny,MPI_REAL,rankR,
!     $                          tag,cartcomm,ierr)
 
!         call MPI_WAIT(req1,stat,ierr)

!         M(:,:,:,:,1) = EMwait(:,:,:,:)


C	RECEIVE DATA FROM THE RIGHT NEIGHBOR

!         tag = 3
!         Call MPI_IRECV(EMwait(:,:,:,:),3*3*nx*ny,MPI_REAL,rankR,
!     $                          tag,cartcomm,req1,ierr)

C	SEND DATA TO THE LEFT NEIGHBOR
!         tag = 3
!         Call MPI_SEND(M(:,:,:,:,2),3*3*nx*ny,MPI_REAL,rankL,
!     $                          tag,cartcomm,ierr)

!         call MPI_WAIT(req1,stat,ierr)

!         M(:,:,:,:,nz) = EMwait(:,:,:,:)

      return
      end


C**********************************************************************
C Subroutine stressfield:
C     Obtain stressfield (RHS in Stokes equations)
C**********************************************************************

      subroutine stressfield
 
      use Glob

      implicit none

      integer ix,iy,iz,iph,ith,i,j,k
      integer ixp,ixm,iyp,iym,izp,izm
C      integer stat(MPI_STATUS_SIZE)
      real p(3),del1(3),del2(3),del3(3),dot1,dot2,dot3

      dy2 = dy
c      dy2 = dy*sqrt(1.d0+(Shear*Stime)**2.d0)

		do iz = 1,nz
                  do iy = 1,ny
                    do ix = 1,nx

			do j = 1,3
			   do i = 1,3
			
			      MM(i,j,ix,iy,iz) = 0.d0     

			    do k = 1,3
 
			   MM(i,j,ix,iy,iz) = MM(i,j,ix,iy,iz) + 
     $			   M(i,k,ix,iy,iz)*M(k,j,ix,iy,iz)

			    enddo

		 	   enddo
		        enddo

		    enddo
		  enddo
		enddo

     
                do iz = 2,nz-1

		    izm = iz-1
	            izp = iz+1

                  do iy = 1,ny

                    iym = iy-1
                    if (iym.lt.1)  iym = ny
                    iyp = iy+1
                    if (iyp.gt.ny) iyp = 1

                    do ix = 1,nx
                     ixm = ix-1
                     if (ixm.lt.1) ixm = nx
                     ixp = ix+1
                     if (ixp.gt.nx) ixp = 1

                           do i = 1,3

			divMM(i,ix,iy,iz) = 
     $                     (MM(i,1,ixp,iy,iz)-MM(i,1,ixm,iy,iz))/(2*dx)
     $                    +(MM(i,2,ix,iyp,iz)-MM(i,2,ix,iym,iz))/(2*dy2)
     $                    +(MM(i,3,ix,iy,izp)-MM(i,3,ix,iy,izm))/(2*dz)

			   enddo

                    enddo
                  enddo
                enddo

		fx = 0.d0
		fy = 0.d0
		fz = 0.d0

      do ith = 1,nth
         do iph = 1,nph

            p(1) = sin(theta(ith))*cos(phi(iph))
            p(2) = sin(theta(ith))*sin(phi(iph))
            p(3) = cos(theta(ith))

              do iz = 1,nz
                 do iy = 1,ny
                    do ix = 1,nx

                        ppE(ix,iy,iz) = 0.d0

                        do j = 1,3
                           do i = 1,3

              ppE(ix,iy,iz) = ppE(ix,iy,iz) + p(i)*p(j)*E(i,j,ix,iy,iz)

                           enddo
                        enddo
                    enddo
                enddo
              enddo

              do iz = 1,nz
                 do iy = 1,ny
                    do ix = 1,nx

                        ppM(ix,iy,iz) = 0.d0

                        do j = 1,3
                           do i = 1,3

              ppM(ix,iy,iz) = ppM(ix,iy,iz) + p(i)*p(j)*M(i,j,ix,iy,iz)

                           enddo
                        enddo
                    enddo
                 enddo
              enddo


            do iz = 2,nz-1

               izm = iz-1
               izp = iz+1

               do iy = 1,ny

                  iym = iy-1
                  if (iym.lt.1)  iym = ny
                  iyp = iy+1
                  if (iyp.gt.ny) iyp = 1

                  do ix = 1,nx
                     ixm = ix-1
                     if (ixm.lt.1) ixm = nx
                     ixp = ix+1
                     if (ixp.gt.nx) ixp = 1

C Del psi
               del1(1) = (c(ixp,iy,iz,iph,ith)
     $                   -c(ixm,iy,iz,iph,ith))/2.d0/dx
            
               del1(2) = (c(ix,iyp,iz,iph,ith)
     $                   -c(ix,iym,iz,iph,ith))/2.d0/dy2

               del1(3) = (c(ix,iy,izp,iph,ith)
     $                   -c(ix,iy,izm,iph,ith))/2.d0/dz


C Del ppE*psi

	       del2(1) = (c(ixp,iy,iz,iph,ith)*ppE(ixp,iy,iz)
     $                   -c(ixm,iy,iz,iph,ith)*ppE(ixm,iy,iz))/2.d0/dx

               del2(2) = (c(ix,iyp,iz,iph,ith)*ppE(ix,iyp,iz)
     $                   -c(ix,iym,iz,iph,ith)*ppE(ix,iym,iz))/2.d0/dy2
               del2(3) = (c(ix,iy,izp,iph,ith)*ppE(ix,iy,izp)
     $                   -c(ix,iy,izm,iph,ith)*ppE(ix,iy,izm))/2.d0/dz

C Del ppM*psi

               del3(1) = (c(ixp,iy,iz,iph,ith)*ppM(ixp,iy,iz)
     $                   -c(ixm,iy,iz,iph,ith)*ppM(ixm,iy,iz))/2.d0/dx

               del3(2) = (c(ix,iyp,iz,iph,ith)*ppM(ix,iyp,iz)
     $                   -c(ix,iym,iz,iph,ith)*ppM(ix,iym,iz))/2.d0/dy2

               del3(3) = (c(ix,iy,izp,iph,ith)*ppM(ix,iy,izp)
     $                   -c(ix,iy,izm,iph,ith)*ppM(ix,iy,izm))/2.d0/dz


C Do matrix vector multiply
C Remember that D is the trace-free version of M
C Recast equation (22) in terms of M, and then the calculations below
C are straight forward to follow.
               dot1 = p(1)*del1(1)+p(2)*del1(2)+p(3)*del1(3)
               dot2 = p(1)*del2(1)+p(2)*del2(2)+p(3)*del2(3)
	           dot3 = p(1)*del3(1)+p(2)*del3(2)+p(3)*del3(3)

                    fx(ix,iy,iz) = fx(ix,iy,iz) +
     $		     ((dot1*p(1)-del1(1)/3.d0)*alpha
     $             + (dot2*p(1)-del2(1)/3.d0)*beta*nu
     $             + (dot3*p(1))*2*U0*beta*nu)
     $                  *sin(theta(ith))*dph*dth

                    fy(ix,iy,iz) = fy(ix,iy,iz) +
     $               ((dot1*p(2)-del1(2)/3.d0)*alpha
     $             + (dot2*p(2)-del2(2)/3.d0)*beta*nu
     $             + (dot3*p(2))*2*U0*beta*nu)
     $                  *sin(theta(ith))*dph*dth

                    fz(ix,iy,iz) = fz(ix,iy,iz) +
     $               ((dot1*p(3)-del1(3)/3.d0)*alpha
     $             + (dot2*p(3)-del2(3)/3.d0)*beta*nu
     $             + (dot3*p(3))*2*U0*beta*nu)
     $                  *sin(theta(ith))*dph*dth
 
               enddo
             enddo
           enddo

         enddo
      enddo

	do iz =2,nz-1
	   do iy = 1,ny
	      do ix = 1,nx

	fx(ix,iy,iz) = fx(ix,iy,iz) - 2*U0*beta*nu*divMM(1,ix,iy,iz)
	fy(ix,iy,iz) = fy(ix,iy,iz) - 2*U0*beta*nu*divMM(2,ix,iy,iz)
	fz(ix,iy,iz) = fz(ix,iy,iz) - 2*U0*beta*nu*divMM(3,ix,iy,iz)

	      enddo
	      enddo
          enddo


      return
      end


C**********************************************************************
C Subroutine velfield:
C     Obtain velocity field spectrally (Hasimoto solution)
C**********************************************************************

      subroutine velfield
      
      use Glob

      implicit none
      include 'mpif.h'



      integer ix,iy,iz,nn(3),is
      integer ixp,ixm,iyp,iym,izp,izm,i,j,k
C      integer stat(MPI_STATUS_SIZE)


      real kx,ky,kz,ksq,fk(2,3),dot(2)
      real cst,kys
      real W(3,3),dU(3,3)
      complex img

      img = (0,1) 

       dy2 = dy
c      dy2 = dy*sqrt(1.d0+(Shear*Stime)**2.d0)

C*********************************************************************************
C*********************************************************************************
C*********************************************************************************

       
       do iz = 2,nz-1
          do iy = 1,ny
             do ix = 1,nx

                fxft(((iz-2)*ny+(iy-1))*nx+ix) = fx(ix,iy,iz)
                fyft(((iz-2)*ny+(iy-1))*nx+ix) = fy(ix,iy,iz)
                fzft(((iz-2)*ny+(iy-1))*nx+ix) = fz(ix,iy,iz)
         
             end do
          end do
       end do
       

       work = 0
       Call fftwnd_f77_mpi(plan,1,fxft,work,0,fftw_norm_order)
       work = 0
       Call fftwnd_f77_mpi(plan,1,fyft,work,0,fftw_norm_order)
       work = 0
       Call fftwnd_f77_mpi(plan,1,fzft,work,0,fftw_norm_order)
 
       
 
       do iz = 2,nz-1
          do iy = 1,ny
             do ix = 1,nx

              fcompx(1,ix,iy,iz) = real(fxft(((iz-2)*ny+(iy-1))*nx+ix)) 
              fcompy(1,ix,iy,iz) = real(fyft(((iz-2)*ny+(iy-1))*nx+ix))
              fcompz(1,ix,iy,iz) = real(fzft(((iz-2)*ny+(iy-1))*nx+ix))

              fcompx(2,ix,iy,iz) = aimag(fxft(((iz-2)*ny+(iy-1))*nx+ix))
              fcompy(2,ix,iy,iz) = aimag(fyft(((iz-2)*ny+(iy-1))*nx+ix))
              fcompz(2,ix,iy,iz) = aimag(fzft(((iz-2)*ny+(iy-1))*nx+ix))
        
             end do
          end do
       end do
 
       


C Get Fourier coefficients of the velocity field (Hasimoto)

      call MPI_CART_COORDS(cartcomm,rank,1,Scoords,ierr)

      Kmin = Scoords(1)*(nz-2) + 1
      Kmax = Scoords(1)*(nz-2) + nz-2



      do iz = Kmin,Kmax
         kz = iz-1
         if (iz.gt.(totnz/2)) kz = iz-1-totnz
         kz = kz/Lz

         do iy = 1,ny
            ky = iy-1
            if (iy.gt.(totny/2)) ky = iy-1-totny
            ky = ky/Ly

            do ix = 1,nx
               kx = ix-1
               if (ix.gt.(totnx/2)) kx = ix-1-totnx
               kx = kx/Lx
           
              
              kys = ky

              i = ix
              j = iy
              k = iz + 1 - Scoords(1)*(nz-2)          
 
              ksq = kx*kx+kys*kys+kz*kz

              fk(1,1) = fcompx(1,i,j,k)
              fk(1,2) = fcompy(1,i,j,k)
              fk(1,3) = fcompz(1,i,j,k)
              fk(2,1) = fcompx(2,i,j,k)
              fk(2,2) = fcompy(2,i,j,k)
              fk(2,3) = fcompz(2,i,j,k)
                 
              dot(1) = fk(1,1)*kx+fk(1,2)*kys+fk(1,3)*kz
              dot(2) = fk(2,1)*kx+fk(2,2)*kys+fk(2,3)*kz
                  
              ukx(1,i,j,k) = (fk(1,1)-dot(1)*kx/ksq)/ksq/(4*pi**2.00)
              ukx(2,i,j,k) = (fk(2,1)-dot(2)*kx/ksq)/ksq/(4*pi**2.00)
              uky(1,i,j,k) = (fk(1,2)-dot(1)*kys/ksq)/ksq/(4*pi**2.00)
              uky(2,i,j,k) = (fk(2,2)-dot(2)*kys/ksq)/ksq/(4*pi**2.00)
              ukz(1,i,j,k) = (fk(1,3)-dot(1)*kz/ksq)/ksq/(4*pi**2.00)
              ukz(2,i,j,k) = (fk(2,3)-dot(2)*kz/ksq)/ksq/(4*pi**2.00)
               


            enddo
         enddo
      enddo


      ix = 1
      iy = 1
      iz = 1

      if (iz.ge.Kmin.AND.iz.le.Kmax) then
 
          k = iz + 1 - Scoords(1)*(nz-2)
               
          ukx(1,ix,iy,k) = 0.00
          ukx(2,ix,iy,k) = 0.00
          uky(1,ix,iy,k) = 0.00
          uky(2,ix,iy,k) = 0.00
          ukz(1,ix,iy,k) = 0.00
          ukz(2,ix,iy,k) = 0.00

       end if
      
       Call MPI_Barrier(cartcomm,ierr)


       do iz = 2,nz-1
          do iy = 1,ny
             do ix = 1,nx

                uxft(((iz-2)*ny+(iy-1))*nx+ix) =
     $          (ukx(1,ix,iy,iz)+img*ukx(2,ix,iy,iz))

                uyft(((iz-2)*ny+(iy-1))*nx+ix) =
     $          (uky(1,ix,iy,iz)+img*uky(2,ix,iy,iz))

                uzft(((iz-2)*ny+(iy-1))*nx+ix) =
     $          (ukz(1,ix,iy,iz)+img*ukz(2,ix,iy,iz))

        
             end do
          end do
       end do

       work = 0

       Call fftwnd_f77_mpi(iplan,1,uxft,work,0,fftw_norm_order)
       work = 0

       Call fftwnd_f77_mpi(iplan,1,uyft,work,0,fftw_norm_order)
       work = 0

       Call fftwnd_f77_mpi(iplan,1,uzft,work,0,fftw_norm_order)
       work = 0


       cst = real(totnx*totny*totnz)

       do iz = 2,nz-1
          do iy = 1,ny
             do ix = 1,nx

              ux(ix,iy,iz) = real(uxft(((iz-2)*ny+(iy-1))*nx+ix))/cst 
              uy(ix,iy,iz) = real(uyft(((iz-2)*ny+(iy-1))*nx+ix))/cst
              uz(ix,iy,iz) = real(uzft(((iz-2)*ny+(iy-1))*nx+ix))/cst
       
             end do
          end do
       end do
 
C COPY TO GHOST CELLS FROM NEIGHBORS

         call MPI_CART_COORDS(cartcomm,rank,1,Scoords,ierr)

C	SEND DATA TO THE RIGHT NEIGHBOR 
         coordsR(:) = Scoords(:)
         coordsR(1) = Scoords(1) + 1

         call MPI_CART_RANK(cartcomm, coordsR, rankR, ierr)

         coordsL(:) = Scoords(:)
         coordsL(1) = Scoords(1) - 1

         call MPI_CART_RANK(cartcomm, coordsL, rankL, ierr)

        
C     	RECEIVE DATA FROM THE LEFT NEIGHBOR 

         tag = 0
         Call MPI_IRECV(uxwait(:,:),nx*ny,MPI_REAL,rankL,
     $                          tag,cartcomm,req1,ierr)
         tag = 1
         Call MPI_IRECV(uywait(:,:),nx*ny,MPI_REAL,rankL,
     $                          tag,cartcomm,req2,ierr)
         tag = 2
         Call MPI_IRECV(uzwait(:,:),nx*ny,MPI_REAL,rankL,
     $                          tag,cartcomm,req3,ierr)



         tag = 0
         Call MPI_SEND(ux(:,:,nz-1),nx*ny,MPI_REAL,rankR,
     $                          tag,cartcomm,ierr)
         tag = 1
         Call MPI_SEND(uy(:,:,nz-1),nx*ny,MPI_REAL,rankR,
     $                          tag,cartcomm,ierr)
         tag = 2
         Call MPI_SEND(uz(:,:,nz-1),nx*ny,MPI_REAL,rankR,
     $                          tag,cartcomm,ierr)

 
         call MPI_WAIT(req1,stat,ierr)
         call MPI_WAIT(req2,stat,ierr)
         call MPI_WAIT(req3,stat,ierr)

         ux(:,:,1) = uxwait(:,:)
         uy(:,:,1) = uywait(:,:)
         uz(:,:,1) = uzwait(:,:)


C	RECEIVE DATA FROM THE RIGHT NEIGHBOR

         tag = 3
         Call MPI_IRECV(uxwait(:,:),nx*ny,MPI_REAL,rankR,
     $                          tag,cartcomm,req1,ierr)
         tag = 4
         Call MPI_IRECV(uywait(:,:),nx*ny,MPI_REAL,rankR,
     $                          tag,cartcomm,req2,ierr)
         tag = 5
         Call MPI_IRECV(uzwait(:,:),nx*ny,MPI_REAL,rankR,
     $                          tag,cartcomm,req3,ierr)



C	SEND DATA TO THE LEFT NEIGHBOR
         tag = 3
         Call MPI_SEND(ux(:,:,2),nx*ny,MPI_REAL,rankL,
     $                          tag,cartcomm,ierr)
         tag = 4
         Call MPI_SEND(uy(:,:,2),nx*ny,MPI_REAL,rankL,
     $                          tag,cartcomm,ierr)
         tag = 5
         Call MPI_SEND(uz(:,:,2),nx*ny,MPI_REAL,rankL,
     $                          tag,cartcomm,ierr)


         call MPI_WAIT(req1,stat,ierr)
         call MPI_WAIT(req2,stat,ierr)
         call MPI_WAIT(req3,stat,ierr)

         ux(:,:,nz) = uxwait(:,:)
         uy(:,:,nz) = uywait(:,:)
         uz(:,:,nz) = uzwait(:,:)



C	RECEIVING DATA FROM NEIGHBORS AND STORING IT IN THE GHOST CELLS


C Get tensor cbEW for angle dynamics

      do iz = 2,nz-1

         izm = iz-1
         izp = iz+1

         do iy = 1,ny

            iym = iy-1
            if (iym.lt.1)  iym = ny
            iyp = iy+1
            if (iyp.gt.ny) iyp = 1

            do ix = 1,nx
               ixm = ix-1
               if (ixm.lt.1)  ixm = nx
               ixp = ix+1
               if (ixp.gt.nx) ixp = 1

                  dU(1,1) = (ux(ixp,iy,iz)-ux(ixm,iy,iz))/2.00/dx

                  dU(1,2) = (ux(ix,iyp,iz)-ux(ix,iym,iz))/2.00/dy2

                  dU(1,3) = (ux(ix,iy,izp)-ux(ix,iy,izm))/2.00/dz


                  dU(2,1) = (uy(ixp,iy,iz)-uy(ixm,iy,iz))/2.00/dx

                  dU(2,2) = (uy(ix,iyp,iz)-uy(ix,iym,iz))/2.00/dy2

                  dU(2,3) = (uy(ix,iy,izp)-uy(ix,iy,izm))/2.00/dz


                  dU(3,1) = (uz(ixp,iy,iz)-uz(ixm,iy,iz))/2.00/dx

                  dU(3,2) = (uz(ix,iyp,iz)-uz(ix,iym,iz))/2.00/dy2


                  dU(3,3) = (uz(ix,iy,izp)-uz(ix,iy,izm))/2.00/dz

 
            do j = 1,3
            do i = 1,3
               E(i,j,ix,iy,iz) = (dU(i,j)+dU(j,i))/2.00
               W(i,j) = (dU(i,j)-dU(j,i))/2.00
               bEW(i,j,ix,iy,iz) = gamma*E(i,j,ix,iy,iz)+W(i,j)
            enddo
            enddo

            enddo
         enddo
      enddo


C COPY TO GHOST CELLS FROM NEIGHBORS

         call MPI_CART_COORDS(cartcomm,rank,1,Scoords,ierr)

C	SEND DATA TO THE RIGHT NEIGHBOR 
         coordsR(:) = Scoords(:)
         coordsR(1) = Scoords(1) + 1

         call MPI_CART_RANK(cartcomm, coordsR, rankR, ierr)

         coordsL(:) = Scoords(:)
         coordsL(1) = Scoords(1) - 1

         call MPI_CART_RANK(cartcomm, coordsL, rankL, ierr)

        
C     	RECEIVE DATA FROM THE LEFT NEIGHBOR 

         tag = 0
         Call MPI_IRECV(EMwait(:,:,:,:),3*3*nx*ny,MPI_REAL,rankL,
     $                          tag,cartcomm,req1,ierr)

         tag = 0
         Call MPI_SEND(E(:,:,:,:,nz-1),3*3*nx*ny,MPI_REAL,rankR,
     $                          tag,cartcomm,ierr)
 
         call MPI_WAIT(req1,stat,ierr)

         E(:,:,:,:,1) = EMwait(:,:,:,:)


C	RECEIVE DATA FROM THE RIGHT NEIGHBOR

         tag = 3
         Call MPI_IRECV(EMwait(:,:,:,:),3*3*nx*ny,MPI_REAL,rankR,
     $                          tag,cartcomm,req1,ierr)

C	SEND DATA TO THE LEFT NEIGHBOR
         tag = 3
         Call MPI_SEND(E(:,:,:,:,2),3*3*nx*ny,MPI_REAL,rankL,
     $                          tag,cartcomm,ierr)

         call MPI_WAIT(req1,stat,ierr)

         E(:,:,:,:,nz) = EMwait(:,:,:,:)

      return
      end



C**********************************************************************
C Subroutine fluxes:
C     Calculate fluxes in evolution equation for the 
C     active particle configuration field
C**********************************************************************

      subroutine flux
 
      use Glob

      implicit none
      include 'mpif.h'

      integer ix,iy,iz,iph,ith,iphth
      integer ixm,iym,izm,iphm,ithm,ixp,iyp,izp,iphp,ithp
      integer i,j,k,is
C      integer stat(MPI_STATUS_SIZE)

      real Vx,Vy,Vz,Vxm,Vym,Vzm,Vxp,Vyp,Vzp
      real delta(3,3),pdot(3),p(3),Vymx,Vypx
       dy2 = dy

      do j = 1,3
         do i = 1,3
            delta(i,j) = 0.00
         enddo
         delta(j,j) = 1.00
      enddo

C Initialize flux

      F = 0.d0

C Calculate convective flux in x-y-z
      do ith = 1,nth
         do iph = 1,nph
            do iz = 2,nz-1

               izm = iz-1
               izp = iz+1

               do iy = 1,ny

                  iym = iy-1
                  if (iym.lt.1)  iym = ny
                  iyp = iy+1
                  if (iyp.gt.ny) iyp = 1

                  do ix = 1,nx

                     ixm = ix-1
                     if (ixm.lt.1)  ixm = nx
                     ixp = ix+1
                     if (ixp.gt.nx) ixp = 1

                     Vxm = sin(theta(ith))*cos(phi(iph)) +
     $                     ux(ixm,iy,iz)
                     Vxp = sin(theta(ith))*cos(phi(iph)) +
     $                     ux(ixp,iy,iz)
                     F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) +
     $               (Vxp*c(ixp,iy,iz,iph,ith)-
     $               Vxm*c(ixm,iy,iz,iph,ith))/2.00/dx

                     Vym = sin(theta(ith))*sin(phi(iph)) +
     $                     uy(ix,iym,iz)
                     Vyp = sin(theta(ith))*sin(phi(iph)) +
     $                     uy(ix,iyp,iz)

                     Vymx = sin(theta(ith))*sin(phi(iph)) +
     $                     uy(ixm,iy,iz)
                     Vypx = sin(theta(ith))*sin(phi(iph)) +
     $                     uy(ixp,iy,iz)

                     F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) +
     $               (Vyp*c(ix,iyp,iz,iph,ith)-
     $               Vym*c(ix,iym,iz,iph,ith))/2.00/dy2 

                     Vzm = cos(theta(ith)) + uz(ix,iy,izm)
                     Vzp = cos(theta(ith)) + uz(ix,iy,izp)
                     F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) +
     $               (Vzp*c(ix,iy,izp,iph,ith)-
     $               Vzm*c(ix,iy,izm,iph,ith))/2.00/dz

                enddo
              enddo
            enddo
         enddo
      enddo

C Calculate convective flux in theta

      do ith = 1,nth
         do iph = 1,nph

            p(1) = sin(theta(ith))*cos(phi(iph))
            p(2) = sin(theta(ith))*sin(phi(iph))
            p(3) = cos(theta(ith))

            do iz = 2,nz-1
               do iy = 1,ny
                  do ix = 1,nx
	 
C Calculate thetadot
               do i = 1,3
               pdot(i) = 0.00
               do k = 1,3
               do j = 1,3
                  pdot(i) = pdot(i) + (delta(i,j)-p(i)*p(j))
     $            *(bEW(j,k,ix,iy,iz)+2.d0*U0*M(j,k,ix,iy,iz))*p(k)
               enddo
               enddo
               enddo

C	f1 = pthetadot, f2 = pphidot or projections of \dot(p) on \hat{theta} and \hat{phi}

               f1(ix,iy,iz,iph,ith) =
     $         cos(theta(ith))*cos(phi(iph))*pdot(1) +
     $         cos(theta(ith))*sin(phi(iph))*pdot(2) -
     $         sin(theta(ith))*pdot(3)

            enddo
            enddo 
            enddo

         enddo
      enddo

      ith = 1
   
      ithp = ith+1
      do iph = 1,nph

         iphm = iph-1
         if (iphm.lt.1) iphm = nph
         iphp = iph+1
         if (iphp.gt.nph) iphp = 1

         iphth = iph + nph/2
         if (iphth.gt.nph) iphth = iphth - nph

         do iz = 2,nz-1
            do iy =1,ny
               do ix =1,nx


               F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) +
     $        (f1(ix,iy,iz,iph,ithp)*c(ix,iy,iz,iph,ithp)*
     $         sin(theta(ithp))
     $        +f1(ix,iy,iz,iphth,ith)*c(ix,iy,iz,iphth,ith)*
     $         sin(theta(ith))) / 2.d0/dth/Sin(theta(ith))
 

 
              end do
           end do
        end do
      end do 


      do ith = 2,nth-1

         ithm = ith-1
         ithp = ith+1
         do iph = 1,nph
            iphm = iph-1
            if (iphm.lt.1) iphm = nph
            iphp = iph+1
            if (iphp.gt.nph) iphp = 1
            do iz = 2,nz-1
               do iy =1,ny
                  do ix =1,nx


             F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) +
     $        ((f1(ix,iy,iz,iph,ithp)*c(ix,iy,iz,iph,ithp)
     *         *sin(theta(ithp))
     $        -f1(ix,iy,iz,iph,ithm)*c(ix,iy,iz,iph,ithm)
     $         *sin(theta(ithm)))/sin(theta(ith))/2.d0/dth)
              
                 enddo 
               enddo               
            enddo
         enddo
      enddo

 
      ith = nth
      
      ithm = ith-1
      do iph = 1,nph
         iphm = iph-1
         if (iphm.lt.1) iphm = nph
         iphp = iph+1
         if (iphp.gt.nph) iphp = 1

         iphth = iph + nph/2
         if (iphth.gt.nph) iphth = iphth - nph

         do iz = 2,nz-1
            do iy =1,ny
               do ix =1,nx

               F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) +
     $        (-f1(ix,iy,iz,iphth,ith)*c(ix,iy,iz,iphth,ith)*
     $         sin(theta(ith))
     $        -f1(ix,iy,iz,iph,ithm)*c(ix,iy,iz,iph,ithm)*
     $         sin(theta(ithm)))/2.d0/dth/sin(theta(ith))

               end do
            end do
          end do
       end do



      do ith = 1,nth
         do iph = 1,nph
            p(1) = sin(theta(ith))*cos(phi(iph))
            p(2) = sin(theta(ith))*sin(phi(iph))
            p(3) = cos(theta(ith))

            do iz = 2,nz-1
               do iy = 1,ny
                  do ix = 1,nx

C Calculate thetadot
               do i = 1,3
               pdot(i) = 0.00
               do k = 1,3
               do j = 1,3
                  pdot(i) = pdot(i) + (delta(i,j)-p(i)*p(j))
     $            *(bEW(j,k,ix,iy,iz)+2.d0*U0*M(j,k,ix,iy,iz))*p(k)
               enddo
               enddo
               enddo

C	f1 = pphidot

               f1(ix,iy,iz,iph,ith) =
     $         -sin(phi(iph))*pdot(1)+cos(phi(iph))*pdot(2)


            enddo
            enddo 
            enddo
         enddo
      enddo


      do ith = 1,nth

         ithm = ith-1
         ithp = ith+1

         do iph = 1,nph
            iphm = iph-1
            if (iphm.lt.1) iphm = nph
            iphp = iph+1
            if (iphp.gt.nph) iphp = 1
            do iz = 2,nz-1
               do iy =1,ny
                  do ix =1,nx


             F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) +
     $        (f1(ix,iy,iz,iphp,ith)*c(ix,iy,iz,iphp,ith)
     $        -f1(ix,iy,iz,iphm,ith)*c(ix,iy,iz,iphm,ith))
     $        /2.d0/dph/sin(theta(ith))

              
                 enddo 
               enddo               
            enddo
         enddo
      enddo


C Calculate diffusive flux in x-y, and diffusive flux in theta

      ith = 1
      
      ithp = ith+1
      do iph = 1,nph
         iphm = iph-1
         if (iphm.lt.1) iphm = nph
         iphp = iph+1
         if (iphp.gt.nph) iphp = 1

         iphth = iph + nph/2
         if (iphth.gt.nph) iphth = iphth - nph

            do iz = 2,nz-1

               izm = iz-1
               izp = iz+1

               do iy = 1,ny

                  iym = iy-1
                  if (iym.lt.1)  iym = ny
                  iyp = iy+1
                  if (iyp.gt.ny) iyp = 1

                  do ix = 1,nx

                     ixm = ix-1
                     if (ixm.lt.1)  ixm = nx
                     ixp = ix+1
                     if (ixp.gt.nx) ixp = 1


               F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) - nu*dfx*(
     $         (c(ixp,iy,iz,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ixm,iy,iz,iph,ith))/dx**2.00
     $         +(c(ix,iyp,iz,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iym,iz,iph,ith))/dy2**2.00
     $         +(c(ix,iy,izp,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iy,izm,iph,ith))/dz**2.00)



               F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) - dfr/nu*(
     $         (1.d0/((sin(theta(ith)))**2.d0))*
     $         (c(ix,iy,iz,iphp,ith)
     $         -2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iy,iz,iphm,ith))/(dph*dph)
     $         +
     $          (sin(theta(ith)+pi/2.d0/nth)*(c(ix,iy,iz,iph,ithp)-
     $         c(ix,iy,iz,iph,ith))
     $         - sin(theta(ith)-pi/2.d0/nth)*(c(ix,iy,iz,iph,ith)-
     $         c(ix,iy,iz,iphth,ith)))
     $         /(dth*dth*sin(theta(ith))))




              end do
           end do
        end do
      end do 


      do ith = 2,nth-1

         ithm = ith-1
         ithp = ith+1
         do iph = 1,nph
            iphm = iph-1
            if (iphm.lt.1) iphm = nph
            iphp = iph+1
            if (iphp.gt.nph) iphp = 1

            do iz = 2,nz-1

               izm = iz-1
               izp = iz+1

               do iy = 1,ny

                  iym = iy-1
                  if (iym.lt.1)  iym = ny
                  iyp = iy+1
                  if (iyp.gt.ny) iyp = 1

                  do ix = 1,nx

                     ixm = ix-1
                     if (ixm.lt.1)  ixm = nx
                     ixp = ix+1
                     if (ixp.gt.nx) ixp = 1



               F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) - nu*dfx*(
     $         (c(ixp,iy,iz,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ixm,iy,iz,iph,ith))/dx**2.00
     $         +(c(ix,iyp,iz,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iym,iz,iph,ith))/dy2**2.00
     $         +(c(ix,iy,izp,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iy,izm,iph,ith))/dz**2.00)



              F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) - dfr/nu*(
     $         (1.d0/(sin(theta(ith))**2.d0))*
     $         (c(ix,iy,iz,iphp,ith)
     $         -2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iy,iz,iphm,ith))/(dph*dph)
     $         +
     $          (sin(theta(ith)+pi/2.d0/nth)*(c(ix,iy,iz,iph,ithp)-
     $         c(ix,iy,iz,iph,ith))
     $         - sin(theta(ith)-pi/2.d0/nth)*(c(ix,iy,iz,iph,ith)-
     $         c(ix,iy,iz,iph,ithm)))
     $         /(dth*dth*sin(theta(ith))))



               
                 enddo 
               enddo               
            enddo
         enddo
      enddo


      ith = nth
      
      ithm = ith-1
      do iph = 1,nph
         iphm = iph-1
         if (iphm.lt.1) iphm = nph
         iphp = iph+1
         if (iphp.gt.nph) iphp = 1

         iphth = iph + nph/2
         if (iphth.gt.nph) iphth = iphth - nph

            do iz = 2,nz-1

               izm = iz-1
               izp = iz+1

               do iy = 1,ny

                  iym = iy-1
                  if (iym.lt.1)  iym = ny
                  iyp = iy+1
                  if (iyp.gt.ny) iyp = 1

                  do ix = 1,nx

                     ixm = ix-1
                     if (ixm.lt.1)  ixm = nx
                     ixp = ix+1
                     if (ixp.gt.nx) ixp = 1




               F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) - nu*dfx*(
     $         (c(ixp,iy,iz,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ixm,iy,iz,iph,ith))/dx**2.00
     $         +(c(ix,iyp,iz,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iym,iz,iph,ith))/dy2**2.00
     $         +(c(ix,iy,izp,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iy,izm,iph,ith))/dz**2.00)



              F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) - dfr/nu*(
     $         (1.d0/(sin(theta(ith))**2.d0))*
     $         (c(ix,iy,iz,iphp,ith)
     $         -2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iy,iz,iphm,ith))/(dph*dph)
     $         +
     $          (sin(theta(ith)+pi/2.d0/nth)*(c(ix,iy,iz,iphth,ith)-
     $         c(ix,iy,iz,iph,ith))
     $         - sin(theta(ith)-pi/2.d0/nth)*(c(ix,iy,iz,iph,ith)-
     $         c(ix,iy,iz,iph,ithm)))
     $         /(dth*dth*sin(theta(ith))))

               end do
            end do
          end do
       end do

      return
      end


C**********************************************************************
C Subroutine euler:
C     Perform initial step using Euler algorithm
C**********************************************************************

      subroutine euler

      use Glob

      implicit none

      integer ix,iy,iz,iph,ith
       dy2 = dy


      call secondmoment
      call MPI_Barrier(cartcomm,ierr)
C Get stress field for Stokes equations 
      call stressfield
      call MPI_Barrier(cartcomm,ierr)
C Get velocity field spectrally
      call velfield
      call MPI_Barrier(cartcomm,ierr)
C
C Get fluxes
      call flux
      call MPI_Barrier(cartcomm,ierr)
C
      c = c - dt * F

      Fo = F

C Update time
      time = time + dt
c      Stime = Stime + dt

      return
      end


C**********************************************************************
C Subroutine update:
C     Advance configuration field using 2nd order Adams-Bashforth
C**********************************************************************

      subroutine update
 
      use Glob

      implicit none
      include 'mpif.h'



      integer ix,iy,iz,iph,ith
      integer i,j,k,is,iremesh
      real imp,imp1
       dy2 = dy
c      dy2 = dy*sqrt(1.d0+(Shear*Stime)**2.d0)


C COPY TO GHOST CELLS FROM NEIGHBORS


C******************************************************************
C	        COORDINATES AND RANKS OF NEIGHBORS
 
         call MPI_CART_COORDS(cartcomm,rank,1,Scoords,ierr)

         coordsR(:) = Scoords(:)
         coordsR(1) = Scoords(1) + 1

         call MPI_CART_RANK(cartcomm, coordsR, rankR, ierr)


         coordsL(:) = Scoords(:)
         coordsL(1) = Scoords(1) - 1

         call MPI_CART_RANK(cartcomm, coordsL, rankL, ierr)


C********************************************************************

         tag = 0
         Call MPI_IRECV(cwait(:,:,:,:),nx*ny*nph*nth,MPI_REAL,rankL,
     $                          tag,cartcomm,req,ierr)
 
         csend(:,:,:,:) = c(:,:,nz-1,:,:)

         tag = 0
         Call MPI_SEND(csend(:,:,:,:),nx*ny*nph*nth,MPI_REAL,rankR,
     $                          tag,cartcomm,ierr)

         call MPI_WAIT(req,stat,ierr)

         c(:,:,1,:,:) = cwait(:,:,:,:)

C     	RECEIVE DATA FROM THE LEFT NEIGHBOR 

C*********************************************************************


C********************************************************************

         tag = 1
         Call MPI_IRECV(cwait(:,:,:,:),nx*ny*nph*nth,MPI_REAL,rankR,
     $                          tag,cartcomm,req,ierr)
 
 
         csend(:,:,:,:) = c(:,:,2,:,:)


         tag = 1
         Call MPI_SEND(csend(:,:,:,:),nx*ny*nph*nth,MPI_REAL,rankL,
     $                          tag,cartcomm,ierr)

         call MPI_WAIT(req,stat,ierr)

         c(:,:,nz,:,:) = cwait(:,:,:,:)

C**********************************************************************
C********************************************************************
C***********************************************************************

      call secondmoment
      Call MPI_Barrier(cartcomm,ierr)

C Get stress field for Stokes equations 

      call stressfield
      Call MPI_Barrier(cartcomm,ierr)

C Get velocity field spectrally

      call velfield
      Call MPI_Barrier(cartcomm,ierr)

C Get fluxes

      call flux
      Call MPI_Barrier(cartcomm,ierr)

C Update configuration field using explicit Euler scheme
C and copy flux
C      omg = 0.0100

c      imp1 = dt*dfx*((1+(Shear*Stime)**2.d0)/(dx*dx)
c     $                + 1/(dy2*dy2)+1/(dz*dz))

c      Do ith = 1,nth

c      imp = imp1
c     $    + dt*dfr/(dph*Sin(theta(ith)))**2.d0


c      c(:,:,:,:,ith)  = ( c(:,:,:,:,ith)*(1+imp) - 0.5d0*dt*
c     $     (  3.d0*F(:,:,:,:,ith)-Fo(:,:,:,:,ith))) / (1+3*imp)
c     $       
c
c      End Do

      c  = c - 0.5d0*dt*(3.d0*F - Fo)

      Fo = F

C Update time
      time = time + dt
c      Stime = Stime + dt


      return
      end



C*******************************************************************
C Function fourn:
C     Computes the fast Fourier transform of a multidimensional
C     array of complex numbers.
C     From Numerical Recipes, Press et al., 2201.
C*******************************************************************

      SUBROUTINE fourn(data,nn,ndim,isign) 

      INTEGER isign,ndim,nn(ndim) 
      real data(*) 
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2 
      INTEGER ip3,k1,k2,n,nprev,nrem,ntot 
      real tempi,tempr 
      REAL theta,wi,wpi,wpr,wr,wtemp 

      ntot=1 
      do idim=1,ndim 
         ntot=ntot*nn(idim) 
      enddo 
      nprev=1 
      do idim=1,ndim 
         n=nn(idim) 
         nrem=ntot/(n*nprev) 
         ip1=2*nprev 
         ip2=ip1*n 
         ip3=ip2*nrem 
         i2rev=1 
         do i2=1,ip2,ip1 
            if(i2.lt.i2rev) then 
               do i1=i2,i2+ip1-2,2 
                  do i3=i1,ip3,ip2 
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
            theta=isign*6.2831852071795900/(ifp2/ip1) 
            wpr=-2.00*sin(0.500*theta)**2 
            wpi=sin(theta) 
            wr=1.00 
            wi=0.00 
            do i3=1,ifp1,ip1 
               do i1=i3,i3+ip1-2,2 
                  do i2=i1,ip3,ifp2 
                     k1=i2 
                     k2=k1+ifp1 
                     tempr=(wr)*data(k2)-(wi)*data(k2+1) 
                     tempi=(wr)*data(k2+1)+(wi)*data(k2) 
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
      END


C*******************************************************************
C Function ran:   
C     Random number generator.
C*******************************************************************

      function ran(idum)

      implicit none
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran,am,eps,rnmx
      parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     &     ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     &     ir2=3791,ntab=16,ndiv=1+imm1/ntab,eps=1.2E-7,rnmx=1.-eps)
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
