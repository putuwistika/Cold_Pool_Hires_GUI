c - DTDM2 adding moisture
c  namelist moisture input ignored, do WK rh in soundiing
c  NOT in anelastic sub YET

c fv should be added to i=2 and nx...

c========================================================================
c  DTDM = Dynamics and Thermodynamics Demonstration Model
c  History
c      Version 1.0 June, 2006
c      Version 1.1 July, 2006
c      Version 1.2 April, 2007
c      Version 1.2.2 included linear friction on U' at surface
c      Version 2.0 May, 2007 
c         - started including moisture
c         - atmospheric heat source modified
c      Version 2.0.1 May, 2007
c
c  Requirements:
c    A Fortran 77/90 compiler and GrADS.  g77 and g95 suffice.
c
c  Author:
c    Robert Fovell
c    Department of Atmospheric and Oceanic Sciences
c    University of California, Los Angeles
c    rfovell@ucla.edu
c
c  For demonstration purposes only.  No guarantees or warranties.
c  Some physical processes or constants may be exaggerated for
c    computational efficiency and/or demonstration purposes.
c  Not for research-quality applications.
c
c  Please send bug reports to rfovell@ucla.edu
c
c========================================================================


c dtdm.f - based on thermal3.f
c    ** if iaccel=1 compute DU and DWDT total forcing and
c        write out in place of U and W
c    ** if iaccel=2 write out buoy forcing Fbx, Fbz
c    ** if iaccel=3 write out dynamic forcing Fdx, Fdz
c open the lateral boundaries

c======================================================================
      program dtdm
      include 'storage.txt'

c
c get input parameters and define some constants
c

      call readnl   ! read in input in namelist format
c      call readin  ! this read in the old input scripts
c
c open files for grads output.  existing files are overwritten.
c

c      open(44,file=ctlfile,status='unknown')
      open(66,file=datfile,status='unknown',
     1 form='unformatted')
      igradscnt=0

c
c set up base state
c
 
      call sounding

c
c set up cases
c

      call zeroit
      if(ithermal.eq.1) call thermal
      if(ishflux.eq.1) call surface
      if(istrfcn.eq.1) call streamfcn
      if(ihsrc.eq.1) call heatsrc
      if(iseabreeze.eq.1) call seabreeze
      if(icoolzone.ge.1) call coolzone
      
c
c set up anelastic pressure solver (blktri)
c

      dxr2=1./(dx*dx)
      do  i=2,nx-1
       am(i-1)=dxr2
       cm(i-1)=dxr2
       bm(i-1)=-am(i-1)-cm(i-1)
      enddo
c commenting next 4 lines produces periodic lateral bcs
       bm(1)=-cm(1)
       bm(nx-2)=-am(nx-2)
       am(1)=0.
       cm(nx-2)=0.
      do  k=2,nz-1
       temp1=1./(rhou(k)*tinit(k)*dz)
        if(k.ne.2) then
         an(k-1)=rhow(k)*0.5*(tinit(k)+tinit(k-1))*temp1/dz
        else
         an(k-1)=0.
        end if
        if(k.ne.nz-1) then
         cn(k-1)=rhow(k+1)*0.5*(tinit(k+1)+tinit(k))*temp1/dz
        else
         cn(k-1)=0.
        end if
      bn(k-1)=-an(k-1)-cn(k-1)
      enddo
      cn(nz-2)=0.

      call blktri(0,1,nz-2,an,bn,cn,1,nx-2,am,bm,cm,nxm,utemp,ier,wa)
      write(6,17002) ier,wa(1)
17002 format(1x,'blktri initial ier = ',i1,' wa(1) is ',e15.6)
      if(linear.eq.1.and.ianelastic.eq.1) stop ' anelastic - nonlinear'
      if(linear.eq.1) ipressure=0 ! overwrite, if needed
      if(ipressure.eq.1.or.ianelastic.eq.1) call pressure
c 
      pinit=0
     
      if(pinit.eq.1)then
c if initial perturbation was applied, hydrostatically rebalance

      do i=2,nx-1
       pi(i,nz-1)=0.
       do k=nz-2,2,-1
        tup=th(i,k+1)/(tinit(k+1)**2)
        tdn=th(i, k )/(tinit( k )**2)
        pi(i,k)=pi(i,k+1)-0.5*(g/cp)*(tup+tdn)*dz
        pim(i,k)=pi(i,k)
       enddo
       pi(i,1)=pi(i,2)
       pim(i,1)=pim(i,2)
      enddo
      else
      print *,' using blktri to get initial pressure field'
      do i=2,nx-1
       do k=2,nz-1
        pi(i,k)=ptot(i,k)
        pim(i,k)=pi(i,k)
       enddo
      enddo
      endif
c
c force lateral boundaries, tho should't be necessary here
c

      call scalarbc(pi)
      call scalarbc(pim)
c
c dimensional pressure perturbation, for plotting purposes
c

      do i=2,nx-1
      do k=1,nz
       pprt(i,k)=rhou(k)*cp*tinit(k)*pi(i,k)
      enddo
      enddo
      
      
c initial time
      time=0
      call prepare_grads
c write out temporary control file, so we can monitor model as it runs
      call write_to_grads_ctl(nx,nz,ipltint,0)
c start integration
      d2t=dt
      nmax=timend/dt
      n=0
      iplt=0
      plot_tmp=plot/dt
      ipltint=ifix(plot_tmp)  ! plotting interval in time steps
      print *,' plot interval ',plot,' sec or ',ipltint,' time steps'
1000  n=n+1
      time=time+dt
      iplt=iplt+1
      dtx=d2t/dx
      dtz=d2t/dz
      if(time.gt.timend) go to 1001
c time factor for oscillating heat source
      timfac=1.0+0.25*sin(h_freq*time)
      sb_timfac=sin(time*sb_freq)
      s_timfac=cos(time*s_freq)

c      if(time.gt.5000) s_timfac=0.

      if(ianelastic.eq.1)then
        call anelastic ! anelastic model, leapfrog
      else
        call leapfrog  ! compressible model, leapfrog
      endif

      thmax=0
       do k=1,nz
       do i=1,nx
        if(abs(th(i,k)).gt.thmax)then
         thmax=abs(th(i,k))
         thmi=i
         kthm=k
        endif
       enddo
       enddo
       print *,' time ',time,' thmax ',thmax

c history calls
      if(iplt.eq.ipltint)then
        iplt=0
        call prepare_grads
c       endif
c statistics
       thmax=0.
       wmax=0.
       umax=0.
       thmi=0
       kthm=0
       wmi=0
       wmk=0
       thsum=0.
       do k=1,nz
       do i=1,nx
        if(abs(th(i,k)).gt.thmax)then
         thmax=abs(th(i,k))
         thmi=i
         kthm=k
        endif
        if(abs(w(i,k)).gt.wmax)then
         wmax=abs(w(i,k))
         wmi=i
         wmk=k
        endif
        if(abs(u(i,k)).gt.umax) umax=abs(u(i,k))
       enddo
       enddo
       do i=2,nx-1
       do k=2,nz-1
        thsum=thsum+rhou(k)*th(i,k)*dx*dz
       enddo
       enddo
       print *,' time ',time,' thmax ',thmax,thmi,kthm
       print *,' wmax ',wmax,wmi,wmk,' thsum ',thsum
       print *,' umax ',umax
       endif ! iplt
             
c double the time step
      d2t=dt+dt
      go to 1000
c end of integration
1001  continue


c write out grads ctl file at end of model run
      call write_to_grads_ctl(nx,nz,ipltint,1)
      print *,' =========================================== '
      print *,'  normal model stop'
      print *,' =========================================== '
      stop   
      end

c
c subroutine prepare_grads
c
      subroutine prepare_grads
      include 'storage.txt'
      ngradsvars=10+istrfcn+ihsrc+iseabreeze+ienable_v
     1 +9*ipressure+imoist*4
      do i=1,nx
       do k=1,nz
        temp(i,k)=rhou(k)
       enddo
      enddo
      call write_to_grads_dat(temp,nxm,nzm,nx,nz,zero,0.,0)   ! mean density at U levels
      do i=1,nx
       do k=1,nz
        temp(i,k)=(float(k)-1.5)*dz
       enddo
      enddo
      call write_to_grads_dat(temp,nxm,nzm,nx,nz,zero,0.,0)   ! height at U levels
c some calls pass a non-zero base state array to augment to the 2D field
      call write_to_grads_dat(u,nxm,nzm,nx,nz,zero,0.,1)    ! full U is predicted in model
      call write_to_grads_dat(u,nxm,nzm,nx,nz,uinit,-1.,1)  ! pert u made by subtracting uinit
      if(ienable_v.eq.1)then
       call write_to_grads_dat(v,nxm,nzm,nx,nz,zero,0.,1)   ! v component
      endif
      call write_to_grads_dat(w,nxm,nzm,nx,nz,zero,0.,2)
      call write_to_grads_dat(th,nxm,nzm,nx,nz,tinit,1.,0)  ! full pot temp
      call write_to_grads_dat(th,nxm,nzm,nx,nz,zero,0.,0)   ! pert pot temp
      if(imoist.eq.1)then
       call write_to_grads_dat(qv,nxm,nzm,nx,nz,qinit,1.,0)  ! full q
       call write_to_grads_dat(qv,nxm,nzm,nx,nz,zero,0.,0)   ! pert q
       call calc_rh ! rel hum returns in array "temp"
       call write_to_grads_dat(temp,nxm,nzm,nx,nz,zero,0.,0) ! rh 
       call write_to_grads_dat(accum_cond,nxm,nzm,nx,nz,zero,0.,0) ! conden sum
       do i=1,nx
        do k=1,nz
         accum_cond(i,k)=0. ! zero out cumulated conden array
        enddo
       enddo
      endif
      call write_to_grads_dat(pi,nxm,nzm,nx,nz,pk,1.,0)     ! full ndim pressure        
      call write_to_grads_dat(pi,nxm,nzm,nx,nz,zero,0.,0)   ! pert ndim pressure
      do i=1,nx
       do k=1,nz
        temp(i,k)=0.01*rhou(k)*cp*tinit(k)*pi(i,k)
       enddo
      enddo
      call write_to_grads_dat(temp,nxm,nzm,nx,nz,zero,0.,0)   ! pert pressure in millibars
      if(ipressure.eq.1)then
       do i=1,nx
        do k=1,nz
         pbyc(i,k)=0.01*rhou(k)*cp*tinit(k)*pbyc(i,k) ! convert to mb
         pdyn(i,k)=0.01*rhou(k)*cp*tinit(k)*pdyn(i,k) ! convert to mb
         ptot(i,k)=0.01*rhou(k)*cp*tinit(k)*ptot(i,k) ! convert to mb
        enddo
       enddo
       call write_to_grads_dat(pbyc,nxm,nzm,nx,nz,zero,0.,0)   ! buoyancy pressure
       call write_to_grads_dat(pdyn,nxm,nzm,nx,nz,zero,0.,0)   ! dynamic pressure
       call write_to_grads_dat(ptot,nxm,nzm,nx,nz,zero,0.,0)   ! total pressure
       call write_to_grads_dat(dudtd,nxm,nzm,nx,nz,zero,0.,0)   ! U accel - dynamic
       call write_to_grads_dat(dudtb,nxm,nzm,nx,nz,zero,0.,0)   ! U accel - buoyancy
       call write_to_grads_dat(dudtt,nxm,nzm,nx,nz,zero,0.,0)   ! U accel - total
       call write_to_grads_dat(dwdtd,nxm,nzm,nx,nz,zero,0.,0)   ! W accel - dynamic
       call write_to_grads_dat(dwdtb,nxm,nzm,nx,nz,zero,0.,0)   ! W accel - buoyancy
       call write_to_grads_dat(dwdtt,nxm,nzm,nx,nz,zero,0.,0)   ! W accel - total
      endif
      if(istrfcn.eq.1)then
       do i=1,nx
        do k=1,nz
         temp(i,k)=strfcn(i,k)*s_timfac
        enddo
       enddo
       call write_to_grads_dat(temp,nxm,nzm,nx,nz,zero,0.,3) ! streamfunction
       endif
      if(ihsrc.eq.1)then
       do i=1,nx
        do k=1,nz
         temp(i,k)=hsrc(i,k)*timfac
        enddo
       enddo
       call write_to_grads_dat(temp,nxm,nzm,nx,nz,zero,0.,0) ! heat source
      endif
      if(iseabreeze.eq.1)then
       do i=1,nx
        do k=1,nz
         temp(i,k)=sb_hsrc(i,k)*sb_timfac
        enddo
       enddo
       call write_to_grads_dat(temp,nxm,nzm,nx,nz,zero,0.,0) ! seabreeze heat source
      endif
      igradscnt=igradscnt+1
      print *,' done writing grads write number ',igradscnt
      return
      end
c 
c subroutine write_to_grads_dat
c
      subroutine write_to_grads_dat(array,nxx,nzx,nxg,nzg,zmean,zfactor,
     1 iavg)
      dimension array(nxx,nzx),zmean(nzx)
      dimension dummy(nxg,1,nzg)
      do i=1,nxg
       do k=1,nzg
        dummy(i,1,k)=0.
       enddo
      enddo
      if(iavg.eq.0)then ! no averaging needed
       do i=1,nxg
        do k=1,nzg
         dummy(i,1,k)=array(i,k)+zfactor*zmean(k)
        enddo
       enddo
      else if(iavg.eq.1)then ! average in x-direction
       do k=1,nzg
        do i=1,nxg-1
        dummy(i,1,k)=0.5*(array(i+1,k)+array(i,k))+zfactor*zmean(k)
        enddo
       enddo
      else if(iavg.eq.2)then ! average in z-direction
       do i=1,nxg
        do k=1,nzg-1
         dummy(i,1,k)=0.5*(array(i,k+1)+zfactor*zmean(k+1)
     1                    +array(i, k )+zfactor*zmean( k ))
        enddo
       enddo
      else if(iavg.eq.3)then ! average from strfcn point
       do i=2,nxg-1
        do k=2,nzg-1
         top_pair=0.5*(array(i+1,k+1)+array(i,k+1))+zfactor*zmean(k+1)
         bot_pair=0.5*(array(i+1, k )+array(i, k ))+zfactor*zmean( k )
         dummy(i,1,k)=0.5*(top_pair+bot_pair)
        enddo
       enddo
      endif ! iavg
      do k=2,nzg-1
       write(66) ((dummy(i,j,k),i=2,nxg-1),j=1,1)
      enddo
      return
      end

c
c subroutine write_to_grads_ctl
c
      subroutine write_to_grads_ctl(nxg,nzg,ipltint,ifinal)
      include 'storage.txt'
      close(44)
      open(44,file=ctlfile,status='unknown')
      print *,' writing grads control file'

      gradstinc=float(ipltint)*dt/60.  ! minutes
      igradstinc=ifix(gradstinc) ! if not integer number of min, time is counted wrong
      print *,' nxg= ',nxg,' nzg= ',nzg,' gradstinc ',gradstinc
      
      write(44,900) datfile ! grads does not like underscores in ctl, dat file names
 900  format('DSET ^',a80)
      write(44,901)
 901  format('TITLE DTDM demo simulation',/,
     1 'OPTIONS sequential',/,'UNDEF -9.99E33')
      if(byteswap.eq.1) write(44,902)
 902  format('OPTIONS byteswapped')
      write(44,100) nxg-2,dx/2/1000.,dx/1000. ! dx in km
 100  format('XDEF ',i3,' LINEAR ',2f10.5,/,'YDEF 1 LINEAR 1.0 1.0')
      write(44,101) nzg-2
 101  format('ZDEF ',i3,' levels')
      do k=2,nzg-1
       write(44,102) (float(k)-1.5)*dz/1000.
 102   format(f10.3)
      enddo
      if(igradstinc.lt.1)then
       print *,' *** WARNING - since plotting interval less ',
     1  'than one minute, GrADS will not report time ',
     1  'correctly.'
       print *,' *** Plot interval recorded as 1 min in ',
     1  'control file'
       igradstinc=1
      endif
      if(ifinal.eq.0)then
       itemp=ifix(timend/plot)+1 ! expected number of writes
      else
       itemp=igradscnt
      endif
      write(44,103) itemp,igradstinc
 103  format('TDEF ',i5,' LINEAR 00:00Z01JAN2000 ',i3,'mn')
      write(44,104) ngradsvars
 104  format('VARS ',i4)
      write(44,205) nzg-2
 205  format('rho ',i3,' 00   mean density')
      write(44,206) nzg-2
 206  format('hgt ',i3,' 00   height')
      write(44,105) nzg-2
 105  format('u  ',i3,' 00   horizontal velocity')
      write(44,106) nzg-2
 106  format('up ',i3,' 00   pert horizontal velocity')
      if(ienable_v.eq.1)then
      write(44,199) nzg-2
 199  format('v  ',i3,' 00   north-south velocity')
      endif
      write(44,107) nzg-2
 107  format('w  ',i3,' 00   vertical velocity')
      write(44,108) nzg-2
 108  format('th ',i3,' 00   potential temperature')
      write(44,109) nzg-2
 109  format('thp ',i3,' 00   pert potential temperature')
      if(imoist.eq.1)then
      write(44,120) nzg-2
 120  format('qv ',i3,' 00   vapor specific humidity')
      write(44,121) nzg-2
 121  format('qvp ',i3,' 00   pert vapor specific humidity')
      write(44,122) nzg-2
 122  format('rh  ',i3,' 00   relative humidity')
      write(44,123) nzg-2
 123  format('cond',i3,' 00   condensation in last interval')
      endif
      write(44,110) nzg-2
 110  format('pi ',i3,' 00   ndim pressure')
      write(44,111) nzg-2
 111  format('pip ',i3,' 00  pert ndim pressure')
      write(44,112) nzg-2
 112  format('ppmb ',i3,' 00  pert pressure in millibars')
      if(ipressure.eq.1)then
      write(44,113) nzg-2
 113  format('pbyc ',i3,' 00  buoyancy pressure in millibars')
      write(44,114) nzg-2
 114  format('pdyn ',i3,' 00  dynamic pressure in millibars')
      write(44,115) nzg-2
 115  format('ptot ',i3,' 00  total pressure in millibars')
      write(44,140) nzg-2
 140  format('dudtd ',i3,' 00  U acceleration - dynamic')
      write(44,141) nzg-2
 141  format('dudtb ',i3,' 00  U acceleration - buoyancy')
      write(44,142) nzg-2
 142  format('dudtt ',i3,' 00  U acceleration - total')
      write(44,143) nzg-2
 143  format('dwdtd ',i3,' 00  W acceleration - dynamic')
      write(44,144) nzg-2
 144  format('dwdtb ',i3,' 00  W acceleration - buoyancy')
      write(44,145) nzg-2
 145  format('dwdtt ',i3,' 00  W acceleration - total')
      endif
      if(istrfcn.eq.1) write(44,116) nzg-2
 116  format('str ',i3,' 00  streamfunction')
      if(ihsrc.eq.1) write(44,117) nzg-2
 117  format('hsrc ',i3,' 00  heat source')
      if(iseabreeze.eq.1) write(44,118) nzg-2
 118  format('sb_hsrc ',i3,' 00  seabreeze heat source')
      write(44,999)
 999  format('ENDVARS')
      close(44) 

      return
      end

c
c subroutine pressure
c

      subroutine pressure
      include 'storage.txt'
      do k=2,nz-1
          advu(2,k)=
     1      -amin1((u(2 ,k)-cstar),0.)*(um(3,k)-um(2,k))*dtx
          advu(nx,k)=
     1      -amax1((u(nx,k)+cstar),0.)*(um(nx,k)-um(nx-1,k))*dtx
      do i=3,nx-1
         advu(i,k)=-.25*((u(i+1,k)+u(i,k))**2
     1                          -(u(i-1,k)+u(i,k))**2)/dx
     2                  -.25*(rhow(k+1)*(w(i,k+1)+w(i-1,k+1))
     2                                     *(u(i,k+1)+u(i,k))
     3                         -rhow( k )*(w(i,k  )+w(i-1,k  ))
     3                                   *(u(i,k-1)+u(i,k)))/
     4      (dz*rhou(k))
         if(fcoef.ne.0..and.k.eq.2.and.i.gt.icoast)then
           if(ienable_v.eq.1) stop 'fric-no_coriolis'
            fric=-fcoef*(um(i,k)-uinit(k))  ! linear friction on perturbations
         else
           fric=0.
         endif
         advu(i,k)=advu(i,k)+fric
       enddo     
       enddo 
       call scalarbc(advu)
      
c w eqn
c w is zero at k=2 and nz
       do k=3,nz-1
       do i=2,nx-1
        buoy(i,k)=g*0.5*(th(i,k)/tinit(k) + th(i,k-1)/tinit(k-1))
        if(moist.eq.1) 
     1  buoy(i,k)=buoy(i,k)+g*0.5*0.61*(qv(i,k)+qv(i,k-1))
        advw(i,k)=-.25*((u(i+1,k)+u(i+1,k-1))
     1                   *(w(i+1,k)+w(i,k))
     2             -(u(i  ,k)+u(i  ,k-1))
     2                   *(w(i-1,k)+w(i,k)))/dx
     3    -.25*(rhou( k )*(w(i,k+1)+w(i,k))**2
     4             -rhou(k-1)*(w(i,k-1)+w(i,k))**2)/(rhow(k)*dz)
       enddo
       enddo

c boundary conditions on buoy
       call wbc(buoy)
       call wbc(advw)    
     
c pressure decomposition
      do k=2,nz-1
       a1=1./(rhou(k)*dz)
       a2=1./(cp*tinit(k))
       do i=2,nx-1
        byctrm(i-1,k-1)=(buoy(i,k+1)*rhow(k+1)-buoy(i,k)*rhow(k))*a1*a2
        dyntrm(i-1,k-1)=((advu(i+1,k)-advu(i,k))/dx
     1               +(advw(i,k+1)*rhow(k+1)-advw(i,k)*rhow(k))*a1)*a2
        alltrm(i-1,k-1)=dyntrm(i-1,k-1)+byctrm(i-1,k-1)
       enddo
      enddo

      call blktri(1,1,nz-2,an,bn,cn,1,nx-2,am,bm,cm,nxm,byctrm,ier,wa)
      if(ier.ne.0) then
       print *,' pbyc ier ',ier
       stop 'blktri'
      endif
      call blktri(1,1,nz-2,an,bn,cn,1,nx-2,am,bm,cm,nxm,dyntrm,ier,wa)
      if(ier.ne.0) then
       print *,' pdyn ier ',ier
       stop 'blktri'
      endif
      call blktri(1,1,nz-2,an,bn,cn,1,nx-2,am,bm,cm,nxm,alltrm,ier,wa)
      if(ier.ne.0) then
       print *,' ptot ier ',ier
       stop 'blktri'
      endif
      do k=2,nz-1
       do i=2,nx-1
        pbyc(i,k)=byctrm(i-1,k-1)-byctrm(1,1)
        pdyn(i,k)=dyntrm(i-1,k-1)-dyntrm(1,1)
        ptot(i,k)=alltrm(i-1,k-1)-alltrm(1,1)
       enddo
      enddo
      call scalarbc(pbyc)
      call scalarbc(pdyn)   
      call scalarbc(ptot)  
      do k=2,nz-1
         do i=2,nx-1
            dudtd(i,k)=-cp*tinit(k)*(pdyn(i,k)-pdyn(i-1,k))/dx
c            if(fcoef.ne.0..and.k.eq.2.and.i.gt.icoast)then ! removed since it was double counting
c             dudtd(i,k)=dudtd(i,k)-fcoef*(um(i,k)-uinit(k))  ! linear friction on perturbations
c            endif
            dudtb(i,k)=-cp*tinit(k)*(pbyc(i,k)-pbyc(i-1,k))/dx
            dudtt(i,k)=dudtb(i,k)+dudtd(i,k)
         enddo
      enddo
      do i=2,nx-1
          do k=3,nz-1
            dwdtd(i,k)=-0.5*cp*(tinit(k)+tinit(k-1))*
     1         (pdyn(i,k)-pdyn(i,k-1))/dz
            dwdtb(i,k)=-0.5*cp*(tinit(k)+tinit(k-1))*
     1         (pbyc(i,k)-pbyc(i,k-1))/dz
     1         +g*0.5*(th(i,k)/tinit(k) + th(i,k-1)/tinit(k-1))
            if(moist.eq.1) 
     1      dwdtb(i,k)=dwdtb(i,k)+g*0.5*0.61*(qv(i,k)+qv(i,k-1))
            dwdtt(i,k)=dwdtd(i,k)+dwdtb(i,k)
          enddo
      enddo
      call wbc(dwdtb)
      call wbc(dwdtd)
      call wbc(dwdtt)
      return
      end


c subroutine scalarbc
c
      subroutine scalarbc(x)
      include 'storage.txt'
      dimension x(nxm,nzm)
c zero gradient top and bottom
      do i=2,nx-1
       x(i,1)=x(i,2)
       x(i,nz)=x(i,nz-1)
      enddo
c now k=1,nz has been done
c zero gradient lateral boundaries
      do k=1,nz
       x(1,k)=x(2,k)
       x(nx,k)=x(nx-1,k)
      enddo
      return
      end
c
c subroutine wbc
c
      subroutine wbc(x)
      include 'storage.txt'
      dimension x(nxm,nzm)
c w is zero at k=2 and nz. also k=1
      do i=2,nx-1
       x(i,2)=0.
       x(i,1)=0.
       x(i,nz)=0.
      enddo
c now k=1,nz has been done
c zero gradient lateral boundaries
      do k=1,nz
       x(1,k)=x(2,k)
       x(nx,k)=x(nx-1,k)
      enddo    
      return
      end

c
c subroutine calc_rh
c

      subroutine calc_rh ! rh returns in array "temp"
      include 'storage.txt'
      do k=1,nz
      do i=1,nx
       temp(i,k)=0.
      enddo
      enddo
      do k=2,nz-1
       p0=(pk(k)**xki)*psl
       do i=2,nx-1
        tempx=(th(i,k)+tinit(k))*pk(k)
        temp(i,k)=(380./p0)*exp(17.27*(tempx-273.)
     1                   /(tempx-36.))
        temp(i,k)=(qv(i,k)+qinit(k))/temp(i,k)
       enddo
      enddo
      return
      end

c
c subroutine satadj
c

      subroutine satadj ! do saturation adjustment for NCCM
      include 'storage.txt'

c in the NCCM, cloud water generated is ignored.  cloudy areas
c  are diagnosed by where RH = 100%.  There is no precip.

      do k=2,nz-1
       p0=(pk(k)**xki)*psl 
       do i=2,nx-1
        qvp(i,k) = amax1(qvp(i,k),-qinit(k)) ! remove negative vapor
        tempqv = qvp(i,k)+qinit(k)           ! vapor mixing ratio
        tempt  = (thp(i,k)+tinit(k))*pk(k)   ! temp from pot. temp
        qvs = (380./p0)*exp(17.27*(tempt-273.)/(tempt-36.))
        rsub1 = qvs*(17.27*237.*hlvcp)
     1            /((tempt-36.)*(tempt-36.))
        conden = (tempqv-qvs)/(1.+rsub1)     ! condensation production
        conden = amax1(conden,0.)            ! no evaporation in NCCM
        conden = amin1(conden,tempqv)        ! limit by available qv
        thp(i,k) = thp(i,k)+conden*hlvcp/pk(k) ! latent heating
        qvp(i,k) = qvp(i,k)-conden           ! vapor sink
        accum_cond(i,k)=accum_cond(i,k)+conden ! conden summed over
                                                   ! plotting interval
c        tempt1=(thp(i,k)+tinit(k))*pk(k) 
c        tempqv1=qvp(i,k)+qinit(k)
c        qvs = (380./p0)*exp(17.27*(tempt1-273.)/(tempt1-36.))
c        rh=tempqv1/qvs
c        if(rh.gt.1.0001)then
c         print *,' ik ',i,k,' rh ',rh,' thp ',thp(i,k),' qvp ',
c     1    qvp(i,k)
c        endif
       enddo
      enddo
      call scalarbc(thp)
      call scalarbc(qvp)
      return
      end

c 
c subroutine leapfrog
c
      subroutine leapfrog
      include 'storage.txt'
c u eqn
       if(linear.ne.1)then
       do k=2,nz-1
          up(2 ,k)=um(2 ,k)
     1      -amin1((u(2 ,k)-cstar),0.)*(um(3,k)-um(2,k))*dtx
          up(nx,k)=um(nx,k)
     1      -amax1((u(nx,k)+cstar),0.)*(um(nx,k)-um(nx-1,k))*dtx
        do i=3,nx-1
         up(i,k)=um(i,k)-.25*dtx*((u(i+1,k)+u(i,k))**2
     1                           -(u(i-1,k)+u(i,k))**2)
     2                  -.25*dtz*(rhow(k+1)*(w(i,k+1)+w(i-1,k+1))
     2                                     *(u(i,k+1)+u(i,k))
     3                           -rhow( k )*(w(i,k  )+w(i-1,k  ))
     3                                   *(u(i,k-1)+u(i,k)))/rhou(k)
     4                 -dtx*cp*tinit(k)*(pi(i,k)-pi(i-1,k))
     5    -dtz*s_timfac*(strfcn(i,k+1)-strfcn(i,k))/rhou(k)
     6    +dkx*(d2t/dx**2)*(um(i+1,k)-2.*um(i,k)+um(i-1,k))
c          temp(i,k)=-dtz*(strfcn(i,k+1)-strfcn(i,k))/rhou(k)
c Coriolis
          if(ienable_v.eq.1) up(i,k)=up(i,k)+d2t*corf*v(i,k)
          up(i,k)=up(i,k)
     1     +dkz*(d2t/dz**2)*(um(i,k+1)-2.*um(i,k)+um(i,k-1)
     2                      -uinit(k+1)+2.*uinit(k)-uinit(k-1))
         if(iaccel.eq.1) dudt(i,k)=-cp*tinit(k)*(pi(i,k)-pi(i-1,k))/dx
        enddo     
       enddo 
       else ! linear model
       do k=2,nz-1
          up(2 ,k)=um(2 ,k)
     1      -amin1((uinit(k)-cstar),0.)*(um(3,k)-um(2,k))*dtx
          up(nx,k)=um(nx,k)
     1      -amax1((uinit(k)+cstar),0.)*(um(nx,k)-um(nx-1,k))*dtx
        do i=3,nx-1
         up(i,k)=um(i,k)-.25*dtx*uinit(k)*(u(i+1,k)-u(i-1,k))
     2                  -.25*dtz*(rhow(k+1)*(w(i,k+1)+w(i-1,k+1))
     2                                     *(uinit(k+1)+uinit(k))
     3                           -rhow( k )*(w(i,k  )+w(i-1,k  ))
     3                                   *(uinit(k-1)+uinit(k)))/rhou(k)
     4                 -dtx*cp*tinit(k)*(pi(i,k)-pi(i-1,k))
     5    -dtz*s_timfac*(strfcn(i,k+1)-strfcn(i,k))/rhou(k)
     6    +dkx*(d2t/dx**2)*(um(i+1,k)-2.*um(i,k)+um(i-1,k))
c Coriolis
          if(ienable_v.eq.1) up(i,k)=up(i,k)+d2t*corf*v(i,k)
          up(i,k)=up(i,k)
     1     +dkz*(d2t/dz**2)*(um(i,k+1)-2.*um(i,k)+um(i,k-1)
     2                      -uinit(k+1)+2.*uinit(k)-uinit(k-1))
c experimental - add friction... not for coriolis yet, not for linear
         if(fcoef.ne.0..and.k.eq.2.and.i.gt.icoast)then
           if(ienable_v.eq.1) stop 'fric-no_coriolis'
            fric=-fcoef*(um(i,k)-uinit(k))  ! linear friction on perturbations
         else
           fric=0.
         endif
         up(i,k)=up(i,k)+fric*d2t
         if(iaccel.eq.1) dudt(i,k)=-cp*tinit(k)*(pi(i,k)-pi(i-1,k))/dx
     1                             +fric
        enddo     
       enddo 
       endif
       do i=2,nx
        up(i,nz)=up(i,nz-1)
        up(i,1)=up(i,2)
       enddo

       if(iaccel.eq.1) call scalarbc(dudt)

c v eqn
       if(ienable_v.eq.1)then
       if(linear.eq.0)then
       do k=2,nz-1
          vp(2 ,k)=vm(2 ,k)
     1      -amin1((v(2 ,k)-cstar),0.)*(vm(3,k)-vm(2,k))*dtx
          vp(nx,k)=vm(nx,k)
     1      -amax1((v(nx,k)+cstar),0.)*(vm(nx,k)-vm(nx-1,k))*dtx
        do i=3,nx-1
         vp(i,k)=vm(i,k)-.25*dtx*((u(i+1,k)+u(i,k))*(v(i+1,k)+v(i,k))
     1                          -(u(i-1,k)+u(i,k))*(v(i-1,k)+v(i,k)))
     2                  -.25*dtz*(rhow(k+1)*(w(i,k+1)+w(i-1,k+1))
     2                                     *(v(i,k+1)+v(i,k))
     3                           -rhow( k )*(w(i,k  )+w(i-1,k  ))
     3                                   *(v(i,k-1)+v(i,k)))/rhou(k)
     6    +dkx*(d2t/dx**2)*(vm(i+1,k)-2.*vm(i,k)+vm(i-1,k))
c Coriolis force
          if(ienable_v.eq.1) vp(i,k)=vp(i,k)-d2t*corf*(u(i,k)-uinit(k))
          vp(i,k)=vp(i,k)
     1     +dkz*(d2t/dz**2)*(vm(i,k+1)-2.*vm(i,k)+vm(i,k-1))
        enddo     
       enddo 
       else ! linear
       do k=2,nz-1
          vp(2 ,k)=vm(2 ,k)
     1      -amin1((-1.*cstar),0.)*(vm(3,k)-vm(2,k))*dtx
          vp(nx,k)=vm(nx,k)
     1      -amax1(cstar,0.)*(vm(nx,k)-vm(nx-1,k))*dtx
        do i=3,nx-1
         vp(i,k)=vm(i,k)-.25*dtx*uinit(k)*(v(i+1,k)-v(i-1,k))
     6    +dkx*(d2t/dx**2)*(vm(i+1,k)-2.*vm(i,k)+vm(i-1,k))
c Coriolis force
          if(ienable_v.eq.1) vp(i,k)=vp(i,k)-d2t*corf*(u(i,k)-uinit(k))
          vp(i,k)=vp(i,k)
     1     +dkz*(d2t/dz**2)*(vm(i,k+1)-2.*vm(i,k)+vm(i,k-1))
        enddo     
       enddo 
       endif ! linear
       do i=2,nx
        vp(i,nz)=vp(i,nz-1)
        vp(i,1)=vp(i,2)
       enddo

       endif ! ienable_v
      
c w eqn
c w is zero at k=2 and nz
       do k=3,nz-1
       do i=2,nx-1
c needs u(nx)
c needs w(nx)=w(2); w(1)=w(nx-1)
        if(linear.eq.0)then
        wp(i,k)=wm(i,k)
     1    -.25*dtx*((u(i+1,k)+u(i+1,k-1))*(w(i+1,k)+w(i,k))
     2             -(u(i  ,k)+u(i  ,k-1))*(w(i-1,k)+w(i,k)))
     3    -.25*dtz*(rhou( k )*(w(i,k+1)+w(i,k))**2
     4             -rhou(k-1)*(w(i,k-1)+w(i,k))**2)/rhow(k)
     5    -0.5*dtz*cp*(tinit(k)+tinit(k-1))*(pi(i,k)-pi(i,k-1))
     7    +dtx*s_timfac*(strfcn(i+1,k)-strfcn(i,k))/rhow(k)
     8    +dkx*(d2t/dx**2)*(wm(i+1,k)-2.*wm(i,k)+wm(i-1,k))
     9    +dkz*(d2t/dz**2)*(wm(i,k+1)-2.*wm(i,k)+wm(i,k-1))
        if(imoist.eq.1)then
         wp(i,k)=wp(i,k)+d2t*g*0.5*
     1    (th(i, k )/tinit( k )+0.61*qv(i, k ) 
     2    +th(i,k-1)/tinit(k-1)+0.61*qv(i,k-1))
        else
         wp(i,k)=wp(i,k)+d2t*g*0.5*
     1    (th(i,k)/tinit(k) + th(i,k-1)/tinit(k-1))
        endif
        if(iaccel.eq.1) ! NOT mod for moisture
     1   dwdt(i,k)=-0.5*cp*(tinit(k)+tinit(k-1))*(pi(i,k)-pi(i,k-1))/dz
     1    +g*0.5*(th(i,k)/tinit(k) + th(i,k-1)/tinit(k-1))
         xlwdt(i,k)=dwdt(i,k)
     1    -.25*((u(i+1,k)+u(i+1,k-1))*(w(i+1,k)+w(i,k))
     2             -(u(i  ,k)+u(i  ,k-1))*(w(i-1,k)+w(i,k)))/dx
     3    -.25*(rhou( k )*(w(i,k+1)+w(i,k))**2
     4             -rhou(k-1)*(w(i,k-1)+w(i,k))**2)/rhow(k)/dz
       else ! linear
        wp(i,k)=wm(i,k)
     1    -.25*dtx*(uinit(k)+uinit(k-1))*(w(i+1,k)-w(i-1,k))
     5    -0.5*dtz*cp*(tinit(k)+tinit(k-1))*(pi(i,k)-pi(i,k-1))
     6    +d2t*g*0.5*(th(i,k)/tinit(k) + th(i,k-1)/tinit(k-1))
     7    +dtx*s_timfac*(strfcn(i+1,k)-strfcn(i,k))/rhow(k)
     8    +dkx*(d2t/dx**2)*(wm(i+1,k)-2.*wm(i,k)+wm(i-1,k))
     9    +dkz*(d2t/dz**2)*(wm(i,k+1)-2.*wm(i,k)+wm(i,k-1))
       endif ! linear
       enddo
       enddo
c boundary conditions on w
       call wbc(wp)
       if(iaccel.eq.1) call wbc(dwdt)
       if(iaccel.eq.1) call wbc(xlwdt)

c th eqn.  mean state vert advec can't be in flux form
c needs tinit(nz), tinit(1), but w is zero there
c needs th(1,k), th(nx,k)
       do k=2,nz-1
        ztmp=dz*(float(k)-1.5)
        do i=2,nx-1
         if(linear.eq.0)then
         if(tflux.eq.1)then
          thp(i,k)=thm(i,k)
     1     -0.5*dtx*(u(i+1,k)*(th(i+1,k)+th(i,k))
     2              -u( i ,k)*(th(i-1,k)+th(i,k)))
     3     -0.5*dtz*(rhow(k+1)*w(i,k+1)*(th(i,k+1)+th(i,k))
     4              -rhow( k )*w(i, k )*(th(i,k-1)+th(i,k)))/rhou(k)
     5     -0.5*dtz*(rhow(k+1)*w(i,k+1)*(tinit(k+1)-tinit( k ))
     6              +rhow( k )*w(i, k )*(tinit( k )-tinit(k-1)))/rhou(k) 
         else
c here t in advective form. Minus sign is used on lines 2,4 because
c   gradients are reversed (i-1)-(i) instead of opposite
          thp(i,k)=thm(i,k)
     1     -0.5*dtx*(u(i+1,k)*(th(i+1,k)-th(i,k))
     2              -u( i ,k)*(th(i-1,k)-th(i,k)))
     3     -0.5*dtz*(rhow(k+1)*w(i,k+1)*(th(i,k+1)-th(i,k))
     4              -rhow( k )*w(i, k )*(th(i,k-1)-th(i,k)))/rhou(k)
     5     -0.5*dtz*(rhow(k+1)*w(i,k+1)*(tinit(k+1)-tinit( k ))
     6              +rhow( k )*w(i, k )*(tinit( k )-tinit(k-1)))/rhou(k)
         endif
         else ! linear
          thp(i,k)=thm(i,k)
     1     -0.5*dtx*uinit(k)*(th(i+1,k)-th(i-1,k))
     5     -0.5*dtz*(rhow(k+1)*w(i,k+1)*(tinit(k+1)-tinit( k ))
     6              +rhow( k )*w(i, k )*(tinit( k )-tinit(k-1)))/rhou(k) 
         endif ! linear
c heat source(s)
         thp(i,k)=thp(i,k)+hsrc(i,k)*d2t*timfac
     1       +sb_hsrc(i,k)*d2t*sb_timfac
c cooling zone - if icoolzone = 1, it is maintained
         if(icoolzone.eq.1)then 
          if(i.ge.icooll.and.i.le.icoolr.and.ztmp.le.cz_depth)then
            thp(i,k)=thp(i,k)-d2t*(thm(i,k)+cz_ampl)*rcoolrate
          endif ! i
         endif ! icoolzone

c surface heat flux, if ishflux=1 and k=2.  a very crude convective velocity approach.
c   denom should be dz not dz/2.
         if(ishflux.eq.1.and.k.eq.2)then
          tdif=amax1(tground(i)-(thm(i,k)+tinit(k)),0.)
          avgu=0.5*abs(u(i+1,k)+u(i,k))
          avgu=amax1(avgu,5.0)
          wnetc=2.*sqrt(tdif) ! convec vel adj, from sbmodel
          vel=sqrt(avgu*avgu+wnetc*wnetc)
          thp(i,k)=thp(i,k)+d2t*cdh*vel*addflx(i)*tdif/dz
         endif
c diffusion         
         thp(i,k)=thp(i,k)
     1    +dkx*(d2t/dx**2)*(thm(i+1,k)-2.*thm(i,k)+thm(i-1,k))
     6    +dkz*(d2t/dz**2)*(thm(i,k+1)-2.*thm(i,k)+thm(i,k-1))
       enddo
       enddo
c boundary conditions
       call scalarbc(thp)


c qv eqn.  mean state vert advec can't be in flux form
      if(imoist.eq.1)then
       if(linear.eq.1) stop "turn off moisture"
       do k=2,nz-1
        ztmp=dz*(float(k)-1.5)
        do i=2,nx-1
          qvp(i,k)=qvm(i,k)
     1     -0.5*dtx*(u(i+1,k)*(qv(i+1,k)+qv(i,k))
     2              -u( i ,k)*(qv(i-1,k)+qv(i,k)))
     3     -0.5*dtz*(rhow(k+1)*w(i,k+1)*(qv(i,k+1)+qv(i,k))
     4              -rhow( k )*w(i, k )*(qv(i,k-1)+qv(i,k)))/rhou(k)
     5     -0.5*dtz*(rhow(k+1)*w(i,k+1)*(qinit(k+1)-qinit( k ))
     6              +rhow( k )*w(i, k )*(qinit( k )-qinit(k-1)))/rhou(k) 

c diffusion         
         qvp(i,k)=qvp(i,k)
     1    +dkx*(d2t/dx**2)*(qvm(i+1,k)-2.*qvm(i,k)+qvm(i-1,k))
     6    +dkz*(d2t/dz**2)*(qvm(i,k+1)-2.*qvm(i,k)+qvm(i,k-1))
       enddo
       enddo
c boundary conditions
       call scalarbc(qvp)
      endif ! imoist
c --------
c pi eqn
       do k=2,nz-1
        do i=2,nx-1
         pcof=csq/(rhou(k)*cp*tinit(k)*tinit(k))
         pip(i,k)=pim(i,k)
     1    -dtx*pcof*rhou(k)*tinit(k)*(u(i+1,k)-u(i,k))
     2    -0.5*dtz*pcof*(rhow(k+1)*w(i,k+1)*(tinit(k+1)+tinit(k))
     3                  -rhow( k )*w(i, k )*(tinit(k-1)+tinit(k)))
     4   +dkx*(d2t/dx**2)*(pim(i+1,k)-2.*pim(i,k)+pim(i-1,k))
     1   +dkz*(d2t/dz**2)*(pim(i,k+1)-2.*pim(i,k)+pim(i,k-1))
       enddo
       enddo
c boundary conditions
       call scalarbc(pip)
       if(ipressure.eq.1) call pressure

c time filter
       do i=1,nx
        do k=1,nz
         u (i,k)=u (i,k)+eps*(up (i,k)-2.*u (i,k)+um (i,k))
         if(ienable_v.eq.1)
     1   v (i,k)=v (i,k)+eps*(vp (i,k)-2.*v (i,k)+vm (i,k))
         w (i,k)=w (i,k)+eps*(wp (i,k)-2.*w (i,k)+wm (i,k)) 
         th(i,k)=th(i,k)+eps*(thp(i,k)-2.*th(i,k)+thm(i,k))
         if(imoist.eq.1)
     1   qv(i,k)=qv(i,k)+eps*(qvp(i,k)-2.*qv(i,k)+qvm(i,k)) 
         pi(i,k)=pi(i,k)+eps*(pip(i,k)-2.*pi(i,k)+pim(i,k))
        enddo
       enddo

c saturation adjustment, if moisture included

      if(imoist.eq.1) call satadj

c set for next time step
       do k=1,nz
        do i=1,nx
         um(i,k)=u(i,k)
         u(i,k)=up(i,k)
         up(i,k)=0.
         if(ienable_v.eq.1)then
          vm(i,k)=v(i,k)
          v(i,k)=vp(i,k)
          vp(i,k)=0.
         endif
         wm(i,k)=w(i,k)
         w(i,k)=wp(i,k)
         wp(i,k)=0.
         thm(i,k)=th(i,k)
         th(i,k)=thp(i,k)
         thp(i,k)=0.
         if(imoist.eq.1)then
          qvm(i,k)=qv(i,k)
          qv(i,k)=qvp(i,k)
          qvp(i,k)=0.
         endif
         pim(i,k)=pi(i,k)
         pi(i,k)=pip(i,k)
         pip(i,k)=0.
        enddo
       enddo 

c if maintained cooling zone is active, reconsider location for time > cz_coolrate/3 (was cz_coolrate/3)
      if(icoolzone.eq.1.and.time.gt.cz_coolrate/3)then
c  -- find gust front (east edge of cooling zone)
        icrit=5
        locgf_old=locgf
        locgf=9999
        do i=nx-1,2,-1
         if(th(i,2).le.-1.0)then
          locgf=i
          go to 101
         endif
        enddo
101    continue 
       if(locgf.ne.locgf_old)
     1  print *,' time ',time,' GF at i= ',locgf,' icoolr is ',icoolr,
     2 ' icooll is ',icooll
       if(locgf.lt.nx-10)then
        locdif=locgf-icoolr
        icoolr_old=icoolr
c        if(locgf.gt.icoolr)then
c         icoolr=icoolr+1
c        else if(locgf.lt.icoolr)then
c         icoolr=icoolr-1
c        endif
        if(locdif.gt.icrit)then
         icoolr=icoolr+1
        else if(locdif.lt.icrit)then
         icoolr=icoolr-1
        endif
        icooll=icoolr-icoolwd+1
        if(icoolr.ne.icoolr_old) 
     1 print *,' new cooling zone right edge at grid point ',icoolr,
     2 ' and icooll= ',icooll
       endif ! locgf
      endif ! icoolzone
      return
      end

c 
c subroutine anelastic
c
      subroutine anelastic
      include 'storage.txt'

      do i=1,nx
      do k=1,nz
       advu(i,k)=0.
       advw(i,k)=0.
       buoy(i,k)=0.
       byctrm(i,k)=0.
       dyntrm(i,k)=0.
       alltrm(i,k)=0.
      enddo
      enddo
c u eqn
       if(linear.ne.1)then
       do k=2,nz-1
          advu(2 ,k)=
     1      -amin1((u(2 ,k)-cstar),0.)*(um(3,k)-um(2,k))*dtx
          if(ienable_v.eq.1) advu(2,k)=advu(2,k)+corf*v(2,k)
          advu(nx,k)=
     1      -amax1((u(nx,k)+cstar),0.)*(um(nx,k)-um(nx-1,k))*dtx
          if(ienable_v.eq.1) advu(nx,k)=advu(nx,k)+corf*v(nx,k)
        do i=3,nx-1
         advu(i,k)=-.25*((u(i+1,k)+u(i,k))**2
     1                           -(u(i-1,k)+u(i,k))**2)/dx
     2                  -.25*(rhow(k+1)*(w(i,k+1)+w(i-1,k+1))
     2                                     *(u(i,k+1)+u(i,k))
     3                       -rhow( k )*(w(i,k  )+w(i-1,k  ))
     3                              *(u(i,k-1)+u(i,k)))/rhou(k)/dz
c     4                 -dtx*cp*tinit(k)*(pi(i,k)-pi(i-1,k))
     5    -s_timfac*(strfcn(i,k+1)-strfcn(i,k))/rhou(k)/dz
     6    +dkx*(1/dx**2)*(um(i+1,k)-2.*um(i,k)+um(i-1,k))
c          temp(i,k)=-dtz*(strfcn(i,k+1)-strfcn(i,k))/rhou(k)
c Coriolis
          if(ienable_v.eq.1) advu(i,k)=advu(i,k)+corf*v(i,k)
          advu(i,k)=advu(i,k)
     1     +dkz*(1/dz**2)*(um(i,k+1)-2.*um(i,k)+um(i,k-1)
     2                      -uinit(k+1)+2.*uinit(k)-uinit(k-1))
c experimental - add friction... not for coriolis yet
         if(fcoef.ne.0..and.k.eq.2.and.i.gt.icoast)then
           if(ienable_v.eq.1) stop 'fric-no_coriolis'
            fric=-fcoef*(um(i,k)-uinit(k))  ! linear friction on perturbations
         else
           fric=0.
         endif
         advu(i,k)=advu(i,k)+fric

c
        enddo     
       enddo 
       else ! linear model
        stop 'no linear'
c was here
       endif
c       do i=2,nx
c        advu(i,nz)=advu(i,nz-1)
c        advu(i,1)=advu(i,2)
c       enddo


       call scalarbc(advu)


c v eqn
       if(ienable_v.eq.1)then
        stop 'anelastic-not-set-up-for-V'
       if(linear.eq.0)then
       do k=2,nz-1
          vp(2 ,k)=vm(2 ,k)
     1      -amin1((v(2 ,k)-cstar),0.)*(vm(3,k)-vm(2,k))*dtx
          if(ienable_v.eq.1) vp(2,k)=vp(2,k)
     1 -d2t*corf*(u(2,k)-uinit(k))
          vp(nx,k)=vm(nx,k)
     1      -amax1((v(nx,k)+cstar),0.)*(vm(nx,k)-vm(nx-1,k))*dtx
          if(ienable_v.eq.1) vp(nx,k)=vp(nx,k)
     1 -d2t*corf*(u(nx,k)-uinit(k))
        do i=3,nx-1
         vp(i,k)=vm(i,k)-.25*dtx*((u(i+1,k)+u(i,k))*(v(i+1,k)+v(i,k))
     1                          -(u(i-1,k)+u(i,k))*(v(i-1,k)+v(i,k)))
     2                  -.25*dtz*(rhow(k+1)*(w(i,k+1)+w(i-1,k+1))
     2                                     *(v(i,k+1)+v(i,k))
     3                           -rhow( k )*(w(i,k  )+w(i-1,k  ))
     3                                   *(v(i,k-1)+v(i,k)))/rhou(k)
     6    +dkx*(d2t/dx**2)*(vm(i+1,k)-2.*vm(i,k)+vm(i-1,k))
c Coriolis force
          if(ienable_v.eq.1) vp(i,k)=vp(i,k)-d2t*corf*(u(i,k)-uinit(k))
          vp(i,k)=vp(i,k)
     1     +dkz*(d2t/dz**2)*(vm(i,k+1)-2.*vm(i,k)+vm(i,k-1))
        enddo     
       enddo 
       else ! linear
c was here
       endif ! linear
       do i=2,nx
        vp(i,nz)=vp(i,nz-1)
        vp(i,1)=vp(i,2)
       enddo

       endif ! ienable_v
      
c w eqn
c w is zero at k=2 and nz
       do k=3,nz-1
       do i=2,nx-1
c needs u(nx)
c needs w(nx)=w(2); w(1)=w(nx-1)
        if(linear.eq.0)then
        advw(i,k)=
     1    -.25*((u(i+1,k)+u(i+1,k-1))*(w(i+1,k)+w(i,k))
     2             -(u(i  ,k)+u(i  ,k-1))*(w(i-1,k)+w(i,k)))/dx
     3    -.25*(rhou( k )*(w(i,k+1)+w(i,k))**2
     4             -rhou(k-1)*(w(i,k-1)+w(i,k))**2)/rhow(k)/dz
c     5    -0.5*dtz*cp*(tinit(k)+tinit(k-1))*(pi(i,k)-pi(i,k-1))
     7    +s_timfac*(strfcn(i+1,k)-strfcn(i,k))/rhow(k)/dx
     8    +dkx*(1/dx**2)*(wm(i+1,k)-2.*wm(i,k)+wm(i-1,k))
     9    +dkz*(1/dz**2)*(wm(i,k+1)-2.*wm(i,k)+wm(i,k-1))
c        buoy(i,k)=g*0.5*(th(i,k)/tinit(k) + th(i,k-1)/tinit(k-1))        
        if(imoist.eq.1)then
         buoy(i,k)=g*0.5*
     1    (th(i, k )/tinit( k )+0.61*qv(i, k ) 
     2    +th(i,k-1)/tinit(k-1)+0.61*qv(i,k-1))
        else
         buoy(i,k)=g*0.5*
     1    (th(i,k)/tinit(k) + th(i,k-1)/tinit(k-1))
        endif

       else ! linear
c was

       endif ! linear
       enddo
       enddo
c boundary conditions on w
       call wbc(advw)
       call wbc(buoy)


c th eqn.  mean state vert advec can't be in flux form
c needs tinit(nz), tinit(1), but w is zero there
c needs th(1,k), th(nx,k)
       do k=2,nz-1
        ztmp=dz*(float(k)-1.5)
        do i=2,nx-1
         if(linear.eq.0)then
         if(tflux.eq.1)then
          thp(i,k)=thm(i,k)
     1     -0.5*dtx*(u(i+1,k)*(th(i+1,k)+th(i,k))
     2              -u( i ,k)*(th(i-1,k)+th(i,k)))
     3     -0.5*dtz*(rhow(k+1)*w(i,k+1)*(th(i,k+1)+th(i,k))
     4              -rhow( k )*w(i, k )*(th(i,k-1)+th(i,k)))/rhou(k)
     5     -0.5*dtz*(rhow(k+1)*w(i,k+1)*(tinit(k+1)-tinit( k ))
     6              +rhow( k )*w(i, k )*(tinit( k )-tinit(k-1)))/rhou(k) 
         else
c here t in advective form. Minus sign is used on lines 2,4 because
c   gradients are reversed (i-1)-(i) instead of opposite
          thp(i,k)=thm(i,k)
     1     -0.5*dtx*(u(i+1,k)*(th(i+1,k)-th(i,k))
     2              -u( i ,k)*(th(i-1,k)-th(i,k)))
     3     -0.5*dtz*(rhow(k+1)*w(i,k+1)*(th(i,k+1)-th(i,k))
     4              -rhow( k )*w(i, k )*(th(i,k-1)-th(i,k)))/rhou(k)
     5     -0.5*dtz*(rhow(k+1)*w(i,k+1)*(tinit(k+1)-tinit( k ))
     6              +rhow( k )*w(i, k )*(tinit( k )-tinit(k-1)))/rhou(k)
         endif
         else ! linear
          thp(i,k)=thm(i,k)
     1     -0.5*dtx*uinit(k)*(th(i+1,k)-th(i-1,k))
     5     -0.5*dtz*(rhow(k+1)*w(i,k+1)*(tinit(k+1)-tinit( k ))
     6              +rhow( k )*w(i, k )*(tinit( k )-tinit(k-1)))/rhou(k) 
         endif ! linear
c heat source(s)
         thp(i,k)=thp(i,k)+hsrc(i,k)*d2t*timfac
     1       +sb_hsrc(i,k)*d2t*sb_timfac
c cooling zone - if icoolzone = 1, it is maintained
         if(icoolzone.eq.1)then 
          if(i.ge.icooll.and.i.le.icoolr.and.ztmp.le.cz_depth)then
            thp(i,k)=thp(i,k)-d2t*(thm(i,k)+cz_ampl)*rcoolrate
          endif ! i
         endif ! icoolzone

c surface heat flux, if ishflux=1 and k=2.  a very crude convective velocity approach.
c   denom should be dz not dz/2.
         if(ishflux.eq.1.and.k.eq.2)then
          tdif=amax1(tground(i)-(thm(i,k)+tinit(k)),0.)
          avgu=0.5*abs(u(i+1,k)+u(i,k))
          avgu=amax1(avgu,5.0)
          wnetc=2.*sqrt(tdif) ! convec vel adj, from sbmodel
          vel=sqrt(avgu*avgu+wnetc*wnetc)
          thp(i,k)=thp(i,k)+d2t*cdh*vel*addflx(i)*tdif/dz
         endif
c diffusion         
         thp(i,k)=thp(i,k)
     1    +dkx*(d2t/dx**2)*(thm(i+1,k)-2.*thm(i,k)+thm(i-1,k))
     6    +dkz*(d2t/dz**2)*(thm(i,k+1)-2.*thm(i,k)+thm(i,k-1))
       enddo
       enddo
c boundary conditions
       call scalarbc(thp)

c qv eqn.  mean state vert advec can't be in flux form
      if(imoist.eq.1)then
       if(linear.eq.1) stop "turn off moisture"
       do k=2,nz-1
        ztmp=dz*(float(k)-1.5)
        do i=2,nx-1
          qvp(i,k)=qvm(i,k)
     1     -0.5*dtx*(u(i+1,k)*(qv(i+1,k)+qv(i,k))
     2              -u( i ,k)*(qv(i-1,k)+qv(i,k)))
     3     -0.5*dtz*(rhow(k+1)*w(i,k+1)*(qv(i,k+1)+qv(i,k))
     4              -rhow( k )*w(i, k )*(qv(i,k-1)+qv(i,k)))/rhou(k)
     5     -0.5*dtz*(rhow(k+1)*w(i,k+1)*(qinit(k+1)-qinit( k ))
     6              +rhow( k )*w(i, k )*(qinit( k )-qinit(k-1)))/rhou(k) 

c diffusion         
         qvp(i,k)=qvp(i,k)
     1    +dkx*(d2t/dx**2)*(qvm(i+1,k)-2.*qvm(i,k)+qvm(i-1,k))
     6    +dkz*(d2t/dz**2)*(qvm(i,k+1)-2.*qvm(i,k)+qvm(i,k-1))
       enddo
       enddo
c boundary conditions
       call scalarbc(qvp)
      endif ! imoist


c blktri
      do k=2,nz-1
       a1=1./(rhou(k)*dz)
       a2=1./(cp*tinit(k))
       do i=2,nx-1
        byctrm(i-1,k-1)=(buoy(i,k+1)*rhow(k+1)-buoy(i,k)*rhow(k))*a1*a2
        dyntrm(i-1,k-1)=((advu(i+1,k)-advu(i,k))/dx
     1               +(advw(i,k+1)*rhow(k+1)-advw(i,k)*rhow(k))*a1)*a2
        alltrm(i-1,k-1)=dyntrm(i-1,k-1)+byctrm(i-1,k-1)
     2              +((um(i+1,k)-um(i,k))/dx
     3                +(wm(i,k+1)*rhow(k+1)-wm(i,k)*rhow(k))*a1)*a2/d2t

       enddo
      enddo

      call blktri(1,1,nz-2,an,bn,cn,1,nx-2,am,bm,cm,nxm,alltrm,ier,wa)
      if(ier.ne.0) then
       print *,' anelastic routine ier ',ier
       stop 'blktri'
      endif
c      print *,' pi at corner ',alltrm(nx-2,1)
      call blktri(1,1,nz-2,an,bn,cn,1,nx-2,am,bm,cm,nxm,byctrm,ier,wa)
      if(ier.ne.0) then
       print *,' pbyc ier ',ier
       stop 'blktri'
      endif
      call blktri(1,1,nz-2,an,bn,cn,1,nx-2,am,bm,cm,nxm,dyntrm,ier,wa)
      if(ier.ne.0) then
       print *,' pdyn ier ',ier
       stop 'blktri'
      endif
      do k=2,nz-1
       do i=2,nx-1
        pi(i,k)=alltrm(i-1,k-1)-alltrm(1,1)
        pbyc(i,k)=byctrm(i-1,k-1)-byctrm(1,1)
        pdyn(i,k)=dyntrm(i-1,k-1)-dyntrm(1,1)
        ptot(i,k)=pbyc(i,k)+pdyn(i,k)
       enddo
      enddo
      call scalarbc(pi)
      call scalarbc(pbyc)
      call scalarbc(pdyn)   
      call scalarbc(ptot)  

      do k=2,nz-1
         do i=2,nx-1
            dudtd(i,k)=-cp*tinit(k)*(pdyn(i,k)-pdyn(i-1,k))/dx
c friction added to dynamic accel
            if(fcoef.ne.0..and.k.eq.2.and.i.gt.icoast)then
             dudtd(i,k)=dudtd(i,k)-fcoef*(um(i,k)-uinit(k))  ! linear friction on perturbations
            endif

            dudtb(i,k)=-cp*tinit(k)*(pbyc(i,k)-pbyc(i-1,k))/dx
            dudtt(i,k)=dudtb(i,k)+dudtd(i,k)
c            if(i.eq.50) print *,' k ',k,dudtd(i,k),dudtb(i,k),dudtt(i,k),pdyn(i,k)
         enddo
      enddo
      call scalarbc(dudtb)
      call scalarbc(dudtd)
      call scalarbc(dudtt)
      do i=2,nx-1
          do k=3,nz-1
            dwdtd(i,k)=-0.5*cp*(tinit(k)+tinit(k-1))*
     1         (pdyn(i,k)-pdyn(i,k-1))/dz
            dwdtb(i,k)=-0.5*cp*(tinit(k)+tinit(k-1))*
     1         (pbyc(i,k)-pbyc(i,k-1))/dz
     1         +g*0.5*(th(i,k)/tinit(k) + th(i,k-1)/tinit(k-1))
            dwdtt(i,k)=dwdtd(i,k)+dwdtb(i,k)
          enddo
      enddo
      call wbc(dwdtb)
      call wbc(dwdtd)
      call wbc(dwdtt)



      do k=2,nz-1
       do i=2,nx ! pgf does not affect i=2 or i=nx owing to BCs applied to pi
        up(i,k)=um(i,k)+advu(i,k)*d2t 
     1    -dtx*cp*tinit(k)*(pi(i,k)-pi(i-1,k))
       enddo
      enddo
      do k=3,nz-1
       do i=2,nx-1
        wp(i,k)=wm(i,k)+(advw(i,k)+buoy(i,k))*d2t
     1     -0.5*dtz*cp*(tinit(k)+tinit(k-1))*(pi(i,k)-pi(i,k-1))
       enddo
      enddo
       do i=2,nx
        up(i,nz)=up(i,nz-1)
        up(i,1)=up(i,2)
       enddo
       call wbc(wp)

c pi eqn
c       do k=2,nz-1
c        do i=2,nx-1
c         pcof=csq/(rhou(k)*cp*tinit(k)*tinit(k))
c         pip(i,k)=pim(i,k)
c     1    -dtx*pcof*rhou(k)*tinit(k)*(u(i+1,k)-u(i,k))
c     2    -0.5*dtz*pcof*(rhow(k+1)*w(i,k+1)*(tinit(k+1)+tinit(k))
c     3                  -rhow( k )*w(i, k )*(tinit(k-1)+tinit(k)))
c     4   +dkx*(d2t/dx**2)*(pim(i+1,k)-2.*pim(i,k)+pim(i-1,k))
c     1   +dkz*(d2t/dz**2)*(pim(i,k+1)-2.*pim(i,k)+pim(i,k-1))
c       enddo
c       enddo
c boundary conditions
c       call scalarbc(pip)
c       if(ipressure.eq.1) call pressure

c time filter
       do i=1,nx
        do k=1,nz
         u (i,k)=u (i,k)+eps*(up (i,k)-2.*u (i,k)+um (i,k))
         if(ienable_v.eq.1)
     1   v (i,k)=v (i,k)+eps*(vp (i,k)-2.*v (i,k)+vm (i,k))
         w (i,k)=w (i,k)+eps*(wp (i,k)-2.*w (i,k)+wm (i,k)) 
         th(i,k)=th(i,k)+eps*(thp(i,k)-2.*th(i,k)+thm(i,k))
         if(imoist.eq.1)
     1   qv(i,k)=qv(i,k)+eps*(qvp(i,k)-2.*qv(i,k)+qvm(i,k)) 
         pi(i,k)=pi(i,k)+eps*(pip(i,k)-2.*pi(i,k)+pim(i,k))
        enddo
       enddo

c saturation adjustment, if moisture included

      if(imoist.eq.1) call satadj

c set for next time step
       do k=1,nz
        do i=1,nx
         um(i,k)=u(i,k)
         u(i,k)=up(i,k)
         up(i,k)=0.
         if(ienable_v.eq.1)then
          vm(i,k)=v(i,k)
          v(i,k)=vp(i,k)
          vp(i,k)=0.
         endif
         wm(i,k)=w(i,k)
         w(i,k)=wp(i,k)
         wp(i,k)=0.
         thm(i,k)=th(i,k)
         th(i,k)=thp(i,k)
         thp(i,k)=0.
c         pim(i,k)=pi(i,k)
c         pi(i,k)=pip(i,k)
c         pip(i,k)=0.
        enddo
       enddo 

c if maintained cooling zone is active, reconsider location for time > cz_coolrate/3 (was cz_coolrate/3)
      if(icoolzone.eq.1.and.time.gt.cz_coolrate/3)then
c  -- find gust front (east edge of cooling zone)
        icrit=5
        locgf_old=locgf
        locgf=9999
        do i=nx-1,2,-1
         if(th(i,2).le.-1.0)then
          locgf=i
          go to 101
         endif
        enddo
101    continue 
       if(locgf.ne.locgf_old)
     1  print *,' time ',time,' GF at i= ',locgf,' icoolr is ',icoolr,
     2 ' icooll is ',icooll
       if(locgf.lt.nx-10)then
        locdif=locgf-icoolr
        icoolr_old=icoolr
c        if(locgf.gt.icoolr)then
c         icoolr=icoolr+1
c        else if(locgf.lt.icoolr)then
c         icoolr=icoolr-1
c        endif
        if(locdif.gt.icrit)then
         icoolr=icoolr+1
        else if(locdif.lt.icrit)then
         icoolr=icoolr-1
        endif
        icooll=icoolr-icoolwd+1
        if(icoolr.ne.icoolr_old) 
     1 print *,' new cooling zone right edge at grid point ',icoolr,
     2 ' and icooll= ',icooll
       endif ! locgf
      endif ! icoolzone
      return
      end

c 
c subroutine sounding
c
      subroutine sounding
      include 'storage.txt'
c zero arrays to start
      do k=1,nz
       qinit(k)=0.
       uinit(k)=0.
       vinit(k)=0.
       qinit(k)=0.
       zero (k)=0.
      enddo
      do i=1,nx
       do k=1,nz
        u(i,k)=0.
        um(i,k)=0.
        w(i,k)=0.
        wm(i,k)=0.
        th(i,k)=0.
        thm(i,k)=0.
	qv(i,k)=0.
	qvm(i,k)=0.
        pi(i,k)=0.
        pim(i,k)=0.
        strfcn(i,k)=0.
       enddo
      enddo

c set up base state horizontal wind

      uinit(2)=usurf
      uinit(1)=uinit(2)
      zusfc=0.5*dz        ! height of lowest U point above model surface
      do k=3,nz-1
       zu=(k-1.5)*dz
       if(zu.le.sdepth1)then
         uinit(k)=usurf+shear1*(zu-zusfc)
       else if(zu.le.sdepth1+sdepth2)then
         uinit(k)=usurf+shear1*sdepth1+shear2*(zu-sdepth1)
       else
         uinit(k)=usurf+shear1*sdepth1+shear2*sdepth2
     1              +shear3*(zu-sdepth2)
       endif
       print *,' zu ',zu,' u ',uinit(k)
      enddo
      uinit(nz)=uinit(nz-1)


c set up potential temperature of base state
      surface_theta=thsurf ! 300.
      tinit(2)=surface_theta
      tinit(1)=tinit(2)
      do k=3,nz-1
       zu=(k-1.5)*dz
       if(zu.lt.pbld)then           ! PBL
        tinit(k)=alog(tinit(k-1))+xnsq_pbl*dz/g
       else if(zu.lt.tropo)then     ! troposphere
        tinit(k)=alog(tinit(k-1))+xnsq_trop*dz/g
       else                         ! stratosphere
        tinit(k)=alog(tinit(k-1))+xnsq_strat*dz/g
       endif
       tinit(k)=exp(tinit(k))
      enddo
      tinit(nz)=tinit(nz-1)

c set up initial specific humidity, if imoist=1
      if(imoist.eq.1)then
       qsurf  =rhsurf*0.01   ! initial values are RH in pct
       qvalue1=rhvalue1*0.01
       qvalue2=rhvalue2*0.01
       qvalue3=rhvalue3*0.01
       qheight1=rhheight1
       qheight2=rhheight2
       qheight3=rhheight3
       qinit(2)=qsurf
       qinit(1)=qinit(2)
       qrate_1=(qvalue1-qsurf)/qheight1 ! lapse rates
       qrate_2=(qvalue2-qvalue1)/(qheight2-qheight1)
       qrate_3=(qvalue3-qvalue2)/(qheight3-qheight2)
       do k=3,nz-1
        zu=(k-1.5)*dz
        if(zu.lt.qheight1)then
         qinit(k)=qsurf+qrate_1*zu
        else if(zu.lt.qheight2)then
         qinit(k)=qvalue1+qrate_2*(zu-qheight1)
        else if(zu.lt.qheight3)then
         qinit(k)=qvalue2+qrate_3*(zu-qheight2)
        else ! above qheight4
         qinit(k)=qvalue3
        endif ! zu
       enddo
       qinit(nz)=qinit(nz-1)

        do k=1,nz
         rhtemp(k)=qinit(k)
        enddo

      endif ! imoist
      do k=2,nz-1
       zu=(k-1.5)*dz
       print *,' z ',zu,' tinit ',tinit(k),' uinit ',uinit(k),
     1 ' qinit ',qinit(k)
      enddo

c iterate to get mean pressure, density.... if moist=1
      maxiter=0
      if(imoist.eq.1) maxiter=1
      do iter=0,maxiter
       print *,'iter ',iter, ' maxiter ',maxiter
       xiter=float(iter)
c PIK is mean pi at theta locations - SKIPPING QINIT since it is RH - "here"
       pk(2)=(psurf/psl)**xk
     1            -g*dz/(2.*tinit(2)*(1.+.61*qinit(2)*xiter)*cp) 
       do k=3,nz-1
        pk(k)=pk(k-1)-2.*g*dz/((tinit(k-1)*(1.+.61*qinit(k-1)*xiter)       
     1                    +tinit( k )*(1.+.61*qinit( k )*xiter))*cp) 
       enddo
c rhou is mean density at U heights
      do k=2,nz-1
       den=pk(k)**cvrd
       rhou(k)=den*psl/(rd*tinit(k)*(1.+.61*qinit(k)*xiter)) ! added q  ! here
      enddo
c fake points for rhou
      rhou(1)=rhou(2)
      rhou(nz)=rhou(nz-1)
c rhow is density at W heights
c rhow at surface:
      sfcpi=(psurf/psl)**xk
      rhow(2)=(sfcpi**cvrd)*psl/(rd*tinit(2)*(1.+.61*qinit(2)*xiter)) ! added q  ! here
      rhow(1)=rhow(2)
      do k=3,nz
        rhow(k)=0.5*(rhou(k-1)+rhou(k))
      enddo
c write out base state
      do k=2,nz
       zw=(k-2)*dz
       zu=(k-1.5)*dz
       print *,' k',k,' zw',zw,' zu ',zu,' rou ',rhou(k),' row ',
     1  rhow(k),' pib ',pk(k)
      enddo
      if(imoist.eq.1)then
      do k=2,nz-1
       zu=(k-1.5)*dz
       p0=(pk(k)**xki)*psl ! two ways of calc p, same to w/in roundoff
       p01=rhou(k)*rd*tinit(k)*(1.+0.61*qinit(k)*xiter)*pk(k)  ! here
       qvs=(380./p0)*exp(17.27*(tinit(k)*pk(k)-273.)
     1                        /(tinit(k)*pk(k)-36.))
c       rh=qinit(k)/qvs
       
       rh=qinit(k)
       qinit(k)=rhtemp(k)*qvs  ! change RH to vapor mixing ratio
c OVERWRITE MOISTURE INPUT if code activated =======
       IF(1.GT.2)THEN ! ******************
       if(zu.lt.12000.)then
        rh=1.-.8*(zu/12000.)**1.25  ! was mult by .75
       else
        rh=0.25
       endif
       qinit(k)=rh*qvs
       qinit(k)=amin1(qinit(k),0.012)

       ENDIF ! 1 vs 2
       print *,' z ',zu,p0,qinit(k),' qvs ',qvs,' rh ',rhtemp(k)
      enddo
      qinit(1)=qinit(2)
      qinit(nz)=qinit(nz-1)
      endif ! imoist
      enddo ! iter
c populate U arrays with base state wind
      do i=1,nx
      do k=1,nz
       u(i,k)=uinit(k)
       um(i,k)=uinit(k)
      enddo
      enddo
      return
      end

c
c subroutine zeroit
c
      subroutine zeroit
      include 'storage.txt'
      ienable_v=0
      corf=0.
      do i=1,nxm
       do k=1,nzm
        hsrc(i,k)=0.
        sb_hsrc(i,k)=0.
        advu(i,k)=0.
        advw(i,k)=0.
        buoy(i,k)=0.
        byctrm(i,k)=0.
        dyntrm(i,k)=0.
        alltrm(i,k)=0.
        pbyc(i,k)=0.
        pdyn(i,k)=0.
        ptot(i,k)=0.
       enddo
      enddo
      return
      end

c
c subroutine thermal
c
      subroutine thermal
      include 'storage.txt'

c initial thermal
      imid=(nx+1)/2
c      imid=120
      do k=2,nz-1
      do i=2,nx-1
       argz=((dz*(float(k)-1.5)-zcnt)/radz)**2
       argx=(dx*(i-imid)/radx)**2
       rad=sqrt(argz+argx)
       if(rad.le.1.)then
        th(i,k)=0.5*delt*(cos(trigpi*rad)+1)
        print *,' ik ',i,k,' rad ',rad,' th ',th(i,k)
       else
        th(i,k)=0.
       endif
c test - inf. wide
c      if(k.lt.7) th(i,k)=delt
       thm(i,k)=th(i,k)
      enddo
      enddo
c force lateral boundaries, tho should't be necessary here
      call scalarbc(th)
      call scalarbc(thm)
      return
      end

c
c subroutine surface
c
      subroutine surface
      include 'storage.txt'
      do i=1,nx
       if(i.ge.icoast)then
        tground(i)=tinit(2)+tdelt
       else
        tground(i)=tinit(2)
       endif
      enddo
c random perts added to surface heat flux?
      if(irand.eq.1)then
       call randnum
      else
       do i=1,nx
        addflx(i)=1.0
       enddo
      endif
      return
      end

c
c subroutine streamfcn
c

      subroutine streamfcn
      include 'storage.txt'
      print *,' s_ampl ',s_ampl
      print *,' psik ',s_psik,' psim ',s_psim
      imid=(nx+1)/2
c      imid=(nx+1)/2+100
       if(s_repeat.eq.1)then
         argu1_left=-999
         argu1_right=999
c         argu1_left=-20*trigpi ! -999
c         argu1_right=999 ! 10*trigpi ! 999
       else
         argu1_left=-trigpi    ! considered .25 
         argu1_right=trigpi
       endif
      print *,' argu1_left= ', argu1_left
      print *,' argu1_right= ', argu1_right
      do i=2,nx-1
       do k=2,nz-1
         right_taper=1.0 ! amin1(float(nx-i)/120,1.0)
         argu1=s_psik*dx*(float(i-imid)-0.5)
         argu2=s_psim*(dz*float(k-1)-s_znaught)
         argu=argu1+argu2
         if(argu1.gt.argu1_left.and.argu1.lt.argu1_right.and.
     1     argu2.gt.-0.5*trigpi.and.argu2.lt.0.5*trigpi)then
c         if(argu2.gt.-0.5*trigpi.and.argu2.lt.0.5*trigpi)then
           strfcn(i,k)=s_ampl*sin(argu1)*cos(argu2)
     1      +0.0*sin(1.2*argu1)*cos(argu2)             ! off
           strfcn(i,k)=strfcn(i,k)*right_taper
           if(k.eq.28) print *,' i ',i,' str ',strfcn(i,k)
         else
           strfcn(i,k)=0.
c         endif ! argu2
         endif ! argu1
       enddo
      enddo
      return
      end

c
c subroutine heatsrc
c

      subroutine heatsrc
      include 'storage.txt'
      print *,' called heatsrc'
      imid=(nx+1)/2
      xmode2=0.
      
      if(h_modes.eq.2) xmode2=1.0
      do k=2,nz-1
       argz=((dz*(float(k)-0.0)-h_center_z)/h_radius_z)**2 ! was -1.5
       vert=1.0*(trigpi/h_radius_z)*(dz*(float(k)-0.0)) ! was 0.5*, was -1.5
       do i=2,nx-1
        argx=(dx*(i-imid)/h_radius_x)**2
        rad=sqrt(argz+argx)
        if(rad.le.1.)then
         hsrc(i,k)=0.5*h_ampl*(cos(trigpi*rad)+1)
     1                        *(sin(vert)-0.5*xmode2*sin(2*vert))
c         hsrc(i,k)=amax1(hsrc(i,k),0.)
c         if(hsrc(i,k).lt.h_ampl/10.) hsrc(i,k)=0.
         if(i.eq.imid) print *,' k ',k,' hsrc ',hsrc(i,k),' z ',
     1 dz*float(k)-1.5,' vert ',vert, ' rad ',rad, ' argz ',argz
        endif
       enddo
      enddo
      return
      end

c
c subroutine seabreeze
c

      subroutine seabreeze
      include 'storage.txt'
      print *,' called seabreeze'
      if(sb_latitude.ne.0)then
       ienable_v=1
       corf=2.*7.292e-5*sin(sb_latitude*trigpi/180.)
       print *,' coriolis factor ',corf
      endif
      imid=(nx+1)/2
      tp2=0.5*trigpi
      if(sb_period.ne.0)then
       sb_freq=(2.*trigpi/86400.)/sb_period ! 7.292e-5/sb_period
      else
       stop ' sb_period cannot be 0'
      endif
      do k=2,nz-1
       zloc=dz*(float(k)-1.5)
       zterm=exp(-1.*zloc/sb_z0)
       do i=2,nx-1
        xterm=float(i-imid)*dx/sb_x0
        sb_hsrc(i,k)=sb_ampl*(tp2+atan(xterm))*zterm
        sb_hsrc(i,k)=amax1(sb_hsrc(i,k),0.)
c        print *,' ik ',i,k,' sbhsrc ',sb_hsrc(i,k),sb_ampl
       enddo
      enddo
      return
      end

c
c subroutine coolzone
c

      subroutine coolzone
      include 'storage.txt'
      locgf=9999
      if(cz_coolrate.ne.0)then
       rcoolrate=1./cz_coolrate
      else if(icoolzone.ne.2)then
       stop ' cooling rate is zero with icoolzone=1.  cannot continue.'
      endif
c right edge of cooling zone
      xtmp=cz_rightedge*float(nx)
      icoolr=ifix(xtmp)
c width determines left edge
      icoolwd=cz_width/dx
      print *, ' cooling zone width is ',icoolwd,' grid points '
      icooll=icoolr-icoolwd+1
      icooll=max(icooll,2)
      if(icoolzone.eq.2)then
       do k=2,nz-1
        ztmp=dz*(float(k)-1.5)
        if(ztmp.le.cz_depth)then
         do i=icooll,icoolr
          th(i,k)=-1.*cz_ampl
          thm(i,k)=th(i,k)
         enddo
        endif
       enddo ! k
       call scalarbc(th)
       call scalarbc(thm)
      endif ! icoolzone = 2
        
      return
      end

c
c subroutine randnum
c
      subroutine randnum
      include 'storage.txt'
      real randx,randrng

c the following generates an array addflx dimensioned nx x ny which
c contains a set of values over the range 1.0 plus/minus randrng
c such that when the array addflx at point (x,y) is multiplied by
c the heat flux, it will be increased or decreased by a random value
c from 0 to ranrng percent
c note: random flux only takes place over land source (= 1.0 over sea)

c (iseed is the pseudorandom number generator seed)
      iseed=12357.

      randrng=0.05
      idiffx=0

      do irnd=2,nx-1
          iseed=2045*iseed+1
          iseed=iseed-(iseed/1048576)*1048576
          randx=float(iseed+1)/1048577.0

          if (irnd.ge.icoast) then
            addflx(irnd)=1.0+(randrng*2.*(randx-0.5))
          else
            addflx(irnd)=1.0
          end if


        print *,' i ',irnd,' addflx ',addflx(irnd)
14099   format(10e13.6)
14046   continue

      enddo

      return
      end


c
c subroutine readnl - reads namelist
c
      subroutine readnl
      include 'storage.txt'

c - local declarations section -
      character*80 casename
      character*4 ctl,dat


c namelist declarations
      namelist /experiment/ casename
      namelist /grid_run/ nx,nz,dx,dz,dt,timend,plot
      namelist /framework/ ipressure,ianelastic,csnd,imoist,
     1 fcoef
      namelist /numerics/ byteswap,cstar,dkx,dkz,eps

      namelist /environ/ thsurf,bvpbl,pbld,bvtropo,tropo,
     1 bvstrat,psurf,usurf,shear1,sdepth1,shear2,sdepth2,shear3,
     1 rhsurf,rhvalue1,rhheight1,rhvalue2,rhheight2,rhvalue3,
     1 rhheight3
      namelist /thermal/ ithermal,delt,radx,radz,zcnt
      namelist /surface_flux/ ishflux,tdelt,icoast,cdh,irand
      namelist /streamfunction/ istrfcn,s_repeat,s_ampl,
     1 s_znaught,s_hwavel,s_vwavel,s_period
      namelist /atmos_heat_source/ ihsrc,h_ampl,h_radius_x,
     1 h_radius_z,h_center_z,h_freq,h_modes
      namelist /rotunno_seabreeze/ iseabreeze,sb_ampl,
     1 sb_x0,sb_z0,sb_period,sb_latitude,sb_linear
      namelist /cooling_zone/ icoolzone,cz_ampl,cz_rightedge,
     1 cz_depth,cz_width,cz_coolrate

c - end local declarations -

c
c executable code starts here
c
      trigpi=4.0*atan(1.0)       ! trigonometric pi


c-----------------------------------------------------------------
c declare default values for all items in the namelist
c-----------------------------------------------------------------
c
c - experiment -
      casename = 'testcase'
c - grid_run
      nx = 101
      nz = 84
      dx = 1000.
      dz = 250.
      dt = 1.0
      timend = 7200.
      plot = 300.
c - framework -
      ipressure = 0
      ianelastic = 0
      csnd = 50.
      imoist = 0
       fcoef = 0
c - numerics -
      byteswap = 0
      cstar = 50.
      dkx = 500.
      dkz = 100.
      eps = 0.005
c - environ -
        thsurf=300.
      	bvpbl = 0.01
	pbld = 2000.
	bvtropo = 0.01
	tropo = 12000.
	bvstrat = 0.02
	psurf = 965.
	usurf = 0.
	shear1 = 0.
	sdepth1 = 3000.
	shear2 = 0.
	sdepth2 = 1500.
	shear3 = 0.
        rhsurf = 0.
        rhvalue1 = 0.
        rhheight = 2000.
        rhvalue2 = 0.
        rhheight2 = 8000.
        rhvalue3 = 0.
        rhheight3 = 0.
c - thermal -
	ithermal = 0
	delt = 3.
	radx = 8000.
	radz = 4000.
	zcnt = 3000.
c - surface_flux -
	ishflux = 0
	tdelt = 12.
	icoast = 30
	cdh = 7.2e-3
	irand = 0
c - streamfunction -
	istrfcn = 0
	s_repeat = 0
	s_ampl = 40.
	s_znaught = 6000.
	s_hwavel = 40000.
	s_vwavel = 18000.
	s_period = 1200.
c - atmos_heat-source -
	ihsrc = 0
	h_ampl = 0.075
	h_radius_x = 3000.
	h_radius_z = 3000.
	h_center_z = 3000.
	h_freq = 0.005
	h_modes = 2
c - rotunno_seabreeze -
	iseabreeze = 1
	sb_ampl = 0.000175
	sb_x0 = 1000.
	sb_z0 = 1000.
	sb_period = 1.0
	sb_latitude = 60.
	sb_linear = 1
c - cooling_zone -
	icoolzone = 0
	cz_ampl = 4.5
	cz_rightedge = 0.60
	cz_depth = 2500.
	cz_width = 46000.
	cz_coolrate = 600.


c
c read in the namelist
c

      read(5,experiment,end=100)

      knx=index(casename,' ')
      write(ctl,'(a4)') '.ctl'
      write(dat,'(a4)') '.dat'
      ctlfile=casename(1:knx-1)//ctl
      datfile=casename(1:knx-1)//dat
      write(6,*) 'casename is ',casename
      write(6,*) 'ctlfile is ',ctlfile
      write(6,*) 'datfile is ',datfile

      read(5,grid_run,end=100)
      write(6,'(a)')'Namelist grid_run sucessfully read.'

      read(5,framework,end=100)
      write(6,'(a)')'Namelist framework sucessfully read.'

      read(5,numerics,end=100)
      write(6,'(a)')'Namelist numerics sucessfully read.'

      read(5,environ,end=100)
      write(6,'(a)')'Namelist environ sucessfully read.'

      read(5,thermal,end=100)
      write(6,'(a)')'Namelist thermal sucessfully read.'

      read(5,surface_flux,end=100)
      write(6,'(a)')'Namelist surface_flux sucessfully read.'

      read(5,streamfunction,end=100)
      write(6,'(a)')'Namelist streamfunction sucessfully read.'

      read(5,atmos_heat_source,end=100)
      write(6,'(a)')'Namelist atmos_heat_source sucessfully read.'

      read(5,rotunno_seabreeze,end=100)
      write(6,'(a)')'Namelist rotunno_seabreeze sucessfully read.'

      read(5,cooling_zone,end=100)
      write(6,'(a)')'Namelist cooling_zone sucessfully read.'

c
c print info and perform manipulations based on parameter choices
c

      print *,' ========================================= '
      print *,' input parameters to DTDM model '
      print *,' ========================================= '
      print *,' BYTESWAP = ',byteswap,' =1 if byteswapping needed'
      print *,' NX = ',nx,' number of horizontal grid points'
      print *,' NZ = ',nz,' number of vertical grid points'
      print *,' DX = ',dx,' horizontal grid spacing (m)'
      print *,' DZ = ',dz,' vertical grid spacing (m)'
      print *,' DT = ',dt,' time step (s)'
      print *,' TIMEND = ',timend,' model runtime (s)'
      print *,' PLOT = ',plot,' plotting interval (s)'
      print *,' IPRESSURE = ',ipressure,' output pressure decomps '
      print *,' IANELASTIC = ',ianelastic,' anelastic=1 compressible=0'
      print *,' CSND = ',csnd,' sound speed (m/s)'
      print *,' IMOIST = ',imoist,' moist run = 1, dry = 0'
      print *,'  FCOEF = ',fcoef,' surface linear friction coef'
      print *,' CSTAR = ',cstar,' gravity wave speed for boundary (m/s)'
      print *,' DKX = ',dkx,' horizontal diffusion (m*m/s)'
      print *,' DKZ = ',dkz,' vertical diffusion (m*m/s)'
      print *,' EPS = ',eps,' leapfrog time filter (ndim)'
      print *,' THSURF = ',thsurf,' sfc potential temperature (K)'
      print *,' BVPBL = ',bvpbl,' PBL BV freq (1/s)'
      print *,' PBLD = ',pbld,' PBL depth (m)'
      print *,' BVTROPO = ',bvtropo,' free tropospheric BV freq (1/s)'
      print *,' TROPO = ',tropo,' tropopause height (m)'
      print *,' BVSTRAT = ',bvstrat,' stratospheric BV freq (1/s)'
      print *,' PSURF = ',psurf,' surface pressure (mb)'
      print *,' USURF = ',usurf,' surface wind speed (m/s)'
      print *,'   1st layer shear ',shear1,' (1/s); depth ',sdepth1,' m'
      print *,'   2nd layer shear ',shear2,' (1/s); depth ',sdepth2,' m'
      print *,'   3rd layer shear ',shear3,' (1/s)'
      if(imoist.eq.1)then
       print *,' MOISTURE INFORMATION (rel hum %, heights m)'
       print *,'   RHSURF = ',rhsurf,' surface relative humidity '
       print *,'   RHVALUE1 = ',rhvalue1,' at RHHEIGHT1 = ',rhheight1
       print *,'   RHVALUE2 = ',rhvalue2,' at RHHEIGHT2 = ',rhheight2
       print *,'   RHVALUE3 = ',rhvalue3,' at RHHEIGHT3 = ',rhheight3
      endif
      print *,' ========================================= '

      if(ithermal.eq.1)then
        print *,' Initial thermal is activated '
        print *,'  DELT = ',delt,' amplitude (K)'
        print *,'  RADX = ',radx,' horizontal radius (m)'
        print *,'  RADZ = ',radz,' vertical radius (m)'
        print *,'  ZCNT = ',zcnt,' thermal center height (m)'
        print *,' ========================================= '
      else
        delt=0
        radx=0
        radz=0
        zcnt=0
      endif
      if(ishflux.eq.1)then
        print *,' Surface heat flux activated '
        print *,'  TDELT = ',tdelt,' initial ground-air T diff (K)'
        print *,'  ICOAST = ',icoast,' coastline location (gridpt)'
        print *,'  CDH = ',cdh,' heat flux coefficient (ndim)'
        print *,'  IRAND = ',irand,' randomize land surface heat flux'
        print *,' ========================================= '
      else
        tdelt=0
        icoast=0
        cdh=0
        irand=0
      endif
      if(istrfcn.eq.1)then
        print *,' Momentum streamfunction source activated '
        print *,'  S_REPEAT = ',s_repeat,' (1 for repeating source)'
        print *,'  S_AMPL = ',s_ampl,' source amplitude '
        print *,'  S_ZNAUGHT = ',s_znaught,' source center height (m)'
        print *,'  S_HWAVEL = ',s_hwavel,' source horiz wavelength (m)'
        print *,'  S_VWAVEL = ',s_vwavel,' source vert wavelength (m)'
        print *,'  S_PERIOD = ',s_period,' oscillation period (s)'
        print *,' ========================================= '
c manipulations for momentum source
        s_psik=2.*trigpi/s_hwavel
        s_psim=2.*trigpi/s_vwavel
        if(s_period.ne.0.)then 
         s_freq=2.*trigpi/s_period
        else
         s_freq=0.  ! period=0 turns oscillation off
        endif
      else
        hsdelt=0
        s_ampl=0
        s_znaught=0
        s_hwavel=0
        s_vwavel=0
        s_freq=0
      endif
      if(ihsrc.eq.1)then
        print *,' Heat source activated '
        print *,'  H_AMPL = ',h_ampl,' heat source amplitude (K/s)'
        print *,'  H_RADIUS_X = ',h_radius_x,' horizontal radius (m)'
        print *,'  H_RADIUS_Z = ',h_radius_z,' vertical radius (m)'
        print *,'  H_CENTER_Z = ',h_center_z,' source center height (m)'
        print *,'  H_FREQ = ',h_freq,' oscillation frequency (1/s)'
        print *,'  H_MODES = ',h_modes,' vertical modes (1 or 2)'
        print *,' ========================================= '
      else
        h_ampl=0
        h_radius_x=0
        h_radius_z=0
        h_center_z=0
        h_freq=0
        h_modes=0
      endif
      if(iseabreeze.eq.1)then
        print *,' Rotunno seabreeze heat source activated '
        print *,'  SB_AMPL = ',sb_ampl,' heat source amplitude (K/s)'
        print *,'  SB_X0 = ',sb_x0,' horizontal factor (m)'
        print *,'  SB_Z0 = ',sb_z0,' vertical depth (m)'
        print *,'  SB_PERIOD = ',sb_period,' source period (days)'
        print *,'  SB_LATITUDE = ',sb_latitude,' latitude (degrees)'
        print *,'  SB_LINEAR = ',sb_linear,' 1=linearize model '
        print *,' ========================================= '
      else
        sb_ampl=0
        sb_x0=0
        sb_z0=0
        sb_period=0
        sb_latitude=0
        sb_linear=0
      endif
      if(icoolzone.eq.1)then
        print *,' Storm-adaptive cooling zone activated '
        print *,'  CZ_AMPL = ',cz_ampl,' cooling zone amplitude (K)'
        print *,'  CZ_RIGHTEDGE = ',cz_rightedge, ' zone right edge (%)'
        print *,'  CZ_DEPTH = ',cz_depth,' zone depth (m)'
        print *,'  CZ_WIDTH = ',cz_width,' zone width (m)'
        print *,'  CZ_COOLRATE = ',cz_coolrate,' cooling rate (K/s)'
        print *,' ========================================= '
      else if(icoolzone.eq.2)then
        print *,' Impulsive cold block implemented '
        print *,'  CZ_AMPL = ',cz_ampl,' perturbation in cold block (K)'
        print *,'  CZ_RIGHTEDGE = ',cz_rightedge,' block right edge (%)'
        print *,'  CZ_DEPTH = ',cz_depth,' block depth (m)'
        print *,'  CZ_WIDTH = ',cz_width,' block width (m)'
        cz_coolrate=0.
      else
        cz_ampl=0.
        cz_rightedge=0.
        cz_depth=0.
        cz_width=0.
        cz_coolrate=0.
      endif
c other stuff
      cc1=1.0
      cc2=0.0
      i2d=1
      dy=dx
      tflux=1.
      iaccel=0
c constants
      g=9.8
      cp=1004.
      rd=287.
      psl=1000.e2
      hlv=2.5e6
c general manipulations
      xnsq_pbl=bvpbl*bvpbl
      xnsq_trop=bvtropo*bvtropo
      xnsq_strat=bvstrat*bvstrat
      psurf=psurf*100.   ! mb to Pa
      xk=rd/cp
      xki=1./xk
      cpv=1./(1.-xk)
      rcp=rd/cp
      cv=cp-rd
      cvrd=cv/rd
      csq=csnd*csnd
      hlvcp=hlv/cp 
      linear=0
      if(sb_linear.eq.1) linear=1

c are we running a case or not?
      isum=ithermal+ishflux+istrfcn+ihsrc+iseabreeze+icoolzone
      if(isum.eq.0)then
       print *,' ***** NO CASE SPECIFIED IN NAMELIST ***** '
       print *,' Please activate a case. '
       stop 'no_case'
      endif
      return
100   continue
      print *,' *********************************************** '
      print *,' THERE WAS AN ERROR READING THE INPUT NAMELIST.'
      print *,' The model will not run.'
      print *,' *********************************************** '
      stop 'error'
      end
