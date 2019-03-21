c
c     testing code for FMM - tests charges and dipoles against
c     O(N^2) direct method 
c
c
        implicit real *8 (a-h,o-z)
        parameter(lw=40 000 000)
        dimension w(lw)
c       
        call test_part(w,lw)
c
        stop
        end
c
c
c
c
c
        subroutine test_part(w,lw)
        implicit real *8 (a-h,o-z)
        dimension source(3,1 000 000)
        complex *16 charge(1 000 000)
        complex *16 dipstr(1 000 000)
        dimension dipvec(3,1 000 000)
        complex *16 cjvec(3,1 000 000)
        complex *16 cmvec(3,1 000 000)

        complex *16 pot(1 000 000)
        complex *16 fld(3,1 000 000)
c       
        complex *16 pot2(1 000 000)
        complex *16 fld2(3,1 000 000)
c       
        complex *16 evec(3,1 000 000)
        complex *16 hvec(3,1 000 000)
c
        complex *16 evec2(3,1 000 000)
        complex *16 hvec2(3,1 000 000)
c
        dimension target(3,2 000 000)
        complex *16 pottarg(2 000 000)
        complex *16 fldtarg(3,2 000 000)
c
        complex *16 ptemp,ftemp(3)
c       
ccc        parameter(lw=120 000 000)
        dimension w(1)
c       
        dimension xyz(3)
c
        real *8, allocatable :: sphererad(:), center(:,:)
        integer, allocatable :: icenter(:)
c
        complex *16, allocatable :: ampoleout(:,:,:)
        complex *16, allocatable :: ampoleinc(:,:,:)
        complex *16, allocatable :: ampoleout2(:,:,:)
        complex *16, allocatable :: ampoleinc2(:,:,:)
c
        complex *16, allocatable :: bmpoleout(:,:,:)
        complex *16, allocatable :: bmpoleinc(:,:,:)
        complex *16, allocatable :: bmpoleout2(:,:,:)
        complex *16, allocatable :: bmpoleinc2(:,:,:)
c
        dimension xnodes(10000),wts(10000)
        dimension rnodes(3,1000000),rwts(1000000)
c
        complex *16 ima
        complex *16 zk
        data ima/(0.0d0,1.0d0)/
c
c
        nterms0 = 3
        maxspheres = 1000000        
c
c
        allocate(sphererad(maxspheres))
        allocate(center(3,maxspheres))
        allocate(icenter(maxspheres))

        allocate(ampoleout(0:nterms0,-nterms0:nterms0,maxspheres))
        allocate(ampoleinc(0:nterms0,-nterms0:nterms0,maxspheres))
        allocate(ampoleout2(0:nterms0,-nterms0:nterms0,maxspheres))
        allocate(ampoleinc2(0:nterms0,-nterms0:nterms0,maxspheres))
c
        allocate(bmpoleout(0:nterms0,-nterms0:nterms0,maxspheres))
        allocate(bmpoleinc(0:nterms0,-nterms0:nterms0,maxspheres))
        allocate(bmpoleout2(0:nterms0,-nterms0:nterms0,maxspheres))
        allocate(bmpoleinc2(0:nterms0,-nterms0:nterms0,maxspheres))
c
c
        done=1
        pi=4*atan(done)
c
c
c       SET ALL PARAMETERS
c        
        call prini(6,13)
ccc        call prini(0,13)
c
        print *, 'ENTER n'
        read *, nsource
c
c
        call prinf('nsource=*',nsource,1)
c
        zk = 1.2d0 + ima*0.1d0
        zk = 1d0
c
c
        do i=1,lw
        w(i)=0
        enddo
c
        idist=4
c
        if( idist .eq. 2 ) then
c
c       ... construct charge distribution on a curve in R^3
c       
        do i=1,nsource
        a=2*pi*dble(i)/nsource
        hh=1.0d0/16.d0/4.0d0 / nsource
        source(1,i)=sin(a)
        source(2,i)=cos(a)
        source(3,i)=cos(a)
        sphererad(i) = hh
        center(1,i)=source(1,i)
        center(2,i)=source(2,i)
        center(3,i)=source(3,i)
        enddo
c
        endif
c       
        if( idist .eq. 4 ) then
c
c       ... construct the grid of charges
c
        ngrid=nsource-1
        hh=1.0d0/16.d0/4.0d0
        kk=0
        do i=1,ngrid+1
        do j=1,ngrid+1
        do k=1,ngrid+1
        kk=kk+1
        source(1,kk)=(i-1.0)/ngrid
        source(2,kk)=(j-1.0)/ngrid
        source(3,kk)=(k-1.0)/ngrid
        source(1,kk)=source(1,kk)+hkrand(0)*.01
        source(2,kk)=source(2,kk)+hkrand(0)*.01
        source(3,kk)=source(3,kk)+hkrand(0)*.01
ccc        call prin2('source=*',source(1,kk),3)
        sphererad(kk) = hh
        center(1,kk)=source(1,kk)
        center(2,kk)=source(2,kk)
        center(3,kk)=source(3,kk)
        enddo
        enddo
        enddo
        nsource=kk
c       
        call prinf('after grid build, nsource=*',nsource,1)
c
        endif
c
c
        if( idist .eq. 5 ) then
c       ... construct sphere distribution on a curve in R^3
c       
        rbig = 1.0d0
	hh = rbig/4.0d0
	x1 = -rbig/2+hh/2
	y1 = -rbig/2+hh/2
	z1 = -rbig/2+hh/2
        do i=1,nsource
           sphererad(i) = 0.25*(2*pi*rbig/(2*nsource))
           sphererad(i) = hh/16.0d0 
        enddo
        do i=1,nsource
        a=2*pi*dble(i)/nsource
        center(1,i)=x1
        center(2,i)=y1
        center(3,i)=z1
	x1 = x1 + hh
	y1 = y1 + hh
	z1 = z1 + hh
c        center(1,i)=rbig*sin(a)
c        center(2,i)=rbig*cos(a)
c        center(3,i)=cos(3.1*a)
c        center(1,i)=125*sin(a)
c        center(2,i)=125*cos(a)
c        center(3,i)=125*cos(3.1*a)
c        center(1,i)=250*rand(0)
c        center(2,i)=250*rand(0)
c        center(3,i)=0

        source(1,i)=center(1,i)
        source(2,i)=center(2,i)
        source(3,i)=center(3,i)

        enddo
c
ccc        call prin2('sphererad=*',sphererad,nsource)
ccc        call prin2('center=*',center,3*nsource)
c
        endif
c
c
c       ... set up the targets
c
        if( idist .eq. 2 .or. 
     $     idist .eq. 4 .or. idist .eq. 5 ) then
        do i=1,nsource
        target(1,i)=source(1,i) + 0.1*sphererad(i) 
        target(2,i)=source(2,i) - 0.2*sphererad(i)
        target(3,i)=source(3,i) + 0.3*sphererad(i)
        enddo 
        ntarget=nsource
        endif
c
        ntarget=0
        call prinf('ntarget=*',ntarget,1)
c       
ccc        iw=71
ccc        read(iw,*) nsource
ccc        do i=1,nsource
ccc        read(iw,*) source(1,i),source(2,i),source(3,i)
ccc        enddo
c
ccc        call prinf('read data from file iw=*',iw,1)
ccc        call prinf('nsource=*',n,1)
c
c
        do i=1,nsource
        icenter(i)=i
        enddo
c
        do i = 1,nsource
        call em3dzero(ampoleout(0,-nterms0,i),nterms0)
        call em3dzero(ampoleout2(0,-nterms0,i),nterms0)
        call em3dzero(ampoleinc(0,-nterms0,i),nterms0)
        call em3dzero(ampoleinc2(0,-nterms0,i),nterms0)
        call em3dzero(bmpoleout(0,-nterms0,i),nterms0)
        call em3dzero(bmpoleout2(0,-nterms0,i),nterms0)
        call em3dzero(bmpoleinc(0,-nterms0,i),nterms0)
        call em3dzero(bmpoleinc2(0,-nterms0,i),nterms0)
        enddo
c
c
        iprec=1
c       
        call prinf('iprec=*',iprec,1)
c       
        ifpot=1
        iffld=1
c
        ifcharge=1
        ifdipole=0
c
        ifpottarg=1
        iffldtarg=0
c
        ifcjvec=1
        ifcmvec=0
c
        itest=2
c
c       
c
        if (itest .eq. 2) then
c
        if (ifcharge .eq. 1 ) then
c
        do i=1,nsource
        charge(i)=0
        enddo
c
        endif
c
        if (ifdipole .eq. 1) then
c
        do i=1,nsource
           dipstr(i)=0
           dipvec(1,i)=0
           dipvec(2,i)=0
           dipvec(3,i)=0
        enddo
c
        endif
c
        if (ifcjvec .eq. 1 ) then
c
        do i=1,nsource
        cjvec(1,i)=0
        cjvec(2,i)=0
        cjvec(3,i)=0
        enddo
c
        endif
c
        if (ifcmvec .eq. 1 ) then
c
        do i=1,nsource
        cjvec(1,i)=0
        cjvec(2,i)=0
        cjvec(3,i)=0
        enddo
c
        endif
c

        if (ifcharge .eq. 1 ) then
c
        do i=1,1
        charge(i)=1
        enddo
c
        endif
c
        if (ifdipole .eq. 1) then
c
        do i=1,1
           dipstr(i)=1
           dipvec(1,i)=1
           dipvec(2,i)=2
           dipvec(3,i)=3
        enddo
c
        endif
c
        if (ifcjvec .eq. 1 ) then
c
        do i=1,1
        cjvec(1,i)=1
        cjvec(2,i)=2
        cjvec(3,i)=3
        enddo
c
        endif
c
        if (ifcmvec .eq. 1 ) then
c
        do i=1,1
        cjvec(1,i)=1
        cjvec(2,i)=2
        cjvec(3,i)=3
        enddo
c
        endif
c
c
        do i=1,1
c
        xyz(1) = center(1,i)+0.2*sphererad(i)
        xyz(2) = center(2,i)+0.1*sphererad(i)
        xyz(3) = center(3,i)-0.2*sphererad(i)
        npts = 1
c
        scale0 = 1
        call em3formmp
     $     (zk,xyz,cjvec(1,i),cmvec(1,i),npts,center(1,i),
     1     ampoleout(0,-nterms0,i),bmpoleout(0,-nterms0,i),nterms0)
        call em3formmp
     $     (zk,xyz,cjvec(1,i),cmvec(1,i),npts,center(1,i),
     1     ampoleout2(0,-nterms0,i),bmpoleout2(0,-nterms0,i),nterms0)

        enddo
c
        endif
c
c
        t1=second()
C$        t1=omp_get_wtime()
c       
        call emfmm3dsph(ier,iprec,zk,
     $     nsource,source,sphererad,ampoleout,bmpoleout,
     $     ampoleinc,bmpoleinc,nterms0)
c       
        t2=second()
C$        t2=omp_get_wtime()
c       
c       
        call prinf('nsource=*',nsource,1)
        call prinf('ntarget=*',ntarget,1)
        call prin2('after fmm, time (sec)=*',t2-t1,1)
ccc        call prin2('after fmm, speed (points/sec)=*',nsource/(t2-t1),1)
        call prin2('after fmm, speed (spheres/sec)=*',
     $     (nsource+ntarget)/(t2-t1),1)
c       
c
ccc        call prin2('after fmm, mpoleout=*',ampoleout,2*12)
c        call prin2('after fmm, mpoleinc=*',ampoleinc,2*12)
c        call prin2('after fmm, mpoleinc=*',ampoleinc,
c     $     2*(nterms0+1)*(2*nterms0+1)*nsource)

c
ccc        m=nsource
ccc        m=min(nsource,100)
        m=min(nsource,20)
c
c
        ifprint=1
        if (ifprint .eq. 1) then
ccc        call prin2('source=*',source,3*nsource)
        endif

        t1=second()
C$        t1=omp_get_wtime()
c
        nquad=2*nterms0
        ifinit=1
        call legewhts(nquad,xnodes,wts,ifinit)
c
        itype=1
        nphi=nquad+1
        call fftnextn(nphi,nphifft)
        nphi=nphifft
        ntheta=nquad+1
        call prinf('nphi=*',nphi,1)
        call prinf('ntheta=*',ntheta,1)
        call e3fgrid(itype,nquad,nphi,ntheta,rnodes,rwts,nnodes)        
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,ftemp) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4) 
        do 6215 j=1,m
        do 6214 i=1,nsource
        if( i .eq. j ) goto 6214

        call em3mpta3_add(zk,source(1,i),
     $     ampoleout2(0,-nterms0,i),bmpoleout2(0,-nterms0,i),nterms0,
     $     source(1,j),
     $     ampoleinc2(0,-nterms0,j),bmpoleinc2(0,-nterms0,j),nterms0,
     $     sphererad(j),rnodes,rwts,nphi,ntheta)

 6214   continue
 6215   continue
C$OMP END PARALLEL DO
c
        t2=second()
C$        t2=omp_get_wtime()
c
        call prin2('directly, estimated time (sec)=*',
     $     (t2-t1)*dble(nsource)/dble(m),1)
        call prin2('directly, estimated speed (spheres/sec)=*',
     $     m/(t2-t1),1)
c       
        call h3derror_mpole(ampoleinc,ampoleinc2,nterms0,m,aerr,rerr)
ccc        call prin2('absolute L2 error in mpole=*',aerr,1)
        call prin2('relative L2 error in A-mpole=*',rerr,1)
c
        call h3derror_mpole(bmpoleinc,bmpoleinc2,nterms0,m,aerr,rerr)
ccc        call prin2('absolute L2 error in mpole=*',aerr,1)
        call prin2('relative L2 error in B-mpole=*',rerr,1)
c
c       
c       now, check the induced E and H fields inside each sphere
c       

        do i=1,nsource
        evec(1,i)=0
        evec(2,i)=0
        evec(3,i)=0
        hvec(1,i)=0
        hvec(2,i)=0
        hvec(3,i)=0
        enddo

        do i=1,nsource
        evec2(1,i)=0
        evec2(2,i)=0
        evec2(3,i)=0
        hvec2(1,i)=0
        hvec2(2,i)=0
        hvec2(3,i)=0
        enddo


C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,ftemp) 
C$OMP$SCHEDULE(DYNAMIC)
C$OMP$NUM_THREADS(1) 
        do i=2,m

        call dipole3emt(zk,xyz,target(1,i),
     $     cjvec(1,1),cmvec(1,1),evec(1,i),hvec(1,i))

        call em3taeval(zk,center(1,i),
     $     ampoleinc(0,-nterms0,i),bmpoleinc(0,-nterms0,i),nterms0,
     $     target(1,i),evec2(1,i),hvec2(1,i))

        enddo
C$OMP END PARALLEL DO
c
        m1=m-1

        if( ifprint .eq. 1 ) then
        call prin2('directly, E=*',evec(1,2),3*2*m1)
        call prin2('directly, H=*',hvec(1,2),3*2*m1)

        call prin2('via fmm, E=*',evec2(1,2),3*2*m1)
        call prin2('via fmm, H=*',hvec2(1,2),3*2*m1)
        endif

c
        call h3derror(evec,evec2,3*m1,aerr,rerr)
ccc         call prin2('absolute L2 error in E-field=*',aerr,1)
        call prin2('relative L2 error in E-field=*',rerr,1)

        call h3derror(hvec,hvec2,3*m1,aerr,rerr)
ccc         call prin2('absolute L2 error in H-field=*',aerr,1)
        call prin2('relative L2 error in H-field=*',rerr,1)

        return
        end
c
c
c
c
c
        subroutine h3dmperr(mpole1,mpole2,nterms,d)
        implicit real *8 (a-h,o-z)
c       
        complex *16 mpole1(0:nterms,-nterms:nterms)
        complex *16 mpole2(0:nterms,-nterms:nterms)
c       
        d=0
c       
        do n=0,nterms
        do m=-n,n
        d=d+abs(mpole1(n,m)-mpole2(n,m))**2
        enddo
        enddo
c       
        d=d/(nterms+1)
        d=d/(2*nterms+1)
        d=sqrt(d)
c       
        return
        end
c
c
c
c
c
        subroutine h3dmpnorm(mpole,nterms,d)
        implicit real *8 (a-h,o-z)
c
        complex *16 mpole(0:nterms,-nterms:nterms)
c
        d=0
c
        do n=0,nterms
        do m=-n,n
        d=d+abs(mpole(n,m))**2
        enddo
        enddo
c
        d=d/(nterms+1)
        d=d/(2*nterms+1)
        d=sqrt(d)
c
        return
        end
c
c
c
c
c
        subroutine h3derror(pot1,pot2,n,ae,re)
        implicit real *8 (a-h,o-z)
c
c       evaluate absolute and relative errors
c
        complex *16 pot1(n),pot2(n)
c
        d=0
        a=0
c       
        do i=1,n
        d=d+abs(pot1(i)-pot2(i))**2
        a=a+abs(pot1(i))**2
        enddo
c       
        d=d/n
        d=sqrt(d)
        a=a/n
        a=sqrt(a)
c       
        ae=d
        re=d/a
c       
        return
        end
c
c
c
c
c
        subroutine h3derror_mpole
     $     (ampole,ampole2,nterms,m,ae,re)
        implicit real *8 (a-h,o-z)
c
c       evaluate absolute and relative errors
c
        complex *16 ampole(0:nterms,-nterms:nterms,m)
        complex *16 ampole2(0:nterms,-nterms:nterms,m)
c
        d=0
        a=0
c
        do i = 1,m
           do j=0,nterms
           do k=-j,j
           d=d+abs((ampole(j,k,i)-ampole2(j,k,i)))**2 
           a=a+abs(ampole2(j,k,i))**2 
           enddo
           enddo
        enddo
c
        n = (nterms+1)*(2*nterms+1)
        d=d/n
        d=sqrt(d)
        a=a/n
        a=sqrt(a)
c
        ae=d
        re=d/a
c
        return
        end
c
c
c
c
c
c
        subroutine fftnextn(n,m)
        implicit real *8 (a-h,o-z)
        dimension ipow2(101)
        data ifinit/0/
c
c       ... number of fft points
c
        if( ifinit .eq. 0 ) then
        ifinit=1
        ipow2(1)=1
        do i=2,100
        ipow2(i)=ipow2(i-1)*2
        enddo
        endif
c       
        m=n
        if( mod(n,2) .eq. 1 ) m=m+1
c
        return
        end
c
c
