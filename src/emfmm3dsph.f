cc Copyright (C) 2009-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c       
c     This file contains the main FMM routines and some related
c     subroutines for evaluating Maxwell sphere (multipole) interactions.
c     (FORTRAN 90 VERSION)
c
c     emfmm3dsph - Maxwell FMM in R^3: evaluate all pairwise sphere
c         interactions (ignoring self-interaction)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the Maxwell sphere FMM in R^3
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine emfmm3dsph(ier,iprec,zk,
     $     nsource,source,sphererad,ampoleout,bmpoleout,
     $     ampoleinc,bmpoleinc,nterms0)
        implicit real *8 (a-h,o-z)
c              
c       Maxwell sphere FMM in R^3: evaluate all pairwise sphere
c       interactions (ignoring self-interaction) 
c
c       We use (exp(ikr)/r) for the Green's function,
c       without the (1/4 pi) scaling.  Self-interactions are not included.
c   
c       This is primarily a memory management code. 
c       The actual work is carried out in subroutine emfmm3dspheremain.
c
c       INPUT PARAMETERS:
c
c       iprec:  FMM precision flag
c
c                 -2 => tolerance =.5d0
c                 -1 => tolerance =.5d-1
c                  0 => tolerance =.5d-2
c                  1 => tolerance =.5d-3
c                  2 => tolerance =.5d-6
c                  3 => tolerance =.5d-9
c                  4 => tolerance =.5d-12
c                  5 => tolerance =.5d-15
c
c       zk: complex *16: Helmholtz parameter
c       nsource: integer:  number of source spheres
c       source: real *8 (3,nsource):  sphere locations
c       sphererad: real *8 (3,nsource):  sphere radii
c       ampoleout: complex *16 (0:nterms,-nterms:nterms,nsource): 
c                    outgoing multipole expansions
c       nterms: number of terms in multipole expansions
c
c       OUTPUT PARAMETERS:
c
c       ier   =  error return code
c                ier=0     =>  normal execution
c                ier=4     =>  cannot allocate tree workspace
c                ier=8     =>  cannot allocate bulk FMM  workspace
c                ier=16    =>  cannot allocate mpole expansion
c                              workspace in FMM
c
c       ampoleinc: complex *16 (0:nterms,-nterms:nterms,nsource): 
c                    incoming multipole expansions
c
c       DEVELOPMENT NOTES:
c  
c     10/2/08 In forming multipole expansions, spheres may be near box
c     corner.  The radius used for the box must be larger than for mpmp
c     shifts - otherwise the field sampled is too nonsmooth.  using
c     radius = radius*1.5 for now, but this must be handled more
c     carefully.
c
c
        dimension source(3,nsource)
        dimension sphererad(nsource)
        complex *16 ampoleout(0:nterms0,-nterms0:nterms0,nsource)
        complex *16 bmpoleout(0:nterms0,-nterms0:nterms0,nsource)
        complex *16 ampoleinc(0:nterms0,-nterms0:nterms0,nsource)
        complex *16 bmpoleinc(0:nterms0,-nterms0:nterms0,nsource)
        complex *16 zk
        complex *16 ima
        dimension timeinfo(10)
c
c     Note: various arrays dimensioned here to 200.
c     That allows for 200 evels of refinment, which is 
c     more than enough for any non-pathological case.
c
 
        dimension laddr(2,200)
        dimension bsize(0:200)
        dimension nterms(0:200)
        integer box(20)
        integer box1(20)
        dimension scale(0:200)
        dimension center(3)
        dimension center0(3),corners0(3,8)
        dimension center1(3),corners1(3,8)
        real *8, allocatable :: w(:)
        real *8, allocatable :: wlists(:)
        real *8, allocatable :: wrmlexp(:)
        complex *16 ptemp,ftemp(3)
c       
        dimension target(3)
        data ima/(0.0d0,1.0d0)/
c       
        ntarget=0
c
        ier=0
        lused7 = 0
c       
        done=1
        pi=4*atan(done)
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c       
        ifprint=1
c
c     set fmm tolerance based on iprec flag.
c
        if( iprec .eq. -2 ) epsfmm=.5d-0 
        if( iprec .eq. -1 ) epsfmm=.5d-1
        if( iprec .eq. 0 ) epsfmm=.5d-2
        if( iprec .eq. 1 ) epsfmm=.5d-3
        if( iprec .eq. 2 ) epsfmm=.5d-6
        if( iprec .eq. 3 ) epsfmm=.5d-9
        if( iprec .eq. 4 ) epsfmm=.5d-12
        if( iprec .eq. 5 ) epsfmm=.5d-15
        if( iprec .eq. 6 ) epsfmm=0
c       
        if (ifprint .ge. 1) call prin2('epsfmm=*',epsfmm,1)
c
c
c     set criterion for box subdivision (number of sources per box)
c
        if( iprec .eq. -2 ) nbox=40/4
        if( iprec .eq. -1 ) nbox=50/4
        if( iprec .eq. 0 ) nbox=80/4
        if( iprec .eq. 1 ) nbox=160/4
        if( iprec .eq. 2 ) nbox=400/4
        if( iprec .eq. 3 ) nbox=800/4
        if( iprec .eq. 4 ) nbox=1200/4
        if( iprec .eq. 5 ) nbox=1400/4
        if( iprec .eq. 6 ) nbox=nsource+ntarget
c
        if (ifprint .ge. 1) call prinf('nbox=*',nbox,1)
c
c
c     create oct-tree data structure
c
        ntot = 100*(nsource+ntarget)+10000
        do ii = 1,10
           allocate (wlists(ntot))
           call emfmm3dparttree(ier,iprec,zk,
     $        nsource,source,ntarget,target,
     $        nbox,epsfmm,iisource,iitarget,iwlists,lwlists,
     $        nboxes,laddr,nlev,center,size,
     $        wlists,ntot,lused7)
           if (ier.ne.0) then
              deallocate(wlists)
              ntot = ntot*1.5
              call prinf(' increasing allocation, ntot is *',ntot,1)
           else
             goto 1200
           endif
        enddo
1200    continue
        if (ier.ne.0) then
           call prinf(' exceeded max allocation, ntot is *',ntot,1)
           ier = 4
           return
        endif
c
c     lused7 is counter that steps through workspace,
c     keeping track of total memory used.
c
        lused7=1
c
c       ... prepare data structures 
c
        do i = 0,nlev
        scale(i) = 1.0d0
        boxsize = abs((size/2.0**i)*zk)
cccccccccccccccc        if (boxsize .lt. 1) scale(i) = boxsize
        enddo
c       
        if (ifprint .ge. 1) call prin2('scale=*',scale,nlev+1)
c       
c
c       carve up workspace further
c
c     isourcesort is pointer for sorted source coordinates
c     itargetsort is pointer for sorted target locations
c     iampoleoutsort is pointer for sorted incoming multipoles
c     iampoleincsort is pointer for sorted outgoing multipoles
c
c
        isourcesort = lused7 + 5
        lsourcesort = 3*nsource
        isphereradsort = isourcesort+lsourcesort
        lsphereradsort = nsource
        iampoleoutsort = isphereradsort+lsphereradsort
        lampoleoutsort = 2*(nterms0+1)*(2*nterms0+1)*nsource
        ibmpoleoutsort = iampoleoutsort+lampoleoutsort
        lbmpoleoutsort = 2*(nterms0+1)*(2*nterms0+1)*nsource
        iampoleincsort = ibmpoleoutsort+lbmpoleoutsort
        lampoleincsort = 2*(nterms0+1)*(2*nterms0+1)*nsource
        ibmpoleincsort = iampoleincsort+lampoleincsort
        lbmpoleincsort = 2*(nterms0+1)*(2*nterms0+1)*nsource
        lused7 = ibmpoleincsort+lbmpoleincsort
c
        if (ifprint .ge. 1) call prinf(' lused7 is *',lused7,1)
c
c       based on FMM tolerance, compute expansion lengths nterms(i)
c 
        nmax = 0
        do i = 0,nlev
           bsize(i)=size/2.0d0**i
           call h3dterms(bsize(i),zk,epsfmm, nterms(i), ier)
           if (nterms(i).gt. nmax) nmax = nterms(i)
ccc           if (nterms(i).gt. nmax .and. i.ge. 2) nmax = nterms(i)
        enddo
c
        if (ifprint.eq.1) 
     $     call prin2('in emfmm3dpart, bsize(0) zk/2 pi=*',
     $     abs(bsize(0)*zk)/2/pi,1)
c
        if (ifprint.eq.1) call prin2('zk=*',zk,2)
        if (ifprint.eq.1) call prin2('bsize=*',bsize,nlev+1)
c
        nquad=2*nmax        
c       
c     ixnodes is pointer for quadrature nodes
c     iwhts is pointer for quadrature weights
c
        ixnodes = lused7 
        iwts = ixnodes + nquad
        lused7 = iwts + nquad
c
        if (ifprint .ge. 1) call prinf('nterms=*',nterms,nlev+1)
        if (ifprint .ge. 1) call prinf('nmax=*',nmax,1)
c
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(4,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c   
c       ... allocate iaddr and temporary arrays
c
        iiaddr = lused7 
        imptemp = iiaddr + 4*nboxes
        lmptemp = (nmax+1)*(2*nmax+1)*2
        lused7 = imptemp + lmptemp
        allocate(w(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate bulk FMM workspace,
     1                  lused7 is *',lused7,1)
           ier = 8
           return
        endif
c
c     reorder sources, targets so that each box holds
c     contiguous list of source/target numbers.
c
        call em3dsphreorder(nsource,source,sphererad,
     $     ampoleout,bmpoleout,wlists(iisource),
     1     w(isourcesort),w(isphereradsort),
     $     w(iampoleoutsort),w(ibmpoleoutsort),nterms0) 
c       
        if (ifprint .ge. 1) call prinf('finished reordering=*',ier,1)
        if (ifprint .ge. 1) call prinf('ier=*',ier,1)
        if (ifprint .ge. 1) call prinf('nboxes=*',nboxes,1)
        if (ifprint .ge. 1) call prinf('nlev=*',nlev,1)
        if (ifprint .ge. 1) call prinf('nboxes=*',nboxes,1)
        if (ifprint .ge. 1) call prinf('lused7=*',lused7,1)
c
        ifinit=1
        call legewhts(nquad,w(ixnodes),w(iwts),ifinit)
c
ccc        call prin2('xnodes=*',xnodes,nquad)
ccc        call prin2('wts=*',wts,nquad)

c     allocate memory need by multipole, local expansions at all
c     levels
c     irmlexp is pointer for workspace need by various fmm routines,
c
        call em3dsphmpalloc
     $     (wlists(iwlists),w(iiaddr),nboxes,lmptot,nterms)
c
        if (ifprint .ge. 1) call prinf(' lmptot is *',lmptot,1)
c       
        irmlexp = 1
        lused7 = irmlexp + lmptot 
        if (ifprint .ge. 1) call prinf(' lused7 is *',lused7,1)
        allocate(wrmlexp(lused7),stat=ier)
        if (ier.ne.0) then
           call prinf(' cannot allocate mpole expansion workspace,
     1                  lused7 is *',lused7,1)
           ier = 16
           return
        endif
c
c       
ccc        do i=lused7+1,lused7+1+100
ccc        w(i)=777
ccc        enddo
c
c     Memory allocation is complete. 
c     Call main fmm routine. There are, unfortunately, a lot
c     of parameters here. ifevalfar and ifevalloc determine
c     whether far field and local fields (respectively) are to 
c     be evaluated. Setting both to 1 means that both will be
c     computed (which is the normal scenario).
c
        ifevalfar=1
        ifevalloc=1
c
        call emfmm3dsphmain(ier,iprec,zk,
     $     ifevalfar,ifevalloc,
     $     nsource,w(isourcesort),wlists(iisource),
     $     w(isphereradsort),
     $     w(iampoleoutsort),w(ibmpoleoutsort),
     $     w(iampoleincsort),w(ibmpoleincsort),nterms0,
     $     epsfmm,w(iiaddr),wrmlexp(irmlexp),w(imptemp),lmptemp,
     $     w(ixnodes),w(iwts),nquad,
     $     nboxes,laddr,nlev,scale,bsize,nterms,
     $     wlists(iwlists),lwlists)
c
c       parameter ier from targmain routine is currently meaningless, reset to 0
        if( ier .ne. 0 ) ier = 0
c
        if (ifprint .ge. 1) call prinf('lwlists=*',lused,1)
        if (ifprint .ge. 1) call prinf('lused total =*',lused7,1)       
c
        if (ifprint .ge. 1) 
     $      call prin2('memory / point = *',(lused7)/dble(nsource),1)
c       
ccc        call prin2('after w=*', w(1+lused7-100), 2*100)
c
        call em3dmpolesort(nsource,wlists(iisource),
     $     w(iampoleincsort),ampoleinc,nterms0)
c       
        call em3dmpolesort(nsource,wlists(iisource),
     $     w(ibmpoleincsort),bmpoleinc,nterms0)
c       
        return
        end
c
c
c
c
c
        subroutine emfmm3dsphmain(ier,iprec,zk,
     $     ifevalfar,ifevalloc,
     $     nsource,sourcesort,isource,
     $     sphereradsort,ampoleoutsort,bmpoleoutsort,
     $     ampoleincsort,bmpoleincsort,nterms0,
     $     epsfmm,iaddr,rmlexp,mptemp,lmptemp,xnodes,wts,nquad,
     $     nboxes,laddr,nlev,scale,bsize,nterms,
     $     wlists,lwlists)
        implicit real *8 (a-h,o-z)
        dimension sourcesort(3,1), isource(1)
        dimension sphereradsort(nsource)
        complex *16 ampoleoutsort(0:nterms0,-nterms0:nterms0,nsource)
        complex *16 bmpoleoutsort(0:nterms0,-nterms0:nterms0,nsource)
        complex *16 ampoleincsort(0:nterms0,-nterms0:nterms0,nsource)
        complex *16 bmpoleincsort(0:nterms0,-nterms0:nterms0,nsource)
        complex *16 zk
        complex *16 ima
        dimension wlists(1)
        dimension iaddr(4,nboxes)
        real *8 rmlexp(1)
        complex *16 mptemp(lmptemp)
        dimension xnodes(nquad),wts(nquad)
        dimension timeinfo(10)
        dimension center(3)
        dimension laddr(2,200)
        dimension scale(0:200)
        dimension bsize(0:200)
        dimension nterms(0:200)
        dimension list(10 000)
        complex *16 ptemp,ftemp(3)
        integer box(20)
        dimension center0(3),corners0(3,8)
        integer box1(20)
        dimension center1(3),corners1(3,8)
        dimension itable(-3:3,-3:3,-3:3)
        dimension wlege(100 000)
        dimension nterms_eval(4,0:200)
c
        real *8, allocatable :: xnodes2(:), wts2(:)
        real *8, allocatable :: rnodes(:,:), rwts(:)
c
        data ima/(0.0d0,1.0d0)/

c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c       
        ifprint=1
c
        max_nodes = 10000
        allocate( xnodes2(max_nodes) )
        allocate( wts2(max_nodes) )
c
        max_nodes = 100000
        allocate( rnodes(3,max_nodes) )
        allocate( rwts(max_nodes) )
c
c
c       ... set the incoming multipoles to zero
c
        do i=1,nsource
        call em3dzero(ampoleincsort(0,-nterms0,i),nterms0)
        call em3dzero(bmpoleincsort(0,-nterms0,i),nterms0)
        enddo
c       
        do i=1,10
        timeinfo(i)=0
        enddo
c
c
        if( ifevalfar .eq. 0 ) goto 8000
c       
c
c       ... initialize Legendre function evaluation routines
c
        nlege=200
        lw7=100 000
        call ylgndrfwini(nlege,wlege,lw7,lused7)
ccc        write(*,*)' lused7 from  ylgndrfwini is',lused7
c
        do i=0,nlev
        do itype=1,4
        call h3dterms_eval(itype,bsize(i),zk,epsfmm,
     1       nterms_eval(itype,i),ier)
        enddo
        enddo
c
        if (ifprint .ge. 2) 
     $     call prinf('nterms_eval=*',nterms_eval,4*(nlev+1))
c
c       ... set all multipole and local expansions to zero
c
        do ibox = 1,nboxes
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        level=box(1)
        call em3dzero(rmlexp(iaddr(1,ibox)),nterms(level))
        call em3dzero(rmlexp(iaddr(2,ibox)),nterms(level))
        call em3dzero(rmlexp(iaddr(3,ibox)),nterms(level))
        call em3dzero(rmlexp(iaddr(4,ibox)),nterms(level))
        enddo
c
c
        if(ifprint .ge. 1) 
     $     call prinf('=== STEP 1 (form mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions
c
ccc        do 1200 ibox=1,nboxes
        do 1300 ilev=3,nlev+1
        nquad2=nterms(ilev-1)*1.2
        nquad2=max(6,nquad2)
        ifinit2=1
        call legewhts(nquad2,xnodes2,wts2,ifinit2)

        itype=1
        nphi=nquad2+1
        call em3fftnextn(nphi,nphifft)
        nphi=nphifft
        ntheta=nquad2+1
cc        call prinf('nphi=*',nphi,1)
cc        call prinf('ntheta=*',ntheta,1)
        call e3fgrid(itype,nquad2,nphi,ntheta,rnodes,rwts,nnodes)        

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,radius)
C$OMP$PRIVATE(lused,ier,i,j,ptemp,ftemp,cd) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 1200 ibox=laddr(1,ilev),laddr(1,ilev)+laddr(2,ilev)-1
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        level=box(1)
c
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
c        ipts=box(14)
c        npts=box(15)
c        call prinf('ipts=*',ipts,1)
c        call prinf('npts=*',npts,1)
        npts=box(15)
        if (ifprint .ge. 2) then
           call prinf('npts=*',npts,1)
           call prinf('isource=*',isource(box(14)),box(15))
        endif
        endif
c
c       ... prune all sourceless boxes
c
        if( box(15) .eq. 0 ) goto 1200
c
        if (nkids .eq. 0) then
c
c       ... form multipole expansions
c
	    radius = (corners0(1,1) - center0(1))**2
	    radius = radius + (corners0(2,1) - center0(2))**2
	    radius = radius + (corners0(3,1) - center0(3))**2
	    radius = sqrt(radius)
c
c       to avoid blow up if a source is exactly in the box corner,
c       make radius larger...
c
            radius = radius*1.5d0
c
            do i=box(14),box(14)+box(15)-1

        call em3mpmp3_add(zk,
     $     sourcesort(1,i),
     $     ampoleoutsort(0,-nterms0,i),bmpoleoutsort(0,-nterms0,i),
     $     nterms0,
     $     center0,rmlexp(iaddr(1,ibox)),rmlexp(iaddr(2,ibox)),
     $     nterms(level),
     $     radius,rnodes,rwts,nphi,ntheta)

            enddo

         endif
c
 1200    continue
C$OMP END PARALLEL DO
 1300    continue
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(1)=t2-t1
c       
         if(ifprint .ge. 1)
     $      call prinf('=== STEP 2 (form lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 2, adaptive part, form local expansions, 
c           or evaluate the potentials and fields directly
c 
         do 3251 ibox=1,nboxes
c
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
c
         itype=3
         call d3tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then 
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list3=*',list,nlist)
            endif
         endif
c
c       ... prune all sourceless boxes
c
         if( box(15) .eq. 0 ) nlist=0
c
c
c       ... note that lists 3 and 4 are dual
c
c       ... form local expansions for all boxes in list 3
c       ... if target is childless, evaluate directly (if cheaper)
c        
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(level,npts,nkids)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect3,radius)
C$OMP$PRIVATE(lused,ier,i,j,ptemp,ftemp,cd,ilist) 
cccC$OMP$SCHEDULE(DYNAMIC)
C$OMP$NUM_THREADS(1) 
         do ilist=1,nlist
            jbox=list(ilist)
            call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c        
            level1=box1(1)
c
            ifdirect3 = 0
            if( box1(15) .lt. (nterms(level1)+1)**2/4 .and.
     $          box(15) .lt. (nterms(level1)+1)**2/4 ) ifdirect3 = 1
c
            ifdirect3 = 0
c
            if( ifdirect3 .eq. 0 ) then
               npts=box(15)
               if_use_trunc = 1
c
	    radius = (corners1(1,1) - center1(1))**2
	    radius = radius + (corners1(2,1) - center1(2))**2
	    radius = radius + (corners1(3,1) - center1(3))**2
	    radius = sqrt(radius)
c
c       to avoid blow up if a source is exactly in the box corner,
c       make radius larger...
c
            radius = radius*1.5d0
c
            nquad2=nterms(level1-1)*1.2
            nquad2=max(6,nquad2)
            ifinit2=1
            call legewhts(nquad2,xnodes2,wts2,ifinit2)
c
        itype=1
        nphi=nquad2+1
        call em3fftnextn(nphi,nphifft)
        nphi=nphifft
        ntheta=nquad2+1
cc        call prinf('nphi=*',nphi,1)
cc        call prinf('ntheta=*',ntheta,1)
        call e3fgrid(itype,nquad2,nphi,ntheta,rnodes,rwts,nnodes)        

            do i=box(14),box(14)+box(15)-1

        call em3mpta3_add(zk,
     $     sourcesort(1,i),
     $     ampoleoutsort(0,-nterms0,i),bmpoleoutsort(0,-nterms0,i),
     $     nterms0,
     $     center1,rmlexp(iaddr(3,jbox)),rmlexp(iaddr(4,jbox)),
     $     nterms(level1),
     $     radius,rnodes,rwts,nphi,ntheta)

            enddo

            else

            write(*,*) '***not implemented***'
            stop
c            call emfmm3dpart_direct(zk,box,box1,sourcesort,
c     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
c     $         ifpot,pot,iffld,fld,
c     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)

            endif
         enddo
C$OMP END PARALLEL DO
c
 3251    continue
c
         t2=second()
C$        t2=omp_get_wtime()
ccc        call prin2('time=*',t2-t1,1)
         timeinfo(2)=t2-t1
c
        if(ifprint .ge. 1)
     $      call prinf('=== STEPS 3,4,5 ====*',i,0)
        ifprune_list2 = 1
ccc        if (ifpot.eq.1) ifprune_list2 = 0
ccc        if (iffld.eq.1) ifprune_list2 = 0
        ifprune_list2 = 0
        call emfmm3d_list2
     $     (zk,bsize,nlev,laddr,scale,nterms,rmlexp,iaddr,epsfmm,
     $     timeinfo,wlists,mptemp,lmptemp,xnodes,wts,nquad,
     $     ifprune_list2)
c
c
        nquad2=nterms0*1.2
        nquad2=max(6,nquad2)
        ifinit2=1
        call legewhts(nquad2,xnodes2,wts2,ifinit2)
c
        itype=1
        nphi=nquad2+1
        call em3fftnextn(nphi,nphifft)
        nphi=nphifft
        ntheta=nquad2+1
cc        call prinf('nphi=*',nphi,1)
cc        call prinf('ntheta=*',ntheta,1)
        call e3fgrid(itype,nquad2,nphi,ntheta,rnodes,rwts,nnodes)        
c
c
c
        if(ifprint .ge. 1)
     $     call prinf('=== STEP 6 (eval mp) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 6, adaptive part, evaluate multipole expansions, 
c           or evaluate the potentials and fields directly
c
         do 3252 ibox=1,nboxes
         call d3tgetb(ier,ibox,box,center0,corners0,wlists)
c
         itype=4
         call d3tgetl(ier,ibox,itype,list,nlist,wlists)
         if (nlist .gt. 0) then 
            if (ifprint .ge. 2) then
               call prinf('ibox=*',ibox,1)
               call prinf('list4=*',list,nlist)
            endif
         endif
c
c       ... prune all sourceless boxes
c
         if( box(15) .eq. 0 ) nlist=0
c
c       ... note that lists 3 and 4 are dual
c
c       ... evaluate multipole expansions for all boxes in list 4 
c       ... if source is childless, evaluate directly (if cheaper)
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(jbox,box1,center1,corners1,level1,ifdirect4,level)
C$OMP$PRIVATE(ier,i,j,ptemp,ftemp,cd,ilist) 
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
         do ilist=1,nlist
            jbox=list(ilist)
            call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
            level=box(1)
c
            ifdirect4 = 0
c
            if (box1(15) .lt. (nterms(level)+1)**2/4 .and.
     $         box(15) .lt. (nterms(level)+1)**2/4 ) ifdirect4 = 1
c
            ifdirect4 = 0
c
            if (ifdirect4 .eq. 0) then

            do i=box1(14),box1(14)+box1(15)-1

        call em3mpta3_add(zk,center0,
     $     rmlexp(iaddr(1,ibox)),rmlexp(iaddr(2,ibox)),
     $     nterms(level),
     $     sourcesort(1,i),
     $     ampoleincsort(0,-nterms0,i),bmpoleincsort(0,-nterms0,i),
     $     nterms0,
     $     sphereradsort(i),rnodes,rwts,nphi,ntheta)

            enddo

            else
            
            write(*,*) '***not implemented***'
            stop
c            call emfmm3dpart_direct(zk,box,box1,sourcesort,
c     $         ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
c     $         ifpot,pot,iffld,fld,
c     $         targetsort,ifpottarg,pottarg,iffldtarg,fldtarg)

            endif
        enddo
C$OMP END PARALLEL DO
 3252   continue
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(6)=t2-t1
c

        if(ifprint .ge. 1)
     $     call prinf('=== STEP 7 (eval lo) ====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 7, evaluate local expansions
c       and all fields directly
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,level,npts,nkids,ier)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 6201 ibox=1,nboxes
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
            npts=box(15)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(14)),box(15))
            endif
        endif
c
        if (nkids .eq. 0) then
c
c       ... evaluate local expansions
c       
        level=box(1)
        npts=box(15)
c       
        if (level .ge. 2) then

        do i=box(14),box(14)+box(15)-1

        call em3tata3_add(zk,center0,
     $     rmlexp(iaddr(3,ibox)),rmlexp(iaddr(4,ibox)),
     $     nterms(level),
     $     sourcesort(1,i),
     $     ampoleincsort(0,-nterms0,i),bmpoleincsort(0,-nterms0,i),
     $     nterms0,
     $     sphereradsort(i),rnodes,rwts,nphi,ntheta)

        enddo
        
        endif
c
        endif
c
 6201   continue
C$OMP END PARALLEL DO
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(7)=t2-t1
c
c
 8000   continue
c
c
        if( ifevalloc .eq. 0 ) goto 9000
c 
c
        if(ifprint .ge. 1)
     $     call prinf('=== STEP 8 (direct) =====*',i,0)
        t1=second()
C$        t1=omp_get_wtime()
c
c       ... step 8, evaluate direct interactions 
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,box,center0,corners0,nkids,list,nlist,npts)
C$OMP$PRIVATE(jbox,box1,center1,corners1)
C$OMP$PRIVATE(ier,ilist,itype) 
C$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(1) 
        do 6202 ibox=1,nboxes
c
        call d3tgetb(ier,ibox,box,center0,corners0,wlists)
        call d3tnkids(box,nkids)
c
        if (ifprint .ge. 2) then
           call prinf('ibox=*',ibox,1)
           call prinf('box=*',box,20)
           call prinf('nkids=*',nkids,1)
        endif
c
        if (nkids .eq. 0) then
            npts=box(15)
            if (ifprint .ge. 2) then
               call prinf('npts=*',npts,1)
               call prinf('isource=*',isource(box(14)),box(15))
            endif
        endif
c
        level=box(1)
c
        if (nkids .eq. 0) then
c
c       ... evaluate self interactions
c
        do 6212 j=box(14),box(14)+box(15)-1
        do 6211 i=box(14),box(14)+box(15)-1
        if( i .eq. j ) goto 6211

        call em3mpta3_add(zk,sourcesort(1,i),
     $     ampoleoutsort(0,-nterms0,i),bmpoleoutsort(0,-nterms0,i),
     $     nterms0,
     $     sourcesort(1,j),
     $     ampoleincsort(0,-nterms0,j),bmpoleincsort(0,-nterms0,j),
     $     nterms0,
     $     sphereradsort(j),rnodes,rwts,nphi,ntheta)

 6211   continue
 6212   continue

c
c
c       ... retrieve list #1
c
c       ... evaluate interactions with the nearest neighbours
c
        itype=1
        call d3tgetl(ier,ibox,itype,list,nlist,wlists)
        if (ifprint .ge. 2) call prinf('list1=*',list,nlist)
c
c       ... for all pairs in list #1, 
c       evaluate the potentials and fields directly
c
            do 6203 ilist=1,nlist
               jbox=list(ilist)
               call d3tgetb(ier,jbox,box1,center1,corners1,wlists)
c
c       ... prune all sourceless boxes
c
         if( box1(15) .eq. 0 ) goto 6203
c    
        do 6215 j=box(14),box(14)+box(15)-1
        do 6214 i=box1(14),box1(14)+box1(15)-1

        call em3mpta3_add(zk,sourcesort(1,i),
     $     ampoleoutsort(0,-nterms0,i),bmpoleoutsort(0,-nterms0,i),
     $     nterms0,
     $     sourcesort(1,j),
     $     ampoleincsort(0,-nterms0,j),bmpoleincsort(0,-nterms0,j),
     $     nterms0,
     $     sphereradsort(j),rnodes,rwts,nphi,ntheta)

 6214   continue
 6215   continue

 6203       continue
        endif
c
 6202   continue
C$OMP END PARALLEL DO
c
ccc        call prin2('inside fmm, pot=*',pot,2*nsource)
c
c
        t2=second()
C$        t2=omp_get_wtime()
ccc     call prin2('time=*',t2-t1,1)
        timeinfo(8)=t2-t1
c
 9000   continue
c
ccc        call prinf('=== DOWNWARD PASS COMPLETE ===*',i,0)
c
        if (ifprint .ge. 1) call prin2('timeinfo=*',timeinfo,8)
c       
        if (ifprint .ge. 1) then
        call prinf('nboxes=*',nboxes,1)
        call prinf('nsource=*',nsource,1)
        call prinf('ntarget=*',ntarget,1)
        endif
c       
        return
        end
c
c
c
c
c
        subroutine em3dsphreorder
     $     (nsource,source,sphererad,ampole,bmpole,isource,
     1     sourcesort,sphereradsort,
     $     ampolesort,bmpolesort,nterms)
        implicit real *8 (a-h,o-z)
        dimension source(3,*),sourcesort(3,*),isource(*)
        complex *16 ampole(0:nterms,-nterms:nterms,*)
        complex *16 bmpole(0:nterms,-nterms:nterms,*)
        complex *16 ampolesort(0:nterms,-nterms:nterms,*)
        complex *16 bmpolesort(0:nterms,-nterms:nterms,*)
        dimension sphererad(*),sphereradsort(*)
c       
        do i = 1,nsource
        sourcesort(1,i) = source(1,isource(i))
        sourcesort(2,i) = source(2,isource(i))
        sourcesort(3,i) = source(3,isource(i))
        do n=0,nterms
        do m=-nterms,nterms
        ampolesort(n,m,i) = ampole(n,m,isource(i))
        bmpolesort(n,m,i) = bmpole(n,m,isource(i))
        enddo
        enddo
        sphereradsort(i) = sphererad(isource(i))
        enddo
c
        return
        end
c
c
c
c
c
        subroutine em3dmpolesort
     $     (nsource,isource,mpolesort,mpole,nterms)
        implicit real *8 (a-h,o-z)
        dimension isource(*)
        complex *16 mpole(0:nterms,-nterms:nterms,*)
        complex *16 mpolesort(0:nterms,-nterms:nterms,*)
c
cccc        call prinf('isource=*',isource,nsource)
c        
        do i = 1,nsource
        do n=0,nterms
        do m=-nterms,nterms
        mpole(n,m,isource(i)) = mpolesort(n,m,i)
        enddo
        enddo
        enddo
c
        return
        end
c
c
c
c
