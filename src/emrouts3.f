cc Copyright (C) 2010: Leslie Greengard and Zydrunas Gimbutas
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
c    $Date: 2013-12-20 21:59:08 -0500 (Fri, 20 Dec 2013) $
c    $Revision: 3502 $
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging routines and the beginning of
c       the code for the evaluation of electromagnetic multipoles and
c       local expansions, translation operators and far field signatures
c
c       Maxwell's equations in R^3
c
c
c       Fortran F90 version, allocatable arrays
c
c       NOTE #1: heavy use of allocatable arrays to simplify programming 
c       We use automatic deallocation feature of Fortran F95
c
c       NOTE #2: Translation operators are of order O(p^3),
c       Point and shoot algorithm, via stabilized projection
c
c       NOTE #3: Legendre nodes in theta direction are used in all routines
c       
c       NOTE #4: Spherical Hankel and Bessel functions are NOT scaled
c       in this version, which may cause instability at low frequencies
c
c       NOTE #5: the code has not been extensively tested and may fail for
c       certain point locations (x,y,z)=0, z=0, etc. (but it should be ok)
c
c       NOTE #6:
c
c         There are many definitions of the spherical harmonics,
c         which differ in terms of normalization constants. We
c         adopt the following convention:
c
c         For m>0, we define Y_n^m according to 
c
c         Y_n^m = \sqrt{2n+1} \sqrt{\frac{ (n-m)!}{(n+m)!}} \cdot
c                 P_n^m(\cos \theta)  e^{i m phi} 
c         and
c 
c         Y_n^-m = dconjg( Y_n^m )
c    
c         We omit the Condon-Shortley phase factor (-1)^m in the 
c         definition of Y_n^m for m<0. (This is standard in several
c         communities.)
c
c         We also omit the factor \sqrt{\frac{1}{4 \pi}}, so that
c         the Y_n^m are orthogonal on the unit sphere but not 
c         orthonormal.  (This is also standard in several communities.)
c         More precisely, 
c
c                 \int_S Y_n^m Y_n^m d\Omega = 4 \pi. 
c
c         Using our standard definition, the addition theorem takes 
c         the simple form 
c
c         e^( i k r}/(ikr) = 
c         \sum_n \sum_m  j_n(k|S|) Ylm*(S) h_n(k|T|) Ylm(T)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       This file contains a set of subroutines for the handling of
c       electromagnetic field multipoles. It contains 21 subroutines that
c       are user-callable. Following is a brief description of these
c       subroutines.
c
c
c       em3formmp - form the electromagnetic multipole (ampole,bmpole) 
c           located at the location center due to electric dipole cjvec 
c           and magnetic dipole cmvec located at an arbitrary source 
c
c       em3formta - form the electromagnetic local expansion (ampole,bmpole) 
c           located at the location center due to electric dipole cjvec 
c           and magnetic dipole cmvec located at an arbitrary source 
c
c       em3mpeval - evaluate E and H fields at the location target due
c           to the electromagnetic multipole (ampole,bmpole) 
c           located at the location center
c
c       em3taeval - evaluate E and H fields at the location target due
c           to the electromagnetic local expansion (ampole,bmpole) 
c           located at the location center
c
c
c
c       em3fgrid - form the nodes and weights for quadrature on the unit sphere
c
c       em3sgrid - compressed form for quadrature on the unit sphere
c
c       em3ehgrid - form the E and H field grids at the location center
c           due to electric dipole cjvec and magnetic dipole cmvec located
c           at an arbitrary source
c
c       em3ehformmp - form the outgoing EM-multipole expansion
c           centered at the center due to the E and H field grids
c       
c       em3ehformta - form the incoming EM-multipole expansion
c           centered at the center due to the E and H field grids
c
c       em3mpevaleh - evaluate E and H field grid at due
c           to the outgoing electromagnetic multipole (ampole,bmpole) 
c           located at the location center
c
c       em3taevaleh - evaluate E and H field grid at due
c           to the incoming electromagnetic multipole (ampole,bmpole) 
c           located at the location center
c
c
c       em3mpmp - shift the electromagnetic multipole (ampole,bmpole) 
c           located at the location center to the new location 
c
c       em3mpta - convert the electromagnetic multipole (ampole,bmpole)
c           located at the location center to the electromagnetic local
c           expansion multipole at the new location
c
c       em3tata - shift the electromagnetic local expansion (ampole,bmpole) 
c           located at the location center to the new location
c       
c
c
c       em3mpfar - evaluate the far field signature of the
c           electromagnetic multipole (ampole,bmpole) located at the
c           location center 
c
c       em3mpfareh - evaluate the far field signature grid of the
c           electromagnetic multipole (ampole,bmpole) located at the
c           location center
c
c
c       em3ehfar - form the far field signature due to electric dipole cjvec 
c           and magnetic dipole cmvec located at an arbitrary source 
c
c       em3farfar - shift the center of the outgoing far field signature 
c
c       em3farpll - convert the outgoing far field signature into the incoming
c           far field signature aroun a sufficiently shifted center
c
c       em3pllpll - shift the center of the incoming far field signature 
c
c       em3plleval - evaluate E and H fields at the location xyz due to
c           incoming far field signature located at the location center
c
c
c
c       TO DO...
c
c       em3pllta - convert the incoming far field signature 
c           to the electromagnetic local expansion (ampole,bmpole) 
c           located at the location center
c
c       em3farta - convert the far field signature 
c           to the electromagnetic local expansion (ampole,bmpole) 
c           located at the location center
c       
c       em3fareval - evaluate E and H fields at the location xyz due to
c           outgoing far field signature located at the location center
c       
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Form and evaluate multipole and local expansions 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine em3formmp
     $     (rk,xyz,cjvec,cmvec,npts,center,ampole,bmpole,nterms)
        implicit real *8 (a-h,o-z)
c
c       This subroutine forms the outgoing EM-multipole expansion
c       centered at the origin due to the monochromatic electric dipole
c       cjvec and the monochromatic magnetic dipole cmvec located at the
c       location xyz
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of dipoles in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of outgoing EM-multipole expansion
c        
        complex *16 rk
        dimension xyz(3,1),center(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 cjvec(3,1),cmvec(3,1)
c
        do n=0,nterms
        do m=-n,n
        ampole(n,m)=0
        bmpole(n,m)=0
        enddo
        enddo
c
        do 1200 i=1,npts
        itype=1
        call em3formex
     $     (itype,rk,xyz(1,i),cjvec(1,i),cmvec(1,i),
     $     center,ampole,bmpole,nterms)
c
 1200   continue
c
        do n=1,nterms
        do m=-n,n
c
        ampole(n,m)=-ampole(n,m)*(rk**2) 
        bmpole(n,m)=-bmpole(n,m)*(rk**2) 
c
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine em3formmp_one
     $     (rk,xyz,cjvec,cmvec,center,ampole,bmpole,nterms)
        implicit real *8 (a-h,o-z)
c
c       This subroutine forms the outgoing EM-multipole expansion
c       centered at the origin due to the monochromatic electric dipole
c       cjvec and the monochromatic magnetic dipole cmvec located at the
c       location xyz
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of dipoles in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of outgoing EM-multipole expansion
c        
        complex *16 rk
        dimension xyz(3),center(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 cjvec(3),cmvec(3)
c
        do n=0,nterms
        do m=-n,n
        ampole(n,m)=0
        bmpole(n,m)=0
        enddo
        enddo
c
        itype=1
        call em3formex
     $     (itype,rk,xyz,cjvec,cmvec,center,ampole,bmpole,nterms)
c
        do n=1,nterms
        do m=-n,n
c
        ampole(n,m)=-ampole(n,m)*(rk**2) 
        bmpole(n,m)=-bmpole(n,m)*(rk**2) 
c
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine em3formta
     $     (rk,xyz,cjvec,cmvec,npts,center,ampole,bmpole,nterms)
        implicit real *8 (a-h,o-z)
c
c       This subroutine forms the incoming EM-multipole expansion
c       centered at the center due to the monochromatic electric dipole
c       cjvec and the monochromatic magnetic dipole cmvec located at the
c       location xyz
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of dipoles in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of incoming EM-multipole expansion
c        
        complex *16 rk
        dimension xyz(3,1),center(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 cjvec(3,1),cmvec(3,1)
c
        do n=0,nterms
        do m=-n,n
        ampole(n,m)=0
        bmpole(n,m)=0
        enddo
        enddo
c
        do 1200 i=1,npts
        itype=2
        call em3formex
     $     (itype,rk,xyz(1,i),cjvec(1,i),cmvec(1,i),
     $     center,ampole,bmpole,nterms)
c
 1200   continue
c
        do n=1,nterms
        do m=-n,n
c
        ampole(n,m)=-ampole(n,m)*(rk**2) 
        bmpole(n,m)=-bmpole(n,m)*(rk**2) 
c
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine em3formta_one
     $     (rk,xyz,cjvec,cmvec,center,ampole,bmpole,nterms)
        implicit real *8 (a-h,o-z)
c
c       This subroutine forms the incoming EM-multipole expansion
c       centered at the center due to the monochromatic electric dipole
c       cjvec and the monochromatic magnetic dipole cmvec located at the
c       location xyz
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of dipoles in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of incoming EM-multipole expansion
c        
        complex *16 rk
        dimension xyz(3),center(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 cjvec(3),cmvec(3)
c
        do n=0,nterms
        do m=-n,n
        ampole(n,m)=0
        bmpole(n,m)=0
        enddo
        enddo
c
        itype=2
        call em3formex
     $     (itype,rk,xyz,cjvec,cmvec,center,ampole,bmpole,nterms)
c
        do n=1,nterms
        do m=-n,n
c
        ampole(n,m)=-ampole(n,m)*(rk**2) 
        bmpole(n,m)=-bmpole(n,m)*(rk**2) 
c
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine em3formex
     $     (itype,rk,xyz,cjvec,cmvec,center,ampole,bmpole,nterms)
        implicit real *8 (a-h,o-z)
c
c       This subroutine forms the outgoing/incoming EM-multipole expansion
c       centered at the origin due to the monochromatic electric dipole
c       cjvec and the monochromatic magnetic dipole cmvec located at the
c       location xyz
c
c          Input parameters:
c
c       
c       itype - the type of expansion
c          the expansions are outgoing, for itype=1
c          the expansions are incoming, for itype=2
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of dipoles in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of outgoing EM-multipole expansion
c        
        complex *16 rk,z
        dimension xyz(3),center(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 cjvec(3),cmvec(3)
c
        complex *16, allocatable :: qvals(:), qders(:)
        complex *16, allocatable :: emtqvals(:), emrqvals(:)
c
        real *8, allocatable :: pnm(:,:), dnm(:,:)
        complex *16, allocatable :: pxnm2(:,:,:)
        complex *16, allocatable :: punm2(:,:,:)
        complex *16, allocatable :: pxnm3(:,:,:)
        complex *16, allocatable :: punm3(:,:,:)
c
        complex *16, allocatable :: ynm(:,:)
        complex *16, allocatable :: xnm3(:,:,:)
        complex *16, allocatable :: unm3(:,:,:)
c
        complex *16, allocatable :: ephi(:)
c
        complex *16 cd,cp,cr
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        allocate( qvals(0:nterms) )
        allocate( qders(0:nterms) )
        allocate( emtqvals(0:nterms) )
        allocate( emrqvals(0:nterms) )
        allocate( pnm(0:nterms,0:nterms) )
        allocate( dnm(0:nterms,0:nterms) )
        allocate( pxnm2(2,0:nterms,0:nterms) )
        allocate( punm2(2,0:nterms,0:nterms) ) 
        allocate( pxnm3(3,0:nterms,0:nterms) )
        allocate( punm3(3,0:nterms,0:nterms) ) 
        allocate( ynm(0:nterms,-nterms:nterms) )
        allocate( xnm3(3,0:nterms,-nterms:nterms) )
        allocate( unm3(3,0:nterms,-nterms:nterms) )
        allocate( ephi(-nterms:nterms) )
c
c
        dx=xyz(1)-center(1)
        dy=xyz(2)-center(2)
        dz=xyz(3)-center(3)
c
        r=sqrt(dx*dx+dy*dy+dz*dz)        
c
        if( r .eq. 0 ) then
          rx=1
          ry=0
          rz=0
        else
          rx=dx/r
          ry=dy/r
          rz=dz/r
        endif
c
        z=rk*r
c
        if( itype .eq. 1 ) then
c
c       ... evaluate jvals, emtjvals, emrjvals
c
        call emjevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c
        if( itype .eq. 2 ) then
c
c       ... evaluate hvals, emthvals, emrhvals
c
        call emhevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c
c
c       ... evaluate xnm3, unm3, ynm, ephi
c        
        call em3phi(dx,dy,phi)
        do m=-nterms,nterms
        ephi(m)=exp(ima*m*phi)
        enddo
c
        call em3ctheta(dz,r,costheta)
        call xulgndr(nterms,costheta,pnm,dnm,pxnm2,punm2)
c
        theta=acos(costheta)
        call sph2cart(rx,ry,rz,phi,theta,(nterms+1)**2,pxnm2,pxnm3)
        call sph2cart(rx,ry,rz,phi,theta,(nterms+1)**2,punm2,punm3)
c
        do n=0,nterms
        do m=0,n
        ynm(n,m)=pnm(n,m)*ephi(m) 
        xnm3(1,n,m)=pxnm3(1,n,m)*ephi(m) 
        xnm3(2,n,m)=pxnm3(2,n,m)*ephi(m) 
        xnm3(3,n,m)=pxnm3(3,n,m)*ephi(m) 
        unm3(1,n,m)=punm3(1,n,m)*ephi(m) 
        unm3(2,n,m)=punm3(2,n,m)*ephi(m) 
        unm3(3,n,m)=punm3(3,n,m)*ephi(m) 
        enddo
        do m=-n,-1
        ynm(n,m)=conjg(ynm(n,-m)) 
        xnm3(1,n,m)=conjg(xnm3(1,n,-m)) 
        xnm3(2,n,m)=conjg(xnm3(2,n,-m)) 
        xnm3(3,n,m)=conjg(xnm3(3,n,-m)) 
        unm3(1,n,m)=conjg(unm3(1,n,-m)) 
        unm3(2,n,m)=conjg(unm3(2,n,-m)) 
        unm3(3,n,m)=conjg(unm3(3,n,-m)) 
        enddo
        enddo
c
c
c
c       ... form outgoing a and b multipole expansions
c
c
        if( 1 .eq. 2 ) then
c
        do n=0,nterms
        do m=-n,n
        ampole(n,m)=0
        bmpole(n,m)=0
        enddo
        enddo
c
        endif
c
c
        do n=1,nterms
        do m=-n,n
        ampole(n,m)=ampole(n,m)+
     $     (
     $     cjvec(1)*conjg(xnm3(1,n,m))+
     $     cjvec(2)*conjg(xnm3(2,n,m))+
     $     cjvec(3)*conjg(xnm3(3,n,m))
     $     )*qvals(n)
        bmpole(n,m)=bmpole(n,m)+
     $     (
     $     cjvec(1)*conjg(unm3(1,n,m))+
     $     cjvec(2)*conjg(unm3(2,n,m))+
     $     cjvec(3)*conjg(unm3(3,n,m))
     $     )*emtqvals(n)*(ima)
        bmpole(n,m)=bmpole(n,m)+
     $     (
     $     cjvec(1)*rx*conjg(ynm(n,m))+
     $     cjvec(2)*ry*conjg(ynm(n,m))+
     $     cjvec(3)*rz*conjg(ynm(n,m))
     $     )*emrqvals(n)*(ima) 
        enddo
        enddo
c
c
        do n=1,nterms
        do m=-n,n
        bmpole(n,m)=bmpole(n,m)+
     $     (
     $     cmvec(1)*conjg(xnm3(1,n,m))+
     $     cmvec(2)*conjg(xnm3(2,n,m))+
     $     cmvec(3)*conjg(xnm3(3,n,m))
     $     )*qvals(n) *(-1)
        ampole(n,m)=ampole(n,m)+
     $     (
     $     cmvec(1)*conjg(unm3(1,n,m))+
     $     cmvec(2)*conjg(unm3(2,n,m))+
     $     cmvec(3)*conjg(unm3(3,n,m))
     $     )*emtqvals(n)*(ima)
        ampole(n,m)=ampole(n,m)+
     $     (
     $     cmvec(1)*rx*conjg(ynm(n,m))+
     $     cmvec(2)*ry*conjg(ynm(n,m))+
     $     cmvec(3)*rz*conjg(ynm(n,m))
     $     )*emrqvals(n)*(ima) 
        enddo
        enddo
c
c
        if( 1 .eq. 2 ) then
c
        do n=1,nterms
        do m=-n,n
c
        ampole(n,m)=-ampole(n,m)*(rk**2) 
        bmpole(n,m)=-bmpole(n,m)*(rk**2) 
c
        enddo
        enddo
c
        endif
c
ccc        call prin2('ampole=*',ampole,(nterms+1)*(2*nterms+1)*2)
ccc        call prin2('bmpole=*',bmpole,(nterms+1)*(2*nterms+1)*2)
c
ccc        call prinm1(ampole,nterms,2)
ccc        call prinm1(bmpole,nterms,2)
c
ccc        call prinml2('ampole=*',ampole,nterms,nterms)
ccc        call prinml2('bmpole=*',bmpole,nterms,nterms)
c
        return
        end
c
c
c
c
c
        subroutine em3mpeval
     $     (rk,center,ampole,bmpole,nterms,xyz,evec,hvec)
        implicit real *8 (a-h,o-z)
c        
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the outgoing EM-multipole expansion (ampole,bmpole) located
c       at the center
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       center - the location of EM-multipole expansion in R^3
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of outgoing EM-multipole expansion
c       nterms - the number of terms in EM-multipole
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        complex *16 rk,z
        dimension xyz(3),center(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 evec(3),hvec(3)
c
        itype=1
c
        call em3exeval
     $     (itype,rk,center,ampole,bmpole,nterms,xyz,evec,hvec)
 1200   continue
c
c
        return
        end
c
c
c
c
c
        subroutine em3taeval
     $     (rk,center,ampole,bmpole,nterms,xyz,evec,hvec)
        implicit real *8 (a-h,o-z)
c        
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the incoming EM-multipole expansion (ampole,bmpole) located
c       at the center
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       center - the location of EM-multipole expansion in R^3
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of incoming EM-multipole expansion
c       nterms - the number of terms in EM-multipole
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        complex *16 rk,z
        dimension xyz(3),center(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 evec(3),hvec(3)
c
        itype=2

        call em3exeval
     $     (itype,rk,center,ampole,bmpole,nterms,xyz,evec,hvec)
c
c
        return
        end
c
c
c
c
c
        subroutine em3exeval
     $     (itype,rk,center,ampole,bmpole,nterms,xyz,evec,hvec)
        implicit real *8 (a-h,o-z)
c        
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the EM-multipole expansion (ampole,bmpole) located at the center
c
c          Input parameters:
c
c       itype - the type of expansion
c          the expansions are outgoing, for itype=1
c          the expansions are incoming, for itype=2
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       center - the location of EM-multipole expansion in R^3
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of EM-multipole expansion
c       nterms - the number of terms in EM-multipole
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        complex *16 rk,z
        dimension xyz(3),center(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 evec(3),hvec(3)
c
        complex *16, allocatable :: qvals(:), qders(:)
        complex *16, allocatable :: emtqvals(:), emrqvals(:)
c
        real *8, allocatable :: pnm(:,:), dnm(:,:)
        complex *16, allocatable :: pxnm2(:,:,:)
        complex *16, allocatable :: punm2(:,:,:)
        complex *16, allocatable :: pxnm3(:,:,:)
        complex *16, allocatable :: punm3(:,:,:)
c
        complex *16, allocatable :: ynm(:,:)
        complex *16, allocatable :: xnm3(:,:,:)
        complex *16, allocatable :: unm3(:,:,:)
c
        complex *16, allocatable :: ephi(:)
c
        complex *16 cd,cp,cr
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        allocate( qvals(0:nterms) )
        allocate( qders(0:nterms) )
        allocate( emtqvals(0:nterms) )
        allocate( emrqvals(0:nterms) )
        allocate( pnm(0:nterms,0:nterms) )
        allocate( dnm(0:nterms,0:nterms) )
        allocate( pxnm2(2,0:nterms,0:nterms) )
        allocate( punm2(2,0:nterms,0:nterms) ) 
        allocate( pxnm3(3,0:nterms,0:nterms) )
        allocate( punm3(3,0:nterms,0:nterms) ) 
        allocate( ynm(0:nterms,-nterms:nterms) )
        allocate( xnm3(3,0:nterms,-nterms:nterms) )
        allocate( unm3(3,0:nterms,-nterms:nterms) )
        allocate( ephi(-nterms:nterms) )
c
        dx=xyz(1)-center(1)
        dy=xyz(2)-center(2)
        dz=xyz(3)-center(3)
c
        r=sqrt(dx*dx+dy*dy+dz*dz)        
c
        if( r .eq. 0 ) then
          rx=1
          ry=0
          rz=0
        else
          rx=dx/r
          ry=dy/r
          rz=dz/r
        endif
c
        z=rk*r
c
        if( itype .eq. 1 ) then
c
c       ... evaluate hvals, emthvals, emrhvals
c
        call emhevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c
        if( itype .eq. 2 ) then
c
c       ... evaluate jvals, emtjvals, emrjvals
c
        call emjevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c
c       ... evaluate xnm3, unm3, ynm, ephi
c
        call em3phi(dx,dy,phi)
        do m=-nterms,nterms
        ephi(m)=exp(ima*m*phi)
        enddo
c
        call em3ctheta(dz,r,costheta)
        call xulgndr(nterms,costheta,pnm,dnm,pxnm2,punm2)
c
        theta=acos(costheta)
        call sph2cart(rx,ry,rz,phi,theta,(nterms+1)**2,pxnm2,pxnm3)
        call sph2cart(rx,ry,rz,phi,theta,(nterms+1)**2,punm2,punm3)
c
        do n=0,nterms
        do m=-nterms,nterms
        ynm(n,m)=0
        xnm3(1,n,m)=0
        xnm3(2,n,m)=0
        xnm3(3,n,m)=0
        unm3(1,n,m)=0
        unm3(2,n,m)=0
        unm3(3,n,m)=0
        enddo
        enddo
c
        do n=0,nterms
        do m=0,n
        ynm(n,m)=pnm(n,m)*ephi(m) 
        xnm3(1,n,m)=pxnm3(1,n,m)*ephi(m) 
        xnm3(2,n,m)=pxnm3(2,n,m)*ephi(m) 
        xnm3(3,n,m)=pxnm3(3,n,m)*ephi(m) 
        unm3(1,n,m)=punm3(1,n,m)*ephi(m) 
        unm3(2,n,m)=punm3(2,n,m)*ephi(m) 
        unm3(3,n,m)=punm3(3,n,m)*ephi(m) 
        enddo
        do m=-n,-1
        ynm(n,m)=conjg(ynm(n,-m)) 
        xnm3(1,n,m)=conjg(xnm3(1,n,-m)) 
        xnm3(2,n,m)=conjg(xnm3(2,n,-m)) 
        xnm3(3,n,m)=conjg(xnm3(3,n,-m)) 
        unm3(1,n,m)=conjg(unm3(1,n,-m)) 
        unm3(2,n,m)=conjg(unm3(2,n,-m)) 
        unm3(3,n,m)=conjg(unm3(3,n,-m)) 
        enddo
        enddo
c
        if( itype .eq. 2 ) then
c        write(*,*) dx,dy,dz,r
c        write(*,*) rx,ry,rz
c        pause
c        call prin2('qvals=*',qvals,2*(nterms+1))
c        call prin2('emrqvals=*',emrqvals,2*(nterms+1))
c        call prin2('emtqvals=*',emtqvals,2*(nterms+1))
c        pause
c        call prin2('pxnm2=*',pxnm2,2*2*(nterms+1)*(nterms+1))
c        call prin2('punm2=*',punm2,2*2*(nterms+1)*(nterms+1))
c        pause
c        call prin2('xnm3=*',xnm3,2*3*(nterms+1)*(2*nterms+1))
c        call prin2('unm3=*',unm3,2*3*(nterms+1)*(2*nterms+1))
c        pause
        endif
c
        evec(1)=0
        evec(2)=0
        evec(3)=0
c
        hvec(1)=0
        hvec(2)=0
        hvec(3)=0
c
c
c       ... evaluate E and H fields via outgoing a and b multipole expansions
c
c
        do n=1,nterms
        do m=-n,n
c
        cd=qvals(n)
        cp=emtqvals(n)
        cr=emrqvals(n)
c
        evec(1)=evec(1)+ampole(n,m)*cd*xnm3(1,n,m)
        evec(2)=evec(2)+ampole(n,m)*cd*xnm3(2,n,m)
        evec(3)=evec(3)+ampole(n,m)*cd*xnm3(3,n,m)
c
        evec(1)=evec(1)+bmpole(n,m)*cp*unm3(1,n,m)*(-ima)
        evec(2)=evec(2)+bmpole(n,m)*cp*unm3(2,n,m)*(-ima)
        evec(3)=evec(3)+bmpole(n,m)*cp*unm3(3,n,m)*(-ima)
c
        evec(1)=evec(1)+bmpole(n,m)*cr*ynm(n,m)*rx*(-ima)
        evec(2)=evec(2)+bmpole(n,m)*cr*ynm(n,m)*ry*(-ima)
        evec(3)=evec(3)+bmpole(n,m)*cr*ynm(n,m)*rz*(-ima)
c
        hvec(1)=hvec(1)+bmpole(n,m)*cd*xnm3(1,n,m)
        hvec(2)=hvec(2)+bmpole(n,m)*cd*xnm3(2,n,m)
        hvec(3)=hvec(3)+bmpole(n,m)*cd*xnm3(3,n,m)
c
        hvec(1)=hvec(1)+ampole(n,m)*cp*unm3(1,n,m)*(+ima)
        hvec(2)=hvec(2)+ampole(n,m)*cp*unm3(2,n,m)*(+ima)
        hvec(3)=hvec(3)+ampole(n,m)*cp*unm3(3,n,m)*(+ima)
c
        hvec(1)=hvec(1)+ampole(n,m)*cr*ynm(n,m)*rx*(+ima)
        hvec(2)=hvec(2)+ampole(n,m)*cr*ynm(n,m)*ry*(+ima)
        hvec(3)=hvec(3)+ampole(n,m)*cr*ynm(n,m)*rz*(+ima)
c
        enddo
        enddo
c
        return
        end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Far field signature of multipole expansion
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine em3mpfar
     $     (rk,center,ampole,bmpole,nterms,dir,efar,hfar)
        implicit real *8 (a-h,o-z)
c        
        complex *16 rk,z
        dimension center(3),dir(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 evec(3),hvec(3)
        complex *16 efar(3),hfar(3)
c
        complex *16, allocatable :: hvals(:), hders(:)
        complex *16, allocatable :: emthvals(:), emrhvals(:)
c
        real *8, allocatable :: pnm(:,:), dnm(:,:)
        complex *16, allocatable :: pxnm2(:,:,:)
        complex *16, allocatable :: punm2(:,:,:)
        complex *16, allocatable :: pxnm3(:,:,:)
        complex *16, allocatable :: punm3(:,:,:)
c
        complex *16, allocatable :: ynm(:,:)
        complex *16, allocatable :: xnm3(:,:,:)
        complex *16, allocatable :: unm3(:,:,:)
c
        complex *16, allocatable :: ephi(:)
c
        dimension far(3)
        complex *16 cfar
c
        complex *16 cd,cp,cr
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        allocate( hvals(0:nterms) )
        allocate( hders(0:nterms) )
        allocate( emthvals(0:nterms) )
        allocate( emrhvals(0:nterms) )
        allocate( pnm(0:nterms,0:nterms) )
        allocate( dnm(0:nterms,0:nterms) )
        allocate( pxnm2(2,0:nterms,0:nterms) )
        allocate( punm2(2,0:nterms,0:nterms) ) 
        allocate( pxnm3(3,0:nterms,0:nterms) )
        allocate( punm3(3,0:nterms,0:nterms) ) 
        allocate( ynm(0:nterms,-nterms:nterms) )
        allocate( xnm3(3,0:nterms,-nterms:nterms) )
        allocate( unm3(3,0:nterms,-nterms:nterms) )
        allocate( ephi(-nterms:nterms) )
c
c
        dx=dir(1)
        dy=dir(2)
        dz=dir(3)
c
        r=sqrt(dx*dx+dy*dy+dz*dz)        
c
        if( r .eq. 0 ) then
          rx=1
          ry=0
          rz=0
        else
          rx=dx/r
          ry=dy/r
          rz=dz/r
        endif
c
        far(1)=rx
        far(2)=ry
        far(3)=rz
c
c       ... evaluate hvals, emthvals, emrhvals
c
        z=rk*r
        call emhevalrt(nterms,z,hvals,hders,emthvals,emrhvals)
c
c       ... evaluate xnm3, unm3, ynm, ephi
c       
        call em3phi(dx,dy,phi)
        do m=-nterms,nterms
        ephi(m)=exp(ima*m*phi)
        enddo
c
        call em3ctheta(dz,r,costheta)
        call xulgndr(nterms,costheta,pnm,dnm,pxnm2,punm2)
c
        theta=acos(costheta)
        call sph2cart(rx,ry,rz,phi,theta,(nterms+1)**2,pxnm2,pxnm3)
        call sph2cart(rx,ry,rz,phi,theta,(nterms+1)**2,punm2,punm3)
c
        do n=0,nterms
        do m=-nterms,nterms
        ynm(n,m)=0
        xnm3(1,n,m)=0
        xnm3(2,n,m)=0
        xnm3(3,n,m)=0
        unm3(1,n,m)=0
        unm3(2,n,m)=0
        unm3(3,n,m)=0
        enddo
        enddo
c
        do n=0,nterms
        do m=0,n
        ynm(n,m)=pnm(n,m)*ephi(m) 
        xnm3(1,n,m)=pxnm3(1,n,m)*ephi(m) 
        xnm3(2,n,m)=pxnm3(2,n,m)*ephi(m) 
        xnm3(3,n,m)=pxnm3(3,n,m)*ephi(m) 
        unm3(1,n,m)=punm3(1,n,m)*ephi(m) 
        unm3(2,n,m)=punm3(2,n,m)*ephi(m) 
        unm3(3,n,m)=punm3(3,n,m)*ephi(m) 
        enddo
        do m=-n,-1
        ynm(n,m)=conjg(ynm(n,-m)) 
        xnm3(1,n,m)=conjg(xnm3(1,n,-m)) 
        xnm3(2,n,m)=conjg(xnm3(2,n,-m)) 
        xnm3(3,n,m)=conjg(xnm3(3,n,-m)) 
        unm3(1,n,m)=conjg(unm3(1,n,-m)) 
        unm3(2,n,m)=conjg(unm3(2,n,-m)) 
        unm3(3,n,m)=conjg(unm3(3,n,-m)) 
        enddo
        enddo
c
c
        evec(1)=0
        evec(2)=0
        evec(3)=0
c
        hvec(1)=0
        hvec(2)=0
        hvec(3)=0
c
c
c       ... evaluate E and H fields via outgoing a and b multipole expansions
c
c
        do n=1,nterms
        do m=-n,n
c
ccc        cd=hvals(n)
ccc        cp=emthvals(n)
ccc        cr=emrhvals(n)
c
        cd=(-ima)**(n+1) / rk
        cp=ima*(-ima)**(n+1) / rk 
        cr=0
c
        evec(1)=evec(1)+ampole(n,m)*cd*xnm3(1,n,m)
        evec(2)=evec(2)+ampole(n,m)*cd*xnm3(2,n,m)
        evec(3)=evec(3)+ampole(n,m)*cd*xnm3(3,n,m)
c
        evec(1)=evec(1)+bmpole(n,m)*cp*unm3(1,n,m)*(-ima)
        evec(2)=evec(2)+bmpole(n,m)*cp*unm3(2,n,m)*(-ima)
        evec(3)=evec(3)+bmpole(n,m)*cp*unm3(3,n,m)*(-ima)
c
        hvec(1)=hvec(1)+bmpole(n,m)*cd*xnm3(1,n,m)
        hvec(2)=hvec(2)+bmpole(n,m)*cd*xnm3(2,n,m)
        hvec(3)=hvec(3)+bmpole(n,m)*cd*xnm3(3,n,m)
c
        hvec(1)=hvec(1)+ampole(n,m)*cp*unm3(1,n,m)*(+ima)
        hvec(2)=hvec(2)+ampole(n,m)*cp*unm3(2,n,m)*(+ima)
        hvec(3)=hvec(3)+ampole(n,m)*cp*unm3(3,n,m)*(+ima)
c
        enddo
        enddo
c
c
        d=center(1)*far(1)+center(2)*far(2)+center(3)*far(3)
        cfar=exp(-ima*rk*d)
c
        efar(1)=evec(1)*cfar
        efar(2)=evec(2)*cfar
        efar(3)=evec(3)*cfar
c
        hfar(1)=hvec(1)*cfar
        hfar(2)=hvec(2)*cfar
        hfar(3)=hvec(3)*cfar
c
        return
        end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       E,H grid evaluation routines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine e3fgrid
     $     (itype,nterms,nphi,ntheta,rnodes,weights,nnodes)
        implicit real *8 (a-h,o-z)
        dimension rnodes(3,1),weights(1)
ccc        dimension ts(10 000), ws(10 000)
        real *8, allocatable :: ts(:), ws(:)
c
        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )

        done=1
        pi=4*atan(done)
c
c       ... construct the Gaussian nodes and weights on the interval [-1,1]
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
        kk=0
        do 1400 i=1,ntheta
c
        z=ts(i)
        r=sqrt(1-z*z)
c
c       ... construct FFT compatible angular quadratures
c
        do 1200 j=1,nphi
c
        phi=(j-1)*(2*pi/nphi)
        x=r*cos(phi)
        y=r*sin(phi)
        kk=kk+1
        rnodes(1,kk)=x
        rnodes(2,kk)=y
        rnodes(3,kk)=z
        weights(kk)=ws(i)*(2*pi/nphi)
 1200   continue
 1400   continue
c
        nnodes=kk
c       
        return
c
c       ... test the quadrature nodes and weights
c
ccc        do kk=1,nnodes
ccc        write(23,*) rnodes(1,kk),rnodes(2,kk),rnodes(3,kk)
ccc        enddo
c
        d=0
c
        do kk=1,nnodes
        d=d+weights(kk)
        enddo
c
        call prin2('sum of weights=*',d,1)
c
        done=1
        pi=4*atan(done)
c
        call prin2('sum of weights/4/pi=*',d/4/pi,1)
c
        return
        end
c
c
c
c
c        
        subroutine e3sgrid(itype,nterms,nphi,ntheta,rtheta,wtheta)
        implicit real *8 (a-h,o-z)
        dimension rtheta(3,1),rwhts(1)
c
c       ... construct the Gaussian nodes and weights on the interval [-1,1]
c
        ifwhts=1
        call legewhts(ntheta,rtheta,wtheta,ifwhts)
c
        return
        end
c
c
c
c
c        
        subroutine em3ehgrid(rk,xyz,cjvec,cmvec,
     $     gridc,radius,egrid,hgrid,rnodes,weights,nphi,ntheta)
        implicit real *8 (a-h,o-z)
c        
        complex *16 rk,z
        dimension xyz(3)
        dimension gridc(3)
c
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c
        dimension target(3)
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
c
        nnodes=nphi*ntheta
c
        do 1400 itheta=1,ntheta
        do 1300 iphi=1,nphi
c       
        target(1)=gridc(1)+radius*rnodes(1,iphi,itheta)
        target(2)=gridc(2)+radius*rnodes(2,iphi,itheta)
        target(3)=gridc(3)+radius*rnodes(3,iphi,itheta)
c
        call dipole3emt(rk,xyz,target,cjvec,cmvec,
     $     egrid(1,iphi,itheta),hgrid(1,iphi,itheta))
c
 1300   continue
 1400   continue
c
        return
        end
c
c
c
c
c
        subroutine em3ehformmp
     $     (rk,gridc,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center,ampole,bmpole,nterms)
        implicit real *8 (a-h,o-z)
c
c       This subroutine forms the outgoing EM-multipole expansion
c       centered at the center due to the E and H field grids
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of dipoles in R^3
c       radius - the radius of E and H grids
c       egrid (complex *16)(3,nphi,ntheta) - E field grid
c       hgrid (complex *16)(3,nphi,ntheta) - H field grid
c       rnodes (real *8)(3,nphi,ntheta) - the grid node coordinates in R^3
c       weights (real *8)(nphi,ntheta) - the grid weights
c       nphi, ntheta - the number of grid points in each direction
c       center - the center of EM-multipole expansion
c
c          Output parameters:
c
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of EM-multipole expansion
c        
        complex *16 rk,z
        dimension center(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        dimension gridc(3)
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c
        itype=1
        call em3ehformex_fast
     $     (itype,rk,gridc,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center,ampole,bmpole,nterms)
c
        return
        end
c
c
c
c
c
        subroutine em3ehformmp2
     $     (rk,gridc,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center,ampole,bmpole,nterms)
        implicit real *8 (a-h,o-z)
c
c       This subroutine forms the outgoing EM-multipole expansion
c       centered at the center due to the E and H field grids
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of dipoles in R^3
c       radius - the radius of E and H grids
c       egrid (complex *16)(3,nphi,ntheta) - E field grid
c       hgrid (complex *16)(3,nphi,ntheta) - H field grid
c       rnodes (real *8)(3,nphi,ntheta) - the grid node coordinates in R^3
c       weights (real *8)(nphi,ntheta) - the grid weights
c       nphi, ntheta - the number of grid points in each direction
c       center - the center of EM-multipole expansion
c
c          Output parameters:
c
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of EM-multipole expansion
c        
        complex *16 rk,z
        dimension center(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        dimension gridc(3)
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c
        itype=1
        call em3ehformex_fast2
     $     (itype,rk,gridc,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center,ampole,bmpole,nterms)
c
        return
        end
c
c
c
c
c
        subroutine em3ehformmp2_trunc
     $     (rk,gridc,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center,ampole,bmpole,nterms,nterms1)
        implicit real *8 (a-h,o-z)
c
c       This subroutine forms the outgoing EM-multipole expansion
c       centered at the center due to the E and H field grids
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of dipoles in R^3
c       radius - the radius of E and H grids
c       egrid (complex *16)(3,nphi,ntheta) - E field grid
c       hgrid (complex *16)(3,nphi,ntheta) - H field grid
c       rnodes (real *8)(3,nphi,ntheta) - the grid node coordinates in R^3
c       weights (real *8)(nphi,ntheta) - the grid weights
c       nphi, ntheta - the number of grid points in each direction
c       center - the center of EM-multipole expansion
c
c          Output parameters:
c
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of EM-multipole expansion
c        
        complex *16 rk,z
        dimension center(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        dimension gridc(3)
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c
        itype=1
        call em3ehformex_fast2_trunc
     $     (itype,rk,gridc,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center,ampole,bmpole,nterms,nterms1)
c
        return
        end
c
c
c
c
c
        subroutine em3ehformta
     $     (rk,gridc,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center,ampole,bmpole,nterms)
        implicit real *8 (a-h,o-z)
c
c       This subroutine forms the incoming EM-multipole expansion
c       centered at the center due to the E and H field grids
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of dipoles in R^3
c       radius - the radius of E and H grids
c       egrid (complex *16)(3,nphi,ntheta) - E field grid
c       hgrid (complex *16)(3,nphi,ntheta) - H field grid
c       rnodes (real *8)(3,nphi,ntheta) - the grid node coordinates in R^3
c       weights (real *8)(nphi,ntheta) - the grid weights
c       nphi, ntheta - the number of grid points in each direction
c       center - the center of EM-multipole expansion
c
c          Output parameters:
c
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of EM-multipole expansion
c        
        complex *16 rk,z
        dimension center(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        dimension gridc(3)
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c
        itype=2
        call em3ehformex_fast
     $     (itype,rk,gridc,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center,ampole,bmpole,nterms)

        return
        end
c
c
c
c
c
        subroutine em3ehformta2
     $     (rk,gridc,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center,ampole,bmpole,nterms)
        implicit real *8 (a-h,o-z)
c
c       This subroutine forms the incoming EM-multipole expansion
c       centered at the center due to the E and H field grids
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of dipoles in R^3
c       radius - the radius of E and H grids
c       egrid (complex *16)(3,nphi,ntheta) - E field grid
c       hgrid (complex *16)(3,nphi,ntheta) - H field grid
c       rnodes (real *8)(3,nphi,ntheta) - the grid node coordinates in R^3
c       weights (real *8)(nphi,ntheta) - the grid weights
c       nphi, ntheta - the number of grid points in each direction
c       center - the center of EM-multipole expansion
c
c          Output parameters:
c
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of EM-multipole expansion
c        
        complex *16 rk,z
        dimension center(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        dimension gridc(3)
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c
        itype=2
        call em3ehformex_fast2
     $     (itype,rk,gridc,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center,ampole,bmpole,nterms)

        return
        end
c
c
c
c
c
        subroutine em3ehformex_slow
     $     (itype,rk,gridc,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center,ampole,bmpole,nterms)
        implicit real *8 (a-h,o-z)
c
c       This subroutine forms the outgoing/incoming EM-multipole expansion
c       centered at the center due to the E and H field grids
c
c       slow algorithm, do not separate phi and theta directions, accumulate 
c       all sums directly
c
c          Input parameters:
c
c       
c       itype - the type of expansion
c          the expansions are outgoing, for itype=1
c          the expansions are incoming, for itype=2
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of dipoles in R^3
c       radius - the radius of E and H grids
c       egrid (complex *16)(3,nphi,ntheta) - E field grid
c       hgrid (complex *16)(3,nphi,ntheta) - H field grid
c       rnodes (real *8)(3,nphi,ntheta) - the grid node coordinates in R^3
c       weights (real *8)(nphi,ntheta) - the grid weights
c       nphi, ntheta - the number of grid points in each direction
c       center - the center of EM-multipole expansion
c
c          Output parameters:
c
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of EM-multipole expansion
c        
        complex *16 rk,z
        dimension xyz(3),center(3)
        dimension gridc(3),target(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c
        complex *16, allocatable :: qvals(:), qders(:)
        complex *16, allocatable :: emtqvals(:), emrqvals(:)
c
        real *8, allocatable :: pnm(:,:), dnm(:,:)
        complex *16, allocatable :: pxnm2(:,:,:)
        complex *16, allocatable :: punm2(:,:,:)
        complex *16, allocatable :: pxnm3(:,:,:)
        complex *16, allocatable :: punm3(:,:,:)
c
        complex *16, allocatable :: ynm(:,:)
        complex *16, allocatable :: xnm3(:,:,:)
        complex *16, allocatable :: unm3(:,:,:)
c
        complex *16, allocatable :: ephi(:)
c
        complex *16 cd,cp,cr
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        allocate( qvals(0:nterms) )
        allocate( qders(0:nterms) )
        allocate( emtqvals(0:nterms) )
        allocate( emrqvals(0:nterms) )
        allocate( pnm(0:nterms,0:nterms) )
        allocate( dnm(0:nterms,0:nterms) )
        allocate( pxnm2(2,0:nterms,0:nterms) )
        allocate( punm2(2,0:nterms,0:nterms) ) 
        allocate( pxnm3(3,0:nterms,0:nterms) )
        allocate( punm3(3,0:nterms,0:nterms) ) 
        allocate( ynm(0:nterms,-nterms:nterms) )
        allocate( xnm3(3,0:nterms,-nterms:nterms) )
        allocate( unm3(3,0:nterms,-nterms:nterms) )
        allocate( ephi(-nterms:nterms) )
c
c
        do n=0,nterms
        do m=-n,n
        ampole(n,m)=0
        bmpole(n,m)=0
        enddo
        enddo
c
c
        done=1
        pi=4*atan(done)
        nnodes=nphi*ntheta
c
c       
c       Note, in this routine, gridc should be the same as center
c
        do 1400 itheta=1,ntheta
        do 1300 iphi=1,nphi
c       
        target(1)=gridc(1)+radius*rnodes(1,iphi,itheta)
        target(2)=gridc(2)+radius*rnodes(2,iphi,itheta)
        target(3)=gridc(3)+radius*rnodes(3,iphi,itheta)
c
        dx=target(1)-center(1)
        dy=target(2)-center(2)
        dz=target(3)-center(3)
c
        r=sqrt(dx*dx+dy*dy+dz*dz)        
c
        if( r .eq. 0 ) then
          rx=1
          ry=0
          rz=0
        else
          rx=dx/r
          ry=dy/r
          rz=dz/r
        endif

c
c       ... evaluate xnm3, unm3, ynm, ephi
c
        call em3phi(dx,dy,phi)
        do m=-nterms,nterms
        ephi(m)=exp(ima*m*phi)
        enddo
c
        call em3ctheta(dz,r,costheta)
        call xulgndr(nterms,costheta,pnm,dnm,pxnm2,punm2)
c
        theta=acos(costheta)
        call sph2cart(rx,ry,rz,phi,theta,(nterms+1)**2,pxnm2,pxnm3)
        call sph2cart(rx,ry,rz,phi,theta,(nterms+1)**2,punm2,punm3)
c
        do n=0,nterms
        do m=0,n
        ynm(n,m)=pnm(n,m)*ephi(m) 
        xnm3(1,n,m)=pxnm3(1,n,m)*ephi(m) 
        xnm3(2,n,m)=pxnm3(2,n,m)*ephi(m) 
        xnm3(3,n,m)=pxnm3(3,n,m)*ephi(m) 
        unm3(1,n,m)=punm3(1,n,m)*ephi(m) 
        unm3(2,n,m)=punm3(2,n,m)*ephi(m) 
        unm3(3,n,m)=punm3(3,n,m)*ephi(m) 
        enddo
        do m=-n,-1
        ynm(n,m)=conjg(ynm(n,-m)) 
        xnm3(1,n,m)=conjg(xnm3(1,n,-m)) 
        xnm3(2,n,m)=conjg(xnm3(2,n,-m)) 
        xnm3(3,n,m)=conjg(xnm3(3,n,-m)) 
        unm3(1,n,m)=conjg(unm3(1,n,-m)) 
        unm3(2,n,m)=conjg(unm3(2,n,-m)) 
        unm3(3,n,m)=conjg(unm3(3,n,-m)) 
        enddo
        enddo
c
c
c
c       ... accumulate a and b multipole expansions
c
c
        do n=1,nterms
        do m=-n,n
        ampole(n,m)=ampole(n,m)+weights(iphi,itheta)*
     $     (
     $     egrid(1,iphi,itheta)*conjg(xnm3(1,n,m))+
     $     egrid(2,iphi,itheta)*conjg(xnm3(2,n,m))+
     $     egrid(3,iphi,itheta)*conjg(xnm3(3,n,m))
     $     )
        bmpole(n,m)=bmpole(n,m)+weights(iphi,itheta)*
     $     (
     $     egrid(1,iphi,itheta)*conjg(unm3(1,n,m))+
     $     egrid(2,iphi,itheta)*conjg(unm3(2,n,m))+
     $     egrid(3,iphi,itheta)*conjg(unm3(3,n,m))
     $     )
        enddo
        enddo
c
        do n=1,nterms
        do m=-n,n
        ampole(n,m)=ampole(n,m)+weights(iphi,itheta)*
     $     (
     $     hgrid(1,iphi,itheta)*conjg(unm3(1,n,m))+
     $     hgrid(2,iphi,itheta)*conjg(unm3(2,n,m))+
     $     hgrid(3,iphi,itheta)*conjg(unm3(3,n,m))
     $     )
        bmpole(n,m)=bmpole(n,m)+weights(iphi,itheta)*
     $     (
     $     hgrid(1,iphi,itheta)*conjg(xnm3(1,n,m))+
     $     hgrid(2,iphi,itheta)*conjg(xnm3(2,n,m))+
     $     hgrid(3,iphi,itheta)*conjg(xnm3(3,n,m))
     $     )
        enddo
        enddo
c
 1300   continue
 1400   continue
c
c
c       ... now, rescale the multipole expansions
c
        z=rk*radius
c
        if( itype .eq. 1 ) then
c
c       ... evaluate hvals, emthvals, emrhvals
c
        call emhevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c
        if( itype .eq. 2 ) then
c
c       ... evaluate jvals, emtjvals, emrjvals
c
        call emjevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c

        do n=1,nterms
        do m=-n,n
c
        cd=qvals(n)
        cp=emtqvals(n)
        cr=emrqvals(n)
        ampole(n,m)=ampole(n,m)/(cd+cp*ima)/(4*pi)
        bmpole(n,m)=bmpole(n,m)/(cd-cp*ima)/(4*pi)
c
        enddo
        enddo
c
c
ccc        call prin2('ampole=*',ampole,(nterms+1)*(2*nterms+1)*2)
ccc        call prin2('bmpole=*',bmpole,(nterms+1)*(2*nterms+1)*2)
c

        return
        end
c
c
c
c
c
        subroutine em3mpevaleh(rk,center,ampole,bmpole,nterms,
     $     gridc,radius,egrid,hgrid,rnodes,weights,nphi,ntheta)
        implicit real *8 (a-h,o-z)
c        
        complex *16 rk,z
        dimension center(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        dimension gridc(3)
c
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c       
        itype=1
        call em3exevaleh_slow(itype,rk,center,ampole,bmpole,nterms,
     $     gridc,radius,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        return
        end
c
c
c
c
c
        subroutine em3taevaleh(rk,center,ampole,bmpole,nterms,
     $     gridc,radius,egrid,hgrid,rnodes,weights,nphi,ntheta)
        implicit real *8 (a-h,o-z)
c        
        complex *16 rk,z
        dimension center(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        dimension gridc(3)
c
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c       
        itype=2
        call em3exevaleh_slow(itype,rk,center,ampole,bmpole,nterms,
     $     gridc,radius,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        return
        end
c
c
c
c
c
        subroutine em3exevaleh_slow
     $     (itype,rk,center,ampole,bmpole,nterms,
     $     gridc,radius,egrid,hgrid,rnodes,weights,nphi,ntheta)
        implicit real *8 (a-h,o-z)
c        
        complex *16 rk
        dimension center(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        dimension gridc(3)
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c
        dimension target(3)
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
c
        nnodes=nphi*ntheta
c
        do 1400 itheta=1,ntheta
        do 1300 iphi=1,nphi
c       
        target(1)=gridc(1)+radius*rnodes(1,iphi,itheta)
        target(2)=gridc(2)+radius*rnodes(2,iphi,itheta)
        target(3)=gridc(3)+radius*rnodes(3,iphi,itheta)
c
        call em3exeval(itype,rk,center,ampole,bmpole,nterms,target,
     $     egrid(1,iphi,itheta),hgrid(1,iphi,itheta))
c
 1300   continue
 1400   continue
c
        return
        end
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Diagonal far field translation operators
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine em3ehfar(rk,xyz,cjvec,cmvec,
     $     gridc,radius,egrid,hgrid,rnodes,weights,nphi,ntheta)
        implicit real *8 (a-h,o-z)
c        
        complex *16 rk,z
        dimension xyz(3)
        dimension gridc(3)
c
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c
        dimension dir(3)
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
c
        nnodes=nphi*ntheta
c
        do 1400 itheta=1,ntheta
        do 1300 iphi=1,nphi
c       
        dir(1)=rnodes(1,iphi,itheta)
        dir(2)=rnodes(2,iphi,itheta)
        dir(3)=rnodes(3,iphi,itheta)
c
        call dipole3emfar(rk,xyz,cjvec,cmvec,dir,
     $     egrid(1,iphi,itheta),hgrid(1,iphi,itheta))
c
 1300   continue
 1400   continue
c
        return
        end
c
c
c
c
c
        subroutine em3farfar
     $     (rk,center0,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,nterms)
        implicit real *8 (a-h,o-z)
c       
        dimension center0(3)
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c      
        dimension center1(3)
c
        complex *16 z,trans,ima
        data ima/(0.0d0,1.0d0)/
c
c
        dx=center1(1)-center0(1)
        dy=center1(2)-center0(2)
        dz=center1(3)-center0(3)
c
c
        do 1400 itheta=1,ntheta
        do 1300 iphi=1,nphi
c
        proj=
     $     dx*rnodes(1,iphi,itheta)+
     $     dy*rnodes(2,iphi,itheta)+
     $     dz*rnodes(3,iphi,itheta)
c
        trans=exp(ima*rk*proj)
c
        egrid(1,iphi,itheta)=egrid(1,iphi,itheta)*trans
        egrid(2,iphi,itheta)=egrid(2,iphi,itheta)*trans
        egrid(3,iphi,itheta)=egrid(3,iphi,itheta)*trans
c
        hgrid(1,iphi,itheta)=hgrid(1,iphi,itheta)*trans
        hgrid(2,iphi,itheta)=hgrid(2,iphi,itheta)*trans
        hgrid(3,iphi,itheta)=hgrid(3,iphi,itheta)*trans
c
 1300   continue
 1400   continue
c
        return
        end
c
c
c
c
c
        subroutine em3pllpll
     $     (rk,center0,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,nterms)
        implicit real *8 (a-h,o-z)
c       
        dimension center0(3)
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c      
        dimension center1(3)
c
        complex *16 z,trans,ima
        data ima/(0.0d0,1.0d0)/
c
c
        dx=center1(1)-center0(1)
        dy=center1(2)-center0(2)
        dz=center1(3)-center0(3)
c
        do 1400 itheta=1,ntheta
        do 1300 iphi=1,nphi
c
        proj=
     $     dx*rnodes(1,iphi,itheta)+
     $     dy*rnodes(2,iphi,itheta)+
     $     dz*rnodes(3,iphi,itheta)
c
        trans=exp(ima*rk*proj)
c
        egrid(1,iphi,itheta)=egrid(1,iphi,itheta)*trans
        egrid(2,iphi,itheta)=egrid(2,iphi,itheta)*trans
        egrid(3,iphi,itheta)=egrid(3,iphi,itheta)*trans
c
        hgrid(1,iphi,itheta)=hgrid(1,iphi,itheta)*trans
        hgrid(2,iphi,itheta)=hgrid(2,iphi,itheta)*trans
        hgrid(3,iphi,itheta)=hgrid(3,iphi,itheta)*trans
c
 1300   continue
 1400   continue
c
        return
        end
c
c
c
c
c
        subroutine em3farpll
     $     (rk,center0,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,nterms)
        implicit real *8 (a-h,o-z)
c       
        dimension center0(3)
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c      
        dimension center1(3)
c
        real *8, allocatable :: pols(:)
c
        complex *16, allocatable :: hvals(:)
        complex *16, allocatable :: hders(:)
        complex *16, allocatable :: emthvals(:)
        complex *16, allocatable :: emrhvals(:)
c
        complex *16 z,trans,ima
        data ima/(0.0d0,1.0d0)/
c
        allocate( pols(0:nterms) )
        allocate( hvals(0:nterms) )
        allocate( hders(0:nterms) )
        allocate( emthvals(0:nterms) )
        allocate( emrhvals(0:nterms) )
c
        dx=center1(1)-center0(1)
        dy=center1(2)-center0(2)
        dz=center1(3)-center0(3)
c
        r=sqrt(dx*dx+dy*dy+dz*dz)
c
        rx=dx/r
        ry=dy/r
        rz=dz/r
c
        z=rk*r
        call emhevalrt(nterms,z,hvals,hders,emthvals,emrhvals)
c
c
        do 1400 itheta=1,ntheta
        do 1300 iphi=1,nphi
c
        cosphi=
     $     rx*rnodes(1,iphi,itheta)+
     $     ry*rnodes(2,iphi,itheta)+
     $     rz*rnodes(3,iphi,itheta)
c
        call em3legpols(cosphi,nterms,pols)
c
        trans=0
        do 1200 n=0,nterms
        trans=trans+(2*n+1)*hvals(n)*pols(n) *ima**n
 1200   continue
c
        egrid(1,iphi,itheta)=egrid(1,iphi,itheta)*trans
        egrid(2,iphi,itheta)=egrid(2,iphi,itheta)*trans
        egrid(3,iphi,itheta)=egrid(3,iphi,itheta)*trans
c
        hgrid(1,iphi,itheta)=hgrid(1,iphi,itheta)*trans
        hgrid(2,iphi,itheta)=hgrid(2,iphi,itheta)*trans
        hgrid(3,iphi,itheta)=hgrid(3,iphi,itheta)*trans
c
 1300   continue
 1400   continue
c
        return
        end
c
c
c
c
c
       subroutine em3legpols(x,n,p)
       implicit real *8 (a-h,o-z)
       dimension p(0:n)
c
c       this subroutine calculates legendre polynomials
c       functions of orders 0 through n of the real
c       argument x.
c
c                  input parameters:
c
c  x - the argument
c  n - the highest order of the legendre polynomial
c        to be returned
c
c                  output parameters:
c
c  p - the legendre polynomials of x.
c        note that this array is assumed to have element
c        number zero.
c
c        . . . construct p0(x), p1(x)
c
       p(0)=1
       p(1)=x
c
c       conduct the recursion
c
       do 1200 i=1,n-1
       p(i+1)=( (2*i+1)*x*p(i)-i*p(i-1) )/(i+1)
 1200 continue
       return
       end
c
c
c
c
c
        subroutine em3plleval
     $     (rk,xyz,center,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,evec,hvec)
        implicit real *8 (a-h,o-z)
c
        dimension xyz(3)
        complex *16 evec(3),hvec(3)
c
        dimension center(3)
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c      
        dimension center1(3)
c
        complex *16 trans,ima
        data ima/(0.0d0,1.0d0)/
c
c
        dx=xyz(1)-center(1)
        dy=xyz(2)-center(2)
        dz=xyz(3)-center(3)
c
        evec(1)=0
        evec(2)=0
        evec(3)=0
c
        hvec(1)=0
        hvec(2)=0
        hvec(3)=0
c
c
        do 1400 itheta=1,ntheta
        do 1300 iphi=1,nphi
c
        proj=
     $     dx*rnodes(1,iphi,itheta)+
     $     dy*rnodes(2,iphi,itheta)+
     $     dz*rnodes(3,iphi,itheta)
c
        trans=exp(ima*rk*proj)
c
        evec(1)=evec(1)+egrid(1,iphi,itheta)*trans*weights(iphi,itheta)
        evec(2)=evec(2)+egrid(2,iphi,itheta)*trans*weights(iphi,itheta)
        evec(3)=evec(3)+egrid(3,iphi,itheta)*trans*weights(iphi,itheta)
c
        hvec(1)=hvec(1)+hgrid(1,iphi,itheta)*trans*weights(iphi,itheta)
        hvec(2)=hvec(2)+hgrid(2,iphi,itheta)*trans*weights(iphi,itheta)
        hvec(3)=hvec(3)+hgrid(3,iphi,itheta)*trans*weights(iphi,itheta)
c
 1300   continue
 1400   continue
c
c
        done=1
        pi=4*atan(done)
c
        evec(1)=evec(1)*ima/(4*pi)
        evec(2)=evec(2)*ima/(4*pi)
        evec(3)=evec(3)*ima/(4*pi)
c
        hvec(1)=hvec(1)*ima/(4*pi)
        hvec(2)=hvec(2)*ima/(4*pi)
        hvec(3)=hvec(3)*ima/(4*pi)
c
        return
        end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Translation operators
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine em3mpmp_slow(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: multipole to multipole
c
c
c       Shift the electromagnetic multipole (ampole,bmpole) 
c       located at the location center to the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        call em3mpevaleh(rk,center,ampole,bmpole,nterms,
     $     center1,radius1,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        call em3ehformmp
     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,ampole1,bmpole1,nterms1)
c       
        return
        end
c
c
c
        subroutine em3mpmp_fast(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: multipole to multipole
c
c
c       Shift the electromagnetic multipole (ampole,bmpole) 
c       located at the location center to the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        call em3mpmp(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta)
        return
        end
c
c
c
        subroutine em3mpmp(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: multipole to multipole
c
c
c       Shift the electromagnetic multipole (ampole,bmpole) 
c       located at the location center to the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        dimension dir(3)
        dimension center(3),center1(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 ampole1(0:nterms1,-nterms1:nterms1)
        complex *16 bmpole1(0:nterms1,-nterms1:nterms1)
c
        complex *16, allocatable :: ampole2(:,:)
        complex *16, allocatable :: bmpole2(:,:)
c
        complex *16, allocatable :: ampole3(:,:)
        complex *16, allocatable :: bmpole3(:,:)
c
        allocate( ampole2(0:nterms,-nterms:nterms) )
        allocate( bmpole2(0:nterms,-nterms:nterms) )
        allocate( ampole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( bmpole3(0:nterms1,-nterms1:nterms1) ) 
c
c       ... translation operator: multipole to multipole
c
        dir(1)=center1(1)-center(1)
        dir(2)=center1(2)-center(2)
        dir(3)=center1(3)-center(3)
c
        zshift=sqrt(dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3))
c
        call emabrotyzf(nterms,dir,ampole,ampole2)
        call emabrotyzf(nterms,dir,bmpole,bmpole2)
c
        itype=1
        call em3exevalehz_fast(itype,rk,ampole2,bmpole2,nterms,
     $     zshift,radius1,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        call em3ehformmp
     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,ampole1,bmpole1,nterms1)
c       
        call emabrotyzb(nterms1,dir,ampole1,ampole3)
        call emabrotyzb(nterms1,dir,bmpole1,bmpole3)
c
        do n=0,nterms1
        do m=-nterms1,nterms1
        ampole1(n,m)=0
        bmpole1(n,m)=0
        enddo
        enddo
c
        do n=0,nterms1
        do m=-n,n
        ampole1(n,m)=ampole3(n,m)
        bmpole1(n,m)=bmpole3(n,m)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine em3mpmp2(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: multipole to multipole
c
c
c       Shift the electromagnetic multipole (ampole,bmpole) 
c       located at the location center to the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        dimension dir(3)
        dimension center(3),center1(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 ampole1(0:nterms1,-nterms1:nterms1)
        complex *16 bmpole1(0:nterms1,-nterms1:nterms1)
c
        complex *16, allocatable :: ampole2(:,:)
        complex *16, allocatable :: bmpole2(:,:)
c
        complex *16, allocatable :: ampole3(:,:)
        complex *16, allocatable :: bmpole3(:,:)
c
        allocate( ampole2(0:nterms,-nterms:nterms) )
        allocate( bmpole2(0:nterms,-nterms:nterms) )
        allocate( ampole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( bmpole3(0:nterms1,-nterms1:nterms1) ) 
c
c       ... translation operator: multipole to multipole
c
        dir(1)=center1(1)-center(1)
        dir(2)=center1(2)-center(2)
        dir(3)=center1(3)-center(3)
c
        zshift=sqrt(dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3))
c
        call emabrotyzf(nterms,dir,ampole,ampole2)
        call emabrotyzf(nterms,dir,bmpole,bmpole2)
c
        itype=1
        call em3exevalehz_fast2(itype,rk,ampole2,bmpole2,nterms,
     $     zshift,radius1,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        call em3ehformmp2
     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,ampole1,bmpole1,nterms1)
c       
        call emabrotyzb(nterms1,dir,ampole1,ampole3)
        call emabrotyzb(nterms1,dir,bmpole1,bmpole3)
c
        do n=0,nterms1
        do m=-nterms1,nterms1
        ampole1(n,m)=0
        bmpole1(n,m)=0
        enddo
        enddo
c
        do n=0,nterms1
        do m=-n,n
        ampole1(n,m)=ampole3(n,m)
        bmpole1(n,m)=bmpole3(n,m)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine em3mpmp3(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,rnodes,weights,
     $     nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: multipole to multipole
c
c
c       Shift the electromagnetic multipole (ampole,bmpole) 
c       located at the location center to the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        dimension dir(3)
        dimension center(3),center1(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 ampole1(0:nterms1,-nterms1:nterms1)
        complex *16 bmpole1(0:nterms1,-nterms1:nterms1)
c
        complex *16, allocatable :: ampole2(:,:)
        complex *16, allocatable :: bmpole2(:,:)
c
        complex *16, allocatable :: ampole3(:,:)
        complex *16, allocatable :: bmpole3(:,:)
c
        complex *16, allocatable :: egrid(:,:,:)
        complex *16, allocatable :: hgrid(:,:,:)
c
        allocate( ampole2(0:nterms,-nterms:nterms) )
        allocate( bmpole2(0:nterms,-nterms:nterms) )
        allocate( ampole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( bmpole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( egrid(3,nphi,ntheta) )
        allocate( hgrid(3,nphi,ntheta) )
c
c
c       ... translation operator: multipole to multipole
c
        dir(1)=center1(1)-center(1)
        dir(2)=center1(2)-center(2)
        dir(3)=center1(3)-center(3)
c
        zshift=sqrt(dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3))
c
        call emabrotyzf(nterms,dir,ampole,ampole2)
        call emabrotyzf(nterms,dir,bmpole,bmpole2)
c
        itype=1
        call em3exevalehz_fast2(itype,rk,ampole2,bmpole2,nterms,
     $     zshift,radius1,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        call em3ehformmp2
     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,ampole1,bmpole1,nterms1)
c       
        call emabrotyzb(nterms1,dir,ampole1,ampole3)
        call emabrotyzb(nterms1,dir,bmpole1,bmpole3)
c
        do n=0,nterms1
        do m=-nterms1,nterms1
        ampole1(n,m)=0
        bmpole1(n,m)=0
        enddo
        enddo
c
        do n=0,nterms1
        do m=-n,n
        ampole1(n,m)=ampole3(n,m)
        bmpole1(n,m)=bmpole3(n,m)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine em3mpmp3trunc(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,rnodes,weights,
     $     nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: multipole to multipole
c
c
c       Shift the electromagnetic multipole (ampole,bmpole) 
c       located at the location center to the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        dimension dir(3)
        dimension center(3),center1(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 ampole1(0:nterms1,-nterms1:nterms1)
        complex *16 bmpole1(0:nterms1,-nterms1:nterms1)
c
        complex *16, allocatable :: ampole2(:,:)
        complex *16, allocatable :: bmpole2(:,:)
c
        complex *16, allocatable :: ampole3(:,:)
        complex *16, allocatable :: bmpole3(:,:)
c
        complex *16, allocatable :: egrid(:,:,:)
        complex *16, allocatable :: hgrid(:,:,:)
c
        allocate( ampole2(0:nterms,-nterms:nterms) )
        allocate( bmpole2(0:nterms,-nterms:nterms) )
        allocate( ampole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( bmpole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( egrid(3,nphi,ntheta) )
        allocate( hgrid(3,nphi,ntheta) )
c
c
c       ... translation operator: multipole to multipole
c
        dir(1)=center1(1)-center(1)
        dir(2)=center1(2)-center(2)
        dir(3)=center1(3)-center(3)
c
        zshift=sqrt(dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3))
c
        call emabrotyzf(nterms,dir,ampole,ampole2)
        call emabrotyzf(nterms,dir,bmpole,bmpole2)
c
        itype=1
        call em3exevalehz_fast2(itype,rk,ampole2,bmpole2,nterms,
     $     zshift,radius1,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
c        call em3ehformmp2
c     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
c     $     nphi,ntheta,center1,ampole1,bmpole1,nterms1)
c       
c        call emabrotyzb(nterms1,dir,ampole1,ampole3)
c        call emabrotyzb(nterms1,dir,bmpole1,bmpole3)
c
        call em3ehformmp2_trunc
     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,ampole1,bmpole1,nterms1,nterms)
c       
        call emabrotyzb_trunc(nterms1,nterms,dir,ampole1,ampole3)
        call emabrotyzb_trunc(nterms1,nterms,dir,bmpole1,bmpole3)
c
        do n=0,nterms1
        do m=-nterms1,nterms1
        ampole1(n,m)=0
        bmpole1(n,m)=0
        enddo
        enddo
c
        do n=0,nterms1
        do m=-n,n
        ampole1(n,m)=ampole3(n,m)
        bmpole1(n,m)=bmpole3(n,m)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine em3mpmp3_add(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,rnodes,weights,nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: multipole to multipole
c
c
c       Shift the electromagnetic multipole (ampole,bmpole) 
c       located at the location center to the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        dimension dir(3)
        dimension center(3),center1(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 ampole1(0:nterms1,-nterms1:nterms1)
        complex *16 bmpole1(0:nterms1,-nterms1:nterms1)
c
        complex *16, allocatable :: ampole2(:,:)
        complex *16, allocatable :: bmpole2(:,:)
c
        complex *16, allocatable :: ampole3(:,:)
        complex *16, allocatable :: bmpole3(:,:)
c
        complex *16, allocatable :: ampole4(:,:)
        complex *16, allocatable :: bmpole4(:,:)
c
        complex *16, allocatable :: egrid(:,:,:)
        complex *16, allocatable :: hgrid(:,:,:)
c
        allocate( ampole2(0:nterms,-nterms:nterms) )
        allocate( bmpole2(0:nterms,-nterms:nterms) )
        allocate( ampole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( bmpole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( ampole4(0:nterms1,-nterms1:nterms1) ) 
        allocate( bmpole4(0:nterms1,-nterms1:nterms1) ) 
        allocate( egrid(3,nphi,ntheta) )
        allocate( hgrid(3,nphi,ntheta) )
c
c
c       ... translation operator: multipole to multipole
c
        dir(1)=center1(1)-center(1)
        dir(2)=center1(2)-center(2)
        dir(3)=center1(3)-center(3)
c
        zshift=sqrt(dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3))
c
        call emabrotyzf(nterms,dir,ampole,ampole2)
        call emabrotyzf(nterms,dir,bmpole,bmpole2)
c
        itype=1
        call em3exevalehz_fast2(itype,rk,ampole2,bmpole2,nterms,
     $     zshift,radius1,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        call em3ehformmp2
     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,ampole4,bmpole4,nterms1)
c       
        call emabrotyzb(nterms1,dir,ampole4,ampole3)
        call emabrotyzb(nterms1,dir,bmpole4,bmpole3)
c
        do n=0,nterms1
        do m=-n,n
        ampole1(n,m)=ampole1(n,m)+ampole3(n,m)
        bmpole1(n,m)=bmpole1(n,m)+bmpole3(n,m)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine em3mpta_slow(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: multipole to local
c
c
c       Convert the electromagnetic multipole (ampole,bmpole) 
c       located at the location center to the electromagnetic 
c       local expansion multipole at the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        call em3mpevaleh(rk,center,ampole,bmpole,nterms,
     $     center1,radius1,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        call em3ehformta
     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,ampole1,bmpole1,nterms1)
c       
        return
        end
c
c
c
        subroutine em3mpta_fast(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta)
        implicit real *8 (a-h,o-z)
        call em3mpta(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta)
        return
        end
c
c
c
        subroutine em3mpta(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: multipole to local
c
c
c       Convert the electromagnetic multipole (ampole,bmpole) 
c       located at the location center to the electromagnetic 
c       local expansion multipole at the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        dimension dir(3)
        dimension center(3),center1(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 ampole1(0:nterms1,-nterms1:nterms1)
        complex *16 bmpole1(0:nterms1,-nterms1:nterms1)
c
        complex *16, allocatable :: ampole2(:,:)
        complex *16, allocatable :: bmpole2(:,:)
c
        complex *16, allocatable :: ampole3(:,:)
        complex *16, allocatable :: bmpole3(:,:)
c
        allocate( ampole2(0:nterms,-nterms:nterms) )
        allocate( bmpole2(0:nterms,-nterms:nterms) )
        allocate( ampole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( bmpole3(0:nterms1,-nterms1:nterms1) ) 
c
c       ... translation operator: multipole to local
c
        dir(1)=center1(1)-center(1)
        dir(2)=center1(2)-center(2)
        dir(3)=center1(3)-center(3)
c
        zshift=sqrt(dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3))
c
        call emabrotyzf(nterms,dir,ampole,ampole2)
        call emabrotyzf(nterms,dir,bmpole,bmpole2)
c
        itype=1
        call em3exevalehz_fast(itype,rk,ampole2,bmpole2,nterms,
     $     zshift,radius1,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        call em3ehformta
     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,ampole1,bmpole1,nterms1)
c       
        call emabrotyzb(nterms1,dir,ampole1,ampole3)
        call emabrotyzb(nterms1,dir,bmpole1,bmpole3)
c
        do n=0,nterms1
        do m=-nterms1,nterms1
        ampole1(n,m)=0
        bmpole1(n,m)=0
        enddo
        enddo
c
        do n=0,nterms1
        do m=-n,n
        ampole1(n,m)=ampole3(n,m)
        bmpole1(n,m)=bmpole3(n,m)
        enddo
        enddo
c
        return
        end
c
c
c
        subroutine em3mpta2(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: multipole to local
c
c
c       Convert the electromagnetic multipole (ampole,bmpole) 
c       located at the location center to the electromagnetic 
c       local expansion multipole at the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        dimension dir(3)
        dimension center(3),center1(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 ampole1(0:nterms1,-nterms1:nterms1)
        complex *16 bmpole1(0:nterms1,-nterms1:nterms1)
c
        complex *16, allocatable :: ampole2(:,:)
        complex *16, allocatable :: bmpole2(:,:)
c
        complex *16, allocatable :: ampole3(:,:)
        complex *16, allocatable :: bmpole3(:,:)
c
        allocate( ampole2(0:nterms,-nterms:nterms) )
        allocate( bmpole2(0:nterms,-nterms:nterms) )
        allocate( ampole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( bmpole3(0:nterms1,-nterms1:nterms1) ) 
c
c       ... translation operator: multipole to multipole
c
        dir(1)=center1(1)-center(1)
        dir(2)=center1(2)-center(2)
        dir(3)=center1(3)-center(3)
c
        zshift=sqrt(dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3))
c
        call emabrotyzf(nterms,dir,ampole,ampole2)
        call emabrotyzf(nterms,dir,bmpole,bmpole2)
c
        itype=1
        call em3exevalehz_fast2(itype,rk,ampole2,bmpole2,nterms,
     $     zshift,radius1,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        call em3ehformta2
     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,ampole1,bmpole1,nterms1)
c       
        call emabrotyzb(nterms1,dir,ampole1,ampole3)
        call emabrotyzb(nterms1,dir,bmpole1,bmpole3)
c
        do n=0,nterms1
        do m=-nterms1,nterms1
        ampole1(n,m)=0
        bmpole1(n,m)=0
        enddo
        enddo
c
        do n=0,nterms1
        do m=-n,n
        ampole1(n,m)=ampole3(n,m)
        bmpole1(n,m)=bmpole3(n,m)
        enddo
        enddo
c
        return
        end
c
c
c
        subroutine em3mpta3(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,rnodes,weights,nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: multipole to local
c
c
c       Convert the electromagnetic multipole (ampole,bmpole) 
c       located at the location center to the electromagnetic 
c       local expansion multipole at the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        dimension dir(3)
        dimension center(3),center1(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 ampole1(0:nterms1,-nterms1:nterms1)
        complex *16 bmpole1(0:nterms1,-nterms1:nterms1)
c
        complex *16, allocatable :: ampole2(:,:)
        complex *16, allocatable :: bmpole2(:,:)
c
        complex *16, allocatable :: ampole3(:,:)
        complex *16, allocatable :: bmpole3(:,:)
c
        complex *16, allocatable :: egrid(:,:,:)
        complex *16, allocatable :: hgrid(:,:,:)
c
        allocate( ampole2(0:nterms,-nterms:nterms) )
        allocate( bmpole2(0:nterms,-nterms:nterms) )
        allocate( ampole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( bmpole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( egrid(3,nphi,ntheta) )
        allocate( hgrid(3,nphi,ntheta) )
c
c       ... translation operator: multipole to local
c
        dir(1)=center1(1)-center(1)
        dir(2)=center1(2)-center(2)
        dir(3)=center1(3)-center(3)
c
        zshift=sqrt(dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3))
c
        call emabrotyzf(nterms,dir,ampole,ampole2)
        call emabrotyzf(nterms,dir,bmpole,bmpole2)
c
        itype=1
        call em3exevalehz_fast2(itype,rk,ampole2,bmpole2,nterms,
     $     zshift,radius1,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        call em3ehformta2
     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,ampole1,bmpole1,nterms1)
c       
        call emabrotyzb(nterms1,dir,ampole1,ampole3)
        call emabrotyzb(nterms1,dir,bmpole1,bmpole3)
c
        do n=0,nterms1
        do m=-nterms1,nterms1
        ampole1(n,m)=0
        bmpole1(n,m)=0
        enddo
        enddo
c
        do n=0,nterms1
        do m=-n,n
        ampole1(n,m)=ampole3(n,m)
        bmpole1(n,m)=bmpole3(n,m)
        enddo
        enddo
c
        return
        end
c
c
c
        subroutine em3mpta3_add(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,rnodes,weights,nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: multipole to local
c
c
c       Convert the electromagnetic multipole (ampole,bmpole) 
c       located at the location center to the electromagnetic 
c       local expansion multipole at the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        dimension dir(3)
        dimension center(3),center1(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 ampole1(0:nterms1,-nterms1:nterms1)
        complex *16 bmpole1(0:nterms1,-nterms1:nterms1)
c
        complex *16, allocatable :: ampole2(:,:)
        complex *16, allocatable :: bmpole2(:,:)
c
        complex *16, allocatable :: ampole3(:,:)
        complex *16, allocatable :: bmpole3(:,:)
c
        complex *16, allocatable :: ampole4(:,:)
        complex *16, allocatable :: bmpole4(:,:)
c
        complex *16, allocatable :: egrid(:,:,:)
        complex *16, allocatable :: hgrid(:,:,:)
c
        allocate( ampole2(0:nterms,-nterms:nterms) )
        allocate( bmpole2(0:nterms,-nterms:nterms) )
        allocate( ampole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( bmpole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( ampole4(0:nterms1,-nterms1:nterms1) ) 
        allocate( bmpole4(0:nterms1,-nterms1:nterms1) ) 
        allocate( egrid(3,nphi,ntheta) )
        allocate( hgrid(3,nphi,ntheta) )
c
c       ... translation operator: multipole to local
c
        dir(1)=center1(1)-center(1)
        dir(2)=center1(2)-center(2)
        dir(3)=center1(3)-center(3)
c
        zshift=sqrt(dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3))
c
        call emabrotyzf(nterms,dir,ampole,ampole2)
        call emabrotyzf(nterms,dir,bmpole,bmpole2)
c
        itype=1
        call em3exevalehz_fast2(itype,rk,ampole2,bmpole2,nterms,
     $     zshift,radius1,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        call em3ehformta2
     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,ampole4,bmpole4,nterms1)
c       
        call emabrotyzb(nterms1,dir,ampole4,ampole3)
        call emabrotyzb(nterms1,dir,bmpole4,bmpole3)
c
        do n=0,nterms1
        do m=-n,n
        ampole1(n,m)=ampole1(n,m)+ampole3(n,m)
        bmpole1(n,m)=bmpole1(n,m)+bmpole3(n,m)
        enddo
        enddo
c
        return
        end
c
c
c
        subroutine em3tata_slow(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: local to local
c
c
c       Shift the electromagnetic local expansion (ampole,bmpole) 
c       located at the location center to the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        call em3taevaleh(rk,center,ampole,bmpole,nterms,
     $     center1,radius1,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        call em3ehformta
     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,ampole1,bmpole1,nterms1)
c       
        return
        end
c
c
c
        subroutine em3tata_fast(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: local to local
c
c
c       Shift the electromagnetic local expansion (ampole,bmpole) 
c       located at the location center to the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        call em3tata(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta)
        return
        end
c
c
c
        subroutine em3tata(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: local to local
c
c
c       Shift the electromagnetic local expansion (ampole,bmpole) 
c       located at the location center to the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        dimension dir(3)
        dimension center(3),center1(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 ampole1(0:nterms1,-nterms1:nterms1)
        complex *16 bmpole1(0:nterms1,-nterms1:nterms1)
c
        complex *16, allocatable :: ampole2(:,:)
        complex *16, allocatable :: bmpole2(:,:)
c
        complex *16, allocatable :: ampole3(:,:)
        complex *16, allocatable :: bmpole3(:,:)
c
        allocate( ampole2(0:nterms,-nterms:nterms) )
        allocate( bmpole2(0:nterms,-nterms:nterms) )
        allocate( ampole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( bmpole3(0:nterms1,-nterms1:nterms1) ) 
c
c       ... translation operator: local to local
c
        dir(1)=center1(1)-center(1)
        dir(2)=center1(2)-center(2)
        dir(3)=center1(3)-center(3)
c
        zshift=sqrt(dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3))
c
        call emabrotyzf(nterms,dir,ampole,ampole2)
        call emabrotyzf(nterms,dir,bmpole,bmpole2)
c
        itype=2
        call em3exevalehz_fast(itype,rk,ampole2,bmpole2,nterms,
     $     zshift,radius1,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        call em3ehformta
     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,ampole1,bmpole1,nterms1)
c       
        call emabrotyzb(nterms1,dir,ampole1,ampole3)
        call emabrotyzb(nterms1,dir,bmpole1,bmpole3)
c
        do n=0,nterms1
        do m=-nterms1,nterms1
        ampole1(n,m)=0
        bmpole1(n,m)=0
        enddo
        enddo
c
        do n=0,nterms1
        do m=-n,n
        ampole1(n,m)=ampole3(n,m)
        bmpole1(n,m)=bmpole3(n,m)
        enddo
        enddo
c
        return
        end
c
c
c
        subroutine em3tata2(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: local to local
c
c
c       Shift the electromagnetic local expansion (ampole,bmpole) 
c       located at the location center to the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        dimension dir(3)
        dimension center(3),center1(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 ampole1(0:nterms1,-nterms1:nterms1)
        complex *16 bmpole1(0:nterms1,-nterms1:nterms1)
c
        complex *16, allocatable :: ampole2(:,:)
        complex *16, allocatable :: bmpole2(:,:)
c
        complex *16, allocatable :: ampole3(:,:)
        complex *16, allocatable :: bmpole3(:,:)
c
        allocate( ampole2(0:nterms,-nterms:nterms) )
        allocate( bmpole2(0:nterms,-nterms:nterms) )
        allocate( ampole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( bmpole3(0:nterms1,-nterms1:nterms1) ) 
c
c       ... translation operator: local to local
c
        dir(1)=center1(1)-center(1)
        dir(2)=center1(2)-center(2)
        dir(3)=center1(3)-center(3)
c
        zshift=sqrt(dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3))
c
        call emabrotyzf(nterms,dir,ampole,ampole2)
        call emabrotyzf(nterms,dir,bmpole,bmpole2)
c
        itype=2
        call em3exevalehz_fast2(itype,rk,ampole2,bmpole2,nterms,
     $     zshift,radius1,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        call em3ehformta2
     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,ampole1,bmpole1,nterms1)
c       
        call emabrotyzb(nterms1,dir,ampole1,ampole3)
        call emabrotyzb(nterms1,dir,bmpole1,bmpole3)
c
        do n=0,nterms1
        do m=-nterms1,nterms1
        ampole1(n,m)=0
        bmpole1(n,m)=0
        enddo
        enddo
c
        do n=0,nterms1
        do m=-n,n
        ampole1(n,m)=ampole3(n,m)
        bmpole1(n,m)=bmpole3(n,m)
        enddo
        enddo
c
        return
        end
c
c
c
        subroutine em3tata3(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,rnodes,weights,
     $     nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: local to local
c
c
c       Shift the electromagnetic local expansion (ampole,bmpole) 
c       located at the location center to the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        dimension dir(3)
        dimension center(3),center1(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 ampole1(0:nterms1,-nterms1:nterms1)
        complex *16 bmpole1(0:nterms1,-nterms1:nterms1)
c
        complex *16, allocatable :: ampole2(:,:)
        complex *16, allocatable :: bmpole2(:,:)
c
        complex *16, allocatable :: ampole3(:,:)
        complex *16, allocatable :: bmpole3(:,:)
c
        complex *16, allocatable :: egrid(:,:,:)
        complex *16, allocatable :: hgrid(:,:,:)
c
        allocate( ampole2(0:nterms,-nterms:nterms) )
        allocate( bmpole2(0:nterms,-nterms:nterms) )
        allocate( ampole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( bmpole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( egrid(3,nphi,ntheta) )
        allocate( hgrid(3,nphi,ntheta) )
c
c       ... translation operator: local to local
c
        dir(1)=center1(1)-center(1)
        dir(2)=center1(2)-center(2)
        dir(3)=center1(3)-center(3)
c
        zshift=sqrt(dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3))
c
        call emabrotyzf(nterms,dir,ampole,ampole2)
        call emabrotyzf(nterms,dir,bmpole,bmpole2)
c
        itype=2
        call em3exevalehz_fast2(itype,rk,ampole2,bmpole2,nterms,
     $     zshift,radius1,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        call em3ehformta2
     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,ampole1,bmpole1,nterms1)
c       
        call emabrotyzb(nterms1,dir,ampole1,ampole3)
        call emabrotyzb(nterms1,dir,bmpole1,bmpole3)
c
        do n=0,nterms1
        do m=-nterms1,nterms1
        ampole1(n,m)=0
        bmpole1(n,m)=0
        enddo
        enddo
c
        do n=0,nterms1
        do m=-n,n
        ampole1(n,m)=ampole3(n,m)
        bmpole1(n,m)=bmpole3(n,m)
        enddo
        enddo
c
        return
        end
c
c
c
        subroutine em3tata3trunc(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,rnodes,weights,
     $     nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: local to local
c
c
c       Shift the electromagnetic local expansion (ampole,bmpole) 
c       located at the location center to the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        dimension dir(3)
        dimension center(3),center1(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 ampole1(0:nterms1,-nterms1:nterms1)
        complex *16 bmpole1(0:nterms1,-nterms1:nterms1)
c
        complex *16, allocatable :: ampole2(:,:)
        complex *16, allocatable :: bmpole2(:,:)
c
        complex *16, allocatable :: ampole3(:,:)
        complex *16, allocatable :: bmpole3(:,:)
c
        complex *16, allocatable :: egrid(:,:,:)
        complex *16, allocatable :: hgrid(:,:,:)
c
        allocate( ampole2(0:nterms,-nterms:nterms) )
        allocate( bmpole2(0:nterms,-nterms:nterms) )
        allocate( ampole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( bmpole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( egrid(3,nphi,ntheta) )
        allocate( hgrid(3,nphi,ntheta) )
c
c       ... translation operator: local to local
c
        dir(1)=center1(1)-center(1)
        dir(2)=center1(2)-center(2)
        dir(3)=center1(3)-center(3)
c
        zshift=sqrt(dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3))
c
ccc        call emabrotyzf(nterms,dir,ampole,ampole2)
ccc        call emabrotyzf(nterms,dir,bmpole,bmpole2)
c
        call emabrotyzf_trunc(nterms,nterms1,dir,ampole,ampole2)
        call emabrotyzf_trunc(nterms,nterms1,dir,bmpole,bmpole2)
c
        itype=2
        call em3exevalehz_fast2_trunc
     $     (itype,rk,ampole2,bmpole2,nterms,nterms1,
     $     zshift,radius1,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        call em3ehformta2
     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,ampole1,bmpole1,nterms1)
c       
        call emabrotyzb(nterms1,dir,ampole1,ampole3)
        call emabrotyzb(nterms1,dir,bmpole1,bmpole3)
c
        do n=0,nterms1
        do m=-nterms1,nterms1
        ampole1(n,m)=0
        bmpole1(n,m)=0
        enddo
        enddo
c
        do n=0,nterms1
        do m=-n,n
        ampole1(n,m)=ampole3(n,m)
        bmpole1(n,m)=bmpole3(n,m)
        enddo
        enddo
c
        return
        end
c
c
c
c
        subroutine em3tata3_add(rk,center,ampole,bmpole,nterms,
     $     center1,ampole1,bmpole1,nterms1,
     $     radius1,rnodes,weights,nphi,ntheta)
        implicit real *8 (a-h,o-z)
c
c       ... translation operator: local to local
c
c
c       Shift the electromagnetic local expansion (ampole,bmpole) 
c       located at the location center to the new location
c
c       This is a reasonably fast "point and shoot" version which first
c       rotates the coordinate system, then shifts along the Z-axis, and
c       then rotates back to the original coordinates.
c
        dimension dir(3)
        dimension center(3),center1(3)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 ampole1(0:nterms1,-nterms1:nterms1)
        complex *16 bmpole1(0:nterms1,-nterms1:nterms1)
c
        complex *16, allocatable :: ampole2(:,:)
        complex *16, allocatable :: bmpole2(:,:)
c
        complex *16, allocatable :: ampole3(:,:)
        complex *16, allocatable :: bmpole3(:,:)
c
        complex *16, allocatable :: ampole4(:,:)
        complex *16, allocatable :: bmpole4(:,:)
c
        complex *16, allocatable :: egrid(:,:,:)
        complex *16, allocatable :: hgrid(:,:,:)
c
        allocate( ampole2(0:nterms,-nterms:nterms) )
        allocate( bmpole2(0:nterms,-nterms:nterms) )
        allocate( ampole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( bmpole3(0:nterms1,-nterms1:nterms1) ) 
        allocate( ampole4(0:nterms1,-nterms1:nterms1) ) 
        allocate( bmpole4(0:nterms1,-nterms1:nterms1) ) 
        allocate( egrid(3,nphi,ntheta) )
        allocate( hgrid(3,nphi,ntheta) )
c
c       ... translation operator: local to local
c
        dir(1)=center1(1)-center(1)
        dir(2)=center1(2)-center(2)
        dir(3)=center1(3)-center(3)
c
        zshift=sqrt(dir(1)*dir(1)+dir(2)*dir(2)+dir(3)*dir(3))
c
        call emabrotyzf(nterms,dir,ampole,ampole2)
        call emabrotyzf(nterms,dir,bmpole,bmpole2)
c
        itype=2
        call em3exevalehz_fast2(itype,rk,ampole2,bmpole2,nterms,
     $     zshift,radius1,egrid,hgrid,rnodes,weights,nphi,ntheta)
c
        call em3ehformta2
     $     (rk,center1,radius1,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center1,ampole4,bmpole4,nterms1)
c       
        call emabrotyzb(nterms1,dir,ampole4,ampole3)
        call emabrotyzb(nterms1,dir,bmpole4,bmpole3)
c
        do n=0,nterms1
        do m=-n,n
        ampole1(n,m)=ampole1(n,m)+ampole3(n,m)
        bmpole1(n,m)=bmpole1(n,m)+bmpole3(n,m)
        enddo
        enddo
c
        return
        end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Spherical and cartesian coordinates
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine em3phi(x,y,phi)
        implicit real *8 (a-h,o-z)
c
        if( abs(x) .eq. 0 .and. abs(y) .eq. 0 ) then
        phi = 0
        else
        phi = atan2(y,x)
        endif
c
        return
        end
c
c
c
c
c
        subroutine em3ctheta(z,r,ctheta)
        implicit real *8 (a-h,o-z)
c
        if( abs(r) .gt. 0 ) then
        ctheta = z/r
        else
        ctheta = 0.0d0
        endif
c
        return
        end
c
c
c
c
c
        subroutine em3polar(x,y,z,r,ctheta,ephi)
        implicit real *8 (a-h,o-z)
        complex *16 ephi,ima
        data ima/(0.0d0,1.0d0)/
c
        r = sqrt(x*x+y*y+z*z)
c
        if( abs(r) .gt. 0 ) then
        ctheta = z/r
        else
        ctheta = 0.0d0
        endif
c
        proj = sqrt(x*x+y*y)
c
        if( abs(proj) .gt. 0 ) then
        ephi = x/proj + ima*y/proj
        else
        ephi = 0.0d0
        endif
c
        return
        end
c
c
c
c
c
        subroutine sph2cart(dx,dy,dz,phi,theta,nnout,unm2,unm3)
        implicit real *8 (a-h,o-z)
        dimension c(3,3)
        complex *16 unm2(2,1),unm3(3,1)
c
c       ... convert to cartesian coordinates
c
        r=sqrt(dx*dx+dy*dy+dz*dz)
c
        x=dx/r
        y=dy/r
        z=dz/r
c
        c(1,1)=-sin(phi)
        c(2,1)=cos(phi)
        c(3,1)=0
c
        c(1,2)=+cos(phi)*z
        c(2,2)=+sin(phi)*z
        c(3,2)=-sqrt(1-z**2) 
c
        c(1,3)=x
        c(2,3)=y
        c(3,3)=z
c
c        call prin2('x=*',x,1)
c        call prin2('y=*',y,1)
c        call prin2('z=*',z,1)
c        call prin2('phi=*',phi,1)
c        call prin2('c=*',c,3*3)
c
        do i=1,nnout
        unm3(1,i)=unm2(1,i)*c(1,1)+unm2(2,i)*c(1,2)
        unm3(2,i)=unm2(1,i)*c(2,1)+unm2(2,i)*c(2,2)
        unm3(3,i)=unm2(1,i)*c(3,1)+unm2(2,i)*c(3,2)
        enddo
c       
        return
        end
c
c
c
c
        subroutine cart2sphi(dx,dy,dz,phi,theta,cvals)
        implicit real *8 (a-h,o-z)
        dimension c(3,3)
        complex *16 cvals(3),c1,c2,c3
c
c       ... convert to spherical coordinates, in place
c
        r=sqrt(dx*dx+dy*dy+dz*dz)
c
        x=dx/r
        y=dy/r
        z=dz/r
c
        c(1,1)=-sin(phi)
        c(2,1)=cos(phi)
        c(3,1)=0
c       
        c(1,2)=+cos(phi)*z
        c(2,2)=+sin(phi)*z
        c(3,2)=-sqrt(1-z**2) 
c       
        c(1,3)=x
        c(2,3)=y
        c(3,3)=z
c
cccc        call prin2('c=*',c,3*3)
c
        c1=cvals(1)*c(1,1)+cvals(2)*c(2,1)+cvals(3)*c(3,1)
        c2=cvals(1)*c(1,2)+cvals(2)*c(2,2)+cvals(3)*c(3,2)
        c3=cvals(1)*c(1,3)+cvals(2)*c(2,3)+cvals(3)*c(3,3)
c
        cvals(1)=c1
        cvals(2)=c2
        cvals(3)=c3
c
        return
        end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Spherical Bessel j and h functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine emjheval(nterms,z,jvals,hvals)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the  j, and h values, for
c       the complex parameter z.
c
c       Input parameters:
c
c       nterms - the number of terms in em-multipole expansion
c
c       Output parameters:
c
c       hvals - the j_n(k r) values, (nterms+1) 
c       hvals - the h_n(k r) values, (nterms+1) 
c
c
        complex *16 z
c
        complex *16 jvals(0:nterms),hvals(0:nterms)
        complex *16 jders(0:nterms),hders(0:nterms)
c
        complex *16, allocatable :: work(:)
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        lwork=500000
        allocate( work(1:lwork) )
c
c       ... construct the spherical bessel j, h functions 
c
c       
        if( abs(z) < 1 ) then
        scale=abs(z)
        else
        scale=1
        endif
c       
        ifder=0
c
c
c       ... construct the spherical bessel h functions 
c
        call chfuns3d(ier,nterms,z,scale,hvals,ifder,hders,work,lwork)
ccc        call prin2('scale=*',scale,1)
c       
        do i=0,nterms
        hvals(i)=hvals(i)/scale**(i+1)
        if( ifder .eq. 1 ) hders(i)=hders(i)/scale**(i+1)
        enddo
c       
c
c       ... construct the spherical bessel j functions 
c
        call cjfuns3d(ier,nterms,z,scale,jvals,ifder,jders,work,lwork)
ccc        call prin2('scale=*',scale,1)
c       
        do i=0,nterms
        jvals(i)=jvals(i)*scale**(i)
        if( ifder .eq. 1 ) jders(i)=jders(i)*scale**(i)
        enddo
c       
        return
        end
c
c
c
c
c
        subroutine emhevalrt(nterms,z,hvals,hders,emthvals,emrhvals)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the  h values, derivatives
c       together with the tangential and radial em-scaling factors for
c       the complex parameter z.
c
c       Input parameters:
c
c       nterms - the number of terms in em-multipole expansion
c
c       Output parameters:
c
c       hvals - the h_n(k r) values (nterms+1) of them
c       hders - the h_n(k r) derivatives (nterms+1) of them
c       emthvals - the tangential outgoing em-scaling factor (nterms+1) of them
c       emrhvals - the radial outgoing em-scaling factor (nterms+1) of them
c
c
        complex *16 z
c
        complex *16 hvals(0:nterms),hders(0:nterms)
        complex *16 emthvals(0:nterms),emrhvals(0:nterms)
c
        complex *16, allocatable :: work(:)
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        lwork=500000
        allocate( work(1:lwork) )
c       
c
c       ... construct the spherical bessel h functions 
c
c       
        if( abs(z) < 1 ) then
        scale=abs(z)
        else
        scale=1
        endif
c       
        ifder=1
c
c
c       ... construct the spherical bessel h functions 
c
        call chfuns3d(ier,nterms,z,scale,hvals,ifder,hders,work,lwork)
ccc        call prin2('scale=*',scale,1)
c       
        do i=0,nterms
        hvals(i)=hvals(i)/scale**(i+1)
        hders(i)=hders(i)/scale**(i+1)
        enddo
c       
c
c       ... construct the scaling factors for the tangential and radial
c
c
        if( nterms .ge. 0 ) emthvals(0)=0
        if( nterms .ge. 0 ) emrhvals(0)=0
        do i=1,nterms
        emthvals(i)=(1/z*hvals(i)+hders(i)) 
        emrhvals(i)=1/z*hvals(i)*sqrt(i*(i+1.0d0))
        enddo
c
c
ccc        call prin2('bessel h vals=*',hvals,2*(nterms+1))
ccc        call prin2('bessel h ders=*',hders,2*(nterms+1))
c       
ccc        call prin2('emthvals=*',emthvals,2*(nterms+1))
ccc        call prin2('emrhvals=*',emrhvals,2*(nterms+1))
c
        return
        end
c
c
c
c
c
        subroutine emjevalrt(nterms,z,jvals,jders,emtjvals,emrjvals)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the  j values, derivatives
c       together with the tangential and radial em-scaling factors for
c       the complex parameter z.
c
c       Input parameters:
c
c       nterms - the number of terms in em-multipole expansion
c
c       Output parameters:
c
c       jvals - the j_n(k r) values (nterms+1) of them
c       jders - the j_n(k r) derivatives (nterms+1) of them
c       emtjvals - the tangential incoming em-scaling factor (nterms+1) of them
c       emrjvals - the radial incoming em-scaling factor (nterms+1) of them
c
c
        complex *16 z
c
        complex *16 jvals(0:nterms),jders(0:nterms)
        complex *16 emtjvals(0:nterms),emrjvals(0:nterms)
c
        complex *16, allocatable :: work(:)
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        lwork=500000
        allocate( work(1:lwork) )
c       
c
c       ... construct the spherical bessel j functions 
c
c
        if( abs(z) < 1 ) then
        scale=abs(z)
        else
        scale=1
        endif
c
c
        ifder=1
c
c
c       ... construct the spherical bessel j functions 
c
        call cjfuns3d(ier,nterms,z,scale,jvals,ifder,jders,work,lwork)
ccc        call prin2('scale=*',scale,1)
c       
        do i=0,nterms
        jvals(i)=jvals(i)*scale**(i)
        jders(i)=jders(i)*scale**(i)
        enddo
c       
c
c       ... construct the scaling factors for the tangential and radial
c
c
        if( abs(z) .gt. 0 ) then
c
        if( nterms .ge. 0 ) emtjvals(0)=0
        if( nterms .ge. 0 ) emrjvals(0)=0
        do i=1,nterms
        emtjvals(i)=(1/z*jvals(i)+jders(i)) 
        emrjvals(i)=1/z*jvals(i)*sqrt(i*(i+1.0d0))
        enddo
c
        else
c
        if( nterms .ge. 0 ) jvals(0)=1
        if( nterms .ge. 0 ) jders(0)=0
        if( nterms .ge. 0 ) emtjvals(0)=0
        if( nterms .ge. 0 ) emrjvals(0)=0
        if( nterms .ge. 1 ) jvals(1)=0
        if( nterms .ge. 1 ) jders(1)=1/3.0d0
        if( nterms .ge. 1 ) emtjvals(1)=2/3.0d0
        if( nterms .ge. 1 ) emrjvals(1)=1/3.0d0*sqrt(2.0d0)
        do i=2,nterms
        emtjvals(i)=0
        emrjvals(i)=0
        enddo
c
        endif
c
c
ccc        call prin2('bessel j vals=*',jvals,2*(nterms+1))
ccc        call prin2('bessel j ders=*',jders,2*(nterms+1))
c       
ccc        call prin2('emtjvals=*',emtjvals,2*(nterms+1))
ccc        call prin2('emrjvals=*',emrjvals,2*(nterms+1))
c
        return
        end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine em3zero(mpole,nterms)
        implicit real *8 (a-h,o-z)
        complex *16 mpole(0:nterms,-nterms:nterms)
c
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=0
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine em3prinm(mess,mpole,nterms)
        implicit real *8 (a-h,o-z)
        character *1 mess(1)
        complex *16 mpole(0:nterms,-nterms:nterms)
c
        call prin2(mess,i,0)
        do m=-nterms,nterms
        call prin2(' *',mpole(0,m),2*(nterms+1))
        enddo
c
        return
        end
c
c
c
c
c
        subroutine em3prinm1(mess,mpole,nterms,nterms1)
        implicit real *8 (a-h,o-z)
        character *1 mess(1)
        complex *16 mpole(0:nterms,-nterms:nterms)
c
        call prin2(mess,i,0)
        do m=-nterms1,nterms1
        call prin2(' *',mpole(0,m),2*(nterms1+1))
        enddo
c
        return
        end
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Maxwell's equations in R^3
c
c       Fast O(p^3) translation routines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       
c
c
c
        subroutine em3sphfft(cvals,nphi,ntheta,ts,ws)
        implicit real *8 (a-h,o-z)
c
        dimension ts(1),ws(1)
        complex *16 work(nphi)
        complex *16 cvals(3,nphi,ntheta)
        dimension c(3,3)
        complex *16 c1,c2,c3
c
        real *8, allocatable :: wsave(:)
c
        allocate( wsave(10*nphi+15) )
c
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
        call zffti(nphi,wsave)
c
        done=1
        pi=4*atan(done)
c
c
c       ... project and convert to spherical coordinates
c
        zshift=0
        zradius=1
        itype=1
        call em3sphcar(itype,zshift,zradius,cvals,nphi,ntheta,ts,ws)
c
c       ... perform the Fourier transform along each parallel
c
        do 1200 i=1,ntheta
c
        do j=1,nphi
        work(j)=cvals(1,j,i)
        enddo
        call zfftf(nphi,work,wsave)
        do j=1,nphi
        cvals(1,j,i)=work(j)
        enddo
c
        do j=1,nphi
        work(j)=cvals(2,j,i)
        enddo
        call zfftf(nphi,work,wsave)
        do j=1,nphi
        cvals(2,j,i)=work(j)
        enddo
c
        do j=1,nphi
        work(j)=cvals(3,j,i)
        enddo
        call zfftf(nphi,work,wsave)
        do j=1,nphi
        cvals(3,j,i)=work(j)
        enddo
c
 1200   continue
c
c
c       ... multiply by quadrature weights
c
        scale=(2*pi)/dble(nphi)
c
        do 1500 i=1,ntheta
        do 1400 j=1,nphi
        cvals(1,j,i)=cvals(1,j,i)*scale
        cvals(2,j,i)=cvals(2,j,i)*scale
        cvals(3,j,i)=cvals(3,j,i)*scale
 1400   continue
 1500   continue
c
c
        return
        end
c
c
c
c
c
        subroutine em3sph(cvals,nphi,ntheta,ts,ws)
        implicit real *8 (a-h,o-z)
c
        dimension ts(1),ws(1)
        complex *16 cvals(3,nphi,ntheta)
        dimension c(3,3)
        complex *16 c1,c2,c3
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
        done=1
        pi=4*atan(done)
c
c
c       ... project and convert to spherical coordinates
c
        zshift=0
        zradius=1
        itype=1
        call em3sphcar(itype,zshift,zradius,cvals,nphi,ntheta,ts,ws)
c
        return
        end
c
c
c
c
c
        subroutine em3fftcar(zshift,zradius,cvals,nphi,ntheta,ts,ws)
        implicit real *8 (a-h,o-z)
c
        dimension ts(1),ws(1)
        complex *16 cvals(3,nphi,ntheta)
        dimension c(3,3)
        complex *16 c1,c2,c3
c
        real *8, allocatable :: wsave(:)
        complex *16, allocatable :: work(:)
c
        allocate( wsave(10*nphi+15) )
        allocate( work(nphi) )
c
c
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
        call zffti(nphi,wsave)
c
        done=1
        pi=4*atan(done)
c
c
c       ... perform the Fourier transform along each parallel
c
        do 1200 i=1,ntheta
c
        do j=1,nphi
        work(j)=cvals(1,j,i)
        enddo
        call zfftb(nphi,work,wsave)
        do j=1,nphi
        cvals(1,j,i)=work(j)
        enddo
c
        do j=1,nphi
        work(j)=cvals(2,j,i)
        enddo
        call zfftb(nphi,work,wsave)
        do j=1,nphi
        cvals(2,j,i)=work(j)
        enddo
c
        do j=1,nphi
        work(j)=cvals(3,j,i)
        enddo
        call zfftb(nphi,work,wsave)
        do j=1,nphi
        cvals(3,j,i)=work(j)
        enddo
c
 1200   continue
c
c
c       ... multiply by quadrature weights
c
        scale=1
c
        do 1500 i=1,ntheta
        do 1400 j=1,nphi
        cvals(1,j,i)=cvals(1,j,i)*scale
        cvals(2,j,i)=cvals(2,j,i)*scale
        cvals(3,j,i)=cvals(3,j,i)*scale
 1400   continue
 1500   continue
c
c
c
c       ... convert the spherical coordinates to cartesian
c
        itype=2
        call em3sphcar(itype,zshift,zradius,cvals,nphi,ntheta,ts,ws)
c
        return
        end
c
c
c
c
c
        subroutine em3car(zshift,zradius,cvals,nphi,ntheta,ts,ws)
        implicit real *8 (a-h,o-z)
        dimension ts(1),ws(1)
        complex *16 cvals(3,nphi,ntheta)
        dimension c(3,3)
        complex *16 c1,c2,c3
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
        done=1
        pi=4*atan(done)
c
c
c       ... multiply by quadrature weights
c
        scale=2*pi
c
        do 1500 i=1,ntheta
        do 1400 j=1,nphi
        cvals(1,j,i)=cvals(1,j,i)*scale
        cvals(2,j,i)=cvals(2,j,i)*scale
        cvals(3,j,i)=cvals(3,j,i)*scale
 1400   continue
 1500   continue
c
c
c
c       ... convert the spherical coordinates to cartesian
c
        itype=2
        call em3sphcar(itype,zshift,zradius,cvals,nphi,ntheta,ts,ws)
c
        return
        end
c
c
c
c
c
        subroutine em3sphcar
     $     (itype,zshift,zradius,cvals,nphi,ntheta,ts,ws)
        implicit real *8 (a-h,o-z)
        dimension ts(1),ws(1)
        complex *16 cvals(3,nphi,ntheta)
        dimension c(3,3)
        complex *16 c1,c2,c3
c
        done=1
        pi=4*atan(done)
c
        do i=1,ntheta
c
        z=ts(i)
c
        z=zshift+ts(i)*zradius
        x=sqrt(1-ts(i)*ts(i))*zradius
        r=sqrt(x*x+z*z)
        z=z/r
c
        do j=1,nphi
c
           phi=(2*pi)*(j-1)/dble(nphi)
           x=cos(phi)*sqrt(1-z**2)
           y=sin(phi)*sqrt(1-z**2)
c
           if( itype .eq. 1 ) then
c
c       ... convert the cartesian coordinates to spherical
c
           c(1,1)=-sin(phi)
           c(2,1)=cos(phi)
           c(3,1)=0
c
           c(1,2)=+cos(phi)*z
           c(2,2)=+sin(phi)*z
           c(3,2)=-sqrt(1-z**2) 
c
           c(1,3)=x
           c(2,3)=y
           c(3,3)=z
c
           endif
c
           if( itype .eq. 2 ) then
c
c       ... convert the spherical coordinates to cartesian
c
c       ... this is a transpose (inverse) matrix
c
           c(1,1)=-sin(phi)
           c(1,2)=cos(phi)
           c(1,3)=0
c
           c(2,1)=+cos(phi)*z
           c(2,2)=+sin(phi)*z
           c(2,3)=-sqrt(1-z**2) 
c
           c(3,1)=x
           c(3,2)=y
           c(3,3)=z
c
           endif
c
           c1=
     $        cvals(1,j,i)*c(1,1)+
     $        cvals(2,j,i)*c(2,1)+
     $        cvals(3,j,i)*c(3,1)
c
           c2=
     $        cvals(1,j,i)*c(1,2)+
     $        cvals(2,j,i)*c(2,2)+
     $        cvals(3,j,i)*c(3,2)
c
           c3=
     $        cvals(1,j,i)*c(1,3)+
     $        cvals(2,j,i)*c(2,3)+
     $        cvals(3,j,i)*c(3,3)
c
           cvals(1,j,i)=c1
           cvals(2,j,i)=c2
           cvals(3,j,i)=c3
c
        enddo
        enddo
c
c
        return
        end
c
c
c
c
c
        subroutine em3ehformex_fast
     $     (itype,rk,gridc,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center,ampole,bmpole,nterms)
        implicit real *8 (a-h,o-z)
c
c       This subroutine forms the outgoing/incoming EM-multipole expansion
c       centered at the center due to the E and H field grids
c
c       Fast O(p^3) routine
c
c       Fast algorithm, perform FFT along all parallels
c
c          Input parameters:
c
c       
c       itype - the type of expansion
c          the expansions are outgoing, for itype=1
c          the expansions are incoming, for itype=2
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of dipoles in R^3
c       radius - the radius of E and H grids
c       egrid (complex *16)(3,nphi,ntheta) - E field grid
c       hgrid (complex *16)(3,nphi,ntheta) - H field grid
c       rnodes (real *8)(3,nphi,ntheta) - the grid node coordinates in R^3
c       weights (real *8)(nphi,ntheta) - the grid weights
c       nphi, ntheta - the number of grid points in each direction
c       center - the center of EM-multipole expansion
c
c          Output parameters:
c
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of EM-multipole expansion
c        
        complex *16 rk,z
        dimension xyz(3),center(3)
        dimension gridc(3),target(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c
        complex *16, allocatable :: qvals(:), qders(:)
        complex *16, allocatable :: emtqvals(:), emrqvals(:)
c
        real *8, allocatable :: pnm(:,:), dnm(:,:)
        complex *16, allocatable :: pxnm2(:,:,:)
        complex *16, allocatable :: punm2(:,:,:)
        complex *16, allocatable :: pxnm3(:,:,:)
        complex *16, allocatable :: punm3(:,:,:)
c
        complex *16, allocatable :: ynm(:,:)
        complex *16, allocatable :: xnm3(:,:,:)
        complex *16, allocatable :: unm3(:,:,:)
c
        complex *16, allocatable :: ephi(:)
c
        complex *16, allocatable :: xnm2(:,:,:)
        complex *16, allocatable :: unm2(:,:,:)
c
        real *8, allocatable :: ts(:), ws(:)
c
        complex *16, allocatable :: emgrid(:,:,:)
        complex *16, allocatable :: hmgrid(:,:,:)
c       
        complex *16 cd,cp,cr
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        allocate( qvals(0:nterms) )
        allocate( qders(0:nterms) )
        allocate( emtqvals(0:nterms) )
        allocate( emrqvals(0:nterms) )
        allocate( pnm(0:nterms,0:nterms) )
        allocate( dnm(0:nterms,0:nterms) )
        allocate( pxnm2(2,0:nterms,0:nterms) )
        allocate( punm2(2,0:nterms,0:nterms) ) 
        allocate( pxnm3(3,0:nterms,0:nterms) )
        allocate( punm3(3,0:nterms,0:nterms) ) 
        allocate( ynm(0:nterms,-nterms:nterms) )
        allocate( xnm3(3,0:nterms,-nterms:nterms) )
        allocate( unm3(3,0:nterms,-nterms:nterms) )
        allocate( ephi(-nterms:nterms) )
c
        allocate( xnm2(2,0:nterms,-nterms:nterms) )
        allocate( unm2(2,0:nterms,-nterms:nterms) )       
        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )
c
        allocate( emgrid(3,-nterms:nterms,ntheta) )
        allocate( hmgrid(3,-nterms:nterms,ntheta) )
c
c
        do n=0,nterms
        do m=-n,n
        ampole(n,m)=0
        bmpole(n,m)=0
        enddo
        enddo
c
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
        done=1
        pi=4*atan(done)
        nnodes=nphi*ntheta
c
c       
c       Note, in this routine, gridc should be the same as center
c
c       ... first, convert E and H grids to spherical coordinates
c       and perform FFT along each parallel
c
        call em3sphfft(egrid,nphi,ntheta,ts,ws)
        call em3sphfft(hgrid,nphi,ntheta,ts,ws)
c
c
        do 1400 itheta=1,ntheta
c
        dz=radius*rnodes(3,1,itheta)
        r=radius        
c
        costheta=dz/r
        call xulgndr(nterms,costheta,pnm,dnm,pxnm2,punm2)
        theta=acos(costheta)
c
        do n=0,nterms
        do m=0,n
        ynm(n,m)=pnm(n,m)
        xnm2(1,n,m)=pxnm2(1,n,m)
        xnm2(2,n,m)=pxnm2(2,n,m)
        unm2(1,n,m)=punm2(1,n,m)
        unm2(2,n,m)=punm2(2,n,m)
        enddo
        do m=-n,-1
        ynm(n,m)=conjg(ynm(n,-m))
        xnm2(1,n,m)=conjg(pxnm2(1,n,-m))
        xnm2(2,n,m)=conjg(pxnm2(2,n,-m))
        unm2(1,n,m)=conjg(punm2(1,n,-m))
        unm2(2,n,m)=conjg(punm2(2,n,-m))
        enddo
        enddo
c
c
c       ... accumulate a and b multipole expansions
c
c
        do n=1,nterms
        do m=-n,n
        if( m .ge. 0 ) mf=m+1
        if( m .lt. 0 ) mf=nphi+m+1
        ampole(n,m)=ampole(n,m)+ws(itheta)*
     $     (
     $     egrid(1,mf,itheta)*conjg(xnm2(1,n,m))+
     $     egrid(2,mf,itheta)*conjg(xnm2(2,n,m))
     $     )
        bmpole(n,m)=bmpole(n,m)+ws(itheta)*
     $     (
     $     egrid(1,mf,itheta)*conjg(unm2(1,n,m))+
     $     egrid(2,mf,itheta)*conjg(unm2(2,n,m))
     $     )
        enddo
        enddo
c
        do n=1,nterms
        do m=-n,n
        if( m .ge. 0 ) mf=m+1
        if( m .lt. 0 ) mf=nphi+m+1
        ampole(n,m)=ampole(n,m)+ws(itheta)*
     $     (
     $     hgrid(1,mf,itheta)*conjg(unm2(1,n,m))+
     $     hgrid(2,mf,itheta)*conjg(unm2(2,n,m))
     $     )
        bmpole(n,m)=bmpole(n,m)+ws(itheta)*
     $     (
     $     hgrid(1,mf,itheta)*conjg(xnm2(1,n,m))+
     $     hgrid(2,mf,itheta)*conjg(xnm2(2,n,m))
     $     )
        enddo
        enddo
c
 1300   continue
 1400   continue
c
c
c
c       ... now, rescale the multipole expansions
c
        z=rk*radius
c
        if( itype .eq. 1 ) then
c
c       ... evaluate hvals, emthvals, emrhvals
c
        call emhevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c
        if( itype .eq. 2 ) then
c
c       ... evaluate jvals, emtjvals, emrjvals
c
        call emjevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c
cc        call prin2('qvals=*',qvals,2*(nterms+1))
cc        call prin2('emtqvals=*',emtqvals,2*(nterms+1))
c
cc        call em3prinm('ampole, raw=*',ampole,nterms)
cc        call em3prinm('bmpole, raw=*',bmpole,nterms)
c
cc        call prinm1(ampole,nterms,2)
cc        call prinm1(bmpole,nterms,2)
c
cc        call prinml2('ampole, raw=*',ampole,nterms,nterms)
cc        call prinml2('bmpole, raw=*',bmpole,nterms,nterms)
c
        do n=1,nterms
        do m=-n,n
c
        cd=qvals(n)
        cp=emtqvals(n)
        cr=emrqvals(n)
        ampole(n,m)=ampole(n,m)/(cd+cp*ima)/(4*pi)
        bmpole(n,m)=bmpole(n,m)/(cd-cp*ima)/(4*pi)
c
        enddo
        enddo
c
c
ccc        call prin2('ampole=*',ampole,(nterms+1)*(2*nterms+1)*2)
ccc        call prin2('bmpole=*',bmpole,(nterms+1)*(2*nterms+1)*2)
c

        return
        end
c
c
c
c
c
        subroutine em3ehformex_fast2
     $     (itype,rk,gridc,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center,ampole,bmpole,nterms)
        implicit real *8 (a-h,o-z)
c
c       This subroutine forms the outgoing/incoming EM-multipole expansion
c       centered at the center due to the E and H field grids
c
c       Fast O(p^3) routine
c
c       No FFT, evaluation is done in fourier domain
c
c          Input parameters:
c
c       
c       itype - the type of expansion
c          the expansions are outgoing, for itype=1
c          the expansions are incoming, for itype=2
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of dipoles in R^3
c       radius - the radius of E and H grids
c       egrid (complex *16)(3,nphi,ntheta) - E field grid
c       hgrid (complex *16)(3,nphi,ntheta) - H field grid
c       rnodes (real *8)(3,nphi,ntheta) - the grid node coordinates in R^3
c       weights (real *8)(nphi,ntheta) - the grid weights
c       nphi, ntheta - the number of grid points in each direction
c       center - the center of EM-multipole expansion
c
c          Output parameters:
c
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of EM-multipole expansion
c        
        complex *16 rk,z
        dimension xyz(3),center(3)
        dimension gridc(3),target(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c
        complex *16, allocatable :: qvals(:), qders(:)
        complex *16, allocatable :: emtqvals(:), emrqvals(:)
c
        real *8, allocatable :: pnm(:,:), dnm(:,:)
        complex *16, allocatable :: pxnm2(:,:,:)
        complex *16, allocatable :: punm2(:,:,:)
        complex *16, allocatable :: pxnm3(:,:,:)
        complex *16, allocatable :: punm3(:,:,:)
c
        complex *16, allocatable :: ynm(:,:)
        complex *16, allocatable :: xnm3(:,:,:)
        complex *16, allocatable :: unm3(:,:,:)
c
        complex *16, allocatable :: ephi(:)
c
        complex *16, allocatable :: xnm2(:,:,:)
        complex *16, allocatable :: unm2(:,:,:)
c
        real *8, allocatable :: ts(:), ws(:)
c
        complex *16, allocatable :: emgrid(:,:,:)
        complex *16, allocatable :: hmgrid(:,:,:)
c       
        complex *16 cd,cp,cr
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        allocate( qvals(0:nterms) )
        allocate( qders(0:nterms) )
        allocate( emtqvals(0:nterms) )
        allocate( emrqvals(0:nterms) )
        allocate( pnm(0:nterms,0:nterms) )
        allocate( dnm(0:nterms,0:nterms) )
        allocate( pxnm2(2,0:nterms,0:nterms) )
        allocate( punm2(2,0:nterms,0:nterms) ) 
        allocate( pxnm3(3,0:nterms,0:nterms) )
        allocate( punm3(3,0:nterms,0:nterms) ) 
        allocate( ynm(0:nterms,-nterms:nterms) )
        allocate( xnm3(3,0:nterms,-nterms:nterms) )
        allocate( unm3(3,0:nterms,-nterms:nterms) )
        allocate( ephi(-nterms:nterms) )
c
        allocate( xnm2(2,0:nterms,-nterms:nterms) )
        allocate( unm2(2,0:nterms,-nterms:nterms) )       
        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )
c
        allocate( emgrid(3,-nterms:nterms,ntheta) )
        allocate( hmgrid(3,-nterms:nterms,ntheta) )
c
c
        do n=0,nterms
        do m=-n,n
        ampole(n,m)=0
        bmpole(n,m)=0
        enddo
        enddo
c
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
        done=1
        pi=4*atan(done)
        nnodes=nphi*ntheta
c
c       
c       Note, in this routine, gridc should be the same as center
c
c       ... first, convert E and H grids to spherical coordinates
c
        call em3sph(egrid,nphi,ntheta,ts,ws)
        call em3sph(hgrid,nphi,ntheta,ts,ws)
c
c
        do 1400 itheta=1,ntheta
c
        dz=radius*rnodes(3,1,itheta)
        r=radius        
c
        costheta=dz/r
        call xulgndr(nterms,costheta,pnm,dnm,pxnm2,punm2)
        theta=acos(costheta)
c
        do n=0,nterms
        do m=0,n
        ynm(n,m)=pnm(n,m)
        xnm2(1,n,m)=pxnm2(1,n,m)
        xnm2(2,n,m)=pxnm2(2,n,m)
        unm2(1,n,m)=punm2(1,n,m)
        unm2(2,n,m)=punm2(2,n,m)
        enddo
        do m=-n,-1
        ynm(n,m)=conjg(ynm(n,-m))
        xnm2(1,n,m)=conjg(pxnm2(1,n,-m))
        xnm2(2,n,m)=conjg(pxnm2(2,n,-m))
        unm2(1,n,m)=conjg(punm2(1,n,-m))
        unm2(2,n,m)=conjg(punm2(2,n,-m))
        enddo
        enddo
c
c
c       ... accumulate a and b multipole expansions
c
c
        do n=1,nterms
        do m=-n,n
        if( m .ge. 0 ) mf=m+1
        if( m .lt. 0 ) mf=nphi+m+1
        ampole(n,m)=ampole(n,m)+ws(itheta)*
     $     (
     $     egrid(1,mf,itheta)*conjg(xnm2(1,n,m))+
     $     egrid(2,mf,itheta)*conjg(xnm2(2,n,m))
     $     )
        bmpole(n,m)=bmpole(n,m)+ws(itheta)*
     $     (
     $     egrid(1,mf,itheta)*conjg(unm2(1,n,m))+
     $     egrid(2,mf,itheta)*conjg(unm2(2,n,m))
     $     )
        enddo
        enddo
c
        do n=1,nterms
        do m=-n,n
        if( m .ge. 0 ) mf=m+1
        if( m .lt. 0 ) mf=nphi+m+1
        ampole(n,m)=ampole(n,m)+ws(itheta)*
     $     (
     $     hgrid(1,mf,itheta)*conjg(unm2(1,n,m))+
     $     hgrid(2,mf,itheta)*conjg(unm2(2,n,m))
     $     )
        bmpole(n,m)=bmpole(n,m)+ws(itheta)*
     $     (
     $     hgrid(1,mf,itheta)*conjg(xnm2(1,n,m))+
     $     hgrid(2,mf,itheta)*conjg(xnm2(2,n,m))
     $     )
        enddo
        enddo
c
 1300   continue
 1400   continue
c
c
c
c       ... now, rescale the multipole expansions
c
        z=rk*radius
c
        if( itype .eq. 1 ) then
c
c       ... evaluate hvals, emthvals, emrhvals
c
        call emhevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c
        if( itype .eq. 2 ) then
c
c       ... evaluate jvals, emtjvals, emrjvals
c
        call emjevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c
cc        call prin2('qvals=*',qvals,2*(nterms+1))
cc        call prin2('emtqvals=*',emtqvals,2*(nterms+1))
c
cc        call em3prinm('ampole, raw=*',ampole,nterms)
cc        call em3prinm('bmpole, raw=*',bmpole,nterms)
c
c
        do n=1,nterms
        do m=-n,n
c
        cd=qvals(n)
        cp=emtqvals(n)
        cr=emrqvals(n)
        ampole(n,m)=ampole(n,m)/(cd+cp*ima)/(4*pi)
        bmpole(n,m)=bmpole(n,m)/(cd-cp*ima)/(4*pi)
c
        enddo
        enddo
c
c
ccc        call prin2('ampole=*',ampole,(nterms+1)*(2*nterms+1)*2)
ccc        call prin2('bmpole=*',bmpole,(nterms+1)*(2*nterms+1)*2)
c

        return
        end
c
c
c
c
c
        subroutine em3ehformex_fast2_trunc
     $     (itype,rk,gridc,radius,egrid,hgrid,rnodes,weights,
     $     nphi,ntheta,center,ampole,bmpole,nterms,nterms1)
        implicit real *8 (a-h,o-z)
c
c       This subroutine forms the outgoing/incoming EM-multipole expansion
c       centered at the center due to the E and H field grids
c
c       Fast O(p^3) routine
c
c       No FFT, evaluation is done in fourier domain
c
c          Input parameters:
c
c       
c       itype - the type of expansion
c          the expansions are outgoing, for itype=1
c          the expansions are incoming, for itype=2
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of dipoles in R^3
c       radius - the radius of E and H grids
c       egrid (complex *16)(3,nphi,ntheta) - E field grid
c       hgrid (complex *16)(3,nphi,ntheta) - H field grid
c       rnodes (real *8)(3,nphi,ntheta) - the grid node coordinates in R^3
c       weights (real *8)(nphi,ntheta) - the grid weights
c       nphi, ntheta - the number of grid points in each direction
c       center - the center of EM-multipole expansion
c
c          Output parameters:
c
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of EM-multipole expansion
c        
        complex *16 rk,z
        dimension xyz(3),center(3)
        dimension gridc(3),target(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 egrid(3,nphi,ntheta)
        complex *16 hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta)
        dimension weights(nphi,ntheta)
c
        complex *16, allocatable :: qvals(:), qders(:)
        complex *16, allocatable :: emtqvals(:), emrqvals(:)
c
        real *8, allocatable :: pnm(:,:), dnm(:,:)
        complex *16, allocatable :: pxnm2(:,:,:)
        complex *16, allocatable :: punm2(:,:,:)
        complex *16, allocatable :: pxnm3(:,:,:)
        complex *16, allocatable :: punm3(:,:,:)
c
        complex *16, allocatable :: ynm(:,:)
        complex *16, allocatable :: xnm3(:,:,:)
        complex *16, allocatable :: unm3(:,:,:)
c
        complex *16, allocatable :: ephi(:)
c
        complex *16, allocatable :: xnm2(:,:,:)
        complex *16, allocatable :: unm2(:,:,:)
c
        real *8, allocatable :: ts(:), ws(:)
c
        complex *16, allocatable :: emgrid(:,:,:)
        complex *16, allocatable :: hmgrid(:,:,:)
c       
        complex *16 cd,cp,cr
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        allocate( qvals(0:nterms) )
        allocate( qders(0:nterms) )
        allocate( emtqvals(0:nterms) )
        allocate( emrqvals(0:nterms) )
        allocate( pnm(0:nterms,0:nterms) )
        allocate( dnm(0:nterms,0:nterms) )
        allocate( pxnm2(2,0:nterms,0:nterms) )
        allocate( punm2(2,0:nterms,0:nterms) ) 
        allocate( pxnm3(3,0:nterms,0:nterms) )
        allocate( punm3(3,0:nterms,0:nterms) ) 
        allocate( ynm(0:nterms,-nterms:nterms) )
        allocate( xnm3(3,0:nterms,-nterms:nterms) )
        allocate( unm3(3,0:nterms,-nterms:nterms) )
        allocate( ephi(-nterms:nterms) )
c
        allocate( xnm2(2,0:nterms,-nterms:nterms) )
        allocate( unm2(2,0:nterms,-nterms:nterms) )       
        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )
c
        allocate( emgrid(3,-nterms:nterms,ntheta) )
        allocate( hmgrid(3,-nterms:nterms,ntheta) )
c
c
        do n=0,nterms
        do m=-n,n
        ampole(n,m)=0
        bmpole(n,m)=0
        enddo
        enddo
c
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
        done=1
        pi=4*atan(done)
        nnodes=nphi*ntheta
c
c       
c       Note, in this routine, gridc should be the same as center
c
c       ... first, convert E and H grids to spherical coordinates
c
        call em3sph(egrid,nphi,ntheta,ts,ws)
        call em3sph(hgrid,nphi,ntheta,ts,ws)
c
c
        do 1400 itheta=1,ntheta
c
        dz=radius*rnodes(3,1,itheta)
        r=radius        
c
        costheta=dz/r
        call xulgndr_trunc(nterms,nterms1,costheta,pnm,dnm,pxnm2,punm2)
        theta=acos(costheta)
c
        do n=0,nterms
        do m=0,min(n,nterms1)
        ynm(n,m)=pnm(n,m)
        xnm2(1,n,m)=pxnm2(1,n,m)
        xnm2(2,n,m)=pxnm2(2,n,m)
        unm2(1,n,m)=punm2(1,n,m)
        unm2(2,n,m)=punm2(2,n,m)
        enddo
        do m=-min(n,nterms1),-1
        ynm(n,m)=conjg(ynm(n,-m))
        xnm2(1,n,m)=conjg(pxnm2(1,n,-m))
        xnm2(2,n,m)=conjg(pxnm2(2,n,-m))
        unm2(1,n,m)=conjg(punm2(1,n,-m))
        unm2(2,n,m)=conjg(punm2(2,n,-m))
        enddo
        enddo
c
c
c       ... accumulate a and b multipole expansions
c
c
        do n=1,nterms
        do m=-min(n,nterms1),min(n,nterms1)
        if( m .ge. 0 ) mf=m+1
        if( m .lt. 0 ) mf=nphi+m+1
        ampole(n,m)=ampole(n,m)+ws(itheta)*
     $     (
     $     egrid(1,mf,itheta)*conjg(xnm2(1,n,m))+
     $     egrid(2,mf,itheta)*conjg(xnm2(2,n,m))
     $     )
        bmpole(n,m)=bmpole(n,m)+ws(itheta)*
     $     (
     $     egrid(1,mf,itheta)*conjg(unm2(1,n,m))+
     $     egrid(2,mf,itheta)*conjg(unm2(2,n,m))
     $     )
        enddo
        enddo
c
        do n=1,nterms
        do m=-min(n,nterms1),min(n,nterms1)
        if( m .ge. 0 ) mf=m+1
        if( m .lt. 0 ) mf=nphi+m+1
        ampole(n,m)=ampole(n,m)+ws(itheta)*
     $     (
     $     hgrid(1,mf,itheta)*conjg(unm2(1,n,m))+
     $     hgrid(2,mf,itheta)*conjg(unm2(2,n,m))
     $     )
        bmpole(n,m)=bmpole(n,m)+ws(itheta)*
     $     (
     $     hgrid(1,mf,itheta)*conjg(xnm2(1,n,m))+
     $     hgrid(2,mf,itheta)*conjg(xnm2(2,n,m))
     $     )
        enddo
        enddo
c
 1300   continue
 1400   continue
c
c
c
c       ... now, rescale the multipole expansions
c
        z=rk*radius
c
        if( itype .eq. 1 ) then
c
c       ... evaluate hvals, emthvals, emrhvals
c
        call emhevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c
        if( itype .eq. 2 ) then
c
c       ... evaluate jvals, emtjvals, emrjvals
c
        call emjevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c
cc        call prin2('qvals=*',qvals,2*(nterms+1))
cc        call prin2('emtqvals=*',emtqvals,2*(nterms+1))
c
cc        call em3prinm('ampole, raw=*',ampole,nterms)
cc        call em3prinm('bmpole, raw=*',bmpole,nterms)
c
c
        do n=1,nterms
        do m=-min(n,nterms1),min(n,nterms1)
c
        cd=qvals(n)
        cp=emtqvals(n)
        cr=emrqvals(n)
        ampole(n,m)=ampole(n,m)/(cd+cp*ima)/(4*pi)
        bmpole(n,m)=bmpole(n,m)/(cd-cp*ima)/(4*pi)
c
        enddo
        enddo
c
c
ccc        call prin2('ampole=*',ampole,(nterms+1)*(2*nterms+1)*2)
ccc        call prin2('bmpole=*',bmpole,(nterms+1)*(2*nterms+1)*2)
c

        return
        end
c
c
c
c
c
        subroutine em3exevalehz_fast
     $     (itype,rk,ampole,bmpole,nterms,
     $     zshift,zradius,egrid,hgrid,rnodes,weights,nphi,ntheta)
        implicit real *8 (a-h,o-z)
c        
c       This subroutine evaluates E and H field grids 
c       to the EM-multipole expansion (ampole,bmpole) located at the origin
c
c       Fast O(p^3) routine
c
c          Input parameters:
c
c       itype - the type of expansion
c          the expansions are outgoing, for itype=1
c          the expansions are incoming, for itype=2
c       rk (complex *16)  - the frequency parameter
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of EM-multipole expansion
c       nterms - the number of terms in EM-multipole
c       zshift (real *8) - location of the target sphere center in z direction
c       zradius - the radius of target sphere
c       
c       
c
c          Output parameters:
c
c       egrid (complex*16) - the electric field grid at the target
c       hgrid (complex*16) - the magnetic field grid at the target
c
c
        complex *16 rk,z
        dimension xyz(3),center(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 evec(3),hvec(3)
        complex *16 egrid(3,nphi,ntheta),hgrid(3,nphi,ntheta)
c
        complex *16, allocatable :: qvals(:), qders(:)
        complex *16, allocatable :: emtqvals(:), emrqvals(:)
c
        real *8, allocatable :: pnm(:,:), dnm(:,:)
        complex *16, allocatable :: pxnm2(:,:,:)
        complex *16, allocatable :: punm2(:,:,:)
        complex *16, allocatable :: pxnm3(:,:,:)
        complex *16, allocatable :: punm3(:,:,:)
c
        complex *16, allocatable :: ynm(:,:)
        complex *16, allocatable :: xnm3(:,:,:)
        complex *16, allocatable :: unm3(:,:,:)
c
        complex *16, allocatable :: ephi(:)
c
        complex *16, allocatable :: xnm2(:,:,:)
        complex *16, allocatable :: unm2(:,:,:)
c
        real *8, allocatable :: ts(:), ws(:)
c
c
        complex *16 cd,cp,cr
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        allocate( qvals(0:nterms) )
        allocate( qders(0:nterms) )
        allocate( emtqvals(0:nterms) )
        allocate( emrqvals(0:nterms) )
        allocate( pnm(0:nterms,0:nterms) )
        allocate( dnm(0:nterms,0:nterms) )
        allocate( pxnm2(2,0:nterms,0:nterms) )
        allocate( punm2(2,0:nterms,0:nterms) ) 
        allocate( pxnm3(3,0:nterms,0:nterms) )
        allocate( punm3(3,0:nterms,0:nterms) ) 
        allocate( ynm(0:nterms,-nterms:nterms) )
        allocate( xnm3(3,0:nterms,-nterms:nterms) )
        allocate( unm3(3,0:nterms,-nterms:nterms) )
        allocate( ephi(-nterms:nterms) )
c
        allocate( xnm2(2,0:nterms,-nterms:nterms) )
        allocate( unm2(2,0:nterms,-nterms:nterms) )       
        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )
c
c
        do itheta=1,ntheta
        do iphi=1,nphi
        egrid(1,iphi,itheta)=0
        egrid(2,iphi,itheta)=0
        egrid(3,iphi,itheta)=0
        hgrid(1,iphi,itheta)=0
        hgrid(2,iphi,itheta)=0
        hgrid(3,iphi,itheta)=0
        enddo
        enddo
c
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
c
        do 2000 k=1,ntheta
c
c
        z0=zshift+ts(k)*zradius
        x0=sqrt(1-ts(k)*ts(k))*zradius
        r=sqrt(x0*x0+z0*z0)
c
        dx=x0
        dy=0
        dz=z0
c
        if( r .eq. 0 ) then
          rx=1
          ry=0
          rz=0
        else
          rx=dx/r
          ry=dy/r
          rz=dz/r
        endif
c

        z=rk*r
c
        if( itype .eq. 1 ) then
c
c       ... evaluate hvals, emthvals, emrhvals
c
        call emhevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c
        if( itype .eq. 2 ) then
c
c       ... evaluate jvals, emtjvals, emrjvals
c
        call emjevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c
c       ... evaluate xnm3, unm3, ynm, ephi
c
        call em3ctheta(dz,r,costheta)
        call xulgndr(nterms,costheta,pnm,dnm,pxnm2,punm2)
c
        do n=0,nterms
        do m=0,n
        ynm(n,m)=pnm(n,m)
        xnm2(1,n,m)=pxnm2(1,n,m)
        xnm2(2,n,m)=pxnm2(2,n,m)
        unm2(1,n,m)=punm2(1,n,m)
        unm2(2,n,m)=punm2(2,n,m)
        enddo
        do m=-n,-1
        ynm(n,m)=conjg(ynm(n,-m)) 
        xnm2(1,n,m)=conjg(pxnm2(1,n,-m)) 
        xnm2(2,n,m)=conjg(pxnm2(2,n,-m)) 
        unm2(1,n,m)=conjg(punm2(1,n,-m)) 
        unm2(2,n,m)=conjg(punm2(2,n,-m)) 
        enddo
        enddo
c
c
c
c       ... evaluate E and H fields via outgoing a and b multipole expansions
c
c
        do n=1,nterms
        do m=-n,n
c
        cd=qvals(n)
        cp=emtqvals(n)
        cr=emrqvals(n)
c
        if( m .ge. 0 ) mf=m+1
        if( m .lt. 0 ) mf=nphi+m+1
c
        if( mf .ge. 1 .and. mf .le. nphi ) then
        egrid(1,mf,k)=egrid(1,mf,k)+ampole(n,m)*cd*xnm2(1,n,m)
        egrid(2,mf,k)=egrid(2,mf,k)+ampole(n,m)*cd*xnm2(2,n,m)
c
        egrid(1,mf,k)=egrid(1,mf,k)+bmpole(n,m)*cp*unm2(1,n,m)*(-ima)
        egrid(2,mf,k)=egrid(2,mf,k)+bmpole(n,m)*cp*unm2(2,n,m)*(-ima)
c
        egrid(3,mf,k)=egrid(3,mf,k)+bmpole(n,m)*cr*ynm(n,m)*(-ima)
c
        hgrid(1,mf,k)=hgrid(1,mf,k)+bmpole(n,m)*cd*xnm2(1,n,m)
        hgrid(2,mf,k)=hgrid(2,mf,k)+bmpole(n,m)*cd*xnm2(2,n,m)
c
        hgrid(1,mf,k)=hgrid(1,mf,k)+ampole(n,m)*cp*unm2(1,n,m)*(+ima)
        hgrid(2,mf,k)=hgrid(2,mf,k)+ampole(n,m)*cp*unm2(2,n,m)*(+ima)
c
        hgrid(3,mf,k)=hgrid(3,mf,k)+ampole(n,m)*cr*ynm(n,m)*(+ima)
        endif
c
        enddo
        enddo
c
c
 2000   continue
c       
c
        call em3fftcar(zshift,zradius,egrid,nphi,ntheta,ts,ws)
        call em3fftcar(zshift,zradius,hgrid,nphi,ntheta,ts,ws)
c
        return
        end
c
c
c
c
c
        subroutine em3exevalehz_fast2
     $     (itype,rk,ampole,bmpole,nterms,
     $     zshift,zradius,egrid,hgrid,rnodes,weights,nphi,ntheta)
        implicit real *8 (a-h,o-z)
c        
c
c       This subroutine evaluates E and H field grids 
c       to the EM-multipole expansion (ampole,bmpole) located at the origin
c
c       Fast O(p^3) routine
c
c       No FFT, evaluation is done in fourier domain
c
c          Input parameters:
c
c       itype - the type of expansion
c          the expansions are outgoing, for itype=1
c          the expansions are incoming, for itype=2
c       rk (complex *16)  - the frequency parameter
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of EM-multipole expansion
c       nterms - the number of terms in EM-multipole
c       zshift (real *8) - location of the target sphere center in z direction
c       zradius - the radius of target sphere
c       
c       
c
c          Output parameters:
c
c       egrid (complex*16) - the electric field grid at the target
c       hgrid (complex*16) - the magnetic field grid at the target
c
c
        complex *16 rk,z
        dimension xyz(3),center(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 evec(3),hvec(3)
        complex *16 egrid(3,nphi,ntheta),hgrid(3,nphi,ntheta)
c
        complex *16, allocatable :: qvals(:), qders(:)
        complex *16, allocatable :: emtqvals(:), emrqvals(:)
c
        real *8, allocatable :: pnm(:,:), dnm(:,:)
        complex *16, allocatable :: pxnm2(:,:,:)
        complex *16, allocatable :: punm2(:,:,:)
        complex *16, allocatable :: pxnm3(:,:,:)
        complex *16, allocatable :: punm3(:,:,:)
c
        complex *16, allocatable :: ynm(:,:)
        complex *16, allocatable :: xnm3(:,:,:)
        complex *16, allocatable :: unm3(:,:,:)
c
        complex *16, allocatable :: ephi(:)
c
        complex *16, allocatable :: xnm2(:,:,:)
        complex *16, allocatable :: unm2(:,:,:)
c
        real *8, allocatable :: ts(:), ws(:)
c
c
        complex *16 cd,cp,cr
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        allocate( qvals(0:nterms) )
        allocate( qders(0:nterms) )
        allocate( emtqvals(0:nterms) )
        allocate( emrqvals(0:nterms) )
        allocate( pnm(0:nterms,0:nterms) )
        allocate( dnm(0:nterms,0:nterms) )
        allocate( pxnm2(2,0:nterms,0:nterms) )
        allocate( punm2(2,0:nterms,0:nterms) ) 
        allocate( pxnm3(3,0:nterms,0:nterms) )
        allocate( punm3(3,0:nterms,0:nterms) ) 
        allocate( ynm(0:nterms,-nterms:nterms) )
        allocate( xnm3(3,0:nterms,-nterms:nterms) )
        allocate( unm3(3,0:nterms,-nterms:nterms) )
        allocate( ephi(-nterms:nterms) )
c
        allocate( xnm2(2,0:nterms,-nterms:nterms) )
        allocate( unm2(2,0:nterms,-nterms:nterms) )       
        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )
c
c
        do itheta=1,ntheta
        do iphi=1,nphi
        egrid(1,iphi,itheta)=0
        egrid(2,iphi,itheta)=0
        egrid(3,iphi,itheta)=0
        hgrid(1,iphi,itheta)=0
        hgrid(2,iphi,itheta)=0
        hgrid(3,iphi,itheta)=0
        enddo
        enddo
c
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
c
        do 2000 k=1,ntheta
c
c
        z0=zshift+ts(k)*zradius
        x0=sqrt(1-ts(k)*ts(k))*zradius
        r=sqrt(x0*x0+z0*z0)
c
        dx=x0
        dy=0
        dz=z0
c
        if( r .eq. 0 ) then
          rx=1
          ry=0
          rz=0
        else
          rx=dx/r
          ry=dy/r
          rz=dz/r
        endif
c

        z=rk*r
c
        if( itype .eq. 1 ) then
c
c       ... evaluate hvals, emthvals, emrhvals
c
        call emhevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c
        if( itype .eq. 2 ) then
c
c       ... evaluate jvals, emtjvals, emrjvals
c
        call emjevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c
c       ... evaluate xnm3, unm3, ynm, ephi
c
        call em3ctheta(dz,r,costheta)
        call xulgndr(nterms,costheta,pnm,dnm,pxnm2,punm2)
c
        do n=0,nterms
        do m=0,n
        ynm(n,m)=pnm(n,m)
        xnm2(1,n,m)=pxnm2(1,n,m)
        xnm2(2,n,m)=pxnm2(2,n,m)
        unm2(1,n,m)=punm2(1,n,m)
        unm2(2,n,m)=punm2(2,n,m)
        enddo
        do m=-n,-1
        ynm(n,m)=conjg(ynm(n,-m)) 
        xnm2(1,n,m)=conjg(pxnm2(1,n,-m)) 
        xnm2(2,n,m)=conjg(pxnm2(2,n,-m)) 
        unm2(1,n,m)=conjg(punm2(1,n,-m)) 
        unm2(2,n,m)=conjg(punm2(2,n,-m)) 
        enddo
        enddo
c
c
c
c       ... evaluate E and H fields via outgoing a and b multipole expansions
c
c
        do n=1,nterms
        do m=-n,n
c
        cd=qvals(n)
        cp=emtqvals(n)
        cr=emrqvals(n)
c
        if( m .ge. 0 ) mf=m+1
        if( m .lt. 0 ) mf=nphi+m+1
c
        if( mf .ge. 1 .and. mf .le. nphi ) then
        egrid(1,mf,k)=egrid(1,mf,k)+ampole(n,m)*cd*xnm2(1,n,m)
        egrid(2,mf,k)=egrid(2,mf,k)+ampole(n,m)*cd*xnm2(2,n,m)
c
        egrid(1,mf,k)=egrid(1,mf,k)+bmpole(n,m)*cp*unm2(1,n,m)*(-ima)
        egrid(2,mf,k)=egrid(2,mf,k)+bmpole(n,m)*cp*unm2(2,n,m)*(-ima)
c
        egrid(3,mf,k)=egrid(3,mf,k)+bmpole(n,m)*cr*ynm(n,m)*(-ima)
c
        hgrid(1,mf,k)=hgrid(1,mf,k)+bmpole(n,m)*cd*xnm2(1,n,m)
        hgrid(2,mf,k)=hgrid(2,mf,k)+bmpole(n,m)*cd*xnm2(2,n,m)
c
        hgrid(1,mf,k)=hgrid(1,mf,k)+ampole(n,m)*cp*unm2(1,n,m)*(+ima)
        hgrid(2,mf,k)=hgrid(2,mf,k)+ampole(n,m)*cp*unm2(2,n,m)*(+ima)
c
        hgrid(3,mf,k)=hgrid(3,mf,k)+ampole(n,m)*cr*ynm(n,m)*(+ima)
        endif
c
        enddo
        enddo
c
c
 2000   continue
c
c
c       ... finally, convert to cartesian coordinates
c       
        call em3car(zshift,zradius,egrid,nphi,ntheta,ts,ws)
        call em3car(zshift,zradius,hgrid,nphi,ntheta,ts,ws)
c
        return
        end
c
c
c
c
c
        subroutine em3exevalehz_fast2_trunc
     $     (itype,rk,ampole,bmpole,nterms,nterms1,
     $     zshift,zradius,egrid,hgrid,rnodes,weights,nphi,ntheta)
        implicit real *8 (a-h,o-z)
c        
c
c       This subroutine evaluates E and H field grids 
c       to the EM-multipole expansion (ampole,bmpole) located at the origin
c
c       Fast O(p^3) routine
c
c       No FFT, evaluation is done in fourier domain
c
c          Input parameters:
c
c       itype - the type of expansion
c          the expansions are outgoing, for itype=1
c          the expansions are incoming, for itype=2
c       rk (complex *16)  - the frequency parameter
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of EM-multipole expansion
c       nterms - the number of terms in EM-multipole
c       zshift (real *8) - location of the target sphere center in z direction
c       zradius - the radius of target sphere
c       
c       
c
c          Output parameters:
c
c       egrid (complex*16) - the electric field grid at the target
c       hgrid (complex*16) - the magnetic field grid at the target
c
c
        complex *16 rk,z
        dimension xyz(3),center(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 evec(3),hvec(3)
        complex *16 egrid(3,nphi,ntheta),hgrid(3,nphi,ntheta)
c
        complex *16, allocatable :: qvals(:), qders(:)
        complex *16, allocatable :: emtqvals(:), emrqvals(:)
c
        real *8, allocatable :: pnm(:,:), dnm(:,:)
        complex *16, allocatable :: pxnm2(:,:,:)
        complex *16, allocatable :: punm2(:,:,:)
        complex *16, allocatable :: pxnm3(:,:,:)
        complex *16, allocatable :: punm3(:,:,:)
c
        complex *16, allocatable :: ynm(:,:)
        complex *16, allocatable :: xnm3(:,:,:)
        complex *16, allocatable :: unm3(:,:,:)
c
        complex *16, allocatable :: ephi(:)
c
        complex *16, allocatable :: xnm2(:,:,:)
        complex *16, allocatable :: unm2(:,:,:)
c
        real *8, allocatable :: ts(:), ws(:)
c
c
        complex *16 cd,cp,cr
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        allocate( qvals(0:nterms) )
        allocate( qders(0:nterms) )
        allocate( emtqvals(0:nterms) )
        allocate( emrqvals(0:nterms) )
        allocate( pnm(0:nterms,0:nterms) )
        allocate( dnm(0:nterms,0:nterms) )
        allocate( pxnm2(2,0:nterms,0:nterms) )
        allocate( punm2(2,0:nterms,0:nterms) ) 
        allocate( pxnm3(3,0:nterms,0:nterms) )
        allocate( punm3(3,0:nterms,0:nterms) ) 
        allocate( ynm(0:nterms,-nterms:nterms) )
        allocate( xnm3(3,0:nterms,-nterms:nterms) )
        allocate( unm3(3,0:nterms,-nterms:nterms) )
        allocate( ephi(-nterms:nterms) )
c
        allocate( xnm2(2,0:nterms,-nterms:nterms) )
        allocate( unm2(2,0:nterms,-nterms:nterms) )       
        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )
c
c
        do itheta=1,ntheta
        do iphi=1,nphi
        egrid(1,iphi,itheta)=0
        egrid(2,iphi,itheta)=0
        egrid(3,iphi,itheta)=0
        hgrid(1,iphi,itheta)=0
        hgrid(2,iphi,itheta)=0
        hgrid(3,iphi,itheta)=0
        enddo
        enddo
c
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
c
        do 2000 k=1,ntheta
c
c
        z0=zshift+ts(k)*zradius
        x0=sqrt(1-ts(k)*ts(k))*zradius
        r=sqrt(x0*x0+z0*z0)
c
        dx=x0
        dy=0
        dz=z0
c
        if( r .eq. 0 ) then
          rx=1
          ry=0
          rz=0
        else
          rx=dx/r
          ry=dy/r
          rz=dz/r
        endif
c

        z=rk*r
c
        if( itype .eq. 1 ) then
c
c       ... evaluate hvals, emthvals, emrhvals
c
        call emhevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c
        if( itype .eq. 2 ) then
c
c       ... evaluate jvals, emtjvals, emrjvals
c
        call emjevalrt(nterms,z,qvals,qders,emtqvals,emrqvals)
c
        endif
c
c       ... evaluate xnm3, unm3, ynm, ephi
c
        call em3ctheta(dz,r,costheta)
        call xulgndr_trunc(nterms,nterms1,costheta,pnm,dnm,pxnm2,punm2)
c
        do n=0,nterms
        do m=0,min(n,nterms1)
        ynm(n,m)=pnm(n,m)
        xnm2(1,n,m)=pxnm2(1,n,m)
        xnm2(2,n,m)=pxnm2(2,n,m)
        unm2(1,n,m)=punm2(1,n,m)
        unm2(2,n,m)=punm2(2,n,m)
        enddo
        do m=-min(n,nterms1),-1
        ynm(n,m)=conjg(ynm(n,-m)) 
        xnm2(1,n,m)=conjg(pxnm2(1,n,-m)) 
        xnm2(2,n,m)=conjg(pxnm2(2,n,-m)) 
        unm2(1,n,m)=conjg(punm2(1,n,-m)) 
        unm2(2,n,m)=conjg(punm2(2,n,-m)) 
        enddo
        enddo
c
c
c
c       ... evaluate E and H fields via outgoing a and b multipole expansions
c
c
        do n=1,nterms
        do m=-min(n,nterms1),min(n,nterms1)
c
        cd=qvals(n)
        cp=emtqvals(n)
        cr=emrqvals(n)
c
        if( m .ge. 0 ) mf=m+1
        if( m .lt. 0 ) mf=nphi+m+1
c
        if( mf .ge. 1 .and. mf .le. nphi ) then
        egrid(1,mf,k)=egrid(1,mf,k)+ampole(n,m)*cd*xnm2(1,n,m)
        egrid(2,mf,k)=egrid(2,mf,k)+ampole(n,m)*cd*xnm2(2,n,m)
c
        egrid(1,mf,k)=egrid(1,mf,k)+bmpole(n,m)*cp*unm2(1,n,m)*(-ima)
        egrid(2,mf,k)=egrid(2,mf,k)+bmpole(n,m)*cp*unm2(2,n,m)*(-ima)
c
        egrid(3,mf,k)=egrid(3,mf,k)+bmpole(n,m)*cr*ynm(n,m)*(-ima)
c
        hgrid(1,mf,k)=hgrid(1,mf,k)+bmpole(n,m)*cd*xnm2(1,n,m)
        hgrid(2,mf,k)=hgrid(2,mf,k)+bmpole(n,m)*cd*xnm2(2,n,m)
c
        hgrid(1,mf,k)=hgrid(1,mf,k)+ampole(n,m)*cp*unm2(1,n,m)*(+ima)
        hgrid(2,mf,k)=hgrid(2,mf,k)+ampole(n,m)*cp*unm2(2,n,m)*(+ima)
c
        hgrid(3,mf,k)=hgrid(3,mf,k)+ampole(n,m)*cr*ynm(n,m)*(+ima)
        endif
c
        enddo
        enddo
c
c
 2000   continue
c
c
c       ... finally, convert to cartesian coordinates
c       
        call em3car(zshift,zradius,egrid,nphi,ntheta,ts,ws)
        call em3car(zshift,zradius,hgrid,nphi,ntheta,ts,ws)
c
        return
        end
c
c
c
c
c
        subroutine em3mpfareh_fast
     $     (rk,center,ampole,bmpole,nterms,
     $     egrid,hgrid,rnodes,weights,nphi,ntheta)
        entry em3mpfareh
     $     (rk,center,ampole,bmpole,nterms,
     $     egrid,hgrid,rnodes,weights,nphi,ntheta)
        entry em3mpfarehfast
     $     (rk,center,ampole,bmpole,nterms,
     $     egrid,hgrid,rnodes,weights,nphi,ntheta)
        implicit real *8 (a-h,o-z)
c        
c
c       This subroutine evaluates the far field signature E and H field grids 
c       of the EM-multipole expansion (ampole,bmpole) located at the center
c
c       Fast O(p^3) routine
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       center - the location of EM-multipole expansion in R^3
c       ampole (complex*16)(0:nterms,-nterms:nterms)
c       bmpole (complex*16)(0:nterms,-nterms:nterms) 
c           - the coefficients of EM-multipole expansion
c       nterms - the number of terms in EM-multipole
c       
c
c          Output parameters:
c
c       egrid (complex*16) - the electric field grid at the target
c       hgrid (complex*16) - the magnetic field grid at the target
c
c
        complex *16 rk,z
        dimension xyz(3),center(3)
c
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 egrid(3,nphi,ntheta),hgrid(3,nphi,ntheta)
        dimension rnodes(3,nphi,ntheta),weights(nphi,ntheta)
c
        real *8, allocatable :: pnm(:,:), dnm(:,:)
        complex *16, allocatable :: pxnm2(:,:,:)
        complex *16, allocatable :: punm2(:,:,:)
        complex *16, allocatable :: pxnm3(:,:,:)
        complex *16, allocatable :: punm3(:,:,:)
c
        complex *16, allocatable :: ynm(:,:)
        complex *16, allocatable :: xnm3(:,:,:)
        complex *16, allocatable :: unm3(:,:,:)
c
        complex *16, allocatable :: ephi(:)
c
        complex *16, allocatable :: xnm2(:,:,:)
        complex *16, allocatable :: unm2(:,:,:)
c
        real *8, allocatable :: ts(:), ws(:)
c
        dimension far(3)
        complex *16 cfar
c
        complex *16 cd,cp,cr
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        allocate( pnm(0:nterms,0:nterms) )
        allocate( dnm(0:nterms,0:nterms) )
        allocate( pxnm2(2,0:nterms,0:nterms) )
        allocate( punm2(2,0:nterms,0:nterms) ) 
        allocate( pxnm3(3,0:nterms,0:nterms) )
        allocate( punm3(3,0:nterms,0:nterms) ) 
        allocate( ynm(0:nterms,-nterms:nterms) )
        allocate( xnm3(3,0:nterms,-nterms:nterms) )
        allocate( unm3(3,0:nterms,-nterms:nterms) )
        allocate( ephi(-nterms:nterms) )
c
        allocate( xnm2(2,0:nterms,-nterms:nterms) )
        allocate( unm2(2,0:nterms,-nterms:nterms) )       
        allocate( ts(1:ntheta) )
        allocate( ws(1:ntheta) )
c
c
        do itheta=1,ntheta
        do iphi=1,nphi
        egrid(1,iphi,itheta)=0
        egrid(2,iphi,itheta)=0
        egrid(3,iphi,itheta)=0
        hgrid(1,iphi,itheta)=0
        hgrid(2,iphi,itheta)=0
        hgrid(3,iphi,itheta)=0
        enddo
        enddo
c
c
        ifwhts=1
        call legewhts(ntheta,ts,ws,ifwhts)
c
c
        do 2000 k=1,ntheta
c
c
        z0=ts(k)
        x0=sqrt(1-z0*z0)
        r=1
c
        dx=x0
        dy=0
        dz=z0
c
        if( r .eq. 0 ) then
          rx=1
          ry=0
          rz=0
        else
          rx=dx/r
          ry=dy/r
          rz=dz/r
        endif
c
c
c       ... evaluate xnm3, unm3, ynm, ephi
c
        call em3ctheta(dz,r,costheta)
        call xulgndr(nterms,costheta,pnm,dnm,pxnm2,punm2)
c
        do n=0,nterms
        do m=0,n
        ynm(n,m)=pnm(n,m)
        xnm2(1,n,m)=pxnm2(1,n,m)
        xnm2(2,n,m)=pxnm2(2,n,m)
        unm2(1,n,m)=punm2(1,n,m)
        unm2(2,n,m)=punm2(2,n,m)
        enddo
        do m=-n,-1
        ynm(n,m)=conjg(ynm(n,-m)) 
        xnm2(1,n,m)=conjg(pxnm2(1,n,-m)) 
        xnm2(2,n,m)=conjg(pxnm2(2,n,-m)) 
        unm2(1,n,m)=conjg(punm2(1,n,-m)) 
        unm2(2,n,m)=conjg(punm2(2,n,-m)) 
        enddo
        enddo
c
c
c
c       ... evaluate E and H fields via outgoing a and b multipole expansions
c
c
        do n=1,nterms
        do m=-n,n
c
        cd=(-ima)**(n+1) / rk
        cp=ima*(-ima)**(n+1) / rk 
        cr=0
c
        if( m .ge. 0 ) mf=m+1
        if( m .lt. 0 ) mf=nphi+m+1
c
        if( mf .ge. 1 .and. mf .le. nphi ) then
        egrid(1,mf,k)=egrid(1,mf,k)+ampole(n,m)*cd*xnm2(1,n,m)
        egrid(2,mf,k)=egrid(2,mf,k)+ampole(n,m)*cd*xnm2(2,n,m)
c
        egrid(1,mf,k)=egrid(1,mf,k)+bmpole(n,m)*cp*unm2(1,n,m)*(-ima)
        egrid(2,mf,k)=egrid(2,mf,k)+bmpole(n,m)*cp*unm2(2,n,m)*(-ima)
c
        egrid(3,mf,k)=egrid(3,mf,k)+bmpole(n,m)*cr*ynm(n,m)*(-ima)
c
        hgrid(1,mf,k)=hgrid(1,mf,k)+bmpole(n,m)*cd*xnm2(1,n,m)
        hgrid(2,mf,k)=hgrid(2,mf,k)+bmpole(n,m)*cd*xnm2(2,n,m)
c
        hgrid(1,mf,k)=hgrid(1,mf,k)+ampole(n,m)*cp*unm2(1,n,m)*(+ima)
        hgrid(2,mf,k)=hgrid(2,mf,k)+ampole(n,m)*cp*unm2(2,n,m)*(+ima)
c
        hgrid(3,mf,k)=hgrid(3,mf,k)+ampole(n,m)*cr*ynm(n,m)*(+ima)
        endif
c
        enddo
        enddo
c
c
 2000   continue
c
c
        zshift=0
        zradius=1
c
        call em3fftcar(zshift,zradius,egrid,nphi,ntheta,ts,ws)
        call em3fftcar(zshift,zradius,hgrid,nphi,ntheta,ts,ws)
c        
c
        do itheta=1,ntheta
        do iphi=1,nphi
c
        far(1)=rnodes(1,iphi,itheta)
        far(2)=rnodes(2,iphi,itheta)
        far(3)=rnodes(3,iphi,itheta)
c
        d=center(1)*far(1)+center(2)*far(2)+center(3)*far(3)
        cfar=exp(-ima*rk*d)
c
        egrid(1,iphi,itheta)=egrid(1,iphi,itheta)*cfar
        egrid(2,iphi,itheta)=egrid(2,iphi,itheta)*cfar
        egrid(3,iphi,itheta)=egrid(3,iphi,itheta)*cfar
        hgrid(1,iphi,itheta)=hgrid(1,iphi,itheta)*cfar
        hgrid(2,iphi,itheta)=hgrid(2,iphi,itheta)*cfar
        hgrid(3,iphi,itheta)=hgrid(3,iphi,itheta)*cfar
        enddo
        enddo
c
        return
        end
c
c
c
c
c
