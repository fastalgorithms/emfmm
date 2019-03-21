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
c    $Date: 2011-10-17 13:56:46 -0400 (Mon, 17 Oct 2011) $
c    $Revision: 2386 $
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c       This is the end of the debugging code and the beginning 
c       of the electromagnetic multipole expansion rotation routines
c       
c       This version is compatible with the notation _without_ (-1)**m
c       scaling factor
c
c       Maxwell's equations in R^3
c
c       Fortran 90 version, allocatable arrays
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       This file contains a set of subroutines for the handling of
c       electromagnetic field multipole expansions. It contains 6
c       use-callable subroutines. Following is a brief description of
c       these subroutines.
c
c
c   emabrotyzf - rotate the z-axis oriented em-multipole
c       expansion to point into the direction dir
c
c   emabrotyzb - rotates the dir-oriented em-multipole
c       expansion to point into the direction z
c
c   emabrotaf - perform the forward Euler rotation of the em-multipole 
c
c   emabrotab - perform the backward Euler rotation of the em-multipole 
c
c
c
        subroutine emabrotyzf(nterms,dir,coefs,couts)
        implicit real *8 (a-h,o-z)
c
c       This subroutine rotates the z-axis oriented em-multipole
c       expansion to point into the direction dir. Such a transformation
c       has one degree of freedom and therefore is intended to be used
c       together with emabrotyzb
c
c       Input parameters:
c
c       nterms - the number of terms in em-multipole expansion
c       dir - the new orientation direction, real *8 (3)
c       coefs - the coefficients of em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
c       Output parameters:
c
c       couts - the coefficients of rotated em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
c
        dimension dir(3)
        complex *16 coefs(0:nterms,-nterms:nterms)
        complex *16 couts(0:nterms,-nterms:nterms)
        complex *16, allocatable :: mpole2(:,:)
        allocatable w(:)
c
c       em-fields
c       rotate em-multipole expansion 
c
c
        x=dir(1)
        y=dir(2)
        z=dir(3)
c
        r=sqrt(x*x+y*y+z*z)
c
        x=x/r
        y=y/r
        z=z/r
c       
        if( r .eq. 0 ) then
          phi=0
          theta=0
        else
          if( abs(x) .eq. 0 .and. abs(y) .eq. 0 ) then
          phi = 0
          else
          phi = atan2(y,x)
          endif
          theta=acos(z)
        endif
c
        allocate( mpole2(0:nterms,-nterms:nterms) ) 
        lw = 2*(nterms+1)*(2*nterms+1) + 2*(2*nterms+1) + 7
        allocate(w(lw),stat=ier)
        if( ier .ne. 0 ) then
        call prinf('in emabrotyzf, allocate, lw=*',lw,1)
        endif
c
        call mpolerotz(nterms,coefs,phi,mpole2)
        call mpoleroty(nterms,mpole2,theta,couts,w,lw,lused)
c
        return
        end
c
c
c
c
c
        subroutine emabrotyzb(nterms,dir,coefs,couts)
        implicit real *8 (a-h,o-z)
c
c       This subroutine rotates the dir-oriented em-multipole
c       expansion to point into the direction z. Such a transformation has
c       one degree of freedom and therefore is intended to be used
c       together with emabrotyzf
c
c       Input parameters:
c
c       nterms - the number of terms in em-multipole expansion
c       dir - the new orientation direction, real *8 (3)
c       coefs - the coefficients of em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
c       Output parameters:
c
c       couts - the coefficients of rotated em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
        dimension dir(3)
        complex *16 coefs(0:nterms,-nterms:nterms)
        complex *16 couts(0:nterms,-nterms:nterms)
        complex *16, allocatable :: mpole2(:,:)
        allocatable w(:)
c
c       em-fields
c       rotate back
c       
        x=dir(1)
        y=dir(2)
        z=dir(3)
c
        r=sqrt(x*x+y*y+z*z)
c
        x=x/r
        y=y/r
        z=z/r
c        
        if( r .eq. 0 ) then
          phi=0
          theta=0
        else
          if( abs(x) .eq. 0 .and. abs(y) .eq. 0 ) then
          phi = 0
          else
          phi = atan2(y,x)
          endif
          theta=acos(z)
        endif
c
        allocate( mpole2(0:nterms,-nterms:nterms) ) 
        lw = 2*(nterms+1)*(2*nterms+1) + 2*(2*nterms+1) + 7
        allocate(w(lw),stat=ier)
        if( ier .ne. 0 ) then
        call prinf('in emabrotyzb, allocate lw=*',lw,1)
        endif
c
        call mpoleroty(nterms,coefs,-theta,mpole2,w,lw,lused)
        call mpolerotz(nterms,mpole2,-phi,couts)
c
        return
        end
c
c
c
c
c
        subroutine emabrotab(nterms,rota,coefs,couts)
        implicit real *8 (a-h,o-z)
c
c       This subroutine performs the forward 
c       Euler rotations the em-multipole expansion.
c
c       Input parameters:
c
c       nterms - the number of terms in em-multipole expansion
c       rota - the Euler angles, real *8 (3)
c       coefs - the coefficients of em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
c       Output parameters:
c
c       couts - the coefficients of rotated em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
        dimension rota(3)
        complex *16 coefs(0:nterms,-nterms:nterms)
        complex *16 couts(0:nterms,-nterms:nterms)
        complex *16, allocatable :: mpole1(:,:)
        complex *16, allocatable :: mpole2(:,:)
        allocatable w(:)
c
c       em-fields
c       rotate em-multipole expansion 
c
        allocate( mpole1(0:nterms,-nterms:nterms) )
        allocate( mpole2(0:nterms,-nterms:nterms) ) 
        lw = 2*(nterms+1)*(2*nterms+1) + 2*(2*nterms+1) + 7
        allocate(w(lw),stat=ier)
        if( ier .ne. 0 ) then
        call prinf('in emabrotab, allocate, lw=*',lw,1)        
        endif
c
        call mpolerotz(nterms,coefs,rota(3),mpole2)
        call mpoleroty(nterms,mpole2,rota(2),mpole1,w,lw,lused)
        call mpolerotz(nterms,mpole1,rota(1),couts)
c
        return
        end
c
c
c
c
c
        subroutine emabrotaf(nterms,rota,coefs,couts)
        implicit real *8 (a-h,o-z)
c
c       This subroutine performs the backward
c       Euler rotations the em-multipole expansion.
c
c       Input parameters:
c
c       nterms - the number of terms in em-multipole expansion
c       rota - the Euler angles, real *8 (3)
c       coefs - the coefficients of em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
c       Output parameters:
c
c       couts - the coefficients of rotated em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
        dimension rota(3)
        complex *16 coefs(0:nterms,-nterms:nterms)
        complex *16 couts(0:nterms,-nterms:nterms)
        complex *16, allocatable :: mpole1(:,:)
        complex *16, allocatable :: mpole2(:,:)
        allocatable w(:)
c
c       em-fields
c       rotate em-multipole expansion 
c
        allocate( mpole1(0:nterms,-nterms:nterms) )
        allocate( mpole2(0:nterms,-nterms:nterms) ) 
        lw = 2*(nterms+1)*(2*nterms+1) + 2*(2*nterms+1) + 7
        allocate(w(lw),stat=ier)
        if( ier .ne. 0 ) then
        call prinf('in emabrotaf, allocate, lw=*',lw,1)
        endif
c
        call mpolerotz(nterms,coefs,-rota(1),mpole2)
        call mpoleroty(nterms,mpole2,-rota(2),mpole1,w,lw,lused)
        call mpolerotz(nterms,mpole1,-rota(3),couts)
c       
        return
        end
c
c
c
c
c
        subroutine mpolerotz(nterms,mpole,phi,mpout)
        implicit real *8 (a-h,o-z)
c
c       This subroutine performs the rotation of the multipole expansion
c       about the z-axis
c
c     Input parameters:
c
c     nterms  (integer   )  : number of terms in the multipole expansion
c     mpole  (complex *16)  : the multipole expansion
c     phi        (real *8)  : the angle of rotation
c
c     Output parameters:
c 
c     mpout  (complex *16)  : the rotated multipole expansion
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 mpout(0:nterms,-nterms:nterms)
        complex *16 ephi,ima
        data ima/(0.0d0,1.0d0)/
c
c
c----- a rotation of PHI radians about the Z-axis 
c
c
        do l=0,nterms
        do m=-l,l
c
        ephi=exp(ima*m*phi)
        mpout(l,m)=ephi*mpole(l,m)
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
        subroutine mpoleroty(nterms,mpole,theta,mpout,w,lw,lused)
        implicit real *8 (a-h,o-z)
c
c       This subroutine performs the rotation of the multipole expansion
c       about the y-axis
c
c     Input parameters:
c
c     nterms  (integer   )  : number of terms in the multipole expansion
c     mpole  (complex *16)  : the multipole expansion
c     theta
c
c     Output parameters:
c 
c     mpout  (complex *16)  : the rotated multipole expansion
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 mpout(0:nterms,-nterms:nterms)
c
c
c----- a rotation of THETA radians about the Y-axis 
c
        call rotviarecur3(theta,nterms,nterms,nterms,mpole,nterms,
     $     mpout,nterms,w,lw,lused)
c
        return
        end
c
c
c
c
c       
        subroutine emabrotyzf_trunc(nterms,m2,dir,coefs,couts)
        implicit real *8 (a-h,o-z)
c
c       This subroutine rotates the z-axis oriented em-multipole
c       expansion to point into the direction dir. Such a transformation
c       has one degree of freedom and therefore is intended to be used
c       together with emabrotyzb
c
c       Input parameters:
c
c       nterms - the number of terms in em-multipole expansion
c       dir - the new orientation direction, real *8 (3)
c       coefs - the coefficients of em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
c       Output parameters:
c
c       couts - the coefficients of rotated em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
c
        dimension dir(3)
        complex *16 coefs(0:nterms,-nterms:nterms)
        complex *16 couts(0:nterms,-nterms:nterms)
        complex *16, allocatable :: mpole2(:,:)
        allocatable w(:)
c
c       em-fields
c       rotate em-multipole expansion 
c
c
        x=dir(1)
        y=dir(2)
        z=dir(3)
c
        r=sqrt(x*x+y*y+z*z)
c
        x=x/r
        y=y/r
        z=z/r
c       
        if( r .eq. 0 ) then
          phi=0
          theta=0
        else
          if( abs(x) .eq. 0 .and. abs(y) .eq. 0 ) then
          phi = 0
          else
          phi = atan2(y,x)
          endif
          theta=acos(z)
        endif
c
        allocate( mpole2(0:nterms,-nterms:nterms) ) 
        lw = 2*(nterms+1)*(2*nterms+1) + 2*(2*nterms+1) + 7
        allocate(w(lw),stat=ier)
        if( ier .ne. 0 ) then
        call prinf('in emabrotyzf, allocate, lw=*',lw,1)
        endif
c
        call mpolerotz(nterms,coefs,phi,mpole2)
        call mpolerotyvar2(nterms,m2,mpole2,theta,couts,w,lw,lused)
c
        return
        end
c
c
c
c
c
        subroutine emabrotyzb_trunc(nterms,m1,dir,coefs,couts)
        implicit real *8 (a-h,o-z)
c
c       This subroutine rotates the dir-oriented em-multipole
c       expansion to point into the direction z. Such a transformation has
c       one degree of freedom and therefore is intended to be used
c       together with emabrotyzf
c
c       Input parameters:
c
c       nterms - the number of terms in em-multipole expansion
c       dir - the new orientation direction, real *8 (3)
c       coefs - the coefficients of em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
c       Output parameters:
c
c       couts - the coefficients of rotated em-multipole expansion
c                complex *16(0:nterms,-nterms:nterms)
c
        dimension dir(3)
        complex *16 coefs(0:nterms,-nterms:nterms)
        complex *16 couts(0:nterms,-nterms:nterms)
        complex *16, allocatable :: mpole2(:,:)
        allocatable w(:)
c
c       em-fields
c       rotate back
c       
        x=dir(1)
        y=dir(2)
        z=dir(3)
c
        r=sqrt(x*x+y*y+z*z)
c
        x=x/r
        y=y/r
        z=z/r
c        
        if( r .eq. 0 ) then
          phi=0
          theta=0
        else
          if( abs(x) .eq. 0 .and. abs(y) .eq. 0 ) then
          phi = 0
          else
          phi = atan2(y,x)
          endif
          theta=acos(z)
        endif
c
        allocate( mpole2(0:nterms,-nterms:nterms) ) 
        lw = 2*(nterms+1)*(2*nterms+1) + 2*(2*nterms+1) + 7
        allocate(w(lw),stat=ier)
        if( ier .ne. 0 ) then
        call prinf('in emabrotyzb, allocate lw=*',lw,1)
        endif
c
        call mpolerotyvar1(nterms,m1,coefs,-theta,mpole2,w,lw,lused)
        call mpolerotz(nterms,mpole2,-phi,couts)
c
        return
        end
c
c
c
c
c
        subroutine mpolerotyvar1(nterms,m1,mpole,theta,mpout,w,lw,lused)
        implicit real *8 (a-h,o-z)
c
c       This subroutine performs the rotation of the multipole expansion
c       about the y-axis
c
c     Input parameters:
c
c     nterms  (integer   )  : number of terms in the multipole expansion
c     mpole  (complex *16)  : the multipole expansion
c     theta
c
c     Output parameters:
c 
c     mpout  (complex *16)  : the rotated multipole expansion
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 mpout(0:nterms,-nterms:nterms)
c
c
c----- a rotation of THETA radians about the Y-axis 
c
c        call rotviarecur3(theta,nterms,m1,nterms,mpole,nterms,
c     $     mpout,nterms,w,lw,lused)
c
        call rotviarecur3f90(theta,nterms,m1,nterms,mpole,nterms,
     $     mpout,nterms)
c
c        call rotviaprojvarf90(theta,nterms,m1,nterms,mpole,nterms,
c     $     mpout,nterms)
c
        return
        end
c
c
c
c
c       
        subroutine mpolerotyvar2(nterms,m2,mpole,theta,mpout,w,lw,lused)
        implicit real *8 (a-h,o-z)
c
c       This subroutine performs the rotation of the multipole expansion
c       about the y-axis
c
c     Input parameters:
c
c     nterms  (integer   )  : number of terms in the multipole expansion
c     mpole  (complex *16)  : the multipole expansion
c     theta
c
c     Output parameters:
c 
c     mpout  (complex *16)  : the rotated multipole expansion
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 mpout(0:nterms,-nterms:nterms)
c
c
c----- a rotation of THETA radians about the Y-axis 
c
c        call rotviarecur3(theta,nterms,nterms,m2,mpole,nterms,
c     $     mpout,nterms,w,lw,lused)
c
        call rotviarecur3f90(theta,nterms,nterms,m2,mpole,nterms,
     $     mpout,nterms)
c
c        call rotviaprojvarf90(theta,nterms,nterms,m2,mpole,nterms,
c     $     mpout,nterms)
c
        return
        end
c
c
c
c
c       
