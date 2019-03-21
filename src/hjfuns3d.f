c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       This is the end of the debugging routines and the beginning of
c       the code for the evaluation of the spherical Bessel functions
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c        This file contains a set of subroutines for the handling 
c        of spherical Bessel functions. It contains 2 subroutines that
c        are user-callable. Following is a brief description of these 
c        subroutines.
c
c     cjfuns3d - evaluates first n+1 spherical j Bessel functions of 
c          the complex argument z
c
c     chfuns3d - evaluates first n+1 spherical h Bessel functions of
c          the complex argument z
c
c
c
        subroutine cjfuns3d(ier,n,z,scale,fjs,ifder,fjder,work,lwork)
        implicit real *8 (a-h,o-z)
c       
        complex *16 fjs(0:1),fjder(0:1),z
        complex *16 work(1)
c       
        ier=0
c       
        lwfjs=2*(n+1)
c       
c       
c       ... estimate lwfjs and ntop in a loop (maybe large in some cases)
c       
        do 1000 k=1,10
c       
        ifjs=1
        iiscale=1+lwfjs
        lused=iiscale+lwfjs
c
        if( lused .gt. lwork ) then
c
c     ... out of memory, return
c
        ier=8
        return
        endif
c
        nterms=n
        call jfuns3d(ier2,nterms,z,scale,work(ifjs),ifder,fjder,
     $     lwfjs,work(iiscale),ntop)
c
ccc      call prinf('ier2=*',ier2,1)
c
ccc      write(*,*) ier2
        if( ier2 .eq. 0 ) then
c
c     ... success, exit from the loop
c
        ier=0
        goto 1100
        else
c
c     ... lwfjs is too small, increase it by a factor 1.5
c
        ier=1
        lwfjs=lwfjs*1.5
        endif
c
 1000   continue
c     
c     ...
c
 1100   continue
c       
        do 1200 i=0,n
        fjs(i)=work(i+1)
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
        subroutine chfuns3d(ier,n,z,scale,fhs,ifder,fhder,work,lwork)
        implicit real *8 (a-h,o-z)
c       
        complex *16 fjs(0:1),fjder(0:1),z
        complex *16 work(1)
c       
        ier=0
c       
        call hfuns3d(n,z,scale,fhs,ifder,fhder)
c       
        return
        end
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine jfuns3d(ier,nterms,z,scale,fjs,ifder,fjder,
     $     lwfjs,iscale,ntop)
        implicit real *8 (a-h,o-z)
c
c PURPOSE:
c
c	This subroutine evaluates the first NTERMS spherical Bessel 
c	functions and if required, their derivatives.
c	It incorporates a scaling parameter SCALE so that
c       
c		fjs_n(z)=j_n(z)/SCALE^n
c		fjder_n(z)=\frac{\partial fjs_n(z)}{\partial z}
c
c	NOTE: The scaling parameter SCALE is meant to be used when
c             abs(z) < 1, in which case we recommend setting
c	      SCALE = abs(z). This prevents the fjs_n from 
c             underflowing too rapidly.
c	      Otherwise, set SCALE=1.
c
c INPUT:
c
c    nterms (integer *4): order of expansion of output array fjs 
c    z     (complex *16): argument of the spherical Bessel functions
c    scale    (real *8) : scaling factor (discussed above)
c    ifder  (integer *4): flag indicating whether to calculate "fjder"
c		          0	NO
c		          1	YES
c    lwfjs  (integer *4): upper limit of input arrays 
c                         fjs(0:1) and iscale(0:1)
c    iscale (integer *4): integer workspace used to keep track of 
c                         internal scaling
c
c OUTPUT:
c
c    ier    (integer *4): error return code 
c                         ier=0 normal return;
c                         ier=8 insufficient array dimension lwfjs
c    fjs   (complex *16): array of scaled Bessel functions.
c    fjder (complex *16): array of derivs of scaled Bessel functions.
c    ntop  (integer *4) : highest index in arrays fjs that is nonzero
c
c
      integer iscale(0:1)
      complex *16 wavek,fjs(0:1),fjder(0:1)
      complex *16 z,zinv,com,fj0,fj1,zscale,ztmp
c
      DATA UPBOUND/1.0d+32/, UPBOUND2/1.0d+40/, UPBOUND2inv/1.0d-40/
      DATA TINY/1.0d-14/,done/1.0d0/,zero/0.0d0/
c
c ... Initializing ...
c
      ier=0
c
c       set to asymptotic values if argument is sufficiently small
c
      if (abs(z).lt.TINY) then
         fjs(0) = done
         do i = 1, nterms
            fjs(i) = zero
	 enddo
c
	 if (ifder.eq.1) then
	    do i=0,nterms
	       fjder(i)=zero
	    enddo
	    fjder(1)=done/(3*scale)
	 endif
c
         RETURN
      endif
c
c ... Step 1: recursion up to find ntop, starting from nterms
c
      ntop=0
      zinv=done/z
      fjs(nterms)=done
      fjs(nterms-1)=zero
c
c     fix potential memory over-run error
c
cc      do 1200 i=nterms,lwfjs
      do 1200 i=nterms,lwfjs-1
         dcoef=2*i+done
         ztmp=dcoef*zinv*fjs(i)-fjs(i-1)
         fjs(i+1)=ztmp
c
         dd = dreal(ztmp)**2 + dimag(ztmp)**2
ccc	   call prin2('dd=*',dd,1)
         if (dd .gt. UPBOUND) then
            ntop=i+1
            goto 1300
         endif
 1200 continue
 1300 continue
      if (ntop.eq.0) then
         ier=8
         return
      endif
cccc        call prinf('ntop=*',ntop,1)
c
c ... Step 2: Recursion back down to generate the unscaled jfuns:
c             if magnitude exceeds UPBOUND2, rescale and continue the 
c	      recursion (saving the order at which rescaling occurred 
c	      in array iscale.
c
      do i=0,ntop
         iscale(i)=0
      enddo
c
      fjs(ntop)=zero
      fjs(ntop-1)=done
      do 2200 i=ntop-1,1,-1
	 dcoef=2*i+done
         ztmp=dcoef*zinv*fjs(i)-fjs(i+1)
         fjs(i-1)=ztmp
c
         dd = dreal(ztmp)**2 + dimag(ztmp)**2
         if (dd.gt.UPBOUND2) then
            fjs(i) = fjs(i)*UPBOUND2inv
            fjs(i-1) = fjs(i-1)*UPBOUND2inv
            iscale(i) = 1
         endif
 2200 continue
ccc	call prin2('fjs before scaling*',fjs(0),ntop*2+2)
cccc	call prinf(' iscale is *',iscale,ntop+1)
c
c ...  Step 3: go back up to the top and make sure that all
c              Bessel functions are scaled by the same factor
c              (i.e. the net total of times rescaling was invoked
c              on the way down in the previous loop).
c              At the same time, add scaling to fjs array.
c
      ncntr=0
      scalinv=done/scale
      sctot = 1.0d0
      do i=1,ntop
         sctot = sctot*scalinv
         if(iscale(i-1).eq.1) sctot=sctot*UPBOUND2inv
         fjs(i)=fjs(i)*sctot
      enddo
ccc      call prin2('fjs after unif scaling*',fjs(0),ntop*2+2)
c
c ... Determine the normalization parameter:
c
      fj0=sin(z)*zinv
      fj1=fj0*zinv-cos(z)*zinv
c
      d0=abs(fj0)
      d1=abs(fj1)
      if (d1 .gt. d0) then
ccc         zscale=fj1/fjs(1)
         zscale=fj1/(fjs(1)*scale)
      else
         zscale=fj0/fjs(0)
      endif
cccc	call prin2('zscale=*',zscale,2)
c
c ... Scale the jfuns by zscale:
c
      ztmp=zscale
      do i=0,nterms
         fjs(i)=fjs(i)*ztmp
      enddo
ccc	call prin2('fjs after scaling*',fjs(0),nterms*2+2)
c
c ... Finally, calculate the derivatives if desired:
c
      if (ifder.eq.1) then
         fjs(nterms+1)=fjs(nterms+1)*ztmp
c
         fjder(0)=-fjs(1)*scale
         do i=1,nterms
            dc1=i/(2*i+done)
            dc2=done-dc1
            dc1=dc1*scalinv
            dc2=dc2*scale
            fjder(i)=dc1*fjs(i-1)-dc2*fjs(i+1)
         enddo
ccc	   call prin2('fjder after scaling*',fjder(0),n*2+2)
      endif
      return
      end
c
c
c
c
c
********************************************************************
      subroutine hfuns3d(n,z,scale,fhs,ifder,fhder)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c Purpose:
c	This subroutine calculates spherical hankel
c	functions of orders 0 through n of the complex
c	argument z. And if required, the derivative of the hankel 
c	functions. NOTE: actually, it can return a scaled version of
c	the hankel functions when |z| < 1:
c		fhs_n(z)=h_n(z)*scale^(n+1)
c		fhder_n(z)=\frac{\partial fhs_n(z)}{\partial z}
c
c	NOTE: the scale only meant to be used when
c			abs(z) < 1, 
c	      in which case set (most cases ||wavek||),
c			scale \aprox abs(z),
c		   (this is to prevent h blows up too fast)
c	      otherwise, set
c			scale=1.
c
c Input:
c	ifder - flag indcating whether to calculate fhder
c		0---NO
c		1---YES
c	n - the number of the hankel functions to be returned
c	z - the complex number where to generate the 
c		spherical h-functions.
c
c Output:
c	fhs - the spherical hankel functions at z.
c	fhder - the derivatives of spherical hankel functions at z.
c          WARNING: BOTH ARRAYS MUST BE DIMENSIONED AS fhs(0:n)!!!
c
      complex *16 fhs(0:1),fhder(0:1)
      complex *16 zk2,z,zinv,ztmp,fhextra
c
      data thresh/1.0d-16/,done/1.0d0/
c
c ... Initializing, checking for error: if d is too small, error!
c
      if (abs(z).lt.THRESH) then
         do i=0,n
            fhs(i)=0
            fhder(i)=0
         enddo
         RETURN
      endif
c
c ... Otherwise, proceed accordingly ...
c
      scal2=scale*scale
      call hank3d01(z,fhs(0),fhs(1))
      fhs(0)=fhs(0)*scale
      fhs(1)=fhs(1)*scal2
c
c ... Calculate the rest of the scaled hankel functions: the formula is
c	(10.1.19)  \hat{h}_{n+1}(k*r)=(2n+1)/r * \hat{h}_n(k*r) - 
c			k**2 * \hat{h}_{n-1}(k*r)
c
      zinv=done/z
      do 1200 i=1,n-1
	 dtmp=scale*(2*i+done)
	 ztmp=zinv*dtmp
	 fhs(i+1)=ztmp*fhs(i)-scal2*fhs(i-1)
 1200 continue
cccc	call prin2('In hfuns3d, fhs=*',fhs(0),n*2+2)
c
c ... Now calculate the derivatives: the used formula is (10.1.21)
c	\hat{fhder}_{n}(k*r)=k**2 * \hat{h}_{n-1}(k*r)
c		- (n+1)/r * \hat{h}_n(k*r)
c
      if (ifder.eq.1) then
c
c	compute fhs(n+1) first:
c
	 dtmp=scale*(2*n+done)
	 ztmp=zinv*dtmp
	 fhextra=ztmp*fhs(n)-scal2*fhs(n-1)
c
c	handle the two ends:
c
	 fhder(0)=-fhs(1)/scale
c
	 scalinv=done/scale
	 cc1=n/(2*n+done)
	 cc2=done-cc1
	 cc1=cc1*scale
	 cc2=cc2*scalinv
	 fhder(n)=cc1*fhs(n-1)-cc2*fhextra
c
c	now the bunch
c
         do i=1,n-1
	    cc1=i/(2*i+done)
	    cc2=done-cc1
	    cc1=cc1*scale
	    cc2=cc2*scalinv
	    fhder(i)=cc1*fhs(i-1)-cc2*fhs(i+1)
	 enddo
      endif
c
      return
      end
c
c
c
c
c
      subroutine hank3d01(z,h0,h1)
      implicit real *8 (a-h,o-z)
c
c Purpose:
c	This subroutine calculates h0 and h1 - the spherical
c	hankel functions of orders 0, 1:
c		h0 = exp(i*z)/(i*z),
c		h1 = exp(i*z)*(-i/z**2-1/z)
c
c Input:
c	z: complex number where to calculate h0
c	  NOTE: when abs(z) < 1.0d-300, this subroutine returns zero.
c
c Output:
c	h0: the hankel function of order 0.
c	h1: the derivative of the hankel function.
c
        complex *16 z,zinv,ima,cd,h0,h1
        data ima/(0.0d0,1.0d0)/
c
c       ... if z is too small, simply return 0:
c
ccc	call prin2(' in hank3d01, z=*',z,2)
c
        done=1
c
        if( abs(z) .lt. 1d-300 ) then
c        
        h0=0.0d0
        h1=0.0d0
c
        return
        endif
c
c
c       ... Regular computation for h0 and h1:
c
        zinv=done/z
        cd=exp(z*ima)
        h0=-ima*cd*zinv
        h1=h0*(zinv-ima)
c
ccc        call prin2(' in hank3d01, h0=*',h0,2)
ccc        call prin2(' in hank3d01, h1=*',h1,2)
c
        return
        end
c
c
c
c
c
        subroutine hfuns3d012(z,h0,h1,h2)
        implicit real *8 (a-h,o-z)
c
c Purpose:
c	This subroutine calculates h0, h1,and h2 - the spherical
c	hankel functions of orders 0, 1, and 2:
c		h0 = exp(i*z)*(-i/z),
c		h1 = exp(i*z)*(-1/z-i/z**2)
c		h2 = exp(i*z)*(+i/z-3/z**2-3*i/z**3)
c
c Input:
c	z: complex number where to calculate h0
c	  NOTE: when abs(z) < 1.0d-300, this subroutine returns zero.
c
c Output:
c	h0: the spherical hankel function of order 0.
c	h1: the spherical hankel function of order 1.
c	h2: the spherical hankel function of order 2
c
        complex *16 z,zinv,ima,cd,h0,h1,h2
        data ima/(0.0d0,1.0d0)/
c       
c       ... if z is too small, simply return 0:
c
ccc	call prin2(' in hfuns3d012, z=*',z,2)
c
        done=1
c
        if( abs(z) .lt. 1d-300 ) then
c 
        h0=0.0d0
        h1=0.0d0
        h2=0.0d0
c
        return
        endif
c
c       ... Regular computation for h0, h1 and h2:
c
        zinv=done/z
        cd=exp(z*ima)
c
        h0=cd*(-ima*zinv)
        h1=cd*(-zinv-ima*zinv**2)
        h2=cd*(+ima*zinv-3*zinv**2-3*ima*zinv**3)
c       
        h0=cd*(-ima*zinv)
        h1=cd*((-1-ima*zinv)*zinv)
        h2=cd*((+ima-(3+3*ima*zinv)*zinv)*zinv)
c       
ccc        call prin2(' in hfuns3d012, h0=*',h0,2)
ccc        call prin2(' in hfuns3d012, h1=*',h1,2)
ccc        call prin2(' in hfuns3d012, h2=*',h2,2)
c       
        return
        end
c
c
c
c
