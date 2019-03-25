c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       This is the end of the debugging routines and the beginning of
c       the code for the evaluation of monochromatic plane waves for the
c       Maxwell's equations in R^3
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine emplanew(rk,xyz,cjvec,cmvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields of the sum of three
c       simple polarized monochromatic plane waves
c
c       The waves are assumed to be of simple form: a sum of three
c       linearly polarized plane waves coming from either positive x, y,
c       or z directions.
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric vector
c       cmvec (complex *16) - the strength of the magnetic vector   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
c
        dimension xyz(3)
        complex *16 cjvec(3),cmvec(3),evec(3),hvec(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
c
        x=xyz(1)
        y=xyz(2)
        z=xyz(3)
c
c
c       (sum of of three linearly polarized waves from positive x, y, z
c       directions, electric part)
c
c       rkvec = (0,-1,0)*rk, epol = (cjvec(2),0,0)
c       rkvec = (0,0,-1)*rk, epol = (0,cjvec(3),0)
c       rkvec = (-1,0,0)*rk, epol = (0,0,cjvec(1))
c
        evec(1)=exp(-ima*rk*y)*cjvec(2)
        evec(2)=exp(-ima*rk*z)*cjvec(3)
        evec(3)=exp(-ima*rk*x)*cjvec(1)
c
        hvec(1)=exp(-ima*rk*z)*cjvec(3)
        hvec(2)=exp(-ima*rk*x)*cjvec(1)
        hvec(3)=exp(-ima*rk*y)*cjvec(2)
c
c
c       (sum of three linearly polarized waves from positive x, y, z
c       directions magnetic part)
c
c       rkvec = (0,0,-1)*rk, epol = (0,0,cmvec(3))
c       rkvec = (-1,0,0)*rk, epol = (cmvec(1),0,0)
c       rkvec = (0,-1,0)*rk, epol = (0,cmvec(2),0)
c
        evec(1)=evec(1)+exp(-ima*rk*z)*cmvec(3)
        evec(2)=evec(2)+exp(-ima*rk*x)*cmvec(1)
        evec(3)=evec(3)+exp(-ima*rk*y)*cmvec(2)
c
        hvec(1)=hvec(1)-exp(-ima*rk*y)*cmvec(2)
        hvec(2)=hvec(2)-exp(-ima*rk*z)*cmvec(3)
        hvec(3)=hvec(3)-exp(-ima*rk*x)*cmvec(1)
c
        return
        end
c
c
c
c
c
        subroutine emplanenw(rk,xyz,cjvec,cmvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields of the sum of three
c       simple polarized monochromatic plane waves
c
c       The waves are assumed to be of simple form: a sum of three
c       linearly polarized plane waves coming from either negative x, y,
c       or z directions.
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric vector
c       cmvec (complex *16) - the strength of the magnetic vector   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
c
        dimension xyz(3)
        complex *16 cjvec(3),cmvec(3),evec(3),hvec(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
c
        x=xyz(1)
        y=xyz(2)
        z=xyz(3)
c
c
c       (sum of of three linearly polarized waves from negative x, y, z
c       directions, electric part)
c
c       rkvec = (0,+1,0)*rk, epol = (cjvec(2),0,0)
c       rkvec = (0,0,+1)*rk, epol = (0,cjvec(3),0)
c       rkvec = (+1,0,0)*rk, epol = (0,0,cjvec(1))
c
c
        evec(1)=exp(ima*rk*y)*cjvec(2)
        evec(2)=exp(ima*rk*z)*cjvec(3)
        evec(3)=exp(ima*rk*x)*cjvec(1)
c
        hvec(1)=-exp(ima*rk*z)*cjvec(3)
        hvec(2)=-exp(ima*rk*x)*cjvec(1)
        hvec(3)=-exp(ima*rk*y)*cjvec(2)
c
c
c       (sum of three linearly polarized waves from positive x, y, z
c       directions magnetic part)
c
c       rkvec = (0,0,+1)*rk, epol = (0,0,cmvec(3))
c       rkvec = (+1,0,0)*rk, epol = (cmvec(1),0,0)
c       rkvec = (0,+1,0)*rk, epol = (0,cmvec(2),0)
c
c
        evec(1)=evec(1)+exp(ima*rk*z)*cmvec(3)
        evec(2)=evec(2)+exp(ima*rk*x)*cmvec(1)
        evec(3)=evec(3)+exp(ima*rk*y)*cmvec(2)
c
        hvec(1)=hvec(1)+exp(ima*rk*y)*cmvec(2)
        hvec(2)=hvec(2)+exp(ima*rk*z)*cmvec(3)
        hvec(3)=hvec(3)+exp(ima*rk*x)*cmvec(1)
c

        return
        end
c
c
c
c
c
        subroutine emplanearb(rkvec,epol,xyz,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields of an arbitrary
c       polarized monochromatic plane wave. Note that epol must be
c       orthogonal to rkvec, otherwise result will not a valid EM field.
c
c          Input parameters:
c
c       rkvec (complex *16) - the k vector
c       epol (complex *16) - the strength of the E polarization, 
c          
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
c
        dimension xyz(3)
        complex *16 rkvec(3),epol(3),hpol(3),evec(3),hvec(3),rk,ima
        complex *16 cd
c
        data ima/(0.0d0,1.0d0)/
c
        x=xyz(1)
        y=xyz(2)
        z=xyz(3)
c
c       rkvec \cdot E = 0, rkvec \cdot H = 0
c
        cd = rkvec(1)*x+rkvec(2)*y+rkvec(3)*z
c
        evec(1)=exp(ima*cd)*epol(1)
        evec(2)=exp(ima*cd)*epol(2)
        evec(3)=exp(ima*cd)*epol(3)
c 
c        hpol = rkvec x epol 
c
        hpol(1)=rkvec(2)*epol(3)-rkvec(3)*epol(2)
        hpol(2)=rkvec(3)*epol(1)-rkvec(1)*epol(3)
        hpol(3)=rkvec(1)*epol(2)-rkvec(2)*epol(1)
c
c       normalize hpol vector, warning, complex valued rk can be
c       difficult to interpret
c
        rk=sqrt(rkvec(1)*rkvec(1)+rkvec(2)*rkvec(2)+rkvec(3)*rkvec(3))
c
        hvec(1)=exp(ima*cd)*hpol(1)/rk
        hvec(2)=exp(ima*cd)*hpol(2)/rk
        hvec(3)=exp(ima*cd)*hpol(3)/rk
c 
        return
        end
c
c
c
c
c
        subroutine emplanelw(rk,xyz,itype,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields of the linearly
c       polarized monochromatic plane waves.  The waves are assumed to
c       be of simple form (coming from either x, y, or z directions with
c       no phase shift).
c       
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       itype (integer) - the type of polarization
c         itype=1 linearly polarized from x-direction (z-plane polarization)
c         itype=2 linearly polarized from y-direction (x-plane polarization)
c         itype=3 linearly polarized from z-direction (y-plane polarization)
c         itype=4 linearly polarized from x-direction (y-plane polarization)
c         itype=5 linearly polarized from y-direction (z-plane polarization)
c         itype=6 linearly polarized from z-direction (x-plane polarization)
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3)
        complex *16 cjvec(3),cmvec(3),evec(3),hvec(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
c
        x=xyz(1)
        y=xyz(2)
        z=xyz(3)
c
c       
        if( itype .eq. 1 ) then
c
c       (linearly polarized wave from x-direction)
c       rkvec = (-1,0,0)*rk, epol = (0,0,1)
c
        evec(1)=0
        evec(2)=0
        evec(3)=exp(-ima*rk*x)
c
        hvec(1)=0
        hvec(2)=exp(-ima*rk*x)  
        hvec(3)=0
c
        endif
c
c
        if( itype .eq. 2 ) then
c
c       (linearly polarized wave from y-direction)
c       rkvec = (0,-1,0)*rk, epol = (1,0,0)
c
        evec(1)=exp(-ima*rk*y)
        evec(2)=0
        evec(3)=0
c
        hvec(1)=0
        hvec(2)=0
        hvec(3)=exp(-ima*rk*y)  
c
        endif
c
c
        if( itype .eq. 3 ) then
c
c       (linearly polarized wave from z-direction)
c       rkvec = (0,0,-1)*rk, epol = (0,1,0)
c
        evec(1)=0
        evec(2)=exp(-ima*rk*z)
        evec(3)=0
c
        hvec(1)=exp(-ima*rk*z)  
        hvec(2)=0
        hvec(3)=0
c
        endif
c
c
c
c       
        if( itype .eq. 4 ) then
c
c       (linearly polarized wave from x-direction)
c       rkvec = (-1,0,0)*rk, epol = (0,1,0)
c
        evec(1)=0
        evec(2)=exp(-ima*rk*x)
        evec(3)=0
c
        hvec(1)=0
        hvec(2)=0
        hvec(3)=-exp(-ima*rk*x)  
c
        endif
c
c
        if( itype .eq. 5 ) then
c
c       (linearly polarized wave from y-direction)
c       rkvec = (0,-1,0)*rk, epol = (0,0,1)
c
        evec(1)=0
        evec(2)=0
        evec(3)=exp(-ima*rk*y)
c
        hvec(1)=-exp(-ima*rk*y)  
        hvec(2)=0
        hvec(3)=0
c
        endif
c
c
        if( itype .eq. 6 ) then
c
c       (linearly polarized wave from z-direction)
c       rkvec = (0,0,-1)*rk, epol = (1,0,0)
c
        evec(1)=exp(-ima*rk*z)
        evec(2)=0
        evec(3)=0
c
        hvec(1)=0
        hvec(2)=-exp(-ima*rk*z)  
        hvec(3)=0
c
        endif
c
c
        return
        end
c
c
c
c
c
        subroutine emplanecw(rk,xyz,itype,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields of the circularly
c       polarized monochromatic plane waves.  The waves are assumed to
c       be of simple form (coming from either x, y, or z directions with
c       no phase shift).
c       
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       itype (integer) - the type of polarization
c         itype=1 circularly polarized from x-direction 
c         itype=2 circularly polarized from y-direction 
c         itype=3 circularly polarized from z-direction 
c         itype=4 circularly polarized from x-direction 
c         itype=5 circularly polarized from y-direction 
c         itype=6 circularly polarized from z-direction 
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3)
        complex *16 cjvec(3),cmvec(3),evec(3),hvec(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
c
        x=xyz(1)
        y=xyz(2)
        z=xyz(3)
c
c       
        if( itype .eq. 1 ) then
c
c       (circularly polarized wave from x-direction)
c       rkvec = (-1,0,0)*rk, epol = (0,1,ima)
c
        evec(1)=0
        evec(2)=exp(-ima*rk*x)
        evec(3)=ima*exp(-ima*rk*x)
c
        hvec(1)=0
        hvec(2)=ima*exp(-ima*rk*x)  
        hvec(3)=-exp(-ima*rk*x)
c
        endif
c
c
        if( itype .eq. 2 ) then
c
c       (circularly polarized wave from y-direction)
c       rkvec = (0,-1,0)*rk, epol = (ima,0,1)
c
        evec(1)=ima*exp(-ima*rk*y)
        evec(2)=exp(-ima*rk*y)
        evec(3)=0
c
        hvec(1)=-exp(-ima*rk*y)  
        hvec(2)=0
        hvec(3)=ima*exp(-ima*rk*y)
c
        endif
c
c
        if( itype .eq. 3 ) then
c
c       (circularly polarized wave from z-direction)
c       rkvec = (0,0,-1)*rk, epol = (1,ima,0)
c
        evec(1)=exp(-ima*rk*z)
        evec(2)=ima*exp(-ima*rk*z)
        evec(3)=0
c
        hvec(1)=ima*exp(-ima*rk*z)  
        hvec(2)=-exp(-ima*rk*z)
        hvec(3)=0
c
        endif
c
c
c
c       
        if( itype .eq. 4 ) then
c
c       (circularly polarized wave from x-direction)
c       rkvec = (-1,0,0)*rk, epol = (0,1,-ima)
c
        evec(1)=0
        evec(2)=exp(-ima*rk*x)
        evec(3)=-ima*exp(-ima*rk*x)
c
        hvec(1)=0
        hvec(2)=-ima*exp(-ima*rk*x)
        hvec(3)=-exp(-ima*rk*x)  
c
        endif
c
c
        if( itype .eq. 5 ) then
c
c       (circularly polarized wave from y-direction)
c       rkvec = (0,-1,0)*rk, epol = (-ima,0,1)
c
        evec(1)=-ima*exp(-ima*rk*y)
        evec(2)=0
        evec(3)=exp(-ima*rk*y)
c
        hvec(1)=-exp(-ima*rk*y)  
        hvec(2)=0
        hvec(3)=-ima*exp(-ima*rk*y)
c
        endif
c
c
        if( itype .eq. 6 ) then
c
c       (circularly polarized wave from z-direction)
c       rkvec = (0,0,-1)*rk, epol = (1,-ima,0)
c
        evec(1)=exp(-ima*rk*z)
        evec(2)=-ima*exp(-ima*rk*z)
        evec(3)=0
c
        hvec(1)=-ima*exp(-ima*rk*z)  
        hvec(2)=-exp(-ima*rk*z)  
        hvec(3)=0
c
        endif
c
c
        return
        end
c
c
c
c
c
        subroutine emplanene(evec,hvec)
        implicit real *8 (a-h,o-z)
        complex *16 evec(3),hvec(3)
c
c       ... reverse the direction of a plane wave
c       
        evec(1)=conjg(evec(1))
        evec(2)=conjg(evec(2))
        evec(3)=conjg(evec(3))
c
        hvec(1)=-conjg(hvec(1))
        hvec(2)=-conjg(hvec(2))
        hvec(3)=-conjg(hvec(3))
c
        return
        end
