
        implicit real *8 (a-h,o-z)
        complex *16 rk,ima,cd
        dimension source(3),target(3),xyz(3)
        complex *16 cjvec(3),cmvec(3)
        complex *16 evec(3),hvec(3),pvec(3)
        dimension dir(3),far(3)
        complex *16 efar(3),hfar(3),pfar(3)
        complex *16 rkvec(3),epol(3)
        data ima/(0.0d0,1.0d0)/

c
c       SET ALL PARAMETERS
c        
        call prini(6,13)
        
        done=1
        pi=4*atan(done)

        rk=1
        call prin2('rk=*',rk,2)

        source(1)=0
        source(2)=0
        source(3)=0

        target(1)=1
        target(2)=2
        target(3)=3

        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)


        call prin2('source=*',source,3)
        call prin2('target=*',target,3)

        call prin2('= electric and magnetic dipoles =*',cjvec,0)

        cjvec(1)=0
        cjvec(2)=0
        cjvec(3)=1

        call prin2('cjvec=*',cjvec,2*3)

        call dipole3e(rk,xyz,cjvec,evec,hvec)
        call empoynting(evec,hvec,pvec)
        call prin2('evec=*',evec,2*3)
        call prin2('hvec=*',hvec,2*3)
        call prin2('pvec=*',pvec,2*3)


        cmvec(1)=0
        cmvec(2)=0
        cmvec(3)=1

        call prin2('cmvec=*',cmvec,2*3)

        call dipole3m(rk,xyz,cmvec,evec,hvec)
        call empoynting(evec,hvec,pvec)

        call prin2('evec=*',evec,2*3)
        call prin2('hvec=*',hvec,2*3)
        call prin2('pvec=*',pvec,2*3)

        call prin2('= far field signatures =*',cjvec,0)

        call prin2('cjvec=*',cjvec,2*3)

        source(1)=1
        source(2)=0
        source(3)=0

        dir(1)=1
        dir(2)=1
        dir(3)=1
        d=sqrt(dir(1)**2+dir(2)**2+dir(3)**2)
        far(1)=dir(1)/d
        far(2)=dir(2)/d
        far(3)=dir(3)/d
        
        target(1)=1d6*far(1)
        target(2)=1d6*far(2)
        target(3)=1d6*far(3)

        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)

        d=sqrt(target(1)**2+target(2)**2+target(3)**2)
        cd=exp(ima*rk*d)/(d)

        call prin2('source=*',source,3)
        call prin2('dir=*',dir,3)

        call dipole3et(rk,source,target,cjvec,evec,hvec)
        evec(1)=evec(1)/cd
        evec(2)=evec(2)/cd
        evec(3)=evec(3)/cd
        call prin2('evec / (exp(ikr)/r)=*',evec,2*3)
        hvec(1)=hvec(1)/cd
        hvec(2)=hvec(2)/cd
        hvec(3)=hvec(3)/cd
        call prin2('hvec / (exp(ikr)/r)=*',hvec,2*3)

        
        call dipole3efar(rk,source,cjvec,dir,efar,hfar)
        call empoynting(efar,hfar,pfar)
        call prin2('efar=*',efar,2*3)
        call prin2('hfar=*',hfar,2*3)
        call prin2('pfar=*',pfar,2*3)


        call prin2('cmvec=*',cmvec,2*3)

        call prin2('source=*',source,3)
        call prin2('dir=*',dir,3)

        call dipole3mt(rk,source,target,cmvec,evec,hvec)
        evec(1)=evec(1)/cd
        evec(2)=evec(2)/cd
        evec(3)=evec(3)/cd
        call prin2('evec / (exp(ikr)/r)=*',evec,2*3)
        hvec(1)=hvec(1)/cd
        hvec(2)=hvec(2)/cd
        hvec(3)=hvec(3)/cd
        call prin2('hvec / (exp(ikr)/r)=*',hvec,2*3)

        
        call dipole3mfar(rk,source,cmvec,dir,efar,hfar)
        call empoynting(efar,hfar,pfar)
        call prin2('efar=*',efar,2*3)
        call prin2('hfar=*',hfar,2*3)
        call prin2('pfar=*',pfar,2*3)

        call prin2('= divergence =*',cjvec,0)

        source(1)=0
        source(2)=0
        source(3)=0

        target(1)=1
        target(2)=2
        target(3)=3

        call prin2('cjvec=*',cjvec,2*3)
        call dipole3e_check_div(rk,source,target,cjvec)

        call prin2('cmvec=*',cjvec,2*3)
        call dipole3m_check_div(rk,source,target,cjvec)

        call prin2('= curl =*',cjvec,0)

        source(1)=0
        source(2)=0
        source(3)=0

        target(1)=1
        target(2)=2
        target(3)=3

        call prin2('cjvec=*',cjvec,2*3)
        call dipole3e_check_curl(rk,source,target,cjvec)

        call prin2('cmvec=*',cjvec,2*3)
        call dipole3m_check_curl(rk,source,target,cjvec)

        itype=6
        call prinf('emplanelw, itype=*',itype,1)
        call emplanelw_check_curl(rk,target,itype)

        rkvec(1)=1
        rkvec(2)=1
        rkvec(3)=1
        epol(1)=1
        epol(2)=-2
        epol(3)=1
        call prin2('emplanearb, rkvec=*',rkvec,6)        
        call prin2('emplanearb, epol=*',epol,6)        
        cd = rkvec(1)*epol(1)+rkvec(2)*epol(2)+rkvec(3)*epol(3)
        call prin2('emplanearb, rkvec dot epol=*',cd,2)        
        call emplanearb_check_curl(rkvec,epol,target)

        stop
        end

c
c
c
c
c
        subroutine dipole3e_check_div(rk,source,target,cjvec)
        implicit real *8 (a-h,o-z)
        dimension source(3),target(3),target_new(3)
        complex *16 rk,cjvec(3)
        complex *16 evec1(3),hvec1(3)
        complex *16 evec2(3),hvec2(3)
        complex *16 eder(3),hder(3),ediv,hdiv

        h=1d-6
        ediv=0
        hdiv=0

        target_new(1)=target(1)+h
        target_new(2)=target(2)
        target_new(3)=target(3)        
        call dipole3et(rk,source,target_new,cjvec,evec1,hvec1)
        target_new(1)=target(1)-h
        target_new(2)=target(2)
        target_new(3)=target(3)        
        call dipole3et(rk,source,target_new,cjvec,evec2,hvec2)
        eder(1)=(evec1(1)-evec2(1))/(2*h)
        eder(2)=(evec1(2)-evec2(2))/(2*h)
        eder(3)=(evec1(3)-evec2(3))/(2*h)
        hder(1)=(hvec1(1)-hvec2(1))/(2*h)
        hder(2)=(hvec1(2)-hvec2(2))/(2*h)
        hder(3)=(hvec1(3)-hvec2(3))/(2*h)

        ediv=ediv+eder(1)
        hdiv=hdiv+hder(1)

        target_new(1)=target(1)
        target_new(2)=target(2)+h
        target_new(3)=target(3)        
        call dipole3et(rk,source,target_new,cjvec,evec1,hvec1)
        target_new(1)=target(1)
        target_new(2)=target(2)-h
        target_new(3)=target(3)        
        call dipole3et(rk,source,target_new,cjvec,evec2,hvec2)
        eder(1)=(evec1(1)-evec2(1))/(2*h)
        eder(2)=(evec1(2)-evec2(2))/(2*h)
        eder(3)=(evec1(3)-evec2(3))/(2*h)
        hder(1)=(hvec1(1)-hvec2(1))/(2*h)
        hder(2)=(hvec1(2)-hvec2(2))/(2*h)
        hder(3)=(hvec1(3)-hvec2(3))/(2*h)

        ediv=ediv+eder(2)
        hdiv=hdiv+hder(2)

        target_new(1)=target(1)
        target_new(2)=target(2)
        target_new(3)=target(3)+h        
        call dipole3et(rk,source,target_new,cjvec,evec1,hvec1)
        target_new(1)=target(1)
        target_new(2)=target(2)
        target_new(3)=target(3)-h        
        call dipole3et(rk,source,target_new,cjvec,evec2,hvec2)
        eder(1)=(evec1(1)-evec2(1))/(2*h)
        eder(2)=(evec1(2)-evec2(2))/(2*h)
        eder(3)=(evec1(3)-evec2(3))/(2*h)
        hder(1)=(hvec1(1)-hvec2(1))/(2*h)
        hder(2)=(hvec1(2)-hvec2(2))/(2*h)
        hder(3)=(hvec1(3)-hvec2(3))/(2*h)

        ediv=ediv+eder(3)
        hdiv=hdiv+hder(3)

        call prin2('ediv=*',ediv,2)
        call prin2('hdiv=*',hdiv,2)

        return
        end
c
c
c
c
c
        subroutine dipole3m_check_div(rk,source,target,cmvec)
        implicit real *8 (a-h,o-z)
        dimension source(3),target(3),target_new(3)
        complex *16 rk,cmvec(3)
        complex *16 evec1(3),hvec1(3)
        complex *16 evec2(3),hvec2(3)
        complex *16 eder(3),hder(3),ediv,hdiv

        h=1d-6
        ediv=0
        hdiv=0

        target_new(1)=target(1)+h
        target_new(2)=target(2)
        target_new(3)=target(3)        
        call dipole3mt(rk,source,target_new,cmvec,evec1,hvec1)
        target_new(1)=target(1)-h
        target_new(2)=target(2)
        target_new(3)=target(3)        
        call dipole3mt(rk,source,target_new,cmvec,evec2,hvec2)
        eder(1)=(evec1(1)-evec2(1))/(2*h)
        eder(2)=(evec1(2)-evec2(2))/(2*h)
        eder(3)=(evec1(3)-evec2(3))/(2*h)
        hder(1)=(hvec1(1)-hvec2(1))/(2*h)
        hder(2)=(hvec1(2)-hvec2(2))/(2*h)
        hder(3)=(hvec1(3)-hvec2(3))/(2*h)

        ediv=ediv+eder(1)
        hdiv=hdiv+hder(1)

        target_new(1)=target(1)
        target_new(2)=target(2)+h
        target_new(3)=target(3)        
        call dipole3mt(rk,source,target_new,cmvec,evec1,hvec1)
        target_new(1)=target(1)
        target_new(2)=target(2)-h
        target_new(3)=target(3)        
        call dipole3mt(rk,source,target_new,cmvec,evec2,hvec2)
        eder(1)=(evec1(1)-evec2(1))/(2*h)
        eder(2)=(evec1(2)-evec2(2))/(2*h)
        eder(3)=(evec1(3)-evec2(3))/(2*h)
        hder(1)=(hvec1(1)-hvec2(1))/(2*h)
        hder(2)=(hvec1(2)-hvec2(2))/(2*h)
        hder(3)=(hvec1(3)-hvec2(3))/(2*h)

        ediv=ediv+eder(2)
        hdiv=hdiv+hder(2)

        target_new(1)=target(1)
        target_new(2)=target(2)
        target_new(3)=target(3)+h        
        call dipole3mt(rk,source,target_new,cmvec,evec1,hvec1)
        target_new(1)=target(1)
        target_new(2)=target(2)
        target_new(3)=target(3)-h        
        call dipole3mt(rk,source,target_new,cmvec,evec2,hvec2)
        eder(1)=(evec1(1)-evec2(1))/(2*h)
        eder(2)=(evec1(2)-evec2(2))/(2*h)
        eder(3)=(evec1(3)-evec2(3))/(2*h)
        hder(1)=(hvec1(1)-hvec2(1))/(2*h)
        hder(2)=(hvec1(2)-hvec2(2))/(2*h)
        hder(3)=(hvec1(3)-hvec2(3))/(2*h)

        ediv=ediv+eder(3)
        hdiv=hdiv+hder(3)

        call prin2('ediv=*',ediv,2)
        call prin2('hdiv=*',hdiv,2)

        return
        end
c
c
c
c
c
        subroutine dipole3e_check_curl(rk,source,target,cjvec)
        implicit real *8 (a-h,o-z)
        dimension source(3),target(3),target_new(3)
        complex *16 rk,ima,cjvec(3)
        complex *16 evec(3),hvec(3)
        complex *16 evec1(3),hvec1(3)
        complex *16 evec2(3),hvec2(3)
        complex *16 ederx(3),hderx(3)
        complex *16 edery(3),hdery(3)
        complex *16 ederz(3),hderz(3)
        complex *16 ecurl(3),hcurl(3)
        data ima/(0.0d0,1.0d0)/

        h=1d-6
        ediv=0
        hdiv=0

        target_new(1)=target(1)+h
        target_new(2)=target(2)
        target_new(3)=target(3)        
        call dipole3et(rk,source,target_new,cjvec,evec1,hvec1)
        target_new(1)=target(1)-h
        target_new(2)=target(2)
        target_new(3)=target(3)        
        call dipole3et(rk,source,target_new,cjvec,evec2,hvec2)
        ederx(1)=(evec1(1)-evec2(1))/(2*h)
        ederx(2)=(evec1(2)-evec2(2))/(2*h)
        ederx(3)=(evec1(3)-evec2(3))/(2*h)
        hderx(1)=(hvec1(1)-hvec2(1))/(2*h)
        hderx(2)=(hvec1(2)-hvec2(2))/(2*h)
        hderx(3)=(hvec1(3)-hvec2(3))/(2*h)

        target_new(1)=target(1)
        target_new(2)=target(2)+h
        target_new(3)=target(3)        
        call dipole3et(rk,source,target_new,cjvec,evec1,hvec1)
        target_new(1)=target(1)
        target_new(2)=target(2)-h
        target_new(3)=target(3)        
        call dipole3et(rk,source,target_new,cjvec,evec2,hvec2)
        edery(1)=(evec1(1)-evec2(1))/(2*h)
        edery(2)=(evec1(2)-evec2(2))/(2*h)
        edery(3)=(evec1(3)-evec2(3))/(2*h)
        hdery(1)=(hvec1(1)-hvec2(1))/(2*h)
        hdery(2)=(hvec1(2)-hvec2(2))/(2*h)
        hdery(3)=(hvec1(3)-hvec2(3))/(2*h)

        target_new(1)=target(1)
        target_new(2)=target(2)
        target_new(3)=target(3)+h        
        call dipole3et(rk,source,target_new,cjvec,evec1,hvec1)
        target_new(1)=target(1)
        target_new(2)=target(2)
        target_new(3)=target(3)-h        
        call dipole3et(rk,source,target_new,cjvec,evec2,hvec2)
        ederz(1)=(evec1(1)-evec2(1))/(2*h)
        ederz(2)=(evec1(2)-evec2(2))/(2*h)
        ederz(3)=(evec1(3)-evec2(3))/(2*h)
        hderz(1)=(hvec1(1)-hvec2(1))/(2*h)
        hderz(2)=(hvec1(2)-hvec2(2))/(2*h)
        hderz(3)=(hvec1(3)-hvec2(3))/(2*h)
c
c   see http://en.wikipedia.org/wiki/Curl_(mathematics)
c
        ecurl(1)=edery(3)-ederz(2)
        ecurl(2)=ederz(1)-ederx(3)
        ecurl(3)=ederx(2)-edery(1)

        hcurl(1)=hdery(3)-hderz(2)
        hcurl(2)=hderz(1)-hderx(3)
        hcurl(3)=hderx(2)-hdery(1)

        call dipole3et(rk,source,target,cjvec,evec,hvec)
        evec(1)=evec(1)*(-ima*rk)
        evec(2)=evec(2)*(-ima*rk)
        evec(3)=evec(3)*(-ima*rk)
        hvec(1)=hvec(1)*(+ima*rk)
        hvec(2)=hvec(2)*(+ima*rk)
        hvec(3)=hvec(3)*(+ima*rk)

        call prin2('ecurl=*',ecurl,2*3)
        call prin2('hvec \times (+ima rk)=*',hvec,2*3)
        call prin2('hcurl=*',hcurl,2*3)
        call prin2('evec \times (-ima rk)=*',evec,2*3)


        return
        end
c
c
c
c
c
        subroutine dipole3m_check_curl(rk,source,target,cmvec)
        implicit real *8 (a-h,o-z)
        dimension source(3),target(3),target_new(3)
        complex *16 rk,ima,cmvec(3)
        complex *16 evec(3),hvec(3)
        complex *16 evec1(3),hvec1(3)
        complex *16 evec2(3),hvec2(3)
        complex *16 ederx(3),hderx(3)
        complex *16 edery(3),hdery(3)
        complex *16 ederz(3),hderz(3)
        complex *16 ecurl(3),hcurl(3)
        data ima/(0.0d0,1.0d0)/

        h=1d-6
        ediv=0
        hdiv=0

        target_new(1)=target(1)+h
        target_new(2)=target(2)
        target_new(3)=target(3)        
        call dipole3mt(rk,source,target_new,cmvec,evec1,hvec1)
        target_new(1)=target(1)-h
        target_new(2)=target(2)
        target_new(3)=target(3)        
        call dipole3mt(rk,source,target_new,cmvec,evec2,hvec2)
        ederx(1)=(evec1(1)-evec2(1))/(2*h)
        ederx(2)=(evec1(2)-evec2(2))/(2*h)
        ederx(3)=(evec1(3)-evec2(3))/(2*h)
        hderx(1)=(hvec1(1)-hvec2(1))/(2*h)
        hderx(2)=(hvec1(2)-hvec2(2))/(2*h)
        hderx(3)=(hvec1(3)-hvec2(3))/(2*h)

        target_new(1)=target(1)
        target_new(2)=target(2)+h
        target_new(3)=target(3)        
        call dipole3mt(rk,source,target_new,cmvec,evec1,hvec1)
        target_new(1)=target(1)
        target_new(2)=target(2)-h
        target_new(3)=target(3)        
        call dipole3mt(rk,source,target_new,cmvec,evec2,hvec2)
        edery(1)=(evec1(1)-evec2(1))/(2*h)
        edery(2)=(evec1(2)-evec2(2))/(2*h)
        edery(3)=(evec1(3)-evec2(3))/(2*h)
        hdery(1)=(hvec1(1)-hvec2(1))/(2*h)
        hdery(2)=(hvec1(2)-hvec2(2))/(2*h)
        hdery(3)=(hvec1(3)-hvec2(3))/(2*h)

        target_new(1)=target(1)
        target_new(2)=target(2)
        target_new(3)=target(3)+h        
        call dipole3mt(rk,source,target_new,cmvec,evec1,hvec1)
        target_new(1)=target(1)
        target_new(2)=target(2)
        target_new(3)=target(3)-h        
        call dipole3mt(rk,source,target_new,cmvec,evec2,hvec2)
        ederz(1)=(evec1(1)-evec2(1))/(2*h)
        ederz(2)=(evec1(2)-evec2(2))/(2*h)
        ederz(3)=(evec1(3)-evec2(3))/(2*h)
        hderz(1)=(hvec1(1)-hvec2(1))/(2*h)
        hderz(2)=(hvec1(2)-hvec2(2))/(2*h)
        hderz(3)=(hvec1(3)-hvec2(3))/(2*h)
c
c   see http://en.wikipedia.org/wiki/Curl_(mathematics)
c
        ecurl(1)=edery(3)-ederz(2)
        ecurl(2)=ederz(1)-ederx(3)
        ecurl(3)=ederx(2)-edery(1)

        hcurl(1)=hdery(3)-hderz(2)
        hcurl(2)=hderz(1)-hderx(3)
        hcurl(3)=hderx(2)-hdery(1)

        call dipole3mt(rk,source,target,cmvec,evec,hvec)
        evec(1)=evec(1)*(-ima*rk)
        evec(2)=evec(2)*(-ima*rk)
        evec(3)=evec(3)*(-ima*rk)
        hvec(1)=hvec(1)*(+ima*rk)
        hvec(2)=hvec(2)*(+ima*rk)
        hvec(3)=hvec(3)*(+ima*rk)

        call prin2('ecurl=*',ecurl,2*3)
        call prin2('hvec \times (+ima rk)=*',hvec,2*3)
        call prin2('hcurl=*',hcurl,2*3)
        call prin2('evec \times (-ima rk)=*',evec,2*3)


        return
        end
c
c
c
c
c
        subroutine emplanelw_check_curl(rk,target,itype)
        implicit real *8 (a-h,o-z)
        dimension source(3),target(3),target_new(3)
        complex *16 rk,ima,cjvec(3)
        complex *16 evec(3),hvec(3)
        complex *16 evec1(3),hvec1(3)
        complex *16 evec2(3),hvec2(3)
        complex *16 ederx(3),hderx(3)
        complex *16 edery(3),hdery(3)
        complex *16 ederz(3),hderz(3)
        complex *16 ecurl(3),hcurl(3)
        data ima/(0.0d0,1.0d0)/

        h=1d-6
        ediv=0
        hdiv=0

        target_new(1)=target(1)+h
        target_new(2)=target(2)
        target_new(3)=target(3)        
        call emplanelw(rk,target_new,itype,evec1,hvec1)
        target_new(1)=target(1)-h
        target_new(2)=target(2)
        target_new(3)=target(3)        
        call emplanelw(rk,target_new,itype,evec2,hvec2)
        ederx(1)=(evec1(1)-evec2(1))/(2*h)
        ederx(2)=(evec1(2)-evec2(2))/(2*h)
        ederx(3)=(evec1(3)-evec2(3))/(2*h)
        hderx(1)=(hvec1(1)-hvec2(1))/(2*h)
        hderx(2)=(hvec1(2)-hvec2(2))/(2*h)
        hderx(3)=(hvec1(3)-hvec2(3))/(2*h)

        target_new(1)=target(1)
        target_new(2)=target(2)+h
        target_new(3)=target(3)        
        call emplanelw(rk,target_new,itype,evec1,hvec1)
        target_new(1)=target(1)
        target_new(2)=target(2)-h
        target_new(3)=target(3)        
        call emplanelw(rk,target_new,itype,evec2,hvec2)
        edery(1)=(evec1(1)-evec2(1))/(2*h)
        edery(2)=(evec1(2)-evec2(2))/(2*h)
        edery(3)=(evec1(3)-evec2(3))/(2*h)
        hdery(1)=(hvec1(1)-hvec2(1))/(2*h)
        hdery(2)=(hvec1(2)-hvec2(2))/(2*h)
        hdery(3)=(hvec1(3)-hvec2(3))/(2*h)

        target_new(1)=target(1)
        target_new(2)=target(2)
        target_new(3)=target(3)+h        
        call emplanelw(rk,target_new,itype,evec1,hvec1)
        target_new(1)=target(1)
        target_new(2)=target(2)
        target_new(3)=target(3)-h        
        call emplanelw(rk,target_new,itype,evec2,hvec2)
        ederz(1)=(evec1(1)-evec2(1))/(2*h)
        ederz(2)=(evec1(2)-evec2(2))/(2*h)
        ederz(3)=(evec1(3)-evec2(3))/(2*h)
        hderz(1)=(hvec1(1)-hvec2(1))/(2*h)
        hderz(2)=(hvec1(2)-hvec2(2))/(2*h)
        hderz(3)=(hvec1(3)-hvec2(3))/(2*h)
c
c   see http://en.wikipedia.org/wiki/Curl_(mathematics)
c
        ecurl(1)=edery(3)-ederz(2)
        ecurl(2)=ederz(1)-ederx(3)
        ecurl(3)=ederx(2)-edery(1)

        hcurl(1)=hdery(3)-hderz(2)
        hcurl(2)=hderz(1)-hderx(3)
        hcurl(3)=hderx(2)-hdery(1)

        call emplanelw(rk,target,itype,evec,hvec)
        evec(1)=evec(1)*(-ima*rk)
        evec(2)=evec(2)*(-ima*rk)
        evec(3)=evec(3)*(-ima*rk)
        hvec(1)=hvec(1)*(+ima*rk)
        hvec(2)=hvec(2)*(+ima*rk)
        hvec(3)=hvec(3)*(+ima*rk)

        call prin2('ecurl=*',ecurl,2*3)
        call prin2('hvec \times (+ima rk)=*',hvec,2*3)
        call prin2('hcurl=*',hcurl,2*3)
        call prin2('evec \times (-ima rk)=*',evec,2*3)

        return
        end
c
c
c
c
c
        subroutine emplanearb_check_curl(rkvec,epol,target)
        implicit real *8 (a-h,o-z)
        dimension source(3),target(3),target_new(3)
        complex *16 rk,ima,cjvec(3),rkvec(3),epol(3)
        complex *16 evec(3),hvec(3)
        complex *16 evec1(3),hvec1(3)
        complex *16 evec2(3),hvec2(3)
        complex *16 ederx(3),hderx(3)
        complex *16 edery(3),hdery(3)
        complex *16 ederz(3),hderz(3)
        complex *16 ecurl(3),hcurl(3)
        data ima/(0.0d0,1.0d0)/

        h=1d-6
        ediv=0
        hdiv=0

        target_new(1)=target(1)+h
        target_new(2)=target(2)
        target_new(3)=target(3)        
        call emplanearb(rkvec,epol,target_new,evec1,hvec1)
        target_new(1)=target(1)-h
        target_new(2)=target(2)
        target_new(3)=target(3)        
        call emplanearb(rkvec,epol,target_new,evec2,hvec2)
        ederx(1)=(evec1(1)-evec2(1))/(2*h)
        ederx(2)=(evec1(2)-evec2(2))/(2*h)
        ederx(3)=(evec1(3)-evec2(3))/(2*h)
        hderx(1)=(hvec1(1)-hvec2(1))/(2*h)
        hderx(2)=(hvec1(2)-hvec2(2))/(2*h)
        hderx(3)=(hvec1(3)-hvec2(3))/(2*h)

        target_new(1)=target(1)
        target_new(2)=target(2)+h
        target_new(3)=target(3)        
        call emplanearb(rkvec,epol,target_new,evec1,hvec1)
        target_new(1)=target(1)
        target_new(2)=target(2)-h
        target_new(3)=target(3)        
        call emplanearb(rkvec,epol,target_new,evec2,hvec2)
        edery(1)=(evec1(1)-evec2(1))/(2*h)
        edery(2)=(evec1(2)-evec2(2))/(2*h)
        edery(3)=(evec1(3)-evec2(3))/(2*h)
        hdery(1)=(hvec1(1)-hvec2(1))/(2*h)
        hdery(2)=(hvec1(2)-hvec2(2))/(2*h)
        hdery(3)=(hvec1(3)-hvec2(3))/(2*h)

        target_new(1)=target(1)
        target_new(2)=target(2)
        target_new(3)=target(3)+h        
        call emplanearb(rkvec,epol,target_new,evec1,hvec1)
        target_new(1)=target(1)
        target_new(2)=target(2)
        target_new(3)=target(3)-h        
        call emplanearb(rkvec,epol,target_new,evec2,hvec2)
        ederz(1)=(evec1(1)-evec2(1))/(2*h)
        ederz(2)=(evec1(2)-evec2(2))/(2*h)
        ederz(3)=(evec1(3)-evec2(3))/(2*h)
        hderz(1)=(hvec1(1)-hvec2(1))/(2*h)
        hderz(2)=(hvec1(2)-hvec2(2))/(2*h)
        hderz(3)=(hvec1(3)-hvec2(3))/(2*h)
c
c   see http://en.wikipedia.org/wiki/Curl_(mathematics)
c
        ecurl(1)=edery(3)-ederz(2)
        ecurl(2)=ederz(1)-ederx(3)
        ecurl(3)=ederx(2)-edery(1)

        hcurl(1)=hdery(3)-hderz(2)
        hcurl(2)=hderz(1)-hderx(3)
        hcurl(3)=hderx(2)-hdery(1)

        call emplanearb(rkvec,epol,target,evec,hvec)
        rk=sqrt(rkvec(1)*rkvec(1)+rkvec(2)*rkvec(2)+rkvec(3)*rkvec(3))        
        evec(1)=evec(1)*(-ima*rk)
        evec(2)=evec(2)*(-ima*rk)
        evec(3)=evec(3)*(-ima*rk)
        hvec(1)=hvec(1)*(+ima*rk)
        hvec(2)=hvec(2)*(+ima*rk)
        hvec(3)=hvec(3)*(+ima*rk)

        call prin2('ecurl=*',ecurl,2*3)
        call prin2('hvec \times (+ima rk)=*',hvec,2*3)
        call prin2('hcurl=*',hcurl,2*3)
        call prin2('evec \times (-ima rk)=*',evec,2*3)

        return
        end
c
c
c
c
c
