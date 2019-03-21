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
c    $Date: 2011-05-15 12:07:04 -0400 (Sun, 15 May 2011) $
c    $Revision: 1973 $
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging routines and the beginning of
c       the code for the evaluation of dyadic Green's functions for the
c       Maxwell's equations in R^3
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c       This file contains a set of subroutines for the handling of
c       electromagnetic field dipoles. It contains 18 subroutines that
c       are user-callable. Following is a brief description of these
c       subroutines.
c
c     dipole3e - evaluate E and H fields at the location xyz due
c           to the monochromatic electric dipole cjvec 
c           located at the origin
c
c     dipole3m - evaluate E and H fields at the location xyz due
c           to the monochromatic magnetic dipole cmvec 
c           located at the origin
c
c     dipole3et - evaluate E and H fields at the location xyz due
c           to the monochromatic electric dipole cjvec 
c           located at the arbitrary source location
c
c     dipole3mt - evaluate E and H fields at the location xyz due
c           to the monochromatic magnetic dipole cmvec 
c           located at the arbitrary source location
c
c     dipole3em - evaluate E and H fields at the location xyz due to the
c           monochromatic electric dipole cjvec and the monochromatic
c           magnetic dipole cmvec located at the origin
c
c     dipole3emt - evaluate E and H fields at the location xyz due to the
c           monochromatic electric dipole cjvec and the monochromatic
c           magnetic dipole cmvec located at the arbitrary source location
c
c     dipole3epot - evaluate electric vector and scalar potentials at
c          the location xyz due to the monochromatic electric dipole cjvec
c          located at the origin
c
c     dipole3mpot - evaluate magnetic vector and scalar potentials at
c          the location xyz due to the monochromatic magnetic dipole cmvec
c          located at the origin
c
c     dipole3etpot - evaluate electric vector and scalar potentials at
c          the location xyz due to the monochromatic electric dipole cjvec
c          located at the arbitrary source location
c
c     dipole3mtpot - evaluate magnetic vector and scalar potentials at
c          the location xyz due to the monochromatic magnetic dipole cmvec
c          located at the arbitrary source location
c
c     green3e - dyadic electric Green's function for EM-field
c
c     green3ez - dyadic electric Green's function for EM-field, 
c                  extract the singular part
c
c     green3ez0 - dyadic electric Green's function for EM-field, 
c                   singular part only
c
c     green3m - dyadic magnetic Green's function for EM-field
c
c     green3mz - dyadic magnetic Green's function for EM-field,
c                  extract the singular part
c
c     green3mz0 - dyadic magnetic Green's function for EM-field,
c                   the singular part only
c
c     dipole3efar - evaluate the electric and magnetic far field
c          signatures in the direction dir due to the monochromatic electric
c          dipole cjvec located at the at the location xyz
c
c     dipole3mfar - evaluate the electric and magnetic far field
c          signatures in the direction dir due to the monochromatic magnetic
c          dipole cmvec located at the at the location xyz
c
c     dipole3epotfar - evaluate the vector and scalar potential far field
c          signatures in the direction dir due to the monochromatic electric
c          dipole cjvec located at the at the location xyz
c
c     dipole3mpotfar - evaluate the vector and scalar potential far field
c          signatures in the direction dir due to the monochromatic magnetic
c          dipole cmvec located at the at the location xyz
c
c     empoynting - evaluate the complex Poynting vector of EM-field
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     charge3pot - evaluates scalar potential at the location xyz due to 
c          the charge located at the origin
c
c     charge3tpot - evaluates scalar potential at the location xyz due to 
c          the charge located at the arbitrary source location
c
c     charge3far - evaluate the far field signature of the scalar field
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ... extract the singularity, use green3ez and green3mz 
c
c     dipole3ez - evaluate E and H fields at the location xyz due
c           to the monochromatic electric dipole cjvec 
c           located at the origin
c
c     dipole3mz - evaluate E and H fields at the location xyz due
c           to the monochromatic magnetic dipole cmvec 
c           located at the origin
c
c     dipole3ezt - evaluate E and H fields at the location xyz due
c           to the monochromatic electric dipole cjvec 
c           located at the arbitrary source location
c
c     dipole3mzt - evaluate E and H fields at the location xyz due
c           to the monochromatic magnetic dipole cmvec 
c           located at the arbitrary source location
c
c     dipole3emz - evaluate E and H fields at the location xyz due to the
c           monochromatic electric dipole cjvec and the monochromatic
c           magnetic dipole cmvec located at the origin
c
c     dipole3emtz - evaluate E and H fields at the location xyz due to the
c           monochromatic electric dipole cjvec and the monochromatic
c           magnetic dipole cmvec located at the arbitrary source location
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ... singular part only, use green3ez0 and green3mz0
c
c     dipole3ez0 - evaluate E and H fields at the location xyz due
c           to the monochromatic electric dipole cjvec 
c           located at the origin
c
c     dipole3mz0 - evaluate E and H fields at the location xyz due
c           to the monochromatic magnetic dipole cmvec 
c           located at the origin
c
c     dipole3ezt0 - evaluate E and H fields at the location xyz due
c           to the monochromatic electric dipole cjvec 
c           located at the arbitrary source location
c
c     dipole3mzt0 - evaluate E and H fields at the location xyz due
c           to the monochromatic magnetic dipole cmvec 
c           located at the arbitrary source location
c
c     dipole3emz0 - evaluate E and H fields at the location xyz due to the
c           monochromatic electric dipole cjvec and the monochromatic
c           magnetic dipole cmvec located at the origin
c
c     dipole3emtz0 - evaluate E and H fields at the location xyz due to the
c           monochromatic electric dipole cjvec and the monochromatic
c           magnetic dipole cmvec located at the arbitrary source location
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine dipole3e(rk,xyz,cjvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic electric dipole cjvec located at the origin
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3)
        complex *16 cjvec(3),evec(3),hvec(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        call green3e(rk,xyz,cjvec,evec)
        evec(1)=evec(1)*(ima*rk)
        evec(2)=evec(2)*(ima*rk)
        evec(3)=evec(3)*(ima*rk)
c       
        call green3m(rk,xyz,cjvec,hvec)
        hvec(1)=hvec(1)
        hvec(2)=hvec(2)
        hvec(3)=hvec(3)
c
        return
        end
c
c
c
c
c
        subroutine dipole3m(rk,xyz,cmvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic magnetic dipole cmvec located at the origin
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3)
        complex *16 cmvec(3),evec(3),hvec(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        call green3m(rk,xyz,cmvec,evec)
        evec(1)=evec(1)
        evec(2)=evec(2)
        evec(3)=evec(3)
c
        call green3e(rk,xyz,cmvec,hvec)
        hvec(1)=hvec(1)*(-ima*rk)
        hvec(2)=hvec(2)*(-ima*rk)
        hvec(3)=hvec(3)*(-ima*rk)
c       
        return
        end
c
c
c
c
c
        subroutine dipole3et(rk,source,target,cjvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic electric dipole cjvec located at the arbitrary 
c       source location
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       source (real *8 ) - the source point in R^3
c       target (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3),source(3),target(3)
        complex *16 cjvec(3),evec(3),hvec(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        call green3e(rk,xyz,cjvec,evec)
        evec(1)=evec(1)*(ima*rk)
        evec(2)=evec(2)*(ima*rk)
        evec(3)=evec(3)*(ima*rk)
c       
        call green3m(rk,xyz,cjvec,hvec)
        hvec(1)=hvec(1)
        hvec(2)=hvec(2)
        hvec(3)=hvec(3)
c
        return
        end
c
c
c
c
c
        subroutine dipole3mt(rk,source,target,cmvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic magnetic dipole cmvec located at the
c       arbitrary source location
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       source (real *8 ) - the source point in R^3
c       target (real *8 ) - the target point in R^3
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3),source(3),target(3)
        complex *16 cmvec(3),evec(3),hvec(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        call green3m(rk,xyz,cmvec,evec)
        evec(1)=evec(1)
        evec(2)=evec(2)
        evec(3)=evec(3)
c
        call green3e(rk,xyz,cmvec,hvec)
        hvec(1)=hvec(1)*(-ima*rk)
        hvec(2)=hvec(2)*(-ima*rk)
        hvec(3)=hvec(3)*(-ima*rk)
c       
        return
        end
c
c
c
c
c
        subroutine dipole3em(rk,xyz,cjvec,cmvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic electric dipole cjvec and the monochromatic
c       magnetic dipole cmvec located at the origin
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
        complex *16 rk
        dimension xyz(3)
        complex *16 cjvec(3),cmvec(3)
        complex *16 evec(3),hvec(3)
        complex *16 evec1(3),hvec1(3)
        complex *16 evec2(3),hvec2(3)
c
        call dipole3e(rk,xyz,cjvec,evec1,hvec1)
        call dipole3m(rk,xyz,cmvec,evec2,hvec2)
c
        evec(1)=evec1(1)+evec2(1)
        evec(2)=evec1(2)+evec2(2)
        evec(3)=evec1(3)+evec2(3)
c
        hvec(1)=hvec1(1)+hvec2(1)
        hvec(2)=hvec1(2)+hvec2(2)
        hvec(3)=hvec1(3)+hvec2(3)
c
        return
        end
c
c
c
c
c
        subroutine dipole3emt(rk,source,target,cjvec,cmvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic electric dipole cjvec and the monochromatic
c       magnetic dipole cmvec located at the arbitray location source
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       source (real *8 ) - the source point in R^3
c       target (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
        complex *16 rk
        dimension xyz(3),source(3),target(3)
        complex *16 cjvec(3),cmvec(3)
        complex *16 evec(3),hvec(3)
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        call dipole3em(rk,xyz,cjvec,cmvec,evec,hvec)
c
        return
        end
c
c
c
c
c
        subroutine dipole3ez(rk,xyz,cjvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic electric dipole cjvec located at the origin
c
c       extract the singularity
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3)
        complex *16 cjvec(3),evec(3),hvec(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        call green3ez(rk,xyz,cjvec,evec)
        evec(1)=evec(1)*(ima*rk)
        evec(2)=evec(2)*(ima*rk)
        evec(3)=evec(3)*(ima*rk)
c       
        call green3mz(rk,xyz,cjvec,hvec)
        hvec(1)=hvec(1)
        hvec(2)=hvec(2)
        hvec(3)=hvec(3)
c
        return
        end
c
c
c
c
c
        subroutine dipole3mz(rk,xyz,cmvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic magnetic dipole cmvec located at the origin
c
c       extract the singularity
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3)
        complex *16 cmvec(3),evec(3),hvec(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        call green3mz(rk,xyz,cmvec,evec)
        evec(1)=evec(1)
        evec(2)=evec(2)
        evec(3)=evec(3)
c
        call green3ez(rk,xyz,cmvec,hvec)
        hvec(1)=hvec(1)*(-ima*rk)
        hvec(2)=hvec(2)*(-ima*rk)
        hvec(3)=hvec(3)*(-ima*rk)
c       
        return
        end
c
c
c
c
c
        subroutine dipole3ezt(rk,source,target,cjvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic electric dipole cjvec located at the arbitrary 
c       source location
c
c       extract the singularity
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       source (real *8 ) - the source point in R^3
c       target (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3),source(3),target(3)
        complex *16 cjvec(3),evec(3),hvec(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        call green3ez(rk,xyz,cjvec,evec)
        evec(1)=evec(1)*(ima*rk)
        evec(2)=evec(2)*(ima*rk)
        evec(3)=evec(3)*(ima*rk)
c       
        call green3mz(rk,xyz,cjvec,hvec)
        hvec(1)=hvec(1)
        hvec(2)=hvec(2)
        hvec(3)=hvec(3)
c
        return
        end
c
c
c
c
c
        subroutine dipole3mzt(rk,source,target,cmvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic magnetic dipole cmvec located at the
c       arbitrary source location
c
c       extract the singularity
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       source (real *8 ) - the source point in R^3
c       target (real *8 ) - the target point in R^3
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3),source(3),target(3)
        complex *16 cmvec(3),evec(3),hvec(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        call green3mz(rk,xyz,cmvec,evec)
        evec(1)=evec(1)
        evec(2)=evec(2)
        evec(3)=evec(3)
c
        call green3ez(rk,xyz,cmvec,hvec)
        hvec(1)=hvec(1)*(-ima*rk)
        hvec(2)=hvec(2)*(-ima*rk)
        hvec(3)=hvec(3)*(-ima*rk)
c       
        return
        end
c
c
c
c
c
        subroutine dipole3emz(rk,xyz,cjvec,cmvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic electric dipole cjvec and the monochromatic
c       magnetic dipole cmvec located at the origin
c
c       extract the singularity
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
        complex *16 rk
        dimension xyz(3)
        complex *16 cjvec(3),cmvec(3)
        complex *16 evec(3),hvec(3)
        complex *16 evec1(3),hvec1(3)
        complex *16 evec2(3),hvec2(3)
c
        call dipole3ez(rk,xyz,cjvec,evec1,hvec1)
        call dipole3mz(rk,xyz,cmvec,evec2,hvec2)
c
        evec(1)=evec1(1)+evec2(1)
        evec(2)=evec1(2)+evec2(2)
        evec(3)=evec1(3)+evec2(3)
c
        hvec(1)=hvec1(1)+hvec2(1)
        hvec(2)=hvec1(2)+hvec2(2)
        hvec(3)=hvec1(3)+hvec2(3)
c
        return
        end
c
c
c
c
c
        subroutine dipole3emzt(rk,source,target,cjvec,cmvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic electric dipole cjvec and the monochromatic
c       magnetic dipole cmvec located at the arbitray location source
c
c       extract the singularity
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       source (real *8 ) - the source point in R^3
c       target (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
        complex *16 rk
        dimension xyz(3),source(3),target(3)
        complex *16 cjvec(3),cmvec(3)
        complex *16 evec(3),hvec(3)
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        call dipole3emz(rk,xyz,cjvec,cmvec,evec,hvec)
c
        return
        end
c
c
c
c
c
        subroutine dipole3ez0(rk,xyz,cjvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic electric dipole cjvec located at the origin
c
c       the singular part
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3)
        complex *16 cjvec(3),evec(3),hvec(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        call green3ez0(rk,xyz,cjvec,evec)
        evec(1)=evec(1)*(ima*rk)
        evec(2)=evec(2)*(ima*rk)
        evec(3)=evec(3)*(ima*rk)
c       
        call green3mz0(rk,xyz,cjvec,hvec)
        hvec(1)=hvec(1)
        hvec(2)=hvec(2)
        hvec(3)=hvec(3)
c
        return
        end
c
c
c
c
c
        subroutine dipole3mz0(rk,xyz,cmvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic magnetic dipole cmvec located at the origin
c
c       the singular part
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3)
        complex *16 cmvec(3),evec(3),hvec(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        call green3mz0(rk,xyz,cmvec,evec)
        evec(1)=evec(1)
        evec(2)=evec(2)
        evec(3)=evec(3)
c
        call green3ez0(rk,xyz,cmvec,hvec)
        hvec(1)=hvec(1)*(-ima*rk)
        hvec(2)=hvec(2)*(-ima*rk)
        hvec(3)=hvec(3)*(-ima*rk)
c       
        return
        end
c
c
c
c
c
        subroutine dipole3ezt0(rk,source,target,cjvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic electric dipole cjvec located at the arbitrary 
c       source location
c
c       the singular part
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       source (real *8 ) - the source point in R^3
c       target (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3),source(3),target(3)
        complex *16 cjvec(3),evec(3),hvec(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        call green3ez0(rk,xyz,cjvec,evec)
        evec(1)=evec(1)*(ima*rk)
        evec(2)=evec(2)*(ima*rk)
        evec(3)=evec(3)*(ima*rk)
c       
        call green3mz0(rk,xyz,cjvec,hvec)
        hvec(1)=hvec(1)
        hvec(2)=hvec(2)
        hvec(3)=hvec(3)
c
        return
        end
c
c
c
c
c
        subroutine dipole3mzt0(rk,source,target,cmvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic magnetic dipole cmvec located at the
c       arbitrary source location
c
c       the singular part
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       source (real *8 ) - the source point in R^3
c       target (real *8 ) - the target point in R^3
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
c
        dimension xyz(3),source(3),target(3)
        complex *16 cmvec(3),evec(3),hvec(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        call green3mz0(rk,xyz,cmvec,evec)
        evec(1)=evec(1)
        evec(2)=evec(2)
        evec(3)=evec(3)
c
        call green3ez0(rk,xyz,cmvec,hvec)
        hvec(1)=hvec(1)*(-ima*rk)
        hvec(2)=hvec(2)*(-ima*rk)
        hvec(3)=hvec(3)*(-ima*rk)
c       
        return
        end
c
c
c
c
c
        subroutine dipole3emz0(rk,xyz,cjvec,cmvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic electric dipole cjvec and the monochromatic
c       magnetic dipole cmvec located at the origin
c
c       the singular part
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
        complex *16 rk
        dimension xyz(3)
        complex *16 cjvec(3),cmvec(3)
        complex *16 evec(3),hvec(3)
        complex *16 evec1(3),hvec1(3)
        complex *16 evec2(3),hvec2(3)
c
        call dipole3ez0(rk,xyz,cjvec,evec1,hvec1)
        call dipole3mz0(rk,xyz,cmvec,evec2,hvec2)
c
        evec(1)=evec1(1)+evec2(1)
        evec(2)=evec1(2)+evec2(2)
        evec(3)=evec1(3)+evec2(3)
c
        hvec(1)=hvec1(1)+hvec2(1)
        hvec(2)=hvec1(2)+hvec2(2)
        hvec(3)=hvec1(3)+hvec2(3)
c
        return
        end
c
c
c
c
c
        subroutine dipole3emzt0(rk,source,target,cjvec,cmvec,evec,hvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates E and H fields at the location xyz due
c       to the monochromatic electric dipole cjvec and the monochromatic
c       magnetic dipole cmvec located at the arbitray location source
c
c       the singular part
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       source (real *8 ) - the source point in R^3
c       target (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       evec (complex*16) - the electric field at the target
c       hvec (complex*16) - the magnetic field at the target
c
        complex *16 rk
        dimension xyz(3),source(3),target(3)
        complex *16 cjvec(3),cmvec(3)
        complex *16 evec(3),hvec(3)
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        call dipole3emz0(rk,xyz,cjvec,cmvec,evec,hvec)
c
        return
        end
c
c
c
c
c
        subroutine dipole3epot(rk,xyz,cjvec,vpot,spot)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates electric vector and scalar potentials
c       at the location xyz due to the monochromatic electric dipole
c       cjvec located at the origin
c
c       (free space, Lorenz gauge)
c
c                       [     ]   exp(I*rk*R)
c       A   = vpot    = [ I_3 ]   -----------
c                       [     ]        R
c
c       where I_3 is a 3x3 identity matrix
c
c       phi = spot    = 1/(I*rk) div A
c
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       vpot (complex*16) - the vector potential at the target
c       spot (complex*16) - the scalar potential at the target
c
c
        dimension xyz(3)
        complex *16 cjvec(3),vpot(3),spot,fout,rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx**2+dy**2+dz**2)
c
c
        fout=exp(ima*rk*cd)/cd
c
        vpot(1)=+fout*cjvec(1)
        vpot(2)=+fout*cjvec(2)
        vpot(3)=+fout*cjvec(3)
c
c
        fout=(dx*cjvec(1)+dy*cjvec(2)+dz*cjvec(3))/(ima*rk)
c
        spot=-fout*(1-ima*rk*cd)*exp(ima*rk*cd)/cd**3
c
        return
        end
c
c
c
c
c
        subroutine dipole3mpot(rk,xyz,cmvec,vpot,spot)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates magnetic vector and scalar potentials
c       at the location xyz due to the monochromatic magnetic dipole
c       cmvec located at the origin
c
c       (free space, Lorenz gauge)
c
c                       [     ]   exp(I*rk*R)
c       A   = vpot    = [ I_3 ]   -----------
c                       [     ]        R
c
c       where I_3 is a 3x3 identity matrix
c
c       phi = spot    = 1/(I*rk) div A
c
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       vpot (complex*16) - the vector potential at the target
c       spot (complex*16) - the scalar potential at the target
c
c
        dimension xyz(3)
        complex *16 cmvec(3),vpot(3),spot,fout,rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx**2+dy**2+dz**2)
c
c
        fout=exp(ima*rk*cd)/cd
c
        vpot(1)=-fout*cmvec(1)
        vpot(2)=-fout*cmvec(2)
        vpot(3)=-fout*cmvec(3)
c
c
        fout=(dx*cmvec(1)+dy*cmvec(2)+dz*cmvec(3))/(ima*rk)
c
        spot=+fout*(1-ima*rk*cd)*exp(ima*rk*cd)/cd**3
c
        return
        end
c
c
c
c
c
        subroutine dipole3etpot(rk,source,target,cjvec,vpot,spot)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates electric vector and scalar potentials
c       at the location target due to the monochromatic electric dipole
c       cjvec located at the location source
c
c       (free space, Lorenz gauge)
c
c                       [     ]   exp(I*rk*R)
c       A   = vpot    = [ I_3 ]   -----------
c                       [     ]        R
c
c       where I_3 is a 3x3 identity matrix
c
c       phi = spot    = 1/(I*rk) div A
c
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       source (real *8 ) - the source point in R^3
c       target (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       vpot (complex*16) - the vector potential at the target
c       spot (complex*16) - the scalar potential at the target
c
c
        dimension xyz(3),source(3),target(3)
        complex *16 cjvec(3),vpot(3),spot,fout,rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx**2+dy**2+dz**2)
c
c
        fout=exp(ima*rk*cd)/cd
c
        vpot(1)=+fout*cjvec(1)
        vpot(2)=+fout*cjvec(2)
        vpot(3)=+fout*cjvec(3)
c
c
        fout=(dx*cjvec(1)+dy*cjvec(2)+dz*cjvec(3))/(ima*rk)
c
        spot=-fout*(1-ima*rk*cd)*exp(ima*rk*cd)/cd**3
c
        return
        end
c
c
c
c
c
        subroutine dipole3mtpot(rk,source,target,cmvec,vpot,spot)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates magnetic vector and scalar potentials
c       at the location target due to the monochromatic magnetic dipole
c       cmvec located at the location source
c
c       (free space, Lorenz gauge)
c
c                       [     ]   exp(I*rk*R)
c       A   = vpot    = [ I_3 ]   -----------
c                       [     ]        R
c
c       where I_3 is a 3x3 identity matrix
c
c       phi = spot    = 1/(I*rk) div A
c
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       source (real *8 ) - the source point in R^3
c       target (real *8 ) - the target point in R^3
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       vpot (complex*16) - the vector potential at the target
c       spot (complex*16) - the scalar potential at the target
c
c
        dimension xyz(3),source(3),target(3)
        complex *16 cmvec(3),vpot(3),spot,fout,rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx**2+dy**2+dz**2)
c
c
        fout=exp(ima*rk*cd)/cd
c
        vpot(1)=-fout*cmvec(1)
        vpot(2)=-fout*cmvec(2)
        vpot(3)=-fout*cmvec(3)
c
c
        fout=(dx*cmvec(1)+dy*cmvec(2)+dz*cmvec(3))/(ima*rk)
c
        spot=+fout*(1-ima*rk*cd)*exp(ima*rk*cd)/cd**3
c
        return
        end
c
c
c
c
c
        subroutine green3e(rk,xyz,cjvec,fvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the electric dyadic Green's function
c       at the location xyz due to the monochromatic electric dipole
c       cjvec located at the origin
c
c       ... dyadic Green's function for the electric field 
c       (free space, Lorenz gauge)
c
c                 [        \grad  \grad  ]   exp(I*rk*R)
c       green3e = [ I_3 +  ------------  ]   -----------
c                 [            rk^2      ]        R
c
c       where I_3 is a 3x3 identity matrix
c
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       fvec (complex*16) - the value of Green's function at the target
c
c
        dimension xyz(3)
        complex *16 cjvec(3),fvec(3),fout,qmat(3,3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx**2+dy**2+dz**2)
c
        fout=exp(ima*rk*cd)/cd
c
        fvec(1)=fout*cjvec(1)
        fvec(2)=fout*cjvec(2)
        fvec(3)=fout*cjvec(3)
c
        qmat(1,1)=(2*dx**2-dy**2-dz**2)*(1-ima*rk*cd)
        qmat(2,2)=(2*dy**2-dz**2-dx**2)*(1-ima*rk*cd)
        qmat(3,3)=(2*dz**2-dx**2-dy**2)*(1-ima*rk*cd)
c
        qmat(1,1)=qmat(1,1)+(-rk**2*dx**2*cd**2)
        qmat(2,2)=qmat(2,2)+(-rk**2*dy**2*cd**2)
        qmat(3,3)=qmat(3,3)+(-rk**2*dz**2*cd**2)
c
        qmat(1,2)=(3-rk**2*cd**2-3*ima*rk*cd)*(dx*dy)
        qmat(2,3)=(3-rk**2*cd**2-3*ima*rk*cd)*(dy*dz)
        qmat(3,1)=(3-rk**2*cd**2-3*ima*rk*cd)*(dz*dx)
c
        qmat(2,1)=qmat(1,2)
        qmat(3,2)=qmat(2,3)
        qmat(1,3)=qmat(3,1)
c
        fout=exp(ima*rk*cd)/cd**5/rk**2
c
        fvec(1)=fvec(1) + fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) + fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) + fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c
        return
        end
c
c
c
c
c
        subroutine green3ez(rk,xyz,cjvec,fvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the electric dyadic Green's function
c       at the location xyz due to the monochromatic electric dipole
c       cjvec located at the origin
c
c       ... dyadic Green's function for the electric field 
c       (free space, Lorenz gauge), extract the singular part
c
c                  [     ]  [ exp(I*rk*R)     1 ]
c       green3ez = [ I_3 ]  [ ------------ -  - ]
c                  [     ]  [      R          R ]
c
c                  [  \grad  \grad  ]  [ exp(I*rk*R)     1   rk**2   ]
c                + [  ------------  ]  [ ------------ -  - +  ---  R ]
c                  [      rk^2      ]  [      R          R     2     ]
c
c       where I_3 is a 3x3 identity matrix
c
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       fvec (complex*16) - the value of Green's function at the target
c
c
        dimension xyz(3)
        complex *16 cjvec(3),fvec(3),fout,qmat(3,3),rk,ima
        complex *16 ca,cb,cc
c
        data ima/(0.0d0,1.0d0)/
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx**2+dy**2+dz**2)
c
        fout=exp(ima*rk*cd)/cd
c
        fout=fout-1/cd
ccc        fout=fout+rk**2*cd/2
c
c                      2               3              4                 5
c 1     I k           k     2       I k    3         k     4         I k    5
c---- + --- a~ - 1/2 ---- a~  - 1/6 ---- a~  + 1/24 ---- a~  + 1/120 ---- a~
c a~    a~            a~             a~              a~               a~
c
c           6
c     + O(a~ )
c
        if( abs(rk*cd) .lt. 1d-4 ) then
ccc        call prin2('in green3ez, before taylor, fout=*',fout,2)
        fout = rk * (+ima - ima*(rk*cd)**2/6 + 
     $     (rk*cd)**3/24 + ima*(rk*cd)**4/120 )
        fout = fout - rk**2*cd/2
ccc        call prin2('in green3ez, after taylor, fout=*',fout,2)
        endif
c
        fvec(1)=fout*cjvec(1)
        fvec(2)=fout*cjvec(2)
        fvec(3)=fout*cjvec(3)
c

        if( abs(rk*cd) .lt. 1d-4 ) then
c
c      4  2   -1            5  2          3             6  2        4
c 1/8 k  x  a~   + (1/15 I k  x  - 1/3 I k ) + (- 1/48 k  x  + 1/8 k ) a~ +
c
c             5            7  2    2             6           8  2    3
c    (1/30 I k  - 1/210 I k  x ) a~  + (- 1/144 k  + 1/1152 k  x ) a~  +
c
c                7             9  2    4       5
c    (- 1/840 I k  + 1/7560 I k  x ) a~  + O(a~ )
c

c
c      4       -1           5             6                   7       2
c 1/8 k  x y a~   + 1/15 I k  x y - 1/48 k  x y a~ - 1/210 I k  x y a~  +
c
c            8       3             9       4       5
c    1/1152 k  x y a~  + 1/7560 I k  x y a~  + O(a~ )
c
c
        ca = rk**3 *( -ima/3 + (rk*cd)/8 + ima*(rk*cd)**2/30 
     $     - (rk*cd)**3/144 - ima*(rk*cd)**4/840 )

        cb = rk**4/cd/8 +
     $     rk**5 * (ima/15 - (rk*cd)/48 - ima*(rk*cd)**2/210 
     $     + (rk*cd)**3/1152 + ima*(rk*cd)**4/7560 )

        qmat(1,1) = ca+cb*(dx*dx)
        qmat(2,2) = ca+cb*(dy*dy)
        qmat(3,3) = ca+cb*(dz*dz)

        qmat(1,2) = cb*(dx*dy)
        qmat(2,3) = cb*(dy*dz)
        qmat(3,1) = cb*(dz*dx)

        qmat(2,1) = qmat(1,2)
        qmat(3,2) = qmat(2,3)
        qmat(1,3) = qmat(3,1)
c
        fout = 1/rk**2
c
        fvec(1)=fvec(1) + fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) + fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) + fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c

        else

c
        qmat(1,1)=(2*dx**2-dy**2-dz**2)*(1-ima*rk*cd)
        qmat(2,2)=(2*dy**2-dz**2-dx**2)*(1-ima*rk*cd)
        qmat(3,3)=(2*dz**2-dx**2-dy**2)*(1-ima*rk*cd)
c
        qmat(1,1)=qmat(1,1)+(-rk**2*dx**2*cd**2)
        qmat(2,2)=qmat(2,2)+(-rk**2*dy**2*cd**2)
        qmat(3,3)=qmat(3,3)+(-rk**2*dz**2*cd**2)
c
        qmat(1,2)=(3-rk**2*cd**2-3*ima*rk*cd)*(dx*dy)
        qmat(2,3)=(3-rk**2*cd**2-3*ima*rk*cd)*(dy*dz)
        qmat(3,1)=(3-rk**2*cd**2-3*ima*rk*cd)*(dz*dx)
c
        qmat(2,1)=qmat(1,2)
        qmat(3,2)=qmat(2,3)
        qmat(1,3)=qmat(3,1)
c
        fout=exp(ima*rk*cd)/cd**5/rk**2
c
        fvec(1)=fvec(1) + fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) + fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) + fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c
c
        qmat(1,1)=(2*dx**2-dy**2-dz**2)
        qmat(2,2)=(2*dy**2-dz**2-dx**2)
        qmat(3,3)=(2*dz**2-dx**2-dy**2)
c
        qmat(1,2)=3*(dx*dy)
        qmat(2,3)=3*(dy*dz)
        qmat(3,1)=3*(dz*dx)
c
        qmat(2,1)=qmat(1,2)
        qmat(3,2)=qmat(2,3)
        qmat(1,3)=qmat(3,1)
c
        fout=1/cd**5/rk**2
c
        fvec(1)=fvec(1) - fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) - fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) - fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c
c       
        qmat(1,1)=dy**2+dz**2
        qmat(2,2)=dz**2+dx**2
        qmat(3,3)=dx**2+dy**2
c
        qmat(1,2)=-(dx*dy)
        qmat(2,3)=-(dy*dz)
        qmat(3,1)=-(dz*dx)
c
        qmat(2,1)=qmat(1,2)
        qmat(3,2)=qmat(2,3)
        qmat(1,3)=qmat(3,1)
c
        fout=1/cd**3/2
c
        fvec(1)=fvec(1) + fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) + fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) + fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c
        endif
c        
        return
        end
c
c
c
c
c
        subroutine green3ez0(rk,xyz,cjvec,fvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the electric dyadic Green's function
c       at the location xyz due to the monochromatic electric dipole
c       cjvec located at the origin
c
c       ... dyadic Green's function for the electric field 
c       (free space, Lorenz gauge), the singular part
c
c                   [     ]  [   1  ]
c       green3ez0 = [ I_3 ]  [ - -  ]
c                   [     ]  [   R  ]
c
c                   [ \grad  \grad  ]  [   1   rk**2   ]
c                 + [ ------------  ]  [ - - +  ---  R ]
c                   [     rk^2      ]  [   R     2     ]
c
c       where I_3 is a 3x3 identity matrix
c
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       fvec (complex*16) - the value of Green's function at the target
c
c
        dimension xyz(3)
        complex *16 cjvec(3),fvec(3),fout,qmat(3,3),rk,ima
        complex *16 ca,cb,cc
c
        data ima/(0.0d0,1.0d0)/
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx**2+dy**2+dz**2)
c
        fout=0
c
        fout=fout-1/cd
ccc        fout=fout+rk**2*cd/2
c
c                      2               3              4                 5
c 1     I k           k     2       I k    3         k     4         I k    5
c---- + --- a~ - 1/2 ---- a~  - 1/6 ---- a~  + 1/24 ---- a~  + 1/120 ---- a~
c a~    a~            a~             a~              a~               a~
c
c           6
c     + O(a~ )
c
        if( abs(rk*cd) .lt. 1d-4 ) then
ccc        call prin2('in green3ez, before taylor, fout=*',fout,2)        
        fout = rk * (+ima - ima*(rk*cd)**2/6 + 
     $     (rk*cd)**3/24 + ima*(rk*cd)**4/120 )
        fout = fout - rk**2*cd/2
ccc        call prin2('in green3ez, after taylor, fout=*',fout,2)
        endif
c
        fvec(1)=fout*cjvec(1)
        fvec(2)=fout*cjvec(2)
        fvec(3)=fout*cjvec(3)
c
        if( abs(rk*cd) .lt. 1d-4 ) then
c
c      4  2   -1            5  2          3             6  2        4
c 1/8 k  x  a~   + (1/15 I k  x  - 1/3 I k ) + (- 1/48 k  x  + 1/8 k ) a~ +
c
c             5            7  2    2             6           8  2    3
c    (1/30 I k  - 1/210 I k  x ) a~  + (- 1/144 k  + 1/1152 k  x ) a~  +
c
c                7             9  2    4       5
c    (- 1/840 I k  + 1/7560 I k  x ) a~  + O(a~ )
c

c
c      4       -1           5             6                   7       2
c 1/8 k  x y a~   + 1/15 I k  x y - 1/48 k  x y a~ - 1/210 I k  x y a~  +
c
c            8       3             9       4       5
c    1/1152 k  x y a~  + 1/7560 I k  x y a~  + O(a~ )
c
c
        ca = rk**3 *( -ima/3 + (rk*cd)/8 + ima*(rk*cd)**2/30 
     $     - (rk*cd)**3/144 - ima*(rk*cd)**4/840 )

        cb = rk**4/cd/8 +
     $     rk**5 * (ima/15 - (rk*cd)/48 - ima*(rk*cd)**2/210 
     $     + (rk*cd)**3/1152 + ima*(rk*cd)**4/7560 )

        qmat(1,1) = ca+cb*(dx*dx)
        qmat(2,2) = ca+cb*(dy*dy)
        qmat(3,3) = ca+cb*(dz*dz)

        qmat(1,2) = cb*(dx*dy)
        qmat(2,3) = cb*(dy*dz)
        qmat(3,1) = cb*(dz*dx)

        qmat(2,1) = qmat(1,2)
        qmat(3,2) = qmat(2,3)
        qmat(1,3) = qmat(3,1)
c
        fout = 1/rk**2
c
        fvec(1)=fvec(1) + fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) + fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) + fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c

        else

c
        qmat(1,1)=(2*dx**2-dy**2-dz**2)
        qmat(2,2)=(2*dy**2-dz**2-dx**2)
        qmat(3,3)=(2*dz**2-dx**2-dy**2)
c
        qmat(1,2)=3*(dx*dy)
        qmat(2,3)=3*(dy*dz)
        qmat(3,1)=3*(dz*dx)
c
        qmat(2,1)=qmat(1,2)
        qmat(3,2)=qmat(2,3)
        qmat(1,3)=qmat(3,1)
c
        fout=1/cd**5/rk**2
c
        fvec(1)=fvec(1) - fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) - fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) - fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c
c       
        qmat(1,1)=dy**2+dz**2
        qmat(2,2)=dz**2+dx**2
        qmat(3,3)=dx**2+dy**2
c
        qmat(1,2)=-(dx*dy)
        qmat(2,3)=-(dy*dz)
        qmat(3,1)=-(dz*dx)
c
        qmat(2,1)=qmat(1,2)
        qmat(3,2)=qmat(2,3)
        qmat(1,3)=qmat(3,1)
c
        fout=1/cd**3/2
c
        fvec(1)=fvec(1) + fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) + fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) + fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c
        endif
c        
        return
        end
c
c
c
c
c
        subroutine green3ez1(rk,xyz,cjvec,fvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the electric dyadic Green's function
c       at the location xyz due to the monochromatic electric dipole
c       cjvec located at the origin
c
c       ... dyadic Green's function for the electric field 
c       (free space, Lorenz gauge), extract the singular part
c
c                   [        \grad  \grad  ]  [ exp(I*rk*R)     1   rk**2   ]
c       green3ez1 = [ I_3 +  ------------  ]  [ ------------ -  - +  ---  R ]
c                   [            rk^2      ]  [      R          R     2     ]
c
c                  [        \grad  \grad  ]   rk**4     
c                - [        ------------  ]   ----- R^3 
c                  [            rk^2      ]     24       
c
c       where I_3 is a 3x3 identity matrix
c
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       fvec (complex*16) - the value of Green's function at the target
c
c
        dimension xyz(3)
        complex *16 cjvec(3),fvec(3),fout,qmat(3,3),rk,ima
        complex *16 ca,cb,cc
c
        data ima/(0.0d0,1.0d0)/
c
        fvec(1)=0
        fvec(2)=0
        fvec(3)=0
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx**2+dy**2+dz**2)
c
        fout=exp(ima*rk*cd)/cd
c
        fout=fout-1/cd
        fout=fout+rk**2*cd/2
c
c                      2               3              4                 5
c 1     I k           k     2       I k    3         k     4         I k    5
c---- + --- a~ - 1/2 ---- a~  - 1/6 ---- a~  + 1/24 ---- a~  + 1/120 ---- a~
c a~    a~            a~             a~              a~               a~
c
c           6
c     + O(a~ )
c        
        if( abs(rk*cd) .lt. 1d-4 ) then
c        call prin2('in green3ez, before taylor, fout=*',fout,2)
c        call prin2('abs(rk \times cd)=*',abs(rk*cd),1)
        fout = rk * (+ima - ima*(rk*cd)**2/6 + 
     $     (rk*cd)**3/24 + ima*(rk*cd)**4/120 )
c        call prin2('in green3ez, after taylor, fout=*',fout,2)
        endif
c
        fvec(1)=fout*cjvec(1)
        fvec(2)=fout*cjvec(2)
        fvec(3)=fout*cjvec(3)
c
c
        if( abs(rk*cd) .lt. 1d-4 ) then
c
c         5  2          3          6  2               5            7  2    2
c(1/15 I k  x  - 1/3 I k ) - 1/48 k  x  a~ + (1/30 I k  - 1/210 I k  x ) a~  +
c
c              6           8  2    3               7             9  2    4
c    (- 1/144 k  + 1/1152 k  x ) a~  + (- 1/840 I k  + 1/7560 I k  x ) a~  +
c
c        5
c    O(a~ )
c
c        5             6                   7       2           8       3
c1/15 I k  x y - 1/48 k  x y a~ - 1/210 I k  x y a~  + 1/1152 k  x y a~  +
c
c              9       4       5
c    1/7560 I k  x y a~  + O(a~ )
c
c
        ca = rk**3 *( -ima/3 + ima*(rk*cd)**2/30 
     $     - (rk*cd)**3/144 - ima*(rk*cd)**4/840 )

        cb = 
     $     rk**5 * (ima/15 - (rk*cd)/48 - ima*(rk*cd)**2/210 
     $     + (rk*cd)**3/1152 + ima*(rk*cd)**4/7560 )

        qmat(1,1) = ca+cb*(dx*dx)
        qmat(2,2) = ca+cb*(dy*dy)
        qmat(3,3) = ca+cb*(dz*dz)

        qmat(1,2) = cb*(dx*dy)
        qmat(2,3) = cb*(dy*dz)
        qmat(3,1) = cb*(dz*dx)

        qmat(2,1) = qmat(1,2)
        qmat(3,2) = qmat(2,3)
        qmat(1,3) = qmat(3,1)
c
        fout = 1/rk**2
c
        fvec(1)=fvec(1) + fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) + fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) + fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c

        else

c
        qmat(1,1)=(2*dx**2-dy**2-dz**2)*(1-ima*rk*cd)
        qmat(2,2)=(2*dy**2-dz**2-dx**2)*(1-ima*rk*cd)
        qmat(3,3)=(2*dz**2-dx**2-dy**2)*(1-ima*rk*cd)
c
        qmat(1,1)=qmat(1,1)+(-rk**2*dx**2*cd**2)
        qmat(2,2)=qmat(2,2)+(-rk**2*dy**2*cd**2)
        qmat(3,3)=qmat(3,3)+(-rk**2*dz**2*cd**2)
c
        qmat(1,2)=(3-rk**2*cd**2-3*ima*rk*cd)*(dx*dy)
        qmat(2,3)=(3-rk**2*cd**2-3*ima*rk*cd)*(dy*dz)
        qmat(3,1)=(3-rk**2*cd**2-3*ima*rk*cd)*(dz*dx)
c
        qmat(2,1)=qmat(1,2)
        qmat(3,2)=qmat(2,3)
        qmat(1,3)=qmat(3,1)
c
        fout=exp(ima*rk*cd)/cd**5/rk**2
c
        fvec(1)=fvec(1) + fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) + fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) + fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c
c
        qmat(1,1)=(2*dx**2-dy**2-dz**2)
        qmat(2,2)=(2*dy**2-dz**2-dx**2)
        qmat(3,3)=(2*dz**2-dx**2-dy**2)
c
        qmat(1,2)=3*(dx*dy)
        qmat(2,3)=3*(dy*dz)
        qmat(3,1)=3*(dz*dx)
c
        qmat(2,1)=qmat(1,2)
        qmat(3,2)=qmat(2,3)
        qmat(1,3)=qmat(3,1)
c
        fout=1/cd**5/rk**2
c
        fvec(1)=fvec(1) - fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) - fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) - fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c
c       
        qmat(1,1)=dy**2+dz**2
        qmat(2,2)=dz**2+dx**2
        qmat(3,3)=dx**2+dy**2
c
        qmat(1,2)=-(dx*dy)
        qmat(2,3)=-(dy*dz)
        qmat(3,1)=-(dz*dx)
c
        qmat(2,1)=qmat(1,2)
        qmat(3,2)=qmat(2,3)
        qmat(1,3)=qmat(3,1)
c
        fout=1/cd**3/2
c
        fvec(1)=fvec(1) + fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) + fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) + fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c
c      
        qmat(1,1)=3*(2*dx**2+dy**2+dz**2)
        qmat(2,2)=3*(2*dy**2+dz**2+dx**2)
        qmat(3,3)=3*(2*dz**2+dx**2+dy**2)
c
        qmat(1,2)=3*(dx*dy)
        qmat(2,3)=3*(dy*dz)
        qmat(3,1)=3*(dz*dx)
c
        qmat(2,1)=qmat(1,2)
        qmat(3,2)=qmat(2,3)
        qmat(1,3)=qmat(3,1)
c
        fout=(1/cd)*(rk**2/24)
c
        fvec(1)=fvec(1) - fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) - fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) - fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c
        endif
c
        return
        end
c
c
c
c
c
        subroutine green3ez2(rk,xyz,cjvec,fvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the electric dyadic Green's function
c       at the location xyz due to the monochromatic electric dipole
c       cjvec located at the origin
c
c       ... dyadic Green's function for the electric field 
c       (free space, Lorenz gauge), extract the singular part
c
c                   [        \grad  \grad  ]  [ exp(I*rk*R)     1   rk**2   ]
c       green3ez2 = [ I_3 +  ------------  ]  [ ------------ -  - +  ---  R ]
c                   [            rk^2      ]  [      R          R     2     ]
c
c       where I_3 is a 3x3 identity matrix
c
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       fvec (complex*16) - the value of Green's function at the target
c
c
        dimension xyz(3)
        complex *16 cjvec(3),fvec(3),fout,qmat(3,3),rk,ima
        complex *16 ca,cb,cc
c
        data ima/(0.0d0,1.0d0)/
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx**2+dy**2+dz**2)
c
        fout=exp(ima*rk*cd)/cd
c
        fout=fout-1/cd
        fout=fout+rk**2*cd/2
c
c                      2               3              4                 5
c 1     I k           k     2       I k    3         k     4         I k    5
c---- + --- a~ - 1/2 ---- a~  - 1/6 ---- a~  + 1/24 ---- a~  + 1/120 ---- a~
c a~    a~            a~             a~              a~               a~
c
c           6
c     + O(a~ )
c
        if( abs(rk*cd) .lt. 1d-4 ) then
ccc        call prin2('in green3ez, before taylor, fout=*',fout,2)
        fout = rk * (+ima - ima*(rk*cd)**2/6 + 
     $     (rk*cd)**3/24 + ima*(rk*cd)**4/120 )
ccc        call prin2('in green3ez, after taylor, fout=*',fout,2)
        endif
c
        fvec(1)=fout*cjvec(1)
        fvec(2)=fout*cjvec(2)
        fvec(3)=fout*cjvec(3)
c

        if( abs(rk*cd) .lt. 1d-4 ) then
c
c      4  2   -1            5  2          3             6  2        4
c 1/8 k  x  a~   + (1/15 I k  x  - 1/3 I k ) + (- 1/48 k  x  + 1/8 k ) a~ +
c
c             5            7  2    2             6           8  2    3
c    (1/30 I k  - 1/210 I k  x ) a~  + (- 1/144 k  + 1/1152 k  x ) a~  +
c
c                7             9  2    4       5
c    (- 1/840 I k  + 1/7560 I k  x ) a~  + O(a~ )
c

c
c      4       -1           5             6                   7       2
c 1/8 k  x y a~   + 1/15 I k  x y - 1/48 k  x y a~ - 1/210 I k  x y a~  +
c
c            8       3             9       4       5
c    1/1152 k  x y a~  + 1/7560 I k  x y a~  + O(a~ )
c
c
        ca = rk**3 *( -ima/3 + (rk*cd)/8 + ima*(rk*cd)**2/30 
     $     - (rk*cd)**3/144 - ima*(rk*cd)**4/840 )

        cb = rk**4/cd/8 +
     $     rk**5 * (ima/15 - (rk*cd)/48 - ima*(rk*cd)**2/210 
     $     + (rk*cd)**3/1152 + ima*(rk*cd)**4/7560 )

        qmat(1,1) = ca+cb*(dx*dx)
        qmat(2,2) = ca+cb*(dy*dy)
        qmat(3,3) = ca+cb*(dz*dz)

        qmat(1,2) = cb*(dx*dy)
        qmat(2,3) = cb*(dy*dz)
        qmat(3,1) = cb*(dz*dx)

        qmat(2,1) = qmat(1,2)
        qmat(3,2) = qmat(2,3)
        qmat(1,3) = qmat(3,1)
c
        fout = 1/rk**2
c
        fvec(1)=fvec(1) + fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) + fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) + fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c

        else

c
        qmat(1,1)=(2*dx**2-dy**2-dz**2)*(1-ima*rk*cd)
        qmat(2,2)=(2*dy**2-dz**2-dx**2)*(1-ima*rk*cd)
        qmat(3,3)=(2*dz**2-dx**2-dy**2)*(1-ima*rk*cd)
c
        qmat(1,1)=qmat(1,1)+(-rk**2*dx**2*cd**2)
        qmat(2,2)=qmat(2,2)+(-rk**2*dy**2*cd**2)
        qmat(3,3)=qmat(3,3)+(-rk**2*dz**2*cd**2)
c
        qmat(1,2)=(3-rk**2*cd**2-3*ima*rk*cd)*(dx*dy)
        qmat(2,3)=(3-rk**2*cd**2-3*ima*rk*cd)*(dy*dz)
        qmat(3,1)=(3-rk**2*cd**2-3*ima*rk*cd)*(dz*dx)
c
        qmat(2,1)=qmat(1,2)
        qmat(3,2)=qmat(2,3)
        qmat(1,3)=qmat(3,1)
c
        fout=exp(ima*rk*cd)/cd**5/rk**2
c
        fvec(1)=fvec(1) + fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) + fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) + fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c
c
        qmat(1,1)=(2*dx**2-dy**2-dz**2)
        qmat(2,2)=(2*dy**2-dz**2-dx**2)
        qmat(3,3)=(2*dz**2-dx**2-dy**2)
c
        qmat(1,2)=3*(dx*dy)
        qmat(2,3)=3*(dy*dz)
        qmat(3,1)=3*(dz*dx)
c
        qmat(2,1)=qmat(1,2)
        qmat(3,2)=qmat(2,3)
        qmat(1,3)=qmat(3,1)
c
        fout=1/cd**5/rk**2
c
        fvec(1)=fvec(1) - fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) - fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) - fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c
c       
        qmat(1,1)=dy**2+dz**2
        qmat(2,2)=dz**2+dx**2
        qmat(3,3)=dx**2+dy**2
c
        qmat(1,2)=-(dx*dy)
        qmat(2,3)=-(dy*dz)
        qmat(3,1)=-(dz*dx)
c
        qmat(2,1)=qmat(1,2)
        qmat(3,2)=qmat(2,3)
        qmat(1,3)=qmat(3,1)
c
        fout=1/cd**3/2
c
        fvec(1)=fvec(1) + fout*
     $     (qmat(1,1)*cjvec(1)+qmat(1,2)*cjvec(2)+qmat(1,3)*cjvec(3))
c
        fvec(2)=fvec(2) + fout*
     $     (qmat(2,1)*cjvec(1)+qmat(2,2)*cjvec(2)+qmat(2,3)*cjvec(3))
c
        fvec(3)=fvec(3) + fout*
     $     (qmat(3,1)*cjvec(1)+qmat(3,2)*cjvec(2)+qmat(3,3)*cjvec(3))
c
        endif
c        
        return
        end
c
c
c
c
c
        subroutine green3m(rk,xyz,cjvec,fvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the magnetic dyadic Green's function
c       at the location xyz due to the monochromatic electric dipole
c       cjvec located at the origin
c
c       ... dyadic Green's function for the magnetic field 
c       (free space, Lorenz gauge)
c
c       green3m =  curl(green3e) = \grad \times green3e
c
c
c                 [        exp(I*rk*R) ]
c       green3m = [ \grad  ----------- ] \times
c                 [             R      ]
c      
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       fvec (complex*16) - the value of Green's function at the target
c
c
c
        dimension xyz(3)
        complex *16 cjvec(3),fvec(3),fout,fder(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx**2+dy**2+dz**2)
c
        fout=(1-ima*rk*cd)*exp(ima*rk*cd)/cd**3
c
        fder(1)=dx*fout
        fder(2)=dy*fout
        fder(3)=dz*fout
c
        call green3cvecp(fder,cjvec,fvec)
c
        fvec(1)=-fvec(1)
        fvec(2)=-fvec(2)
        fvec(3)=-fvec(3)
c
        return
        end
c
c
c
c
c
        subroutine green3mz(rk,xyz,cjvec,fvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the magnetic dyadic Green's function
c       at the location xyz due to the monochromatic electric dipole
c       cjvec located at the origin
c
c       ... dyadic Green's function for the magnetic field 
c       (free space, Lorenz gauge), extract the singular part
c
c       green3mz =  curl(green3ez) = \grad \times green3ez
c
c
c                  [       [ exp(I*rk*R)   1   rk^2   ]]
c       green3mz = [ \grad [ ----------- - - + ---- R ]] \times
c                  [       [      R        R     2    ]]
c      
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       fvec (complex*16) - the value of Green's function at the target
c
c
c
        dimension xyz(3)
        complex *16 cjvec(3),fvec(3),fout,fder(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx**2+dy**2+dz**2)
c
        fout=(1-ima*rk*cd)*exp(ima*rk*cd)/cd**3
c
c
c                    2
c     x             k  x            3          4    2 1/2             5   2
c- -------- - 1/2 -------- - 1/3 I k  x + 1/8 k  (a~ )    x + 1/30 I k  a~  x
c     2 3/2          2 1/2
c  (a~ )          (a~ )
c
c              6    2 3/2              7   4         6
c     - 1/144 k  (a~ )    x - 1/840 I k  a~  x + O(a~  x)
c
c
        fout=fout-1/cd**3
        fout=fout-rk**2/cd/2
c
        if( abs(rk*cd) .lt. 1d-4 ) then
ccc        call prin2('in green3mz, before taylor, fout=*',fout,2)
        fout = rk**3 * ( +ima/3 - (rk*cd)/8 - ima*(rk*cd)**2/30 
     $     + (rk*cd)**3/144 + ima*(rk*cd)**4/840 ) 
ccc        call prin2('in green3mz, after taylor, fout=*',fout,2)
        endif
c
        fder(1)=dx*fout
        fder(2)=dy*fout
        fder(3)=dz*fout
c
        call green3cvecp(fder,cjvec,fvec)
c
        fvec(1)=-fvec(1)
        fvec(2)=-fvec(2)
        fvec(3)=-fvec(3)
c
        return
        end
c
c
c
c
c
        subroutine green3mz0(rk,xyz,cjvec,fvec)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the magnetic dyadic Green's function
c       at the location xyz due to the monochromatic electric dipole
c       cjvec located at the origin
c
c       ... dyadic Green's function for the magnetic field 
c       (free space, Lorenz gauge), the singular part
c
c       green3mz0 =  curl(green3ez0) = \grad \times green3ez0
c
c
c                   [       [   1   rk^2   ]]
c       green3mz0 = [ \grad [ - - + ---- R ]] \times
c                   [       [   R     2    ]]
c      
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       fvec (complex*16) - the value of Green's function at the target
c
c
c
        dimension xyz(3)
        complex *16 cjvec(3),fvec(3),fout,fder(3),rk,ima
c
        data ima/(0.0d0,1.0d0)/

c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx**2+dy**2+dz**2)
c
        fout=0
c
c
c                    2
c     x             k  x            3          4    2 1/2             5   2
c- -------- - 1/2 -------- - 1/3 I k  x + 1/8 k  (a~ )    x + 1/30 I k  a~  x
c     2 3/2          2 1/2
c  (a~ )          (a~ )
c
c              6    2 3/2              7   4         6
c     - 1/144 k  (a~ )    x - 1/840 I k  a~  x + O(a~  x)
c
c
        fout=fout-1/cd**3
        fout=fout-rk**2/cd/2
c
        if( abs(rk*cd) .lt. 1d-3 ) then
ccc        call prin2('in green3mz, before taylor, fout=*',fout,2)
        fout = rk**3 * ( +ima/3 - (rk*cd)/8 - ima*(rk*cd)**2/30 
     $     + (rk*cd)**3/144 + ima*(rk*cd)**4/840 ) 
ccc        call prin2('in green3mz, after taylor, fout=*',fout,2)
        endif
c
        fder(1)=dx*fout
        fder(2)=dy*fout
        fder(3)=dz*fout
c
        call green3cvecp(fder,cjvec,fvec)
c
        fvec(1)=-fvec(1)
        fvec(2)=-fvec(2)
        fvec(3)=-fvec(3)
c
        return
        end
c
c
c
c
c
        subroutine dipole3efar(rk,xyz,cjvec,dir,efar,hfar)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the electric and magnetic far field
c       signatures in the direction dir due to the monochromatic electric
c       dipole cjvec located at the at the location xyz
c
c       (free space, Lorenz gauge)
c
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of the electric dipole in R^3
c       dir (real *8 ) - the far field direction vector in R^3 
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       efar (complex*16) - the far field signature of the E-field
c       hfar (complex*16) - the far field signature of the H-field
c
c
        dimension xyz(3),dir(3),far(3)
        complex *16 cjvec(3),efar(3),hfar(3),fout,rk,ima,cfvec(3)
        complex *16 cfar
c
        data ima/(0.0d0,1.0d0)/
c
c       ... first, find the point on the unit sphere
c        
        d=sqrt(dir(1)**2+dir(2)**2+dir(3)**2)
        far(1)=dir(1)/d
        far(2)=dir(2)/d
        far(3)=dir(3)/d
c
c       ... find the far-field signature of the unit charge located at
c       the location xyz
c
        d=xyz(1)*far(1)+xyz(2)*far(2)+xyz(3)*far(3)
        cfar=exp(-ima*rk*d)
c       
c
c       ... and the far-field signature of the E-field
c
        efar(1)=cjvec(1)
        efar(2)=cjvec(2)
        efar(3)=cjvec(3)
c       
        fout=far(1)*cjvec(1)+far(2)*cjvec(2)+far(3)*cjvec(3)
c
        efar(1)=efar(1)-far(1)*fout
        efar(2)=efar(2)-far(2)*fout
        efar(3)=efar(3)-far(3)*fout
c
        efar(1)=(ima*rk)*cfar*efar(1)
        efar(2)=(ima*rk)*cfar*efar(2)
        efar(3)=(ima*rk)*cfar*efar(3)
c
c       ... and the far-field signature of the H-field
c
        cfvec(1)=far(1)
        cfvec(2)=far(2)
        cfvec(3)=far(3)
        call green3cvecp(cfvec,cjvec,hfar)     
        hfar(1)=(ima*rk)*cfar*hfar(1)
        hfar(2)=(ima*rk)*cfar*hfar(2)
        hfar(3)=(ima*rk)*cfar*hfar(3)        
c
        return
        end
c
c
c
c
c
        subroutine dipole3epotfar(rk,xyz,cjvec,dir,vfar,sfar)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the vector and scalar potential far field
c       signatures in the direction dir due to the monochromatic electric
c       dipole cjvec located at the at the location xyz
c
c       (free space, Lorenz gauge)
c
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of the electric dipole in R^3
c       dir (real *8 ) - the far field direction vector in R^3 
c       cjvec (complex *16) - the strength of the electric dipole   
c
c          Output parameters:
c
c       vfar (complex*16) - the far field signature of the vector potential
c       sfar (complex*16) - the far field signature of the scalar potential
c
c
        dimension xyz(3),dir(3),far(3)
        complex *16 cjvec(3),vfar(3),sfar,fout,rk,ima,cfvec(3)
        complex *16 cfar
c
        data ima/(0.0d0,1.0d0)/
c
c       ... first, find the point on the unit sphere
c        
        d=sqrt(dir(1)**2+dir(2)**2+dir(3)**2)
        far(1)=dir(1)/d
        far(2)=dir(2)/d
        far(3)=dir(3)/d
c
c       ... find the far-field signature of the unit charge located at
c       the location xyz
c
        d=xyz(1)*far(1)+xyz(2)*far(2)+xyz(3)*far(3)
        cfar=exp(-ima*rk*d)
c       
c
c       ... and the far-field signature of the vector potential
c
        vfar(1)=cjvec(1)
        vfar(2)=cjvec(2)
        vfar(3)=cjvec(3)
c       
        vfar(1)=cfar*vfar(1)
        vfar(2)=cfar*vfar(2)
        vfar(3)=cfar*vfar(3)
c
c       ... and the far-field signature of the scalar potential
c
        fout=cjvec(1)*far(1)+cjvec(2)*far(2)+cjvec(3)*far(3)
c       
        sfar=cfar*fout
c
        return
        end
c
c
c
c
c
        subroutine dipole3mfar(rk,xyz,cmvec,dir,efar,hfar)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the electric and magnetic far field
c       signatures in the direction dir due to the monochromatic magnetic
c       dipole cmvec located at the at the location xyz
c
c       (free space, Lorenz gauge)
c
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of the electric dipole in R^3
c       dir (real *8 ) - the far field direction vector in R^3 
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       efar (complex*16) - the far field signature of the E-field
c       hfar (complex*16) - the far field signature of the H-field
c
c
        dimension xyz(3),dir(3),far(3)
        complex *16 cmvec(3),efar(3),hfar(3),fout,rk,ima,cfvec(3)
        complex *16 cfar
c
        data ima/(0.0d0,1.0d0)/
c
c       ... first, find the point on the unit sphere
c        
        d=sqrt(dir(1)**2+dir(2)**2+dir(3)**2)
        far(1)=dir(1)/d
        far(2)=dir(2)/d
        far(3)=dir(3)/d
c
c       ... find the far-field signature of the unit charge located at
c       the location xyz
c
        d=xyz(1)*far(1)+xyz(2)*far(2)+xyz(3)*far(3)
        cfar=exp(-ima*rk*d)
c
c
c       ... and the far-field signature of the E-field
c
        cfvec(1)=far(1)
        cfvec(2)=far(2)
        cfvec(3)=far(3)
        call green3cvecp(cfvec,cmvec,efar)     
        efar(1)=(ima*rk)*cfar*efar(1)
        efar(2)=(ima*rk)*cfar*efar(2)
        efar(3)=(ima*rk)*cfar*efar(3)        
c
c
c       ... and the far-field signature of the H-field
c
        hfar(1)=cmvec(1)
        hfar(2)=cmvec(2)
        hfar(3)=cmvec(3)
c       
        fout=far(1)*cmvec(1)+far(2)*cmvec(2)+far(3)*cmvec(3)
c
        hfar(1)=hfar(1)-far(1)*fout
        hfar(2)=hfar(2)-far(2)*fout
        hfar(3)=hfar(3)-far(3)*fout
c
        hfar(1)=(-ima*rk)*cfar*hfar(1)
        hfar(2)=(-ima*rk)*cfar*hfar(2)
        hfar(3)=(-ima*rk)*cfar*hfar(3)
c
        return
        end
c
c
c
c
c
        subroutine dipole3mpotfar(rk,xyz,cmvec,dir,vfar,sfar)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the vector and scalar potential far field
c       signatures in the direction dir due to the monochromatic magnetic
c       dipole cmvec located at the at the location xyz
c
c       (free space, Lorenz gauge)
c
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of the electric dipole in R^3
c       dir (real *8 ) - the far field direction vector in R^3 
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       vfar (complex*16) - the far field signature of the vector potential
c       sfar (complex*16) - the far field signature of the scalar potential
c
c
        dimension xyz(3),dir(3),far(3)
        complex *16 cmvec(3),vfar(3),sfar,fout,rk,ima,cfvec(3)
        complex *16 cfar
c
        data ima/(0.0d0,1.0d0)/
c
c       ... first, find the point on the unit sphere
c        
        d=sqrt(dir(1)**2+dir(2)**2+dir(3)**2)
        far(1)=dir(1)/d
        far(2)=dir(2)/d
        far(3)=dir(3)/d
c
c       ... find the far-field signature of the unit charge located at
c       the location xyz
c
        d=xyz(1)*far(1)+xyz(2)*far(2)+xyz(3)*far(3)
        cfar=exp(-ima*rk*d)
c
c
c       ... and the far-field signature of the vector potential
c
        vfar(1)=cmvec(1)
        vfar(2)=cmvec(2)
        vfar(3)=cmvec(3)
c       
        vfar(1)=-cfar*vfar(1)
        vfar(2)=-cfar*vfar(2)
        vfar(3)=-cfar*vfar(3)
c
c
c       ... and the far-field signature of the scalar potential
c
        fout=cmvec(1)*far(1)+cmvec(2)*far(2)+cmvec(3)*far(3)
c       
        sfar=-cfar*fout
c
        return
        end
c
c
c
c
c
        subroutine dipole3emfar(rk,xyz,cjvec,cmvec,dir,efar,hfar)
        implicit real *8 (a-h,o-z)
c
c
c       This subroutine evaluates the electric and magnetic far field
c       signatures in the direction dir due to the monochromatic
c       electric dipole cjvec and the monochromatic magnetic dipole
c       cmvec located at the at the location xyz
c
c       (free space, Lorenz gauge)
c
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of the electric dipole in R^3
c       dir (real *8 ) - the far field direction vector in R^3 
c       cjvec (complex *16) - the strength of the electric dipole   
c       cmvec (complex *16) - the strength of the magnetic dipole   
c
c          Output parameters:
c
c       efar (complex*16) - the far field signature of the E-field
c       hfar (complex*16) - the far field signature of the H-field
c
c
        complex *16 rk
        dimension xyz(3),dir(3)
        complex *16 cjvec(3),cmvec(3)
        complex *16 efar(3),hfar(3)
        complex *16 efar1(3),hfar1(3)
        complex *16 efar2(3),hfar2(3)
c
        call dipole3efar(rk,xyz,cjvec,dir,efar1,hfar1)
        call dipole3mfar(rk,xyz,cmvec,dir,efar2,hfar2)
c
        efar(1)=efar1(1)+efar2(1)
        efar(2)=efar1(2)+efar2(2)
        efar(3)=efar1(3)+efar2(3)
c
        hfar(1)=hfar1(1)+hfar2(1)
        hfar(2)=hfar1(2)+hfar2(2)
        hfar(3)=hfar1(3)+hfar2(3)
c
        return
        end
c
c
c
c
c
        subroutine empoynting(evec,hvec,pvec)
        implicit real *8 (a-h,o-z)
        complex *16 evec(3),hvec(3),pvec(3)
c
c       this subroutine evaluates the complex Poynting vector of EM field
c
c       P = E x conjg(H)
c
c       WARNING: the time averaged complex Poynting vector is usually defined
c       as <P> = 0.5 Re(E x conjg(H)) (as in J. Jackson's "Classical
c       Electrodynamics") USER BEWARE!
c
c         Input parameters:
c
c       evec (complex *16) - E-field 
c       hvec (complex *16) - H-field 
c
c         Output parameters:
c
c       pvec (complex *16) - Poynting vector: E x conjg(H)
c
        pvec(1)=evec(2)*dconjg(hvec(3))-evec(3)*dconjg(hvec(2))
        pvec(2)=evec(3)*dconjg(hvec(1))-evec(1)*dconjg(hvec(3))
        pvec(3)=evec(1)*dconjg(hvec(2))-evec(2)*dconjg(hvec(1))
c
        return
        end
c
c
c
c
c
        subroutine empoynting2(evec,hvec,pvec)
        implicit real *8 (a-h,o-z)
        complex *16 evec(3),hvec(3),pvec(3)
c
c       This subroutine evaluates correctly defined
c       complex Poynting vector of EM field
c
c       P = 1/2 (E x conjg(H))
c       
c       The time averaged Poynting vector is equal to the real part of
c       complex Poynting vector: <P> = 0.5 Re(E x conjg(H))
c
c       (see J. Jackson's "Classical Electrodynamics") 
c
c         Input parameters:
c
c       evec (complex *16) - E-field 
c       hvec (complex *16) - H-field 
c
c         Output parameters:
c
c       pvec (complex *16) - Poynting vector: E x conjg(H)
c
        pvec(1)=evec(2)*dconjg(hvec(3))-evec(3)*dconjg(hvec(2))
        pvec(2)=evec(3)*dconjg(hvec(1))-evec(1)*dconjg(hvec(3))
        pvec(3)=evec(1)*dconjg(hvec(2))-evec(2)*dconjg(hvec(1))
c
        pvec(1)=pvec(1)*0.5d0
        pvec(2)=pvec(2)*0.5d0
        pvec(3)=pvec(3)*0.5d0
c
        return
        end
c
c
c
c
c
        subroutine green3cvecp(x,y,z)     
        implicit real *8 (a-h,o-z)
        complex *16 x(3),y(3),z(3)
c       
c       return the complex-valued vector product in R^3
c
c       z = [x,y] = x \times y
c       
        z(1)=x(2)*y(3)-y(2)*x(3)
        z(2)=x(3)*y(1)-y(3)*x(1)
        z(3)=x(1)*y(2)-y(1)*x(2)
        return
        end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine charge3pot(rk,xyz,charge,cpot)
        implicit real *8 (a-h,o-z)
c
        dimension xyz(3)
        complex *16 rk,charge,cpot,ima
        data ima/(0.0d0,1.0d0)/
c
c       This subroutine evaluates the scalar potential at the location
c       xyz due to the charge located at the origin
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the target point in R^3
c       charge (complex *16) - the strength of the charge
c
c          Output parameters:
c
c       cpot (complex*16) - the potential at the target
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c       
        cpot=exp(ima*rk*cd)/cd
        cpot=cpot*charge
c
        return
        end
c
c
c
c
c
        subroutine charge3tpot(rk,source,target,charge,cpot)
        implicit real *8 (a-h,o-z)
c
        dimension source(3),target(3)
        complex *16 rk,charge,cpot,ima
        data ima/(0.0d0,1.0d0)/
c
c       This subroutine evaluates the scalar potential at the location
c       target due to the charge located at the location source
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       source (real *8 ) - the source point in R^3
c       target (real *8 ) - the target point in R^3
c       charge (complex *16) - the strength of the charge
c
c          Output parameters:
c
c       cpot (complex*16) - the potential at the target
c
        dx=target(1)-source(1)
        dy=target(2)-source(2)
        dz=target(3)-source(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c       
        cpot=exp(ima*rk*cd)/cd
        cpot=cpot*charge
c
        return
        end
c
c
c
c
c
        subroutine charge3far(rk,xyz,charge,dir,cfar)
        implicit real *8 (a-h,o-z)
c
c       This subroutine evaluates the far field
c       signatures in the direction dir due to the monochromatic
c       charge located at the at the location xyz
c
c          Input parameters:
c
c       rk (complex *16)  - the frequency parameter
c       xyz (real *8 ) - the location of the electric dipole in R^3
c       dir (real *8 ) - the far field direction vector in R^3 
c       charge (complex *16) - the strength of the charge
c
c          Output parameters:
c
c       cfar (complex*16) - the far field signature of the potential
c
        dimension xyz(3),dir(3),far(3)
        complex *16 rk,charge,cfar,ima
        data ima/(0.0d0,1.0d0)/
c
c       ... first, find the point on the unit sphere
c        
        d=sqrt(dir(1)**2+dir(2)**2+dir(3)**2)
        far(1)=dir(1)/d
        far(2)=dir(2)/d
        far(3)=dir(3)/d
c
c       ... find the far-field signature of the unit charge located at
c       the location xyz
c
        d=xyz(1)*far(1)+xyz(2)*far(2)+xyz(3)*far(3)
        cfar=exp(-ima*rk*d)
c
        cfar=cfar*charge
c
        return
        end
