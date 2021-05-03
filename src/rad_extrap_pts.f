c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for setting the mask chrbdy to a value different 
c from ex at the outermost grid points that will then be used 
c for radial extrapolation.  
c
c The extrapolation technique goes as follows:
c 1. define a range in which we look for grid points  
c to use for extrapolation. Consider the outermost  
c points in this range. 
c 2. Take one of such points, denoted by p1.  
c This will be the vertex of a cube (the other vertices  
c are other grid points). Consider the radial direction  
c connecting the point to the origin. This direction  
c intersects a face of the cube at a point p2. This will be  
c the second point used for extrapolation. The angular coords  
c of p2 are the same as p1, the radial coord must be obtained  
c and it depends on the quadrant of the 3-dimensional grid in  
c which p1 is.
c 3. To obtain the value of the function at p2, we use bilinear  
c extrapolation from the 4 vertices of the face of the cube (if  
c these are not accessible from the current process, we discard  
c the point p1 and move on to the next one).
c----------------------------------------------------------------------

        subroutine nexttobdypoints_radextrap(
     &                  chrbdy,
     &                  numbdypoints,
     &                  bdy_extrap_order,
     &                  currentres_ratio_Lhighres_Llowres,
     &                  half_steps_from_bdy_ext,
     &                  half_steps_from_bdy_int,
     &                  x,y,z,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width)

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        real*8 chrbdy(Nx,Ny,Nz),chrbdy2(Nx,Ny,Nz)

        integer i,j,k,is,ie,js,je,ks,ke
        integer ipa,jpa,kpa
        integer ipb,jpb,kpb
        integer ipc,jpc,kpc
        integer ipd,jpd,kpd

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz

        real*8 x0,y0,z0,rho0,q

        real*8 dx,dy,dz
        real*8 xp1,yp1,zp1,rhop1,chip1,xip1
        real*8 maxxyzp1
        integer numbdypoints
        integer half_steps_from_bdy_ext
        integer half_steps_from_bdy_int
        real*8 currentres_ratio_Lhighres_Llowres
        integer bdy_extrap_order

        real*8 PI
        parameter (PI=3.141592653589793d0)

!----------------------------------------------------------------------

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        numbdypoints=0

        ! set index bounds for main loop
        is=2
        ie=Nx-1
        js=2
        je=Ny-1
        ks=2
        ke=Nz-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)

         do i=1,Nx
          do j=1,Ny
           do k=1,Nz
           xp1=x(i)
           yp1=y(j)
           zp1=z(k)
           rhop1=sqrt(xp1**2+yp1**2+zp1**2)
           chip1=(1/PI)*acos(xp1/rhop1)
           if (zp1.lt.0) then
              xip1=(1/(2*PI))*(atan2(zp1,yp1)+2*PI)
           else
              xip1=(1/(2*PI))*atan2(zp1,yp1)
           end if

!select only points that have radius smaller than (1-5*dx/2).
!For each resolution, points between rhobdy=1 and 1-dx/2 are excised, points between rhobdy=1 and 1-3*dx/2 use forward/backward stencils (so when we cannot expect convergence at these points because the stencils used are different for different resolutions), points between rhobdy=1 and 1-5*dx/2 use points that are set by forward/backward stencils, so we cannot expect convergence. If we want to check convergence at the grid points used for extrapolation at the boundary, we need to pick points that have rho<1-5*dx/2. However, we see that we get good convergence only if we use points with rho<1-9*dx/2. This should be understood in the future. The condition rhop1>currentres_ratio_Lhighres_Llowres ensures that the first suitable point for extrapolation is not too far from the AdS boundary.
            if ((rhop1.lt.(1-half_steps_from_bdy_ext*dx/2))
     &          .and.(rhop1.gt.
     &                (1-half_steps_from_bdy_int
     &                  *currentres_ratio_Lhighres_Llowres*dx/2))
     &          .and.(chr(i,j,k).ne.ex)
     &         ) then
              chrbdy(i,j,k) =ex-1.0d0
              chrbdy2(i,j,k)=ex-1.0d0
            else
              chrbdy(i,j,k) =ex
              chrbdy2(i,j,k)=ex
            end if


            if ((i.lt.is).or.(i.gt.ie).or.
     &          (j.lt.js).or.(j.gt.je).or.
     &          (k.lt.ks).or.(k.gt.ke)) then

               chrbdy(i,j,k) =ex
               chrbdy2(i,j,k)=ex
            end if

           end do
          end do
         end do


         do i=is,ie
          do j=js,je
           do k=ks,ke

           xp1=x(i)
           yp1=y(j)
           zp1=z(k)
           rhop1=sqrt(xp1**2+yp1**2+zp1**2)
           chip1=(1/PI)*acos(xp1/rhop1)
           if (zp1.lt.0) then
              xip1=(1/(2*PI))*(atan2(zp1,yp1)+2*PI)
           else
              xip1=(1/(2*PI))*atan2(zp1,yp1)
           end if


           if (chrbdy(i,j,k).ne.ex) then


!! eliminate troublesome points
!!removing points with abs(x)=abs(y)=abs(z): even if the grid function quasiset_tracell converges at these points, the extrapolated boundary quantity quasisettrace, obtained using these points, is not converging. It would be interesting to understand why. For now, we just eliminate these points from the set of points used for extrapolation
!            if (
!     &          (abs(abs(xp1)-abs(yp1)).lt.10.0d0**(-10)).and.
!     &          (abs(abs(xp1)-abs(zp1)).lt.10.0d0**(-10))
!     &         ) then
!                  chrbdy(i,j,k)=ex
!            end if




! eliminate troublesome points
!removing z=0 implies that we remove in particular the troublesome points 
!with chi=0,1 (which have y=z=0,x=1,-1) and points with xi=0,1 (which have z=0,y>0,any x). 
!We will fill and impose regularity at these points in post-processing in Mathematica.
            if (
     &          (abs(zp1).lt.10.0d0**(-10)) 
     &         ) then
             chrbdy(i,j,k)=ex
            end if






!If we use derivatives to define near boundary quantities, we will only define them at points between is and ie (js and je, ks and ke). Therefore, for extrapolation, we can only select near boundary points whose neighbors used for extrapolation in the direction of the bulk along the axes (i.e. the direction of extrapolation) are within that range
!We also need to make sure that those neighbours are not excised.
!We also define near boundary quantities at points where y0 and z0 are not both 0. So we need to make sure that we select points such that neighbouring points used for extrapolation don't have such values of y0,z0. Notice, we've already imposed that z(k) is not 0, we only need to impose the condition when extrapolation is along z, so the value of the z-coordinate of the second point used is different from the first one.
!The condition (chrbdy2(i+1,j,k).ne.ex) makes sure that (i,j,k) is the outmost point satisfying the conditions of the previous for-loop, which sets chrbdy2 as well as chrbdy. In other words, if there's an outer point w.r.t. (i,j,k) that satisfies those conditions, then we don't want to use (i,j,k) for extrapolation, but we will use that other point. 


            if (bdy_extrap_order.eq.1) then



              if ((abs(xp1).gt.abs(yp1)).and.
     &            (abs(xp1).gt.abs(zp1))) then !(i.e., |xp1|>|yp1|,|zp1|, so xp1 cannot be 0)
            
                if (xp1.gt.0) then !(i.e., we are in the upper part of 
                                    !the spatial grid, called "a" sector)

                  ipa=i-1
                  ipb=i-1
                  ipc=i-1
                  ipd=i-1

                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((i+1).le.ie) then 
                    if (chrbdy2(i+1,j,k).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(yp1).gt.abs(zp1)) then !(i.e., |yp1|>|zp1|, so yp1 cannot be 0)
                    if (yp1.gt.0) then !(i.e., either quadrant Ia or IVa)

                      jpa=j
                      jpb=j
                      jpc=j-1
                      jpd=j-1

                      if (zp1.gt.0) then !(i.e., quadrant Ia)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IVa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., yp1<=0, so Q.IIa or IIIa)

                      jpa=j
                      jpb=j
                      jpc=j+1
                      jpd=j+1

                      if (zp1.gt.0) then !(i.e., quadrant IIa)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IIIa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |zp1|>=|yp1|, it's possible that zp1=yp1=0)

                    if (zp1.gt.0) then !(i.e., either quadrant Ia or IIa)

                      kpa=k
                      kpb=k
                      kpc=k-1
                      kpd=k-1

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIa)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., zp1<=0, so Q.IIIa or IVa)

                      kpa=k
                      kpb=k
                      kpc=k+1
                      kpd=k+1

                      if (yp1.gt.0) then !(i.e., quadrant IVa)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIIa)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(yp1).gt.abs(zp1))

                else !(i.e., xp1<0, i.e., we are in the lower part of 
                     !the spatial grid, called "b" sector)

                  ipa=i+1
                  ipb=i+1
                  ipc=i+1
                  ipd=i+1


                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((i-1).ge.is) then 
                    if (chrbdy2(i-1,j,k).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(yp1).gt.abs(zp1)) then !(i.e., |yp1|>|zp1|, so yp1 cannot be 0)
                    if (yp1.gt.0) then !(i.e., either quadrant Ia or IVa)

                      jpa=j
                      jpb=j
                      jpc=j-1
                      jpd=j-1

                      if (zp1.gt.0) then !(i.e., quadrant Ia)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IVa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., yp1<=0, so Q.IIa or IIIa)

                      jpa=j
                      jpb=j
                      jpc=j+1
                      jpd=j+1

                      if (zp1.gt.0) then !(i.e., quadrant IIa)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IIIa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |zp1|>=|yp1|, it's possible that zp1=yp1=0)

                    if (zp1.gt.0) then !(i.e., either quadrant Ia or IIa)

                      kpa=k
                      kpb=k
                      kpc=k-1
                      kpd=k-1

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIa)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., zp1<=0, so Q.IIIa or IVa)

                      kpa=k
                      kpb=k
                      kpc=k+1
                      kpd=k+1

                      if (yp1.gt.0) then !(i.e., quadrant IVa)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIIa)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(yp1).gt.abs(zp1))

                end if !closes condition on (xp1.gt.0)











              else if ((abs(yp1).gt.abs(zp1)).and.
     &                 (abs(yp1).gt.abs(xp1))) then !(i.e., |yp1|>|zp1|,|xp1|, so yp1 cannot be 0)


                if (yp1.gt.0) then !(i.e., we are in either Q.Ia or Q.Ib or Q.IVa or Q.IVb)

                  jpa=j-1
                  jpb=j-1
                  jpc=j-1
                  jpd=j-1

                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((j+1).le.je) then 
                    if (chrbdy2(i,j+1,k).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(zp1).gt.abs(xp1)) then !(i.e., |zp1|>|xp1|, so zp1 cannot be 0)
                    if (zp1.gt.0) then !(i.e., either quadrant Ia or Ib)

                      kpa=k
                      kpb=k
                      kpc=k-1
                      kpd=k-1

                      if (xp1.gt.0) then !(i.e., quadrant Ia)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant Ib)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., zp1<=0, so Q.IVa or Q.IVb)

                      kpa=k
                      kpb=k
                      kpc=k+1
                      kpd=k+1

                      if (xp1.gt.0) then !(i.e., quadrant IVa)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IVb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |xp1|>=|zp1| - xp1=zp1=0 is not an issue)

                    if (xp1.gt.0) then !(i.e., either quadrant Ia or IVa)

                      ipa=i
                      ipb=i
                      ipc=i-1
                      ipd=i-1

                      if (zp1.gt.0) then !(i.e., quadrant Ia)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IVa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., xp1<=0, so Q.Ib or IVb)

                      ipa=i
                      ipb=i
                      ipc=i+1
                      ipd=i+1

                      if (zp1.gt.0) then !(i.e., quadrant Ib)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IVb)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(zp1).gt.abs(xp1))

                else !(i.e., yp1<0, so we are in either Q.IIa or Q.IIb or Q.IIIa or Q.IIIb)

                  jpa=j+1
                  jpb=j+1
                  jpc=j+1
                  jpd=j+1


                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((j-1).ge.js) then 
                    if (chrbdy2(i,j-1,k).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(zp1).gt.abs(xp1)) then !(i.e., |zp1|>|xp1|, so zp1 cannot be 0)
                    if (zp1.gt.0) then !(i.e., either quadrant IIa or IIb)

                      kpa=k
                      kpb=k
                      kpc=k-1
                      kpd=k-1

                      if (xp1.gt.0) then !(i.e., quadrant IIa)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IIb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., zp1<=0, so Q.IIIa or IIIb)

                      kpa=k
                      kpb=k
                      kpc=k+1
                      kpd=k+1

                      if (xp1.gt.0) then !(i.e., quadrant IIIa)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IIIb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |xp1|>=|zp1|, it's possible that xp1=zp1=0)

                    if (xp1.gt.0) then !(i.e., either quadrant IIa or IIIa)

                      ipa=i
                      ipb=i
                      ipc=i-1
                      ipd=i-1

                      if (zp1.gt.0) then !(i.e., quadrant IIa)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IIIa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., xp1<=0, so Q.IIb or Q.IIIb)

                      ipa=i
                      ipb=i
                      ipc=i+1
                      ipd=i+1

                      if (zp1.gt.0) then !(i.e., quadrant IIb)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IIIb)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(zp1).gt.abs(xp1))

                end if !closes condition on (yp1.gt.0)











              else !(i.e., |zp1|>=|xp1|,|yp1|, so zp1 cannot be 0)

                if (zp1.gt.0) then !(i.e., we are in either Q.Ia or Q.Ib or Q.IIa or Q.IIb)

                  kpa=k-1
                  kpb=k-1
                  kpc=k-1
                  kpd=k-1

                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((k+1).le.ke) then 
                    if (chrbdy2(i,j,k+1).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(xp1).gt.abs(yp1)) then !(i.e., |xp1|>|yp1|, so xp1 cannot be 0)
                    if (xp1.gt.0) then !(i.e., either quadrant Ia or IIa)

                      ipa=i
                      ipb=i
                      ipc=i-1
                      ipd=i-1

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIa)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., xp1<=0, so Q.Ib or IIb)

                      ipa=i
                      ipb=i
                      ipc=i+1
                      ipd=i+1

                      if (yp1.gt.0) then !(i.e., quadrant Ib)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIb)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |yp1|>=|xp1| - yp1=xp1=0 is not an issue)

                    if (yp1.gt.0) then !(i.e., either quadrant Ia or Ib)

                      jpa=j
                      jpb=j
                      jpc=j-1
                      jpd=j-1

                      if (xp1.gt.0) then !(i.e., quadrant Ia)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant Ib)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., yp1<=0, so Q.IIa or IIb)

                      jpa=j
                      jpb=j
                      jpc=j+1
                      jpd=j+1

                      if (xp1.gt.0) then !(i.e., quadrant IIa)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IIb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(xp1).gt.abs(yp1))


                else !(i.e., zp1<0, so we are in either Q.IIIa or Q.IIIb or Q.IVa or Q.IVb)

                  kpa=k+1
                  kpb=k+1
                  kpc=k+1
                  kpd=k+1


                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((k+1).le.ie) then 
                    if (chrbdy2(i,j,k+1).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(xp1).gt.abs(yp1)) then !(i.e., |xp1|>|yp1|, so xp1 cannot be 0)
                    if (xp1.gt.0) then !(i.e., either quadrant IIIa or IVa)

                      ipa=i
                      ipb=i
                      ipc=i-1
                      ipd=i-1

                      if (yp1.gt.0) then !(i.e., quadrant IVa)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIIa)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., xp1<=0, so Q.IIIb or IVb)

                      ipa=i
                      ipb=i
                      ipc=i+1
                      ipd=i+1

                      if (yp1.gt.0) then !(i.e., quadrant IVb)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIIb)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |yp1|>=|xp1|, it's possible that yp1=xp1=0)

                    if (yp1.gt.0) then !(i.e., either quadrant IVa or IVb)

                      jpa=j
                      jpb=j
                      jpc=j-1
                      jpd=j-1

                      if (xp1.gt.0) then !(i.e., quadrant IVa)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IVb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., yp1<=0, so Q.IIIa or IIIb)

                      jpa=j
                      jpb=j
                      jpc=j+1
                      jpd=j+1

                      if (xp1.gt.0) then !(i.e., quadrant IIIa)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IIIb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(xp1).gt.abs(yp1))

                end if !closes condition on (zp1.gt.0)






              end if !closes condition on ((abs(xp1).gt.abs(yp1)).and.(abs(xp1).gt.abs(zp1)))
           

            end if !closes condition on bdy_extrap_order.eq.1


          end if !closes condition on chrbdy(i,j,k).ne.ex





          if (chrbdy(i,j,k).ne.ex) then
            numbdypoints=numbdypoints+1
          end if


         end do
        end do
      end do

        return
        end


c------------------------------------------------------------------------------------------------------


c------------------------------------------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the Cartesian coordinates of each AdS boundary point pbdy where we want to extrapolate the value of functions 
c AND the coordinates of the outermost point pout used for radial extrapolation for each of the boundary points
c-------------------------------------------------------------------------------------------------------------------------

        subroutine xyz_bdy_out_radextrap(
     &                  xpbdy,ypbdy,zpbdy,
     &                  xpout,ypout,zpout,
     &                  chrbdy,
     &                  numbdypoints,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,ghost_width)

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        real*8 chrbdy(Nx,Ny,Nz)

        integer i,j,k,is,ie,js,je,ks,ke

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz
        integer lind

        real*8 x0,y0,z0,rho0,q

        real*8 dx,dy,dz
        real*8 xp1,yp1,zp1
        real*8 rhop1,chip1,xip1
        real*8 rhopbdy,chipbdy,xipbdy
        real*8 maxxyzp1
        integer numbdypoints
        real*8 xpbdy(numbdypoints)
        real*8 ypbdy(numbdypoints)
        real*8 zpbdy(numbdypoints)
        real*8 xpout(numbdypoints)
        real*8 ypout(numbdypoints)
        real*8 zpout(numbdypoints)

        real*8 PI
        parameter (PI=3.141592653589793d0)

!----------------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

        ! set index bounds for main loop
        is=2
        ie=Nx-1
        js=2
        je=Ny-1
        ks=2
        ke=Nz-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)

        lind=0
        do i=is,ie
         do j=js,je
          do k=ks,ke

            if (chrbdy(i,j,k).ne.ex) then
              lind=lind+1

              xp1=x(i)
              yp1=y(j)
              zp1=z(k)
              rhop1=sqrt(xp1**2+yp1**2+zp1**2)
              chip1=(1/PI)*acos(xp1/rhop1)
              if (zp1.lt.0) then
                xip1=(1/(2*PI))*(atan2(zp1,yp1)+2*PI)
              else
                xip1=(1/(2*PI))*atan2(zp1,yp1)
              end if

              xpout(lind)=xp1
              ypout(lind)=yp1
              zpout(lind)=zp1

              rhopbdy=1
              chipbdy=chip1
              xipbdy=xip1

              xpbdy(lind)=rhopbdy*cos(PI*chipbdy)
              ypbdy(lind)=rhopbdy*sin(PI*chipbdy)
     &              *cos(2*PI*xipbdy)
              zpbdy(lind)=rhopbdy*sin(PI*chipbdy)
     &              *sin(2*PI*xipbdy)


            end if !closes condition on (chrbdy(i,j,k).ne.ex)
          end do
         end do
        end do

        return
        end
c--------------------------------------------------------------------------------------