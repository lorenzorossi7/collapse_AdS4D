c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the scalar field phi of the boundary CFT
c----------------------------------------------------------------------

        subroutine bdyphi_radextrap(bdyphi,
     &                  leadordcoeff_phi1,
     &                  xpbdy,ypbdy,zpbdy,
     &                  chrbdy,numbdypoints,
     &                  bdy_extrap_order,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width)

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        integer numbdypoints
        integer bdy_extrap_order
        real*8 chrbdy(Nx,Ny,Nz)
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        integer i,j,k,is,ie,js,je,ks,ke,lind
        integer ipa,jpa,kpa
        integer ipb,jpb,kpb
        integer ipc,jpc,kpc
        integer ipd,jpd,kpd
        integer a,b,c,d

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz


        real*8 x0,y0,z0,rho0,q,chi0,xi0
        real*8 xpa,xpb,xpc,xpd
        real*8 ypa,ypb,ypc,ypd
        real*8 zpa,zpb,zpc,zpd

        real*8 dx,dy,dz

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 leadordcoeff_phi1(Nx,Ny,Nz)
        real*8 leadordcoeff_phi1_p1

        real*8 leadordcoeff_phi1_pa
        real*8 leadordcoeff_phi1_pb
        real*8 leadordcoeff_phi1_pc
        real*8 leadordcoeff_phi1_pd

        real*8 leadordcoeff_phi1_p2
        real*8 leadordcoeff_phi1_p3
        real*8 leadordcoeff_phi1_p4
        real*8 bdyphi(numbdypoints)

        real*8 xpbdy(numbdypoints)
        real*8 ypbdy(numbdypoints)
        real*8 zpbdy(numbdypoints)

        real*8 xp1,yp1,zp1
        real*8 xp2,yp2,zp2
        real*8 xp3,yp3,zp3
        real*8 xp4,yp4,zp4
        real*8 xex,yex,zex,rhoex,chiex,xiex
        real*8 rhop1,chip1,xip1
        real*8 rhop2,chip2,xip2
        real*8 maxxyzp1

        real*8 bilinear_interp
        real*8 firstord_extrap
        real*8 secondord_extrap
        real*8 thirdord_extrap
        real*8 fourthord_extrap
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
  
              leadordcoeff_phi1_p1=leadordcoeff_phi1(i,j,k)

              xex=xpbdy(lind)
              yex=ypbdy(lind)
              zex=zpbdy(lind)
              rhoex=1
              chiex=(1/PI)*acos(xex/rhoex)
              if (zex.lt.0) then
                xiex=(1/(2*PI))*(atan2(zex,yex)+2*PI)
              else
                xiex=(1/(2*PI))*atan2(zex,yex)
              end if



            if (bdy_extrap_order.eq.1) then

              chip2=chip1
              xip2=xip1


              if ((abs(xp1).gt.abs(yp1)).and.
     &            (abs(xp1).gt.abs(zp1))) then !(i.e., |xp1|>|yp1|,|zp1|, so xp1 cannot be 0)
            
                if (xp1.gt.0) then !(i.e., we are in the upper part of 
                                    !the spatial grid, called "a" sector)

                  ipa=i-1
                  ipb=i-1
                  ipc=i-1
                  ipd=i-1

                  xpa=x(ipa)
                  xpb=x(ipb)
                  xpc=x(ipc)
                  xpd=x(ipd)
                  xp2=x(ipa)

                  if (abs(yp1).gt.abs(zp1)) then !(i.e., |yp1|>|zp1|, so yp1 cannot be 0)
                    if (yp1.gt.0) then !(i.e., either quadrant Ia or IVa)

                      jpa=j
                      jpb=j
                      jpc=j-1
                      jpd=j-1

                      ypa=y(jpa)
                      ypb=y(jpb)
                      ypc=y(jpc)
                      ypd=y(jpd)

                      if (zp1.gt.0) then !(i.e., quadrant Ia)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k

                        zpa=z(kpa)
                        zpb=z(kpb)
                        zpc=z(kpc)
                        zpd=z(kpd)

                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      else !(i.e., zp1<=0: quadrant IVa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        zpa=z(kpa)
                        zpb=z(kpb)
                        zpc=z(kpc)
                        zpd=z(kpd)

                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    else !(i.e., yp1<=0, so Q.IIa or IIIa)

                      jpa=j
                      jpb=j
                      jpc=j+1
                      jpd=j+1

                      ypa=y(jpa)
                      ypb=y(jpb)
                      ypc=y(jpc)
                      ypd=y(jpd)

                      if (zp1.gt.0) then !(i.e., quadrant IIa)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k

                        zpa=z(kpa)
                        zpb=z(kpb)
                        zpc=z(kpc)
                        zpd=z(kpd)

                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      else !(i.e., zp1<=0: quadrant IIIa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        zpa=z(kpa)
                        zpb=z(kpb)
                        zpc=z(kpc)
                        zpd=z(kpd)

                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    end if

                  else !(i.e., |zp1|>=|yp1|, it's possible that zp1=yp1=0)

                    if (zp1.gt.0) then !(i.e., either quadrant Ia or IIa)

                      kpa=k
                      kpb=k
                      kpc=k-1
                      kpd=k-1

                      zpa=z(kpa)
                      zpb=z(kpb)
                      zpc=z(kpc)
                      zpd=z(kpd)

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j

                        ypa=y(jpa)
                        ypb=y(jpb)
                        ypc=y(jpc)
                        ypd=y(jpd)


                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


                      else !(i.e., yp1<=0: quadrant IIa)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        ypa=y(jpa)
                        ypb=y(jpb)
                        ypc=y(jpc)
                        ypd=y(jpd)


                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    else !(i.e., zp1<=0, so Q.IIIa or IVa)

                      kpa=k
                      kpb=k
                      kpc=k+1
                      kpd=k+1

                      zpa=z(kpa)
                      zpb=z(kpb)
                      zpc=z(kpc)
                      zpd=z(kpd)

                      if (yp1.gt.0) then !(i.e., quadrant IVa)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j

                        ypa=y(jpa)
                        ypb=y(jpb)
                        ypc=y(jpc)
                        ypd=y(jpd)


                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


                      else !(i.e., yp1<=0: quadrant IIIa)

                        if ((abs(yp1).lt.10.0d0**(-10)).and.
     &                      (abs(zp1).lt.10.0d0**(-10))) then
                          ypa=0.0d0
                          ypb=0.0d0
                          ypc=0.0d0
                          ypd=0.0d0

                          yp2=0.0d0
                          zp2=0.0d0
                          rhop2=abs(xp2)

                          leadordcoeff_phi1_p2=
     &                        leadordcoeff_phi1(i-1,j,k)


                        else

                          jpa=j
                          jpb=j+1
                          jpc=j+1
                          jpd=j

                          ypa=y(jpa)
                          ypb=y(jpb)
                          ypc=y(jpc)
                          ypd=y(jpd)


                          yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                          zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                          rhop2=abs(xp2/cos(PI*chip2))

                          leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                          leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                          leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                          leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                          leadordcoeff_phi1_p2=
     &                        bilinear_interp(
     &                             leadordcoeff_phi1_pa,
     &                             leadordcoeff_phi1_pb,
     &                             leadordcoeff_phi1_pc,
     &                             leadordcoeff_phi1_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                        end if

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


                      end if

                    end if

                  end if !closes condition on (abs(yp1).gt.abs(zp1))

                else !(i.e., xp1<0, i.e., we are in the lower part of 
                     !the spatial grid, called "b" sector)

                  ipa=i+1
                  ipb=i+1
                  ipc=i+1
                  ipd=i+1

                  xpa=x(ipa)
                  xpb=x(ipb)
                  xpc=x(ipc)
                  xpd=x(ipd)
                  xp2=x(ipa)

                  if (abs(yp1).gt.abs(zp1)) then !(i.e., |yp1|>|zp1|, so yp1 cannot be 0)
                    if (yp1.gt.0) then !(i.e., either quadrant Ia or IVa)
                      jpa=j
                      jpb=j
                      jpc=j-1
                      jpd=j-1

                      ypa=y(jpa)
                      ypb=y(jpb)
                      ypc=y(jpc)
                      ypd=y(jpd)

                      if (zp1.gt.0) then !(i.e., quadrant Ia)
                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k

                        zpa=z(kpa)
                        zpb=z(kpb)
                        zpc=z(kpc)
                        zpd=z(kpd)

                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      else !(i.e., zp1<=0: quadrant IVa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        zpa=z(kpa)
                        zpb=z(kpb)
                        zpc=z(kpc)
                        zpd=z(kpd)

                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    else !(i.e., yp1<=0, so Q.IIa or IIIa)

                      jpa=j
                      jpb=j
                      jpc=j+1
                      jpd=j+1

                      ypa=y(jpa)
                      ypb=y(jpb)
                      ypc=y(jpc)
                      ypd=y(jpd)

                      if (zp1.gt.0) then !(i.e., quadrant IIa)
                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k

                        zpa=z(kpa)
                        zpb=z(kpb)
                        zpc=z(kpc)
                        zpd=z(kpd)

                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      else !(i.e., zp1<=0: quadrant IIIa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        zpa=z(kpa)
                        zpb=z(kpb)
                        zpc=z(kpc)
                        zpd=z(kpd)

                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    end if

                  else !(i.e., |zp1|>=|yp1|, it's possible that zp1=yp1=0)

                    if (zp1.gt.0) then !(i.e., either quadrant Ia or IIa)

                      kpa=k
                      kpb=k
                      kpc=k-1
                      kpd=k-1

                      zpa=z(kpa)
                      zpb=z(kpb)
                      zpc=z(kpc)
                      zpd=z(kpd)

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j

                        ypa=y(jpa)
                        ypb=y(jpb)
                        ypc=y(jpc)
                        ypd=y(jpd)


                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


                      else !(i.e., yp1<=0: quadrant IIa)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        ypa=y(jpa)
                        ypb=y(jpb)
                        ypc=y(jpc)
                        ypd=y(jpd)


                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    else !(i.e., zp1<=0, so Q.IIIa or IVa)

                      kpa=k
                      kpb=k
                      kpc=k+1
                      kpd=k+1

                      zpa=z(kpa)
                      zpb=z(kpb)
                      zpc=z(kpc)
                      zpd=z(kpd)

                      if (yp1.gt.0) then !(i.e., quadrant IVa)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j

                        ypa=y(jpa)
                        ypb=y(jpb)
                        ypc=y(jpc)
                        ypd=y(jpd)


                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


                      else !(i.e., yp1<=0: quadrant IIIa)

                        if ((abs(yp1).lt.10.0d0**(-10)).and.
     &                      (abs(zp1).lt.10.0d0**(-10))) then
                          ypa=0.0d0
                          ypb=0.0d0
                          ypc=0.0d0
                          ypd=0.0d0

                          yp2=0.0d0
                          zp2=0.0d0
                          rhop2=abs(xp2)

                          leadordcoeff_phi1_p2=
     &                        leadordcoeff_phi1(i+1,j,k)


                        else

                          jpa=j
                          jpb=j+1
                          jpc=j+1
                          jpd=j

                          ypa=y(jpa)
                          ypb=y(jpb)
                          ypc=y(jpc)
                          ypd=y(jpd)


                          yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                          zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                          rhop2=abs(xp2/cos(PI*chip2))

                          leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                          leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                          leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                          leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                          leadordcoeff_phi1_p2=
     &                        bilinear_interp(
     &                             leadordcoeff_phi1_pa,
     &                             leadordcoeff_phi1_pb,
     &                             leadordcoeff_phi1_pc,
     &                             leadordcoeff_phi1_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                        end if

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


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

                  ypa=y(jpa)
                  ypb=y(jpb)
                  ypc=y(jpc)
                  ypd=y(jpd)
                  yp2=y(jpa)

                  if (abs(zp1).gt.abs(xp1)) then !(i.e., |zp1|>|xp1|, so zp1 cannot be 0)
                    if (zp1.gt.0) then !(i.e., either quadrant Ia or Ib)

                      kpa=k
                      kpb=k
                      kpc=k-1
                      kpd=k-1

                      zpa=z(kpa)
                      zpb=z(kpb)
                      zpc=z(kpc)
                      zpd=z(kpd)

                      if (xp1.gt.0) then !(i.e., quadrant Ia)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i

                        xpa=x(ipa)
                        xpb=x(ipb)
                        xpc=x(ipc)
                        xpd=x(ipd)

                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      else !(i.e., xp1<=0: quadrant Ib)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        xpa=x(ipa)
                        xpb=x(ipb)
                        xpc=x(ipc)
                        xpd=x(ipd)

                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    else !(i.e., zp1<=0, so Q.IVa or IVb)

                      kpa=k
                      kpb=k
                      kpc=k+1
                      kpd=k+1

                      zpa=z(kpa)
                      zpb=z(kpb)
                      zpc=z(kpc)
                      zpd=z(kpd)

                      if (xp1.gt.0) then !(i.e., quadrant IVa)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i

                        xpa=x(ipa)
                        xpb=x(ipb)
                        xpc=x(ipc)
                        xpd=x(ipd)

                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      else !(i.e., xp1<=0: quadrant IVb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        xpa=x(ipa)
                        xpb=x(ipb)
                        xpc=x(ipc)
                        xpd=x(ipd)

                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    end if

                  else !(i.e., |xp1|>=|zp1| - xp1=zp1=0 is not an issue)

                    if (xp1.gt.0) then !(i.e., either quadrant Ia or IVa)

                      ipa=i
                      ipb=i
                      ipc=i-1
                      ipd=i-1

                      xpa=x(ipa)
                      xpb=x(ipb)
                      xpc=x(ipc)
                      xpd=x(ipd)

                      if (zp1.gt.0) then !(i.e., quadrant Ia)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k

                        zpa=z(kpa)
                        zpb=z(kpb)
                        zpc=z(kpc)
                        zpd=z(kpd)


                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


                      else !(i.e., zp1<=0: quadrant IVa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        zpa=z(kpa)
                        zpb=z(kpb)
                        zpc=z(kpc)
                        zpd=z(kpd)


                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    else !(i.e., xp1<=0, so Q.Ib or IVb)

                      ipa=i
                      ipb=i
                      ipc=i+1
                      ipd=i+1

                      xpa=x(ipa)
                      xpb=x(ipb)
                      xpc=x(ipc)
                      xpd=x(ipd)

                      if (zp1.gt.0) then !(i.e., quadrant Ib)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k

                        zpa=z(kpa)
                        zpb=z(kpb)
                        zpc=z(kpc)
                        zpd=z(kpd)

                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


                      else !(i.e., zp1<=0: quadrant IVb)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        zpa=z(kpa)
                        zpb=z(kpb)
                        zpc=z(kpc)
                        zpd=z(kpd)

                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


                      end if

                    end if

                  end if !closes condition on (abs(zp1).gt.abs(xp1))

                else !(i.e., yp1<0, so we are in either Q.IIa or Q.IIb or Q.IIIa or Q.IIIb)

                  jpa=j+1
                  jpb=j+1
                  jpc=j+1
                  jpd=j+1

                  ypa=y(jpa)
                  ypb=y(jpb)
                  ypc=y(jpc)
                  ypd=y(jpd)
                  yp2=y(jpa)

                  if (abs(zp1).gt.abs(xp1)) then !(i.e., |zp1|>|xp1|, so zp1 cannot be 0)
                    if (zp1.gt.0) then !(i.e., either quadrant IIa or IIb)
                      kpa=k
                      kpb=k
                      kpc=k-1
                      kpd=k-1

                      zpa=z(kpa)
                      zpb=z(kpb)
                      zpc=z(kpc)
                      zpd=z(kpd)

                      if (xp1.gt.0) then !(i.e., quadrant Ia)
                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i

                        xpa=x(ipa)
                        xpb=x(ipb)
                        xpc=x(ipc)
                        xpd=x(ipd)

                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      else !(i.e., xp1<=0: quadrant IIb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        xpa=x(ipa)
                        xpb=x(ipb)
                        xpc=x(ipc)
                        xpd=x(ipd)

                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    else !(i.e., zp1<=0, so Q.IIIa or IIIb)

                      kpa=k
                      kpb=k
                      kpc=k+1
                      kpd=k+1

                      zpa=z(kpa)
                      zpb=z(kpb)
                      zpc=z(kpc)
                      zpd=z(kpd)

                      if (xp1.gt.0) then !(i.e., quadrant IIIa)
                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i

                        xpa=x(ipa)
                        xpb=x(ipb)
                        xpc=x(ipc)
                        xpd=x(ipd)

                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      else !(i.e., xp1<=0: quadrant IIIb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        xpa=x(ipa)
                        xpb=x(ipb)
                        xpc=x(ipc)
                        xpd=x(ipd)

                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    end if

                  else !(i.e., |xp1|>=|zp1| - the case xp1=zp1=0 is not an issue)

                    if (xp1.gt.0) then !(i.e., either quadrant IIa or IIIa)

                      ipa=i
                      ipb=i
                      ipc=i-1
                      ipd=i-1

                      xpa=x(ipa)
                      xpb=x(ipb)
                      xpc=x(ipc)
                      xpd=x(ipd)

                      if (zp1.gt.0) then !(i.e., quadrant IIa)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k

                        zpa=z(kpa)
                        zpb=z(kpb)
                        zpc=z(kpc)
                        zpd=z(kpd)


                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


                      else !(i.e., zp1<=0: quadrant IIIa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        zpa=z(kpa)
                        zpb=z(kpb)
                        zpc=z(kpc)
                        zpd=z(kpd)


                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    else !(i.e., xp1<=0, so Q.IIb or IIIb)

                      ipa=i
                      ipb=i
                      ipc=i+1
                      ipd=i+1

                      xpa=x(ipa)
                      xpb=x(ipb)
                      xpc=x(ipc)
                      xpd=x(ipd)

                      if (zp1.gt.0) then !(i.e., quadrant IIb)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k

                        zpa=z(kpa)
                        zpb=z(kpb)
                        zpc=z(kpc)
                        zpd=z(kpd)


                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


                      else !(i.e., zp1<=0: quadrant IIIb)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        zpa=z(kpa)
                        zpb=z(kpb)
                        zpc=z(kpc)
                        zpd=z(kpd)


                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)


                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


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

                  zpa=z(kpa)
                  zpb=z(kpb)
                  zpc=z(kpc)
                  zpd=z(kpd)
                  zp2=z(kpa)

                  if (abs(xp1).gt.abs(yp1)) then !(i.e., |xp1|>|yp1|, so xp1 cannot be 0)
                    if (xp1.gt.0) then !(i.e., either quadrant Ia or IIa)

                      ipa=i
                      ipb=i
                      ipc=i-1
                      ipd=i-1

                      xpa=x(ipa)
                      xpb=x(ipb)
                      xpc=x(ipc)
                      xpd=x(ipd)

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j

                        ypa=y(jpa)
                        ypb=y(jpb)
                        ypc=y(jpc)
                        ypd=y(jpd)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      else !(i.e., yp1<=0: quadrant IIa)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        ypa=y(jpa)
                        ypb=y(jpb)
                        ypc=y(jpc)
                        ypd=y(jpd)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    else !(i.e., xp1<=0, so Q.Ib or IIb)

                      ipa=i
                      ipb=i
                      ipc=i+1
                      ipd=i+1

                      xpa=x(ipa)
                      xpb=x(ipb)
                      xpc=x(ipc)
                      xpd=x(ipd)

                      if (yp1.gt.0) then !(i.e., quadrant Ib)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j

                        ypa=y(jpa)
                        ypb=y(jpb)
                        ypc=y(jpc)
                        ypd=y(jpd)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      else !(i.e., yp1<=0: quadrant IIb)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        ypa=y(jpa)
                        ypb=y(jpb)
                        ypc=y(jpc)
                        ypd=y(jpd)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    end if

                  else !(i.e., |yp1|>=|xp1| - yp1=xp1=0 is not an issue)

                    if (yp1.gt.0) then !(i.e., either quadrant Ia or Ib)

                      jpa=j
                      jpb=j
                      jpc=j-1
                      jpd=j-1

                      ypa=y(jpa)
                      ypb=y(jpb)
                      ypc=y(jpc)
                      ypd=y(jpd)

                      if (xp1.gt.0) then !(i.e., quadrant Ia)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i

                        xpa=x(ipa)
                        xpb=x(ipb)
                        xpc=x(ipc)
                        xpd=x(ipd)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


                      else !(i.e., xp1<=0: quadrant Ib)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        xpa=x(ipa)
                        xpb=x(ipb)
                        xpc=x(ipc)
                        xpd=x(ipd)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    else !(i.e., yp1<=0, so Q.IIa or IIb)

                      jpa=j
                      jpb=j
                      jpc=j+1
                      jpd=j+1

                      ypa=y(jpa)
                      ypb=y(jpb)
                      ypc=y(jpc)
                      ypd=y(jpd)

                      if (xp1.gt.0) then !(i.e., quadrant IIa)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i

                        xpa=x(ipa)
                        xpb=x(ipb)
                        xpc=x(ipc)
                        xpd=x(ipd)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


                      else !(i.e., xp1<=0: quadrant IIb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        xpa=x(ipa)
                        xpb=x(ipb)
                        xpc=x(ipc)
                        xpd=x(ipd)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


                      end if

                    end if

                  end if !closes condition on (abs(xp1).gt.abs(yp1))

                else !(i.e., zp1<0, so we are in either Q.IIIa or Q.IIIb or Q.IVa or Q.IVb)

                  kpa=k+1
                  kpb=k+1
                  kpc=k+1
                  kpd=k+1

                  zpa=z(kpa)
                  zpb=z(kpb)
                  zpc=z(kpc)
                  zpd=z(kpd)
                  zp2=z(kpa)

                  if (abs(xp1).gt.abs(yp1)) then !(i.e., |xp1|>|yp1|, so xp1 cannot be 0)
                    if (xp1.gt.0) then !(i.e., either quadrant IIIa or IVa)
                      ipa=i
                      ipb=i
                      ipc=i-1
                      ipd=i-1

                      xpa=x(ipa)
                      xpb=x(ipb)
                      xpc=x(ipc)
                      xpd=x(ipd)

                      if (yp1.gt.0) then !(i.e., quadrant IVa)
                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j

                        ypa=y(jpa)
                        ypb=y(jpb)
                        ypc=y(jpc)
                        ypd=y(jpd)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      else !(i.e., yp1<=0: quadrant IIIa)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        ypa=y(jpa)
                        ypb=y(jpb)
                        ypc=y(jpc)
                        ypd=y(jpd)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    else !(i.e., xp1<=0, so Q.IIIb or IVb)

                      ipa=i
                      ipb=i
                      ipc=i+1
                      ipd=i+1

                      xpa=x(ipa)
                      xpb=x(ipb)
                      xpc=x(ipc)
                      xpd=x(ipd)

                      if (yp1.gt.0) then !(i.e., quadrant IVb)
                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j

                        ypa=y(jpa)
                        ypb=y(jpb)
                        ypc=y(jpc)
                        ypd=y(jpd)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      else !(i.e., yp1<=0: quadrant IIIb)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        ypa=y(jpa)
                        ypb=y(jpb)
                        ypc=y(jpc)
                        ypd=y(jpd)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    end if

                  else !(i.e., |yp1|>=|xp1| - yp1=xp1=0 is not an issue)

                    if (yp1.gt.0) then !(i.e., either quadrant IVa or IVb)

                      jpa=j
                      jpb=j
                      jpc=j-1
                      jpd=j-1

                      ypa=y(jpa)
                      ypb=y(jpb)
                      ypc=y(jpc)
                      ypd=y(jpd)

                      if (xp1.gt.0) then !(i.e., quadrant IVa)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i

                        xpa=x(ipa)
                        xpb=x(ipb)
                        xpc=x(ipc)
                        xpd=x(ipd)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


                      else !(i.e., xp1<=0: quadrant IVb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        xpa=x(ipa)
                        xpb=x(ipb)
                        xpc=x(ipc)
                        xpd=x(ipd)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)

                      end if

                    else !(i.e., yp1<=0, so Q.IIIa or IIIb)

                      jpa=j
                      jpb=j
                      jpc=j+1
                      jpd=j+1

                      ypa=y(jpa)
                      ypb=y(jpb)
                      ypc=y(jpc)
                      ypd=y(jpd)

                      if (xp1.gt.0) then !(i.e., quadrant IIIa)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i

                        xpa=x(ipa)
                        xpb=x(ipb)
                        xpc=x(ipc)
                        xpd=x(ipd)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


                      else !(i.e., xp1<=0: quadrant IIIb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        xpa=x(ipa)
                        xpb=x(ipb)
                        xpc=x(ipc)
                        xpd=x(ipd)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        leadordcoeff_phi1_pa=
     &                        leadordcoeff_phi1(ipa,jpa,kpa)
                        leadordcoeff_phi1_pb=
     &                        leadordcoeff_phi1(ipb,jpb,kpb)
                        leadordcoeff_phi1_pc=
     &                        leadordcoeff_phi1(ipc,jpc,kpc)
                        leadordcoeff_phi1_pd=
     &                        leadordcoeff_phi1(ipd,jpd,kpd)

                        leadordcoeff_phi1_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_phi1_pa,
     &                           leadordcoeff_phi1_pb,
     &                           leadordcoeff_phi1_pc,
     &                           leadordcoeff_phi1_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        bdyphi(lind)=
     &                      firstord_extrap(
     &                           leadordcoeff_phi1_p1,
     &                           leadordcoeff_phi1_p2,
     &                           rhop1,rhop2,rhoex)


                      end if

                    end if

                  end if !closes condition on (abs(xp1).gt.abs(yp1))

                end if !closes condition on (zp1.gt.0)






              end if !closes condition on ((abs(xp1).gt.abs(yp1)).and.(abs(xp1).gt.abs(zp1)))



!                  bdyphi(lind)=leadordcoeff_phi1_p1  !TEST
!              write(*,*) "lind-1,xp1,yp1,zp1",lind-1,xp1,yp1,zp1
!             write(*,*) "bdyphi(lind)=",bdyphi(lind)

            end if !closes condition on bdy_extrap_order.eq.1


           end if !closes condition on chrbdy(i,j,k).ne.ex




          end do
         end do
        end do


        return
        end
c--------------------------------------------------------------------------------------


c----------------------------------------------------------------------
c bilinear interpolation on a plane using values 
c Tpa at (xpa,ypa), Tpb at (xpb,ypb), Tpc at (xpc,ypc), Tpd at (xpd,ypd)  
c to obtain the value at (xp0,yp0) 
c----------------------------------------------------------------------
        real*8 function bilinear_interp(Tpa,Tpb,Tpc,Tpd,
     &                                  xpa,xpb,xpc,xpd,
     &                                  ypa,ypb,ypc,ypd,
     &                                  xp0,yp0)
        implicit none
        real*8 Tpa,Tpb,Tpc,Tpd
        real*8 xpa,xpb,xpc,xpd
        real*8 ypa,ypb,ypc,ypd
        real*8 xp0,yp0

        !--------------------------------------------------------------

        bilinear_interp=(1/((xpc-xpa)*(ypc-ypa)))*
     &    (Tpa*(xpc-xp0)*(ypc-yp0)+
     &     Tpb*(xpc-xp0)*(yp0-ypa)+
     &     Tpc*(xp0-xpa)*(yp0-ypa)+
     &     Tpd*(xp0-xpa)*(ypc-yp0))

        return
        end
c--------------------------------------------------------------------------------------