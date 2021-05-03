c-------------------------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the asymptotic quasilocal stress-energy of AdS4D_polar  
c using a 1-rho expansion about rho=1 AT the boundary through extrapolation.
c The tensor components are given in spherical polar coordinates.
c-------------------------------------------------------------------------------------

        subroutine quasiset_radextrap(
     &                  quasiset_tt,quasiset_tchi,quasiset_txi,
     &                  quasiset_chichi,quasiset_chixi,
     &                  quasiset_xixi,
     &                  quasiset_trace,
     &                  quasiset_massdensity,
     &                  quasiset_angmomdensityx,
     &                  quasiset_angmomdensityy,
     &                  quasiset_angmomdensityz,
     &                  quasiset_tt_ll,quasiset_tchi_ll,quasiset_txi_ll,
     &                  quasiset_chichi_ll,quasiset_chixi_ll,
     &                  quasiset_xixi_ll,
     &                  quasiset_tracell,
     &                  quasiset_massdensityll,
     &                  quasiset_angmomdensityxll,
     &                  quasiset_angmomdensityyll,
     &                  quasiset_angmomdensityzll,
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
        real*8 gb_tt_np1(Nx,Ny,Nz),gb_tt_n(Nx,Ny,Nz),gb_tt_nm1(Nx,Ny,Nz)
        real*8 gb_tx_np1(Nx,Ny,Nz),gb_tx_n(Nx,Ny,Nz),gb_tx_nm1(Nx,Ny,Nz)
        real*8 gb_ty_np1(Nx,Ny,Nz),gb_ty_n(Nx,Ny,Nz),gb_ty_nm1(Nx,Ny,Nz)
        real*8 gb_tz_np1(Nx,Ny,Nz),gb_tz_n(Nx,Ny,Nz),gb_tz_nm1(Nx,Ny,Nz)
        real*8 gb_xx_np1(Nx,Ny,Nz),gb_xx_n(Nx,Ny,Nz),gb_xx_nm1(Nx,Ny,Nz)
        real*8 gb_xy_np1(Nx,Ny,Nz),gb_xy_n(Nx,Ny,Nz),gb_xy_nm1(Nx,Ny,Nz)
        real*8 gb_xz_np1(Nx,Ny,Nz),gb_xz_n(Nx,Ny,Nz),gb_xz_nm1(Nx,Ny,Nz)
        real*8 gb_yy_np1(Nx,Ny,Nz),gb_yy_n(Nx,Ny,Nz),gb_yy_nm1(Nx,Ny,Nz)
        real*8 gb_yz_np1(Nx,Ny,Nz),gb_yz_n(Nx,Ny,Nz),gb_yz_nm1(Nx,Ny,Nz)
        real*8 gb_zz_np1(Nx,Ny,Nz),gb_zz_n(Nx,Ny,Nz),gb_zz_nm1(Nx,Ny,Nz)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex
        real*8 dx,dy,dz

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 quasiset_tt_ll(Nx,Ny,Nz),quasiset_tchi_ll(Nx,Ny,Nz)
        real*8 quasiset_txi_ll(Nx,Ny,Nz),quasiset_chichi_ll(Nx,Ny,Nz)
        real*8 quasiset_chixi_ll(Nx,Ny,Nz),quasiset_xixi_ll(Nx,Ny,Nz)
        real*8 quasiset_tracell(Nx,Ny,Nz)
        real*8 quasiset_massdensityll(Nx,Ny,Nz)
        real*8 quasiset_angmomdensityxll(Nx,Ny,Nz)
        real*8 quasiset_angmomdensityyll(Nx,Ny,Nz)
        real*8 quasiset_angmomdensityzll(Nx,Ny,Nz)

        real*8 quasiset_tt_p1,quasiset_tchi_p1
        real*8 quasiset_txi_p1,quasiset_chichi_p1
        real*8 quasiset_chixi_p1,quasiset_xixi_p1
        real*8 quasiset_trace_p1
        real*8 quasiset_massdensity_p1
        real*8 quasiset_angmomdensityx_p1
        real*8 quasiset_angmomdensityy_p1
        real*8 quasiset_angmomdensityz_p1

        real*8 quasiset_tt_pa,quasiset_tchi_pa
        real*8 quasiset_txi_pa,quasiset_chichi_pa
        real*8 quasiset_chixi_pa,quasiset_xixi_pa
        real*8 quasiset_trace_pa
        real*8 quasiset_massdensity_pa
        real*8 quasiset_angmomdensityx_pa
        real*8 quasiset_angmomdensityy_pa
        real*8 quasiset_angmomdensityz_pa

        real*8 quasiset_tt_pb,quasiset_tchi_pb
        real*8 quasiset_txi_pb,quasiset_chichi_pb
        real*8 quasiset_chixi_pb,quasiset_xixi_pb
        real*8 quasiset_trace_pb
        real*8 quasiset_massdensity_pb
        real*8 quasiset_angmomdensityx_pb
        real*8 quasiset_angmomdensityy_pb
        real*8 quasiset_angmomdensityz_pb

        real*8 quasiset_tt_pc,quasiset_tchi_pc
        real*8 quasiset_txi_pc,quasiset_chichi_pc
        real*8 quasiset_chixi_pc,quasiset_xixi_pc
        real*8 quasiset_trace_pc
        real*8 quasiset_massdensity_pc
        real*8 quasiset_angmomdensityx_pc
        real*8 quasiset_angmomdensityy_pc
        real*8 quasiset_angmomdensityz_pc

        real*8 quasiset_tt_pd,quasiset_tchi_pd
        real*8 quasiset_txi_pd,quasiset_chichi_pd
        real*8 quasiset_chixi_pd,quasiset_xixi_pd
        real*8 quasiset_trace_pd
        real*8 quasiset_massdensity_pd
        real*8 quasiset_angmomdensityx_pd
        real*8 quasiset_angmomdensityy_pd
        real*8 quasiset_angmomdensityz_pd

        real*8 quasiset_tt_p2,quasiset_tchi_p2
        real*8 quasiset_txi_p2,quasiset_chichi_p2
        real*8 quasiset_chixi_p2,quasiset_xixi_p2
        real*8 quasiset_trace_p2
        real*8 quasiset_massdensity_p2
        real*8 quasiset_angmomdensityx_p2
        real*8 quasiset_angmomdensityy_p2
        real*8 quasiset_angmomdensityz_p2

        real*8 quasiset_tt_p3,quasiset_tchi_p3
        real*8 quasiset_txi_p3,quasiset_chichi_p3
        real*8 quasiset_chixi_p3,quasiset_xixi_p3
        real*8 quasiset_trace_p3
        real*8 quasiset_massdensity_p3
        real*8 quasiset_angmomdensityx_p3
        real*8 quasiset_angmomdensityy_p3
        real*8 quasiset_angmomdensityz_p3

        real*8 quasiset_tt_p4,quasiset_tchi_p4
        real*8 quasiset_txi_p4,quasiset_chichi_p4
        real*8 quasiset_chixi_p4,quasiset_xixi_p4
        real*8 quasiset_trace_p4
        real*8 quasiset_massdensity_p4
        real*8 quasiset_angmomdensityx_p4
        real*8 quasiset_angmomdensityy_p4
        real*8 quasiset_angmomdensityz_p4


        real*8 gamma0sphbdy_uu_tt
        real*8 gamma0sphbdy_uu_tchi
        real*8 gamma0sphbdy_uu_txi
        real*8 gamma0sphbdy_uu_chichi
        real*8 gamma0sphbdy_uu_chixi
        real*8 gamma0sphbdy_uu_xixi

        real*8 quasiset_tt(numbdypoints),quasiset_tchi(numbdypoints)
        real*8 quasiset_txi(numbdypoints),quasiset_chichi(numbdypoints)
        real*8 quasiset_chixi(numbdypoints),quasiset_xixi(numbdypoints)
        real*8 quasiset_trace(numbdypoints)
        real*8 quasiset_massdensity(numbdypoints)
        real*8 quasiset_angmomdensityx(numbdypoints)
        real*8 quasiset_angmomdensityy(numbdypoints)
        real*8 quasiset_angmomdensityz(numbdypoints)

        real*8 xpbdy(numbdypoints)
        real*8 ypbdy(numbdypoints)
        real*8 zpbdy(numbdypoints)

        real*8 rhoextrap
        real*8 chiextrap(numbdypoints)
        real*8 xiextrap(numbdypoints)

        real*8 chiextrap_min,chiextrap_max
        real*8 xiextrap_min,xiextrap_max

        integer bdy_Nchi,bdy_Nxi

!        real*8 chibdy(bdy_Nchi)
!        real*8 xibdy(bdy_Nxi)

        integer i,j,k,is,ie,js,je,ks,ke,lind,m,e,increase

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz

        integer ipa,jpa,kpa
        integer ipb,jpb,kpb
        integer ipc,jpc,kpc
        integer ipd,jpd,kpd

        real*8 x0,y0,z0,rho0,q
        real*8 xpa,xpb,xpc,xpd
        real*8 ypa,ypb,ypc,ypd
        real*8 zpa,zpb,zpc,zpd

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

        real*8 dp1p2

        real*8 AdS_mass

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
  
              quasiset_tt_p1=
     &                        quasiset_tt_ll(i,j,k)
              quasiset_tchi_p1=
     &                        quasiset_tchi_ll(i,j,k)
              quasiset_txi_p1=
     &                        quasiset_txi_ll(i,j,k)
              quasiset_chichi_p1=
     &                        quasiset_chichi_ll(i,j,k)
              quasiset_chixi_p1=
     &                        quasiset_chixi_ll(i,j,k)
              quasiset_xixi_p1=
     &                        quasiset_xixi_ll(i,j,k)
              quasiset_trace_p1=
     &                        quasiset_tracell(i,j,k)
              quasiset_massdensity_p1=
     &                        quasiset_massdensityll(i,j,k)
              quasiset_angmomdensityx_p1=
     &                        quasiset_angmomdensityxll(i,j,k)
              quasiset_angmomdensityy_p1=
     &                        quasiset_angmomdensityyll(i,j,k)
              quasiset_angmomdensityz_p1=
     &                        quasiset_angmomdensityzll(i,j,k)

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

              !inverse of conformal metric on AdS boundary (needed for trace) at extrapolated point
              gamma0sphbdy_uu_tt=-1
              gamma0sphbdy_uu_tchi=0
              gamma0sphbdy_uu_txi=0
              gamma0sphbdy_uu_chichi=1/(PI**2)
              gamma0sphbdy_uu_chixi=0
              gamma0sphbdy_uu_xixi=1/((sin(PI*chiex))**2)/4/PI**2



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

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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

                          quasiset_tt_p2=
     &                         quasiset_tt_ll(i-1,j,k)
                          quasiset_tchi_p2=
     &                         quasiset_tchi_ll(i-1,j,k)
                          quasiset_txi_p2=
     &                         quasiset_txi_ll(i-1,j,k)
                          quasiset_chichi_p2=
     &                         quasiset_chichi_ll(i-1,j,k)
                          quasiset_chixi_p2=
     &                         quasiset_chixi_ll(i-1,j,k)
                          quasiset_xixi_p2=
     &                         quasiset_xixi_ll(i-1,j,k)
                          quasiset_trace_p2=
     &                         quasiset_tracell(i-1,j,k)
                          quasiset_massdensity_p2=
     &                         quasiset_massdensityll(i-1,j,k)
                          quasiset_angmomdensityx_p2=
     &                         quasiset_angmomdensityxll(i-1,j,k)
                          quasiset_angmomdensityy_p2=
     &                         quasiset_angmomdensityyll(i-1,j,k)
                          quasiset_angmomdensityz_p2=
     &                         quasiset_angmomdensityzll(i-1,j,k)


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

                          quasiset_tt_pa=
     &                         quasiset_tt_ll(ipa,jpa,kpa)
                          quasiset_tt_pb=
     &                         quasiset_tt_ll(ipb,jpb,kpb)
                          quasiset_tt_pc=
     &                         quasiset_tt_ll(ipc,jpc,kpc)
                          quasiset_tt_pd=
     &                         quasiset_tt_ll(ipd,jpd,kpd)

                          quasiset_tt_p2=
     &                        bilinear_interp(
     &                             quasiset_tt_pa,
     &                             quasiset_tt_pb,
     &                             quasiset_tt_pc,
     &                             quasiset_tt_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_tchi_pa=
     &                         quasiset_tchi_ll(ipa,jpa,kpa)
                          quasiset_tchi_pb=
     &                         quasiset_tchi_ll(ipb,jpb,kpb)
                          quasiset_tchi_pc=
     &                         quasiset_tchi_ll(ipc,jpc,kpc)
                          quasiset_tchi_pd=
     &                         quasiset_tchi_ll(ipd,jpd,kpd)

                          quasiset_tchi_p2=
     &                        bilinear_interp(
     &                             quasiset_tchi_pa,
     &                             quasiset_tchi_pb,
     &                             quasiset_tchi_pc,
     &                             quasiset_tchi_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_txi_pa=
     &                         quasiset_txi_ll(ipa,jpa,kpa)
                          quasiset_txi_pb=
     &                         quasiset_txi_ll(ipb,jpb,kpb)
                          quasiset_txi_pc=
     &                         quasiset_txi_ll(ipc,jpc,kpc)
                          quasiset_txi_pd=
     &                         quasiset_txi_ll(ipd,jpd,kpd)

                          quasiset_txi_p2=
     &                        bilinear_interp(
     &                             quasiset_txi_pa,
     &                             quasiset_txi_pb,
     &                             quasiset_txi_pc,
     &                             quasiset_txi_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_chichi_pa=
     &                         quasiset_chichi_ll(ipa,jpa,kpa)
                          quasiset_chichi_pb=
     &                         quasiset_chichi_ll(ipb,jpb,kpb)
                          quasiset_chichi_pc=
     &                         quasiset_chichi_ll(ipc,jpc,kpc)
                          quasiset_chichi_pd=
     &                         quasiset_chichi_ll(ipd,jpd,kpd)

                          quasiset_chichi_p2=
     &                        bilinear_interp(
     &                             quasiset_chichi_pa,
     &                             quasiset_chichi_pb,
     &                             quasiset_chichi_pc,
     &                             quasiset_chichi_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_chixi_pa=
     &                         quasiset_chixi_ll(ipa,jpa,kpa)
                          quasiset_chixi_pb=
     &                         quasiset_chixi_ll(ipb,jpb,kpb)
                          quasiset_chixi_pc=
     &                         quasiset_chixi_ll(ipc,jpc,kpc)
                          quasiset_chixi_pd=
     &                         quasiset_chixi_ll(ipd,jpd,kpd)

                          quasiset_chixi_p2=
     &                        bilinear_interp(
     &                             quasiset_chixi_pa,
     &                             quasiset_chixi_pb,
     &                             quasiset_chixi_pc,
     &                             quasiset_chixi_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_xixi_pa=
     &                         quasiset_xixi_ll(ipa,jpa,kpa)
                          quasiset_xixi_pb=
     &                         quasiset_xixi_ll(ipb,jpb,kpb)
                          quasiset_xixi_pc=
     &                         quasiset_xixi_ll(ipc,jpc,kpc)
                          quasiset_xixi_pd=
     &                         quasiset_xixi_ll(ipd,jpd,kpd)

                          quasiset_xixi_p2=
     &                        bilinear_interp(
     &                             quasiset_xixi_pa,
     &                             quasiset_xixi_pb,
     &                             quasiset_xixi_pc,
     &                             quasiset_xixi_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_trace_pa=
     &                         quasiset_tracell(ipa,jpa,kpa)
                          quasiset_trace_pb=
     &                         quasiset_tracell(ipb,jpb,kpb)
                          quasiset_trace_pc=
     &                         quasiset_tracell(ipc,jpc,kpc)
                          quasiset_trace_pd=
     &                         quasiset_tracell(ipd,jpd,kpd)

                          quasiset_trace_p2=
     &                        bilinear_interp(
     &                             quasiset_trace_pa,
     &                             quasiset_trace_pb,
     &                             quasiset_trace_pc,
     &                             quasiset_trace_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_massdensity_pa=
     &                         quasiset_massdensityll(ipa,jpa,kpa)
                          quasiset_massdensity_pb=
     &                         quasiset_massdensityll(ipb,jpb,kpb)
                          quasiset_massdensity_pc=
     &                         quasiset_massdensityll(ipc,jpc,kpc)
                          quasiset_massdensity_pd=
     &                         quasiset_massdensityll(ipd,jpd,kpd)

                          quasiset_massdensity_p2=
     &                        bilinear_interp(
     &                             quasiset_massdensity_pa,
     &                             quasiset_massdensity_pb,
     &                             quasiset_massdensity_pc,
     &                             quasiset_massdensity_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_angmomdensityx_pa=
     &                         quasiset_angmomdensityxll(ipa,jpa,kpa)
                          quasiset_angmomdensityx_pb=
     &                         quasiset_angmomdensityxll(ipb,jpb,kpb)
                          quasiset_angmomdensityx_pc=
     &                         quasiset_angmomdensityxll(ipc,jpc,kpc)
                          quasiset_angmomdensityx_pd=
     &                         quasiset_angmomdensityxll(ipd,jpd,kpd)

                          quasiset_angmomdensityx_p2=
     &                        bilinear_interp(
     &                             quasiset_angmomdensityx_pa,
     &                             quasiset_angmomdensityx_pb,
     &                             quasiset_angmomdensityx_pc,
     &                             quasiset_angmomdensityx_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_angmomdensityy_pa=
     &                         quasiset_angmomdensityyll(ipa,jpa,kpa)
                          quasiset_angmomdensityy_pb=
     &                         quasiset_angmomdensityyll(ipb,jpb,kpb)
                          quasiset_angmomdensityy_pc=
     &                         quasiset_angmomdensityyll(ipc,jpc,kpc)
                          quasiset_angmomdensityy_pd=
     &                         quasiset_angmomdensityyll(ipd,jpd,kpd)

                          quasiset_angmomdensityy_p2=
     &                        bilinear_interp(
     &                             quasiset_angmomdensityy_pa,
     &                             quasiset_angmomdensityy_pb,
     &                             quasiset_angmomdensityy_pc,
     &                             quasiset_angmomdensityy_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_angmomdensityz_pa=
     &                         quasiset_angmomdensityzll(ipa,jpa,kpa)
                          quasiset_angmomdensityz_pb=
     &                         quasiset_angmomdensityzll(ipb,jpb,kpb)
                          quasiset_angmomdensityz_pc=
     &                         quasiset_angmomdensityzll(ipc,jpc,kpc)
                          quasiset_angmomdensityz_pd=
     &                         quasiset_angmomdensityzll(ipd,jpd,kpd)

                          quasiset_angmomdensityz_p2=
     &                        bilinear_interp(
     &                             quasiset_angmomdensityz_pa,
     &                             quasiset_angmomdensityz_pb,
     &                             quasiset_angmomdensityz_pc,
     &                             quasiset_angmomdensityz_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                        end if

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           ypa,ypb,ypc,ypd,
     &                           zpa,zpb,zpc,zpd,
     &                           yp2,
     &                           zp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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

                          quasiset_tt_p2=
     &                         quasiset_tt_ll(i+1,j,k)
                          quasiset_tchi_p2=
     &                         quasiset_tchi_ll(i+1,j,k)
                          quasiset_txi_p2=
     &                         quasiset_txi_ll(i+1,j,k)
                          quasiset_chichi_p2=
     &                         quasiset_chichi_ll(i+1,j,k)
                          quasiset_chixi_p2=
     &                         quasiset_chixi_ll(i+1,j,k)
                          quasiset_xixi_p2=
     &                         quasiset_xixi_ll(i+1,j,k)
                          quasiset_trace_p2=
     &                         quasiset_tracell(i+1,j,k)
                          quasiset_massdensity_p2=
     &                         quasiset_massdensityll(i+1,j,k)
                          quasiset_angmomdensityx_p2=
     &                         quasiset_angmomdensityxll(i+1,j,k)
                          quasiset_angmomdensityy_p2=
     &                         quasiset_angmomdensityyll(i+1,j,k)
                          quasiset_angmomdensityz_p2=
     &                         quasiset_angmomdensityzll(i+1,j,k)


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

                          quasiset_tt_p2=
     &                        bilinear_interp(
     &                             quasiset_tt_pa,
     &                             quasiset_tt_pb,
     &                             quasiset_tt_pc,
     &                             quasiset_tt_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_tchi_pa=
     &                         quasiset_tchi_ll(ipa,jpa,kpa)
                          quasiset_tchi_pb=
     &                         quasiset_tchi_ll(ipb,jpb,kpb)
                          quasiset_tchi_pc=
     &                         quasiset_tchi_ll(ipc,jpc,kpc)
                          quasiset_tchi_pd=
     &                         quasiset_tchi_ll(ipd,jpd,kpd)

                          quasiset_tchi_p2=
     &                        bilinear_interp(
     &                             quasiset_tchi_pa,
     &                             quasiset_tchi_pb,
     &                             quasiset_tchi_pc,
     &                             quasiset_tchi_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_txi_pa=
     &                         quasiset_txi_ll(ipa,jpa,kpa)
                          quasiset_txi_pb=
     &                         quasiset_txi_ll(ipb,jpb,kpb)
                          quasiset_txi_pc=
     &                         quasiset_txi_ll(ipc,jpc,kpc)
                          quasiset_txi_pd=
     &                         quasiset_txi_ll(ipd,jpd,kpd)

                          quasiset_txi_p2=
     &                        bilinear_interp(
     &                             quasiset_txi_pa,
     &                             quasiset_txi_pb,
     &                             quasiset_txi_pc,
     &                             quasiset_txi_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_chichi_pa=
     &                         quasiset_chichi_ll(ipa,jpa,kpa)
                          quasiset_chichi_pb=
     &                         quasiset_chichi_ll(ipb,jpb,kpb)
                          quasiset_chichi_pc=
     &                         quasiset_chichi_ll(ipc,jpc,kpc)
                          quasiset_chichi_pd=
     &                         quasiset_chichi_ll(ipd,jpd,kpd)

                          quasiset_chichi_p2=
     &                        bilinear_interp(
     &                             quasiset_chichi_pa,
     &                             quasiset_chichi_pb,
     &                             quasiset_chichi_pc,
     &                             quasiset_chichi_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_chixi_pa=
     &                         quasiset_chixi_ll(ipa,jpa,kpa)
                          quasiset_chixi_pb=
     &                         quasiset_chixi_ll(ipb,jpb,kpb)
                          quasiset_chixi_pc=
     &                         quasiset_chixi_ll(ipc,jpc,kpc)
                          quasiset_chixi_pd=
     &                         quasiset_chixi_ll(ipd,jpd,kpd)

                          quasiset_chixi_p2=
     &                        bilinear_interp(
     &                             quasiset_chixi_pa,
     &                             quasiset_chixi_pb,
     &                             quasiset_chixi_pc,
     &                             quasiset_chixi_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_xixi_pa=
     &                         quasiset_xixi_ll(ipa,jpa,kpa)
                          quasiset_xixi_pb=
     &                         quasiset_xixi_ll(ipb,jpb,kpb)
                          quasiset_xixi_pc=
     &                         quasiset_xixi_ll(ipc,jpc,kpc)
                          quasiset_xixi_pd=
     &                         quasiset_xixi_ll(ipd,jpd,kpd)

                          quasiset_xixi_p2=
     &                        bilinear_interp(
     &                             quasiset_xixi_pa,
     &                             quasiset_xixi_pb,
     &                             quasiset_xixi_pc,
     &                             quasiset_xixi_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_trace_pa=
     &                         quasiset_tracell(ipa,jpa,kpa)
                          quasiset_trace_pb=
     &                         quasiset_tracell(ipb,jpb,kpb)
                          quasiset_trace_pc=
     &                         quasiset_tracell(ipc,jpc,kpc)
                          quasiset_trace_pd=
     &                         quasiset_tracell(ipd,jpd,kpd)

                          quasiset_trace_p2=
     &                        bilinear_interp(
     &                             quasiset_trace_pa,
     &                             quasiset_trace_pb,
     &                             quasiset_trace_pc,
     &                             quasiset_trace_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_massdensity_pa=
     &                         quasiset_massdensityll(ipa,jpa,kpa)
                          quasiset_massdensity_pb=
     &                         quasiset_massdensityll(ipb,jpb,kpb)
                          quasiset_massdensity_pc=
     &                         quasiset_massdensityll(ipc,jpc,kpc)
                          quasiset_massdensity_pd=
     &                         quasiset_massdensityll(ipd,jpd,kpd)

                          quasiset_massdensity_p2=
     &                        bilinear_interp(
     &                             quasiset_massdensity_pa,
     &                             quasiset_massdensity_pb,
     &                             quasiset_massdensity_pc,
     &                             quasiset_massdensity_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_angmomdensityx_pa=
     &                         quasiset_angmomdensityxll(ipa,jpa,kpa)
                          quasiset_angmomdensityx_pb=
     &                         quasiset_angmomdensityxll(ipb,jpb,kpb)
                          quasiset_angmomdensityx_pc=
     &                         quasiset_angmomdensityxll(ipc,jpc,kpc)
                          quasiset_angmomdensityx_pd=
     &                         quasiset_angmomdensityxll(ipd,jpd,kpd)

                          quasiset_angmomdensityx_p2=
     &                        bilinear_interp(
     &                             quasiset_angmomdensityx_pa,
     &                             quasiset_angmomdensityx_pb,
     &                             quasiset_angmomdensityx_pc,
     &                             quasiset_angmomdensityx_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_angmomdensityy_pa=
     &                         quasiset_angmomdensityyll(ipa,jpa,kpa)
                          quasiset_angmomdensityy_pb=
     &                         quasiset_angmomdensityyll(ipb,jpb,kpb)
                          quasiset_angmomdensityy_pc=
     &                         quasiset_angmomdensityyll(ipc,jpc,kpc)
                          quasiset_angmomdensityy_pd=
     &                         quasiset_angmomdensityyll(ipd,jpd,kpd)

                          quasiset_angmomdensityy_p2=
     &                        bilinear_interp(
     &                             quasiset_angmomdensityy_pa,
     &                             quasiset_angmomdensityy_pb,
     &                             quasiset_angmomdensityy_pc,
     &                             quasiset_angmomdensityy_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                          quasiset_angmomdensityz_pa=
     &                         quasiset_angmomdensityzll(ipa,jpa,kpa)
                          quasiset_angmomdensityz_pb=
     &                         quasiset_angmomdensityzll(ipb,jpb,kpb)
                          quasiset_angmomdensityz_pc=
     &                         quasiset_angmomdensityzll(ipc,jpc,kpc)
                          quasiset_angmomdensityz_pd=
     &                         quasiset_angmomdensityzll(ipd,jpd,kpd)

                          quasiset_angmomdensityz_p2=
     &                        bilinear_interp(
     &                             quasiset_angmomdensityz_pa,
     &                             quasiset_angmomdensityz_pb,
     &                             quasiset_angmomdensityz_pc,
     &                             quasiset_angmomdensityz_pd,
     &                             ypa,ypb,ypc,ypd,
     &                             zpa,zpb,zpc,zpd,
     &                             yp2,
     &                             zp2)

                        end if

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
     &                        /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
     &                        /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
     &                        /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
     &                        /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
     &                        /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
     &                        /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
     &                        /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
     &                        /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
     &                        /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
     &                        /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
     &                        /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
     &                        /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
     &                        /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
     &                        /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
     &                        /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
     &                        /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           zpa,zpb,zpc,zpd,
     &                           xpa,xpb,xpc,xpd,
     &                           zp2,
     &                           xp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
                        yp2  = abs(zp2/tan(PI*chip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
                        yp2  = -abs(zp2/tan(PI*chip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
                        yp2  = abs(zp2/tan(PI*chip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
                        yp2  = -abs(zp2/tan(PI*chip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
                        yp2  = abs(zp2/tan(PI*chip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
                        yp2  = abs(zp2/tan(PI*chip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
                        yp2  = -abs(zp2/tan(PI*chip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
                        yp2  = -abs(zp2/tan(PI*chip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
                        yp2  = abs(zp2/tan(PI*chip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
                        yp2  = -abs(zp2/tan(PI*chip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
                        yp2  = abs(zp2/tan(PI*chip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
                        yp2  = -abs(zp2/tan(PI*chip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
                        yp2  = abs(zp2/tan(PI*chip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
                        yp2  = abs(zp2/tan(PI*chip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
                        yp2  = -abs(zp2/tan(PI*chip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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
                        yp2  = -abs(zp2/tan(PI*chip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        quasiset_tt_pa=
     &                        quasiset_tt_ll(ipa,jpa,kpa)
                        quasiset_tt_pb=
     &                        quasiset_tt_ll(ipb,jpb,kpb)
                        quasiset_tt_pc=
     &                        quasiset_tt_ll(ipc,jpc,kpc)
                        quasiset_tt_pd=
     &                        quasiset_tt_ll(ipd,jpd,kpd)

                        quasiset_tt_p2=
     &                      bilinear_interp(
     &                           quasiset_tt_pa,
     &                           quasiset_tt_pb,
     &                           quasiset_tt_pc,
     &                           quasiset_tt_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tt(lind)=
     &                      firstord_extrap(
     &                           quasiset_tt_p1,
     &                           quasiset_tt_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_tchi_pa=
     &                        quasiset_tchi_ll(ipa,jpa,kpa)
                        quasiset_tchi_pb=
     &                        quasiset_tchi_ll(ipb,jpb,kpb)
                        quasiset_tchi_pc=
     &                        quasiset_tchi_ll(ipc,jpc,kpc)
                        quasiset_tchi_pd=
     &                        quasiset_tchi_ll(ipd,jpd,kpd)

                        quasiset_tchi_p2=
     &                      bilinear_interp(
     &                           quasiset_tchi_pa,
     &                           quasiset_tchi_pb,
     &                           quasiset_tchi_pc,
     &                           quasiset_tchi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_tchi(lind)=
     &                      firstord_extrap(
     &                           quasiset_tchi_p1,
     &                           quasiset_tchi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_txi_pa=
     &                        quasiset_txi_ll(ipa,jpa,kpa)
                        quasiset_txi_pb=
     &                        quasiset_txi_ll(ipb,jpb,kpb)
                        quasiset_txi_pc=
     &                        quasiset_txi_ll(ipc,jpc,kpc)
                        quasiset_txi_pd=
     &                        quasiset_txi_ll(ipd,jpd,kpd)

                        quasiset_txi_p2=
     &                      bilinear_interp(
     &                           quasiset_txi_pa,
     &                           quasiset_txi_pb,
     &                           quasiset_txi_pc,
     &                           quasiset_txi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_txi(lind)=
     &                      firstord_extrap(
     &                           quasiset_txi_p1,
     &                           quasiset_txi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chichi_pa=
     &                        quasiset_chichi_ll(ipa,jpa,kpa)
                        quasiset_chichi_pb=
     &                        quasiset_chichi_ll(ipb,jpb,kpb)
                        quasiset_chichi_pc=
     &                        quasiset_chichi_ll(ipc,jpc,kpc)
                        quasiset_chichi_pd=
     &                        quasiset_chichi_ll(ipd,jpd,kpd)

                        quasiset_chichi_p2=
     &                      bilinear_interp(
     &                           quasiset_chichi_pa,
     &                           quasiset_chichi_pb,
     &                           quasiset_chichi_pc,
     &                           quasiset_chichi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chichi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chichi_p1,
     &                           quasiset_chichi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_chixi_pa=
     &                        quasiset_chixi_ll(ipa,jpa,kpa)
                        quasiset_chixi_pb=
     &                        quasiset_chixi_ll(ipb,jpb,kpb)
                        quasiset_chixi_pc=
     &                        quasiset_chixi_ll(ipc,jpc,kpc)
                        quasiset_chixi_pd=
     &                        quasiset_chixi_ll(ipd,jpd,kpd)

                        quasiset_chixi_p2=
     &                      bilinear_interp(
     &                           quasiset_chixi_pa,
     &                           quasiset_chixi_pb,
     &                           quasiset_chixi_pc,
     &                           quasiset_chixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_chixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_chixi_p1,
     &                           quasiset_chixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_xixi_pa=
     &                        quasiset_xixi_ll(ipa,jpa,kpa)
                        quasiset_xixi_pb=
     &                        quasiset_xixi_ll(ipb,jpb,kpb)
                        quasiset_xixi_pc=
     &                        quasiset_xixi_ll(ipc,jpc,kpc)
                        quasiset_xixi_pd=
     &                        quasiset_xixi_ll(ipd,jpd,kpd)

                        quasiset_xixi_p2=
     &                      bilinear_interp(
     &                           quasiset_xixi_pa,
     &                           quasiset_xixi_pb,
     &                           quasiset_xixi_pc,
     &                           quasiset_xixi_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_xixi(lind)=
     &                      firstord_extrap(
     &                           quasiset_xixi_p1,
     &                           quasiset_xixi_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_trace_pa=
     &                        quasiset_tracell(ipa,jpa,kpa)
                        quasiset_trace_pb=
     &                        quasiset_tracell(ipb,jpb,kpb)
                        quasiset_trace_pc=
     &                        quasiset_tracell(ipc,jpc,kpc)
                        quasiset_trace_pd=
     &                        quasiset_tracell(ipd,jpd,kpd)

                        quasiset_trace_p2=
     &                      bilinear_interp(
     &                           quasiset_trace_pa,
     &                           quasiset_trace_pb,
     &                           quasiset_trace_pc,
     &                           quasiset_trace_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_trace(lind)=
     &                      firstord_extrap(
     &                           quasiset_trace_p1,
     &                           quasiset_trace_p2,
     &                           rhop1,rhop2,rhoex)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )

                        quasiset_massdensity_pa=
     &                        quasiset_massdensityll(ipa,jpa,kpa)
                        quasiset_massdensity_pb=
     &                        quasiset_massdensityll(ipb,jpb,kpb)
                        quasiset_massdensity_pc=
     &                        quasiset_massdensityll(ipc,jpc,kpc)
                        quasiset_massdensity_pd=
     &                        quasiset_massdensityll(ipd,jpd,kpd)

                        quasiset_massdensity_p2=
     &                      bilinear_interp(
     &                           quasiset_massdensity_pa,
     &                           quasiset_massdensity_pb,
     &                           quasiset_massdensity_pc,
     &                           quasiset_massdensity_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_massdensity(lind)=
     &                      firstord_extrap(
     &                           quasiset_massdensity_p1,
     &                           quasiset_massdensity_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityx_pa=
     &                        quasiset_angmomdensityxll(ipa,jpa,kpa)
                        quasiset_angmomdensityx_pb=
     &                        quasiset_angmomdensityxll(ipb,jpb,kpb)
                        quasiset_angmomdensityx_pc=
     &                        quasiset_angmomdensityxll(ipc,jpc,kpc)
                        quasiset_angmomdensityx_pd=
     &                        quasiset_angmomdensityxll(ipd,jpd,kpd)

                        quasiset_angmomdensityx_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityx_pa,
     &                           quasiset_angmomdensityx_pb,
     &                           quasiset_angmomdensityx_pc,
     &                           quasiset_angmomdensityx_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityx(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityx_p1,
     &                           quasiset_angmomdensityx_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityy_pa=
     &                        quasiset_angmomdensityyll(ipa,jpa,kpa)
                        quasiset_angmomdensityy_pb=
     &                        quasiset_angmomdensityyll(ipb,jpb,kpb)
                        quasiset_angmomdensityy_pc=
     &                        quasiset_angmomdensityyll(ipc,jpc,kpc)
                        quasiset_angmomdensityy_pd=
     &                        quasiset_angmomdensityyll(ipd,jpd,kpd)

                        quasiset_angmomdensityy_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityy_pa,
     &                           quasiset_angmomdensityy_pb,
     &                           quasiset_angmomdensityy_pc,
     &                           quasiset_angmomdensityy_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityy(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityy_p1,
     &                           quasiset_angmomdensityy_p2,
     &                           rhop1,rhop2,rhoex)

                        quasiset_angmomdensityz_pa=
     &                        quasiset_angmomdensityzll(ipa,jpa,kpa)
                        quasiset_angmomdensityz_pb=
     &                        quasiset_angmomdensityzll(ipb,jpb,kpb)
                        quasiset_angmomdensityz_pc=
     &                        quasiset_angmomdensityzll(ipc,jpc,kpc)
                        quasiset_angmomdensityz_pd=
     &                        quasiset_angmomdensityzll(ipd,jpd,kpd)

                        quasiset_angmomdensityz_p2=
     &                      bilinear_interp(
     &                           quasiset_angmomdensityz_pa,
     &                           quasiset_angmomdensityz_pb,
     &                           quasiset_angmomdensityz_pc,
     &                           quasiset_angmomdensityz_pd,
     &                           xpa,xpb,xpc,xpd,
     &                           ypa,ypb,ypc,ypd,
     &                           xp2,
     &                           yp2)

                        quasiset_angmomdensityz(lind)=
     &                      firstord_extrap(
     &                           quasiset_angmomdensityz_p1,
     &                           quasiset_angmomdensityz_p2,
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