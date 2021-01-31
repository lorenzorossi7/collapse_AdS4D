c----------------------------------------------------------------------
c substitutes the scalar field profile with a Lorentz-boosted Gaussian profile 
c centred at (xu0,yu0,zu0) and with boost velocity (boost_vx,boost_vy,boost_vz)
c----------------------------------------------------------------------

        subroutine subs_boost_phi1(
     &                     f_np1,f_n,f_nm1,f_t_n,
     &                     A0,B0,C0,
     &                     r0,delta,xu0,yu0,zu0,
     &                     ex,ey,ez,
     &                     boost_vx,boost_vy,boost_vz,
     &                     L,x,y,z,dt,chr,exc,Nx,Ny,Nz)

        implicit none
        integer Nx,Ny,Nz
        real*8 f_n(Nx,Ny,Nz),f_t_n(Nx,Ny,Nz)
        real*8 f_np1(Nx,Ny,Nz),f_nm1(Nx,Ny,Nz)
        real*8 x(Nx),y(Ny),z(Nz),dt
        real*8 A0,B0,C0,r0,delta,ex,ey,ez,xu0,yu0,zu0,L
        real*8 chr(Nx,Ny,Nz),exc

        logical is_nan

        integer i,j,k
        integer stype
        integer a,b
        real*8 r,x0,y0,z0,rho0,xi0,chi0,csr,xb,yb,zb
        real*8 boost_vx,boost_vy,boost_vz
        real*8 boost_vnorm
        real*8 gamma
        real*8 lambda_boost(4,4),lambdainv_boost(4,4)
        real*8 t0,t0_invboost
        real*8 x0_invboost,y0_invboost,z0_invboost
        real*8 xb_invboost,yb_invboost,zb_invboost
        real*8 r_invboost
        real*8 x0_invboost_t,y0_invboost_t,z0_invboost_t
        real*8 r_invboost_t
        real*8 derexp0

        real*8 rhoc,rhod
        real*8 f1,trans

        real*8 PI
        parameter (PI=3.141592653589793d0)

        ! initialize fixed-size variables
        data i,j,k/0,0,0/
        data r,x0,y0,z0,rho0,xi0,chi0,csr,xb,yb,zb
     &       /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
        data lambda_boost,lambdainv_boost/16*0.0,16*0.0/

        !--------------------------------------------------------------

       if ((abs(boost_vx).gt.10.0d0**(-10)).or.
     &     (abs(boost_vy).gt.10.0d0**(-10)).or.
     &     (abs(boost_vz).gt.10.0d0**(-10)) ) then
        do i=1,Nx
           do j=1,Ny
            do k=1,Nz
             if (chr(i,j,k).ne.exc) then
              f_n(i,j,k)=0
              f_np1(i,j,k)=0
              f_nm1(i,j,k)=0
              f_t_n(i,j,k)=0
              x0=x(i)
              y0=y(j)
              z0=z(k)
              rho0=sqrt(x0**2+y0**2+z0**2)
              if (rho0.ne.0.0d0) then
               xi0=acos(x0/rho0)
              end if
              if ((y0.ne.0.0d0).or.(z0.ne.0.0d0)) then
                chi0=atan2(z0,y0)
                if (chi0.lt.0) chi0=chi0+2*PI
              end if
!            write (*,*) 'DEBUG from gauss3d'
!            write(*,*) 'rho0=',rho0
!            write(*,*) 'chi0,xi0=',chi0,xi0


                boost_vnorm=sqrt(boost_vx**2+boost_vy**2+boost_vz**2)
                gamma=1/sqrt(1-boost_vnorm**2)

                lambda_boost(1,1)=gamma
                lambda_boost(1,2)=-gamma*boost_vx
                lambda_boost(1,3)=-gamma*boost_vy
                lambda_boost(1,4)=-gamma*boost_vz
                lambda_boost(2,2)=1+(gamma-1)*(boost_vx**2)
     &           /boost_vnorm**2
                lambda_boost(2,3)=(gamma-1)*boost_vx*boost_vy
     &           /boost_vnorm**2
                lambda_boost(2,4)=(gamma-1)*boost_vx*boost_vz
     &           /boost_vnorm**2
                lambda_boost(3,3)=1+(gamma-1)*(boost_vy**2)
     &           /boost_vnorm**2
                lambda_boost(3,4)=(gamma-1)*boost_vy*boost_vz
     &           /boost_vnorm**2
                lambda_boost(4,4)=1+(gamma-1)*(boost_vz**2)
     &           /boost_vnorm**2
                !the matrix of Lorentz boosts is symmetric
                do a=1,4
                 do b=a+1,4
                  lambda_boost(b,a)=lambda_boost(a,b)
                 end do
                end do

                lambdainv_boost(1,1)=-1/((-1+boost_vnorm**2)*gamma)
                lambdainv_boost(1,2)=-boost_vx
     &           /((-1+boost_vnorm**2)*gamma)
                lambdainv_boost(1,3)=-boost_vy
     &           /((-1+boost_vnorm**2)*gamma)
                lambdainv_boost(1,4)=-boost_vz
     &           /((-1+boost_vnorm**2)*gamma)
                lambdainv_boost(2,2)=(1/boost_vnorm**2)*
     &           (boost_vy**2+boost_vz**2-(boost_vx**2)
     &           /((-1+boost_vnorm**2)*gamma))
                lambdainv_boost(2,3)=-(boost_vx*boost_vy*
     &           (1+(-1+boost_vnorm**2)*gamma))
     &           /((-1+boost_vnorm**2)*boost_vnorm**2*gamma)
                lambdainv_boost(2,4)=-(boost_vx*boost_vz*
     &           (1+(-1+boost_vnorm**2)*gamma))
     &           /((-1+boost_vnorm**2)*boost_vnorm**2*gamma)
                lambdainv_boost(3,3)=(1/boost_vnorm**2)*
     &           (boost_vx**2+boost_vz**2-(boost_vy**2)
     &           /((-1+boost_vnorm**2)*gamma))
                lambdainv_boost(3,4)=-(boost_vy*boost_vz*
     &           (1+(-1+boost_vnorm**2)*gamma))
     &           /((-1+boost_vnorm**2)*boost_vnorm**2*gamma)
                lambdainv_boost(4,4)=(1/boost_vnorm**2)*
     &           (boost_vx**2+boost_vy**2-(boost_vz**2)
     &           /((-1+boost_vnorm**2)*gamma))
                !the matrix of Lorentz boosts is symmetric
                do a=1,4
                 do b=a+1,4
                  lambdainv_boost(b,a)=lambdainv_boost(a,b)
                 end do
                end do

                !The inverse Lorentz boost is lambdainv_boost^(-1)^\mu_\nu x^\mu in UNCOMPACTIFIED CARTESIAN COORDINATES
                !We compute this in compactified Cartesian coords at t=0
                t0=0
                !computing t0_invboost is not necessary but we do it for completeness

                t0_invboost=
     -  lambdainv_boost(1,1)*t0 - 
     -  (2*(lambdainv_boost(1,2)*x0 + lambdainv_boost(1,3)*y0 + 
     -       lambdainv_boost(1,4)*z0))/
     -   (-1 + x0**2 + y0**2 + z0**2)
                x0_invboost=
     -         ((lambdainv_boost(1,2)*t0 - 
     -      (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -           lambdainv_boost(2,4)*z0))/
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    (-1 + Sqrt(1 + 
     -        (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*t0 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*t0 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((lambdainv_boost(1,2)*t0 - 
     -       (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -            lambdainv_boost(2,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2 + 
     -    (lambdainv_boost(1,3)*t0 - 
     -       (2*(lambdainv_boost(2,3)*x0 + lambdainv_boost(3,3)*y0 + 
     -            lambdainv_boost(3,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2 + 
     -    (lambdainv_boost(1,4)*t0 - 
     -       (2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -            lambdainv_boost(4,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2)
                y0_invboost=
     -          ((-2*(lambdainv_boost(2,3)*x0 
     -      + lambdainv_boost(3,3)*y0 + 
     -         lambdainv_boost(3,4)*z0) + 
     -      lambdainv_boost(1,3)*t0*
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    Sqrt(1 - (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + 
     -               lambdainv_boost(2,4)*z0))/
     -           (-1 + x0**2 + y0**2 + z0**2))**2/
     -       ((lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -         ))*(-1 + 
     -      Sqrt(1 + 
     -        (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*t0 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*t0 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((-1 + x0**2 + y0**2 + z0**2)*
     -    Sqrt((lambdainv_boost(1,3)*t0 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*t0 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2)*
     -    Sqrt((lambdainv_boost(1,2)*t0 - 
     -         (2*(lambdainv_boost(2,2)*x0 + 
     -              lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,3)*t0 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*t0 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2))
                z0_invboost=
     -  ((-2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -         lambdainv_boost(4,4)*z0) + 
     -      lambdainv_boost(1,4)*t0*
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    Sqrt(1 - (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + 
     -               lambdainv_boost(2,4)*z0))/
     -           (-1 + x0**2 + y0**2 + z0**2))**2/
     -       ((lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -         ))*(-1 + 
     -      Sqrt(1 + 
     -        (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*t0 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*t0 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((-1 + x0**2 + y0**2 + z0**2)*
     -    Sqrt((lambdainv_boost(1,3)*t0 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*t0 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2)*
     -    Sqrt((lambdainv_boost(1,2)*t0 - 
     -         (2*(lambdainv_boost(2,2)*x0 + 
     -              lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,3)*t0 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*t0 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2))

              if ((abs(y0).lt.10.0d0**(-10))
     -         .and.(abs(z0).lt.10.0d0**(-10)) ) then
                y0_invboost=0
                z0_invboost=0
              end if

              if ((abs(x0).lt.10.0d0**(-10))
     -         .and.(abs(z0).lt.10.0d0**(-10)) ) then
                x0_invboost=0
                z0_invboost=0
              end if

              if ((abs(x0).lt.10.0d0**(-10))
     -         .and.(abs(y0).lt.10.0d0**(-10)) ) then
                x0_invboost=0
                y0_invboost=0
              end if

              if ((abs(x0).lt.10.0d0**(-10))
     -         .and.(abs(y0).lt.10.0d0**(-10))
     -          .and.(abs(z0).lt.10.0d0**(-10)) ) then
                x0_invboost=0
                y0_invboost=0
                z0_invboost=0
              end if




               xb_invboost=x0_invboost-xu0
               yb_invboost=y0_invboost-yu0
               zb_invboost=z0_invboost-zu0

               r_invboost=sqrt(xb_invboost**2+yb_invboost**2
     &          +zb_invboost**2
     &           -ex**2*xb_invboost**2-ey**2*yb_invboost**2
     &           -ez**2*zb_invboost**2)

               f_n(i,j,k)=
     &             A0*exp(-((r_invboost-r0)/delta)**2)





















               !set derivative of f(i,j,k) w.r.t. t at t=0

               x0_invboost_t=
     -         ((lambdainv_boost(1,2)*t0 - 
     -       (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -            lambdainv_boost(2,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (2*lambdainv_boost(1,2)*
     -        (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,3)*
     -        (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,4)*
     -        (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )))/
     -   (2.*((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt(1 + (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)) - 
     -  ((lambdainv_boost(1,2)*t0 - 
     -       (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -            lambdainv_boost(2,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (2*lambdainv_boost(1,2)*
     -        (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,3)*
     -        (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,4)*
     -        (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   ((lambdainv_boost(1,2)*t0 - 
     -         (2*(lambdainv_boost(2,2)*x0 + 
     -              lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0)
     -            )/(-1 + x0**2 + y0**2 + z0**2))**2
     -        + (lambdainv_boost(1,3)*t0 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0)
     -            )/(-1 + x0**2 + y0**2 + z0**2))**2
     -        + (lambdainv_boost(1,4)*t0 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0)
     -            )/(-1 + x0**2 + y0**2 + z0**2))**2
     -      )**2 + (lambdainv_boost(1,2)*
     -     (-1 + Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   ((lambdainv_boost(1,2)*t0 - 
     -        (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -             lambdainv_boost(2,4)*z0))/
     -         (-1 + x0**2 + y0**2 + z0**2))**2 + 
     -     (lambdainv_boost(1,3)*t0 - 
     -        (2*(lambdainv_boost(2,3)*x0 + lambdainv_boost(3,3)*y0 + 
     -             lambdainv_boost(3,4)*z0))/
     -         (-1 + x0**2 + y0**2 + z0**2))**2 + 
     -     (lambdainv_boost(1,4)*t0 - 
     -        (2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -             lambdainv_boost(4,4)*z0))/
     -         (-1 + x0**2 + y0**2 + z0**2))**2)

               y0_invboost_t=
     -         ((-2*(lambdainv_boost(2,3)*x0 + lambdainv_boost(3,3)*y0 + 
     -          lambdainv_boost(3,4)*z0) + 
     -       lambdainv_boost(1,3)*t0*
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (2*lambdainv_boost(1,2)*
     -        (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,3)*
     -        (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,4)*
     -        (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ))*Sqrt(1 - 
     -       (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          )))/
     -   (2.*(-1 + x0**2 + y0**2 + z0**2)*
     -     Sqrt((lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt(1 + (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)) + 
     -  ((-2*(lambdainv_boost(2,3)*x0 + lambdainv_boost(3,3)*y0 + 
     -          lambdainv_boost(3,4)*z0) + 
     -       lambdainv_boost(1,3)*t0*
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -           *(2*lambdainv_boost(1,2)*
     -             (lambdainv_boost(1,2)*t0 - 
     -               (2*
     -                  (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -                (-1 + x0**2 + y0**2 + z0**2))
     -             + 2*lambdainv_boost(1,3)*
     -             (lambdainv_boost(1,3)*t0 - 
     -               (2*
     -                  (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -                (-1 + x0**2 + y0**2 + z0**2))
     -             + 2*lambdainv_boost(1,4)*
     -             (lambdainv_boost(1,4)*t0 - 
     -               (2*
     -                  (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -                (-1 + x0**2 + y0**2 + z0**2)))
     -          )/
     -        ((lambdainv_boost(1,2)*t0 - 
     -              (2*
     -                 (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -               (-1 + x0**2 + y0**2 + z0**2))**
     -            2 + 
     -           (lambdainv_boost(1,3)*t0 - 
     -              (2*
     -                 (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -               (-1 + x0**2 + y0**2 + z0**2))**
     -            2 + 
     -           (lambdainv_boost(1,4)*t0 - 
     -              (2*
     -                 (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -               (-1 + x0**2 + y0**2 + z0**2))**
     -            2)**2 - 
     -       (2*lambdainv_boost(1,2)*
     -          (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2)))/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   (2.*(-1 + x0**2 + y0**2 + z0**2)*
     -     Sqrt((lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt(1 - (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))) - 
     -  ((-2*(lambdainv_boost(2,3)*x0 + lambdainv_boost(3,3)*y0 + 
     -          lambdainv_boost(3,4)*z0) + 
     -       lambdainv_boost(1,3)*t0*
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (2*lambdainv_boost(1,2)*
     -        (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,3)*
     -        (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,4)*
     -        (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ))*Sqrt(1 - 
     -       (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   (2.*(-1 + x0**2 + y0**2 + z0**2)*
     -     Sqrt((lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     ((lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*t0 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*t0 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -       **1.5) - 
     -  ((-2*(lambdainv_boost(2,3)*x0 + lambdainv_boost(3,3)*y0 + 
     -          lambdainv_boost(3,4)*z0) + 
     -       lambdainv_boost(1,3)*t0*
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (2*lambdainv_boost(1,3)*
     -        (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,4)*
     -        (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ))*Sqrt(1 - 
     -       (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   (2.*(-1 + x0**2 + y0**2 + z0**2)*
     -     ((lambdainv_boost(1,3)*t0 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*t0 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -       **1.5*Sqrt((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)) + 
     -  (lambdainv_boost(1,3)*
     -     Sqrt(1 - (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   (Sqrt((lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2))

               z0_invboost_t=
     -         ((-2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -          lambdainv_boost(4,4)*z0) + 
     -       lambdainv_boost(1,4)*t0*
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (2*lambdainv_boost(1,2)*
     -        (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,3)*
     -        (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,4)*
     -        (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ))*Sqrt(1 - 
     -       (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          )))/
     -   (2.*(-1 + x0**2 + y0**2 + z0**2)*
     -     Sqrt((lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt(1 + (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)) + 
     -  ((-2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -          lambdainv_boost(4,4)*z0) + 
     -       lambdainv_boost(1,4)*t0*
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -           *(2*lambdainv_boost(1,2)*
     -             (lambdainv_boost(1,2)*t0 - 
     -               (2*
     -                  (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -                (-1 + x0**2 + y0**2 + z0**2))
     -             + 2*lambdainv_boost(1,3)*
     -             (lambdainv_boost(1,3)*t0 - 
     -               (2*
     -                  (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -                (-1 + x0**2 + y0**2 + z0**2))
     -             + 2*lambdainv_boost(1,4)*
     -             (lambdainv_boost(1,4)*t0 - 
     -               (2*
     -                  (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -                (-1 + x0**2 + y0**2 + z0**2)))
     -          )/
     -        ((lambdainv_boost(1,2)*t0 - 
     -              (2*
     -                 (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -               (-1 + x0**2 + y0**2 + z0**2))**
     -            2 + 
     -           (lambdainv_boost(1,3)*t0 - 
     -              (2*
     -                 (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -               (-1 + x0**2 + y0**2 + z0**2))**
     -            2 + 
     -           (lambdainv_boost(1,4)*t0 - 
     -              (2*
     -                 (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -               (-1 + x0**2 + y0**2 + z0**2))**
     -            2)**2 - 
     -       (2*lambdainv_boost(1,2)*
     -          (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2)))/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   (2.*(-1 + x0**2 + y0**2 + z0**2)*
     -     Sqrt((lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt(1 - (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))) - 
     -  ((-2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -          lambdainv_boost(4,4)*z0) + 
     -       lambdainv_boost(1,4)*t0*
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (2*lambdainv_boost(1,2)*
     -        (lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,3)*
     -        (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,4)*
     -        (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ))*Sqrt(1 - 
     -       (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   (2.*(-1 + x0**2 + y0**2 + z0**2)*
     -     Sqrt((lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     ((lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*t0 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*t0 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -       **1.5) - 
     -  ((-2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -          lambdainv_boost(4,4)*z0) + 
     -       lambdainv_boost(1,4)*t0*
     -        (-1 + x0**2 + y0**2 + z0**2))*
     -     (2*lambdainv_boost(1,3)*
     -        (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ) + 2*lambdainv_boost(1,4)*
     -        (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          ))*Sqrt(1 - 
     -       (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   (2.*(-1 + x0**2 + y0**2 + z0**2)*
     -     ((lambdainv_boost(1,3)*t0 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*t0 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -       **1.5*Sqrt((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)) + 
     -  (lambdainv_boost(1,4)*
     -     Sqrt(1 - (lambdainv_boost(1,2)*t0 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2/
     -        ((lambdainv_boost(1,2)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,2)*x0 + 
     -                  lambdainv_boost(2,3)*y0 + 
     -                  lambdainv_boost(2,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,3)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,3)*x0 + 
     -                  lambdainv_boost(3,3)*y0 + 
     -                  lambdainv_boost(3,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -            + (lambdainv_boost(1,4)*t0 - 
     -             (2*
     -                (lambdainv_boost(2,4)*x0 + 
     -                  lambdainv_boost(3,4)*y0 + 
     -                  lambdainv_boost(4,4)*z0))/
     -              (-1 + x0**2 + y0**2 + z0**2))**2
     -          ))*(-1 + 
     -       Sqrt(1 + 
     -         (lambdainv_boost(1,2)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,3)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -          + (lambdainv_boost(1,4)*t0 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2)
     -       ))/
     -   (Sqrt((lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2)*
     -     Sqrt((lambdainv_boost(1,2)*t0 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,3)*t0 - 
     -          (2*(lambdainv_boost(2,3)*x0 + 
     -               lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2 + 
     -       (lambdainv_boost(1,4)*t0 - 
     -          (2*(lambdainv_boost(2,4)*x0 + 
     -               lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -               ))/(-1 + x0**2 + y0**2 + z0**2)
     -          )**2))


              if ((abs(y0).lt.10.0d0**(-10))
     -         .and.(abs(z0).lt.10.0d0**(-10)) ) then
                y0_invboost_t=0
                z0_invboost_t=0
              end if

              if ((abs(x0).lt.10.0d0**(-10))
     -         .and.(abs(z0).lt.10.0d0**(-10)) ) then
                x0_invboost_t=0
                z0_invboost_t=0
              end if

              if ((abs(x0).lt.10.0d0**(-10))
     -         .and.(abs(y0).lt.10.0d0**(-10)) ) then
                x0_invboost_t=0
                y0_invboost_t=0
              end if

              if ((abs(x0).lt.10.0d0**(-10))
     -         .and.(abs(y0).lt.10.0d0**(-10))
     -          .and.(abs(z0).lt.10.0d0**(-10)) ) then
                x0_invboost_t=0
                y0_invboost_t=0
                z0_invboost_t=0
              end if


              r_invboost_t=(1/(2*r_invboost))*
     &             ((1-ex**2)*2*xb_invboost*x0_invboost_t+
     &              (1-ey**2)*2*yb_invboost*y0_invboost_t+
     &              (1-ez**2)*2*zb_invboost*z0_invboost_t)

              derexp0=(-1/delta**2)*2*(r_invboost-r0)*r_invboost_t

               f_t_n(i,j,k)=f_n(i,j,k)*derexp0

               if (rho0.lt.10.0d0**(-10)) then
                f_t_n(i,j,k)=0
               end if
               f_nm1(i,j,k)=f_n(i,j,k) - f_t_n(i,j,k)*dt
               f_np1(i,j,k)=f_n(i,j,k) + f_t_n(i,j,k)*dt

!!DEBUG
!              if ((is_nan(f_n(i,j,k)))
!     &         .or.(is_nan(f_t_n(i,j,k)))) then
!                write (*,*) "x0,y0,z0=",x0,y0,z0
!                write (*,*) "xu0,yu0,zu0=",xu0,yu0,zu0
!                write (*,*) "A0,r0,delta=",A0,r0,delta
!               do a=1,4
!                do b=1,4
!                  write (*,*) "a,b,lambda_boost(a,b)="
!     &             ,a,b,lambda_boost(a,b)
!                  write (*,*) "a,b,lambdainv_boost(a,b)="
!     &             ,a,b,lambdainv_boost(a,b)
!                end do
!               end do
!                  write (*,*) "x0_invboost=",x0_invboost
!                  write (*,*) "y0_invboost=",y0_invboost
!                  write (*,*) "z0_invboost=",z0_invboost
!                  write (*,*) "r_invboost=",r_invboost
!                  write (*,*) "x0_invboost_t=",x0_invboost_t
!                  write (*,*) "y0_invboost_t=",y0_invboost_t
!                  write (*,*) "z0_invboost_t=",z0_invboost_t
!                  write (*,*) "r_invboost_t=",r_invboost_t
!                  write (*,*) "f_n(i,j,k)=",f_n(i,j,k)
!                  write (*,*) "f_t_n(i,j,k)=",f_t_n(i,j,k)
!              end if


!              if ((abs(x0-0.375d0).lt.10.0d0**(-10)).and.
!     &            (abs(y0-0.125d0).lt.10.0d0**(-10)).and.
!     &            (abs(z0-0.0d0).lt.10.0d0**(-10))) then
!                write (*,*) "x0,y0,z0=",x0,y0,z0
!                write (*,*) "xu0,yu0,zu0=",xu0,yu0,zu0
!                write (*,*) "A0,r0,delta=",A0,r0,delta
!               do a=1,4
!                do b=1,4
!                  write (*,*) "a,b,lambda_boost(a,b)="
!     &             ,a,b,lambda_boost(a,b)
!                  write (*,*) "a,b,lambdainv_boost(a,b)="
!     &             ,a,b,lambdainv_boost(a,b)
!                end do
!               end do
!                  write (*,*) "xb_invboost=",xb_invboost
!                  write (*,*) "yb_invboost=",yb_invboost
!                  write (*,*) "zb_invboost=",zb_invboost
!                  write (*,*) "r_invboost=",r_invboost
!                  write (*,*) "f(i,j,k)=",f(i,j,k)
!                  write (*,*) "f_t(i,j,k)=",f_t(i,j,k)
!              end if

             end if
            end do
           end do
        end do
       end if

        return
        end