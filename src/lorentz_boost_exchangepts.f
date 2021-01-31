c----------------------------------------------------------------------
c Lorentz-boost phi1_n,phi1_np1,phi1_nm1
c boost velocity (boost_vx,boost_vy,boost_vz)
c----------------------------------------------------------------------

        subroutine boost_phi1(
     &                     phi1_np1,phi1_n,phi1_nm1,
     &                     boost_vx,boost_vy,boost_vz,
     &                     L,x,y,z,dt,chr,ex,Nx,Ny,Nz)

        implicit none
        integer Nx,Ny,Nz
        real*8 f_n(Nx,Ny,Nz),f_t_n(Nx,Ny,Nz)
        real*8 f_np1(Nx,Ny,Nz),f_nm1(Nx,Ny,Nz)
        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz)
        real*8 L
        real dx,dy,dz,dt

        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)

        real*8 phi1_atinvboostpt0

        real*8 phi1_atinvboostptnp1

        real*8 phi1_atinvboostptnm1

        integer i,j,k
        integer stype
        integer a,b
        real*8 r,x0,y0,z0,rho0,xi0,chi0,csr,xb,yb,zb
        real*8 xunc0,yunc0,zunc0
        real*8 boost_vx,boost_vy,boost_vz
        real*8 boost_vnorm
        real*8 gamma
        real*8 lambda_boost(4,4),lambdainv_boost(4,4)
        real*8 t0,tnp1,tnm1
        real*8 t0_invboost
        real*8 x0_invboost,y0_invboost,z0_invboost
        real*8 tnp1_invboost
        real*8 xnp1_invboost,ynp1_invboost,znp1_invboost
        real*8 tnm1_invboost
        real*8 xnm1_invboost,ynm1_invboost,znm1_invboost
        integer i0_invboost,j0_invboost,k0_invboost
        integer inp1_invboost,jnp1_invboost,knp1_invboost
        integer inm1_invboost,jnm1_invboost,knm1_invboost
        real*8 fx0,fy0,fz0
        real*8 fxnp1,fynp1,fznp1
        real*8 fxnm1,fynm1,fznm1

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

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)
 
       if ((abs(boost_vx).gt.10.0d0**(-10)).or.
     &     (abs(boost_vy).gt.10.0d0**(-10)).or.
     &     (abs(boost_vz).gt.10.0d0**(-10)) ) then
        do i=1,Nx
           do j=1,Ny
            do k=1,Nz
             if (chr(i,j,k).ne.ex) then
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
                !We compute this in compactified Cartesian coords at t=0,dt,-dt
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

              if ((y0.lt.10.0d0**(-10))
     -         .and.(z0.lt.10.0d0**(-10)) ) then
                y0_invboost=0
                z0_invboost=0
                if (x0.lt.10.0d0**(-10)) then
                  x0_invboost=0
                end if
              end if

              !find the (i,j,k) indices of the grid point (x0_invboost,y0_invboost,z0_invboost)
              i0_invboost=(x0_invboost-x(1))/dx+1
              j0_invboost=(y0_invboost-y(1))/dy+1
              k0_invboost=(z0_invboost-z(1))/dz+1

              if ((i0_invboost.gt.(Nx-1)).or.
     -            (j0_invboost.gt.(Ny-1)).or.
     -            (k0_invboost.gt.(Nz-1))) then
                write (*,*) "i0_invboost=",i0_invboost
                write (*,*) "j0_invboost=",j0_invboost
                write (*,*) "k0_invboost=",k0_invboost
                write (*,*) "x0_invboost=",x0_invboost
                write (*,*) "y0_invboost=",y0_invboost
                write (*,*) "z0_invboost=",z0_invboost
                write (*,*) "i=",i
                write (*,*) "j=",j
                write (*,*) "k=",k
                write (*,*) "Nx,Ny,Nz=",Nx,Ny,Nz
                write (*,*) "x0=",x0
                write (*,*) "y0=",y0
                write (*,*) "z0=",z0
                stop

              end if
!
!              if ((i0_invboost.lt.(1)).or.
!     -            (j0_invboost.lt.(1)).or.
!     -            (k0_invboost.lt.(1))) then
!                write (*,*) "i0_invboost=",i0_invboost
!                write (*,*) "j0_invboost=",j0_invboost
!                write (*,*) "k0_invboost=",k0_invboost
!                write (*,*) "x0_invboost=",x0_invboost
!                write (*,*) "y0_invboost=",y0_invboost
!                write (*,*) "z0_invboost=",z0_invboost
!                write (*,*) "i=",i
!                write (*,*) "j=",j
!                write (*,*) "k=",k
!                write (*,*) "x0=",x0
!                write (*,*) "y0=",y0
!                write (*,*) "z0=",z0
!                stop
!
!              end if

              ! for bilinear interpolation of theta from 
              ! (i0_invboost,j0_invboost,k0_invboost),(i0_invboost+1,j0_invboost,k0_invboost),(i0_invboost,j0_invboost+1,k0_invboost),(i0_invboost+1,j0_invboost+1,k0_invboost),(i0_invboost,j0_invboost,k0_invboost+1),(i0_invboost+1,j0_invboost,k0_invboost+1),(i0_invboost,j0_invboost+1,k0_invboost+1),(i0_invboost+1,j0_invboost+1,k0_invboost+1)
              !(NOTE: since, x0_invboost,y0_invboost,z0_invboost does not necessarily lie 
              ! on a grid point these fx,fy,fz are *not* identically zero)
              fx0=(x0_invboost-((i0_invboost-1)*dx+x(1)))/dx
              fy0=(y0_invboost-((j0_invboost-1)*dy+y(1)))/dy
              fz0=(z0_invboost-((k0_invboost-1)*dz+z(1)))/dz

              !from bilinear interpolation from neighbouring points
              !obtain the value of gb at (x0_invboost,y0_invboost,z0_invboost)
              phi1_atinvboostpt0=
     &           (1-fx0)*(1-fy0)*(1-fz0)*
     &              phi1_n(i0_invboost,j0_invboost,k0_invboost)+
     &           (  fx0)*(1-fy0)*(1-fz0)*
     &              phi1_n(i0_invboost+1,j0_invboost,k0_invboost)+
     &           (1-fx0)*(  fy0)*(1-fz0)*
     &              phi1_n(i0_invboost,j0_invboost+1,k0_invboost)+
     &           (1-fx0)*(1-fy0)*(  fz0)*
     &              phi1_n(i0_invboost,j0_invboost,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(1-fz0)*
     &              phi1_n(i0_invboost+1,j0_invboost+1,k0_invboost)+
     &           (  fx0)*(1-fy0)*(  fz0)*
     &              phi1_n(i0_invboost+1,j0_invboost,k0_invboost+1)+
     &           (1-fx0)*(  fy0)*(  fz0)*
     &              phi1_n(i0_invboost,j0_invboost+1,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(  fz0)*
     &              phi1_n(i0_invboost+1,j0_invboost+1,k0_invboost+1)

               phi1_n(i,j,k)=phi1_atinvboostpt0

               !set derivative of phi1_np1 and phi1_nm1
               !for this we do what we did above but for t=dt and -dt
                tnp1=dt
                tnm1=-dt
               !computing tnp1_invboost and tnm1_invboost is not necessary but we do it for completeness

                tnp1_invboost=
     -  lambdainv_boost(1,1)*tnp1 - 
     -  (2*(lambdainv_boost(1,2)*x0 + lambdainv_boost(1,3)*y0 + 
     -       lambdainv_boost(1,4)*z0))/
     -   (-1 + x0**2 + y0**2 + z0**2)
                xnp1_invboost=
     -         ((lambdainv_boost(1,2)*tnp1 - 
     -      (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -           lambdainv_boost(2,4)*z0))/
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    (-1 + Sqrt(1 + 
     -        (lambdainv_boost(1,2)*tnp1 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*tnp1 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*tnp1 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((lambdainv_boost(1,2)*tnp1 - 
     -       (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -            lambdainv_boost(2,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2 + 
     -    (lambdainv_boost(1,3)*tnp1 - 
     -       (2*(lambdainv_boost(2,3)*x0 + lambdainv_boost(3,3)*y0 + 
     -            lambdainv_boost(3,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2 + 
     -    (lambdainv_boost(1,4)*tnp1 - 
     -       (2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -            lambdainv_boost(4,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2)
                ynp1_invboost=
     -          ((-2*(lambdainv_boost(2,3)*x0 
     -      + lambdainv_boost(3,3)*y0 + 
     -         lambdainv_boost(3,4)*z0) + 
     -      lambdainv_boost(1,3)*tnp1*
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    Sqrt(1 - (lambdainv_boost(1,2)*tnp1 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + 
     -               lambdainv_boost(2,4)*z0))/
     -           (-1 + x0**2 + y0**2 + z0**2))**2/
     -       ((lambdainv_boost(1,2)*tnp1 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,3)*tnp1 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,4)*tnp1 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -         ))*(-1 + 
     -      Sqrt(1 + 
     -        (lambdainv_boost(1,2)*tnp1 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*tnp1 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*tnp1 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((-1 + x0**2 + y0**2 + z0**2)*
     -    Sqrt((lambdainv_boost(1,3)*tnp1 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*tnp1 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2)*
     -    Sqrt((lambdainv_boost(1,2)*tnp1 - 
     -         (2*(lambdainv_boost(2,2)*x0 + 
     -              lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,3)*tnp1 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*tnp1 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2))
                znp1_invboost=
     -  ((-2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -         lambdainv_boost(4,4)*z0) + 
     -      lambdainv_boost(1,4)*tnp1*
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    Sqrt(1 - (lambdainv_boost(1,2)*tnp1 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + 
     -               lambdainv_boost(2,4)*z0))/
     -           (-1 + x0**2 + y0**2 + z0**2))**2/
     -       ((lambdainv_boost(1,2)*tnp1 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,3)*tnp1 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,4)*tnp1 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -         ))*(-1 + 
     -      Sqrt(1 + 
     -        (lambdainv_boost(1,2)*tnp1 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*tnp1 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*tnp1 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((-1 + x0**2 + y0**2 + z0**2)*
     -    Sqrt((lambdainv_boost(1,3)*tnp1 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*tnp1 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2)*
     -    Sqrt((lambdainv_boost(1,2)*tnp1 - 
     -         (2*(lambdainv_boost(2,2)*x0 + 
     -              lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,3)*tnp1 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*tnp1 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2))


                tnm1_invboost=
     -  lambdainv_boost(1,1)*tnm1 - 
     -  (2*(lambdainv_boost(1,2)*x0 + lambdainv_boost(1,3)*y0 + 
     -       lambdainv_boost(1,4)*z0))/
     -   (-1 + x0**2 + y0**2 + z0**2)
                xnm1_invboost=
     -         ((lambdainv_boost(1,2)*tnm1 - 
     -      (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -           lambdainv_boost(2,4)*z0))/
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    (-1 + Sqrt(1 + 
     -        (lambdainv_boost(1,2)*tnm1 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*tnm1 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*tnm1 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((lambdainv_boost(1,2)*tnm1 - 
     -       (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -            lambdainv_boost(2,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2 + 
     -    (lambdainv_boost(1,3)*tnm1 - 
     -       (2*(lambdainv_boost(2,3)*x0 + lambdainv_boost(3,3)*y0 + 
     -            lambdainv_boost(3,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2 + 
     -    (lambdainv_boost(1,4)*tnm1 - 
     -       (2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -            lambdainv_boost(4,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2)
                ynm1_invboost=
     -          ((-2*(lambdainv_boost(2,3)*x0 
     -      + lambdainv_boost(3,3)*y0 + 
     -         lambdainv_boost(3,4)*z0) + 
     -      lambdainv_boost(1,3)*tnm1*
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    Sqrt(1 - (lambdainv_boost(1,2)*tnm1 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + 
     -               lambdainv_boost(2,4)*z0))/
     -           (-1 + x0**2 + y0**2 + z0**2))**2/
     -       ((lambdainv_boost(1,2)*tnm1 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,3)*tnm1 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,4)*tnm1 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -         ))*(-1 + 
     -      Sqrt(1 + 
     -        (lambdainv_boost(1,2)*tnm1 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*tnm1 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*tnm1 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((-1 + x0**2 + y0**2 + z0**2)*
     -    Sqrt((lambdainv_boost(1,3)*tnm1 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*tnm1 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2)*
     -    Sqrt((lambdainv_boost(1,2)*tnm1 - 
     -         (2*(lambdainv_boost(2,2)*x0 + 
     -              lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,3)*tnm1 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*tnm1 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2))
                znm1_invboost=
     -  ((-2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -         lambdainv_boost(4,4)*z0) + 
     -      lambdainv_boost(1,4)*tnm1*
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    Sqrt(1 - (lambdainv_boost(1,2)*tnm1 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + 
     -               lambdainv_boost(2,4)*z0))/
     -           (-1 + x0**2 + y0**2 + z0**2))**2/
     -       ((lambdainv_boost(1,2)*tnm1 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,3)*tnm1 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,4)*tnm1 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -         ))*(-1 + 
     -      Sqrt(1 + 
     -        (lambdainv_boost(1,2)*tnm1 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*tnm1 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*tnm1 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((-1 + x0**2 + y0**2 + z0**2)*
     -    Sqrt((lambdainv_boost(1,3)*tnm1 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*tnm1 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2)*
     -    Sqrt((lambdainv_boost(1,2)*tnm1 - 
     -         (2*(lambdainv_boost(2,2)*x0 + 
     -              lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,3)*tnm1 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*tnm1 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2))

              if ((y0.lt.10.0d0**(-10))
     -         .and.(z0.lt.10.0d0**(-10)) ) then
                ynp1_invboost=0
                znp1_invboost=0
                ynm1_invboost=0
                znm1_invboost=0
              end if


              inp1_invboost=(xnp1_invboost-x(1))/dx+1
              jnp1_invboost=(ynp1_invboost-y(1))/dy+1
              knp1_invboost=(znp1_invboost-z(1))/dz+1

              inm1_invboost=(xnm1_invboost-x(1))/dx+1
              jnm1_invboost=(ynm1_invboost-y(1))/dy+1
              knm1_invboost=(znm1_invboost-z(1))/dz+1

              fxnp1=(xnp1_invboost-((inp1_invboost-1)*dx+x(1)))/dx
              fynp1=(ynp1_invboost-((jnp1_invboost-1)*dy+y(1)))/dy
              fznp1=(znp1_invboost-((knp1_invboost-1)*dz+z(1)))/dz

              fxnm1=(xnm1_invboost-((inm1_invboost-1)*dx+x(1)))/dx
              fynm1=(ynm1_invboost-((jnm1_invboost-1)*dy+y(1)))/dy
              fznm1=(znm1_invboost-((knm1_invboost-1)*dz+z(1)))/dz

              phi1_atinvboostptnp1=
     &             (1-fxnp1)*(1-fynp1)*(1-fznp1)*
     &         phi1_np1(inp1_invboost,jnp1_invboost,knp1_invboost)+
     &             (  fxnp1)*(1-fynp1)*(1-fznp1)*
     &         phi1_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost)+
     &             (1-fxnp1)*(  fynp1)*(1-fznp1)*
     &         phi1_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost)+
     &             (1-fxnp1)*(1-fynp1)*(  fznp1)*
     &         phi1_np1(inp1_invboost,jnp1_invboost,knp1_invboost+1)+
     &             (  fxnp1)*(  fynp1)*(1-fznp1)*
     &         phi1_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost)+
     &             (  fxnp1)*(1-fynp1)*(  fznp1)*
     &         phi1_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost+1)+
     &             (1-fxnp1)*(  fynp1)*(  fznp1)*
     &         phi1_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost+1)+
     &             (  fxnp1)*(  fynp1)*(  fznp1)*
     &         phi1_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost+1)




              phi1_atinvboostptnm1=
     &             (1-fxnm1)*(1-fynm1)*(1-fznm1)*
     &         phi1_nm1(inm1_invboost,jnm1_invboost,knm1_invboost)+
     &             (  fxnm1)*(1-fynm1)*(1-fznm1)*
     &         phi1_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost)+
     &             (1-fxnm1)*(  fynm1)*(1-fznm1)*
     &         phi1_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost)+
     &             (1-fxnm1)*(1-fynm1)*(  fznm1)*
     &         phi1_nm1(inm1_invboost,jnm1_invboost,knm1_invboost+1)+
     &             (  fxnm1)*(  fynm1)*(1-fznm1)*
     &         phi1_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost)+
     &             (  fxnm1)*(1-fynm1)*(  fznm1)*
     &         phi1_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost+1)+
     &             (1-fxnm1)*(  fynm1)*(  fznm1)*
     &         phi1_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost+1)+
     &             (  fxnm1)*(  fynm1)*(  fznm1)*
     &         phi1_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost+1)


               phi1_np1(i,j,k)=phi1_atinvboostptnp1

               phi1_nm1(i,j,k)=phi1_atinvboostptnm1


             end if
            end do
           end do
        end do
       end if

        return
        end

c----------------------------------------------------------------------
c Lorentz-boost gb_n,gb_np1,gb_nm1
c boost velocity (boost_vx,boost_vy,boost_vz)
c----------------------------------------------------------------------

        subroutine boost_gb(
     &                     gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                     gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                     gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                     gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &                     gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                     gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                     gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &                     gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                     gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &                     gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &                     boost_vx,boost_vy,boost_vz,
     &                     L,x,y,z,dt,chr,ex,Nx,Ny,Nz)

        implicit none
        integer Nx,Ny,Nz
        real*8 f_n(Nx,Ny,Nz),f_t_n(Nx,Ny,Nz)
        real*8 f_np1(Nx,Ny,Nz),f_nm1(Nx,Ny,Nz)
        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz)
        real*8 L
        real dx,dy,dz,dt

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
        real*8 Hb_t_np1(Nx,Ny,Nz),Hb_t_n(Nx,Ny,Nz),Hb_t_nm1(Nx,Ny,Nz)
        real*8 Hb_x_np1(Nx,Ny,Nz),Hb_x_n(Nx,Ny,Nz),Hb_x_nm1(Nx,Ny,Nz)
        real*8 Hb_y_np1(Nx,Ny,Nz),Hb_y_n(Nx,Ny,Nz),Hb_y_nm1(Nx,Ny,Nz)
        real*8 Hb_z_np1(Nx,Ny,Nz),Hb_z_n(Nx,Ny,Nz),Hb_z_nm1(Nx,Ny,Nz)
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)

        real*8 gb_tt_atinvboostpt0
        real*8 gb_tx_atinvboostpt0
        real*8 gb_ty_atinvboostpt0
        real*8 gb_tz_atinvboostpt0
        real*8 gb_xx_atinvboostpt0
        real*8 gb_xy_atinvboostpt0
        real*8 gb_xz_atinvboostpt0
        real*8 gb_yy_atinvboostpt0
        real*8 gb_yz_atinvboostpt0
        real*8 gb_zz_atinvboostpt0

        real*8 gb_tt_atinvboostptnp1
        real*8 gb_tx_atinvboostptnp1
        real*8 gb_ty_atinvboostptnp1
        real*8 gb_tz_atinvboostptnp1
        real*8 gb_xx_atinvboostptnp1
        real*8 gb_xy_atinvboostptnp1
        real*8 gb_xz_atinvboostptnp1
        real*8 gb_yy_atinvboostptnp1
        real*8 gb_yz_atinvboostptnp1
        real*8 gb_zz_atinvboostptnp1

        real*8 gb_tt_atinvboostptnm1
        real*8 gb_tx_atinvboostptnm1
        real*8 gb_ty_atinvboostptnm1
        real*8 gb_tz_atinvboostptnm1
        real*8 gb_xx_atinvboostptnm1
        real*8 gb_xy_atinvboostptnm1
        real*8 gb_xz_atinvboostptnm1
        real*8 gb_yy_atinvboostptnm1
        real*8 gb_yz_atinvboostptnm1
        real*8 gb_zz_atinvboostptnm1

        integer i,j,k
        integer stype
        integer a,b
        real*8 r,x0,y0,z0,rho0,xi0,chi0,csr,xb,yb,zb
        real*8 boost_vx,boost_vy,boost_vz
        real*8 boost_vnorm
        real*8 gamma
        real*8 lambda_boost(4,4),lambdainv_boost(4,4)
        real*8 t0,tnp1,tnm1
        real*8 t0_invboost
        real*8 x0_invboost,y0_invboost,z0_invboost
        real*8 tnp1_invboost
        real*8 xnp1_invboost,ynp1_invboost,znp1_invboost
        real*8 tnm1_invboost
        real*8 xnm1_invboost,ynm1_invboost,znm1_invboost
        integer i0_invboost,j0_invboost,k0_invboost
        integer inp1_invboost,jnp1_invboost,knp1_invboost
        integer inm1_invboost,jnm1_invboost,knm1_invboost
        real*8 fx0,fy0,fz0
        real*8 fxnp1,fynp1,fznp1
        real*8 fxnm1,fynm1,fznm1

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

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)
 
       if ((abs(boost_vx).gt.10.0d0**(-10)).or.
     &     (abs(boost_vy).gt.10.0d0**(-10)).or.
     &     (abs(boost_vz).gt.10.0d0**(-10)) ) then
        do i=1,Nx
           do j=1,Ny
            do k=1,Nz
             if (chr(i,j,k).ne.ex) then
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

                !lambdainv_boost^(-1)^\mu_\nu x^\mu at t=0,dt,-dt
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

              if ((y0.lt.10.0d0**(-10))
     -         .and.(z0.lt.10.0d0**(-10)) ) then
                y0_invboost=0
                z0_invboost=0
                if (x0.lt.10.0d0**(-10)) then
                  x0_invboost=0
                end if
              end if

              !find the (i,j,k) indices of the grid point (x0_invboost,y0_invboost,z0_invboost)
              i0_invboost=(x0_invboost-x(1))/dx+1
              j0_invboost=(y0_invboost-y(1))/dy+1
              k0_invboost=(z0_invboost-z(1))/dz+1

              ! for bilinear interpolation of theta from 
              ! (i,j,k),(i+1,j,k),(i,j+1,k),(i+1,j+1,k),(i,j,k+1),(i+1,j,k+1),(i,j+1,k+1),(i+1,j+1,k+1)
              !(NOTE: since, x0_invboost,y0_invboost,z0_invboost does not necessarily lie 
              ! on a grid point these fx,fy,fz are *not* identically zero)
              fx0=(x0_invboost-((i0_invboost-1)*dx+x(1)))/dx
              fy0=(y0_invboost-((j0_invboost-1)*dy+y(1)))/dy
              fz0=(z0_invboost-((k0_invboost-1)*dz+z(1)))/dz

              !from bilinear interpolation from neighbouring points
              !obtain the value of gb at (x0_invboost,y0_invboost,z0_invboost)
              gb_tt_atinvboostpt0=
     &           (1-fx0)*(1-fy0)*(1-fz0)*
     &              gb_tt_n(i0_invboost,j0_invboost,k0_invboost)+
     &           (  fx0)*(1-fy0)*(1-fz0)*
     &              gb_tt_n(i0_invboost+1,j0_invboost,k0_invboost)+
     &           (1-fx0)*(  fy0)*(1-fz0)*
     &              gb_tt_n(i0_invboost,j0_invboost+1,k0_invboost)+
     &           (1-fx0)*(1-fy0)*(  fz0)*
     &              gb_tt_n(i0_invboost,j0_invboost,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(1-fz0)*
     &              gb_tt_n(i0_invboost+1,j0_invboost+1,k0_invboost)+
     &           (  fx0)*(1-fy0)*(  fz0)*
     &              gb_tt_n(i0_invboost+1,j0_invboost,k0_invboost+1)+
     &           (1-fx0)*(  fy0)*(  fz0)*
     &              gb_tt_n(i0_invboost,j0_invboost+1,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(  fz0)*
     &              gb_tt_n(i0_invboost+1,j0_invboost+1,k0_invboost+1)

              gb_tx_atinvboostpt0=
     &           (1-fx0)*(1-fy0)*(1-fz0)*
     &              gb_tx_n(i0_invboost,j0_invboost,k0_invboost)+
     &           (  fx0)*(1-fy0)*(1-fz0)*
     &              gb_tx_n(i0_invboost+1,j0_invboost,k0_invboost)+
     &           (1-fx0)*(  fy0)*(1-fz0)*
     &              gb_tx_n(i0_invboost,j0_invboost+1,k0_invboost)+
     &           (1-fx0)*(1-fy0)*(  fz0)*
     &              gb_tx_n(i0_invboost,j0_invboost,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(1-fz0)*
     &              gb_tx_n(i0_invboost+1,j0_invboost+1,k0_invboost)+
     &           (  fx0)*(1-fy0)*(  fz0)*
     &              gb_tx_n(i0_invboost+1,j0_invboost,k0_invboost+1)+
     &           (1-fx0)*(  fy0)*(  fz0)*
     &              gb_tx_n(i0_invboost,j0_invboost+1,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(  fz0)*
     &              gb_tx_n(i0_invboost+1,j0_invboost+1,k0_invboost+1)

              gb_ty_atinvboostpt0=
     &           (1-fx0)*(1-fy0)*(1-fz0)*
     &              gb_ty_n(i0_invboost,j0_invboost,k0_invboost)+
     &           (  fx0)*(1-fy0)*(1-fz0)*
     &              gb_ty_n(i0_invboost+1,j0_invboost,k0_invboost)+
     &           (1-fx0)*(  fy0)*(1-fz0)*
     &              gb_ty_n(i0_invboost,j0_invboost+1,k0_invboost)+
     &           (1-fx0)*(1-fy0)*(  fz0)*
     &              gb_ty_n(i0_invboost,j0_invboost,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(1-fz0)*
     &              gb_ty_n(i0_invboost+1,j0_invboost+1,k0_invboost)+
     &           (  fx0)*(1-fy0)*(  fz0)*
     &              gb_ty_n(i0_invboost+1,j0_invboost,k0_invboost+1)+
     &           (1-fx0)*(  fy0)*(  fz0)*
     &              gb_ty_n(i0_invboost,j0_invboost+1,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(  fz0)*
     &              gb_ty_n(i0_invboost+1,j0_invboost+1,k0_invboost+1)

              gb_tz_atinvboostpt0=
     &           (1-fx0)*(1-fy0)*(1-fz0)*
     &              gb_tz_n(i0_invboost,j0_invboost,k0_invboost)+
     &           (  fx0)*(1-fy0)*(1-fz0)*
     &              gb_tz_n(i0_invboost+1,j0_invboost,k0_invboost)+
     &           (1-fx0)*(  fy0)*(1-fz0)*
     &              gb_tz_n(i0_invboost,j0_invboost+1,k0_invboost)+
     &           (1-fx0)*(1-fy0)*(  fz0)*
     &              gb_tz_n(i0_invboost,j0_invboost,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(1-fz0)*
     &              gb_tz_n(i0_invboost+1,j0_invboost+1,k0_invboost)+
     &           (  fx0)*(1-fy0)*(  fz0)*
     &              gb_tz_n(i0_invboost+1,j0_invboost,k0_invboost+1)+
     &           (1-fx0)*(  fy0)*(  fz0)*
     &              gb_tz_n(i0_invboost,j0_invboost+1,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(  fz0)*
     &              gb_tz_n(i0_invboost+1,j0_invboost+1,k0_invboost+1)

              gb_xx_atinvboostpt0=
     &           (1-fx0)*(1-fy0)*(1-fz0)*
     &              gb_xx_n(i0_invboost,j0_invboost,k0_invboost)+
     &           (  fx0)*(1-fy0)*(1-fz0)*
     &              gb_xx_n(i0_invboost+1,j0_invboost,k0_invboost)+
     &           (1-fx0)*(  fy0)*(1-fz0)*
     &              gb_xx_n(i0_invboost,j0_invboost+1,k0_invboost)+
     &           (1-fx0)*(1-fy0)*(  fz0)*
     &              gb_xx_n(i0_invboost,j0_invboost,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(1-fz0)*
     &              gb_xx_n(i0_invboost+1,j0_invboost+1,k0_invboost)+
     &           (  fx0)*(1-fy0)*(  fz0)*
     &              gb_xx_n(i0_invboost+1,j0_invboost,k0_invboost+1)+
     &           (1-fx0)*(  fy0)*(  fz0)*
     &              gb_xx_n(i0_invboost,j0_invboost+1,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(  fz0)*
     &              gb_xx_n(i0_invboost+1,j0_invboost+1,k0_invboost+1)

              gb_xy_atinvboostpt0=
     &           (1-fx0)*(1-fy0)*(1-fz0)*
     &              gb_xy_n(i0_invboost,j0_invboost,k0_invboost)+
     &           (  fx0)*(1-fy0)*(1-fz0)*
     &              gb_xy_n(i0_invboost+1,j0_invboost,k0_invboost)+
     &           (1-fx0)*(  fy0)*(1-fz0)*
     &              gb_xy_n(i0_invboost,j0_invboost+1,k0_invboost)+
     &           (1-fx0)*(1-fy0)*(  fz0)*
     &              gb_xy_n(i0_invboost,j0_invboost,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(1-fz0)*
     &              gb_xy_n(i0_invboost+1,j0_invboost+1,k0_invboost)+
     &           (  fx0)*(1-fy0)*(  fz0)*
     &              gb_xy_n(i0_invboost+1,j0_invboost,k0_invboost+1)+
     &           (1-fx0)*(  fy0)*(  fz0)*
     &              gb_xy_n(i0_invboost,j0_invboost+1,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(  fz0)*
     &              gb_xy_n(i0_invboost+1,j0_invboost+1,k0_invboost+1)

              gb_xz_atinvboostpt0=
     &           (1-fx0)*(1-fy0)*(1-fz0)*
     &              gb_xz_n(i0_invboost,j0_invboost,k0_invboost)+
     &           (  fx0)*(1-fy0)*(1-fz0)*
     &              gb_xz_n(i0_invboost+1,j0_invboost,k0_invboost)+
     &           (1-fx0)*(  fy0)*(1-fz0)*
     &              gb_xz_n(i0_invboost,j0_invboost+1,k0_invboost)+
     &           (1-fx0)*(1-fy0)*(  fz0)*
     &              gb_xz_n(i0_invboost,j0_invboost,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(1-fz0)*
     &              gb_xz_n(i0_invboost+1,j0_invboost+1,k0_invboost)+
     &           (  fx0)*(1-fy0)*(  fz0)*
     &              gb_xz_n(i0_invboost+1,j0_invboost,k0_invboost+1)+
     &           (1-fx0)*(  fy0)*(  fz0)*
     &              gb_xz_n(i0_invboost,j0_invboost+1,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(  fz0)*
     &              gb_xz_n(i0_invboost+1,j0_invboost+1,k0_invboost+1)

              gb_yy_atinvboostpt0=
     &           (1-fx0)*(1-fy0)*(1-fz0)*
     &              gb_yy_n(i0_invboost,j0_invboost,k0_invboost)+
     &           (  fx0)*(1-fy0)*(1-fz0)*
     &              gb_yy_n(i0_invboost+1,j0_invboost,k0_invboost)+
     &           (1-fx0)*(  fy0)*(1-fz0)*
     &              gb_yy_n(i0_invboost,j0_invboost+1,k0_invboost)+
     &           (1-fx0)*(1-fy0)*(  fz0)*
     &              gb_yy_n(i0_invboost,j0_invboost,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(1-fz0)*
     &              gb_yy_n(i0_invboost+1,j0_invboost+1,k0_invboost)+
     &           (  fx0)*(1-fy0)*(  fz0)*
     &              gb_yy_n(i0_invboost+1,j0_invboost,k0_invboost+1)+
     &           (1-fx0)*(  fy0)*(  fz0)*
     &              gb_yy_n(i0_invboost,j0_invboost+1,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(  fz0)*
     &              gb_yy_n(i0_invboost+1,j0_invboost+1,k0_invboost+1)

              gb_yz_atinvboostpt0=
     &           (1-fx0)*(1-fy0)*(1-fz0)*
     &              gb_yz_n(i0_invboost,j0_invboost,k0_invboost)+
     &           (  fx0)*(1-fy0)*(1-fz0)*
     &              gb_yz_n(i0_invboost+1,j0_invboost,k0_invboost)+
     &           (1-fx0)*(  fy0)*(1-fz0)*
     &              gb_yz_n(i0_invboost,j0_invboost+1,k0_invboost)+
     &           (1-fx0)*(1-fy0)*(  fz0)*
     &              gb_yz_n(i0_invboost,j0_invboost,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(1-fz0)*
     &              gb_yz_n(i0_invboost+1,j0_invboost+1,k0_invboost)+
     &           (  fx0)*(1-fy0)*(  fz0)*
     &              gb_yz_n(i0_invboost+1,j0_invboost,k0_invboost+1)+
     &           (1-fx0)*(  fy0)*(  fz0)*
     &              gb_yz_n(i0_invboost,j0_invboost+1,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(  fz0)*
     &              gb_yz_n(i0_invboost+1,j0_invboost+1,k0_invboost+1)

              gb_zz_atinvboostpt0=
     &           (1-fx0)*(1-fy0)*(1-fz0)*
     &              gb_zz_n(i0_invboost,j0_invboost,k0_invboost)+
     &           (  fx0)*(1-fy0)*(1-fz0)*
     &              gb_zz_n(i0_invboost+1,j0_invboost,k0_invboost)+
     &           (1-fx0)*(  fy0)*(1-fz0)*
     &              gb_zz_n(i0_invboost,j0_invboost+1,k0_invboost)+
     &           (1-fx0)*(1-fy0)*(  fz0)*
     &              gb_zz_n(i0_invboost,j0_invboost,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(1-fz0)*
     &              gb_zz_n(i0_invboost+1,j0_invboost+1,k0_invboost)+
     &           (  fx0)*(1-fy0)*(  fz0)*
     &              gb_zz_n(i0_invboost+1,j0_invboost,k0_invboost+1)+
     &           (1-fx0)*(  fy0)*(  fz0)*
     &              gb_zz_n(i0_invboost,j0_invboost+1,k0_invboost+1)+
     &           (  fx0)*(  fy0)*(  fz0)*
     &              gb_zz_n(i0_invboost+1,j0_invboost+1,k0_invboost+1)




               gb_tt_n(i,j,k)=gb_tt_atinvboostpt0
               gb_tx_n(i,j,k)=gb_tx_atinvboostpt0
               gb_ty_n(i,j,k)=gb_ty_atinvboostpt0
               gb_tz_n(i,j,k)=gb_tz_atinvboostpt0
               gb_xx_n(i,j,k)=gb_xx_atinvboostpt0
               gb_xy_n(i,j,k)=gb_xy_atinvboostpt0
               gb_xz_n(i,j,k)=gb_xz_atinvboostpt0
               gb_yy_n(i,j,k)=gb_yy_atinvboostpt0
               gb_yz_n(i,j,k)=gb_yz_atinvboostpt0
               gb_zz_n(i,j,k)=gb_zz_atinvboostpt0

               !set derivative of gb_np1 and gb_nm1
               !for this we do what we did above but for t=dt and -dt
                tnp1=dt
                tnm1=-dt
               !computing tnp1_invboost and tnm1_invboost is not necessary but we do it for completeness
                tnp1_invboost=
     -  lambdainv_boost(1,1)*tnp1 - 
     -  (2*(lambdainv_boost(1,2)*x0 + lambdainv_boost(1,3)*y0 + 
     -       lambdainv_boost(1,4)*z0))/
     -   (-1 + x0**2 + y0**2 + z0**2)
                xnp1_invboost=
     -         ((lambdainv_boost(1,2)*tnp1 - 
     -      (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -           lambdainv_boost(2,4)*z0))/
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    (-1 + Sqrt(1 + 
     -        (lambdainv_boost(1,2)*tnp1 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*tnp1 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*tnp1 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((lambdainv_boost(1,2)*tnp1 - 
     -       (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -            lambdainv_boost(2,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2 + 
     -    (lambdainv_boost(1,3)*tnp1 - 
     -       (2*(lambdainv_boost(2,3)*x0 + lambdainv_boost(3,3)*y0 + 
     -            lambdainv_boost(3,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2 + 
     -    (lambdainv_boost(1,4)*tnp1 - 
     -       (2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -            lambdainv_boost(4,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2)
                ynp1_invboost=
     -          ((-2*(lambdainv_boost(2,3)*x0 
     -      + lambdainv_boost(3,3)*y0 + 
     -         lambdainv_boost(3,4)*z0) + 
     -      lambdainv_boost(1,3)*tnp1*
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    Sqrt(1 - (lambdainv_boost(1,2)*tnp1 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + 
     -               lambdainv_boost(2,4)*z0))/
     -           (-1 + x0**2 + y0**2 + z0**2))**2/
     -       ((lambdainv_boost(1,2)*tnp1 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,3)*tnp1 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,4)*tnp1 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -         ))*(-1 + 
     -      Sqrt(1 + 
     -        (lambdainv_boost(1,2)*tnp1 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*tnp1 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*tnp1 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((-1 + x0**2 + y0**2 + z0**2)*
     -    Sqrt((lambdainv_boost(1,3)*tnp1 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*tnp1 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2)*
     -    Sqrt((lambdainv_boost(1,2)*tnp1 - 
     -         (2*(lambdainv_boost(2,2)*x0 + 
     -              lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,3)*tnp1 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*tnp1 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2))
                znp1_invboost=
     -  ((-2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -         lambdainv_boost(4,4)*z0) + 
     -      lambdainv_boost(1,4)*tnp1*
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    Sqrt(1 - (lambdainv_boost(1,2)*tnp1 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + 
     -               lambdainv_boost(2,4)*z0))/
     -           (-1 + x0**2 + y0**2 + z0**2))**2/
     -       ((lambdainv_boost(1,2)*tnp1 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,3)*tnp1 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,4)*tnp1 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -         ))*(-1 + 
     -      Sqrt(1 + 
     -        (lambdainv_boost(1,2)*tnp1 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*tnp1 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*tnp1 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((-1 + x0**2 + y0**2 + z0**2)*
     -    Sqrt((lambdainv_boost(1,3)*tnp1 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*tnp1 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2)*
     -    Sqrt((lambdainv_boost(1,2)*tnp1 - 
     -         (2*(lambdainv_boost(2,2)*x0 + 
     -              lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,3)*tnp1 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*tnp1 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2))


                tnm1_invboost=
     -  lambdainv_boost(1,1)*tnm1 - 
     -  (2*(lambdainv_boost(1,2)*x0 + lambdainv_boost(1,3)*y0 + 
     -       lambdainv_boost(1,4)*z0))/
     -   (-1 + x0**2 + y0**2 + z0**2)
                xnm1_invboost=
     -         ((lambdainv_boost(1,2)*tnm1 - 
     -      (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -           lambdainv_boost(2,4)*z0))/
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    (-1 + Sqrt(1 + 
     -        (lambdainv_boost(1,2)*tnm1 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*tnm1 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*tnm1 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((lambdainv_boost(1,2)*tnm1 - 
     -       (2*(lambdainv_boost(2,2)*x0 + lambdainv_boost(2,3)*y0 + 
     -            lambdainv_boost(2,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2 + 
     -    (lambdainv_boost(1,3)*tnm1 - 
     -       (2*(lambdainv_boost(2,3)*x0 + lambdainv_boost(3,3)*y0 + 
     -            lambdainv_boost(3,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2 + 
     -    (lambdainv_boost(1,4)*tnm1 - 
     -       (2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -            lambdainv_boost(4,4)*z0))/
     -        (-1 + x0**2 + y0**2 + z0**2))**2)
                ynm1_invboost=
     -          ((-2*(lambdainv_boost(2,3)*x0 
     -      + lambdainv_boost(3,3)*y0 + 
     -         lambdainv_boost(3,4)*z0) + 
     -      lambdainv_boost(1,3)*tnm1*
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    Sqrt(1 - (lambdainv_boost(1,2)*tnm1 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + 
     -               lambdainv_boost(2,4)*z0))/
     -           (-1 + x0**2 + y0**2 + z0**2))**2/
     -       ((lambdainv_boost(1,2)*tnm1 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,3)*tnm1 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,4)*tnm1 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -         ))*(-1 + 
     -      Sqrt(1 + 
     -        (lambdainv_boost(1,2)*tnm1 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*tnm1 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*tnm1 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((-1 + x0**2 + y0**2 + z0**2)*
     -    Sqrt((lambdainv_boost(1,3)*tnm1 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*tnm1 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2)*
     -    Sqrt((lambdainv_boost(1,2)*tnm1 - 
     -         (2*(lambdainv_boost(2,2)*x0 + 
     -              lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,3)*tnm1 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*tnm1 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2))
                znm1_invboost=
     -  ((-2*(lambdainv_boost(2,4)*x0 + lambdainv_boost(3,4)*y0 + 
     -         lambdainv_boost(4,4)*z0) + 
     -      lambdainv_boost(1,4)*tnm1*
     -       (-1 + x0**2 + y0**2 + z0**2))*
     -    Sqrt(1 - (lambdainv_boost(1,2)*tnm1 - 
     -          (2*(lambdainv_boost(2,2)*x0 + 
     -               lambdainv_boost(2,3)*y0 + 
     -               lambdainv_boost(2,4)*z0))/
     -           (-1 + x0**2 + y0**2 + z0**2))**2/
     -       ((lambdainv_boost(1,2)*tnm1 - 
     -            (2*
     -               (lambdainv_boost(2,2)*x0 + 
     -                 lambdainv_boost(2,3)*y0 + 
     -                 lambdainv_boost(2,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,3)*tnm1 - 
     -            (2*
     -               (lambdainv_boost(2,3)*x0 + 
     -                 lambdainv_boost(3,3)*y0 + 
     -                 lambdainv_boost(3,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -           + (lambdainv_boost(1,4)*tnm1 - 
     -            (2*
     -               (lambdainv_boost(2,4)*x0 + 
     -                 lambdainv_boost(3,4)*y0 + 
     -                 lambdainv_boost(4,4)*z0))/
     -             (-1 + x0**2 + y0**2 + z0**2))**2
     -         ))*(-1 + 
     -      Sqrt(1 + 
     -        (lambdainv_boost(1,2)*tnm1 - 
     -           (2*(lambdainv_boost(2,2)*x0 + 
     -                lambdainv_boost(2,3)*y0 + 
     -                lambdainv_boost(2,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,3)*tnm1 - 
     -           (2*(lambdainv_boost(2,3)*x0 + 
     -                lambdainv_boost(3,3)*y0 + 
     -                lambdainv_boost(3,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2
     -         + (lambdainv_boost(1,4)*tnm1 - 
     -           (2*(lambdainv_boost(2,4)*x0 + 
     -                lambdainv_boost(3,4)*y0 + 
     -                lambdainv_boost(4,4)*z0))/
     -            (-1 + x0**2 + y0**2 + z0**2))**2)
     -      ))/
     -  ((-1 + x0**2 + y0**2 + z0**2)*
     -    Sqrt((lambdainv_boost(1,3)*tnm1 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*tnm1 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2)*
     -    Sqrt((lambdainv_boost(1,2)*tnm1 - 
     -         (2*(lambdainv_boost(2,2)*x0 + 
     -              lambdainv_boost(2,3)*y0 + lambdainv_boost(2,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,3)*tnm1 - 
     -         (2*(lambdainv_boost(2,3)*x0 + 
     -              lambdainv_boost(3,3)*y0 + lambdainv_boost(3,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2 + 
     -      (lambdainv_boost(1,4)*tnm1 - 
     -         (2*(lambdainv_boost(2,4)*x0 + 
     -              lambdainv_boost(3,4)*y0 + lambdainv_boost(4,4)*z0
     -              ))/(-1 + x0**2 + y0**2 + z0**2)
     -         )**2))

              if ((y0.lt.10.0d0**(-10))
     -         .and.(z0.lt.10.0d0**(-10)) ) then
                ynp1_invboost=0
                znp1_invboost=0
                ynm1_invboost=0
                znm1_invboost=0
              end if

              inp1_invboost=(xnp1_invboost-x(1))/dx+1
              jnp1_invboost=(ynp1_invboost-y(1))/dy+1
              knp1_invboost=(znp1_invboost-z(1))/dz+1

              inm1_invboost=(xnm1_invboost-x(1))/dx+1
              jnm1_invboost=(ynm1_invboost-y(1))/dy+1
              knm1_invboost=(znm1_invboost-z(1))/dz+1

              fxnp1=(xnp1_invboost-((inp1_invboost-1)*dx+x(1)))/dx
              fynp1=(ynp1_invboost-((jnp1_invboost-1)*dy+y(1)))/dy
              fznp1=(znp1_invboost-((knp1_invboost-1)*dz+z(1)))/dz

              fxnm1=(xnm1_invboost-((inm1_invboost-1)*dx+x(1)))/dx
              fynm1=(ynm1_invboost-((jnm1_invboost-1)*dy+y(1)))/dy
              fznm1=(znm1_invboost-((knm1_invboost-1)*dz+z(1)))/dz

              gb_tt_atinvboostptnp1=
     &           (1-fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_tt_np1(inp1_invboost,jnp1_invboost,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_tt_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost)+
     &           (1-fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_tt_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost)+
     &           (1-fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_tt_np1(inp1_invboost,jnp1_invboost,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_tt_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_tt_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost+1)+
     &           (1-fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_tt_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_tt_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost+1)

              gb_tx_atinvboostptnp1=
     &           (1-fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_tx_np1(inp1_invboost,jnp1_invboost,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_tx_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost)+
     &           (1-fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_tx_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost)+
     &           (1-fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_tx_np1(inp1_invboost,jnp1_invboost,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_tx_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_tx_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost+1)+
     &           (1-fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_tx_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_tx_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost+1)

              gb_ty_atinvboostptnp1=
     &           (1-fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_ty_np1(inp1_invboost,jnp1_invboost,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_ty_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost)+
     &           (1-fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_ty_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost)+
     &           (1-fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_ty_np1(inp1_invboost,jnp1_invboost,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_ty_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_ty_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost+1)+
     &           (1-fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_ty_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_ty_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost+1)

              gb_tz_atinvboostptnp1=
     &           (1-fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_tz_np1(inp1_invboost,jnp1_invboost,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_tz_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost)+
     &           (1-fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_tz_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost)+
     &           (1-fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_tz_np1(inp1_invboost,jnp1_invboost,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_tz_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_tz_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost+1)+
     &           (1-fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_tz_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_tz_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost+1)

              gb_xx_atinvboostptnp1=
     &           (1-fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_xx_np1(inp1_invboost,jnp1_invboost,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_xx_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost)+
     &           (1-fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_xx_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost)+
     &           (1-fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_xx_np1(inp1_invboost,jnp1_invboost,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_xx_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_xx_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost+1)+
     &           (1-fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_xx_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_xx_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost+1)

              gb_xy_atinvboostptnp1=
     &           (1-fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_xy_np1(inp1_invboost,jnp1_invboost,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_xy_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost)+
     &           (1-fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_xy_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost)+
     &           (1-fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_xy_np1(inp1_invboost,jnp1_invboost,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_xy_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_xy_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost+1)+
     &           (1-fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_xy_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_xy_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost+1)

              gb_xz_atinvboostptnp1=
     &           (1-fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_xz_np1(inp1_invboost,jnp1_invboost,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_xz_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost)+
     &           (1-fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_xz_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost)+
     &           (1-fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_xz_np1(inp1_invboost,jnp1_invboost,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_xz_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_xz_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost+1)+
     &           (1-fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_xz_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_xz_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost+1)

              gb_yy_atinvboostptnp1=
     &           (1-fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_yy_np1(inp1_invboost,jnp1_invboost,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_yy_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost)+
     &           (1-fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_yy_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost)+
     &           (1-fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_yy_np1(inp1_invboost,jnp1_invboost,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_yy_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_yy_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost+1)+
     &           (1-fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_yy_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_yy_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost+1)

              gb_yz_atinvboostptnp1=
     &           (1-fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_yz_np1(inp1_invboost,jnp1_invboost,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_yz_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost)+
     &           (1-fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_yz_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost)+
     &           (1-fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_yz_np1(inp1_invboost,jnp1_invboost,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_yz_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_yz_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost+1)+
     &           (1-fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_yz_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_yz_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost+1)

              gb_zz_atinvboostptnp1=
     &           (1-fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_zz_np1(inp1_invboost,jnp1_invboost,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(1-fznp1)*
     &      gb_zz_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost)+
     &           (1-fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_zz_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost)+
     &           (1-fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_zz_np1(inp1_invboost,jnp1_invboost,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(1-fznp1)*
     &      gb_zz_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost)+
     &           (  fxnp1)*(1-fynp1)*(  fznp1)*
     &      gb_zz_np1(inp1_invboost+1,jnp1_invboost,knp1_invboost+1)+
     &           (1-fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_zz_np1(inp1_invboost,jnp1_invboost+1,knp1_invboost+1)+
     &           (  fxnp1)*(  fynp1)*(  fznp1)*
     &      gb_zz_np1(inp1_invboost+1,jnp1_invboost+1,knp1_invboost+1)





              gb_tt_atinvboostptnm1=
     &           (1-fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_tt_nm1(inm1_invboost,jnm1_invboost,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_tt_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost)+
     &           (1-fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_tt_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost)+
     &           (1-fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_tt_nm1(inm1_invboost,jnm1_invboost,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_tt_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_tt_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost+1)+
     &           (1-fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_tt_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_tt_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost+1)

              gb_tx_atinvboostptnm1=
     &           (1-fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_tx_nm1(inm1_invboost,jnm1_invboost,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_tx_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost)+
     &           (1-fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_tx_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost)+
     &           (1-fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_tx_nm1(inm1_invboost,jnm1_invboost,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_tx_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_tx_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost+1)+
     &           (1-fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_tx_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_tx_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost+1)

              gb_ty_atinvboostptnm1=
     &           (1-fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_ty_nm1(inm1_invboost,jnm1_invboost,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_ty_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost)+
     &           (1-fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_ty_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost)+
     &           (1-fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_ty_nm1(inm1_invboost,jnm1_invboost,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_ty_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_ty_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost+1)+
     &           (1-fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_ty_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_ty_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost+1)

              gb_tz_atinvboostptnm1=
     &           (1-fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_tz_nm1(inm1_invboost,jnm1_invboost,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_tz_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost)+
     &           (1-fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_tz_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost)+
     &           (1-fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_tz_nm1(inm1_invboost,jnm1_invboost,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_tz_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_tz_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost+1)+
     &           (1-fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_tz_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_tz_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost+1)

              gb_xx_atinvboostptnm1=
     &           (1-fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_xx_nm1(inm1_invboost,jnm1_invboost,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_xx_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost)+
     &           (1-fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_xx_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost)+
     &           (1-fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_xx_nm1(inm1_invboost,jnm1_invboost,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_xx_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_xx_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost+1)+
     &           (1-fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_xx_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_xx_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost+1)

              gb_xy_atinvboostptnm1=
     &           (1-fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_xy_nm1(inm1_invboost,jnm1_invboost,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_xy_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost)+
     &           (1-fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_xy_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost)+
     &           (1-fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_xy_nm1(inm1_invboost,jnm1_invboost,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_xy_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_xy_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost+1)+
     &           (1-fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_xy_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_xy_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost+1)

              gb_xz_atinvboostptnm1=
     &           (1-fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_xz_nm1(inm1_invboost,jnm1_invboost,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_xz_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost)+
     &           (1-fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_xz_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost)+
     &           (1-fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_xz_nm1(inm1_invboost,jnm1_invboost,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_xz_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_xz_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost+1)+
     &           (1-fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_xz_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_xz_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost+1)

              gb_yy_atinvboostptnm1=
     &           (1-fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_yy_nm1(inm1_invboost,jnm1_invboost,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_yy_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost)+
     &           (1-fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_yy_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost)+
     &           (1-fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_yy_nm1(inm1_invboost,jnm1_invboost,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_yy_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_yy_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost+1)+
     &           (1-fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_yy_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_yy_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost+1)

              gb_yz_atinvboostptnm1=
     &           (1-fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_yz_nm1(inm1_invboost,jnm1_invboost,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_yz_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost)+
     &           (1-fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_yz_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost)+
     &           (1-fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_yz_nm1(inm1_invboost,jnm1_invboost,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_yz_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_yz_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost+1)+
     &           (1-fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_yz_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_yz_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost+1)

              gb_zz_atinvboostptnm1=
     &           (1-fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_zz_nm1(inm1_invboost,jnm1_invboost,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(1-fznm1)*
     &      gb_zz_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost)+
     &           (1-fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_zz_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost)+
     &           (1-fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_zz_nm1(inm1_invboost,jnm1_invboost,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(1-fznm1)*
     &      gb_zz_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost)+
     &           (  fxnm1)*(1-fynm1)*(  fznm1)*
     &      gb_zz_nm1(inm1_invboost+1,jnm1_invboost,knm1_invboost+1)+
     &           (1-fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_zz_nm1(inm1_invboost,jnm1_invboost+1,knm1_invboost+1)+
     &           (  fxnm1)*(  fynm1)*(  fznm1)*
     &      gb_zz_nm1(inm1_invboost+1,jnm1_invboost+1,knm1_invboost+1)




               gb_tt_np1(i,j,k)=gb_tt_atinvboostptnp1
               gb_tx_np1(i,j,k)=gb_tx_atinvboostptnp1
               gb_ty_np1(i,j,k)=gb_ty_atinvboostptnp1
               gb_tz_np1(i,j,k)=gb_tz_atinvboostptnp1
               gb_xx_np1(i,j,k)=gb_xx_atinvboostptnp1
               gb_xy_np1(i,j,k)=gb_xy_atinvboostptnp1
               gb_xz_np1(i,j,k)=gb_xz_atinvboostptnp1
               gb_yy_np1(i,j,k)=gb_yy_atinvboostptnp1
               gb_yz_np1(i,j,k)=gb_yz_atinvboostptnp1
               gb_zz_np1(i,j,k)=gb_zz_atinvboostptnp1

               gb_tt_nm1(i,j,k)=gb_tt_atinvboostptnm1
               gb_tx_nm1(i,j,k)=gb_tx_atinvboostptnm1
               gb_ty_nm1(i,j,k)=gb_ty_atinvboostptnm1
               gb_tz_nm1(i,j,k)=gb_tz_atinvboostptnm1
               gb_xx_nm1(i,j,k)=gb_xx_atinvboostptnm1
               gb_xy_nm1(i,j,k)=gb_xy_atinvboostptnm1
               gb_xz_nm1(i,j,k)=gb_xz_atinvboostptnm1
               gb_yy_nm1(i,j,k)=gb_yy_atinvboostptnm1
               gb_yz_nm1(i,j,k)=gb_yz_atinvboostptnm1
               gb_zz_nm1(i,j,k)=gb_zz_atinvboostptnm1


             end if
            end do
           end do
        end do
       end if

        return
        end
