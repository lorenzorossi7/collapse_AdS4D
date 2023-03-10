c----------------------------------------------------------------------
c miscellaneous numerical routines
c----------------------------------------------------------------------

c-----------------------------------------------------------------------
c specific x/y first derivative routines, including support
c for excision. (separate AdS boundaries no longer supported)
c
c only "method 4" copied over here from gh3d ... if that's 
c good enough simplifies routines a bit. Also,
c can easily add back y derivatives.
c-----------------------------------------------------------------------
        subroutine df1_int_x(f,f_x,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz,i,j,k
        real*8 f(Nx,Ny,Nz),chr(Nx,Ny,Nz),ex,f_x,x(Nx),y(Ny),z(Nz)

        real*8 dx
        real*8 x0,y0,z0
        real*8 xi,yj,zk
        real*8 xip1,yjp1,zkp1
        real*8 xim1,yjm1,zkm1
        real*8 rhoijk,rhoip1jk,rhoijp1k,rhoip1jp1k,rhoijkp1
        real*8 rhoip1jkp1,rhoijp1kp1,rhoip1jp1kp1
        real*8 rhoim1jk,rhoijm1k,rhoim1jp1k,rhoijkm1
        real*8 rhoim1jkm1,rhoijm1km1,rhoim1jm1km1

        !--------------------------------------------------------------

        logical first
        save first
        data first/.true./

        logical extrap
        data extrap/.true./

        dx=x(2)-x(1)

        xi=x(i)
        yj=y(j)
        zk=z(k)
        rhoijk=sqrt(xi**2+yj**2+zk**2)
        xip1=x(i+1)
        rhoip1jk=sqrt(xip1**2+yj**2+zk**2)
        xim1=x(i-1)
        rhoim1jk=sqrt(xim1**2+yj**2+zk**2)


!!!!!!!!!!!!MYVERSION
        if (i.eq.1) then
               if ((.not.extrap)
     &            .and.(chr(i+1,j,k).ne.ex)
     &            .and.(chr(i+2,j,k).ne.ex)
     &            .and.(chr(i+3,j,k).ne.ex)) then
                   f_x=(-4*f(i,j,k)+7*f(i+1,j,k)
     &                  -4*f(i+2,j,k)+f(i+3,j,k))/2/dx
               else if ((chr(i+1,j,k).ne.ex
     &                 .and.chr(i+2,j,k).ne.ex)) then
                   f_x=(-3*f(i,j,k)+4*f(i+1,j,k)-f(i+2,j,k))/2/dx
               else if (chr(i+1,j,k).ne.ex) then
                   f_x=(-f(i,j,k)+f(i+1,j,k))/dx
!                  write(*,*) 'df1_int_x: warning ... i=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,Nx,Ny,Nz,dx
               else
                 if (first) then
                     first=.false.
                     write(*,*) 'df1_int_x: error in chr stencil (A)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dx,
     &                           xi,yj,zk,rhoijk,
     &                           xip1,rhoip1jk,xim1,rhoim1jk='
     &                           ,i,j,k,Nx,Ny,Nz,dx,
     &                           xi,yj,zk,rhoijk,
     &                           xip1,rhoip1jk,xim1,rhoim1jk
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dx,
     &                           chr(i,j,k),chr(i+1,j,k),chr(i-1,j,k)='
     &                          ,i,j,k,Nx,Ny,Nz,dx,
     &                           chr(i,j,k),chr(i+1,j,k),chr(i-1,j,k)
                     write(*,*) '    (first error only)'
                 end if
                   f_x=0
                 return
               end if

        else if (i.eq.2) then
         if ((chr(i-1,j,k).ne.ex).and.(chr(i+1,j,k).ne.ex)) then
                   f_x=(f(i+1,j,k)-f(i-1,j,k))/2/dx
         else if (chr(i-1,j,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i+1,j,k).ne.ex)
     &            .and.(chr(i+2,j,k).ne.ex)
     &            .and.(chr(i+3,j,k).ne.ex)) then
                   f_x=(-4*f(i,j,k)+7*f(i+1,j,k)
     &                  -4*f(i+2,j,k)+f(i+3,j,k))/2/dx
               else if ((chr(i+1,j,k).ne.ex)
     &                 .and.(chr(i+2,j,k).ne.ex)) then
                   f_x=(-3*f(i,j,k)+4*f(i+1,j,k)-f(i+2,j,k))/2/dx
               else if (chr(i+1,j,k).ne.ex) then
                   f_x=(-f(i,j,k)+f(i+1,j,k))/dx
!                  write(*,*) 'df1_int_x: warning ... i=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,Nx,Ny,Nz,dx
               else
                if (first) then
                     first=.false.
                     first=.false.
                     write(*,*) 'df1_int_x: error in chr stencil (B)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dx,
     &                           xi,yj,zk,rhoijk,
     &                           xip1,rhoip1jk,xim1,rhoim1jk='
     &                           ,i,j,k,Nx,Ny,Nz,dx,
     &                           xi,yj,zk,rhoijk,
     &                           xip1,rhoip1jk,xim1,rhoim1jk
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dx,
     &                           chr(i,j,k),chr(i+1,j,k),chr(i-1,j,k)='
     &                          ,i,j,k,Nx,Ny,Nz,dx,
     &                           chr(i,j,k),chr(i+1,j,k),chr(i-1,j,k)
                     write(*,*) '    (first error only)'
                end if
                   f_x=0
                   return
               end if
         else   !this is the case where (i-1,j,k) is not excised and (i+1,j,k) is excised 
                   f_x=(f(i,j,k)-f(i-1,j,k))/dx
!                  write(*,*) 'df1_int_x: warning ... i=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,Nx,Ny,Nz,dx
         end if

        else if (i.eq.3) then
         if ((chr(i-1,j,k).ne.ex).and.(chr(i+1,j,k).ne.ex)) then
                   f_x=(f(i+1,j,k)-f(i-1,j,k))/2/dx
         else if (chr(i-1,j,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i+1,j,k).ne.ex)
     &            .and.(chr(i+2,j,k).ne.ex)
     &            .and.(chr(i+3,j,k).ne.ex)) then
                   f_x=(-4*f(i,j,k)+7*f(i+1,j,k)
     &                  -4*f(i+2,j,k)+f(i+3,j,k))/2/dx
               else if ((chr(i+1,j,k).ne.ex)
     &                 .and.(chr(i+2,j,k).ne.ex)) then
                   f_x=(-3*f(i,j,k)+4*f(i+1,j,k)-f(i+2,j,k))/2/dx
               else if (chr(i+1,j,k).ne.ex) then
                   f_x=(-f(i,j,k)+f(i+1,j,k))/dx
!                  write(*,*) 'df1_int_x: warning ... i=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,Nx,Ny,Nz,dx
               else
                if (first) then
                     first=.false.
                     first=.false.
                     write(*,*) 'df1_int_x: error in chr stencil (C)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dx,
     &                           xi,yj,zk,rhoijk,
     &                           xip1,rhoip1jk,xim1,rhoim1jk='
     &                           ,i,j,k,Nx,Ny,Nz,dx,
     &                           xi,yj,zk,rhoijk,
     &                           xip1,rhoip1jk,xim1,rhoim1jk
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dx,
     &                           chr(i,j,k),chr(i+1,j,k),chr(i-1,j,k)='
     &                          ,i,j,k,Nx,Ny,Nz,dx,
     &                           chr(i,j,k),chr(i+1,j,k),chr(i-1,j,k)
                     write(*,*) '    (first error only)'
                end if
                   f_x=0
                   return
               end if
         else 
               if (chr(i-2,j,k).ne.ex) then
                   f_x=(3*f(i,j,k)-4*f(i-1,j,k)+f(i-2,j,k))/2/dx
               else
                   f_x=(f(i,j,k)-f(i-1,j,k))/dx
!                  write(*,*) 'df1_int_x: warning ... i=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,Nx,Ny,Nz,dx
               end if
         end if


        else if ((i.ge.4).and.(i.le.(Nx-3))) then
         if ((chr(i-1,j,k).ne.ex).and.(chr(i+1,j,k).ne.ex)) then
                   f_x=(f(i+1,j,k)-f(i-1,j,k))/2/dx
         else if (chr(i-1,j,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i+1,j,k).ne.ex)
     &            .and.(chr(i+2,j,k).ne.ex)
     &            .and.(chr(i+3,j,k).ne.ex)) then
                   f_x=(-4*f(i,j,k)+7*f(i+1,j,k)
     &                  -4*f(i+2,j,k)+f(i+3,j,k))/2/dx
               else if ((chr(i+1,j,k).ne.ex)
     &                 .and.(chr(i+2,j,k).ne.ex)) then
!                 if (z(k).eq.0.0d0) then
!                   write(*,*) "forward stencil at i,j,k,x(i),y(j),z(k)="
!     &                                            ,i,j,k,x(i),y(j),z(k)
!                   write(*,*) "chr(i,j,k),chr(i-1,j,k)="
!     &                        ,chr(i,j,k),chr(i-1,j,k)
!                   write(*,*) "f(i,j,k),f(i+1,j,k),f(i+2,j,k)="
!     &                        ,f(i,j,k),f(i+1,j,k),f(i+2,j,k)
!                 end if
                   f_x=(-3*f(i,j,k)+4*f(i+1,j,k)-f(i+2,j,k))/2/dx
               else if (chr(i+1,j,k).ne.ex) then
                   f_x=(-f(i,j,k)+f(i+1,j,k))/dx
!                  write(*,*) 'df1_int_x: warning ... i=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,Nx,Ny,Nz,dx
               else
                if (first) then
                     first=.false.
                     first=.false.
                     write(*,*) 'df1_int_x: error in chr stencil (D)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dx,
     &                           xi,yj,zk,rhoijk,
     &                           xip1,rhoip1jk,xim1,rhoim1jk='
     &                           ,i,j,k,Nx,Ny,Nz,dx,
     &                           xi,yj,zk,rhoijk,
     &                           xip1,rhoip1jk,xim1,rhoim1jk
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dx,
     &                           chr(i,j,k),chr(i+1,j,k),chr(i-1,j,k)='
     &                          ,i,j,k,Nx,Ny,Nz,dx,
     &                           chr(i,j,k),chr(i+1,j,k),chr(i-1,j,k)
                     write(*,*) '    (first error only)'
                end if
                   f_x=0
                   return
               end if
         else
               if ((.not.extrap)
     &            .and.(chr(i-3,j,k).ne.ex)
     &            .and.(chr(i-2,j,k).ne.ex)) then
                   f_x=(4*f(i,j,k)-7*f(i-1,j,k)
     &                  +4*f(i-2,j,k)-f(i-3,j,k))/2/dx
               else if (chr(i-2,j,k).ne.ex) then
                   f_x=(3*f(i,j,k)-4*f(i-1,j,k)+f(i-2,j,k))/2/dx
               else
                   f_x=(f(i,j,k)-f(i-1,j,k))/dx
!                  write(*,*) 'df1_int_x: warning ... i=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,Nx,Ny,Nz,dx
               end if
         end if


        else if (i.eq.(Nx-2)) then
         if ((chr(i+1,j,k).ne.ex).and.(chr(i-1,j,k).ne.ex)) then
                   f_x=(f(i+1,j,k)-f(i-1,j,k))/2/dx
         else if (chr(i+1,j,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i-1,j,k).ne.ex)
     &            .and.(chr(i-2,j,k).ne.ex)
     &            .and.(chr(i-3,j,k).ne.ex)) then
                   f_x=(4*f(i,j,k)-7*f(i-1,j,k)
     &                  +4*f(i-2,j,k)-f(i-3,j,k))/2/dx
               else if ((chr(i-1,j,k).ne.ex)
     &                 .and.(chr(i-2,j,k).ne.ex)) then
                   f_x=(3*f(i,j,k)-4*f(i-1,j,k)+f(i-2,j,k))/2/dx
               else if (chr(i-1,j,k).ne.ex) then
                   f_x=(f(i,j,k)-f(i-1,j,k))/dx
!                  write(*,*) 'df1_int_x: warning ... i=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,Nx,Ny,Nz,dx
               else
                if (first) then
                     first=.false.
                     first=.false.
                     write(*,*) 'df1_int_x: error in chr stencil (E)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dx,
     &                           xi,yj,zk,rhoijk,
     &                           xip1,rhoip1jk,xim1,rhoim1jk='
     &                           ,i,j,k,Nx,Ny,Nz,dx,
     &                           xi,yj,zk,rhoijk,
     &                           xip1,rhoip1jk,xim1,rhoim1jk
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dx,
     &                           chr(i,j,k),chr(i+1,j,k),chr(i-1,j,k)='
     &                          ,i,j,k,Nx,Ny,Nz,dx,
     &                           chr(i,j,k),chr(i+1,j,k),chr(i-1,j,k)
                     write(*,*) '    (first error only)'
                end if
                   f_x=0
                   return
               end if
         else 
               if (chr(i+2,j,k).ne.ex) then
                   f_x=(-3*f(i,j,k)+4*f(i+1,j,k)-f(i+2,j,k))/2/dx
               else
                   f_x=(-f(i,j,k)+f(i+1,j,k))/dx
!                  write(*,*) 'df1_int_x: warning ... i=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,Nx,Ny,Nz,dx
               end if
         end if

        else if (i.eq.(Nx-1)) then
         if ((chr(i+1,j,k).ne.ex).and.(chr(i-1,j,k).ne.ex)) then
                   f_x=(f(i+1,j,k)-f(i-1,j,k))/2/dx
         else if (chr(i+1,j,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i-1,j,k).ne.ex)
     &            .and.(chr(i-2,j,k).ne.ex)
     &            .and.(chr(i-3,j,k).ne.ex)) then
                   f_x=(4*f(i,j,k)-7*f(i-1,j,k)
     &                  +4*f(i-2,j,k)-f(i-3,j,k))/2/dx
               else if ((chr(i-1,j,k).ne.ex)
     &                 .and.(chr(i-2,j,k).ne.ex)) then
                   f_x=(3*f(i,j,k)-4*f(i-1,j,k)+f(i-2,j,k))/2/dx
               else if (chr(i-1,j,k).ne.ex) then
                   f_x=(f(i,j,k)-f(i-1,j,k))/dx
!                  write(*,*) 'df1_int_x: warning ... i=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,Nx,Ny,Nz,dx
               else
                if (first) then
                     first=.false.
                     first=.false.
                     write(*,*) 'df1_int_x: error in chr stencil (F)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dx,
     &                           xi,yj,zk,rhoijk,
     &                           xip1,rhoip1jk,xim1,rhoim1jk='
     &                           ,i,j,k,Nx,Ny,Nz,dx,
     &                           xi,yj,zk,rhoijk,
     &                           xip1,rhoip1jk,xim1,rhoim1jk
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dx,
     &                           chr(i,j,k),chr(i+1,j,k),chr(i-1,j,k)='
     &                          ,i,j,k,Nx,Ny,Nz,dx,
     &                           chr(i,j,k),chr(i+1,j,k),chr(i-1,j,k)
                     write(*,*) '    (first error only)'
                end if
                   f_x=0
                   return
               end if
         else
                   f_x=(-f(i,j,k)+f(i+1,j,k))/dx
!                  write(*,*) 'df1_int_x: warning ... i=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,Nx,Ny,Nz,dx
         end if


        else if (i.eq.Nx) then
               if ((.not.extrap)
     &            .and.(chr(i-1,j,k).ne.ex)
     &            .and.(chr(i-2,j,k).ne.ex)
     &            .and.(chr(i-3,j,k).ne.ex)) then
                   f_x=(4*f(i,j,k)-7*f(i-1,j,k)
     &                  +4*f(i-2,j,k)-f(i-3,j,k))/2/dx
               else if ((chr(i-1,j,k).ne.ex)
     &                 .and.(chr(i-2,j,k).ne.ex)) then
                   f_x=(3*f(i,j,k)-4*f(i-1,j,k)+f(i-2,j,k))/2/dx
               else if (chr(i-1,j,k).ne.ex) then
                   f_x=(f(i,j,k)-f(i-1,j,k))/dx
!                  write(*,*) 'df1_int_x: warning ... i=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,Nx,Ny,Nz,dx
               else
                 if (first) then
                     first=.false.
                     first=.false.
                     write(*,*) 'df1_int_x: error in chr stencil (G)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dx,
     &                           xi,yj,zk,rhoijk,
     &                           xip1,rhoip1jk,xim1,rhoim1jk='
     &                           ,i,j,k,Nx,Ny,Nz,dx,
     &                           xi,yj,zk,rhoijk,
     &                           xip1,rhoip1jk,xim1,rhoim1jk
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dx,
     &                           chr(i,j,k),chr(i+1,j,k),chr(i-1,j,k)='
     &                          ,i,j,k,Nx,Ny,Nz,dx,
     &                           chr(i,j,k),chr(i+1,j,k),chr(i-1,j,k)
                     write(*,*) '    (first error only)'
                 end if
                   f_x=0
                 return
               end if
        end if

!!!!!!!!!!!!!!!!!


!!!!!!!!OLD VERSION!!!!!!!!!!!!!!!
!        if (i.eq.1.or.(chr(i-1,j,k).eq.ex)) then
!           if (i.le.(Nx-3)
!     &         .and.((chr(i+1,j,k).ne.ex
!     &         .and.chr(i+2,j,k).ne.ex
!     &         .and.chr(i+3,j,k).ne.ex))) then
!             f_x=(-4*f(i,j,k)+7*f(i+1,j,k)-4*f(i+2,j,k)+f(i+3,j,k))/2/dx
!           else if (i.le.(Nx-2)
!     &              .and.((chr(i+1,j,k).ne.ex
!     &              .and.chr(i+2,j,k).ne.ex))) then
!              f_x=(-3*f(i,j,k)+4*f(i+1,j,k)-f(i+2,j,k))/2/dx
!           else if (i.le.(Nx-1).and.chr(i+1,j,k).ne.ex) then
!              f_x=(-f(i,j,k)+f(i+1,j,k))/dx
!!              write(*,*) 'df1_int_x: warning ... i=1 first order'
!!              write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
!           else
!              if (first) then
!                 first=.false.
!                 write(*,*) 'df1_int_x: error in chr stencil (A)'
!                 write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
!                 write(*,*) '    (first error only)'
!              end if
!              f_x=0
!              return
!           end if
!        else if (i.eq.Nx.or.(chr(i+1,j,k).eq.ex)) then
!           if (i.ge.4
!     &         .and.((chr(i-1,j,k).ne.ex
!     &         .and.chr(i-2,j,k).ne.ex
!     &         .and.chr(i-3,j,k).ne.ex))) then
!             f_x=(4*f(i,j,k)-7*f(i-1,j,k)+4*f(i-2,j,k)-f(i-3,j,k))/2/dx
!           else if (i.ge.3
!     &              .and.((chr(i-1,j,k).ne.ex
!     &              .and.chr(i-2,j,k).ne.ex))) then
!              f_x=(3*f(i,j,k)-4*f(i-1,j,k)+f(i-2,j,k))/2/dx
!           else if (i.ge.2.and.chr(i-1,j,k).ne.ex) then
!              f_x=(f(i,j,k)-f(i-1,j,k))/dx
!!              write(*,*) 'df1_int: warning ... i=Nx first order'
!!              write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
!           else
!              if (first) then
!                 first=.false.
!                 write(*,*) 'df1_int: error in chr stencil (B)'
!                 write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
!                 write(*,*) '    (first error only)'
!              end if
!              f_x=0
!              return
!           end if
!        else
!           if ((chr(i+1,j,k).ne.ex.and.chr(i-1,j,k).ne.ex)) then
!              f_x=(f(i+1,j,k)-f(i-1,j,k))/2/dx
!           else
!              if (first) then
!                 first=.false.
!                 write(*,*) 'df1_int: error in chr stencil (C)'
!                 write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
!                 write(*,*) '    (first error only)'
!              end if
!              f_x=0
!              return
!           end if
!        end if
!!!!!!!!!!!!!!!!!!!!!!

        return
        end
!----------------------------------------------------------------------
        subroutine df1_int_y(f,f_y,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz,i,j,k
        real*8 f(Nx,Ny,Nz),chr(Nx,Ny,Nz),ex,f_y,x(Nx),y(Ny),z(Nz)

        real*8 dy
        logical first
        save first
        data first/.true./

        logical extrap
        data extrap/.true./

        real*8 dx
        real*8 x0,y0,z0
        real*8 xi,yj,zk
        real*8 xip1,yjp1,zkp1
        real*8 xim1,yjm1,zkm1
        real*8 rhoijk,rhoip1jk,rhoijp1k,rhoip1jp1k,rhoijkp1
        real*8 rhoip1jkp1,rhoijp1kp1,rhoip1jp1kp1
        real*8 rhoim1jk,rhoijm1k,rhoim1jp1k,rhoijkm1
        real*8 rhoim1jkm1,rhoijm1km1,rhoim1jm1km1

        !--------------------------------------------------------------

        dy=y(2)-y(1)

        xi=x(i)
        yj=y(j)
        zk=z(k)
        rhoijk=sqrt(xi**2+yj**2+zk**2)
        yjp1=y(j+1)
        rhoijp1k=sqrt(xi**2+yjp1**2+zk**2)
        yjm1=y(j-1)
        rhoijm1k=sqrt(xi**2+yjm1**2+zk**2)

!!!!!!!MYVERSION!!!!!!!!!!!!

        if (j.eq.1) then
               if ((.not.extrap)
     &            .and.(chr(i,j+1,k).ne.ex)
     &            .and.(chr(i,j+2,k).ne.ex)
     &            .and.(chr(i,j+3,k).ne.ex)) then
                   f_y=(-4*f(i,j,k)+7*f(i,j+1,k)
     &                  -4*f(i,j+2,k)+f(i,j+3,k))/2/dy
               else if ((chr(i,j+1,k).ne.ex
     &                 .and.chr(i,j+2,k).ne.ex)) then
                   f_y=(-3*f(i,j,k)+4*f(i,j+1,k)-f(i,j+2,k))/2/dy
               else if (chr(i,j+1,k).ne.ex) then
                   f_y=(-f(i,j,k)+f(i,j+1,k))/dy
!                  write(*,*) 'df1_int_y: warning ... j=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
               else
                 if (first) then
                     first=.false.
                     write(*,*) 'df1_int_y: error in chr stencil (A)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dy,
     &                           xi,yj,zk,rhoijk,
     &                           yjp1,rhoijp1k,yjm1,rhoijm1k='
     &                          ,i,j,k,Nx,Ny,Nz,dy,
     &                           xi,yj,zk,rhoijk,
     &                           yjp1,rhoijp1k,yjm1,rhoijm1k
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dy,
     &                           chr(i,j,k),chr(i,j+1,k),chr(i,j-1,k)='
     &                          ,i,j,k,Nx,Ny,Nz,dy,
     &                           chr(i,j,k),chr(i,j+1,k),chr(i,j-1,k)
                     write(*,*) '    (first error only)'
                 end if
                   f_y=0
                 return
               end if

        else if (j.eq.2) then
         if ((chr(i,j-1,k).ne.ex).and.(chr(i,j+1,k).ne.ex)) then
                   f_y=(f(i,j+1,k)-f(i,j-1,k))/2/dy
         else if (chr(i,j-1,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j+1,k).ne.ex)
     &            .and.(chr(i,j+2,k).ne.ex)
     &            .and.(chr(i,j+3,k).ne.ex)) then
                   f_y=(-4*f(i,j,k)+7*f(i,j+1,k)
     &                  -4*f(i,j+2,k)+f(i,j+3,k))/2/dy
               else if ((chr(i,j+1,k).ne.ex)
     &                 .and.(chr(i,j+2,k).ne.ex)) then
                   f_y=(-3*f(i,j,k)+4*f(i,j+1,k)-f(i,j+2,k))/2/dy
               else if (chr(i,j+1,k).ne.ex) then
                   f_y=(-f(i,j,k)+f(i,j+1,k))/dy
!                  write(*,*) 'df1_int_y: warning ... j=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df1_int_y: error in chr stencil (B)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dy,
     &                           xi,yj,zk,rhoijk,
     &                           yjp1,rhoijp1k,yjm1,rhoijm1k='
     &                           ,i,j,k,Nx,Ny,Nz,dy,
     &                           xi,yj,zk,rhoijk,
     &                           yjp1,rhoijp1k,yjm1,rhoijm1k
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dy,
     &                           chr(i,j,k),chr(i,j+1,k),chr(i,j-1,k)='
     &                          ,i,j,k,Nx,Ny,Nz,dy,
     &                           chr(i,j,k),chr(i,j+1,k),chr(i,j-1,k)
                     write(*,*) '    (first error only)'
                end if
                   f_y=0
                   return
               end if
         else
                   f_y=(f(i,j,k)-f(i,j-1,k))/dy
!                  write(*,*) 'df1_int_y: warning ... j=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
         end if

        else if (j.eq.3) then
         if ((chr(i,j-1,k).ne.ex).and.(chr(i,j+1,k).ne.ex)) then
                   f_y=(f(i,j+1,k)-f(i,j-1,k))/2/dy
         else if (chr(i,j-1,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j+1,k).ne.ex)
     &            .and.(chr(i,j+2,k).ne.ex)
     &            .and.(chr(i,j+3,k).ne.ex)) then
                   f_y=(-4*f(i,j,k)+7*f(i,j+1,k)
     &                  -4*f(i,j+2,k)+f(i,j+3,k))/2/dy
               else if ((chr(i,j+1,k).ne.ex)
     &                 .and.(chr(i,j+2,k).ne.ex)) then
                   f_y=(-3*f(i,j,k)+4*f(i,j+1,k)-f(i,j+2,k))/2/dy
               else if (chr(i+1,j,k).ne.ex) then
                   f_y=(-f(i,j,k)+f(i,j+1,k))/dy
!                  write(*,*) 'df1_int_y: warning ... j=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df1_int_y: error in chr stencil (C)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dy,
     &                           xi,yj,zk,rhoijk,
     &                           yjp1,rhoijp1k,yjm1,rhoijm1k='
     &                           ,i,j,k,Nx,Ny,Nz,dy,
     &                           xi,yj,zk,rhoijk,
     &                           yjp1,rhoijp1k,yjm1,rhoijm1k
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dy,
     &                           chr(i,j,k),chr(i,j+1,k),chr(i,j-1,k)='
     &                          ,i,j,k,Nx,Ny,Nz,dy,
     &                           chr(i,j,k),chr(i,j+1,k),chr(i,j-1,k)
                     write(*,*) '    (first error only)'
                end if
                   f_y=0
                   return
               end if
         else
               if (chr(i,j-2,k).ne.ex) then
                   f_y=(3*f(i,j,k)-4*f(i,j-1,k)+f(i,j-2,k))/2/dy
               else
                   f_y=(f(i,j,k)-f(i,j-1,k))/dy
!                  write(*,*) 'df1_int_y: warning ... j=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
               end if
         end if

        else if ((j.ge.4).and.(j.le.(Ny-3))) then
         if ((chr(i,j-1,k).ne.ex).and.(chr(i,j+1,k).ne.ex)) then
                   f_y=(f(i,j+1,k)-f(i,j-1,k))/2/dy
         else if (chr(i,j-1,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j+1,k).ne.ex)
     &            .and.(chr(i,j+2,k).ne.ex)
     &            .and.(chr(i,j+3,k).ne.ex)) then
                   f_y=(-4*f(i,j,k)+7*f(i,j+1,k)
     &                  -4*f(i,j+2,k)+f(i,j+3,k))/2/dy
               else if ((chr(i,j+1,k).ne.ex)
     &                 .and.(chr(i,j+2,k).ne.ex)) then
                   f_y=(-3*f(i,j,k)+4*f(i,j+1,k)-f(i,j+2,k))/2/dy
               else if (chr(i,j+1,k).ne.ex) then
                   f_y=(-f(i,j,k)+f(i,j+1,k))/dy
!                  write(*,*) 'df1_int_y: warning ... j=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df1_int_y: error in chr stencil (D)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dy,
     &                           xi,yj,zk,rhoijk,
     &                           yjp1,rhoijp1k,yjm1,rhoijm1k='
     &                           ,i,j,k,Nx,Ny,Nz,dy,
     &                           xi,yj,zk,rhoijk,
     &                           yjp1,rhoijp1k,yjm1,rhoijm1k
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dy,
     &                           chr(i,j,k),chr(i,j+1,k),chr(i,j-1,k)='
     &                          ,i,j,k,Nx,Ny,Nz,dy,
     &                           chr(i,j,k),chr(i,j+1,k),chr(i,j-1,k)
                     write(*,*) '    (first error only)'
                end if
                   f_y=0
                   return
               end if
         else
               if ((.not.extrap)
     &            .and.(chr(i,j-3,k).ne.ex)
     &            .and.(chr(i,j-2,k).ne.ex)) then
                   f_y=(4*f(i,j,k)-7*f(i,j-1,k)
     &                  +4*f(i,j-2,k)-f(i,j-3,k))/2/dy
               else if (chr(i,j-2,k).ne.ex) then
                   f_y=(3*f(i,j,k)-4*f(i,j-1,k)+f(i,j-2,k))/2/dy
               else
                   f_y=(f(i,j,k)-f(i,j-1,k))/dy
!                  write(*,*) 'df1_int_y: warning ... j=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
               end if
         end if

        else if (j.eq.(Ny-2)) then
         if ((chr(i,j+1,k).ne.ex).and.(chr(i,j-1,k).ne.ex)) then
                   f_y=(f(i,j+1,k)-f(i,j-1,k))/2/dy
         else if (chr(i,j+1,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j-1,k).ne.ex)
     &            .and.(chr(i,j-2,k).ne.ex)
     &            .and.(chr(i,j-3,k).ne.ex)) then
                   f_y=(4*f(i,j,k)-7*f(i,j-1,k)
     &                  +4*f(i,j-2,k)-f(i,j-3,k))/2/dy
               else if ((chr(i,j-1,k).ne.ex)
     &                 .and.(chr(i,j-2,k).ne.ex)) then
                   f_y=(3*f(i,j,k)-4*f(i,j-1,k)+f(i,j-2,k))/2/dy
               else if (chr(i,j-1,k).ne.ex) then
                   f_y=(f(i,j,k)-f(i,j-1,k))/dy
!                  write(*,*) 'df1_int_y: warning ... j=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df1_int_y: error in chr stencil (E)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dy,
     &                           xi,yj,zk,rhoijk,
     &                           yjp1,rhoijp1k,yjm1,rhoijm1k='
     &                           ,i,j,k,Nx,Ny,Nz,dy,
     &                           xi,yj,zk,rhoijk,
     &                           yjp1,rhoijp1k,yjm1,rhoijm1k
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dy,
     &                           chr(i,j,k),chr(i,j+1,k),chr(i,j-1,k)='
     &                          ,i,j,k,Nx,Ny,Nz,dy,
     &                           chr(i,j,k),chr(i,j+1,k),chr(i,j-1,k)
                     write(*,*) '    (first error only)'
                end if
                   f_y=0
                   return
               end if
         else
               if (chr(i,j+2,k).ne.ex) then
                   f_y=(-3*f(i,j,k)+4*f(i,j+1,k)-f(i,j+2,k))/2/dy
               else
                   f_y=(-f(i,j,k)+f(i,j+1,k))/dy
!                  write(*,*) 'df1_int_y: warning ... j=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
               end if
         end if

        else if (j.eq.(Ny-1)) then
         if ((chr(i,j+1,k).ne.ex).and.(chr(i,j-1,k).ne.ex)) then
                   f_y=(f(i,j+1,k)-f(i,j-1,k))/2/dy
         else if (chr(i,j+1,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j-1,k).ne.ex)
     &            .and.(chr(i,j-2,k).ne.ex)
     &            .and.(chr(i,j-3,k).ne.ex)) then
                   f_y=(4*f(i,j,k)-7*f(i,j-1,k)
     &                  +4*f(i,j-2,k)-f(i,j-3,k))/2/dy
               else if ((chr(i,j-1,k).ne.ex)
     &                 .and.(chr(i,j-2,k).ne.ex)) then
                   f_y=(3*f(i,j,k)-4*f(i,j-1,k)+f(i,j-2,k))/2/dy
               else if (chr(i,j-1,k).ne.ex) then
                   f_y=(f(i,j,k)-f(i,j-1,k))/dy
!                  write(*,*) 'df1_int_y: warning ... j=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df1_int_y: error in chr stencil (F)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dy,
     &                           xi,yj,zk,rhoijk,
     &                           yjp1,rhoijp1k,yjm1,rhoijm1k='
     &                           ,i,j,k,Nx,Ny,Nz,dy,
     &                           xi,yj,zk,rhoijk,
     &                           yjp1,rhoijp1k,yjm1,rhoijm1k
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dy,
     &                           chr(i,j,k),chr(i,j+1,k),chr(i,j-1,k)='
     &                          ,i,j,k,Nx,Ny,Nz,dy,
     &                           chr(i,j,k),chr(i,j+1,k),chr(i,j-1,k)
                     write(*,*) '    (first error only)'
                end if
                   f_y=0
                   return
               end if
         else
                   f_y=(-f(i,j,k)+f(i,j+1,k))/dy
!                  write(*,*) 'df1_int_y: warning ... j=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
         end if

        else if (j.eq.Ny) then
               if ((.not.extrap)
     &            .and.(chr(i,j-1,k).ne.ex)
     &            .and.(chr(i,j-2,k).ne.ex)
     &            .and.(chr(i,j-3,k).ne.ex)) then
                   f_y=(4*f(i,j,k)-7*f(i,j-1,k)
     &                  +4*f(i,j-2,k)-f(i,j-3,k))/2/dy
               else if ((chr(i,j-1,k).ne.ex
     &                 .and.chr(i,j-2,k).ne.ex)) then
                   f_y=(3*f(i,j,k)-4*f(i,j-1,k)+f(i,j-2,k))/2/dy
               else if (chr(i,j-1,k).ne.ex) then
                   f_y=(f(i,j,k)-f(i,j-1,k))/dy
!                  write(*,*) 'df1_int_y: warning ... j=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
               else
                 if (first) then
                     first=.false.
                     write(*,*) 'df1_int_y: error in chr stencil (G)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dy,
     &                           xi,yj,zk,rhoijk,
     &                           yjp1,rhoijp1k,yjm1,rhoijm1k='
     &                           ,i,j,k,Nx,Ny,Nz,dy,
     &                           xi,yj,zk,rhoijk,
     &                           yjp1,rhoijp1k,yjm1,rhoijm1k
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dy,
     &                           chr(i,j,k),chr(i,j+1,k),chr(i,j-1,k)='
     &                          ,i,j,k,Nx,Ny,Nz,dy,
     &                           chr(i,j,k),chr(i,j+1,k),chr(i,j-1,k)
                     write(*,*) '    (first error only)'
                 end if
                   f_y=0
                 return
               end if
        end if
!!!!!!!!!!!!!!!!!

!!!!!!!!OLD VERSION!!!!!!!!!!!!!!!!!!
!        if ((j.eq.1).or.(chr(i,j-1,k).eq.ex)) then
!           if (j.le.(Ny-3)
!     &         .and.((chr(i,j+1,k).ne.ex
!     &         .and.chr(i,j+2,k).ne.ex
!     &         .and.chr(i,j+3,k).ne.ex))) then
!             f_y=(-4*f(i,j,k)+7*f(i,j+1,k)-4*f(i,j+2,k)+f(i,j+3,k))/2/dy
!           else if (j.le.(Ny-2).and.((chr(i,j+1,k).ne.ex
!     &              .and.chr(i,j+2,k).ne.ex))) then
!              f_y=(-3*f(i,j,k)+4*f(i,j+1,k)-f(i,j+2,k))/2/dy
!           else if (j.le.(Ny-1).and.chr(i,j+1,k).ne.ex) then
!              f_y=(-f(i,j,k)+f(i,j+1,k))/dy
!!              write(*,*) 'df1_int_y: warning ... j=1 first order'
!!              write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
!           else
!              if (first) then
!                 first=.false.
!                 write(*,*) 'df1_int_y: error in chr stencil (D)'
!                 write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
!                 write(*,*) '    (first error only)'
!              end if
!              f_y=0
!              return
!           end if
!        else if ((j.eq.Ny).or.(chr(i,j+1,k).eq.ex)) then
!           if (j.ge.4
!     &         .and.((chr(i,j-1,k).ne.ex
!     &         .and.chr(i,j-2,k).ne.ex
!     &         .and.chr(i,j-3,k).ne.ex))) then
!             f_y=(4*f(i,j,k)-7*f(i,j-1,k)+4*f(i,j-2,k)-f(i,j-3,k))/2/dy
!           else if (j.ge.3.and.((chr(i,j-1,k).ne.ex
!     &              .and.chr(i,j-2,k).ne.ex))) then
!              f_y=(3*f(i,j,k)-4*f(i,j-1,k)+f(i,j-2,k))/2/dy
!           else if (j.ge.2.and.chr(i,j-1,k).ne.ex) then
!              f_y=(f(i,j,k)-f(i,j-1,k))/dy
!!              write(*,*) 'df1_int_y: warning ... j=Ny first order'
!!              write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
!           else
!              if (first) then
!                 first=.false.
!                 write(*,*) 'df1_int_y: error in chr stencil (E)'
!                 write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
!                 write(*,*) '    (first error only)'
!              end if
!              f_y=0
!              return
!           end if
!        else
!           if ((chr(i,j+1,k).ne.ex.and.chr(i,j-1,k).ne.ex)) then
!              f_y=(f(i,j+1,k)-f(i,j-1,k))/2/dy
!           else
!              if (first) then
!                 first=.false.
!                 write(*,*) 'df1_int_y: error in chr stencil (F)'
!                 write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
!                 write(*,*) '    (first error only)'
!              end if
!              f_y=0
!              return
!           end if
!        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!

        return
        end

!----------------------------------------------------------------------
        subroutine df1_int_z(f,f_z,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz,i,j,k
        real*8 f(Nx,Ny,Nz),chr(Nx,Ny,Nz),ex,f_z,x(Nx),y(Ny),z(Nz)

        real*8 dz
        logical first
        save first
        data first/.true./

        logical extrap
        data extrap/.true./

        real*8 dx
        real*8 x0,y0,z0
        real*8 xi,yj,zk
        real*8 xip1,yjp1,zkp1
        real*8 xim1,yjm1,zkm1
        real*8 rhoijk,rhoip1jk,rhoijp1k,rhoip1jp1k,rhoijkp1
        real*8 rhoip1jkp1,rhoijp1kp1,rhoip1jp1kp1
        real*8 rhoim1jk,rhoijm1k,rhoim1jp1k,rhoijkm1
        real*8 rhoim1jkm1,rhoijm1km1,rhoim1jm1km1

        !--------------------------------------------------------------

        dz=z(2)-z(1)

        xi=x(i)
        yj=y(j)
        zk=z(k)
        rhoijk=sqrt(xi**2+yj**2+zk**2)
        zkp1=z(k+1)
        rhoijkp1=sqrt(xi**2+yj**2+zkp1**2)
        zkm1=z(k-1)
        rhoijkm1=sqrt(xi**2+yj**2+zkm1**2)

!        f_z=0


        if (k.eq.1) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k+1).ne.ex)
     &            .and.(chr(i,j,k+2).ne.ex)
     &            .and.(chr(i,j,k+3).ne.ex)) then
                   f_z=(-4*f(i,j,k)+7*f(i,j,k+1)
     &                  -4*f(i,j,k+2)+f(i,j,k+3))/2/dz
               else if ((chr(i,j,k+1).ne.ex
     &                 .and.chr(i,j,k+2).ne.ex)) then
                   f_z=(-3*f(i,j,k)+4*f(i,j,k+1)-f(i,j,k+2))/2/dz
               else if (chr(i,j,k+1).ne.ex) then
                   f_z=(-f(i,j,k)+f(i,j,k+1))/dz
!                  write(*,*) 'df1_int_z: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               else
                 if (first) then
                     first=.false.
                     write(*,*) 'df1_int_z: error in chr stencil (A)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dz,
     &                           xi,yj,zk,rhoijk,
     &                           zkp1,rhoijkp1,zkm1,rhoijkm1='
     &                           ,i,j,k,Nx,Ny,Nz,dz,
     &                           xi,yj,zk,rhoijk,
     &                           zkp1,rhoijkp1,zkm1,rhoijkm1
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dz,
     &                           chr(i,j,k),chr(i,j,k+1),chr(i,j,k-1)='
     &                          ,i,j,k,Nx,Ny,Nz,dz,
     &                           chr(i,j,k),chr(i,j,k+1),chr(i,j,k-1)
                     write(*,*) '    (first error only)'
                 end if
                   f_z=0
                 return
               end if

        else if (k.eq.2) then
         if ((chr(i,j,k-1).ne.ex).and.(chr(i,j,k+1).ne.ex)) then
                   f_z=(f(i,j,k+1)-f(i,j,k-1))/2/dz
         else if (chr(i,j,k-1).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k+1).ne.ex)
     &            .and.(chr(i,j,k+2).ne.ex)
     &            .and.(chr(i,j,k+3).ne.ex)) then
                   f_z=(-4*f(i,j,k)+7*f(i,j,k+1)
     &                  -4*f(i,j,k+2)+f(i,j,k+3))/2/dz
               else if ((chr(i,j,k+1).ne.ex)
     &                 .and.(chr(i,j,k+2).ne.ex)) then
                   f_z=(-3*f(i,j,k)+4*f(i,j,k+1)-f(i,j,k+2))/2/dz
               else if (chr(i,j,k+1).ne.ex) then
                   f_z=(-f(i,j,k)+f(i,j,k+1))/dz
!                  write(*,*) 'df1_int_z: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df1_int_z: error in chr stencil (B)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dz,
     &                           xi,yj,zk,rhoijk,
     &                           zkp1,rhoijkp1,zkm1,rhoijkm1='
     &                           ,i,j,k,Nx,Ny,Nz,dz,
     &                           xi,yj,zk,rhoijk,
     &                           zkp1,rhoijkp1,zkm1,rhoijkm1
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dz,
     &                           chr(i,j,k),chr(i,j,k+1),chr(i,j,k-1)='
     &                          ,i,j,k,Nx,Ny,Nz,dz,
     &                           chr(i,j,k),chr(i,j,k+1),chr(i,j,k-1)
                     write(*,*) '    (first error only)'
                end if
                   f_z=0
                   return
               end if
         else
                   f_z=(f(i,j,k)-f(i,j,k-1))/dz
!                  write(*,*) 'df1_int_z: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
         end if

        else if (k.eq.3) then
         if ((chr(i,j,k-1).ne.ex).and.(chr(i,j,k+1).ne.ex)) then
                   f_z=(f(i,j,k+1)-f(i,j,k-1))/2/dz
         else if (chr(i,j,k-1).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k+1).ne.ex)
     &            .and.(chr(i,j,k+2).ne.ex)
     &            .and.(chr(i,j,k+3).ne.ex)) then
                   f_z=(-4*f(i,j,k)+7*f(i,j,k+1)
     &                  -4*f(i,j,k+2)+f(i,j,k+3))/2/dz
               else if ((chr(i,j,k+1).ne.ex)
     &                 .and.(chr(i,j,k+2).ne.ex)) then
                   f_z=(-3*f(i,j,k)+4*f(i,j,k+1)-f(i,j,k+2))/2/dz
               else if (chr(i+1,j,k).ne.ex) then
                   f_z=(-f(i,j,k)+f(i,j,k+1))/dz
!                  write(*,*) 'df1_int_z: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df1_int_z: error in chr stencil (C)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dz,
     &                           xi,yj,zk,rhoijk,
     &                           zkp1,rhoijkp1,zkm1,rhoijkm1='
     &                           ,i,j,k,Nx,Ny,Nz,dz,
     &                           xi,yj,zk,rhoijk,
     &                           zkp1,rhoijkp1,zkm1,rhoijkm1
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dz,
     &                           chr(i,j,k),chr(i,j,k+1),chr(i,j,k-1)='
     &                          ,i,j,k,Nx,Ny,Nz,dz,
     &                           chr(i,j,k),chr(i,j,k+1),chr(i,j,k-1)
                     write(*,*) '    (first error only)'
                end if
                   f_z=0
                   return
               end if
         else
               if (chr(i,j,k-2).ne.ex) then
                   f_z=(3*f(i,j,k)-4*f(i,j,k-1)+f(i,j,k-2))/2/dz
               else
                   f_z=(f(i,j,k)-f(i,j,k-1))/dz
!                  write(*,*) 'df1_int_z: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               end if
         end if

        else if ((k.ge.4).and.(k.le.(Nz-3))) then
         if ((chr(i,j,k-1).ne.ex).and.(chr(i,j,k+1).ne.ex)) then
                   f_z=(f(i,j,k+1)-f(i,j,k-1))/2/dz
         else if (chr(i,j,k-1).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k+1).ne.ex)
     &            .and.(chr(i,j,k+2).ne.ex)
     &            .and.(chr(i,j,k+3).ne.ex)) then
                   f_z=(-4*f(i,j,k)+7*f(i,j,k+1)
     &                  -4*f(i,j,k+2)+f(i,j,k+3))/2/dz
               else if ((chr(i,j,k+1).ne.ex)
     &                 .and.(chr(i,j,k+2).ne.ex)) then
                   f_z=(-3*f(i,j,k)+4*f(i,j,k+1)-f(i,j,k+2))/2/dz
               else if (chr(i,j,k+1).ne.ex) then
                   f_z=(-f(i,j,k)+f(i,j,k+1))/dz
!                  write(*,*) 'df1_int_z: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df1_int_z: error in chr stencil (D)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dz,
     &                           xi,yj,zk,rhoijk,
     &                           zkp1,rhoijkp1,zkm1,rhoijkm1='
     &                           ,i,j,k,Nx,Ny,Nz,dz,
     &                           xi,yj,zk,rhoijk,
     &                           zkp1,rhoijkp1,zkm1,rhoijkm1
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dz,
     &                           chr(i,j,k),chr(i,j,k+1),chr(i,j,k-1)='
     &                          ,i,j,k,Nx,Ny,Nz,dz,
     &                           chr(i,j,k),chr(i,j,k+1),chr(i,j,k-1)
                     write(*,*) '    (first error only)'
                end if
                   f_z=0
                   return
               end if
         else
               if ((.not.extrap)
     &            .and.(chr(i,j,k-3).ne.ex)
     &            .and.(chr(i,j,k-2).ne.ex)) then
                   f_z=(4*f(i,j,k)-7*f(i,j,k-1)
     &                  +4*f(i,j,k-2)-f(i,j,k-3))/2/dz
               else if (chr(i,j,k-2).ne.ex) then
                   f_z=(3*f(i,j,k)-4*f(i,j,k-1)+f(i,j,k-2))/2/dz
               else
                   f_z=(f(i,j,k)-f(i,j,k-1))/dz
!                  write(*,*) 'df1_int_z: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               end if
         end if

        else if (k.eq.(Nz-2)) then
         if ((chr(i,j,k+1).ne.ex).and.(chr(i,j,k-1).ne.ex)) then
                   f_z=(f(i,j,k+1)-f(i,j,k-1))/2/dz
         else if (chr(i,j,k+1).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k-1).ne.ex)
     &            .and.(chr(i,j,k-2).ne.ex)
     &            .and.(chr(i,j,k-3).ne.ex)) then
                   f_z=(4*f(i,j,k)-7*f(i,j,k-1)
     &                  +4*f(i,j,k-2)-f(i,j,k-3))/2/dz
               else if ((chr(i,j,k-1).ne.ex)
     &                 .and.(chr(i,j,k-2).ne.ex)) then
                   f_z=(3*f(i,j,k)-4*f(i,j,k-1)+f(i,j,k-2))/2/dz
               else if (chr(i,j,k-1).ne.ex) then
                   f_z=(f(i,j,k)-f(i,j,k-1))/dz
!                  write(*,*) 'df1_int_z: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df1_int_z: error in chr stencil (E)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dz,
     &                           xi,yj,zk,rhoijk,
     &                           zkp1,rhoijkp1,zkm1,rhoijkm1='
     &                           ,i,j,k,Nx,Ny,Nz,dz,
     &                           xi,yj,zk,rhoijk,
     &                           zkp1,rhoijkp1,zkm1,rhoijkm1
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dz,
     &                           chr(i,j,k),chr(i,j,k+1),chr(i,j,k-1)='
     &                          ,i,j,k,Nx,Ny,Nz,dz,
     &                           chr(i,j,k),chr(i,j,k+1),chr(i,j,k-1)
                     write(*,*) '    (first error only)'
                end if
                   f_z=0
                   return
               end if
         else
               if (chr(i,j,k+2).ne.ex) then
                   f_z=(-3*f(i,j,k)+4*f(i,j,k+1)-f(i,j,k+2))/2/dz
               else
                   f_z=(-f(i,j,k)+f(i,j,k+1))/dz
!                  write(*,*) 'df1_int_z: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               end if
         end if

        else if (k.eq.(Nz-1)) then
         if ((chr(i,j,k+1).ne.ex).and.(chr(i,j,k-1).ne.ex)) then
                   f_z=(f(i,j,k+1)-f(i,j,k-1))/2/dz
         else if (chr(i,j,k+1).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k-1).ne.ex)
     &            .and.(chr(i,j,k-2).ne.ex)
     &            .and.(chr(i,j,k-3).ne.ex)) then
                   f_z=(4*f(i,j,k)-7*f(i,j,k-1)
     &                  +4*f(i,j,k-2)-f(i,j,k-3))/2/dz
               else if ((chr(i,j,k-1).ne.ex)
     &                 .and.(chr(i,j,k-2).ne.ex)) then
                   f_z=(3*f(i,j,k)-4*f(i,j,k-1)+f(i,j,k-2))/2/dz
               else if (chr(i,j,k-1).ne.ex) then
                   f_z=(f(i,j,k)-f(i,j,k-1))/dz
!                  write(*,*) 'df1_int_z: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df1_int_z: error in chr stencil (F)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dz,
     &                           xi,yj,zk,rhoijk,
     &                           zkp1,rhoijkp1,zkm1,rhoijkm1='
     &                           ,i,j,k,Nx,Ny,Nz,dz,
     &                           xi,yj,zk,rhoijk,
     &                           zkp1,rhoijkp1,zkm1,rhoijkm1
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dz,
     &                           chr(i,j,k),chr(i,j,k+1),chr(i,j,k-1)='
     &                          ,i,j,k,Nx,Ny,Nz,dz,
     &                           chr(i,j,k),chr(i,j,k+1),chr(i,j,k-1)
                     write(*,*) '    (first error only)'
                end if
                   f_z=0
                   return
               end if
         else
                   f_z=(-f(i,j,k)+f(i,j,k+1))/dz
!                  write(*,*) 'df1_int_z: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
         end if

        else if (k.eq.Nz) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k-1).ne.ex)
     &            .and.(chr(i,j,k-2).ne.ex)
     &            .and.(chr(i,j,k-3).ne.ex)) then
                   f_z=(4*f(i,j,k)-7*f(i,j,k-1)
     &                  +4*f(i,j,k-2)-f(i,j,k-3))/2/dz
               else if ((chr(i,j,k-1).ne.ex
     &                 .and.chr(i,j,k-2).ne.ex)) then
                   f_z=(3*f(i,j,k)-4*f(i,j,k-1)+f(i,j,k-2))/2/dz
               else if (chr(i,j,k-1).ne.ex) then
                   f_z=(f(i,j,k)-f(i,j,k-1))/dz
!                  write(*,*) 'df1_int_z: warning ... k=1 first order'
!                  write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
               else
                 if (first) then
                     first=.false.
                     write(*,*) 'df1_int_z: error in chr stencil (G)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dz,
     &                           xi,yj,zk,rhoijk,
     &                           zkp1,rhoijkp1,zkm1,rhoijkm1='
     &                           ,i,j,k,Nx,Ny,Nz,dz,
     &                           xi,yj,zk,rhoijk,
     &                           zkp1,rhoijkp1,zkm1,rhoijkm1
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dz,
     &                           chr(i,j,k),chr(i,j,k+1),chr(i,j,k-1)='
     &                          ,i,j,k,Nx,Ny,Nz,dz,
     &                           chr(i,j,k),chr(i,j,k+1),chr(i,j,k-1)
                     write(*,*) '    (first error only)'
                 end if
                   f_z=0
                 return
               end if
        end if

        return
        end

c----------------------------------------------------------------------
c the following computes all first derivatives of f,
c at a point i,j, at time level n.
c----------------------------------------------------------------------
        subroutine df1_int(f_np1,f_n,f_nm1,f_t,f_x,f_y,f_z,
     &                     x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,name)
        implicit none
        integer Nx,Ny,Nz,i,j,k
        real*8 f_np1(Nx,Ny,Nz),f_n(Nx,Ny,Nz),f_nm1(Nx,Ny,Nz)
        real*8 f_t,f_x,f_y,f_z,x(Nx),y(Ny),z(Nz),dt,ex,chr(Nx,Ny,Nz)
        character*(*) name

        logical first
        save first
        data first/.true./

        logical ltrace
        parameter (ltrace=.false.)

        !--------------------------------------------------------------

        if (chr(i,j,k).eq.ex) then
          write(*,*) 'df1_int: error ... point excised'
          stop
        end if

        f_t=(f_np1(i,j,k)-f_nm1(i,j,k))/2/dt
  
        call df1_int_x(f_n,f_x,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
        call df1_int_y(f_n,f_y,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
        call df1_int_z(f_n,f_z,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)

        if (ltrace) then
           write(*,*) 'df1_int for ',name
           write(*,*) ' f_t=',f_t
           write(*,*) ' f_x=',f_x
           write(*,*) ' f_y=',f_y
           write(*,*) ' f_z=',f_z
        end if

        return
        end

c----------------------------------------------------------------------
c same as df1_int above, but computes all second derivatives as well 
c
c if (extrap), then we use TE matched first derivative
c operators, plus 2nd order accurate 2nd derivatives ... this
c (with *NO* boundary dissipation) is equivalent to using 
c 4th order extrapolation (along the given direction where the 
c derivative is evaluated) 2-points into the excision surface, and then
c using the standard update scheme (now with interior dissipation)
c applied to the 'old' boundaries.
c
c For AdS: forget now whether extrap was essential for gh3d or not.
c          leaving on for now
c
c CALLING grid function f_n f here, to avoid conflict
c with f_n in the include stuff
c----------------------------------------------------------------------
        subroutine df2_int(f_np1,f,f_nm1,f_t,f_x,f_y,
     &                     f_z,
     &                     f_tt,f_tx,f_ty,
     &                     f_tz,
     &                     f_xx,f_xy,
     &                     f_xz,
     &                     f_yy,
     &                     f_yz,f_zz,
     &                     x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,name)
        implicit none
        integer Nx,Ny,Nz,i,j,k
        real*8 f_np1(Nx,Ny,Nz),f(Nx,Ny,Nz),f_nm1(Nx,Ny,Nz)
        real*8 f_t,f_x,f_y,f_z
        real*8 f_tt,f_tx,f_ty
        real*8 f_tz
        real*8 f_xx,f_xy
        real*8 f_xz
        real*8 f_yy
        real*8 f_yz,f_zz
        real*8 x(Nx),y(Ny),z(Nz),dt,ex,chr(Nx,Ny,Nz)
        character*(*) name

        real*8 dx,dy,dz

        logical first,extrap
        save first
        data first/.true./
        parameter (extrap=.true.)

        logical ltrace
        parameter (ltrace=.false.)

        real*8 f_x_np1,f_x_nm1,f_y_np1,f_y_nm1
        real*8 f_z_np1,f_z_nm1

        real*8 f_y_ip1,f_y_ip2
        real*8 f_z_ip1,f_z_ip2

        real*8 f_y_im1,f_y_im2
        real*8 f_z_im1,f_z_im2

        real*8 f_z_jp1,f_z_jp2
        real*8 f_z_jm1,f_z_jm2




        ! initialize fixed-size variables
        data f_x_np1,f_x_nm1,f_y_np1,f_y_nm1/0.0,0.0,0.0,0.0/
        data f_z_np1,f_z_nm1/0.0,0.0/
        data f_y_ip1,f_y_ip2/0.0,0.0/
        data f_z_ip1,f_z_ip2/0.0,0.0/
        data f_y_im1,f_y_im2/0.0,0.0/
        data f_z_im1,f_z_im2/0.0,0.0/
        data f_z_jp1,f_z_jp2/0.0,0.0/
        data f_z_jm1,f_z_jm2/0.0,0.0/

        !--------------------------------------------------------------
       
        call df1_int(f_np1,f,f_nm1,f_t,f_x,f_y,f_z,
     &               x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,name)

        f_tt=(f_np1(i,j,k)-2*f(i,j,k)+f_nm1(i,j,k))/dt/dt 

        f_xx=0
        f_xy=0
        f_xz=0
        f_yy=0
        f_yz=0
        f_zz=0

        if (chr(i,j,k).eq.ex) then
         write(*,*) 'df2_int: error ... point excised'
         stop
        end if

        call df1_int_x(f_np1,f_x_np1,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
        call df1_int_x(f_nm1,f_x_nm1,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)

        f_tx=(f_x_np1-f_x_nm1)/2/dt

        call df1_int_y(f_np1,f_y_np1,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
        call df1_int_y(f_nm1,f_y_nm1,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)

        f_ty=(f_y_np1-f_y_nm1)/2/dt

        call df1_int_z(f_np1,f_z_np1,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
        call df1_int_z(f_nm1,f_z_nm1,x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)

        f_tz=(f_z_np1-f_z_nm1)/2/dt

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

!!!!!!MYVERSION!!!!!!!!!

        !i

        if (i.eq.1) then
               if ((.not.extrap)
     &            .and.chr(i+1,j,k).ne.ex
     &            .and.chr(i+2,j,k).ne.ex
     &            .and.chr(i+3,j,k).ne.ex
     &            .and.chr(i+4,j,k).ne.ex) then
                   f_xx=(3*f(i,j,k)-9*f(i+1,j,k)+
     &                  10*f(i+2,j,k)-5*f(i+3,j,k)+
     &                  f(i+4,j,k))/dx/dx
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(-3*f_y+4*f_y_ip1-f_y_ip2)/2/dx
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(-3*f_z+4*f_z_ip1-f_z_ip2)/2/dx
               else if (chr(i+1,j,k).ne.ex
     &                 .and.chr(i+2,j,k).ne.ex
     &                 .and.chr(i+3,j,k).ne.ex) then
                   f_xx=(2*f(i,j,k)-5*f(i+1,j,k)+
     &                   4*f(i+2,j,k)-f(i+3,j,k))/dx/dx
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(-3*f_y+4*f_y_ip1-f_y_ip2)/2/dx
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(-3*f_z+4*f_z_ip1-f_z_ip2)/2/dx
               else if (chr(i+1,j,k).ne.ex
     &                  .and.chr(i+2,j,k).ne.ex) then
              !    write(*,*) 'df2_int: warning ... first order i=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
                   f_xx=(f(i+2,j,k)-2*f(i+1,j,k)+f(i,j,k))/dx/dx
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(-3*f_y+4*f_y_ip1-f_y_ip2)/2/dx
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(-3*f_z+4*f_z_ip1-f_z_ip2)/2/dx
               else
                 if (first) then
                   first=.false.
                   write(*,*) 'df2_int: error in chr stencil (A)'
                   write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,Nx,Ny,Nz,dx
                   write(*,*) '    (first error only)'
                 end if
                   return
               end if

        else if (i.eq.2) then
         if ((chr(i-1,j,k).ne.ex).and.(chr(i+1,j,k).ne.ex)) then
                   f_xx=(f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(f_y_ip1-f_y_im1)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(f_z_ip1-f_z_im1)/2/dx
         else if (chr(i-1,j,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i+1,j,k).ne.ex)
     &            .and.(chr(i+2,j,k).ne.ex)
     &            .and.(chr(i+3,j,k).ne.ex)
     &            .and.(chr(i+4,j,k).ne.ex)) then
                   f_xx=(3*f(i,j,k)-9*f(i+1,j,k)+
     &                  10*f(i+2,j,k)-5*f(i+3,j,k)+
     &                  f(i+4,j,k))/dx/dx
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(-3*f_y+4*f_y_ip1-f_y_ip2)/2/dx
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(-3*f_z+4*f_z_ip1-f_z_ip2)/2/dx
               else if (chr(i+1,j,k).ne.ex
     &                 .and.chr(i+2,j,k).ne.ex
     &                 .and.chr(i+3,j,k).ne.ex) then
                   f_xx=(2*f(i,j,k)-5*f(i+1,j,k)+
     &                   4*f(i+2,j,k)-f(i+3,j,k))/dx/dx
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(-3*f_y+4*f_y_ip1-f_y_ip2)/2/dx
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(-3*f_z+4*f_z_ip1-f_z_ip2)/2/dx
               else if (chr(i+1,j,k).ne.ex
     &                 .and.chr(i+2,j,k).ne.ex) then
              !    write(*,*) 'df2_int: warning ... first order i=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
                   f_xx=(f(i+2,j,k)-2*f(i+1,j,k)+f(i,j,k))/dx/dx
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(-3*f_y+4*f_y_ip1-f_y_ip2)/2/dx
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(-3*f_z+4*f_z_ip1-f_z_ip2)/2/dx
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (B)'
                     write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,
     &                                                   Nx,Ny,Nz,dx
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         end if

        else if (i.eq.3) then
         if ((chr(i-1,j,k).ne.ex).and.(chr(i+1,j,k).ne.ex)) then
                   f_xx=(f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(f_y_ip1-f_y_im1)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(f_z_ip1-f_z_im1)/2/dx
         else if (chr(i-1,j,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i+1,j,k).ne.ex)
     &            .and.(chr(i+2,j,k).ne.ex)
     &            .and.(chr(i+3,j,k).ne.ex)
     &            .and.(chr(i+4,j,k).ne.ex)) then
                   f_xx=(3*f(i,j,k)-9*f(i+1,j,k)+
     &                  10*f(i+2,j,k)-5*f(i+3,j,k)+
     &                  f(i+4,j,k))/dx/dx
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(-3*f_y+4*f_y_ip1-f_y_ip2)/2/dx
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(-3*f_z+4*f_z_ip1-f_z_ip2)/2/dx
               else if ((chr(i+1,j,k).ne.ex)
     &                 .and.(chr(i+2,j,k).ne.ex)
     &                 .and.(chr(i+3,j,k).ne.ex)) then
                   f_xx=(2*f(i,j,k)-5*f(i+1,j,k)+
     &                   4*f(i+2,j,k)-f(i+3,j,k))/dx/dx
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(-3*f_y+4*f_y_ip1-f_y_ip2)/2/dx
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(-3*f_z+4*f_z_ip1-f_z_ip2)/2/dx
               else if ((chr(i+1,j,k).ne.ex)
     &                 .and.(chr(i+2,j,k).ne.ex)) then
              !    write(*,*) 'df2_int: warning ... first order i=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
                   f_xx=(f(i+2,j,k)-2*f(i+1,j,k)+f(i,j,k))/dx/dx
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(-3*f_y+4*f_y_ip1-f_y_ip2)/2/dx
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(-3*f_z+4*f_z_ip1-f_z_ip2)/2/dx
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (C)'
                     write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,
     &                                                   Nx,Ny,Nz,dx
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         else  !this is the case where (i-1,j,k) is not excised and (i+1,j,k) is excised
               if (chr(i-2,j,k).ne.ex) then
                   f_xx=(f(i,j,k)-2*f(i-1,j,k)+f(i-2,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(3*f_y-4*f_y_im1+f_y_im2)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(3*f_z-4*f_z_im1+f_z_im2)/2/dx
               end if
         end if

        else if ((i.ge.4).and.(i.le.(Nx-4))) then  !we need to impose i.le.Nx-4 so that chr(i+4,j,k) exists in the condition with .not.extrap. Therefore we will also need to add the case i.eq.Nx-3 below
         if ((chr(i-1,j,k).ne.ex).and.(chr(i+1,j,k).ne.ex)) then
                   f_xx=(f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(f_y_ip1-f_y_im1)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(f_z_ip1-f_z_im1)/2/dx
         else if (chr(i-1,j,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i+1,j,k).ne.ex)
     &            .and.(chr(i+2,j,k).ne.ex)
     &            .and.(chr(i+3,j,k).ne.ex)
     &            .and.(chr(i+4,j,k).ne.ex)) then
                   f_xx=(3*f(i,j,k)-9*f(i+1,j,k)+
     &                  10*f(i+2,j,k)-5*f(i+3,j,k)+
     &                  f(i+4,j,k))/dx/dx
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(-3*f_y+4*f_y_ip1-f_y_ip2)/2/dx
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(-3*f_z+4*f_z_ip1-f_z_ip2)/2/dx
               else if ((chr(i+1,j,k).ne.ex)
     &                 .and.(chr(i+2,j,k).ne.ex)
     &                 .and.(chr(i+3,j,k).ne.ex)) then
                   f_xx=(2*f(i,j,k)-5*f(i+1,j,k)+
     &                   4*f(i+2,j,k)-f(i+3,j,k))/dx/dx
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(-3*f_y+4*f_y_ip1-f_y_ip2)/2/dx
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(-3*f_z+4*f_z_ip1-f_z_ip2)/2/dx
               else if ((chr(i+1,j,k).ne.ex)
     &                 .and.(chr(i+2,j,k).ne.ex)) then
              !    write(*,*) 'df2_int: warning ... first order i=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
                   f_xx=(f(i+2,j,k)-2*f(i+1,j,k)+f(i,j,k))/dx/dx
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(-3*f_y+4*f_y_ip1-f_y_ip2)/2/dx
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(-3*f_z+4*f_z_ip1-f_z_ip2)/2/dx
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (D)'
                     write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,
     &                                                   Nx,Ny,Nz,dx
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         else
               if ((chr(i-3,j,k).ne.ex)
     &            .and.(chr(i-2,j,k).ne.ex)) then
                   f_xx=(2*f(i,j,k)-5*f(i-1,j,k)
     &                   +4*f(i-2,j,k)-f(i-3,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(3*f_y-4*f_y_im1+f_y_im2)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(3*f_z-4*f_z_im1+f_z_im2)/2/dx
               else if (chr(i-2,j,k).ne.ex) then
                   f_xx=(f(i,j,k)-2*f(i-1,j,k)+f(i-2,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(3*f_y-4*f_y_im1+f_y_im2)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(3*f_z-4*f_z_im1+f_z_im2)/2/dx
               end if
         end if



!NEW BIT: added to avoid going out of array boundaries when i=Nx-3 AND we extrap is false
        else if (i.eq.(Nx-3)) then
         if ((chr(i+1,j,k).ne.ex).and.(chr(i-1,j,k).ne.ex)) then
                   f_xx=(f(i-1,j,k)-2*f(i,j,k)+f(i+1,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(f_y_ip1-f_y_im1)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(f_z_ip1-f_z_im1)/2/dx
         else if (chr(i+1,j,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i-1,j,k).ne.ex)
     &            .and.(chr(i-2,j,k).ne.ex)
     &            .and.(chr(i-3,j,k).ne.ex)
     &            .and.(chr(i-4,j,k).ne.ex)) then
                   f_xx=(3*f(i,j,k)-9*f(i-1,j,k)+
     &                  10*f(i-2,j,k)-5*f(i-3,j,k)+
     &                  f(i-4,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(3*f_y-4*f_y_im1+f_y_im2)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(3*f_z-4*f_z_im1+f_z_im2)/2/dx
               else if ((chr(i-1,j,k).ne.ex)
     &                 .and.(chr(i-2,j,k).ne.ex)
     &                 .and.(chr(i-3,j,k).ne.ex)) then
                   f_xx=(2*f(i,j,k)-5*f(i-1,j,k)+
     &                   4*f(i-2,j,k)-f(i-3,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(3*f_y-4*f_y_im1+f_y_im2)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(3*f_z-4*f_z_im1+f_z_im2)/2/dx
               else if ((chr(i-1,j,k).ne.ex)
     &                 .and.(chr(i-2,j,k).ne.ex)) then
              !    write(*,*) 'df2_int: warning ... first order i=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
                   f_xx=(f(i-2,j,k)-2*f(i-1,j,k)+f(i,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(3*f_y-4*f_y_im1+f_y_im2)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(3*f_z-4*f_z_im1+f_z_im2)/2/dx
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (E)'
                     write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,
     &                                                   Nx,Ny,Nz,dx
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         else  !this is the case where (i+1,j,k) is not excised and (i-1,j,k) is excised
               if (chr(i+2,j,k).ne.ex) then
                   f_xx=(f(i,j,k)-2*f(i+1,j,k)+f(i+2,j,k))/dx/dx
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(-3*f_y+4*f_y_ip1-f_y_ip2)/2/dx
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(-3*f_z+4*f_z_ip1-f_z_ip2)/2/dx
               end if
         end if

!!!!!!!!!!!!!!!!!!!!!!!!!






        else if (i.eq.(Nx-2)) then
         if ((chr(i+1,j,k).ne.ex).and.(chr(i-1,j,k).ne.ex)) then
                   f_xx=(f(i-1,j,k)-2*f(i,j,k)+f(i+1,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(f_y_ip1-f_y_im1)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(f_z_ip1-f_z_im1)/2/dx
         else if (chr(i+1,j,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i-1,j,k).ne.ex)
     &            .and.(chr(i-2,j,k).ne.ex)
     &            .and.(chr(i-3,j,k).ne.ex)
     &            .and.(chr(i-4,j,k).ne.ex)) then
                   f_xx=(3*f(i,j,k)-9*f(i-1,j,k)+
     &                  10*f(i-2,j,k)-5*f(i-3,j,k)+
     &                  f(i-4,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(3*f_y-4*f_y_im1+f_y_im2)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(3*f_z-4*f_z_im1+f_z_im2)/2/dx
               else if ((chr(i-1,j,k).ne.ex)
     &                 .and.(chr(i-2,j,k).ne.ex)
     &                 .and.(chr(i-3,j,k).ne.ex)) then
                   f_xx=(2*f(i,j,k)-5*f(i-1,j,k)+
     &                   4*f(i-2,j,k)-f(i-3,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(3*f_y-4*f_y_im1+f_y_im2)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(3*f_z-4*f_z_im1+f_z_im2)/2/dx
               else if ((chr(i-1,j,k).ne.ex)
     &                 .and.(chr(i-2,j,k).ne.ex)) then
              !    write(*,*) 'df2_int: warning ... first order i=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
                   f_xx=(f(i-2,j,k)-2*f(i-1,j,k)+f(i,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(3*f_y-4*f_y_im1+f_y_im2)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(3*f_z-4*f_z_im1+f_z_im2)/2/dx
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (E)'
                     write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,
     &                                                   Nx,Ny,Nz,dx
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         else  !this is the case where (i+1,j,k) is not excised and (i-1,j,k) is excised
               if (chr(i+2,j,k).ne.ex) then
                   f_xx=(f(i,j,k)-2*f(i+1,j,k)+f(i+2,j,k))/dx/dx
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(-3*f_y+4*f_y_ip1-f_y_ip2)/2/dx
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip2,x,y,z,
     &                            i+2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(-3*f_z+4*f_z_ip1-f_z_ip2)/2/dx
               end if
         end if

        else if (i.eq.(Nx-1)) then
         if ((chr(i+1,j,k).ne.ex).and.(chr(i-1,j,k).ne.ex)) then
                   f_xx=(f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(f_y_ip1-f_y_im1)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_ip1,x,y,z,
     &                            i+1,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(f_z_ip1-f_z_im1)/2/dx
         else if (chr(i+1,j,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i-1,j,k).ne.ex)
     &            .and.(chr(i-2,j,k).ne.ex)
     &            .and.(chr(i-3,j,k).ne.ex)
     &            .and.(chr(i-4,j,k).ne.ex)) then
                   f_xx=(3*f(i,j,k)-9*f(i-1,j,k)+
     &                  10*f(i-2,j,k)-5*f(i-3,j,k)+
     &                  f(i-4,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(3*f_y-4*f_y_im1+f_y_im2)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(3*f_z-4*f_z_im1+f_z_im2)/2/dx
               else if (chr(i-1,j,k).ne.ex
     &                 .and.chr(i-2,j,k).ne.ex
     &                 .and.chr(i-3,j,k).ne.ex) then
                   f_xx=(2*f(i,j,k)-5*f(i-1,j,k)+
     &                   4*f(i-2,j,k)-f(i-3,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(3*f_y-4*f_y_im1+f_y_im2)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(3*f_z-4*f_z_im1+f_z_im2)/2/dx
               else if (chr(i-1,j,k).ne.ex
     &                 .and.chr(i-2,j,k).ne.ex) then
              !    write(*,*) 'df2_int: warning ... first order i=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
                   f_xx=(f(i-2,j,k)-2*f(i-1,j,k)+f(i,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(3*f_y-4*f_y_im1+f_y_im2)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(3*f_z-4*f_z_im1+f_z_im2)/2/dx
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (F)'
                     write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,
     &                                                   Nx,Ny,Nz,dx
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         end if

        else if (i.eq.Nx) then
               if ((.not.extrap)
     &            .and.(chr(i-1,j,k).ne.ex)
     &            .and.(chr(i-2,j,k).ne.ex)
     &            .and.(chr(i-3,j,k).ne.ex)
     &            .and.(chr(i-4,j,k).ne.ex)) then
                   f_xx=(3*f(i,j,k)-9*f(i-1,j,k)+
     &                  10*f(i-2,j,k)-5*f(i-3,j,k)+
     &                  f(i-4,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(3*f_y-4*f_y_im1+f_y_im2)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(3*f_z-4*f_z_im1+f_z_im2)/2/dx
               else if (chr(i-1,j,k).ne.ex
     &                 .and.chr(i-2,j,k).ne.ex
     &                 .and.chr(i-3,j,k).ne.ex) then
                   f_xx=(2*f(i,j,k)-5*f(i-1,j,k)+
     &                   4*f(i-2,j,k)-f(i-3,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(3*f_y-4*f_y_im1+f_y_im2)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(3*f_z-4*f_z_im1+f_z_im2)/2/dx
               else if (chr(i-1,j,k).ne.ex
     &                  .and.chr(i-2,j,k).ne.ex) then
              !    write(*,*) 'df2_int: warning ... first order i=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
                   f_xx=(f(i-2,j,k)-2*f(i-1,j,k)+f(i,j,k))/dx/dx
                   call df1_int_y(f,f_y_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_y(f,f_y_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xy=(3*f_y-4*f_y_im1+f_y_im2)/2/dx
                   call df1_int_z(f,f_z_im1,x,y,z,
     &                            i-1,j,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_im2,x,y,z,
     &                            i-2,j,k,chr,ex,Nx,Ny,Nz)
                   f_xz=(3*f_z-4*f_z_im1+f_z_im2)/2/dx
               else
                 if (first) then
                   first=.false.
                   write(*,*) 'df2_int: error in chr stencil (G)'
                   write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,Nx,Ny,Nz,dx
                   write(*,*) '    (first error only)'
                 end if
                   return
               end if
        end if


        !j

        if (j.eq.1) then
               if ((.not.extrap)
     &            .and.(chr(i,j+1,k).ne.ex)
     &            .and.(chr(i,j+2,k).ne.ex)
     &            .and.(chr(i,j+3,k).ne.ex)
     &            .and.(chr(i,j+4,k).ne.ex)) then
                   f_yy=(3*f(i,j,k)-9*f(i,j+1,k)+
     &                  10*f(i,j+2,k)-5*f(i,j+3,k)+
     &                  f(i,j+4,k))/dy/dy
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp2,x,y,z,
     &                            i,j+2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(-3*f_z+4*f_z_jp1-f_z_jp2)/2/dy
               else if (chr(i,j+1,k).ne.ex
     &                 .and.chr(i,j+2,k).ne.ex
     &                 .and.chr(i,j+3,k).ne.ex) then
                   f_yy=(2*f(i,j,k)-5*f(i,j+1,k)+
     &                   4*f(i,j+2,k)-f(i,j+3,k))/dy/dy
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp2,x,y,z,
     &                            i,j+2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(-3*f_z+4*f_z_jp1-f_z_jp2)/2/dy
               else if (chr(i,j+1,k).ne.ex
     &                  .and.chr(i,j+2,k).ne.ex) then
              !    write(*,*) 'df2_int: warning ... first order j=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
                   f_yy=(f(i,j+2,k)-2*f(i,j+1,k)+f(i,j,k))/dy/dy
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp2,x,y,z,
     &                            i,j+2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(-3*f_z+4*f_z_jp1-f_z_jp2)/2/dy
               else
                 if (first) then
                   first=.false.
                   write(*,*) 'df2_int: error in chr stencil (H)'
                   write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,Nx,Ny,Nz,dx
                   write(*,*) '    (first error only)'
                 end if
                   return
               end if

        else if (j.eq.2) then
         if ((chr(i,j-1,k).ne.ex).and.(chr(i,j+1,k).ne.ex)) then
                   f_yy=(f(i,j+1,k)-2*f(i,j,k)+f(i,j-1,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(f_z_jp1-f_z_jm1)/2/dy
         else if (chr(i,j-1,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j+1,k).ne.ex)
     &            .and.(chr(i,j+2,k).ne.ex)
     &            .and.(chr(i,j+3,k).ne.ex)
     &            .and.(chr(i,j+4,k).ne.ex)) then
                   f_yy=(3*f(i,j,k)-9*f(i,j+1,k)+
     &                  10*f(i,j+2,k)-5*f(i,j+3,k)+
     &                  f(i,j+4,k))/dy/dy
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp2,x,y,z,
     &                            i,j+2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(-3*f_z+4*f_z_jp1-f_z_jp2)/2/dy
               else if (chr(i,j+1,k).ne.ex
     &                 .and.chr(i,j+2,k).ne.ex
     &                 .and.chr(i,j+3,k).ne.ex) then
                   f_yy=(2*f(i,j,k)-5*f(i,j+1,k)+
     &                   4*f(i,j+2,k)-f(i,j+3,k))/dy/dy
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp2,x,y,z,
     &                            i,j+2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(-3*f_z+4*f_z_jp1-f_z_jp2)/2/dy
               else if (chr(i,j+1,k).ne.ex
     &                 .and.chr(i,j+2,k).ne.ex) then
              !    write(*,*) 'df2_int: warning ... first order j=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
                   f_yy=(f(i,j+2,k)-2*f(i,j+1,k)+f(i,j,k))/dy/dy
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp2,x,y,z,
     &                            i,j+2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(-3*f_z+4*f_z_jp1-f_z_jp2)/2/dy
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (I)'
                     write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,
     &                                                   Nx,Ny,Nz,dx
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         end if

        else if (j.eq.3) then
         if ((chr(i,j-1,k).ne.ex).and.(chr(i,j+1,k).ne.ex)) then
                   f_yy=(f(i,j+1,k)-2*f(i,j,k)+f(i,j-1,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(f_z_jp1-f_z_jm1)/2/dy
         else if (chr(i,j-1,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j+1,k).ne.ex)
     &            .and.(chr(i,j+2,k).ne.ex)
     &            .and.(chr(i,j+3,k).ne.ex)
     &            .and.(chr(i,j+4,k).ne.ex)) then
                   f_yy=(3*f(i,j,k)-9*f(i,j+1,k)+
     &                  10*f(i,j+2,k)-5*f(i,j+3,k)+
     &                  f(i,j+4,k))/dy/dy
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp2,x,y,z,
     &                            i,j+2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(-3*f_z+4*f_z_jp1-f_z_jp2)/2/dy
               else if ((chr(i,j+1,k).ne.ex)
     &                 .and.(chr(i,j+2,k).ne.ex)
     &                 .and.(chr(i,j+3,k).ne.ex)) then
                   f_yy=(2*f(i,j,k)-5*f(i,j+1,k)+
     &                   4*f(i,j+2,k)-f(i,j+3,k))/dy/dy
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp2,x,y,z,
     &                            i,j+2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(-3*f_z+4*f_z_jp1-f_z_jp2)/2/dy
               else if ((chr(i,j+1,k).ne.ex)
     &                 .and.(chr(i,j+2,k).ne.ex)) then
              !    write(*,*) 'df2_int: warning ... first order j=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
                   f_yy=(f(i,j+2,k)-2*f(i,j+1,k)+f(i,j,k))/dy/dy
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp2,x,y,z,
     &                            i,j+2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(-3*f_z+4*f_z_jp1-f_z_jp2)/2/dy
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (J)'
                     write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,
     &                                                   Nx,Ny,Nz,dx
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         else  !this is the case where (i,j-1,k) is not excised and (i,j+1,k) is excised
               if (chr(i,j-2,k).ne.ex) then
                   f_yy=(f(i,j,k)-2*f(i,j-1,k)+f(i,j-2,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jm2,x,y,z,
     &                            i,j-2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(3*f_z-4*f_z_jm1+f_z_jm2)/2/dy
               end if
         end if

        else if ((j.ge.4).and.(j.le.(Ny-4))) then !we need to impose j.le.Ny-4 so that chr(i,j+4,k) exists in the condition with .not.extrap. Therefore we will also need to add the case j.eq.Ny-3 below
         if ((chr(i,j-1,k).ne.ex).and.(chr(i,j+1,k).ne.ex)) then
                   f_yy=(f(i,j+1,k)-2*f(i,j,k)+f(i,j-1,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(f_z_jp1-f_z_jm1)/2/dy
         else if (chr(i,j-1,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j+1,k).ne.ex)
     &            .and.(chr(i,j+2,k).ne.ex)
     &            .and.(chr(i,j+3,k).ne.ex)
     &            .and.(chr(i,j+4,k).ne.ex)) then
                   f_yy=(3*f(i,j,k)-9*f(i,j+1,k)+
     &                  10*f(i,j+2,k)-5*f(i,j+3,k)+
     &                  f(i,j+4,k))/dy/dy
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp2,x,y,z,
     &                            i,j+2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(-3*f_z+4*f_z_jp1-f_z_jp2)/2/dy
               else if ((chr(i,j+1,k).ne.ex)
     &                 .and.(chr(i,j+2,k).ne.ex)
     &                 .and.(chr(i,j+3,k).ne.ex)) then
                   f_yy=(2*f(i,j,k)-5*f(i,j+1,k)+
     &                   4*f(i,j+2,k)-f(i,j+3,k))/dy/dy
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp2,x,y,z,
     &                            i,j+2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(-3*f_z+4*f_z_jp1-f_z_jp2)/2/dy
               else if ((chr(i,j+1,k).ne.ex)
     &                 .and.(chr(i,j+2,k).ne.ex)) then
              !    write(*,*) 'df2_int: warning ... first order j=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
                   f_yy=(f(i,j+2,k)-2*f(i,j+1,k)+f(i,j,k))/dy/dy
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp2,x,y,z,
     &                            i,j+2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(-3*f_z+4*f_z_jp1-f_z_jp2)/2/dy
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (K)'
                     write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,
     &                                                   Nx,Ny,Nz,dx
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         else
               if ((chr(i,j-3,k).ne.ex)
     &            .and.(chr(i,j-2,k).ne.ex)) then
                   f_yy=(2*f(i,j,k)-5*f(i,j-1,k)
     &                   +4*f(i,j-2,k)-f(i,j-3,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jm2,x,y,z,
     &                            i,j-2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(3*f_z-4*f_z_jm1+f_z_jm2)/2/dy
               else if (chr(i,j-2,k).ne.ex) then
                   f_yy=(f(i,j,k)-2*f(i,j-1,k)+f(i,j-2,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jm2,x,y,z,
     &                            i,j-2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(3*f_z-4*f_z_jm1+f_z_jm2)/2/dy
               end if
         end if

!NEW BIT: added to avoid going out of array boundaries when j=Ny-3 AND we extrap is false

        else if (j.eq.(Ny-3)) then
         if ((chr(i,j+1,k).ne.ex).and.(chr(i,j-1,k).ne.ex)) then
                   f_yy=(f(i,j-1,k)-2*f(i,j,k)+f(i,j+1,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(f_z_jp1-f_z_jm1)/2/dy
         else if (chr(i,j+1,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j-1,k).ne.ex)
     &            .and.(chr(i,j-2,k).ne.ex)
     &            .and.(chr(i,j-3,k).ne.ex)
     &            .and.(chr(i,j-4,k).ne.ex)) then
                   f_yy=(3*f(i,j,k)-9*f(i,j-1,k)+
     &                  10*f(i,j-2,k)-5*f(i,j-3,k)+
     &                  f(i,j-4,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jm2,x,y,z,
     &                            i,j-2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(3*f_z-4*f_z_jm1+f_z_jm2)/2/dy
               else if ((chr(i,j-1,k).ne.ex)
     &                 .and.(chr(i,j-2,k).ne.ex)
     &                 .and.(chr(i,j-3,k).ne.ex)) then
                   f_yy=(2*f(i,j,k)-5*f(i,j-1,k)+
     &                   4*f(i,j-2,k)-f(i,j-3,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jm2,x,y,z,
     &                            i,j-2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(3*f_z-4*f_z_jm1+f_z_jm2)/2/dy
               else if ((chr(i,j-1,k).ne.ex)
     &                 .and.(chr(i,j-2,k).ne.ex)) then
              !    write(*,*) 'df2_int: warning ... first order j=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
                   f_yy=(f(i,j-2,k)-2*f(i,j-1,k)+f(i,j,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jm2,x,y,z,
     &                            i,j-2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(3*f_z-4*f_z_jm1+f_z_jm2)/2/dy
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (L)'
                     write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,
     &                                                   Nx,Ny,Nz,dx
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         else  !this is the case where (i,j+1,k) is not excised and (i,j-1,k) is excised
               if (chr(i,j+2,k).ne.ex) then
                   f_yy=(f(i,j,k)-2*f(i,j+1,k)+f(i,j+2,k))/dy/dy
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp2,x,y,z,
     &                            i,j+2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(-3*f_z+4*f_z_jp1-f_z_jp2)/2/dy
               end if
         end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        else if (j.eq.(Ny-2)) then
         if ((chr(i,j+1,k).ne.ex).and.(chr(i,j-1,k).ne.ex)) then
                   f_yy=(f(i,j-1,k)-2*f(i,j,k)+f(i,j+1,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(f_z_jp1-f_z_jm1)/2/dy
         else if (chr(i,j+1,k).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j-1,k).ne.ex)
     &            .and.(chr(i,j-2,k).ne.ex)
     &            .and.(chr(i,j-3,k).ne.ex)
     &            .and.(chr(i,j-4,k).ne.ex)) then
                   f_yy=(3*f(i,j,k)-9*f(i,j-1,k)+
     &                  10*f(i,j-2,k)-5*f(i,j-3,k)+
     &                  f(i,j-4,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jm2,x,y,z,
     &                            i,j-2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(3*f_z-4*f_z_jm1+f_z_jm2)/2/dy
               else if ((chr(i,j-1,k).ne.ex)
     &                 .and.(chr(i,j-2,k).ne.ex)
     &                 .and.(chr(i,j-3,k).ne.ex)) then
                   f_yy=(2*f(i,j,k)-5*f(i,j-1,k)+
     &                   4*f(i,j-2,k)-f(i,j-3,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jm2,x,y,z,
     &                            i,j-2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(3*f_z-4*f_z_jm1+f_z_jm2)/2/dy
               else if ((chr(i,j-1,k).ne.ex)
     &                 .and.(chr(i,j-2,k).ne.ex)) then
              !    write(*,*) 'df2_int: warning ... first order j=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
                   f_yy=(f(i,j-2,k)-2*f(i,j-1,k)+f(i,j,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jm2,x,y,z,
     &                            i,j-2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(3*f_z-4*f_z_jm1+f_z_jm2)/2/dy
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (L)'
                     write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,
     &                                                   Nx,Ny,Nz,dx
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         else  !this is the case where (i,j+1,k) is not excised and (i,j-1,k) is excised
               if (chr(i,j+2,k).ne.ex) then
                   f_yy=(f(i,j,k)-2*f(i,j+1,k)+f(i,j+2,k))/dy/dy
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp2,x,y,z,
     &                            i,j+2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(-3*f_z+4*f_z_jp1-f_z_jp2)/2/dy
               end if
         end if

        else if (j.eq.(Ny-1)) then
         if ((chr(i,j+1,k).ne.ex).and.(chr(i,j-1,k).ne.ex)) then
                   f_yy=(f(i,j+1,k)-2*f(i,j,k)+f(i,j-1,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jp1,x,y,z,
     &                            i,j+1,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(f_z_jp1-f_z_jm1)/2/dy
         else if (chr(i,j+1,k).eq.ex) then
               if (.not.extrap
     &            .and.(chr(i,j-1,k).ne.ex)
     &            .and.(chr(i,j-2,k).ne.ex)
     &            .and.(chr(i,j-3,k).ne.ex)
     &            .and.(chr(i,j-4,k).ne.ex)) then
                   f_yy=(3*f(i,j,k)-9*f(i,j-1,k)+
     &                  10*f(i,j-2,k)-5*f(i,j-3,k)+
     &                  f(i,j-4,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jm2,x,y,z,
     &                            i,j-2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(3*f_z-4*f_z_jm1+f_z_jm2)/2/dy
               else if (chr(i,j-1,k).ne.ex
     &                 .and.chr(i,j-2,k).ne.ex
     &                 .and.chr(i,j-3,k).ne.ex) then
                   f_yy=(2*f(i,j,k)-5*f(i,j-1,k)+
     &                   4*f(i,j-2,k)-f(i,j-3,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jm2,x,y,z,
     &                            i,j-2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(3*f_z-4*f_z_jm1+f_z_jm2)/2/dy
               else if (chr(i,j-1,k).ne.ex
     &                 .and.chr(i,j-2,k).ne.ex) then
              !    write(*,*) 'df2_int: warning ... first order j=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
                   f_yy=(f(i,j-2,k)-2*f(i,j-1,k)+f(i,j,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jm2,x,y,z,
     &                            i,j-2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(3*f_z-4*f_z_jm1+f_z_jm2)/2/dy
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (M)'
                     write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,
     &                                                   Nx,Ny,Nz,dx
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         end if

        else if (j.eq.Ny) then
               if ((.not.extrap)
     &            .and.(chr(i,j-1,k).ne.ex)
     &            .and.(chr(i,j-2,k).ne.ex)
     &            .and.(chr(i,j-3,k).ne.ex)
     &            .and.(chr(i,j-4,k).ne.ex)) then
                   f_yy=(3*f(i,j,k)-9*f(i,j-1,k)+
     &                  10*f(i,j-2,k)-5*f(i,j-3,k)+
     &                  f(i,j-4,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jm2,x,y,z,
     &                            i,j-2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(3*f_z-4*f_z_jm1+f_z_jm2)/2/dy
               else if (chr(i,j-1,k).ne.ex
     &                 .and.chr(i,j-2,k).ne.ex
     &                 .and.chr(i,j-3,k).ne.ex) then
                   f_yy=(2*f(i,j,k)-5*f(i,j-1,k)+
     &                   4*f(i,j-2,k)-f(i,j-3,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jm2,x,y,z,
     &                            i,j-2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(3*f_z-4*f_z_jm1+f_z_jm2)/2/dy
               else if (chr(i,j-1,k).ne.ex
     &                  .and.chr(i,j-2,k).ne.ex) then
              !    write(*,*) 'df2_int: warning ... first order j=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dy=',i,j,k,Nx,Ny,Nz,dy
                   f_yy=(f(i,j-2,k)-2*f(i,j-1,k)+f(i,j,k))/dy/dy
                   call df1_int_z(f,f_z_jm1,x,y,z,
     &                            i,j-1,k,chr,ex,Nx,Ny,Nz)
                   call df1_int_z(f,f_z_jm2,x,y,z,
     &                            i,j-2,k,chr,ex,Nx,Ny,Nz)
                   f_yz=(3*f_z-4*f_z_jm1+f_z_jm2)/2/dy
               else
                 if (first) then
                   first=.false.
                   write(*,*) 'df2_int: error in chr stencil (N)'
                   write(*,*) '    i,j,k,Nx,Ny,Nz,dx=',i,j,k,Nx,Ny,Nz,dx
                   write(*,*) '    (first error only)'
                 end if
                   return
               end if
        end if


        !k

!        f_zz=0

        if (k.eq.1) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k+1).ne.ex)
     &            .and.(chr(i,j,k+2).ne.ex)
     &            .and.(chr(i,j,k+3).ne.ex)
     &            .and.(chr(i,j,k+4).ne.ex)) then
                   f_zz=(3*f(i,j,k)-9*f(i,j,k+1)+
     &                  10*f(i,j,k+2)-5*f(i,j,k+3)+
     &                  f(i,j,k+4))/dz/dz
               else if (chr(i,j,k+1).ne.ex
     &                 .and.chr(i,j,k+2).ne.ex
     &                 .and.chr(i,j,k+3).ne.ex) then
                   f_zz=(2*f(i,j,k)-5*f(i,j,k+1)+
     &                   4*f(i,j,k+2)-f(i,j,k+3))/dz/dz
               else if (chr(i,j,k+1).ne.ex
     &                  .and.chr(i,j,k+2).ne.ex) then
              !    write(*,*) 'df2_int: warning ... first order k=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
                   f_zz=(f(i,j,k+2)-2*f(i,j,k+1)+f(i,j,k))/dz/dz
               else
                 if (first) then
                   first=.false.
                   write(*,*) 'df2_int: error in chr stencil (O)'
                   write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
                   write(*,*) '    (first error only)'
                 end if
                   return
               end if

        else if (k.eq.2) then
         if ((chr(i,j,k-1).ne.ex).and.(chr(i,j,k+1).ne.ex)) then
                   f_zz=(f(i,j,k+1)-2*f(i,j,k)+f(i,j,k-1))/dz/dz
         else if (chr(i,j,k-1).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k+1).ne.ex)
     &            .and.(chr(i,j,k+2).ne.ex)
     &            .and.(chr(i,j,k+3).ne.ex)
     &            .and.(chr(i,j,k+4).ne.ex)) then
                   f_zz=(3*f(i,j,k)-9*f(i,j,k+1)+
     &                  10*f(i,j,k+2)-5*f(i,j,k+3)+
     &                  f(i,j,k+4))/dz/dz
               else if (chr(i,j,k+1).ne.ex
     &                 .and.chr(i,j,k+2).ne.ex
     &                 .and.chr(i,j,k+3).ne.ex) then
                   f_zz=(2*f(i,j,k)-5*f(i,j,k+1)+
     &                   4*f(i,j,k+2)-f(i,j,k+3))/dz/dz
               else if (chr(i,j,k+1).ne.ex
     &                 .and.chr(i,j,k+2).ne.ex) then
              !    write(*,*) 'df2_int: warning ... first order k=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
                   f_zz=(f(i,j,k+2)-2*f(i,j,k+1)+f(i,j,k))/dz/dz
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (P)'
                     write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,
     &                                                   Nx,Ny,Nz,dz
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         end if

        else if (k.eq.3) then
         if ((chr(i,j,k-1).ne.ex).and.(chr(i,j,k+1).ne.ex)) then
                   f_zz=(f(i,j,k+1)-2*f(i,j,k)+f(i,j,k-1))/dz/dz
         else if (chr(i,j,k-1).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k+1).ne.ex)
     &            .and.(chr(i,j,k+2).ne.ex)
     &            .and.(chr(i,j,k+3).ne.ex)
     &            .and.(chr(i,j,k+4).ne.ex)) then
                   f_zz=(3*f(i,j,k)-9*f(i,j,k+1)+
     &                  10*f(i,j,k+2)-5*f(i,j,k+3)+
     &                  f(i,j,k+4))/dz/dz
               else if ((chr(i,j,k+1).ne.ex)
     &                 .and.(chr(i,j,k+2).ne.ex)
     &                 .and.(chr(i,j,k+3).ne.ex)) then
                   f_zz=(2*f(i,j,k)-5*f(i,j,k+1)+
     &                   4*f(i,j,k+2)-f(i,j,k+3))/dz/dz
               else if ((chr(i,j,k+1).ne.ex)
     &                 .and.(chr(i,j,k+2).ne.ex)) then
              !    write(*,*) 'df2_int: warning ... first order k=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
                   f_zz=(f(i,j,k+2)-2*f(i,j,k+1)+f(i,j,k))/dz/dz
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (Q)'
                     write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,
     &                                                   Nx,Ny,Nz,dz
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         else  !this is the case where (i,j,k-1) is not excised and (i,j,k+1) is excised
               if (chr(i,j,k-2).ne.ex) then
                   f_zz=(f(i,j,k)-2*f(i,j,k-1)+f(i,j,k-2))/dz/dz
               end if
         end if

        else if ((k.ge.4).and.(k.le.(Nz-4))) then !we need to impose k.le.Nz-4 so that chr(i,j,k+4) exists in the condition with .not.extrap. Therefore we will also need to add the case k.eq.Nz-3 below
         if ((chr(i,j,k-1).ne.ex).and.(chr(i,j,k+1).ne.ex)) then
                   f_zz=(f(i,j,k+1)-2*f(i,j,k)+f(i,j,k-1))/dz/dz
         else if (chr(i,j,k-1).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k+1).ne.ex)
     &            .and.(chr(i,j,k+2).ne.ex)
     &            .and.(chr(i,j,k+3).ne.ex)
     &            .and.(chr(i,j,k+4).ne.ex)) then
                   f_zz=(3*f(i,j,k)-9*f(i,j,k+1)+
     &                  10*f(i,j,k+2)-5*f(i,j,k+3)+
     &                  f(i,j,k+4))/dz/dz
               else if ((chr(i,j,k+1).ne.ex)
     &                 .and.(chr(i,j,k+2).ne.ex)
     &                 .and.(chr(i,j,k+3).ne.ex)) then
                   f_zz=(2*f(i,j,k)-5*f(i,j,k+1)+
     &                   4*f(i,j,k+2)-f(i,j,k+3))/dz/dz
               else if ((chr(i,j,k+1).ne.ex)
     &                 .and.(chr(i,j,k+2).ne.ex)) then
              !    write(*,*) 'df2_int: warning ... first order k=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
                   f_zz=(f(i,j,k+2)-2*f(i,j,k+1)+f(i,j,k))/dz/dz
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (R)'
                     write(*,*) 'i,j,k,Nx,Ny,Nz,dx,
     &                           xi,yj,zk,rhoijk,
     &                           zkp1,zkm1,zkp2,zkm2,
     &                           chr(i,j,k),chr(i,j,k+1),chr(i,j,k-1),
     &                           chr(i,j,k+2),chr(i,j,k-2),ex='
     &                           ,i,j,k,Nx,Ny,Nz,dx,
     &                           x(i),y(j),z(k),sqrt(x(i)*x(i)
     &                            +y(j)*y(j)+z(k)*z(k)),
     &                           z(k+1),z(k-1),z(k+2),z(k-2),
     &                           chr(i,j,k),chr(i,j,k+1),chr(i,j,k-1),
     &                           chr(i,j,k+2),chr(i,j,k-2),ex
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         else
               if ((chr(i,j,k-3).ne.ex)
     &            .and.(chr(i,j,k-2).ne.ex)) then
                   f_zz=(2*f(i,j,k)-5*f(i,j,k-1)
     &                   +4*f(i,j,k-2)-f(i,j,k-3))/dz/dz
               else if (chr(i,j,k-2).ne.ex) then
                   f_zz=(f(i,j,k)-2*f(i,j,k-1)+f(i,j,k-2))/dz/dz
               end if
         end if


!NEW BIT: added to avoid going out of array boundaries when k=Nz-3 AND we extrap is false

        else if (k.eq.(Nz-3)) then
         if ((chr(i,j,k+1).ne.ex).and.(chr(i,j,k-1).ne.ex)) then
                   f_zz=(f(i,j,k-1)-2*f(i,j,k)+f(i,j,k+1))/dz/dz
         else if (chr(i,j,k+1).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k-1).ne.ex)
     &            .and.(chr(i,j,k-2).ne.ex)
     &            .and.(chr(i,j,k-3).ne.ex)
     &            .and.(chr(i,j,k-4).ne.ex)) then
                   f_zz=(3*f(i,j,k)-9*f(i,j,k-1)+
     &                  10*f(i,j,k-2)-5*f(i,j,k-3)+
     &                  f(i,j,k-4))/dz/dz
               else if ((chr(i,j,k-1).ne.ex)
     &                 .and.(chr(i,j,k-2).ne.ex)
     &                 .and.(chr(i,j,k-3).ne.ex)) then
                   f_zz=(2*f(i,j,k)-5*f(i,j,k-1)+
     &                   4*f(i,j,k-2)-f(i,j,k-3))/dz/dz
               else if ((chr(i,j,k-1).ne.ex)
     &                 .and.(chr(i,j,k-2).ne.ex)) then
              !    write(*,*) 'df2_int: warning ... first order k=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
                   f_zz=(f(i,j,k-2)-2*f(i,j,k-1)+f(i,j,k))/dz/dz
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (S)'
                     write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,
     &                                                   Nx,Ny,Nz,dz
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         else  !this is the case where (i,j,k+1) is not excised and (i,j,k-1) is excised
               if (chr(i,j,k+2).ne.ex) then
                   f_zz=(f(i,j,k)-2*f(i,j,k+1)+f(i,j,k+2))/dz/dz
               end if
         end if

        else if (k.eq.(Nz-2)) then
         if ((chr(i,j,k+1).ne.ex).and.(chr(i,j,k-1).ne.ex)) then
                   f_zz=(f(i,j,k-1)-2*f(i,j,k)+f(i,j,k+1))/dz/dz
         else if (chr(i,j,k+1).eq.ex) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k-1).ne.ex)
     &            .and.(chr(i,j,k-2).ne.ex)
     &            .and.(chr(i,j,k-3).ne.ex)
     &            .and.(chr(i,j,k-4).ne.ex)) then
                   f_zz=(3*f(i,j,k)-9*f(i,j,k-1)+
     &                  10*f(i,j,k-2)-5*f(i,j,k-3)+
     &                  f(i,j,k-4))/dz/dz
               else if ((chr(i,j,k-1).ne.ex)
     &                 .and.(chr(i,j,k-2).ne.ex)
     &                 .and.(chr(i,j,k-3).ne.ex)) then
                   f_zz=(2*f(i,j,k)-5*f(i,j,k-1)+
     &                   4*f(i,j,k-2)-f(i,j,k-3))/dz/dz
               else if ((chr(i,j,k-1).ne.ex)
     &                 .and.(chr(i,j,k-2).ne.ex)) then
              !    write(*,*) 'df2_int: warning ... first order k=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
                   f_zz=(f(i,j,k-2)-2*f(i,j,k-1)+f(i,j,k))/dz/dz
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (S)'
                     write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,
     &                                                   Nx,Ny,Nz,dz
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         else  !this is the case where (i,j,k+1) is not excised and (i,j,k-1) is excised
               if (chr(i,j,k+2).ne.ex) then
                   f_zz=(f(i,j,k)-2*f(i,j,k+1)+f(i,j,k+2))/dz/dz
               end if
         end if

        else if (k.eq.(Nz-1)) then
         if ((chr(i,j,k+1).ne.ex).and.(chr(i,j,k-1).ne.ex)) then
                   f_zz=(f(i,j,k+1)-2*f(i,j,k)+f(i,j,k-1))/dz/dz
         else if (chr(i,j,k+1).eq.ex) then
               if (.not.extrap
     &            .and.(chr(i,j,k-1).ne.ex)
     &            .and.(chr(i,j,k-2).ne.ex)
     &            .and.(chr(i,j,k-3).ne.ex)
     &            .and.(chr(i,j,k-4).ne.ex)) then
                   f_zz=(3*f(i,j,k)-9*f(i,j,k-1)+
     &                  10*f(i,j,k-2)-5*f(i,j,k-3)+
     &                  f(i,j,k-4))/dz/dz
               else if (chr(i,j,k-1).ne.ex
     &                 .and.chr(i,j,k-2).ne.ex
     &                 .and.chr(i,j,k-3).ne.ex) then
                   f_zz=(2*f(i,j,k)-5*f(i,j,k-1)+
     &                   4*f(i,j,k-2)-f(i,j,k-3))/dz/dz
               else if (chr(i,j,k-1).ne.ex
     &                 .and.chr(i,j,k-2).ne.ex) then
              !    write(*,*) 'df2_int: warning ... first order k=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
                   f_zz=(f(i,j,k-2)-2*f(i,j,k-1)+f(i,j,k))/dz/dz
               else
                if (first) then
                     first=.false.
                     write(*,*) 'df2_int: error in chr stencil (T)'
                     write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,
     &                                                   Nx,Ny,Nz,dz
                     write(*,*) '    (first error only)'
                end if
                   return
               end if
         end if

        else if (k.eq.Nz) then
               if ((.not.extrap)
     &            .and.(chr(i,j,k-1).ne.ex)
     &            .and.(chr(i,j,k-2).ne.ex)
     &            .and.(chr(i,j,k-3).ne.ex)
     &            .and.(chr(i,j,k-4).ne.ex)) then
                   f_zz=(3*f(i,j,k)-9*f(i,j,k-1)+
     &                  10*f(i,j,k-2)-5*f(i,j,k-3)+
     &                  f(i,j,k-4))/dz/dz
               else if (chr(i,j,k-1).ne.ex
     &                 .and.chr(i,j,k-2).ne.ex
     &                 .and.chr(i,j,k-3).ne.ex) then
                   f_zz=(2*f(i,j,k)-5*f(i,j,k-1)+
     &                   4*f(i,j,k-2)-f(i,j,k-3))/dz/dz
               else if (chr(i,j,k-1).ne.ex
     &                  .and.chr(i,j,k-2).ne.ex) then
              !    write(*,*) 'df2_int: warning ... first order k=1'
              !    write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
                   f_zz=(f(i,j,k-2)-2*f(i,j,k-1)+f(i,j,k))/dz/dz
               else
                 if (first) then
                   first=.false.
                   write(*,*) 'df2_int: error in chr stencil (U)'
                   write(*,*) '    i,j,k,Nx,Ny,Nz,dz=',i,j,k,Nx,Ny,Nz,dz
                   write(*,*) '    (first error only)'
                 end if
                   return
               end if
        end if



!!!!!!!!!!!!!!!

!!!!!!!!!!OLDVERSION!!!!!!!!!!!!!!!!
!        !i
!        if (i.eq.1.or.(chr(i-1,j,k).eq.ex)) then
!           if (i.ge.(Nx-1).or.chr(i+1,j,k).eq.ex.or.chr(i+2,j,k).eq.ex) then
!              if (first) then
!                 first=.false.
!                 write(*,*) 'df2_int: error in chr (A)'
!                 write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
!                 write(*,*) '    (first error only)'
!              end if
!              return
!           end if
!           if (i.eq.(Nx-2).or.chr(i+3,j,k).eq.ex) then
!              ! write(*,*) 'df2_int: warning ... first order i=1'
!              ! write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
!              f_xx=(f(i+2,j,k)-2*f(i+1,j,k)+f(i,j,k))/dx/dx 
!           else if ((.not.extrap).and.i.lt.(Nx-3)
!     &               .and.chr(i+4,j,k).ne.ex) then
!              f_xx=(3*f(i,j,k)-9*f(i+1,j,k)+
!     &             10*f(i+2,j,k)-5*f(i+3,j,k)+
!     &             f(i+4,j,k))/dx/dx
!           else
!              f_xx=(2*f(i,j,k)-5*f(i+1,j,k)+
!     &              4*f(i+2,j,k)-f(i+3,j,k))/dx/dx
!           end if
!           call df1_int_y(f,f_y_ip1,x,y,i+1,j,chr,ex,Nx,Ny,Nz)
!           call df1_int_y(f,f_y_ip2,x,y,i+2,j,chr,ex,Nx,Ny,Nz)
!           f_xy=(-3*f_y+4*f_y_ip1-f_y_ip2)/2/dx
!        else if (i.eq.Nx.or.(chr(i+1,j,k).eq.ex)) then
!           if (i.le.2.or.
!     &        chr(i-1,j,k).eq.ex.or.chr(i-2,j,k).eq.ex) then
!              if (first) then
!                 first=.false.
!                 write(*,*) 'df2_int: error in chr (B)'
!                 write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
!                 write(*,*) '    (first error only)'
!              end if
!              return
!           end if
!           if (i.eq.3.or.chr(i-3,j,k).eq.ex) then
!              ! write(*,*) 'df2_int: warning ... first order i=Nx'
!              ! write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
!              f_xx=(f(i,j,k)-2*f(i-1,j,k)+f(i-2,j,k))/dx/dx 
!           else if ((.not.extrap).and.i.gt.4.and.chr(i-4,j,k).ne.ex) then
!              f_xx=(3*f(i,j,k)-9*f(i-1,j,k)+
!     &              10*f(i-2,j,k)-5*f(i-3,j,k)+
!     &              f(i-4,j,k))/dx/dx
!           else
!              f_xx=(2*f(i,j,k)-5*f(i-1,j,k)+
!     &              4*f(i-2,j,k)-f(i-3,j,k))/dx/dx
!           end if
!           call df1_int_y(f,f_y_im1,x,y,i-1,j,chr,ex,Nx,Ny,Nz)
!           call df1_int_y(f,f_y_im2,x,y,i-2,j,chr,ex,Nx,Ny,Nz)
!           f_xy=(3*f_y-4*f_y_im1+f_y_im2)/2/dx
!        else if (chr(i+1,j,k).ne.ex.and.chr(i-1,j,k).ne.ex) then
!           f_xx=(f(i+1,j,k)-2*f(i,j,k)+f(i-1,j,k))/dx/dx 
!           call df1_int_y(f,f_y_im1,x,y,i-1,j,chr,ex,Nx,Ny,Nz)
!           call df1_int_y(f,f_y_ip1,x,y,i+1,j,chr,ex,Nx,Ny,Nz)
!           f_xy=(f_y_ip1-f_y_im1)/2/dx
!        else
!           if (first) then
!              first=.false.
!              write(*,*) 'df2_int: error in chr (C)'
!              write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
!              write(*,*) '    (first error only)'
!           end if
!           return
!        end if
!        !j
!        if (j.eq.1.or.(chr(i,j-1,k).eq.ex)) then
!           if (j.ge.(Ny-1).or.chr(i,j+1,k).eq.ex.or.chr(i,j+2,k).eq.ex) then
!              if (first) then
!                 first=.false.
!                 write(*,*) 'df2_int: error in chr (D)'
!                 write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
!                 write(*,*) '    (first error only)'
!              end if
!              return
!           end if
!           if (j.eq.(Ny-2).or.chr(i,j+3,k).eq.ex) then
!              !write(*,*) 'df2_int: warning ... first order j=1'
!              !write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
!              f_yy=(f(i,j+2,k)-2*f(i,j+1,k)+f(i,j,k))/dy/dy 
!           else if ((.not.extrap).and.j.lt.(Ny-3).and.
!     &               chr(i,j+4,k).ne.ex) then
!              f_yy=(3*f(i,j,k)-9*f(i,j+1,k)+
!     &              10*f(i,j+2,k)-5*f(i,j+3,k)+
!     &              f(i,j+4,k))/dy/dy
!           else
!              f_yy=(2*f(i,j,k)-5*f(i,j+1,k)+
!     &              4*f(i,j+2,k)-f(i,j+3,k))/dy/dy
!           end if
!        else if (j.eq.Ny.or.(chr(i,j+1,k).eq.ex)) then
!           if (j.le.2.or.chr(i,j-1,k).eq.ex.or.chr(i,j-2,k).eq.ex) then
!              if (first) then
!                 first=.false.
!                 write(*,*) 'df2_int: error in chr (E)'
!                 write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
!                 write(*,*) '    name=',name
!                 write(*,*) '    (first error only)'
!              end if
!              return
!           end if
!           if (j.eq.3.or.chr(i,j-3,k).eq.ex) then
!              !write(*,*) 'df2_int: warning ... first order j=Ny'
!              !write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
!              f_yy=(f(i,j,k)-2*f(i,j-1,k)+f(i,j-2,k))/dy/dy 
!           else if ((.not.extrap).and.j.gt.4.and.chr(i,j-4,k).ne.ex) then
!              f_yy=(3*f(i,j,k)-9*f(i,j-1,k)+
!     &              10*f(i,j-2,k)-5*f(i,j-3,k)+
!     &              f(i,j-4,k))/dy/dy
!           else
!              f_yy=(2*f(i,j,k)-5*f(i,j-1,k)+
!     &              4*f(i,j-2,k)-f(i,j-3,k))/dy/dy
!           end if
!        else if (chr(i,j+1,k).ne.ex.and.chr(i,j-1,k).ne.ex) then
!            f_yy=(f(i,j+1,k)-2*f(i,j,k)+f(i,j-1,k))/dy/dy 
!        else
!            if (first) then
!               first=.false.
!               write(*,*) 'df2_int: error in chr (F)'
!               write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
!               write(*,*) '    (first error only)'
!            end if
!            return
!        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        if (ltrace) then
           write(*,*) 'df2_int for ',name
           write(*,*) ' f_tt=',f_tt
           write(*,*) ' f_tx=',f_tx
           write(*,*) ' f_ty=',f_ty
           write(*,*) ' f_tz=',f_tz
           write(*,*) ' f_xx=',f_xx
           write(*,*) ' f_xy=',f_xy
           write(*,*) ' f_xz=',f_xz
           write(*,*) ' f_yy=',f_yy
           write(*,*) ' f_yz=',f_yz
           write(*,*) ' f_zz=',f_zz
        end if

        return
        end

c----------------------------------------------------------------------
c Background AdS values ... with these new variables
c the AdS metric has been factored into the maple already
c----------------------------------------------------------------------
        subroutine init_ghb_ads(gb_tt,gb_tx,gb_ty,
     &                      gb_tz,
     &                      gb_xx,gb_xy,
     &                      gb_xz,
     &                      gb_yy,
     &                      gb_yz,
     &                      gb_zz,Hb_t,Hb_x,Hb_y,
     &                      Hb_z,
     &                      L,x,y,z,chr,ex,Nx,Ny,Nz,regtype)
        implicit none
        integer Nx,Ny,Nz
        integer regtype
        real*8 gb_tt(Nx,Ny,Nz),gb_tx(Nx,Ny,Nz),gb_ty(Nx,Ny,Nz)
        real*8 gb_tz(Nx,Ny,Nz)
        real*8 gb_xx(Nx,Ny,Nz),gb_xy(Nx,Ny,Nz)
        real*8 gb_xz(Nx,Ny,Nz)
        real*8 gb_yy(Nx,Ny,Nz)
        real*8 gb_yz(Nx,Ny,Nz)
        real*8 gb_zz(Nx,Ny,Nz),tfunction(Nx,Ny,Nz),Hb_t(Nx,Ny,Nz)
        real*8 Hb_x(Nx,Ny,Nz),Hb_y(Nx,Ny,Nz)
        real*8 Hb_z(Nx,Ny,Nz)
        real*8 chr(Nx,Ny,Nz),ex,L,x(Nx),y(Ny),z(Nz)

        integer i,j,k

        logical ltrace
        parameter (ltrace=.false.)
 
        ! initialize fixed-size variables
        data i,j,k/0,0,0/
  
        !--------------------------------------------------------------
 
        do i=1,Nx
           do j=1,Ny
            do k=1,Nz
              gb_tt(i,j,k)=0
              gb_tx(i,j,k)=0
              gb_ty(i,j,k)=0
              gb_tz(i,j,k)=0
              gb_xx(i,j,k)=0
              gb_xy(i,j,k)=0
              gb_xz(i,j,k)=0
              gb_yy(i,j,k)=0
              gb_yz(i,j,k)=0
              gb_zz(i,j,k)=0

              Hb_t(i,j,k)=0
              Hb_x(i,j,k)=0
              Hb_y(i,j,k)=0
              Hb_z(i,j,k)=0
            end do
           end do
        end do

!        call axi_reg_g(gb_tt,gb_tx,gb_ty,
!     &                 gb_xx,gb_xy,gb_yy,psi,tfunction,chr,ex,
!     &                 L,x,y,z,Nx,Ny,Nz,regtype)

!        call axi_reg_Hb(Hb_t,Hb_x,Hb_y,chr,ex,L,x,y,z,Nx,Ny,Nz,regtype)

        return
        end

c----------------------------------------------------------------------
c initializes f with a 2d gaussian-like profile ... multiplication
c by (1-rho^2)^3 is for correct asymptotics for AdS 
c
c f = A * (1-rho^2) exp (- (r-r0)^2/delta^2) , r>r0
c   = A , r < r0
c
c where r = sqrt ( (1-ex^2)*(x)^2 + (1-ey^2)*(y)^2 )
c----------------------------------------------------------------------
        subroutine gauss3d(f,A,B,C,r0,delta,xu0,yu0,zu0,ex,ey,ez,
     &                     L,x,y,z,Nx,Ny,Nz,
     &                     rhoc,rhod,stype)
        implicit none
        integer Nx,Ny,Nz
        real*8 f(Nx,Ny,Nz),x(Nx),y(Ny),z(Nz)
        real*8 A,B,C,r0,delta,ex,ey,ez,xu0,yu0,zu0,L

        integer i,j,k
        integer stype
        real*8 r,x0,y0,z0,rho0,xi0,chi0,csr,xb,yb,zb

        real*8 rhoc,rhod
        real*8 f1,trans

        real*8 PI
        parameter (PI=3.141592653589793d0)

        ! initialize fixed-size variables
        data i,j,k/0,0,0/
        data r,x0,y0,z0,rho0,xi0,chi0,csr,xb,yb,zb
     &       /0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
  
        !--------------------------------------------------------------
 

        do i=1,Nx
           do j=1,Ny
            do k=1,Nz
              f(i,j,k)=0
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

              xb=x0-xu0
              yb=y0-yu0
              zb=z0-zu0
              r=sqrt(xb**2+yb**2+zb**2
     &           -ex**2*xb**2-ey**2*yb**2-ez**2*zb**2)

              f1=trans(rho0,rhoc,rhod)

              ! Gaussian phi=amp*exp(-((r-r0)/delta)^2) profile 
              ! remember that phi=phi1*(1-rho^2)^2
              ! NOTE: compactified coordinate here is: rho

              if (stype.eq.0) then

                 if (rho0.ge.1) then
                    f(i,j,k)=0
                 else if (r.gt.r0) then
                   f(i,j,k)=
     &             A*exp(-((r-r0)/delta)**2)
     &            +B*cos(chi0)*4*f1*(1-f1)
     &            +C*cos(xi0)*4*f1*(1-f1)
                 else
                  if (rho0.ne.0.0d0) then
                   f(i,j,k)=
     &             A
     &            +B*cos(chi0)*4*f1*(1-f1)
     &            +C*cos(xi0)*4*f1*(1-f1)
                  else if ((y0.ne.0.0d0).or.(z0.ne.0.0d0)) then
                   f(i,j,k)=
     &             A
     &            +B*cos(chi0)*4*f1*(1-f1)
                  else
                   f(i,j,k)=
     &             A
                  end if
                 end if

              else if (stype.eq.1) then

                 if (rho0.ge.1) then
                    f(i,j,k)=0
                 else if (rho0.gt.r0) then
                    f(i,j,k)=
     &              A*(1-f1)/(1-rho0**2)**2
     &             +B*cos(chi0)*4*f1*(1-f1)/(1-rho0**2)**2
     &             +C*cos(xi0)*4*f1*(1-f1)/(1-rho0**2)**2
                 else
                  if (rho0.ne.0.0d0) then
                    f(i,j,k)=
     &              A*(1-f1)/(1-rho0**2)**2
     &             +B*cos(chi0)*4*f1*(1-f1)/(1-rho0**2)**2
     &             +C*cos(xi0)*4*f1*(1-f1)/(1-rho0**2)**2
                  else if ((y0.ne.0.0d0).or.(z0.ne.0.0d0)) then
                    f(i,j,k)=
     &              A*(1-f1)/(1-rho0**2)**2
     &             +B*cos(chi0)*4*f1*(1-f1)/(1-rho0**2)**2
                  else
                    f(i,j,k)=
     &              A*(1-f1)/(1-rho0**2)**2
                  end if
                 end if

              else if (stype.eq.2) then !annulus ID

                 if (rho0.ge.1) then
                    f(i,j,k)=0
                 else if (rho0.gt.r0) then
                    f(i,j,k)=
     &              A*(1-f1)/(1-rho0**2)**2
     &             +B*cos(chi0)*4*f1*(1-f1)/(1-rho0**2)**2
     &             +C*cos(xi0)*4*f1*(1-f1)/(1-rho0**2)**2
                 else
                    f(i,j,k)=
     &              A*(1-f1)/(1-rho0**2)**2
     &             +B*cos(chi0)*4*f1*(1-f1)/(1-rho0**2)**2
     &             +C*cos(xi0)*4*f1*(1-f1)/(1-rho0**2)**2
                 end if

              end if
            end do
           end do
        end do

        return
        end

c----------------------------------------------------------------------
c initializes complex s.f. phi with approximate q-ball profile
c
c phi = amp * exp (-rho^2/delta^2 + i*omega*t) * (1-rho^2)
c
c phi = d(phi)/dt = i*omega*phi
c
c where rgo = sqrt ((x-x0)^2 + (y-y0)^2)
c
c and profile is multiplied by (1-rho^2)^3 for correct asymptotics
c----------------------------------------------------------------------
!        subroutine approx_qb(phi_r,phi_i,phi_r_dot,phi_i_dot,amp,xu0,
!     &                       delta,omega,v0,L,x,y,z,Nx,Ny,Nz)
!        implicit none
!        integer Nx,Ny,Nz
!        real*8 phi_r(Nx,Ny,Nz),phi_i(Nx,Ny,Nz)
!        real*8 phi_r_dot(Nx,Ny,Nz),phi_i_dot(Nx,Ny,Nz)
!        real*8 x(Nx),y(Ny),z(Nz),L,v0(2)
!        real*8 amp,omega,delta,xu0(2)
!
!        integer i,j,k
!        real*8 r,x0,y0,z0,rho0,xb,yb
!
!        ! initialize fixed-size variables
!        data i,j,k/0,0,0/
!        data r,x0,y0,z0,rho0,xb,yb/0.0,0.0,0.0,0.0,0.0,0.0,0.0/
!  
!        !--------------------------------------------------------------
! 
!        do i=1,Nx
!           do j=1,Ny
!            do k=1,Nz
!              x0=x(i)
!              y0=y(j)
!              z0=z(k)
!              rho0=sqrt(x0**2+y0**2)
!
!              xb=x0-xu0(1)
!              yb=y0-xu0(2)
!              r=sqrt(xb**2+yb**2)
!
!              if (rho0.ge.1) then
!                 phi_r(i,j,k)=0
!                 phi_i(i,j,k)=0
!                 phi_r_dot(i,j,k)=0
!                 phi_i_dot(i,j,k)=0
!              else
!                 phi_r(i,j,k)=amp*exp(-(r/delta)**2)*(1-rho0**2)**2
!                 phi_i(i,j,k)=0
!                 phi_r_dot(i,j,k)=0
!                 phi_i_dot(i,j,k)=omega*phi_r(i,j,k)
!              end if
!           end do
!          end do
!        end do
!
!        return
!        end
!
c-----------------------------------------------------------------------
c for variables with lin_zero_bnd ... zeros residual there
c-----------------------------------------------------------------------
        subroutine lin_zero_bnd_res(f,phys_bdy,all,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz,all
        real*8 f(Nx,Ny,Nz)
        integer phys_bdy(6)

        integer i,j,k,is,ie,js,je,ks,ke

        ! initialize fixed-size variables
        data i,j,k,is,ie,js,je,ks,ke/0,0,0,0,0,0,0,0,0/
  
        !--------------------------------------------------------------
 
        if (phys_bdy(1).eq.1.or.all.eq.1) then
           do j=2,Ny-1
            do k=2,Nz-1
              f(2,j,k)=0
            end do
           end do
        end if

        if (phys_bdy(2).eq.1.or.all.eq.1) then
           do j=2,Ny-1
            do k=2,Nz-1
              f(Nx-1,j,k)=0
            end do
           end do
        end if

        if (phys_bdy(3).eq.1.or.all.eq.1) then
           do i=3,Nx-2
             do k=2,Nz-1
              f(i,2,k)=0
             end do
           end do
        end if

        if (phys_bdy(4).eq.1.or.all.eq.1) then
           do i=3,Nx-2
            do k=2,Nz-1
              f(i,Ny-1,k)=0
            end do
           end do
        end if

        if (phys_bdy(5).eq.1.or.all.eq.1) then
           do i=3,Nx-2
            do j=3,Ny-2
              f(i,j,2)=0
            end do
           end do
        end if

        if (phys_bdy(6).eq.1.or.all.eq.1) then
           do i=3,Nx-2
            do j=3,Ny-2
              f(i,j,Nz-1)=0
            end do
           end do
        end if

        return
        end

c---------------------------------------------------------------------------
c applies a simple KO filter to f and modified
c near excision boundaries as follows:
c (or all boundaries if do_bdy=1, 
c  or neither if do_bdy=0 and internal flag no_uc_ex_bdy=.true.)
c
c with undivided operators Up(f) = f(i+1)-f(i)
c                          Um(f) = f(i)-f(i-1)
c
c Interior: eps*w4*(Up Um)^2
c left+1  : eps*w3*(Up Up Um)
c left    : eps*w2*(Up Um)
c right-1 : eps*w3*(Um Um Up)
c right   : eps*w2*(Up Um)
c
c NOTE: SPECIAL FLAG ... if do_ex>0, eps is treated as a scalar, 
c       else eps is treated as an array.
c
c same as dmdiss3d_ex(), but with 'kmax' flag, applying separate
c sweeps from kmax.eq.1
c------------------------------------------------------------------------
!        subroutine dmdiss3d_ex_gen(f,work,eps,do_bdy,phys_bdy_type,even,
!     &             odd,nx,ny,nz,chr,ex,do_ex,ind_sweeps,kmax)
!        implicit none
!        integer nx,ny,nz,do_bdy,phys_bdy_type(6),even,odd,do_ex,kmax
!        real*8 f(nx,ny,nz),work(nx,ny,nz),eps(nx,ny,nz)
!        real*8 chr(nx,ny,nz),ex
!        logical ind_sweeps
!
!        integer i,j,k,bo1,bo2
!        integer pass,npass,kd
!        real*8 eps_eff,f_hf,norm_f
!        real*8 w4,w3,w2
!        parameter (w4=1.0d0/16.0d0,w3=1.0d0/8.0d0,w2=1.0d0/4.0d0)
!        logical no_uc_ex_diss
!        parameter (no_uc_ex_diss=.true.)
!
!        ! initialize fixed-size variables
!        data i,j,k,bo1,bo2/0,0,0,0,0/
!        data pass,npass,kd/0,0,0/
!        data eps_eff,f_hf,norm_f/0.0,0.0,0.0/
!
!        !--------------------------------------------------------------
!        
!        eps_eff=eps(1,1,1)
!
!        if (do_bdy.eq.0) then
!           bo1=0
!           bo2=0
!           if (no_uc_ex_diss) then
!              bo1=-max(nx,ny,nz)
!              bo2=bo1
!           end if
!        else
!           bo1=1
!           bo2=2
!        end if
!
!        do kd=kmax,1,-1
!
!        npass=1
!        if (ny.gt.1.and.ind_sweeps) npass=2
!        if (nz.gt.1.and.ind_sweeps) npass=3
!
!        do pass=1,npass
!
!         do i=1,nx
!           do j=1,ny
!              do k=1,nz
!                 work(i,j,k)=f(i,j,k)
!              end do
!           end do
!         end do
!
!         do i=1,nx
!          do j=1,ny
!           do k=1,nz
!             if ((chr(i,j,k).ne.ex)) then
!
!               if (do_ex.lt.0) eps_eff=eps(i,j,k)
!
!               if (.not.ind_sweeps.or.pass.eq.1) then
!                f_hf=0
!                if (i.ge.(2*kd+1).and.i.le.(nx-2*kd).and.
!     &           ((chr(i-2*kd,j,k).ne.ex.and.chr(i-kd,j,k).ne.ex
!     &             .and.chr(i+2*kd,j,k).ne.ex.and.chr(i+kd,j,k).ne.ex)))
!     &          then
!                   f_hf=w4*(
!     &                 work(i-2*kd,j,k)+work(i+2*kd,j,k)
!     &                 -4*(work(i-kd,j,k)+work(i+kd,j,k))+6*work(i,j,k))
!                else if (i.ge.(kd+1).and.i.lt.(2*kd+1).and.
!     &                   phys_bdy_type(1).eq.odd.and.
!     &                ((chr(i-kd,j,k).ne.ex.and.
!     &                 chr(i+2*kd,j,k).ne.ex.and.chr(i+kd,j,k).ne.ex)) )
!     &          then
!                   f_hf=w4*(
!     &                (-work(i,j,k))+work(i+2*kd,j,k)
!     &                 -4*(work(i-kd,j,k)+work(i+kd,j,k))+6*work(i,j,k))
!                else if (i.lt.(kd+1).and.
!     &                   phys_bdy_type(1).eq.odd.and.
!     &               ((chr(i+2*kd,j,k).ne.ex.and.chr(i+kd,j,k).ne.ex)) )
!     &          then
!                   f_hf=w4*(
!     &                (-work(i+2*kd,j,k))+work(i+2*kd,j,k)
!     &          -4*((-work(i+1*kd,j,k))+work(i+1*kd,j,k))+6*work(i,j,k))
!                else if (i.gt.(Nx-2*kd).and.i.le.(Nx-kd)
!     &                   .and.phys_bdy_type(2).eq.odd.and.
!     &              ((chr(i+kd,j,k).ne.ex.and.
!     &                chr(i-2*kd,j,k).ne.ex.and.chr(i-kd,j,k).ne.ex)) )
!     &          then
!                   f_hf=w4*(
!     &                 work(i-2*kd,j,k)+(-work(i,j,k))
!     &                 -4*(work(i-kd,j,k)+work(i+kd,j,k))+6*work(i,j,k))
!                else if (i.gt.(Nx-kd).and.phys_bdy_type(2).eq.odd.and.
!     &             ((chr(i-2*kd,j,k).ne.ex.and.chr(i-kd,j,k).ne.ex)) )
!     &          then
!                   f_hf=w4*(
!     &                  work(i-2*kd,j,k)+(-work(i-2*kd,j,k))
!     &              -4*(work(i-kd,j,k)+(-work(i-kd,j,k)))+6*work(i,j,k))
!                else if (i.ge.(kd+1).and.i.lt.(2*kd+1)
!     &                   .and.phys_bdy_type(1).eq.even.and.
!     &              ((chr(i-kd,j,k).ne.ex.and.
!     &                 chr(i+2*kd,j,k).ne.ex.and.chr(i+kd,j,k).ne.ex)) )
!     &          then
!                   f_hf=w4*(
!     &                    (work(i,j,k))+work(i+2*kd,j,k)
!     &                 -4*(work(i-kd,j,k)+work(i+kd,j,k))+6*work(i,j,k))
!                else if (i.lt.(kd+1).and.phys_bdy_type(1).eq.even.and.
!     &               ((chr(i+2*kd,j,k).ne.ex.and.chr(i+kd,j,k).ne.ex)) )
!     &          then
!                   f_hf=w4*(
!     &                   (work(i+2*kd,j,k))+work(i+2*kd,j,k)
!     &               -4*((work(i+kd,j,k))+work(i+kd,j,k))+6*work(i,j,k))
!                else if (i.gt.(Nx-2*kd).and.i.le.(Nx-kd)
!     *                   .and.phys_bdy_type(2).eq.even.and.
!     &              ((chr(i+kd,j,k).ne.ex.and.
!     &                 chr(i-2*kd,j,k).ne.ex.and.chr(i-kd,j,k).ne.ex)) )
!     &          then
!                   f_hf=w4*(
!     &                     work(i-2*kd,j,k)+(work(i,j,k))
!     &                 -4*(work(i-kd,j,k)+work(i+kd,j,k))+6*work(i,j,k))
!                else if (i.gt.(Nx-kd).and.phys_bdy_type(2).eq.even.and.
!     &             ((chr(i-2*kd,j,k).ne.ex.and.chr(i-kd,j,k).ne.ex)) )
!     &          then
!                   f_hf=w4*(
!     &                   work(i-2*kd,j,k)+(work(i-2*kd,j,k))
!     &               -4*(work(i-kd,j,k)+(work(i-kd,j,k)))+6*work(i,j,k))
!                else if (i.gt.(2-bo1)*kd.and.i.le.(nx-2*kd).and.
!     &                 (chr(i-kd,j,k).ne.ex.and.chr(i+kd,j,k).ne.ex.and.
!     &                  chr(i+2*kd,j,k).ne.ex)) 
!     &          then
!                   f_hf=w3*(
!     &                    -work(i-kd,j,k)+3*work(i,j,k)
!     &                  -3*work(i+kd,j,k)+work(i+2*kd,j,k))
!                else if (i.gt.2*kd.and.i.le.(nx-2*kd+bo1*kd).and.
!     &                 (chr(i-kd,j,k).ne.ex.and.chr(i+kd,j,k).ne.ex.and.
!     &                   chr(i-2*kd,j,k).ne.ex) ) 
!     &          then
!                   f_hf=w3*(
!     &                    -work(i+kd,j,k)+3*work(i,j,k)
!     &                  -3*work(i-kd,j,k)+work(i-2*kd,j,k))
!                else if ((i.gt.(2-bo2)*kd.or.   
!     &           (i.le.2*kd.and.chr(max(1,i-kd),j,k).eq.ex))
!     &                       .and.i.le.(nx-2*kd).and.
!     &                (chr(i+kd,j,k).ne.ex.and.chr(i+2*kd,j,k).ne.ex) )
!     &          then
!                   f_hf=w2*(
!     &                   work(i,j,k)-2*work(i+kd,j,k)+work(i+2*kd,j,k))
!                else if (i.gt.2*kd.and.(i.lt.(nx-kd+bo2*kd).or.
!     &              (i.ge.(nx-kd).and.chr(min(nx,i+kd),j,k).eq.ex)).and.
!     &                 (chr(i-kd,j,k).ne.ex.and.chr(i-2*kd,j,k).ne.ex) )
!     &          then
!                   f_hf=w2*(
!     &                    work(i,j,k)-2*work(i-kd,j,k)+work(i-2*kd,j,k))
!                end if
!
!                f(i,j,k)=f(i,j,k)-eps_eff*f_hf
!               end if
!
!               if (.not.ind_sweeps.or.pass.eq.2) then
!                f_hf=0
!                if (j.ge.(2*kd+1).and.j.le.(ny-2*kd).and.
!     &           ((chr(i,j-2*kd,k).ne.ex.and.chr(i,j-kd,k).ne.ex
!     &             .and.chr(i,j+2*kd,k).ne.ex.and.chr(i,j+kd,k).ne.ex)))
!     &          then
!                   f_hf=w4*(
!     &                 work(i,j-2*kd,k)+work(i,j+2*kd,k)
!     &                 -4*(work(i,j-kd,k)+work(i,j+kd,k))+6*work(i,j,k))
!                else if (j.ge.(kd+1).and.j.lt.(2*kd+1).and.
!     &                   phys_bdy_type(3).eq.odd.and.
!     &                ((chr(i,j-kd,k).ne.ex.and.
!     &                 chr(i,j+2*kd,k).ne.ex.and.chr(i,j+kd,k).ne.ex)) )
!     &          then
!                   f_hf=w4*(
!     &                (-work(i,j,k))+work(i,j+2*kd,k)
!     &                 -4*(work(i,j-kd,k)+work(i,j+kd,k))+6*work(i,j,k))
!                else if (j.lt.(kd+1).and.
!     &                   phys_bdy_type(3).eq.odd.and.
!     &               ((chr(i,j+2*kd,k).ne.ex.and.chr(i,j+kd,k).ne.ex)) )
!     &          then
!                   f_hf=w4*(
!     &                (-work(i,j+2*kd,k))+work(i,j+2*kd,k)
!     &          -4*((-work(i,j+1*kd,k))+work(i,j+1*kd,k))+6*work(i,j,k))
!                else if (j.gt.(ny-2*kd).and.j.le.(ny-kd)
!     &                   .and.phys_bdy_type(4).eq.odd.and.
!     &              ((chr(i,j+kd,k).ne.ex.and.
!     &                chr(i,j-2*kd,k).ne.ex.and.chr(i,j-kd,k).ne.ex)) )
!     &          then
!                   f_hf=w4*(
!     &                 work(i,j-2*kd,k)+(-work(i,j,k))
!     &                 -4*(work(i,j-kd,k)+work(i,j+kd,k))+6*work(i,j,k))
!                else if (j.gt.(Ny-kd).and.phys_bdy_type(4).eq.odd.and.
!     &             ((chr(i,j-2*kd,k).ne.ex.and.chr(i,j-kd,k).ne.ex)) )
!     &          then
!                   f_hf=w4*(
!     &                  work(i,j-2*kd,k)+(-work(i,j-2*kd,k))
!     &              -4*(work(i,j-kd,k)+(-work(i,j-kd,k)))+6*work(i,j,k))
!                else if (j.ge.(kd+1).and.j.lt.(2*kd+1)
!     &                   .and.phys_bdy_type(3).eq.even.and.
!     &              ((chr(i,j-kd,k).ne.ex.and.
!     &                 chr(i,j+2*kd,k).ne.ex.and.chr(i,j+kd,k).ne.ex)) )
!     &          then
!                   f_hf=w4*(
!     &                    (work(i,j,k))+work(i,j+2*kd,k)
!     &                 -4*(work(i,j-kd,k)+work(i,j+kd,k))+6*work(i,j,k))
!                else if (j.lt.(kd+1).and.phys_bdy_type(3).eq.even.and.
!     &               ((chr(i,j+2*kd,k).ne.ex.and.chr(i,j+kd,k).ne.ex)) )
!     &          then
!                   f_hf=w4*(
!     &                   (work(i,j+2*kd,k))+work(i,j+2*kd,k)
!     &               -4*((work(i,j+kd,k))+work(i,j+kd,k))+6*work(i,j,k))
!                else if (j.gt.(Ny-2*kd).and.j.le.(Ny-kd)
!     *                   .and.phys_bdy_type(4).eq.even.and.
!     &              ((chr(i,j+kd,k).ne.ex.and.
!     &                 chr(i,j-2*kd,k).ne.ex.and.chr(i,j-kd,k).ne.ex)) )
!     &          then
!                   f_hf=w4*(
!     &                     work(i,j-2*kd,k)+(work(i,j,k))
!     &                 -4*(work(i,j-kd,k)+work(i,j+kd,k))+6*work(i,j,k))
!                else if (j.gt.(Ny-kd).and.phys_bdy_type(4).eq.even.and.
!     &             ((chr(i,j-2*kd,k).ne.ex.and.chr(i,j-kd,k).ne.ex)) )
!     &          then
!                   f_hf=w4*(
!     &                   work(i,j-2*kd,k)+(work(i,j-2*kd,k))
!     &               -4*(work(i,j-kd,k)+(work(i,j-kd,k)))+6*work(i,j,k))
!                else if (j.gt.(2-bo1)*kd.and.j.le.(ny-2*kd).and.
!     &                 (chr(i,j-kd,k).ne.ex.and.chr(i,j+kd,k).ne.ex.and.
!     &                  chr(i,j+2*kd,k).ne.ex)) 
!     &          then
!                   f_hf=w3*(
!     &                    -work(i,j-kd,k)+3*work(i,j,k)
!     &                  -3*work(i,j+kd,k)+work(i,j+2*kd,k))
!                else if (j.gt.2*kd.and.j.le.(ny-2*kd+bo1*kd).and.
!     &                 (chr(i,j-kd,k).ne.ex.and.chr(i,j+kd,k).ne.ex.and.
!     &                   chr(i,j-2*kd,k).ne.ex) ) 
!     &          then
!                   f_hf=w3*(
!     &                    -work(i,j+kd,k)+3*work(i,j,k)
!     &                  -3*work(i,j-kd,k)+work(i,j-2*kd,k))
!                else if ((j.gt.(2-bo2)*kd.or.   
!     &           (j.le.2*kd.and.chr(i,max(1,j-kd),k).eq.ex))
!     &                       .and.j.le.(ny-2*kd).and.
!     &                (chr(i,j+kd,k).ne.ex.and.chr(i,j+2*kd,k).ne.ex) )
!     &          then
!                   f_hf=w2*(
!     &                   work(i,j,k)-2*work(i,j+kd,k)+work(i,j+2*kd,k))
!                else if (j.gt.2*kd.and.(j.lt.(ny-kd+bo2*kd).or.
!     &              (j.ge.(ny-kd).and.chr(i,min(ny,j+kd),k).eq.ex)).and.
!     &                 (chr(i,j-kd,k).ne.ex.and.chr(i,j-2*kd,k).ne.ex) )
!     &          then
!                   f_hf=w2*(
!     &                    work(i,j,k)-2*work(i,j-kd,k)+work(i,j-2*kd,k))
!                end if
!
!                f(i,j,k)=f(i,j,k)-eps_eff*f_hf
!               end if
!
!               if ((.not.ind_sweeps.or.pass.eq.3).and.Nz.gt.1) then
!                f_hf=0
!                write(*,*) 'dmdiss3d_ex_gen: not yet updated for 3D'
!                stop
!               end if
!
!             end if
!            end do
!           end do
!         end do
!        end do
!
!        end do !kd
!
!        return
!        end

c----------------------------------------------------------------------
c initializes the metric to an exact Schwarzschild-AdS black hole solution
c with radius parameter r0=2*BHmass
c----------------------------------------------------------------------
        subroutine init_schwads4d_bh(r0,L,gb_tt,gb_tx,gb_ty,
     &                         gb_tz,
     &                         gb_xx,gb_xy,
     &                         gb_xz,
     &                         gb_yy,
     &                         gb_yz,
     &                         gb_zz,gb_tt_t,gb_tx_t,gb_ty_t,
     &                         gb_tz_t,
     &                         gb_xx_t,gb_xy_t,
     &                         gb_xz_t,
     &                         gb_yy_t,
     &                         gb_yz_t,
     &                         gb_zz_t,Hb_t,Hb_x,Hb_y,
     &                         Hb_z,
     &                         Hb_t_t,Hb_x_t,Hb_y_t,
     &                         Hb_z_t,
     &                         phys_bdy,
     &                         x,y,z,dt,chr,ex,Nx,Ny,Nz,regtype)
        implicit none

        integer Nx,Ny,Nz
        integer regtype
        integer phys_bdy(6)
        real*8 r0
        real*8 dt,ex,L
        real*8 chr(Nx,Ny,Nz)
        real*8 Hb_t(Nx,Ny,Nz),Hb_x(Nx,Ny,Nz)
        real*8 Hb_y(Nx,Ny,Nz)
        real*8 Hb_z(Nx,Ny,Nz)
        real*8 Hb_t_t(Nx,Ny,Nz),Hb_x_t(Nx,Ny,Nz)
        real*8 Hb_y_t(Nx,Ny,Nz)
        real*8 Hb_z_t(Nx,Ny,Nz)
        real*8 gb_tt(Nx,Ny,Nz),gb_tx(Nx,Ny,Nz)
        real*8 gb_ty(Nx,Ny,Nz)
        real*8 gb_tz(Nx,Ny,Nz)
        real*8 gb_xx(Nx,Ny,Nz),gb_xy(Nx,Ny,Nz)
        real*8 gb_xz(Nx,Ny,Nz)
        real*8 gb_yy(Nx,Ny,Nz),gb_zz(Nx,Ny,Nz)
        real*8 gb_yz(Nx,Ny,Nz)
        real*8 gb_tt_t(Nx,Ny,Nz),gb_tx_t(Nx,Ny,Nz)
        real*8 gb_ty_t(Nx,Ny,Nz),gb_zz_t(Nx,Ny,Nz)
        real*8 gb_tz_t(Nx,Ny,Nz)
        real*8 gb_xx_t(Nx,Ny,Nz),gb_xy_t(Nx,Ny,Nz)
        real*8 gb_xz_t(Nx,Ny,Nz)
        real*8 gb_yy_t(Nx,Ny,Nz)
        real*8 gb_yz_t(Nx,Ny,Nz)
        real*8 x(Nx),y(Ny),z(Nz)

        integer n
        parameter (n=3)

        integer i,j,k

        real*8 rho0,f1,f0,cF0,C0,A0,B0,D0
        real*8 x0,y0,z0
        real*8 r_h,rho_h
        real*8 small
        parameter (small=1d-10)

        real*8 tfunction(Nx,Ny,Nz)

        ! initialize fixed-size variables
        data i,j,k/0,0,0/

        data rho0,f1,f0,cF0/0.0,0.0,0.0,0.0/
        data C0,A0,B0,D0/0.0,0.0,0.0,0.0/
        data x0,y0,z0/0.0,0.0,0.0/
        data r_h,rho_h/0.0,0.0/

        !--------------------------------------------------------------

        ! compute horizon global radius r_h and corresponding compactified rho_h
        !(NOTE ... here the x0,y0,rho0 are compactified (CODE)
        ! versions of the coordinates)

        r_h=-L**2
     &   /(3**(1.0d0/3.0d0))
     &   /((9*L**2*(r0/2)+sqrt(3.0d0)*sqrt(L**6+27*L**4*(r0/2)**2))
     &    **(1.0d0/3.0d0))
     &   +((9*L**2*(r0/2)+sqrt(3.0d0)*sqrt(L**6+27*L**4*(r0/2)**2))
     &    **(1.0d0/3.0d0))
     &   /(3**(2.0d0/3.0d0))

        rho_h=(-1 + sqrt(1 + r_h**2))/r_h

        ! initialize metric 
        do i=1,Nx
           do j=1,Ny
            do k=1,Nz
              gb_tt_t(i,j,k)=0
              gb_tx_t(i,j,k)=0
              gb_ty_t(i,j,k)=0
              gb_tz_t(i,j,k)=0
              gb_xx_t(i,j,k)=0
              gb_xy_t(i,j,k)=0
              gb_xz_t(i,j,k)=0
              gb_yy_t(i,j,k)=0
              gb_yz_t(i,j,k)=0
              gb_zz_t(i,j,k)=0
              Hb_t_t(i,j,k)=0
              Hb_x_t(i,j,k)=0
              Hb_y_t(i,j,k)=0
              Hb_z_t(i,j,k)=0
              if (chr(i,j,k).eq.ex) then
                 gb_tt(i,j,k)=0
                 gb_tx(i,j,k)=0
                 gb_ty(i,j,k)=0
                 gb_tz(i,j,k)=0
                 gb_xx(i,j,k)=0
                 gb_xy(i,j,k)=0
                 gb_xz(i,j,k)=0
                 gb_yy(i,j,k)=0
                 gb_yz(i,j,k)=0
                 gb_zz(i,j,k)=0
                 Hb_t(i,j,k)=0
                 Hb_x(i,j,k)=0
                 Hb_y(i,j,k)=0
                 Hb_z(i,j,k)=0
              else
                 x0=x(i)
                 y0=y(j)
                 z0=z(k)
                 rho0=sqrt(x0**2+y0**2+z0**2)

                 ! EF-like-near-horizon Schwarzschild-like-near-bdy coordinates

!!!2+1 version!!!!!!!!                 
!                 gb_tt(i,j,k)=(r0/2)*(1/rho0-rho0)
!
!                 gb_tx(i,j,k)=2*x0*(1+rho0**2)
!     &                      *((1-rho_h)/(1-rho0))**(-n)
!     &                      /rho0/(1-rho0**2)**2
!                 gb_ty(i,j,k)=2*y0*(1+rho0**2)
!     &                      *((1-rho_h)/(1-rho0))**(-n)
!     &                      /rho0/(1-rho0**2)**2
!                 gb_tz(i,j,k)=0
!                 gb_xx(i,j,k)=4*x0**2*(1+rho0**2)**2
!     &                      *(-1/(1+(-2+4/L**2)*rho0**2+rho0**4)
!     &                      +L**2*rho0*(1-((1-rho_h)/(1-rho0))**(-2*n))
!     &                      /(4*rho0**3+L**2*(1-rho0**2)**2
!     &                       *(rho0+(r0/2)*(-1+rho0**2))))
!     &                      /rho0**2/(1-rho0**2)**2
!                 gb_xy(i,j,k)=4*L**2*x0*y0*(1+rho0**2)**2
!     &                      *(-4*rho0**3+L**2*(1-rho0**2)**2
!     &                      *(-rho0-(r0/2)*(-1+rho0**2)
!     &                      *((1-rho_h)/(1-rho0))**(2*n)))
!     &                      *((1-rho_h)/(1-rho0))**(-2*n)
!     &                      /rho0**2/(1-rho0**2)**2
!     &                      /(4*rho0**2+L**2*(1-rho0**2)**2)
!     &                      /(4*rho0**3+L**2*(1-rho0**2)**2
!     &                       *(rho0+(r0/2)*(-1+rho0**2)))
!                 gb_xz(i,j,k)=0
!                 gb_yy(i,j,k)=4*y0**2*(1+rho0**2)**2
!     &                      *(-1/(1+(-2+4/L**2)*rho0**2+rho0**4)
!     &                      +L**2*rho0*(1-((1-rho_h)/(1-rho0))**(-2*n))
!     &                      /(4*rho0**3+L**2*(1-rho0**2)**2
!     &                       *(rho0+(r0/2)*(-1+rho0**2))))
!     &                      /rho0**2/(1-rho0**2)**2
!                 gb_yz(i,j,k)=0
!                 psi(i,j,k)=0
!!!!!!!!!!!!!!!!!!!!!!!

!!!CHECKED WITH MATHEMATICA!!

                 gb_tt(i,j,k)=(r0/2)*(1/rho0-rho0)
                 gb_tx(i,j,k)=2*x0*(1+rho0**2)
     &                      *((-1+rho_h)/(-1+rho0))**(-n)
     &                      /rho0/(-1+rho0**2)**2
                 gb_ty(i,j,k)=2*y0*(1+rho0**2)
     &                      *((-1+rho_h)/(-1+rho0))**(-n)
     &                      /rho0/(-1+rho0**2)**2
                 gb_tz(i,j,k)=2*z0*(1+rho0**2)
     &                      *((-1+rho_h)/(-1+rho0))**(-n)
     &                      /rho0/(-1+rho0**2)**2
                 gb_xx(i,j,k)=(-((8*(-1+L**2)*(x0**2-y0**2-z0**2)
     &                        +8*rho0**2+4*L**2*(1+rho0**4))
     &                        /(4*rho0**2+L**2*(-1+rho0**2)**2))
     &                        +(4*(y0**2+z0**2
     &                         +(L**2*x0**2*rho0*(1+rho0**2)**2
     &                         *(-1+((-1+rho_h)/(-1+rho0))**(2*n))
     &                         *((-1+rho_h)/(-1+rho0))**(-2*n))
     &                        /(4*rho0**3+L**2*(-1+rho0**2)**2
     &                         *(rho0+(1.0d0/2.0d0)*r0*(-1+rho0**2)))))
     &                         /rho0**2)
     &                          /(-1+rho0**2)**2
                 gb_xy(i,j,k)=(4*L**2*x0*y0*(1+rho0**2)**2
     &                        *(-4*rho0**3+L**2*(-1+rho0**2)**2
     &                        *(-rho0-(1.0d0/2.0d0)*r0*(-1+rho0**2)
     &                        *((-1+rho_h)/(-1+rho0))**(2*n)))
     &                        *((-1+rho_h)/(-1+rho0))**(-2*n))
     &                        /(rho0**2*(-1+rho0**2)**2*(4*rho0**2
     &                        +L**2*(-1+rho0**2)**2)*(4*rho0**3
     &                        +L**2*(-1+rho0**2)**2
     &                        *(rho0+(1.0d0/2.0d0)*r0*(-1+rho0**2))))
                 gb_xz(i,j,k)=(4*L**2*x0*z0*(1+rho0**2)**2
     &                        *(-4*rho0**3+L**2*(-1+rho0**2)**2
     &                        *(-rho0-(1.0d0/2.0d0)*r0*(-1+rho0**2)
     &                        *((-1+rho_h)/(-1+rho0))**(2*n)))
     &                        *((-1+rho_h)/(-1+rho0))**(-2*n))
     &                        /(rho0**2*(-1+rho0**2)**2*(4*rho0**2
     &                        +L**2*(-1+rho0**2)**2)*(4*rho0**3
     &                        +L**2*(-1+rho0**2)**2
     &                        *(rho0+(1.0d0/2.0d0)*r0*(-1+rho0**2))))
           if (y0.ne.0.0d0) then
                 gb_yy(i,j,k)=(1/((-1+rho0**2)**2))
     &                        *((8*y0**2*(-x0**2+y0**2+z0**2)
     &                        -8*(y0**2+2*z0**2)*rho0**2
     &                        -4*L**2*(2*y0**4+z0**2*(-1+rho0**2)**2
     &                        +y0**2*(1-2*x0**2+2*z0**2+rho0**4)))
     &                        /((y0**2+z0**2)*(4*rho0**2
     &                         +L**2*(-1+rho0**2)**2))
     &                        +4*(z0**2/(y0**2+z0**2)
     &                        +(x0**2*y0**2)/((y0**2+z0**2)*rho0**2)
     &                        +(L**2*y0**2*(1+rho0**2)**2
     &                         *(-1+((-1+rho_h)/(-1+rho0))**(2*n))
     &                         *((-1+rho_h)/(-1+rho0))**(-2*n))
     &                        /(4*rho0**4+L**2*rho0*(-1+rho0**2)**2
     &                         *(rho0+(1.0d0/2.0d0)*r0*(-1+rho0**2)))))
          else
                 gb_yy(i,j,k)=0.0d0
          end if
                 gb_yz(i,j,k)=(4*L**2*y0*z0*(1+rho0**2)**2
     &                        *(-4*rho0**3+L**2*(-1+rho0**2)**2
     &                        *(-rho0-(1.0d0/2.0d0)*r0*(-1+rho0**2)
     &                        *((-1+rho_h)/(-1+rho0))**(2*n)))
     &                        *((-1+rho_h)/(-1+rho0))**(-2*n))
     &                        /(rho0**2*(-1+rho0**2)**2*(4*rho0**2
     &                        +L**2*(-1+rho0**2)**2)*(4*rho0**3
     &                        +L**2*(-1+rho0**2)**2
     &                        *(rho0+(1.0d0/2.0d0)*r0*(-1+rho0**2))))
           if (z0.ne.0.0d0) then
                 gb_zz(i,j,k)=(1/((-1+rho0**2)**2))
     &                        *(
     &                         (
     &                         8*z0**2*(-x0**2+y0**2+z0**2)
     &                        -8*(2*y0**2+z0**2)*rho0**2
     &                        -4*L**2*(z0**2*(1-2*x0**2
     &                         +2*z0**2+rho0**4)
     &                         +y0**2*(2*z0**2+(-1+rho0**2)**2))
     &                         )
     &                        /((y0**2+z0**2)*(4*rho0**2
     &                         +L**2*(-1+rho0**2)**2))
     &                        +4*(y0**2/(y0**2+z0**2)
     &                        +(x0**2*z0**2)/((y0**2+z0**2)*rho0**2)
     &                        +(L**2*z0**2*(1+rho0**2)**2
     &                         *(-1+((-1+rho_h)/(-1+rho0))**(2*n))
     &                         *((-1+rho_h)/(-1+rho0))**(-2*n))
     &                        /(4*rho0**4+L**2*rho0*(-1+rho0**2)**2
     &                         *(rho0+(1.0d0/2.0d0)*r0*(-1+rho0**2))))
     &                        )
           else
                 gb_zz(i,j,k)=0.0d0
           end if


!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0=',L,i,j,k,x0,y0,z0,rho0
!       write (*,*) ' gb_tt=',gb_tt(i,j,k)
!       write (*,*) ' gb_tx=',gb_tx(i,j,k)
!       write (*,*) ' gb_ty=',gb_ty(i,j,k)
!       write (*,*) ' gb_tz=',gb_tz(i,j,k)
!       write (*,*) ' gb_xx=',gb_xx(i,j,k)
!       write (*,*) ' gb_xy=',gb_xy(i,j,k)
!       write (*,*) ' gb_xz=',gb_xz(i,j,k)
!       write (*,*) ' gb_yy=',gb_yy(i,j,k)
!       write (*,*) ' gb_yz=',gb_yz(i,j,k)
!       write (*,*) ' gb_zz=',gb_zz(i,j,k)
!
!                 ! (Schw coordinates)!
!                 ! TODO: add AdS_L dependence; currently assumes AdS_L=1!
!                 gb_tt(i,j,k)=0
!                 gb_tx(i,j,k)=0
!                 gb_ty(i,j,k)=0
!                 gb_tz(i,j,k)=0
!                 gb_xx(i,j,k)=0
!                 gb_xy(i,j,k)=0
!                 gb_xz(i,j,k)=0
!                 gb_yy(i,j,k)=0
!                 gb_yz(i,j,k)=0
!                 gb_zz(i,j,k)=0

              end if
            end do
           end do
        end do

        ! y=0 axis regularization
!        call axi_reg_g(gb_tt,gb_tx,gb_ty,
!     &                 gb_xx,gb_xy,gb_yy,psi,tfunction,chr,ex,
!     &                 L,x,y,z,Nx,Ny,Nz,regtype)

!        call axi_reg_Hb(Hb_t,Hb_x,Hb_y,chr,ex,L,x,y,z,Nx,Ny,Nz,regtype)

        return
        end

c----------------------------------------------------------------------
c calculates inverse of a 4x4 matrix with components g11,...,g44
c----------------------------------------------------------------------
        subroutine calc_g0uu(g11,g12,g13,g14,g22,g23,g24,
     &                       g33,g34,g44,g0_uu)
        implicit none
        real*8 g11,g12,g13,g14
        real*8 g22,g23,g24
        real*8 g33,g34
        real*8 g44
        real*8 g0_uu(4,4)                        

        real*8 invdenominator

        !--------------------------------------------------------------

        invdenominator=
     &      g14**2*g23**2 - 2*g13*g14*g23*g24 + g13**2*g24**2 - 
     &      g14**2*g22*g33 + 2*g12*g14*g24*g33 - g11*g24**2*g33 + 
     &      2*g13*g14*g22*g34 - 2*g12*g14*g23*g34 - 2*g12*g13*g24*g34 + 
     &      2*g11*g23*g24*g34 + g12**2*g34**2 - g11*g22*g34**2 - 
     &      g13**2*g22*g44 + 2*g12*g13*g23*g44 - g11*g23**2*g44 - 
     &      g12**2*g33*g44 + g11*g22*g33*g44

       if (invdenominator.ne.0.0d0) then
        g0_uu(1,1)=
     &      (
     &      -(g24**2*g33) + 2*g23*g24*g34 - g22*g34**2 - g23**2*g44 
     &      + g22*g33*g44
     &      )
     &      /invdenominator
        g0_uu(1,2)=
     &      (
     &      g14*g24*g33 - g14*g23*g34 - g13*g24*g34 + g12*g34**2 + 
     &      g13*g23*g44 - g12*g33*g44
     &      )
     &      /invdenominator
        g0_uu(1,3)=
     &      (
     &      -(g14*g23*g24) + g13*g24**2 + g14*g22*g34 - g12*g24*g34 - 
     &      g13*g22*g44 + g12*g23*g44
     &      )
     &      /invdenominator
        g0_uu(1,4)=
     &      (
     &      g14*g23**2 - g13*g23*g24 - g14*g22*g33 + g12*g24*g33 + 
     &      g13*g22*g34 - g12*g23*g34
     &      )
     &      /invdenominator
        g0_uu(2,2)=
     &      (
     &      -(g14**2*g33) + 2*g13*g14*g34 - g11*g34**2 - g13**2*g44 
     &      + g11*g33*g44
     &      )
     &      /invdenominator
        g0_uu(2,3)=
     &      (
     &      g14**2*g23 - g13*g14*g24 - g12*g14*g34 + g11*g24*g34 + 
     &      g12*g13*g44 - g11*g23*g44
     &      )
     &      /invdenominator
        g0_uu(2,4)=
     &      (
     &      -(g13*g14*g23) + g13**2*g24 + g12*g14*g33 - g11*g24*g33 - 
     &      g12*g13*g34 + g11*g23*g34
     &      )
     &      /invdenominator
        g0_uu(3,3)=
     &      (
     &      -(g14**2*g22) + 2*g12*g14*g24 - g11*g24**2 - g12**2*g44 
     &      + g11*g22*g44
     &      )
     &      /invdenominator
        g0_uu(3,4)=
     &      (
     &      g13*g14*g22 - g12*g14*g23 - g12*g13*g24 + g11*g23*g24 + 
     &      g12**2*g34 - g11*g22*g34
     &      )
     &      /invdenominator
        g0_uu(4,4)=
     &      (
     &      -(g13**2*g22) + 2*g12*g13*g23 - g11*g23**2 - g12**2*g33 
     &      + g11*g22*g33
     &      )
     &      /invdenominator
       else
        write(*,*) 'the metric g0 is singular, the determinant is 0 
     &              at some point'
!        g0_uu(1,1)=0
!        g0_uu(1,2)=0
!        g0_uu(1,3)=0
!        g0_uu(1,4)=0
!        g0_uu(2,2)=0
!        g0_uu(2,3)=0
!        g0_uu(2,4)=0
!        g0_uu(3,3)=0
!        g0_uu(3,4)=0
!        g0_uu(4,4)=0
         stop
       end if
        
       return
       end  

c----------------------------------------------------------------------
c calculates all the tensorial objects in x,y coordinates, at point i,j
c----------------------------------------------------------------------
        subroutine tensor_init(
     &                  gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                  gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                  gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                  gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &                  gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                  gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                  gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &                  gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                  gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &                  gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &                  Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                  Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                  Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                  Hb_z_np1,Hb_z_n,Hb_z_nm1,
     &                  phi1_np1,phi1_n,phi1_nm1,
     &                  g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &                  gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &                  h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &                  A_l,A_l_x,Hads_l,
     &                  gamma_ull,gamma_ull_x,
     &                  riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                  einstein_ll,set_ll,
     &                  phi10_x,phi10_xx,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,i,j,k)

        implicit none

        integer Nx,Ny,Nz
        integer i,j,k

        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz),dt,L

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

        integer a,b,c,d,e,f,g,h
        real*8 dx,dy,dz
        real*8 x0,y0,z0
        real*8 rho0
        real*8 f0

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 grad_phi1_sq

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,w,y,z)
        !--------------------------------------------------------------
        real*8 g0_ll(4,4),g0_uu(4,4)
        real*8 g0_ll_x(4,4,4),g0_uu_x(4,4,4),g0_ll_xx(4,4,4,4)
        real*8 gads_ll(4,4),gads_uu(4,4)
        real*8 gads_ll_x(4,4,4),gads_uu_x(4,4,4),gads_ll_xx(4,4,4,4)
        real*8 h0_ll(4,4),h0_uu(4,4)
        real*8 h0_ll_x(4,4,4),h0_uu_x(4,4,4),h0_ll_xx(4,4,4,4)
        real*8 gamma_ull(4,4,4),gamma_ull_x(4,4,4,4)
        real*8 riemann_ulll(4,4,4,4)
        real*8 ricci_ll(4,4),ricci_lu(4,4),ricci
        real*8 einstein_ll(4,4),set_ll(4,4)
        real*8 Hads_l(4),A_l(4),A_l_x(4,4)
        real*8 phi10_x(4),phi10_xx(4,4)

        !--------------------------------------------------------------
        ! the following are first and second time derivatives of *n*
        ! level variables, and as these are the only derivatives we
        ! use we drop any _n identifier
        !--------------------------------------------------------------
        real*8 gb_tt_t, gb_tt_x, gb_tt_y
        real*8 gb_tt_z
        real*8 gb_tt_tt,gb_tt_tx,gb_tt_ty
        real*8 gb_tt_tz
        real*8 gb_tt_xx,gb_tt_xy
        real*8 gb_tt_xz
        real*8 gb_tt_yy
        real*8 gb_tt_yz
        real*8 gb_tt_zz

        real*8 gb_tx_t, gb_tx_x, gb_tx_y
        real*8 gb_tx_z
        real*8 gb_tx_tt,gb_tx_tx,gb_tx_ty
        real*8 gb_tx_tz
        real*8 gb_tx_xx,gb_tx_xy
        real*8 gb_tx_xz
        real*8 gb_tx_yy
        real*8 gb_tx_yz
        real*8 gb_tx_zz

        real*8 gb_ty_t, gb_ty_x, gb_ty_y
        real*8 gb_ty_z
        real*8 gb_ty_tt,gb_ty_tx,gb_ty_ty
        real*8 gb_ty_tz
        real*8 gb_ty_xx,gb_ty_xy
        real*8 gb_ty_xz
        real*8 gb_ty_yy
        real*8 gb_ty_yz
        real*8 gb_ty_zz

        real*8 gb_tz_t, gb_tz_x, gb_tz_y
        real*8 gb_tz_z
        real*8 gb_tz_tt,gb_tz_tx,gb_tz_ty
        real*8 gb_tz_tz
        real*8 gb_tz_xx,gb_tz_xy
        real*8 gb_tz_xz
        real*8 gb_tz_yy
        real*8 gb_tz_yz
        real*8 gb_tz_zz

        real*8 gb_xx_t, gb_xx_x, gb_xx_y
        real*8 gb_xx_z
        real*8 gb_xx_tt,gb_xx_tx,gb_xx_ty
        real*8 gb_xx_tz
        real*8 gb_xx_xx,gb_xx_xy
        real*8 gb_xx_xz
        real*8 gb_xx_yy
        real*8 gb_xx_yz
        real*8 gb_xx_zz

        real*8 gb_xy_t, gb_xy_x, gb_xy_y
        real*8 gb_xy_z
        real*8 gb_xy_tt,gb_xy_tx,gb_xy_ty
        real*8 gb_xy_tz
        real*8 gb_xy_xx,gb_xy_xy
        real*8 gb_xy_xz
        real*8 gb_xy_yy
        real*8 gb_xy_yz
        real*8 gb_xy_zz

        real*8 gb_xz_t, gb_xz_x, gb_xz_y
        real*8 gb_xz_z
        real*8 gb_xz_tt,gb_xz_tx,gb_xz_ty
        real*8 gb_xz_tz
        real*8 gb_xz_xx,gb_xz_xy
        real*8 gb_xz_xz
        real*8 gb_xz_yy
        real*8 gb_xz_yz
        real*8 gb_xz_zz

        real*8 gb_yy_t, gb_yy_x, gb_yy_y
        real*8 gb_yy_z
        real*8 gb_yy_tt,gb_yy_tx,gb_yy_ty
        real*8 gb_yy_tz
        real*8 gb_yy_xx,gb_yy_xy
        real*8 gb_yy_xz
        real*8 gb_yy_yy
        real*8 gb_yy_yz
        real*8 gb_yy_zz

        real*8 gb_yz_t, gb_yz_x, gb_yz_y
        real*8 gb_yz_z
        real*8 gb_yz_tt,gb_yz_tx,gb_yz_ty
        real*8 gb_yz_tz
        real*8 gb_yz_xx,gb_yz_xy
        real*8 gb_yz_xz
        real*8 gb_yz_yy
        real*8 gb_yz_yz
        real*8 gb_yz_zz

        real*8 gb_zz_t, gb_zz_x, gb_zz_y
        real*8 gb_zz_z
        real*8 gb_zz_tt,gb_zz_tx,gb_zz_ty
        real*8 gb_zz_tz
        real*8 gb_zz_xx,gb_zz_xy
        real*8 gb_zz_xz
        real*8 gb_zz_yy
        real*8 gb_zz_yz
        real*8 gb_zz_zz

        real*8 phi1_t, phi1_x, phi1_y
        real*8 phi1_z
        real*8 phi1_tt,phi1_tx,phi1_ty
        real*8 phi1_tz
        real*8 phi1_xx,phi1_xy
        real*8 phi1_xz
        real*8 phi1_yy
        real*8 phi1_yz
        real*8 phi1_zz

        real*8 gb_tt0,gb_tx0,gb_ty0
        real*8 gb_tz0
        real*8 gb_xx0,gb_xy0
        real*8 gb_xz0
        real*8 gb_yy0
        real*8 gb_yz0
        real*8 gb_zz0
        real*8 phi10

        real*8 g0_tt_ads_x,g0_tt_ads_xx,g0_tt_ads_xy
        real*8 g0_tt_ads_y,g0_tt_ads_yy
        real*8 g0_tt_ads_z
        real*8 g0_tt_ads_xz
        real*8 g0_tt_ads_yz
        real*8 g0_tt_ads_zz
        real*8 g0_tx_ads_x,g0_tx_ads_xx,g0_tx_ads_xy
        real*8 g0_tx_ads_y,g0_tx_ads_yy
        real*8 g0_tx_ads_z
        real*8 g0_tx_ads_xz
        real*8 g0_tx_ads_yz
        real*8 g0_tx_ads_zz
        real*8 g0_ty_ads_x,g0_ty_ads_xx,g0_ty_ads_xy
        real*8 g0_ty_ads_y,g0_ty_ads_yy
        real*8 g0_ty_ads_z
        real*8 g0_ty_ads_xz
        real*8 g0_ty_ads_yz
        real*8 g0_ty_ads_zz
        real*8 g0_tz_ads_x,g0_tz_ads_xx,g0_tz_ads_xy
        real*8 g0_tz_ads_y,g0_tz_ads_yy
        real*8 g0_tz_ads_z
        real*8 g0_tz_ads_xz
        real*8 g0_tz_ads_yz
        real*8 g0_tz_ads_zz
        real*8 g0_xx_ads_x,g0_xx_ads_xx,g0_xx_ads_xy
        real*8 g0_xx_ads_y,g0_xx_ads_yy
        real*8 g0_xx_ads_z
        real*8 g0_xx_ads_xz
        real*8 g0_xx_ads_yz
        real*8 g0_xx_ads_zz
        real*8 g0_xy_ads_x,g0_xy_ads_xx,g0_xy_ads_xy
        real*8 g0_xy_ads_y,g0_xy_ads_yy
        real*8 g0_xy_ads_z
        real*8 g0_xy_ads_xz
        real*8 g0_xy_ads_yz
        real*8 g0_xy_ads_zz
        real*8 g0_xz_ads_x,g0_xz_ads_xx,g0_xz_ads_xy
        real*8 g0_xz_ads_y,g0_xz_ads_yy
        real*8 g0_xz_ads_z
        real*8 g0_xz_ads_xz
        real*8 g0_xz_ads_yz
        real*8 g0_xz_ads_zz
        real*8 g0_yy_ads_x,g0_yy_ads_xx,g0_yy_ads_xy
        real*8 g0_yy_ads_y,g0_yy_ads_yy
        real*8 g0_yy_ads_z
        real*8 g0_yy_ads_xz
        real*8 g0_yy_ads_yz
        real*8 g0_yy_ads_zz
        real*8 g0_yz_ads_x,g0_yz_ads_xx,g0_yz_ads_xy
        real*8 g0_yz_ads_y,g0_yz_ads_yy
        real*8 g0_yz_ads_z
        real*8 g0_yz_ads_xz
        real*8 g0_yz_ads_yz
        real*8 g0_yz_ads_zz
        real*8 g0_zz_ads_x,g0_zz_ads_xx,g0_zz_ads_xy
        real*8 g0_zz_ads_y,g0_zz_ads_yy
        real*8 g0_zz_ads_z
        real*8 g0_zz_ads_xz
        real*8 g0_zz_ads_yz
        real*8 g0_zz_ads_zz

        real*8 g0_tt_ads0,g0_xx_ads0
        real*8 g0_tx_ads0,g0_ty_ads0,g0_tz_ads0
        real*8 g0_xy_ads0,g0_yy_ads0,g0_zz_ads0
        real*8 g0_xz_ads0,g0_yz_ads0

        real*8 detg0_ads0
        real*8 g0u_tt_ads0,g0u_xx_ads0
        real*8 g0u_tx_ads0,g0u_ty_ads0,g0u_tz_ads0
        real*8 g0u_xy_ads0,g0u_yy_ads0,g0u_zz_ads0
        real*8 g0u_xz_ads0,g0u_yz_ads0

        real*8 Hb_t_t,Hb_t_x,Hb_t_y,Hb_t_z
        real*8 Hb_x_t,Hb_x_x,Hb_x_y,Hb_x_z
        real*8 Hb_y_t,Hb_y_x,Hb_y_y,Hb_y_z
        real*8 Hb_z_t,Hb_z_x,Hb_z_y,Hb_z_z

        real*8 Hb_t0,Hb_x0,Hb_y0
        real*8 Hb_z0

!!!!!!!!!!DEBUG DERIVATIVE STENCILS!!!!!!!!!!!
        real*8 testf1(Nx,Ny,Nz),testf2(Nx,Ny,Nz),testf3(Nx,Ny,Nz)
        real*8 testf1_t,testf1_x,testf1_y,testf1_z
        real*8 testf2_t,testf2_x,testf2_y,testf2_z
        real*8 testf3_t,testf3_x,testf3_y,testf3_z
        real*8 testf1_tt,testf1_tx,testf1_ty
        real*8 testf1_xx,testf1_xy,testf1_yy
        real*8 testf1_tz,testf1_xz,testf1_yz,testf1_zz
        real*8 testf2_tt,testf2_tx,testf2_ty
        real*8 testf2_xx,testf2_xy,testf2_yy
        real*8 testf2_tz,testf2_xz,testf2_yz,testf2_zz
        real*8 testf3_tt,testf3_tx,testf3_ty
        real*8 testf3_xx,testf3_xy,testf3_yy
        real*8 testf3_tz,testf3_xz,testf3_yz,testf3_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------
        
        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        x0=x(i)
        y0=y(j)
        z0=z(k)
        rho0=sqrt(x0**2+y0**2+z0**2)


!!DEBUG!!!
!
!        g0_tt_ads0 =0
!        g0_tx_ads0 =0
!        g0_ty_ads0 =0
!        g0_tz_ads0 =0
!        g0_xx_ads0 =0
!        g0_xy_ads0 =0
!        g0_xz_ads0 =0
!        g0_yy_ads0 =0
!        g0_yz_ads0 =0
!        g0_zz_ads0=0
!
!        g0u_tt_ads0 =0
!        g0u_tx_ads0 =0
!        g0u_ty_ads0 =0
!        g0u_tz_ads0 =0
!        g0u_xx_ads0 =0
!        g0u_xy_ads0 =0
!        g0u_xz_ads0 =0
!        g0u_yy_ads0 =0
!        g0u_yz_ads0 =0
!        g0u_zz_ads0=0
!
!        g0_tt_ads_x  =0
!        g0_tt_ads_y  =0
!        g0_tt_ads_z  =0
!        g0_tt_ads_xx =0
!        g0_tt_ads_xy =0
!        g0_tt_ads_xz =0
!        g0_tt_ads_yy =0
!        g0_tt_ads_yz =0
!        g0_tt_ads_zz =0
!
!        g0_tx_ads_x  =0
!        g0_tx_ads_y  =0
!        g0_tx_ads_z  =0
!        g0_tx_ads_xx =0
!        g0_tx_ads_xy =0
!        g0_tx_ads_xz =0
!        g0_tx_ads_yy =0
!        g0_tx_ads_yz =0
!        g0_tx_ads_zz =0
!
!        g0_ty_ads_x  =0
!        g0_ty_ads_y  =0
!        g0_ty_ads_z  =0
!        g0_ty_ads_xx =0
!        g0_ty_ads_xy =0
!        g0_ty_ads_xz =0
!        g0_ty_ads_yy =0
!        g0_ty_ads_yz =0
!        g0_ty_ads_zz =0
!
!        g0_tz_ads_x  =0
!        g0_tz_ads_y  =0
!        g0_tz_ads_z  =0
!        g0_tz_ads_xx =0
!        g0_tz_ads_xy =0
!        g0_tz_ads_xz =0
!        g0_tz_ads_yy =0
!        g0_tz_ads_yz =0
!        g0_tz_ads_zz =0
!
!        g0_xx_ads_x  =0
!        g0_xx_ads_y  =0
!        g0_xx_ads_z  =0
!        g0_xx_ads_xx =0
!        g0_xx_ads_xy =0
!        g0_xx_ads_xz =0
!        g0_xx_ads_yy =0
!        g0_xx_ads_yz =0
!        g0_xx_ads_zz =0
!
!        g0_xy_ads_x  =0
!        g0_xy_ads_y  =0
!        g0_xy_ads_z  =0
!        g0_xy_ads_xx =0
!        g0_xy_ads_xy =0
!        g0_xy_ads_xz =0
!        g0_xy_ads_yy =0
!        g0_xy_ads_yz =0
!        g0_xy_ads_zz =0
!
!        g0_xz_ads_x  =0
!        g0_xz_ads_y  =0
!        g0_xz_ads_z  =0
!        g0_xz_ads_xx =0
!        g0_xz_ads_xy =0
!        g0_xz_ads_xz =0
!        g0_xz_ads_yy =0
!        g0_xz_ads_yz =0
!        g0_xz_ads_zz =0
!
!        g0_yy_ads_x  =0
!        g0_yy_ads_y  =0
!        g0_yy_ads_z  =0
!        g0_yy_ads_xx =0
!        g0_yy_ads_xy =0
!        g0_yy_ads_xz =0
!        g0_yy_ads_yy =0
!        g0_yy_ads_yz =0
!        g0_yy_ads_zz =0
!
!        g0_yz_ads_x  =0
!        g0_yz_ads_y  =0
!        g0_yz_ads_z  =0
!        g0_yz_ads_xx =0
!        g0_yz_ads_xy =0
!        g0_yz_ads_xz =0
!        g0_yz_ads_yy =0
!        g0_yz_ads_yz =0
!        g0_yz_ads_zz =0
!
!        g0_zz_ads_x  =0
!        g0_zz_ads_y  =0
!        g0_zz_ads_z  =0
!        g0_zz_ads_xx =0
!        g0_zz_ads_xy =0
!        g0_zz_ads_xz =0
!        g0_zz_ads_yy =0
!        g0_zz_ads_yz =0
!        g0_zz_ads_zz =0
!
!        Hads_l(1)=0
!        Hads_l(2)=0
!        Hads_l(3)=0
!        Hads_l(4)=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!2+1 version!!!!!!
!       
!        f0=(1-rho0**2)**2+4*rho0**2/L**2
!
!        ! set gads values
!        g0_tt_ads0 =-f0
!     &             /(1-rho0**2)**2
!        g0_tx_ads0 =0
!        g0_ty_ads0 =0
!        g0_tz_ads0 =0
!        g0_xx_ads0 =(x0**2*(1+rho0**2)**2/f0+y0**2)
!     &             /(1-rho0**2)**2
!     &             /rho0**2
!     &             *4
!        g0_xy_ads0 =((1+rho0**2)**2/f0-1)
!     &             /(1-rho0**2)**2
!     &             /rho0**2
!     &             *x0*y0
!     &             *4
!        g0_xz_ads0 =((1+rho0**2)**2/f0-1)
!     &             /(1-rho0**2)**2
!     &             /rho0**2
!     &             *x0*z0
!     &             *4
!        g0_yy_ads0 =(y0**2*(1+rho0**2)**2/f0+x0**2)
!     &             /(1-rho0**2)**2
!     &             /rho0**2
!     &             *4
!        g0_yz_ads0 =((1+rho0**2)**2/f0-1)
!     &             /(1-rho0**2)**2
!     &             /rho0**2
!     &             *y0*z0
!     &             *4
!        g0_psi_ads0=(y0**2)
!     &             /(1-rho0**2)**2
!     &             *4
!!!!!!!!!!!!!!!!!!!


        g0_tt_ads0 =-(4*rho0**2+L**2*(-1+rho0**2)**2)
     &               /L**2/(-1+rho0**2)**2
        g0_tx_ads0 =0
        g0_ty_ads0 =0
        g0_tz_ads0 =0
        g0_xx_ads0 =(8*(-1+L**2)*(x0**2-y0**2-z0**2)
     &              +8*rho0**2+4*L**2*(1+rho0**4))
     &              /(-1+rho0**2)**2/(4*rho0**2+L**2*(-1+rho0**2)**2)
        g0_xy_ads0 =(16 *(-1 + L**2) *x0* y0)
     &              /((-1 + rho0**2)**2 
     &               *(4 *rho0**2 +L**2 *(-1 +rho0**2)**2))
        g0_xz_ads0 =(16 *(-1 + L**2) *x0* z0)
     &              /((-1 + rho0**2)**2
     &               *(4 *rho0**2 +L**2 *(-1 +rho0**2)**2))
!        g0_yy_ads0 =(-8*y0**2*(-x0**2+y0**2+z0**2)
!     &              +8*(y0**2+2*z0**2)*rho0**2
!     &              +4 *L**2* (2 *y0**4 + z0**2 *(-1 + rho0**2)**2
!     &              + y0**2 *(1 - 2* x0**2 + 2 *z0**2 + rho0**4)))
!     &              /((y0**2 + z0**2)* (-1 + rho0**2)**2
!     &              * (4* rho0**2 + L**2* (-1 + rho0**2)**2))
        g0_yy_ads0 =(4*(4*(x0**2+z0**2)+L**2*(x0**4+(1+y0**2)**2
     &              +2*(-1+y0**2)*z0**2+z0**4
     &              +2*x0**2*(-1+y0**2+z0**2))))
     &              /(L**2*(-1+rho0**2)**4+4*(-1+rho0**2)**2*(rho0**2))
        g0_yz_ads0 =(16 *(-1 + L**2) *y0* z0)
     &              /((-1 + rho0**2)**2
     &               *(4 *rho0**2 +L**2 *(-1 +rho0**2)**2))
!        g0_zz_ads0=(-8*z0**2*(-x0**2+y0**2+z0**2)
!     &              +8*(2*y0**2+z0**2)*rho0**2
!     &              +4*L**2*(z0**2*(1-2*x0**2+2*z0**2+rho0**4)
!     &              +y0**2*(2*z0**2+(-1+rho0**2)**2)))
!     &              /(y0**2+z0**2)/(-1+rho0**2)**2
!     &              /(4*rho0**2+L**2*(-1+rho0**2)**2)
        g0_zz_ads0=(4*(4*(x0**2+y0**2)+L**2*((-1+x0**2+y0**2)**2
     &              +2*(1+x0**2+y0**2)*z0**2+z0**4)))
     &              /(L**2*(-1+rho0**2)**4
     &              +4*(-1+rho0**2)**2*(rho0**2))

!!!!!diagnostic!!!!
!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0=',L,i,j,k,x0,y0,z0,rho0
!       write (*,*) ' g0_tt_ads0=',g0_tt_ads0
!       write (*,*) ' g0_xx_ads0=',g0_xx_ads0
!       write (*,*) ' g0_xy_ads0=',g0_xy_ads0
!       write (*,*) ' g0_xz_ads0=',g0_xz_ads0
!       write (*,*) ' g0_yy_ads0=',g0_yy_ads0 
!       write (*,*) ' g0_yz_ads0=',g0_yz_ads0
!       write (*,*) ' g0_zz_ads0=',g0_zz_ads0
!!!!!!!!!!!!!!!

!!!!!!2+1 version!!!!
!        g0u_tt_ads0 =1/g0_tt_ads0
!        g0u_tx_ads0 =0
!        g0u_ty_ads0 =0
!        g0u_tz_ads0 =0
!        g0u_xx_ads0 =g0_yy_ads0/(g0_xx_ads0*g0_yy_ads0-g0_xy_ads0**2)
!        g0u_xy_ads0 =-g0_xy_ads0/(g0_xx_ads0*g0_yy_ads0-g0_xy_ads0**2)
!        g0u_xz_ads0 =0
!        g0u_yy_ads0 =g0_xx_ads0/(g0_xx_ads0*g0_yy_ads0-g0_xy_ads0**2)
!        g0u_yz_ads0 =0
!        g0u_psi_ads0=1/g0_psi_ads0
!!!!!!!!!!!!!!

        detg0_ads0=-g0_tt_ads0*(g0_xz_ads0**2*g0_yy_ads0
     &       -2*g0_xy_ads0*g0_xz_ads0*g0_yz_ads0
     &       +g0_xy_ads0**2*g0_zz_ads0
     &       +g0_xx_ads0*(g0_yz_ads0**2-g0_yy_ads0*g0_zz_ads0))

      if ((detg0_ads0.ne.0.0d0).and.(g0_tt_ads0.ne.0.0d0)) then        
        g0u_tt_ads0 =1/g0_tt_ads0
        g0u_tx_ads0 =0
        g0u_ty_ads0 =0
        g0u_tz_ads0 =0
        g0u_xx_ads0 =(g0_yz_ads0**2-g0_yy_ads0*g0_zz_ads0)
     &               /(g0_xz_ads0**2*g0_yy_ads0
     &               -2*g0_xy_ads0*g0_xz_ads0*g0_yz_ads0
     &               +g0_xy_ads0**2*g0_zz_ads0
     &               +g0_xx_ads0*(g0_yz_ads0**2-g0_yy_ads0*g0_zz_ads0))
        g0u_xy_ads0 =(-g0_xz_ads0*g0_yz_ads0+g0_xy_ads0*g0_zz_ads0)
     &               /(g0_xz_ads0**2*g0_yy_ads0
     &               -2*g0_xy_ads0*g0_xz_ads0*g0_yz_ads0
     &               +g0_xy_ads0**2*g0_zz_ads0
     &               +g0_xx_ads0*(g0_yz_ads0**2-g0_yy_ads0*g0_zz_ads0))
        g0u_xz_ads0 =(g0_xz_ads0*g0_yy_ads0-g0_xy_ads0*g0_yz_ads0)
     &               /(g0_xz_ads0**2*g0_yy_ads0
     &               -2*g0_xy_ads0*g0_xz_ads0*g0_yz_ads0
     &               +g0_xy_ads0**2*g0_zz_ads0
     &               +g0_xx_ads0*(g0_yz_ads0**2-g0_yy_ads0*g0_zz_ads0))
        g0u_yy_ads0 =(g0_xz_ads0**2-g0_xx_ads0*g0_zz_ads0)
     &               /(g0_xz_ads0**2*g0_yy_ads0
     &               -2*g0_xy_ads0*g0_xz_ads0*g0_yz_ads0
     &               +g0_xy_ads0**2*g0_zz_ads0
     &               +g0_xx_ads0*(g0_yz_ads0**2-g0_yy_ads0*g0_zz_ads0))
        g0u_yz_ads0 =(-g0_xy_ads0*g0_xz_ads0+g0_xx_ads0*g0_yz_ads0)
     &               /(g0_xz_ads0**2*g0_yy_ads0
     &               -2*g0_xy_ads0*g0_xz_ads0*g0_yz_ads0
     &               +g0_xy_ads0**2*g0_zz_ads0
     &               +g0_xx_ads0*(g0_yz_ads0**2-g0_yy_ads0*g0_zz_ads0))
        g0u_zz_ads0 =(g0_xy_ads0**2-g0_xx_ads0*g0_yy_ads0)
     &               /(g0_xz_ads0**2*g0_yy_ads0
     &               -2*g0_xy_ads0*g0_xz_ads0*g0_yz_ads0
     &               +g0_xy_ads0**2*g0_zz_ads0
     &               +g0_xx_ads0*(g0_yz_ads0**2-g0_yy_ads0*g0_zz_ads0))
      else
       write (*,*) 'L,i,j,k,x0,y0,z0,rho0=',L,i,j,k,x0,y0,z0,rho0
       write (*,*) 'detg0_ads0=',detg0_ads0
       write (*,*) 'ERROR: the metric gads0 is singular at this point'
       stop
      end if

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0=',L,i,j,k,x0,y0,z0,rho0
!       write (*,*) ' g0u_tt_ads0=',g0u_tt_ads0
!       write (*,*) ' g0u_xx_ads0=',g0u_xx_ads0
!       write (*,*) ' g0u_xy_ads0=',g0u_xy_ads0
!       write (*,*) ' g0u_xz_ads0=',g0u_xz_ads0
!       write (*,*) ' g0u_yy_ads0=',g0u_yy_ads0
!       write (*,*) ' g0u_yz_ads0=',g0u_yz_ads0
!       write (*,*) ' g0u_zz_ads0=',g0u_zz_ads0

        ! set gbar values
        gb_tt0=gb_tt_n(i,j,k)
        gb_tx0=gb_tx_n(i,j,k)
        gb_ty0=gb_ty_n(i,j,k)
        gb_tz0=gb_tz_n(i,j,k)
        gb_xx0=gb_xx_n(i,j,k)
        gb_xy0=gb_xy_n(i,j,k)
        gb_xz0=gb_xz_n(i,j,k)
        gb_yy0=gb_yy_n(i,j,k)
        gb_yz0=gb_yz_n(i,j,k)
        gb_zz0=gb_zz_n(i,j,k)

        ! set hbar values
        Hb_t0=Hb_t_n(i,j,k)
        Hb_x0=Hb_x_n(i,j,k)
        Hb_y0=Hb_y_n(i,j,k)
        Hb_z0=Hb_z_n(i,j,k)

        ! set phi1 value
        phi10=phi1_n(i,j,k)

!CHECKED WITH MATHEMATICA UP TO HERE

!!!!!!!2+1 version!!!!!
!        ! ASSUMES L=1
!        g0_tt_ads_x  =(8*x0*(1 + rho0**2))/(L**2*(-1 + rho0**2)**3)
!        g0_tt_ads_y  =(8*y0*(1 + rho0**2))/(L**2*(-1 + rho0**2)**3)
!        g0_tt_ads_z  =0
!        g0_tt_ads_xx =(-8*(1 + 3*x0**4 - y0**4 + 2*x0**2*(4 + y0**2)))/
!     &  (L**2*(-1 + x0**2 + y0**2)**4)
!        g0_tt_ads_xy =(-32*x0*y0*(2 + x0**2 + y0**2))/
!     &  (L**2*(-1 + x0**2 + y0**2)**4)
!        g0_tt_ads_xz =0
!        g0_tt_ads_yy =(8*(-1 + x0**4 - 2*(4 + x0**2)*y0**2 - 3*y0**4))/
!     &  (L**2*(-1 + x0**2 + y0**2)**4)
!        g0_tt_ads_yz =0
!        g0_tt_ads_zz =0
!
!        g0_tx_ads_x  =0
!        g0_tx_ads_y  =0
!        g0_tx_ads_z  =0
!        g0_tx_ads_xx =0
!        g0_tx_ads_xy =0
!        g0_tx_ads_xz =0
!        g0_tx_ads_yy =0
!        g0_tx_ads_yz =0
!        g0_tx_ads_zz =0
!
!        g0_ty_ads_x  =0
!        g0_ty_ads_y  =0
!        g0_ty_ads_z  =0
!        g0_ty_ads_xx =0
!        g0_ty_ads_xy =0
!        g0_ty_ads_xz =0
!        g0_ty_ads_yy =0
!        g0_ty_ads_yz =0
!        g0_ty_ads_zz =0
!
!        g0_tz_ads_x  =0
!        g0_tz_ads_y  =0
!        g0_tz_ads_z  =0
!        g0_tz_ads_xx =0
!        g0_tz_ads_xy =0
!        g0_tz_ads_xz =0
!        g0_tz_ads_yy =0
!        g0_tz_ads_yz =0
!        g0_tz_ads_zz =0
!
!        g0_xx_ads_x  =(-16*x0*(8*y0**2*(-1 + 3*x0**2 + 3*y0**2) + 
!     &                L**4*(-1 + rho0**2)**2*
!     &                (3 + x0**4 - 4*y0**2 + y0**4 + 2*x0**2*(2+y0**2))+
!     &                2*L**2*(-1 + x0**6 + 11*y0**2 +5*y0**4*(-3+y0**2)+
!     &                x0**4*(5 + 7*y0**2) + 
!     &                x0**2*(3 - 10*y0**2 + 11*y0**4))))/
!     &                ((-1 + rho0**2)**3*
!     &                (L**2*(-1 + rho0**2)**2 + 4*(rho0**2))**2)
!        g0_xx_ads_y  =(-16*y0*(8*(-x0**4 + 2*y0**4 + x0**2*(1+y0**2))+
!     &                L**4*(-1 + rho0**2)**2*
!     &                (x0**4 + (-1 + y0**2)**2 + 2*x0**2*(3 + y0**2)) + 
!     &                8*L**2*(y0**2*(-1 + y0**2)**2 + x0**4*(3 + y0**2)+
!     &                x0**2*(-1 + y0**2 + 2*y0**4))))/
!     &                ((-1 + rho0**2)**3*
!     &                (L**2*(-1 + rho0**2)**2 + 4*(rho0**2))**2)
!        g0_xx_ads_z  =0
!        g0_xx_ads_xx =(16*(32*y0**2*(3*(x0**2 - 4*x0**4 + 7*x0**6) + 
!     &                (-1 - 8*x0**2 + 39*x0**4)*y0**2 + 
!     &                (4 + 15*x0**2)*y0**4 - 3*y0**6) + 
!     &                L**6*(-1 + x0**2 + y0**2)**4*
!     &                (5*x0**6 - (-3 + y0**2)*(-1 + y0**2)**2 + 
!     &                x0**4*(33 + 9*y0**2) + 
!     &                3*x0**2*(13 - 14*y0**2 + y0**4)) + 
!     &                2*L**4*(-1 + x0**2 + y0**2)**2*
!     &                (9*x0**8 + x0**6*(78 + 60*y0**2) - 
!     &                (-1 + y0**2)**2*(1 - 16*y0**2 + 7*y0**4) + 
!     &                2*x0**4*(40 - 67*y0**2 + 43*y0**4) + 
!     &                2*x0**2*(-11 + 88*y0**2 - 91*y0**4 + 14*y0**6)) + 
!     &                8*L**2*(3*x0**10 + 4*x0**8*(7 + 16*y0**2) - 
!     &                2*(y0 - y0**3)**2*(1 - 7*y0**2 + 4*y0**4) + 
!     &                2*x0**6*(13 - 63*y0**2 + 83*y0**4) + 
!     &                6*x0**4*(-2 + 31*y0**2 - 51*y0**4 + 24*y0**6) + 
!     &                x0**2*(-1 + y0)*(1 + y0)*
!     &                (-3 + 31*y0**2 - 91*y0**4 + 31*y0**6))))/
!     &                ((-1 + x0**2 + y0**2)**4*
!     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)
!        g0_xx_ads_xy =(32*x0*y0*(8*L**2*(1 + 6*x0**2-28*x0**4+50*x0**6 -
!     &                5*x0**8 + 2*(-7 + 5*x0**2*(3 + 3*x0**2 + x0**4))*
!     &                y0**2 + 2*(29 - 45*x0**2 + 30*x0**4)*y0**4 + 
!     &                70*(-1 + x0**2)*y0**6 + 25*y0**8) + 
!     &                L**6*(-1 + x0**2 + y0**2)**4*
!     &                (11 + 3*x0**4 - 14*y0**2 + 3*y0**4 + 
!     &                x0**2*(26 + 6*y0**2)) - 
!     &                32*(3*x0**6 - y0**2 + 4*y0**4 - 9*y0**6 - 
!     &                x0**4*(4 + 3*y0**2) + x0**2*(1 - 15*y0**4)) + 
!     &                4*L**4*(-1 + x0**2 + y0**2)**2*
!     &                (-4 + x0**6 + 31*y0**2 - 38*y0**4 + 11*y0**6 + 
!     &                x0**4*(42 + 13*y0**2) + 
!     &                x0**2*(-3 + 4*y0**2 + 23*y0**4))))/
!     &                ((-1 + x0**2 + y0**2)**4*
!     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)
!        g0_xx_ads_xz =0
!        g0_xx_ads_yy =(-16*(L**6*(-1 + x0**2 + y0**2)**4*
!     &                (-1 - 5*x0**2 + 5*x0**4 + x0**6 - 
!     &                3*(1 + 22*x0**2 + x0**4)*y0**2 - 
!     &                9*(-1 + x0**2)*y0**4 - 5*y0**6) + 
!     &                4*L**4*(-1 + x0**2 + y0**2)**2*
!     &                (x0**8 - 3*y0**2*(-1 + y0**2)**2*(1 + 5*y0**2) + 
!     &                x0**6*(11 + 8*y0**2) - 
!     &                x0**4*(13 + 111*y0**2 + 2*y0**4) + 
!     &                x0**2*(1 + 46*y0**2 - 95*y0**4 - 24*y0**6)) - 
!     &                32*(x0**8 - x0**6*(2 + 11*y0**2) + 
!     &                x0**4*(1 + 14*y0**2 - 15*y0**4) + 
!     &                x0**2*y0**2*(-3 + 18*y0**2 + 7*y0**4) + 
!     &                2*(y0**6 + 5*y0**8)) - 
!     &                8*L**2*(x0**10 + 
!     &                6*y0**4*(-1 + y0**2)**2*(1 + 5*y0**2) - 
!     &                2*x0**8*(8 + 13*y0**2) + 
!     &                x0**6*(22 + 138*y0**2 - 54*y0**4) + 
!     &                2*x0**4*(-4 - 55*y0**2 + 135*y0**4 + 2*y0**6) + 
!     &                x0**2*(1 + 38*y0**2 - 114*y0**4 + 62*y0**6 + 
!     &                61*y0**8))))/
!     &                ((-1 + x0**2 + y0**2)**4*
!     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)
!        g0_xx_ads_yz =0
!        g0_xx_ads_zz =0
!
!        g0_xy_ads_x  =(-16*(-1 + L)*(1 + L)*y0*
!     &                (L**2*(1 + 7*x0**2 - y0**2)*(-1 + rho0**2)**2 + 
!     &                4*(5*x0**4 + y0**2 - y0**4 +x0**2*(-1+4*y0**2))))/
!     &                ((-1 + rho0**2)**3*
!     &                (L**2*(-1 + rho0**2)**2 + 4*(rho0**2))**2)
!        g0_xy_ads_y  =(16*(-1 + L)*(1 + L)*x0*
!     &                (L**2*(-1 + x0**2 - 7*y0**2)*(-1 + rho0**2)**2 + 
!     &                4*(x0**4 + y0**2 - 5*y0**4 - x0**2*(1+4*y0**2))))/
!     &                ((-1 + rho0**2)**3*
!     &                (L**2*(-1 + rho0**2)**2 + 4*(rho0**2))**2)
!        g0_xy_ads_z  =0
!        g0_xy_ads_xx =(-128*x0*y0*(4*(x0**2 - 3*y0**2) - 
!     &                L**6*(3+7*x0**2-3*y0**2)*(-1+ x0**2 + y0**2)**4 + 
!     &                4*(x0**2 + y0**2)*
!     &                (x0**2*(-4 + 15*x0**2) + 6*(2 + x0**2)*y0**2 - 
!     &                9*y0**4) + L**4*(-1 + x0**2 + y0**2)**2*
!     &                (6 + x0**2 - 50*x0**4 + 7*x0**6 + 
!     &                (-33 - 20*x0**2 + 11*x0**4)*y0**2 + 
!     &                (30 + x0**2)*y0**4 - 3*y0**6) + 
!     &                L**2*(-3+2*x0**2+52*x0**4- 138*x0**6 + 39*x0**8 + 
!     &                2*(21 - 34*x0**2 - 87*x0**4 + 48*x0**6)*y0**2 + 
!     &                6*(-20 + 11*x0**2 + 9*x0**4)*y0**4 + 
!     &                6*(17 - 4*x0**2)*y0**6 - 21*y0**8)))/
!     &                ((-1 + x0**2 + y0**2)**4*
!     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)
!        g0_xy_ads_xy =(-16*(-(L**4*(-1 + x0**2 + y0**2)**2*
!     &                (-1 - 4*x0**2 + 66*x0**4 - 68*x0**6 + 7*x0**8 - 
!     &                4*(1 + 35*x0**2 - 109*x0**4 + 13*x0**6)*y0**2 + 
!     &                2*(33 + 218*x0**2 - 59*x0**4)*y0**4 - 
!     &                4*(17 + 13*x0**2)*y0**6 + 7*y0**8)) + 
!     &                L**6*(-1 + x0**2 + y0**2)**4*
!     &                (-1 + 7*x0**4 - 6*y0**2 + 7*y0**4 - 
!     &                6*x0**2*(1 + 11*y0**2)) + 
!     &                16*(-5*x0**8 - y0**4 + 6*y0**6 - 5*y0**8 + 
!     &                x0**6*(6 + 28*y0**2) + 
!     &                2*x0**2*y0**2*(3 - 7*y0**2 + 14*y0**4) + 
!     &                x0**4*(-1 - 14*y0**2 + 66*y0**4)) + 
!     &                16*L**2*(-3*x0**10 + x0**8*(14 + 15*y0**2) + 
!     &                x0**6*(-15 - 64*y0**2 + 60*y0**4) + 
!     &                y0**4*(4 - 15*y0**2 + 14*y0**4 - 3*y0**6) + 
!     &                x0**2*y0**2*
!     &                (-12 + 41*y0**2 - 64*y0**4 + 15*y0**6) + 
!     &                x0**4*(4 + 41*y0**2 - 156*y0**4 + 60*y0**6))))/
!     &                ((-1 + x0**2 + y0**2)**4*
!     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)
!        g0_xy_ads_xz =0
!        g0_xy_ads_yy =(-128*x0*y0*(4*(-3*x0**2 + y0**2) + 
!     &                L**6*(-3 + 3*x0**2 - 7*y0**2)*
!     &                (-1 + x0**2 + y0**2)**4 - 
!     &                4*(x0**2 + y0**2)*
!     &                (9*x0**4+4*y0**2 - 15*y0**4 - 6*x0**2*(2 + y0**2))
!     &                - L**4*(-1 + x0**2 + y0**2)**2*
!     &                (-6 + 3*x0**6 - y0**2 + 50*y0**4 - 7*y0**6 - 
!     &                x0**4*(30 + y0**2) + 
!     &                x0**2*(33 + 20*y0**2 - 11*y0**4)) + 
!     &                L**2*(-3-21*x0**8+2*y0**2+ 52*y0**4 - 138*y0**6 + 
!     &                39*y0**8 + 6*x0**6*(17 - 4*y0**2) + 
!     &                6*x0**4*(-20 + 11*y0**2 + 9*y0**4) + 
!     &                2*x0**2*(21 - 34*y0**2 - 87*y0**4 + 48*y0**6))))/
!     &                ((-1 + x0**2 + y0**2)**4*
!     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)
!        g0_xy_ads_yz =0
!        g0_xy_ads_zz =0
!
!        g0_xz_ads_x  =0
!        g0_xz_ads_y  =0
!        g0_xz_ads_z  =0
!        g0_xz_ads_xx =0
!        g0_xz_ads_xy =0
!        g0_xz_ads_xz =0
!        g0_xz_ads_yy =0
!        g0_xz_ads_yz =0
!        g0_xz_ads_zz =0
!
!        g0_yy_ads_x  =(-16*x0*(16*x0**4 + 8*(1 + x0**2)*y0**2 -8*y0**4+ 
!     &                L**4*(-1 + rho0**2)**2*
!     &                ((-1 + x0**2)**2 + 2*(3 + x0**2)*y0**2 + y0**4) + 
!     &                8*L**2*((x0 - x0**3)**2 + 
!     &                (-1 + x0**2 + 2*x0**4)*y0**2 + (3+x0**2)*y0**4)))/
!     &                ((-1 + rho0**2)**3*
!     &                (L**2*(-1 + rho0**2)**2 + 4*(rho0**2))**2)
!        g0_yy_ads_y  =(-16*y0*(8*x0**2*(-1 + 3*x0**2 + 3*y0**2) + 
!     &                L**4*(-1 + rho0**2)**2*
!     &                (3 + x0**4 + 4*y0**2 + y0**4 +2*x0**2*(-2+y0**2))
!     &                + 2*L**2*(-1 + 5*x0**6 + 3*y0**2 + 5*y0**4 +y0**6+
!     &                x0**4*(-15 + 11*y0**2) + 
!     &                x0**2*(11 - 10*y0**2 + 7*y0**4))))/
!     &                ((-1 + rho0**2)**3*
!     &                (L**2*(-1 + rho0**2)**2 + 4*(rho0**2))**2)
!        g0_yy_ads_z  =0
!        g0_yy_ads_xx =(16*(4*L**4*(-1 + x0**2 + y0**2)**2*
!     &                (3*x0**2*(-1 + x0**2)**2*(1 + 5*x0**2) + 
!     &                (-1 - 46*x0**2 + 95*x0**4 + 24*x0**6)*y0**2 + 
!     &                (13 + 111*x0**2 + 2*x0**4)*y0**4 - 
!     &                (11 + 8*x0**2)*y0**6 - y0**8) + 
!     &                8*L**2*(6*x0**4*(-1 + x0**2)**2*(1 + 5*x0**2) + 
!     &                (1 + 38*x0**2 - 114*x0**4 + 62*x0**6 + 61*x0**8)*
!     &                y0**2 + 2*(-4 - 55*x0**2 + 135*x0**4 + 2*x0**6)*
!     &                y0**4 + 2*(11 + 69*x0**2 - 27*x0**4)*y0**6 - 
!     &                2*(8 + 13*x0**2)*y0**8 + y0**10) + 
!     &                32*(10*x0**8 + y0**4 - 2*y0**6 + y0**8 + 
!     &                3*x0**4*y0**2*(6-5*y0**2) + x0**6*(2 + 7*y0**2) + 
!     &                x0**2*y0**2*(-3 + 14*y0**2 - 11*y0**4)) + 
!     &                L**6*(-1 + x0**2 + y0**2)**4*
!     &                (1 + 5*x0**6 + 5*y0**2 + 9*x0**4*(-1 + y0**2) - 
!     &                y0**4*(5+y0**2) + 3*x0**2*(1 + 22*y0**2 + y0**4)))
!     &                )/((-1 + x0**2 + y0**2)**4*
!     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)
!        g0_yy_ads_xy =(32*x0*y0*(32*(x0**2 - 4*x0**4 + 9*x0**6 + 
!     &                (-1 + 15*x0**4)*y0**2 + (4 + 3*x0**2)*y0**4 - 
!     &                3*y0**6) + L**6*(-1 + x0**2 + y0**2)**4*
!     &                (11 + 3*x0**4 + 26*y0**2 + 3*y0**4 + 
!     &                2*x0**2*(-7 + 3*y0**2)) + 
!     &                4*L**4*(-1 + x0**2 + y0**2)**2*
!     &                (-4 + 11*x0**6 - 3*y0**2 + 42*y0**4 + y0**6 + 
!     &                x0**4*(-38 + 23*y0**2) + 
!     &                x0**2*(31 + 4*y0**2 + 13*y0**4)) + 
!     &                8*L**2*(1+25*x0**8+6*y0**2- 28*y0**4 + 50*y0**6 - 
!     &                5*y0**8 + 70*x0**6*(-1 + y0**2) + 
!     &                x0**4*(58 - 90*y0**2 + 60*y0**4) + 
!     &                2*x0**2*(-7 + 5*y0**2*(3 + 3*y0**2 + y0**4)))))/
!     &                ((-1 + x0**2 + y0**2)**4*
!     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)
!        g0_yy_ads_xz =0
!        g0_yy_ads_yy =(-16*(32*x0**2*(x0**2 - 4*x0**4 + 3*x0**6 + 
!     &                (-3 + 8*x0**2 - 15*x0**4)*y0**2 + 
!     &                3*(4 - 13*x0**2)*y0**4 - 21*y0**6) + 
!     &                L**6*(-1 + x0**2 + y0**2)**4*
!     &                ((-3 + x0**2)*(-1 + x0**2)**2 - 
!     &                3*(13 - 14*x0**2 + x0**4)*y0**2 - 
!     &                3*(11 + 3*x0**2)*y0**4 - 5*y0**6) + 
!     &                2*L**4*(-1 + x0**2 + y0**2)**2*
!     &                ((-1 + x0**2)**2*(1 - 16*x0**2 + 7*x0**4) + 
!     &                2*(11 - 88*x0**2 + 91*x0**4 - 14*x0**6)*y0**2 - 
!     &                2*(40 - 67*x0**2 + 43*x0**4)*y0**4 - 
!     &                6*(13 + 10*x0**2)*y0**6 - 9*y0**8) + 
!     &                8*L**2*(2*(x0-x0**3)**2*(1 - 7*x0**2 + 4*x0**4) + 
!     &                (-3+ 34*x0**2 - 122*x0**4 + 122*x0**6 - 31*x0**8)*
!     &                y0**2 - 6*(-2 + 31*x0**2 - 51*x0**4 + 24*x0**6)*
!     &                y0**4 - 2*(13 - 63*x0**2 + 83*x0**4)*y0**6 - 
!     &                4*(7 + 16*x0**2)*y0**8 - 3*y0**10)))/
!     &                ((-1 + x0**2 + y0**2)**4*
!     &                (L**2*(-1+x0**2+y0**2)**2 + 4*(x0**2 + y0**2))**3)
!        g0_yy_ads_yz =0
!        g0_yy_ads_zz =0
!
!        g0_yz_ads_x  =0
!        g0_yz_ads_y  =0
!        g0_yz_ads_z  =0
!        g0_yz_ads_xx =0
!        g0_yz_ads_xy =0
!        g0_yz_ads_xz =0
!        g0_yz_ads_yy =0
!        g0_yz_ads_yz =0
!        g0_yz_ads_zz =0
!
!        g0_psi_ads_x  =(-16*x0*y0**2)/(-1 + rho0**2)**3
!        g0_psi_ads_y  =(-8*(y0 - x0**2*y0 + y0**3))/(-1 + rho0**2)**3
!        g0_psi_ads_z  =0
!        g0_psi_ads_xx =(-16*y0**2*(-1 - 5*x0**2 + y0**2))/
!     &                 (-1 + x0**2 + y0**2)**4
!        g0_psi_ads_xy =(-32*x0*y0*(-1 + x0**2 - 2*y0**2))/
!     &                 (-1 + x0**2 + y0**2)**4
!        g0_psi_ads_xz =0
!        g0_psi_ads_yy =(8*((-1+x0**2)**2-8*(-1+x0**2)*y0**2+3*y0**4))/
!     &                 (-1 + x0**2 + y0**2)**4
!        g0_psi_ads_yz =0
!        g0_psi_ads_zz =0
!!!!!!!!!!!!!!!!!!!

        g0_tt_ads_x  =(8*x0*(1 + rho0**2))/(L**2*(-1 + rho0**2)**3)
        g0_tt_ads_y  =(8*y0*(1 + rho0**2))/(L**2*(-1 + rho0**2)**3)
        g0_tt_ads_z  =(8*z0*(1 + rho0**2))/(L**2*(-1 + rho0**2)**3)
        g0_tt_ads_xx =-((8*(1+3*x0**4-(y0**2+z0**2)**2
     &                +2*x0**2*(4+y0**2+z0**2)))
     &                /(L**2 *(-1 + rho0**2)**4))
        g0_tt_ads_xy =(-32*x0*y0*(2 + rho0**2))/
     &                (L**2*(-1 + rho0**2)**4)
        g0_tt_ads_xz =(-32*x0*z0*(2 + rho0**2))/
     &                (L**2*(-1 + rho0**2)**4)
        g0_tt_ads_yy =(8*(-1-8*y0**2+(x0**2-3*y0**2+z0**2)*rho0**2))
     &                /(L**2*(-1+rho0**2)**4)
        g0_tt_ads_yz =(-32*y0*z0*(2 + rho0**2))/
     &                (L**2*(-1 + rho0**2)**4)
        g0_tt_ads_zz =(8*(-1-8*z0**2+(x0**2-3*z0**2+y0**2)*rho0**2))
     &                /(L**2*(-1+rho0**2)**4)

        g0_tx_ads_x  =0
        g0_tx_ads_y  =0
        g0_tx_ads_z  =0
        g0_tx_ads_xx =0
        g0_tx_ads_xy =0
        g0_tx_ads_xz =0
        g0_tx_ads_yy =0
        g0_tx_ads_yz =0
        g0_tx_ads_zz =0

        g0_ty_ads_x  =0
        g0_ty_ads_y  =0
        g0_ty_ads_z  =0
        g0_ty_ads_xx =0
        g0_ty_ads_xy =0
        g0_ty_ads_xz =0
        g0_ty_ads_yy =0
        g0_ty_ads_yz =0
        g0_ty_ads_zz =0

        g0_tz_ads_x  =0
        g0_tz_ads_y  =0
        g0_tz_ads_z  =0
        g0_tz_ads_xx =0
        g0_tz_ads_xy =0
        g0_tz_ads_xz =0
        g0_tz_ads_yy =0
        g0_tz_ads_yz =0
        g0_tz_ads_zz =0

        g0_xx_ads_x  =(4*x0* (-32* (y0**2 + z0**2)* (-1 + 3*rho0**2)
     &                 -4* L**4* (-1 +rho0**2)**2
     &                 *(x0**4 + (-3 + y0**2+z0**2)*(-1+y0**2+z0**2)
     &                 +2*x0**2* (2 + y0**2 + z0**2)) 
     &                 +8* L**2* (1 - x0**6 - 11* y0**2 - 11* z0**2
     &                 -5* (-3 + y0**2 + z0**2)* (y0**2 + z0**2)**2 
     &                 - x0**4* (5 + 7* y0**2 + 7 *z0**2) 
     &                 - x0**2* (3 + 11*y0**4 - 10*z0**2+11* z0**4 
     &                 +2*y0**2* (-5 + 11* z0**2)))))
     &                 /((-1 + rho0**2)**3* (L**2* (-1 +rho0**2)**2 
     &                 + 4* (rho0**2))**2)
        g0_xx_ads_y  =(4*y0*(32* (x0**4 - 2* (y0**2 + z0**2)**2 
     &                - x0**2* (1 + y0**2 + z0**2)) 
     &                - 4* L**4* (-1 + rho0**2)**2
     &                * (x0**4 + (-1 + y0**2 + z0**2)**2 
     &                +2* x0**2* (3 + y0**2 + z0**2)) 
     &                - 32* L**2* ((-1+y0**2+z0**2)**2*(y0**2 + z0**2) 
     &                + x0**4* (3 + y0**2 + z0**2) 
     &                +x0**2*(1+y0**2+z0**2)*(-1+2*y0**2+2*z0**2))))
     &                /((-1 +rho0**2)**3
     &                *(L**2* (-1 +rho0**2)**2+4*(rho0**2))**2)
        g0_xx_ads_z  =(4*z0*(32* (x0**4 - 2* (y0**2 + z0**2)**2
     &                - x0**2* (1 + y0**2 + z0**2))
     &                - 4* L**4* (-1 + rho0**2)**2
     &                * (x0**4 + (-1 + y0**2 + z0**2)**2
     &                +2* x0**2* (3 + y0**2 + z0**2))
     &                - 32* L**2* ((-1+y0**2+z0**2)**2*(y0**2 + z0**2)
     &                + x0**4* (3 + y0**2 + z0**2)
     &                +x0**2*(1+y0**2+z0**2)*(-1+2*y0**2+2*z0**2))))
     &                /((-1 +rho0**2)**3
     &                *(L**2* (-1 +rho0**2)**2+4*(rho0**2))**2)
        g0_xx_ads_xx =8 *((12* x0**2* (1 + (-1 + L**2)* x0**2))
     &                    /(-1 +rho0**2)**4 
     &                 + (-2-2*(-1 + L**2)* x0**2*(5+2*x0**2))
     &                   /(-1 + rho0**2)**3 
     &                 + ((-1+L**2)*(1+5* x0**2))
     &                   /(-1 +rho0**2)**2 
     &                 + (1 - L**2)
     &                   /(-1 + rho0**2) 
     &                 - (64*(-1+L**2)**2*x0**4*(4+L**2*(-2+rho0**2)))
     &                   /(L**2*(-1+rho0**2)**2+4*(rho0**2))**3 
     &                 -((-1+L**2)*(-4+L**2*(2+4*x0**2-y0**2-z0**2)))
     &                   /(L**2*(-1 +rho0**2)**2+4*(rho0**2))
     &                 +(2*(-1+L**2)*x0**2
     &                  *(-40+2*L**2*(3*x0**2-5*(-4+y0**2+z0**2))
     &                  +L**4*(2*x0**4+5*(-1+y0**2+z0**2)
     &                  +x0**2*(-3+2*y0**2+2*z0**2))))
     &                  /(L**2*(-1+rho0**2)**2 + 4* (rho0**2))**2)
        g0_xx_ads_xy =(32*x0*y0*(L**6*(-1 +rho0**2)**4
     &                 *(3*x0**4+(-1+y0**2+z0**2)*(-11+3*y0**2+3*z0**2)
     &                 +x0**2*(26+6*y0**2+6*z0**2))
     &                 +32*(-3*x0**6+y0**2+z0**2
     &                 +x0**4*(4+3*y0**2+3*z0**2) 
     &                 +(y0**2+z0**2)**2*(-4+9*y0**2+9*z0**2)
     &                 +x0**2*(-1+15*(y0**2+z0**2)**2))
     &                 +4*L**4* (-1 +rho0**2)**2
     &                 * (-4 + x0**6 + 31* y0**2 +31*z0**2 
     &                 + (y0**2 + z0**2)**2*(-38 + 11* y0**2+11*z0**2)
     &                 +x0**4*(42 + 13*y0**2 + 13*z0**2)
     &                 +x0**2*(-3+23*y0**4+4*z0**2+23*z0**4
     &                  +y0**2*(4+46*z0**2)))
     &                 +8*L**2* (-5* x0**8 + 10* x0**6*(5+y0**2+z0**2)
     &                 +2*x0**4*(-14+30*y0**4+15*z0**2+30*z0**4
     &                  +15*y0**2*(1+4*z0**2)) 
     &                 +(-1+y0**2+z0**2)*(-1+5*y0**2+5*z0**2)
     &                 *(1+5*y0**4-8*z0**2+5*z0**4+2*y0**2*(-4+5*z0**2))
     &                 +2*x0**2*(3+35*y0**6+15*y0**4*(-3+7*z0**2)
     &                 +5*z0**2*(3-9*z0**2+7*z0**4) 
     &                 +15*y0**2*(1-6*z0**2+7*z0**4)))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_xx_ads_xz =(32*x0*z0*(L**6*(-1 +rho0**2)**4
     &                 *(3*x0**4+(-1+y0**2+z0**2)*(-11+3*y0**2+3*z0**2) 
     &                 +x0**2*(26+6*y0**2+6*z0**2))
     &                 +32*(-3*x0**6+y0**2+z0**2
     &                 +x0**4*(4+3*y0**2+3*z0**2) 
     &                 +(y0**2+z0**2)**2*(-4+9*y0**2+9*z0**2)
     &                 +x0**2*(-1+15*(y0**2+z0**2)**2))
     &                 +4*L**4* (-1 +rho0**2)**2
     &                 * (-4 + x0**6 + 31* y0**2 +31*z0**2 
     &                 + (y0**2 + z0**2)**2*(-38 + 11* y0**2+11*z0**2)
     &                 +x0**4*(42 + 13*y0**2 + 13*z0**2)
     &                 +x0**2*(-3+23*y0**4+4*z0**2+23*z0**4
     &                  +y0**2*(4+46*z0**2)))
     &                 +8*L**2* (-5* x0**8 + 10* x0**6*(5+y0**2+z0**2) 
     &                 +2*x0**4*(-14+30*y0**4+15*z0**2+30*z0**4
     &                  +15*y0**2*(1+4*z0**2)) 
     &                 +(-1+y0**2+z0**2)*(-1+5*y0**2+5*z0**2)
     &                 *(1+5*y0**4-8*z0**2+5*z0**4+2*y0**2*(-4+5*z0**2))
     &                 +2*x0**2*(3+35*y0**6+15*y0**4*(-3+7*z0**2)
     &                 +5*z0**2*(3-9*z0**2+7*z0**4) 
     &                 +15*y0**2*(1-6*z0**2+7*z0**4)))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_xx_ads_yy =(16*(-4*L**4*(-1+rho0**2)**2
     &                 *(x0**8-3*(1+5*y0**2-z0**2)
     &                 *(-1+y0**2+z0**2)**2*(y0**2+z0**2)
     &                 +x0**6*(11+8*y0**2+6*z0**2)
     &                 +x0**4*(-13-111*y0**2-2*y0**4
     &                  +(13+10*y0**2)*z0**2+12*z0**4)
     &                 -x0**2*(-1-46*y0**2+95*y0**4+24*y0**6
     &                  +2*(2+51*y0**2+19*y0**4)*z0**2
     &                  +(7+4*y0**2)*z0**4-10*z0**6))
     &                 -L**6*(-1+rho0**2)**4
     &                 *(x0**6-(1+5*y0**2-z0**2)*(-1+y0**2+z0**2)**2
     &                 +x0**4*(5-3*y0**2+3*z0**2)
     &                 -x0**2*(5+9*y0**4-2*z0**2
     &                  -3*z0**4+6*y0**2*(11+z0**2)))
     &                 +32*(x0**8+x0**6*(-2-11*y0**2+z0**2)
     &                 +2*(1+5*y0**2-z0**2)*(y0**2+z0**2)**3
     &                 -x0**4*(-1+15*y0**4+2*z0**2+3*z0**4
     &                  +2*y0**2*(-7+9*z0**2))
     &                 +x0**2*(7*y0**6+z0**2+2*z0**4-5*z0**6
     &                  +9*y0**4*(2+z0**2)+y0**2*(-3+20*z0**2-3*z0**4)))
     &                 +8*L**2*(x0**10+6*(1+5*y0**2-z0**2)
     &                  *(-1+y0**2+z0**2)**2*(y0**2+z0**2)**2
     &                 -2*x0**8*(8+13*y0**2+z0**2)
     &                 -2*x0**6*(-11+27*y0**4+15*z0**2+9*z0**4
     &                  +y0**2*(-69+36*z0**2))
     &                 +2*x0**4*(-4+2*y0**6+13*z0**2+3*z0**4-16*z0**6
     &                  +3*y0**4*(45-4*z0**2)
     &                  +y0**2*(-55+138*z0**2-30*z0**4))
     &                 +x0**2*(1+61*y0**8-2*z0**2-14*z0**4
     &                  +38*z0**6-23*z0**8+2*y0**6*(31+80*z0**2)
     &                  +6*y0**4*(-19+27*z0**2+19*z0**4)
     &                  +2*y0**2*(19-64*z0**2+69*z0**4-4*z0**6)))))
     &                  /((-1 + rho0**2)**4
     &                  *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_xx_ads_yz =(1/((-1 +rho0**2)**4))*32* y0* z0
     &                *(3+(8*(-1 + L)*(1 +L)*x0**2
     &                * (5* L**4*(-1 + rho0**2)**4+48*(rho0**2)**2 
     &                - 8* (-1 + 4*rho0**2) + 6*L**2*(-1 + rho0**2)**2
     &                * (-2 + 5*rho0**2)))
     &                /(L**2* (-1 + rho0**2)**2 + 4* (rho0**2))**3)
        g0_xx_ads_zz =(16*(-L**6*(-1 +rho0**2)**4
     &                *(x0**6+x0**4*(5+3*y0**2-3*z0**2)
     &                +(-1+y0**2-5*z0**2)*(-1+y0**2+z0**2)**2 
     &                +x0**2*(-5+2*y0**2+3*y0**4
     &                 -6*(11+y0**2)*z0**2-9*z0**4))
     &                +32*(x0**8+x0**6*(-2+y0**2-11*z0**2)
     &                 -2*(-1+y0**2-5*z0**2)*(y0**2+z0**2)**3
     &                 -x0**4*(-1+2*y0**2+3*y0**4
     &                  +2*(-7+9*y0**2)*z0**2+15*z0**4)
     &                 +x0**2*(y0**2+2*y0**4-5*y0**6
     &                 +(-3+20*y0**2-3*y0**4)*z0**2
     &                 +9*(2+y0**2)* z0**4 + 7* z0**6))
     &                -4*L**4*(-1 +rho0**2)**2
     &                *(x0**8+3*(-1+y0**2-5*z0**2)*(-1+y0**2+z0**2)**2
     &                 *(y0**2+z0**2)
     &                +x0**6*(11+6*y0**2+8*z0**2)
     &                +x0**4*(-13+12*y0**4-111*z0**2-2*z0**4
     &                 +y0**2*(13+10*z0**2))
     &                +x0**2*(1+10*y0**6+46*z0**2-95*z0**4-24*z0**6
     &                 -y0**4*(7+4*z0**2)
     &                 -2*y0**2*(2+51*z0**2+19*z0**4)))
     &                +8*L**2*(x0**10
     &                -6*(-1+y0**2-5*z0**2)*(-1+y0**2+z0**2)**2
     &                 *(y0**2+z0**2)**2
     &                -2*x0**8*(8+y0**2+13*z0**2)
     &                +x0**4*(-8+26*y0**2+6*y0**4-32*y0**6
     &                 -2*(55-138*y0**2+30*y0**4)*z0**2 
     &                 +6*(45-4*y0**2)*z0**4+4*z0**6)
     &                -2*x0**6*(-11+9*y0**4-69*z0**2+27*z0**4
     &                 +3*y0**2*(5+12*z0**2))
     &                +x0**2*(1-23*y0**8+38*z0**2-114*z0**4+62*z0**6
     &                 +61*z0**8+y0**6*(38-8*z0**2)
     &                 +2*y0**4*(-7+69*z0**2+57*z0**4)
     &                 +2*y0**2*(-1-64*z0**2+81*z0**4+80*z0**6)))))
     &                /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)

        g0_xy_ads_x  =-((16*(-1+L**2)*y0
     &                *(L**2*(1+7*x0**2-y0**2-z0**2)*(-1+rho0**2)**2
     &                +4*(5*x0**4+y0**2+z0**2-(y0**2+z0**2)**2
     &                +x0**2*(-1+4*y0**2+4*z0**2))))
     &                /((-1+rho0**2)**3
     &                *(L**2*(-1+rho0**2)**2+4*(rho0**2))**2))
        g0_xy_ads_y  =(16*(-1+L**2)*x0*(-4*(x0**2-y0**2+z0**2)
     &                +L**2*(-1+x0**2-7*y0**2+z0**2)*(-1+rho0**2)**2
     &                +4*(x0**2-5*y0**2+z0**2)*(rho0**2)))
     &                /((-1+rho0**2)**3
     &                *(L**2*(-1+rho0**2)**2+4*(rho0**2))**2)
        g0_xy_ads_z  =-((128*(-1+L**2)*x0*y0*z0
     &                *(-1+3*rho0**2+L**2*(-1+rho0**2)**2))
     &                /((-1+rho0**2)**3
     &                *(L**2*(-1+rho0**2)**2+4*(rho0**2))**2))
        g0_xy_ads_xx  =(128*(-1+L**2)*x0*y0*(L**4*(-1+rho0**2)**4
     &                 *(7*x0**2-3*(-1+y0**2+z0**2))
     &                 +4*(x0**2-3*(y0**2+z0**2))
     &                 +4*(rho0**2)*(15*x0**4+12*(y0**2+z0**2)
     &                 -9*(y0**2+z0**2)**2
     &                 +x0**2*(-4+6*y0**2+6*z0**2))
     &                 +3*L**2*(-1+rho0**2)**2
     &                  *(-1+8*y0**2+8*z0**2+(rho0**2)
     &                   *(13*x0**2-7*(y0**2+z0**2)))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_xy_ads_xy  =(16*(-1+L**2)*(-L**4*(-1+rho0**2)**4
     &                 *(7*x0**4+6*x0**2*(-1-11*y0**2+z0**2)
     &                  +(1+7*y0**2-z0**2)*(-1+y0**2+z0**2))
     &                 -8*L**2*(-1+rho0**2)**2
     &                 *(6*x0**6+x0**4*(-6-42*y0**2+11*z0**2)
     &                 +(-1+y0**2+z0**2)
     &                  *(6*y0**4+z0**2+5*y0**2*z0**2-z0**4)
     &                 -2*x0**2*(21*y0**4+2*z0**2-2*z0**4
     &                  +y0**2*(-6+19*z0**2)))
     &                 +16*(-5*x0**8+2*x0**6*(3+14*y0**2-7*z0**2)
     &                  -(-1+y0**2+z0**2)*(y0**2+z0**2)
     &                  *(5*y0**4+z0**2-z0**4+y0**2*(-1+4*z0**2))
     &                 +x0**4*(-1+66*y0**4+10*z0**2-12*z0**4
     &                  +2*y0**2*(-7+27*z0**2))
     &                 +2*x0**2*(14*y0**6+z0**4-z0**6
     &                  +y0**4*(-7+27*z0**2)
     &                  +3*y0**2*(1-2*z0**2+4*z0**4)))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_xy_ads_xz  =(128*(-1+L**2)*y0*z0
     &                 *(12*x0**2-4*(y0**2+z0**2)
     &                 +L**4*(1+9*x0**2-y0**2-z0**2)*(-1+rho0**2)**4
     &                 +4*(rho0**2)*(21*x0**4+4*(y0**2+z0**2)
     &                 -3*(y0**2+z0**2)**2
     &                 +6*x0**2*(-2+3*y0**2+3*z0**2))
     &                 +L**2*(-1+rho0**2)**2
     &                  *(-1+53*x0**4+8*y0**2+8*z0**2-7*(y0**2+z0**2)**2
     &                   +2*x0**2*(-8+23*y0**2+23*z0**2))))
     &                 /((-1+rho0**2)**4 
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_xy_ads_yy  =(128*(-1+L**2)*x0*y0
     &                 *(4*(-3*x0**2+y0**2-3*z0**2)-L**4*(-1+rho0**2)**4
     &                  *(-3+3*x0**2-7*y0**2+3*z0**2)
     &                 -3*L**2*(-1+rho0**2)**2
     &                  *(1+7*x0**4-13*y0**4-2*(4+3*y0**2)*z0**2
     &                  +7*z0**4-2*x0**2*(4+3*y0**2-7*z0**2))
     &                 -4*(rho0**2)*(9*x0**4+4*y0**2-15*y0**4
     &                  -6*(2+y0**2)*z0**2+9*z0**4
     &                  -6*x0**2*(2+y0**2-3*z0**2))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_xy_ads_yz  =(32*(-1+L**2)*x0*z0*(-16*(x0**2-3*y0**2+z0**2)
     &                 -4*L**4*(-1+x0**2-9*y0**2+z0**2)*(-1+rho0**2)**4
     &                 -4*L**2*(-1+rho0**2)**2
     &                  *(1+7*x0**4-53*y0**4-8*z0**2+7*z0**4
     &                  +y0**2*(16-46*z0**2)
     &                  -2*x0**2*(4+23*y0**2-7*z0**2))
     &                 -16*(rho0**2)*(3*x0**4-4*z0**2
     &                  -2*x0**2*(2+9*y0**2-3*z0**2)
     &                  +3*(-7*y0**4+z0**4+y0**2*(4-6*z0**2)))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_xy_ads_zz  =(128*(-1+L**2)*x0*y0*(-4*(x0**2+y0**2-3*z0**2)
     &                 -L**4*(-1+x0**2+y0**2-9*z0**2)*(-1+rho0**2)**4
     &                 -4*(rho0**2)*(3*x0**4+3*y0**4
     &                  +2*x0**2*(-2+3*y0**2-9*z0**2)
     &                  +3*z0**2*(4-7*z0**2)-2*y0**2*(2+9*z0**2))
     &                 -L**2*(-1+rho0**2)**2*(1+7*x0**4+7*y0**4+16*z0**2
     &                  -53*z0**4+2*x0**2*(-4+7*y0**2-23*z0**2)
     &                  -2*y0**2*(4+23*z0**2))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)

        g0_xz_ads_x  =-((16*(-1+L**2)*z0
     &                *(L**2*(1+7*x0**2-y0**2-z0**2)*(-1+rho0**2)**2
     &                +4*(5*x0**4+y0**2+z0**2-(y0**2+z0**2)**2
     &                +x0**2*(-1+4*y0**2+4*z0**2))))
     &                /((-1+rho0**2)**3
     &                *(L**2*(-1+rho0**2)**2+4*(rho0**2))**2))
        g0_xz_ads_y  =-((128*(-1+L**2)*x0*y0*z0
     &                *(-1+3*rho0**2+L**2*(-1+rho0**2)**2))
     &                /((-1+rho0**2)**3
     &                *(L**2*(-1+rho0**2)**2+4*(rho0**2))**2))
        g0_xz_ads_z  =(16*(-1+L**2)*x0*(-4*(x0**2+y0**2-z0**2)
     &                +L**2*(-1+x0**2+y0**2-7*z0**2)*(-1+rho0**2)**2
     &                +4*(x0**2+y0**2-5*z0**2)*(rho0**2)))
     &                /((-1+rho0**2)**3
     &                *(L**2*(-1+rho0**2)**2+4*(rho0**2))**2)
        g0_xz_ads_xx  =(128*(-1+L**2)*x0*z0*(L**4*(-1+rho0**2)**4
     &                 *(7*x0**2-3*(-1+y0**2+z0**2))
     &                 +4*(x0**2-3*(y0**2+z0**2))
     &                 +4*(rho0**2)*(15*x0**4+12*(y0**2+z0**2)
     &                 -9*(y0**2+z0**2)**2
     &                 +x0**2*(-4+6*y0**2+6*z0**2))
     &                 +3*L**2*(-1+rho0**2)**2
     &                  *(-1+8*y0**2+8*z0**2+(rho0**2)
     &                   *(13*x0**2-7*(y0**2+z0**2)))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_xz_ads_xy  =(128*(-1+L**2)*y0*z0
     &                 *(12*x0**2-4*(y0**2+z0**2)
     &                 +L**4*(1+9*x0**2-y0**2-z0**2)*(-1+rho0**2)**4
     &                 +4*(rho0**2)*(21*x0**4+4*(y0**2+z0**2)
     &                 -3*(y0**2+z0**2)**2
     &                 +6*x0**2*(-2+3*y0**2+3*z0**2))
     &                 +L**2*(-1+rho0**2)**2
     &                  *(-1+53*x0**4+8*y0**2+8*z0**2-7*(y0**2+z0**2)**2
     &                   +2*x0**2*(-8+23*y0**2+23*z0**2))))
     &                 /((-1+rho0**2)**4 
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_xz_ads_xz  =(16*(-1+L**2)*(-L**4*(-1+rho0**2)**4
     &                 *(7*x0**4+6*x0**2*(-1+y0**2-11*z0**2)
     &                  -(-1+y0**2-7*z0**2)*(-1+y0**2+z0**2))
     &                 -8*L**2*(-1+rho0**2)**2
     &                 *(6*x0**6+x0**4*(-6+11*y0**2-42*z0**2)
     &                 -(-1+y0**2+z0**2)
     &                  *(y0**4-6*z0**4-y0**2*(1+5*z0**2))
     &                 +2*x0**2*(2*y0**4+6*z0**2-21*z0**4
     &                  -y0**2*(2+19*z0**2)))
     &                 +16*(-5*x0**8+x0**6*(6-14*y0**2+28*z0**2)
     &                  +(-1+y0**2+z0**2)*(y0**2+z0**2)
     &                  *(y0**4+z0**2-5*z0**4-y0**2*(1+4*z0**2))
     &                 +x0**4*(-1-12*y0**4-14*z0**2
     &                  +66*z0**4+2*y0**2*(5+27*z0**2))
     &                 +2*x0**2*(3*z0**2-(y0**2+z0**2)
     &                  *(y0**4+7*z0**2-14*z0**4-y0**2*(1+13*z0**2))))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_xz_ads_yy  =(128*(-1+L**2)*x0*z0*(-4*(x0**2-3*y0**2+z0**2)
     &                 -L**4*(-1+x0**2-9*y0**2+z0**2)*(-1+rho0**2)**4
     &                 -L**2*(-1+rho0**2)**2
     &                  *(1+7*x0**4-53*y0**4-8*z0**2+7*z0**4
     &                  +y0**2*(16-46*z0**2)
     &                  -2*x0**2*(4+23*y0**2-7*z0**2))
     &                 -4*(rho0**2)*(3*x0**4-4*z0**2
     &                  -2*x0**2*(2+9*y0**2-3*z0**2)
     &                 +3*(-7*y0**4+z0**4+y0**2*(4-6*z0**2)))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_xz_ads_yz  =(128*(-1+L**2)*x0*y0*(-4*(x0**2+y0**2-3*z0**2)
     &                 -L**4*(-1+x0**2+y0**2-9*z0**2)*(-1+rho0**2)**4
     &                 -4*(rho0**2)*(3*x0**4+3*y0**4
     &                  +2*x0**2*(-2+3*y0**2-9*z0**2)
     &                  +3*z0**2*(4-7*z0**2)-2*y0**2*(2+9*z0**2))
     &                 -L**2*(-1+rho0**2)**2*(1+7*x0**4+7*y0**4+16*z0**2
     &                  -53*z0**4+2*x0**2*(-4+7*y0**2-23*z0**2)
     &                  -2*y0**2*(4+23*z0**2))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_xz_ads_zz  =(128*(-1+L**2)*x0*z0
     &                 *(-L**4*(3*(-1+x0**2+y0**2)
     &                  -7*z0**2)*(-1+rho0**2)**4
     &                  +4*(-3*(x0**2+y0**2)+z0**2)
     &                 -4*(rho0**2)*(3*(x0**2+y0**2)
     &                  *(-4+3*x0**2+3*y0**2)
     &                 -2*(-2+3*x0**2+3*y0**2)*z0**2-15*z0**4)
     &                 -3*L**2*(-1+rho0**2)**2
     &                  *(1+7*x0**4+7*y0**4-13*z0**4
     &                  +2*x0**2*(-4+7*y0**2-3*z0**2)
     &                  -2*y0**2*(4+3*z0**2))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)

!!!CHECKED WITH MATHEMATICA UP TO HERE!!

        g0_yy_ads_x  =-((16*x0*(L**4*(-1+rho0**2)**2
     &                *(1+x0**4+6*y0**2-2*z0**2+2*x0**2*(-1+y0**2+z0**2)
     &                +(y0**2+z0**2)**2)
     &                +8*(y0**2+(rho0**2)*(2*x0**2-y0**2+2*z0**2))
     &                +8*L**2*(x0**6-y0**2+z0**2
     &                 +x0**4*(-2+2*y0**2+3*z0**2)
     &                 +(y0**2+z0**2)*(3*y0**2+(-2+y0**2)*z0**2+z0**4)
     &                +x0**2*(1+y0**2+y0**4
     &                 +4*(-1+y0**2)*z0**2+3*z0**4))))
     &                /((-1 + rho0**2)**3* (L**2* (-1 +rho0**2)**2
     &                + 4* (rho0**2))**2))
        g0_yy_ads_y  =-((16*y0*(8*(x0**2+z0**2)*(-1+3*rho0**2)
     &                +L**4*(-1+rho0**2)**2*(3+4*y0**2-4*z0**2
     &                 +((-2+x0)*x0+y0**2+z0**2)
     &                 *(x0*(2+x0)+y0**2+z0**2))
     &                +2*L**2*(-1+5*x0**6+3*y0**2+11*z0**2
     &                 +x0**4*(11*y0**2+15*(-1+z0**2))
     &                 +(y0**2+z0**2)*(y0**4+5*z0**2*(-3+z0**2)
     &                 +y0**2*(5+6*z0**2))
     &                 +x0**2*(11+7*y0**4+15*z0**2*(-2+z0**2)
     &                 +2*y0**2*(-5+11*z0**2)))))
     &                /((-1 + rho0**2)**3* (L**2* (-1 +rho0**2)**2
     &                + 4* (rho0**2))**2))
        g0_yy_ads_z  =-((16*z0*(L**4*(-1+rho0**2)**2
     &                *(1+x0**4+6*y0**2-2*z0**2+2*x0**2*(-1+y0**2+z0**2)
     &                 +(y0**2+z0**2)**2)
     &                +8*(y0**2+(rho0**2)*(2*x0**2-y0**2+2*z0**2))
     &                +8*L**2*(x0**6-y0**2+z0**2
     &                 +x0**4*(-2+2*y0**2+3*z0**2)
     &                 +(y0**2+z0**2)*(3*y0**2+(-2+y0**2)*z0**2+z0**4)
     &                 +x0**2*(1+y0**2+y0**4
     &                 +4*(-1+y0**2)*z0**2+3*z0**4))))
     &                /((-1 + rho0**2)**3* (L**2* (-1 +rho0**2)**2
     &                + 4* (rho0**2))**2))
        g0_yy_ads_xx  =8*((12*x0**2*(1+(-1+L**2)*y0**2))
     &                 /(-1+rho0**2)**4
     &                 +(-2-2*(-1+L**2)*(1+2*x0**2)*y0**2)
     &                 /(-1+rho0**2)**3
     &                 +((-1+L**2)*y0**2)
     &                 /(-1+rho0**2)**2
     &                 -(64*(-1+L**2)**2*x0**2*y0**2
     &                  *(4+L**2*(-2+rho0**2)))
     &                 /(L**2*(-1+rho0**2)**2
     &                  +4*(rho0**2))**3
     &                 -(L**2*(-1+L**2)*y0**2)
     &                 /(L**2*(-1+rho0**2)**2+4*(rho0**2))
     &                 +(2*(-1+L**2)*y0**2
     &                  *(-8+L**2*(14*x0**2-2*(-4+y0**2+z0**2)
     &                  +L**2*(-1+2*x0**4+y0**2+z0**2
     &                  +x0**2*(-7+2*y0**2+2*z0**2)))))
     &                 /(L**2*(-1+rho0**2)**2+4*(rho0**2))**2)
        g0_yy_ads_xy  =(32*x0*y0*(32*(x0**2-y0**2+z0**2)
     &                 +L**6*(-1+rho0**2)**4
     &                 *(11+3*x0**4+26*y0**2-14*z0**2+3*(y0**2+z0**2)**2
     &                  +2*x0**2*(-7+3*y0**2+3*z0**2))
     &                 +32*(rho0**2)*(9*x0**4-3*y0**4-4*z0**2+9*z0**4
     &                  +y0**2*(4+6*z0**2)+2*x0**2*(-2+3*y0**2+9*z0**2))
     &                 +4*L**4*(-1+rho0**2)**2
     &                 *(-4+11*x0**6-3*y0**2+31*z0**2
     &                  +x0**4*(-38+23*y0**2+33*z0**2)
     &                 +(y0**2+z0**2)*(y0**4-38*z0**2+11*z0**4
     &                  +6*y0**2*(7+2*z0**2))
     &                 +x0**2*(31+13*y0**4-76*z0**2+33*z0**4
     &                  +y0**2*(4+46*z0**2)))
     &                 +8*L**2*(1+25*x0**8+6*y0**2-14*z0**2
     &                  +10*x0**6*(-7+7*y0**2+10*z0**2)
     &                  +2*x0**4*(29-45*y0**2+30*y0**4
     &                  +105*(-1+y0**2)*z0**2+75*z0**4)
     &                 -(y0**2+z0**2)*(5*y0**6-58*z0**2+70*z0**4
     &                  -25*z0**6-5*y0**4*(10+3*z0**2)
     &                  +y0**2*(28+20*z0**2-45*z0**4))
     &                 +2*x0**2*(-7+58*z0**2+5*(y0**6
     &                  +3*y0**4*(1+4*z0**2)
     &                  +z0**4*(-21+10*z0**2)
     &                  +3*y0**2*(1-6*z0**2+7*z0**4))))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_yy_ads_xz  =(1/((-1+rho0**2)**4))
     &                 *32*x0*z0*(3+(8*(-1+L)*(1+L)*y0**2
     &                 *(5*L**4*(-1+rho0**2)**4+48*(rho0**2)**2
     &                  -8*(-1+4*rho0**2)
     &                  +6*L**2*(-1+rho0**2)**2*(-2+5*rho0**2)))
     &                 /(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_yy_ads_yy  =8 *((12* y0**2* (1 + (-1 + L**2)* y0**2))
     &                    /(-1 +rho0**2)**4
     &                 + (-2-2*(-1 + L**2)* y0**2*(5+2*y0**2))
     &                   /(-1 + rho0**2)**3
     &                 + ((-1+L**2)*(1+5* y0**2))
     &                   /(-1 +rho0**2)**2
     &                 + (1 - L**2)
     &                   /(-1 + rho0**2)
     &                 - (64*(-1+L**2)**2*y0**4*(4+L**2*(-2+rho0**2)))
     &                   /(L**2*(-1+rho0**2)**2+4*(rho0**2))**3
     &                 +(-1+L**2)*(4+L**2*(-2+x0**2-4*y0**2+z0**2))
     &                   /(L**2*(-1 +rho0**2)**2+4*(rho0**2))
     &                 +(2*(-1+L**2)*y0**2
     &                  *(-40+2*L**2*(20 - 5*x0**2+3*y0**2-5*z0**2)
     &                  +L**4*(-5-3*y0**2+2*y0**4+x0**2*(5+2*y0**2)
     &                  +(5+2*y0**2)*z0**2)))
     &                  /(L**2*(-1+rho0**2)**2 + 4* (rho0**2))**2)
        g0_yy_ads_yz  =(32*y0*z0*(32*(x0**2-y0**2+z0**2)
     &                 +L**6*(-1+rho0**2)**4
     &                 *(11+3*x0**4+26*y0**2-14*z0**2+3*(y0**2+z0**2)**2
     &                  +2*x0**2*(-7+3*y0**2+3*z0**2))
     &                 +32*(rho0**2)*(9*x0**4-3*y0**4-4*z0**2+9*z0**4
     &                  +y0**2*(4+6*z0**2)+2*x0**2*(-2+3*y0**2+9*z0**2))
     &                 +4*L**4*(-1+rho0**2)**2
     &                 *(-4+11*x0**6-3*y0**2+31*z0**2
     &                  +x0**4*(-38+23*y0**2+33*z0**2)
     &                 +(y0**2+z0**2)*(y0**4-38*z0**2+11*z0**4
     &                  +6*y0**2*(7+2*z0**2))
     &                 +x0**2*(31+13*y0**4-76*z0**2+33*z0**4
     &                  +y0**2*(4+46*z0**2)))
     &                 +8*L**2*(1+25*x0**8+6*y0**2-14*z0**2
     &                  +10*x0**6*(-7+7*y0**2+10*z0**2)
     &                  +2*x0**4*(29-45*y0**2+30*y0**4
     &                  +105*(-1+y0**2)*z0**2+75*z0**4)
     &                 -(y0**2+z0**2)*(5*y0**6-58*z0**2+70*z0**4
     &                  -25*z0**6-5*y0**4*(10+3*z0**2)
     &                  +y0**2*(28+20*z0**2-45*z0**4))
     &                 +2*x0**2*(-7+58*z0**2+5*(y0**6
     &                  +3*y0**4*(1+4*z0**2)
     &                  +z0**4*(-21+10*z0**2)
     &                  +3*y0**2*(1-6*z0**2+7*z0**4))))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_yy_ads_zz  =(16*(8*z0**2*(-1+rho0**2)*(2+L**2*(-1+rho0**2))
     &                 *(L**4*(-1+rho0**2)**2*(1+x0**4+6*y0**2-2*z0**2
     &                  +2*x0**2*(-1+y0**2+z0**2)
     &                  +(y0**2+z0**2)**2)
     &                  +8*(y0**2+(rho0**2)*(2*x0**2-y0**2+2*z0**2))
     &                 +8*L**2*(x0**6-y0**2+z0**2
     &                  +x0**4*(-2+2*y0**2+3*z0**2)
     &                  +(y0**2+z0**2)*(3*y0**2+(-2+y0**2)*z0**2+z0**4)
     &                 +x0**2*(1+y0**2+y0**4+4*(-1+y0**2)*z0**2
     &                  +3*z0**4)))
     &                 +6*z0**2*(L**2*(-1+rho0**2)**2
     &                  +4*(rho0**2))*(L**4*(-1+rho0**2)**2
     &                  *(1+x0**4+6*y0**2-2*z0**2
     &                  +2*x0**2*(-1+y0**2+z0**2)+(y0**2+z0**2)**2)
     &                 +8*(y0**2+(rho0**2)*(2*x0**2-y0**2+2*z0**2))
     &                 +8*L**2*(x0**6-y0**2+z0**2
     &                  +x0**4*(-2+2*y0**2+3*z0**2)
     &                  +(y0**2+z0**2)*(3*y0**2+(-2+y0**2)*z0**2+z0**4)
     &                  +x0**2*(1+y0**2+y0**4+4*(-1+y0**2)*z0**2
     &                  +3*z0**4)))
     &                 -(-1+rho0**2)*(L**2*(-1+rho0**2)**2+4*(rho0**2))
     &                 *(L**4*(-1+rho0**2)**2
     &                 *(1+x0**4+6*y0**2-2*z0**2
     &                  +2*x0**2*(-1+y0**2+z0**2)
     &                  +(y0**2+z0**2)**2)
     &                 +8*(y0**2+(rho0**2)*(2*x0**2-y0**2+2*z0**2))
     &                 +8*L**2*(x0**6-y0**2+z0**2
     &                  +x0**4*(-2+2*y0**2+3*z0**2)
     &                  +(y0**2+z0**2)*(3*y0**2+(-2+y0**2)*z0**2+z0**4)
     &                  +x0**2*(1+y0**2+y0**4+4*(-1+y0**2)*z0**2
     &                  +3*z0**4)))
     &                 -4*z0**2*(-1+rho0**2)*(L**2*(-1+rho0**2)**2
     &                 +4*(rho0**2))*(4*(4*x0**2+y0**2+4*z0**2)
     &                  +2*L**2*(2-8*x0**2+2*y0**2-8*z0**2+2*(rho0**2)
     &                  *(3*x0**2+y0**2+3*z0**2)
     &                 +L**2*(-1+rho0**2)*(x0**4+(1+y0**2)**2
     &                  +2*(-1+y0**2)*z0**2+z0**4
     &                  +2*x0**2*(-1+y0**2+z0**2))))))
     &                 /((-1+rho0**2)**4
     &                  *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)

        g0_yz_ads_x  =-((128*(-1+L**2)*x0*y0*z0
     &                *(-1+3*rho0**2+L**2*(-1+rho0**2)**2))
     &                /((-1+rho0**2)**3
     &                *(L**2*(-1+rho0**2)**2+4*(rho0**2))**2))
        g0_yz_ads_y  =-((16*(-1+L**2)*z0*(4*(x0**2-y0**2+z0**2)
     &                -L**2*(-1+x0**2-7*y0**2+z0**2)*(-1+rho0**2)**2
     &                -4*(x0**2-5*y0**2+z0**2)*(rho0**2)))
     &                /((-1+rho0**2)**3
     &                *(L**2*(-1+rho0**2)**2+4*(rho0**2))**2))
        g0_yz_ads_z  =(16*(-1+L**2)*y0*(-4*(x0**2+y0**2-z0**2)
     &                +L**2*(-1+x0**2+y0**2-7*z0**2)*(-1+rho0**2)**2
     &                +4*(x0**2+y0**2-5*z0**2)*(rho0**2)))
     &                /((-1+rho0**2)**3
     &                *(L**2*(-1+rho0**2)**2+4*(rho0**2))**2)
        g0_yz_ads_xx  =(128*(-1+L**2)*y0*z0
     &                 *(12*x0**2-4*(y0**2+z0**2)
     &                 +L**4*(1+9*x0**2-y0**2-z0**2)*(-1+rho0**2)**4
     &                 +4*(rho0**2)*(21*x0**4+4*(y0**2+z0**2)
     &                 -3*(y0**2+z0**2)**2
     &                 +6*x0**2*(-2+3*y0**2+3*z0**2))
     &                 +L**2*(-1+rho0**2)**2
     &                  *(-1+53*x0**4+8*y0**2+8*z0**2-7*(y0**2+z0**2)**2
     &                   +2*x0**2*(-8+23*y0**2+23*z0**2))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_yz_ads_xy  =(128*(-1+L**2)*x0*z0*(-4*(x0**2-3*y0**2+z0**2)
     &                 -L**4*(-1+x0**2-9*y0**2+z0**2)*(-1+rho0**2)**4
     &                 -L**2*(-1+rho0**2)**2
     &                  *(1+7*x0**4-53*y0**4-8*z0**2+7*z0**4
     &                  +y0**2*(16-46*z0**2)
     &                  -2*x0**2*(4+23*y0**2-7*z0**2))
     &                 -4*(rho0**2)*(3*x0**4-4*z0**2
     &                  -2*x0**2*(2+9*y0**2-3*z0**2)
     &                 +3*(-7*y0**4+z0**4+y0**2*(4-6*z0**2)))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_yz_ads_xz  =(128*(-1+L**2)*x0*y0*(-4*(x0**2+y0**2-3*z0**2)
     &                 -L**4*(-1+x0**2+y0**2-9*z0**2)*(-1+rho0**2)**4
     &                 -4*(rho0**2)*(3*x0**4+3*y0**4
     &                  +2*x0**2*(-2+3*y0**2-9*z0**2)
     &                  +3*z0**2*(4-7*z0**2)-2*y0**2*(2+9*z0**2))
     &                 -L**2*(-1+rho0**2)**2*(1+7*x0**4+7*y0**4+16*z0**2
     &                  -53*z0**4+2*x0**2*(-4+7*y0**2-23*z0**2)
     &                  -2*y0**2*(4+23*z0**2))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_yz_ads_yy  =(128*(-1+L**2)*y0*z0*(4*(-3*x0**2+y0**2-3*z0**2)
     &                 -L**4*(-1+rho0**2)**4
     &                  *(-3+3*x0**2-7*y0**2+3*z0**2)
     &                 -3*L**2*(-1+rho0**2)**2
     &                  *(1+7*x0**4-13*y0**4-2*(4+3*y0**2)*z0**2
     &                  +7*z0**4-2*x0**2*(4+3*y0**2-7*z0**2))
     &                 -4*(rho0**2)*(9*x0**4+4*y0**2-15*y0**4
     &                  -6*(2+y0**2)*z0**2
     &                  +9*z0**4-6*x0**2*(2+y0**2-3*z0**2))))
     &                 /((-1+rho0**2)**4
     &                  *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_yz_ads_yz  =(16*(-1+L**2)*(L**4*(-1+rho0**2)**4
     &                 *(1+x0**4-7*y0**4+6*z0**2-7*z0**4
     &                 -2*x0**2*(1+3*y0**2+3*z0**2)+y0**2*(6+66*z0**2))
     &                 +16*(x0**8-5*y0**8-z0**4+6*z0**6-5*z0**8
     &                  -2*x0**6*(1+y0**2+z0**2)+y0**6*(6+28*z0**2)
     &                  +2*y0**2*z0**2*(3-7*z0**2+14*z0**4)
     &                  +y0**4*(-1-14*z0**2+66*z0**4)
     &                  +x0**4*(1-12*y0**4+2*z0**2-12*z0**4
     &                  +y0**2*(2+24*z0**2))
     &                 -2*x0**2*(7*y0**6+3*y0**2*z0**2*(2-9*z0**2)
     &                  +z0**4*(-5+7*z0**2)-y0**4*(5+27*z0**2)))
     &                 +8*L**2*(-1+rho0**2)**2
     &                 *(x0**6-2*x0**4*(1+2*y0**2+2*z0**2)
     &                  +6*(-y0**6+z0**4-z0**6+y0**2*z0**2*(-2+7*z0**2)
     &                  +y0**4*(1+7*z0**2))
     &                 +x0**2*(1-11*y0**4+4*z0**2-11*z0**4
     &                  +y0**2*(4+38*z0**2)))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_yz_ads_zz  =(128*(-1+L**2)*y0*z0
     &                 *(-L**4*(3*(-1+x0**2+y0**2)
     &                  -7*z0**2)*(-1+rho0**2)**4
     &                  +4*(-3*(x0**2+y0**2)+z0**2)
     &                 -4*(rho0**2)*(3*(x0**2+y0**2)
     &                  *(-4+3*x0**2+3*y0**2)
     &                 -2*(-2+3*x0**2+3*y0**2)*z0**2-15*z0**4)
     &                 -3*L**2*(-1+rho0**2)**2
     &                  *(1+7*x0**4+7*y0**4-13*z0**4
     &                  +2*x0**2*(-4+7*y0**2-3*z0**2)
     &                  -2*y0**2*(4+3*z0**2))))
     &                 /((-1+rho0**2)**4
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)

        g0_zz_ads_x  =-((16*x0*(16*(x0**2+y0**2)**2
     &                 +8*(1+x0**2+y0**2)*z0**2-8*z0**4
     &                 +L**4*(-1+rho0**2)**2*((-1+x0**2+y0**2)**2
     &                  +2*(3+x0**2+y0**2)*z0**2+z0**4)
     &                 +8*L**2*((-1+x0**2+y0**2)**2*(x0**2+y0**2)
     &                 +(1+x0**2+y0**2)*(-1+2*x0**2+2*y0**2)*z0**2
     &                 +(3+x0**2+y0**2)*z0**4)))
     &                 /((-1+rho0**2)**3
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**2))
        g0_zz_ads_y  =-((16*y0*(16*(x0**2+y0**2)**2
     &                 +8*(1+x0**2+y0**2)*z0**2-8*z0**4
     &                 +L**4*(-1+rho0**2)**2*((-1+x0**2+y0**2)**2
     &                  +2*(3+x0**2+y0**2)*z0**2+z0**4)
     &                 +8*L**2*((-1+x0**2+y0**2)**2*(x0**2+y0**2)
     &                 +(1+x0**2+y0**2)*(-1+2*x0**2+2*y0**2)*z0**2
     &                 +(3+x0**2+y0**2)*z0**4)))
     &                 /((-1+rho0**2)**3
     &                 *(L**2*(-1+rho0**2)**2+4*(rho0**2))**2))
        g0_zz_ads_z  =-((16*z0*(8*(x0**2+y0**2)*(-1+3*rho0**2)
     &                +L**4*(-1+rho0**2)**2*(3-4*y0**2+4*z0**2
     &                 +((-2+x0)*x0+y0**2+z0**2)
     &                 *(x0*(2+x0)+y0**2+z0**2))
     &                +2*L**2*(-1+5*x0**6+11*y0**2+3*z0**2
     &                 +x0**4*(-15+15*y0**2+11*z0**2)
     &                 +(y0**2+z0**2)*(5*y0**2*(-3+y0**2)
     &                 +(5+6*y0**2)*z0**2+z0**4)
     &                 +x0**2*(11+15*y0**4-10*z0**2+7*z0**4
     &                 +y0**2*(-30+22*z0**2)))))
     &                /((-1 + rho0**2)**3* (L**2* (-1 +rho0**2)**2
     &                + 4* (rho0**2))**2))
        g0_zz_ads_xx  =(16*(4*L**4*(-1+rho0**2)**2*(3*(1+5*x0**2-y0**2)
     &                   *(-1+x0**2+y0**2)**2*(x0**2+y0**2)
     &                  +(-1+24*x0**6+4*y0**2+7*y0**4-10*y0**6
     &                  +19*x0**4*(5+2*y0**2)
     &                  +2*x0**2*(-23+51*y0**2+2*y0**4))*z0**2
     &                  +(13+2*x0**4-13*y0**2-12*y0**4
     &                   +x0**2*(111-10*y0**2))*z0**4
     &                   -(11+8*x0**2+6*y0**2)*z0**6-z0**8)
     &                  +32*(2*(1+5*x0**2-y0**2)*(x0**2+y0**2)**3
     &                  +(7*x0**6+y0**2+2*y0**4-5*y0**6
     &                   +9*x0**4*(2+y0**2)
     &                  +x0**2*(-3+20*y0**2-3*y0**4))*z0**2
     &                  -(-1+15*x0**4+2*y0**2+3*y0**4
     &                   +2*x0**2*(-7+9*y0**2))*z0**4
     &                   +(-2-11*x0**2+y0**2)*z0**6+z0**8)
     &                  +8*L**2*(6*(1+5*x0**2-y0**2)
     &                  *(-1+x0**2+y0**2)**2*(x0**2+y0**2)**2
     &                  +(1+61*x0**8-2*y0**2-14*y0**4+38*y0**6-23*y0**8
     &                   +2*x0**6*(31+80*y0**2)
     &                   +6*x0**4*(-19+27*y0**2+19*y0**4)
     &                   +2*x0**2*(19-64*y0**2+69*y0**4-4*y0**6))*z0**2
     &                  +2*(-4+2*x0**6+13*y0**2+3*y0**4-16*y0**6
     &                   +3*x0**4*(45-4*y0**2)
     &                   +x0**2*(-55+138*y0**2-30*y0**4))*z0**4
     &                  -2*(-11+27*x0**4+15*y0**2+9*y0**4
     &                   +x0**2*(-69+36*y0**2))*z0**6
     &                   -2*(8+13*x0**2+y0**2)*z0**8+z0**10)
     &                  +L**6*(-1+rho0**2)**4*(5*x0**6
     &                   +9*x0**4*(-1+y0**2+z0**2)
     &                   -(-1+y0**2+z0**2)*((-1+y0**2)**2
     &                   +2*(3+y0**2)*z0**2+z0**4)
     &                   +3*x0**2*((-1+y0**2)**2+2*(11+y0**2)*z0**2
     &                   + z0**4))))
     &                  /((-1+rho0**2)**4
     &                  *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_zz_ads_xy  =(32*x0*y0*(-64*z0**2+64*(rho0**2)
     &                   *(3*(x0**2+y0**2)**2+4*z0**2-3*z0**4)
     &                  +L**6*(-1+rho0**2)**4*(3*(-1+x0**2+y0**2)**2
     &                  +2*(17+3*x0**2+3*y0**2)*z0**2+3*z0**4)
     &                  +4*L**4*(-1+rho0**2)**2
     &                  *(9*(-1+x0**2+y0**2)**2*(x0**2+y0**2)
     &                  +(-25+17*x0**4+44*y0**2+17*y0**4
     &                   +x0**2*(44+34*y0**2))*z0**2
     &                  +(62+7*x0**2+7*y0**2)*z0**4-z0**6)
     &                  +16*L**2*(10*z0**2+(rho0**2)
     &                   *(9*(-1+x0**2+y0**2)**2*(x0**2+y0**2)
     &                  +2*(-17+12*y0**2+6*(x0**4+y0**4
     &                   +2*x0**2*(1+y0**2)))*z0**2
     &                   -3*(-14+x0**2+y0**2)*z0**4-6*z0**6))))
     &                  /((-1+rho0**2)**4
     &                  *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_zz_ads_xz  =(16*x0*(8*z0*(-1+rho0**2)*(2+L**2*(-1+rho0**2))
     &                  *(16*(x0**2+y0**2)**2
     &                  +8*(1+x0**2+y0**2)*z0**2-8*z0**4
     &                  +L**4*(-1+rho0**2)**2*((-1+x0**2+y0**2)**2
     &                   +2*(3+x0**2+y0**2)*z0**2+z0**4)
     &                  +8*L**2*((-1+x0**2+y0**2)**2*(x0**2+y0**2)
     &                  +(1+x0**2+y0**2)*(-1+2*x0**2+2*y0**2)*z0**2
     &                  +(3+x0**2+y0**2)*z0**4))
     &                  +6*z0*(L**2*(-1+rho0**2)**2+4*(rho0**2))
     &                   *(16*(x0**2+y0**2)**2
     &                  +8*(1+x0**2+y0**2)*z0**2-8*z0**4
     &                  +L**4*(-1+rho0**2)**2*((-1+x0**2+y0**2)**2
     &                   +2*(3+x0**2+y0**2)*z0**2+z0**4)
     &                  +8*L**2*((-1+x0**2+y0**2)**2*(x0**2+y0**2)
     &                   +(1+x0**2+y0**2)*(-1+2*x0**2+2*y0**2)*z0**2
     &                   +(3+x0**2+y0**2)*z0**4))
     &                  -8*z0*(-1+rho0**2)*(L**2*(-1+rho0**2)**2
     &                   +4*(rho0**2))*(2*(1+x0**2+y0**2-2*z0**2)
     &                  +L**2*(4*(x0**2+y0**2)*(rho0**2)
     &                  +2*(-1+x0**2+y0**2+6*z0**2)
     &                  +L**2*(-1+rho0**2)*(-1+4*z0**2+(rho0**2)**2)))))
     &                  /((-1+rho0**2)**4
     &                  *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_zz_ads_yy  =(16*(-L**6*(-1+rho0**2)**4*((-1+x0**2-5*y0**2)
     &                   *(-1+x0**2+y0**2)**2
     &                  +(-5+2*x0**2+3*x0**4-6*(11+x0**2)*y0**2
     &                   -9*y0**4)*z0**2
     &                  +(5+3*x0**2-3*y0**2)*z0**4+z0**6)
     &                  +32*(-2*(-1+x0**2-5*y0**2)*(x0**2+y0**2)**3
     &                  +(x0**2+2*x0**4-5*x0**6
     &                  +(-3+20*x0**2-3*x0**4)*y0**2+9*(2+x0**2)*y0**4
     &                   +7*y0**6)*z0**2
     &                  -(-1+2*x0**2+3*x0**4+2*(-7+9*x0**2)*y0**2
     &                   +15*y0**4)*z0**4
     &                   +(-2+x0**2-11*y0**2)*z0**6+z0**8)
     &                  -4*L**4*(-1+rho0**2)**2
     &                   *(3*(-1+x0**2-5*y0**2)*(-1+x0**2+y0**2)**2
     &                   *(x0**2+y0**2)
     &                  +(1+10*x0**6+46*y0**2-95*y0**4-24*y0**6
     &                   -x0**4*(7+4*y0**2)
     &                   -2*x0**2*(2+51*y0**2+19*y0**4))*z0**2
     &                  +(-13+12*x0**4-111*y0**2-2*y0**4
     &                   +x0**2*(13+10*y0**2))*z0**4
     &                  +(11+6*x0**2+8*y0**2)*z0**6+z0**8)
     &                  +8*L**2*(-6*(-1+x0**2-5*y0**2)
     &                   *(-1+x0**2+y0**2)**2*(x0**2+y0**2)**2
     &                   +(1-23*x0**8+38*y0**2-114*y0**4+62*y0**6
     &                   +61*y0**8+x0**6*(38-8*y0**2)
     &                   +2*x0**4*(-7+69*y0**2+57*y0**4)
     &                   +2*x0**2*(-1-64*y0**2+81*y0**4+80*y0**6))*z0**2
     &                  +2*(-4+13*x0**2+3*x0**4-16*x0**6
     &                   +(-55+138*x0**2-30*x0**4)*y0**2
     &                   +3*(45-4*x0**2)*y0**4+2*y0**6)*z0**4
     &                  -2*(-11+9*x0**4-69*y0**2+27*y0**4
     &                   +3*x0**2*(5+12*y0**2))*z0**6
     &                  -2*(8+x0**2+13*y0**2)*z0**8+z0**10)))
     &                  /((-1+rho0**2)**4
     &                  *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_zz_ads_yz  =(16*y0*(8*z0*(-1+rho0**2)*(2+L**2*(-1+rho0**2))
     &                  *(16*(x0**2+y0**2)**2
     &                  +8*(1+x0**2+y0**2)*z0**2-8*z0**4
     &                  +L**4*(-1+rho0**2)**2*((-1+x0**2+y0**2)**2
     &                   +2*(3+x0**2+y0**2)*z0**2+z0**4)
     &                  +8*L**2*((-1+x0**2+y0**2)**2*(x0**2+y0**2)
     &                  +(1+x0**2+y0**2)*(-1+2*x0**2+2*y0**2)*z0**2
     &                  +(3+x0**2+y0**2)*z0**4))
     &                  +6*z0*(L**2*(-1+rho0**2)**2+4*(rho0**2))
     &                   *(16*(x0**2+y0**2)**2
     &                  +8*(1+x0**2+y0**2)*z0**2-8*z0**4
     &                  +L**4*(-1+rho0**2)**2*((-1+x0**2+y0**2)**2
     &                   +2*(3+x0**2+y0**2)*z0**2+z0**4)
     &                  +8*L**2*((-1+x0**2+y0**2)**2*(x0**2+y0**2)
     &                   +(1+x0**2+y0**2)*(-1+2*x0**2+2*y0**2)*z0**2
     &                   +(3+x0**2+y0**2)*z0**4))
     &                  -8*z0*(-1+rho0**2)*(L**2*(-1+rho0**2)**2
     &                   +4*(rho0**2))*(2*(1+x0**2+y0**2-2*z0**2)
     &                  +L**2*(4*(x0**2+y0**2)*(rho0**2)
     &                  +2*(-1+x0**2+y0**2+6*z0**2)
     &                  +L**2*(-1+rho0**2)*(-1+4*z0**2+(rho0**2)**2)))))
     &                  /((-1+rho0**2)**4
     &                  *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)
        g0_zz_ads_zz  =(16*(-L**6*(-1+rho0**2)**4
     &                   *((-3+x0**2+y0**2)*(-1+x0**2+y0**2)**2
     &                  -3*(-13+x0**2+y0**2)*(-1+x0**2+y0**2)*z0**2
     &                  -3*(11+3*x0**2+3*y0**2)*z0**4-5*z0**6)
     &                  +32*(x0**2+y0**2)*(-(-1+x0**2+y0**2)
     &                   *(x0**2+y0**2)*(-1+3*x0**2+3*y0**2)
     &                  +(3+15*x0**4-8*y0**2+15*y0**4
     &                   +x0**2*(-8+30*y0**2))*z0**2
     &                   +3*(-4+13*x0**2+13*y0**2)*z0**4
     &                   +21*z0**6)
     &                  -2*L**4*(-1+rho0**2)**2
     &                  *((-1+x0**2+y0**2)**2
     &                  *(1+7*x0**4-16*y0**2+7*y0**4
     &                   +2*x0**2*(-8+7*y0**2))
     &                  -2*(-1+x0**2+y0**2)
     &                  *(11+14*x0**4-77*y0**2+14*y0**4
     &                   +7*x0**2*(-11+4*y0**2))*z0**2
     &                  -2*(40+43*x0**4-67*y0**2+43*y0**4
     &                   +x0**2*(-67+86*y0**2))*z0**4
     &                  -6*(13+10*x0**2+10*y0**2)*z0**6-9*z0**8)
     &                  +8*L**2*(-2*(-1+x0**2+y0**2)**2*(x0**2+y0**2)
     &                  *(1+4*x0**4-7*y0**2+4*y0**4+x0**2*(-7+8*y0**2))
     &                  +(-1+x0**2+y0**2)
     &                  *(-3+31*x0**6+31*y0**2-91*y0**4+31*y0**6
     &                   +x0**4*(-91+93*y0**2)
     &                   +x0**2*(31-182*y0**2+93*y0**4))*z0**2
     &                   +6*(-2+24*x0**6+31*y0**2-51*y0**4+24*y0**6
     &                   +x0**4*(-51+72*y0**2)
     &                   +x0**2*(31-102*y0**2+72*y0**4))*z0**4
     &                  +2*(13+83*x0**4-63*y0**2+83*y0**4
     &                   +x0**2*(-63+166*y0**2))*z0**6
     &                  +4*(7+16*x0**2+16*y0**2)*z0**8+3*z0**10)))
     &                  /((-1+rho0**2)**4
     &                  *(L**2*(-1+rho0**2)**2+4*(rho0**2))**3)

!!!CHECKED WITH MATHEMATICA UP TO HERE!!

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0=',L,i,j,k,x0,y0,z0,rho0
!
!       write (*,*) ' g0_tt_ads_x=' ,g0_tt_ads_x
!       write (*,*) ' g0_tt_ads_y=' ,g0_tt_ads_y
!       write (*,*) ' g0_tt_ads_z=' ,g0_tt_ads_z
!       write (*,*) ' g0_tt_ads_xx=',g0_tt_ads_xx
!       write (*,*) ' g0_tt_ads_xy=',g0_tt_ads_xy
!       write (*,*) ' g0_tt_ads_xz=',g0_tt_ads_xz
!       write (*,*) ' g0_tt_ads_yy=',g0_tt_ads_yy
!       write (*,*) ' g0_tt_ads_yz=',g0_tt_ads_yz
!       write (*,*) ' g0_tt_ads_zz=',g0_tt_ads_zz
!
!       write (*,*) ' g0_tx_ads_x=' ,g0_tx_ads_x
!       write (*,*) ' g0_tx_ads_y=' ,g0_tx_ads_y
!       write (*,*) ' g0_tx_ads_z=' ,g0_tx_ads_z
!       write (*,*) ' g0_tx_ads_xx=',g0_tx_ads_xx
!       write (*,*) ' g0_tx_ads_xy=',g0_tx_ads_xy
!       write (*,*) ' g0_tx_ads_xz=',g0_tx_ads_xz
!       write (*,*) ' g0_tx_ads_yy=',g0_tx_ads_yy
!       write (*,*) ' g0_tx_ads_yz=',g0_tx_ads_yz
!       write (*,*) ' g0_tx_ads_zz=',g0_tx_ads_zz
!
!       write (*,*) ' g0_ty_ads_x=' ,g0_ty_ads_x
!       write (*,*) ' g0_ty_ads_y=' ,g0_ty_ads_y
!       write (*,*) ' g0_ty_ads_z=' ,g0_ty_ads_z
!       write (*,*) ' g0_ty_ads_xx=',g0_ty_ads_xx
!       write (*,*) ' g0_ty_ads_xy=',g0_ty_ads_xy
!       write (*,*) ' g0_ty_ads_xz=',g0_ty_ads_xz
!       write (*,*) ' g0_ty_ads_yy=',g0_ty_ads_yy
!       write (*,*) ' g0_ty_ads_yz=',g0_ty_ads_yz
!       write (*,*) ' g0_ty_ads_zz=',g0_ty_ads_zz
!
!       write (*,*) ' g0_tz_ads_x=' ,g0_tz_ads_x
!       write (*,*) ' g0_tz_ads_y=' ,g0_tz_ads_y
!       write (*,*) ' g0_tz_ads_z=' ,g0_tz_ads_z
!       write (*,*) ' g0_tz_ads_xx=',g0_tz_ads_xx
!       write (*,*) ' g0_tz_ads_xy=',g0_tz_ads_xy
!       write (*,*) ' g0_tz_ads_xz=',g0_tz_ads_xz
!       write (*,*) ' g0_tz_ads_yy=',g0_tz_ads_yy
!       write (*,*) ' g0_tz_ads_yz=',g0_tz_ads_yz
!       write (*,*) ' g0_tz_ads_zz=',g0_tz_ads_zz
!
!       write (*,*) ' g0_xx_ads_x=' ,g0_xx_ads_x
!       write (*,*) ' g0_xx_ads_y=' ,g0_xx_ads_y
!       write (*,*) ' g0_xx_ads_z=' ,g0_xx_ads_z
!       write (*,*) ' g0_xx_ads_xx=',g0_xx_ads_xx
!       write (*,*) ' g0_xx_ads_xy=',g0_xx_ads_xy
!       write (*,*) ' g0_xx_ads_xz=',g0_xx_ads_xz
!       write (*,*) ' g0_xx_ads_yy=',g0_xx_ads_yy
!       write (*,*) ' g0_xx_ads_yz=',g0_xx_ads_yz
!       write (*,*) ' g0_xx_ads_zz=',g0_xx_ads_zz
!
!       write (*,*) ' g0_xy_ads_x=' ,g0_xy_ads_x
!       write (*,*) ' g0_xy_ads_y=' ,g0_xy_ads_y
!       write (*,*) ' g0_xy_ads_z=' ,g0_xy_ads_z
!       write (*,*) ' g0_xy_ads_xx=',g0_xy_ads_xx
!       write (*,*) ' g0_xy_ads_xy=',g0_xy_ads_xy
!       write (*,*) ' g0_xy_ads_xz=',g0_xy_ads_xz
!       write (*,*) ' g0_xy_ads_yy=',g0_xy_ads_yy
!       write (*,*) ' g0_xy_ads_yz=',g0_xy_ads_yz
!       write (*,*) ' g0_xy_ads_zz=',g0_xy_ads_zz
!
!       write (*,*) ' g0_xz_ads_x=' ,g0_xz_ads_x
!       write (*,*) ' g0_xz_ads_y=' ,g0_xz_ads_y
!       write (*,*) ' g0_xz_ads_z=' ,g0_xz_ads_z
!       write (*,*) ' g0_xz_ads_xx=',g0_xz_ads_xx
!       write (*,*) ' g0_xz_ads_xy=',g0_xz_ads_xy
!       write (*,*) ' g0_xz_ads_xz=',g0_xz_ads_xz
!       write (*,*) ' g0_xz_ads_yy=',g0_xz_ads_yy
!       write (*,*) ' g0_xz_ads_yz=',g0_xz_ads_yz
!       write (*,*) ' g0_xz_ads_zz=',g0_xz_ads_zz
!
!       write (*,*) ' g0_yy_ads_x=' ,g0_yy_ads_x
!       write (*,*) ' g0_yy_ads_y=' ,g0_yy_ads_y
!       write (*,*) ' g0_yy_ads_z=' ,g0_yy_ads_z
!       write (*,*) ' g0_yy_ads_xx=',g0_yy_ads_xx
!       write (*,*) ' g0_yy_ads_xy=',g0_yy_ads_xy
!       write (*,*) ' g0_yy_ads_xz=',g0_yy_ads_xz
!       write (*,*) ' g0_yy_ads_yy=',g0_yy_ads_yy
!       write (*,*) ' g0_yy_ads_yz=',g0_yy_ads_yz
!       write (*,*) ' g0_yy_ads_zz=',g0_yy_ads_zz
!
!       write (*,*) ' g0_yz_ads_x=' ,g0_yz_ads_x
!       write (*,*) ' g0_yz_ads_y=' ,g0_yz_ads_y
!       write (*,*) ' g0_yz_ads_z=' ,g0_yz_ads_z
!       write (*,*) ' g0_yz_ads_xx=',g0_yz_ads_xx
!       write (*,*) ' g0_yz_ads_xy=',g0_yz_ads_xy
!       write (*,*) ' g0_yz_ads_xz=',g0_yz_ads_xz
!       write (*,*) ' g0_yz_ads_yy=',g0_yz_ads_yy
!       write (*,*) ' g0_yz_ads_yz=',g0_yz_ads_yz
!       write (*,*) ' g0_yz_ads_zz=',g0_yz_ads_zz
!
!       write (*,*) ' g0_zz_ads_x=' ,g0_zz_ads_x
!       write (*,*) ' g0_zz_ads_y=' ,g0_zz_ads_y
!       write (*,*) ' g0_zz_ads_z=' ,g0_zz_ads_z
!       write (*,*) ' g0_zz_ads_xx=',g0_zz_ads_xx
!       write (*,*) ' g0_zz_ads_xy=',g0_zz_ads_xy
!       write (*,*) ' g0_zz_ads_xz=',g0_zz_ads_xz
!       write (*,*) ' g0_zz_ads_yy=',g0_zz_ads_yy
!       write (*,*) ' g0_zz_ads_yz=',g0_zz_ads_yz
!       write (*,*) ' g0_zz_ads_zz=',g0_zz_ads_zz

        ! calculate gbar derivatives
        call df2_int(gb_tt_np1,gb_tt_n,gb_tt_nm1,gb_tt_t,
     &       gb_tt_x,gb_tt_y,
     &       gb_tt_z,
     &       gb_tt_tt,gb_tt_tx,gb_tt_ty,
     &       gb_tt_tz,
     &       gb_tt_xx,gb_tt_xy,
     &       gb_tt_xz,
     &       gb_tt_yy,
     &       gb_tt_yz,
     &       gb_tt_zz,
     &       x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'gb_tt')
        call df2_int(gb_tx_np1,gb_tx_n,gb_tx_nm1,gb_tx_t,
     &       gb_tx_x,gb_tx_y,
     &       gb_tx_z,
     &       gb_tx_tt,gb_tx_tx,gb_tx_ty,
     &       gb_tx_tz,
     &       gb_tx_xx,gb_tx_xy,
     &       gb_tx_xz,
     &       gb_tx_yy,
     &       gb_tx_yz,
     &       gb_tx_zz,
     &       x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'gb_tx')
        call df2_int(gb_ty_np1,gb_ty_n,gb_ty_nm1,gb_ty_t,
     &       gb_ty_x,gb_ty_y,
     &       gb_ty_z,
     &       gb_ty_tt,gb_ty_tx,gb_ty_ty,
     &       gb_ty_tz,
     &       gb_ty_xx,gb_ty_xy,
     &       gb_ty_xz,
     &       gb_ty_yy,
     &       gb_ty_yz,
     &       gb_ty_zz,
     &       x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'gb_ty')
        call df2_int(gb_tz_np1,gb_tz_n,gb_tz_nm1,gb_tz_t,
     &       gb_tz_x,gb_tz_y,
     &       gb_tz_z,
     &       gb_tz_tt,gb_tz_tx,gb_tz_ty,
     &       gb_tz_tz,
     &       gb_tz_xx,gb_tz_xy,
     &       gb_tz_xz,
     &       gb_tz_yy,
     &       gb_tz_yz,
     &       gb_tz_zz,
     &       x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'gb_tz')
        call df2_int(gb_xx_np1,gb_xx_n,gb_xx_nm1,gb_xx_t,
     &       gb_xx_x,gb_xx_y,
     &       gb_xx_z,
     &       gb_xx_tt,gb_xx_tx,gb_xx_ty,
     &       gb_xx_tz,
     &       gb_xx_xx,gb_xx_xy,
     &       gb_xx_xz,
     &       gb_xx_yy,
     &       gb_xx_yz,
     &       gb_xx_zz,
     &       x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'gb_xx')
        call df2_int(gb_xy_np1,gb_xy_n,gb_xy_nm1,gb_xy_t,
     &       gb_xy_x,gb_xy_y,
     &       gb_xy_z,
     &       gb_xy_tt,gb_xy_tx,gb_xy_ty,
     &       gb_xy_tz,
     &       gb_xy_xx,gb_xy_xy,
     &       gb_xy_xz,
     &       gb_xy_yy,
     &       gb_xy_yz,
     &       gb_xy_zz,
     &       x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'gb_xy')
        call df2_int(gb_xz_np1,gb_xz_n,gb_xz_nm1,gb_xz_t,
     &       gb_xz_x,gb_xz_y,
     &       gb_xz_z,
     &       gb_xz_tt,gb_xz_tx,gb_xz_ty,
     &       gb_xz_tz,
     &       gb_xz_xx,gb_xz_xy,
     &       gb_xz_xz,
     &       gb_xz_yy,
     &       gb_xz_yz,
     &       gb_xz_zz,
     &       x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'gb_xz')
        call df2_int(gb_yy_np1,gb_yy_n,gb_yy_nm1,gb_yy_t,
     &       gb_yy_x,gb_yy_y,
     &       gb_yy_z,
     &       gb_yy_tt,gb_yy_tx,gb_yy_ty,
     &       gb_yy_tz,
     &       gb_yy_xx,gb_yy_xy,
     &       gb_yy_xz,
     &       gb_yy_yy,
     &       gb_yy_yz,
     &       gb_yy_zz,
     &       x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'gb_yy')
        call df2_int(gb_yz_np1,gb_yz_n,gb_yz_nm1,gb_yz_t,
     &       gb_yz_x,gb_yz_y,
     &       gb_yz_z,
     &       gb_yz_tt,gb_yz_tx,gb_yz_ty,
     &       gb_yz_tz,
     &       gb_yz_xx,gb_yz_xy,
     &       gb_yz_xz,
     &       gb_yz_yy,
     &       gb_yz_yz,
     &       gb_yz_zz,
     &       x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'gb_yz')
        call df2_int(gb_zz_np1,gb_zz_n,gb_zz_nm1,gb_zz_t,gb_zz_x,
     &       gb_zz_y,
     &       gb_zz_z,
     &       gb_zz_tt,gb_zz_tx,gb_zz_ty,
     &       gb_zz_tz,
     &       gb_zz_xx,gb_zz_xy,
     &       gb_zz_xz,
     &       gb_zz_yy,
     &       gb_zz_yz,
     &       gb_zz_zz,
     &       x,y,z,dt,i,j,k,chr,ex,Nx,Ny,Nz,'gb_zz')

        ! calculate hbar derivatives
        call df1_int(Hb_t_np1,Hb_t_n,Hb_t_nm1,Hb_t_t,Hb_t_x,
     &       Hb_t_y,
     &       Hb_t_z,
     &       x,y,z,dt,i,j,k,
     &       chr,ex,Nx,Ny,Nz,'Hb_t')
        call df1_int(Hb_x_np1,Hb_x_n,Hb_x_nm1,Hb_x_t,Hb_x_x,
     &       Hb_x_y,
     &       Hb_x_z,
     &       x,y,z,dt,i,j,k,
     &       chr,ex,Nx,Ny,Nz,'Hb_x')
        call df1_int(Hb_y_np1,Hb_y_n,Hb_y_nm1,Hb_y_t,Hb_y_x,
     &       Hb_y_y,
     &       Hb_y_z,
     &       x,y,z,dt,i,j,k,
     &       chr,ex,Nx,Ny,Nz,'Hb_y')
        call df1_int(Hb_z_np1,Hb_z_n,Hb_z_nm1,Hb_z_t,Hb_z_x,
     &       Hb_z_y,
     &       Hb_z_z,
     &       x,y,z,dt,i,j,k,
     &       chr,ex,Nx,Ny,Nz,'Hb_z')

        ! calculate phi1 derivatives
        call df2_int(phi1_np1,phi1_n,phi1_nm1,phi1_t,phi1_x,
     &       phi1_y,
     &       phi1_z,
     &       phi1_tt,phi1_tx,phi1_ty,
     &       phi1_tz,
     &       phi1_xx,phi1_xy,
     &       phi1_xz,
     &       phi1_yy,
     &       phi1_yz,
     &       phi1_zz,
     &       x,y,z,dt,i,j,k,
     &       chr,ex,Nx,Ny,Nz,'phi1')


!!!!!!!!!!DEBUG DERIVATIVE STENCILS!!!!!!!!!!!
!         if    ((gb_tt_z.ne.0.0d0)
!     &      .or.(gb_tt_tz.ne.0.0d0)
!     &      .or.(gb_tt_xz.ne.0.0d0)
!     &      .or.(gb_tt_yz.ne.0.0d0)
!     &      .or.(gb_tt_zz.ne.0.0d0)) then
!
!          write(*,*) 'DEBUG from misc.f: non zero z-derivative of 
!     &           gb_tt'
!
!          write(*,*) 'gb_tt_nm1(i,j,k),gb_tt_n(i,j,k),gb_tt_np1(i,j,k)='
!     &               ,gb_tt_nm1(i,j,k),gb_tt_n(i,j,k),gb_tt_np1(i,j,k)
!
!          write(*,*) 'gb_tt_z,gb_tt_tz,gb_tt_xz,gb_tt_yz,gb_tt_zz='
!     &               ,gb_tt_z,gb_tt_tz,gb_tt_xz,gb_tt_yz,gb_tt_zz
!
!         end if
!
!         if    ((gb_tx_z.ne.0.0d0)
!     &      .or.(gb_tx_tz.ne.0.0d0)
!     &      .or.(gb_tx_xz.ne.0.0d0)
!     &      .or.(gb_tx_yz.ne.0.0d0)
!     &      .or.(gb_tx_zz.ne.0.0d0)) then
!
!          write(*,*) 'DEBUG from misc.f: non zero z-derivative of 
!     &           gb_tx'
!
!         end if
!
!         if    ((gb_ty_z.ne.0.0d0)
!     &      .or.(gb_ty_tz.ne.0.0d0)
!     &      .or.(gb_ty_xz.ne.0.0d0)
!     &      .or.(gb_ty_yz.ne.0.0d0)
!     &      .or.(gb_ty_zz.ne.0.0d0)) then
!
!          write(*,*) 'DEBUG from misc.f: non zero z-derivative of 
!     &           gb_ty'
!
!         end if
!
!         if    ((gb_tz_n(i,j,k).ne.0.0d0)
!     &      .or.(gb_tz_z.ne.0.0d0)
!     &      .or.(gb_tz_tz.ne.0.0d0)
!     &      .or.(gb_tz_xz.ne.0.0d0)
!     &      .or.(gb_tz_yz.ne.0.0d0)
!     &      .or.(gb_tz_zz.ne.0.0d0)) then
!
!          write(*,*) 'DEBUG from misc.f: non zero
!     &           gb_tz'
!         end if
!
!         if    ((gb_xx_z.ne.0.0d0)
!     &      .or.(gb_xx_tz.ne.0.0d0)
!     &      .or.(gb_xx_xz.ne.0.0d0)
!     &      .or.(gb_xx_yz.ne.0.0d0)
!     &      .or.(gb_xx_zz.ne.0.0d0)) then
!
!          write(*,*) 'DEBUG from misc.f: non zero z-derivative of 
!     &           gb_xx'
!         end if
!
!         if    ((gb_xy_z.ne.0.0d0)
!     &      .or.(gb_xy_tz.ne.0.0d0)
!     &      .or.(gb_xy_xz.ne.0.0d0)
!     &      .or.(gb_xy_yz.ne.0.0d0)
!     &      .or.(gb_xy_zz.ne.0.0d0)) then
!
!          write(*,*) 'DEBUG from misc.f: non zero z-derivative of 
!     &           gb_xy'
!         end if
!
!         if    ((gb_xz_n(i,j,k).ne.0.0d0)
!     &      .or.(gb_xz_z.ne.0.0d0)
!     &      .or.(gb_xz_tz.ne.0.0d0)
!     &      .or.(gb_xz_xz.ne.0.0d0)
!     &      .or.(gb_xz_yz.ne.0.0d0)
!     &      .or.(gb_xz_zz.ne.0.0d0)) then
!
!          write(*,*) 'DEBUG from misc.f: non zero
!     &           gb_xz'
!         end if
!
!         if    ((gb_yy_z.ne.0.0d0)
!     &      .or.(gb_yy_tz.ne.0.0d0)
!     &      .or.(gb_yy_xz.ne.0.0d0)
!     &      .or.(gb_yy_yz.ne.0.0d0)
!     &      .or.(gb_yy_zz.ne.0.0d0)) then
!
!          write(*,*) 'DEBUG from misc.f: non zero z-derivative of 
!     &           gb_yy'
!         end if
!
!         if    ((gb_yz_n(i,j,k).ne.0.0d0)
!     &      .or.(gb_yz_z.ne.0.0d0)
!     &      .or.(gb_yz_tz.ne.0.0d0)
!     &      .or.(gb_yz_xz.ne.0.0d0)
!     &      .or.(gb_yz_yz.ne.0.0d0)
!     &      .or.(gb_yz_zz.ne.0.0d0)) then
!
!          write(*,*) 'DEBUG from misc.f: non zero
!     &           gb_yz'
!         end if
!
!         if    ((gb_zz_z.ne.0.0d0)
!     &      .or.(gb_zz_tz.ne.0.0d0)
!     &      .or.(gb_zz_xz.ne.0.0d0)
!     &      .or.(gb_zz_yz.ne.0.0d0)
!     &      .or.(gb_zz_zz.ne.0.0d0)) then
!
!          write(*,*) 'DEBUG from misc.f: non zero z-derivative of 
!     &           gb_zz'
!         end if
!
!!!TEST
!        do a=1,Nx
!         do b=1,Ny
!          do c=1,Nz
!          testf1(a,b,c)=x(a)**3+2*x(a)
!          testf2(a,b,c)=y(b)**3+2*y(b)*z(c)**7
!          testf3(a,b,c)=x(a)**3+2*y(b)+x(a)**2*y(b)**5*z(c)**5+z(c)
!     &                  +x(a)*z(c)
!          chr(a,b,c)=0
!          chr(2,b,c)=ex
!          chr(4,b,c)=ex
!          chr(7,b,c)=ex
!          chr(11,b,c)=ex
!          chr(a,5,c)=ex
!          chr(a,7,c)=ex
!          chr(a,b,25)=ex
!          chr(a,b,28)=ex
!          end do
!         end do
!        end do
!        
!        do a=1,Nx
!         do b=1,Ny
!          do c=1,Nz
!          if (chr(a,b,c).ne.ex) then
!        call df1_int(testf1,testf1,testf1,testf1_t,testf1_x,
!     &       testf1_y,testf1_z,x,y,z,dt,a,b,c,
!     &       chr,ex,Nx,Ny,Nz,'testf1')
!        call df1_int(testf2,testf2,testf2,testf2_t,testf2_x,
!     &       testf2_y, testf2_z,x,y,z,dt,a,b,c,
!     &       chr,ex,Nx,Ny,Nz,'testf2')
!        call df1_int(testf3,testf3,testf3,testf3_t,testf3_x,
!     &       testf3_y,testf3_z,x,y,z,dt,a,b,c,
!     &       chr,ex,Nx,Ny,Nz,'testf3')
!        call df2_int(testf1,testf1,testf1,testf1_t,testf1_x,
!     &       testf1_y,testf1_z,testf1_tt,testf1_tx,testf1_ty,
!     &       testf1_tz,
!     &       testf1_xx,testf1_xy,
!     &       testf1_xz,
!     &       testf1_yy,
!     &       testf1_yz,
!     &       testf1_zz,
!     &       x,y,z,dt,a,b,c,chr,ex,Nx,Ny,Nz,'testf1')
!        call df2_int(testf2,testf2,testf2,testf2_t,testf2_x,
!     &       testf2_y,testf2_z,testf2_tt,testf2_tx,testf2_ty,
!     &       testf2_tz,
!     &       testf2_xx,testf2_xy,
!     &       testf2_xz,
!     &       testf2_yy,
!     &       testf2_yz,
!     &       testf2_zz,
!     &       x,y,z,dt,a,b,c,chr,ex,Nx,Ny,Nz,'testf2')
!        call df2_int(testf3,testf3,testf3,testf3_t,testf3_x,
!     &       testf3_y,testf3_z,testf3_tt,testf3_tx,testf3_ty,
!     &       testf3_tz,
!     &       testf3_xx,testf3_xy,
!     &       testf3_xz,
!     &       testf3_yy,
!     &       testf3_yz,
!     &       testf3_zz,
!     &       x,y,z,dt,a,b,c,chr,ex,Nx,Ny,Nz,'testf3')
!
!      write(*,*) 'a,b,c,x(a),y(b),z(c),dx,dy,dz='
!     &           ,a,b,c,x(a),y(b),z(c),dx,dy,dz
!      write(*,*) 'testf1_x,testf2_x,testf3_x='
!     &            ,testf1_x,testf2_x,testf3_x
!      write(*,*) 'testf1_z,testf2_z,testf3_z='
!     &            ,testf1_z,testf2_z,testf3_z
!      write(*,*) 'testf1_tz,testf1_xz,testf1_yz,testf1_zz='
!     &            ,testf1_tz,testf1_xz,testf1_yz,testf1_zz
!      write(*,*) 'testf2_tz,testf2_xz,testf2_yz,testf2_zz='
!     &            ,testf2_tz,testf2_xz,testf2_yz,testf2_zz
!      write(*,*) 'testf3_tz,testf3_xz,testf3_yz,testf3_zz='
!     &            ,testf3_tz,testf3_xz,testf3_yz,testf3_zz
!
!         end if
!         end do
!        end do
!       end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!


        ! give values to the metric
        g0_ll(1,1)=g0_tt_ads0+gb_tt0
        g0_ll(1,2)=g0_tx_ads0+gb_tx0
        g0_ll(1,3)=g0_ty_ads0+gb_ty0
        g0_ll(1,4)=g0_tz_ads0+gb_tz0
        g0_ll(2,2)=g0_xx_ads0+gb_xx0
        g0_ll(2,3)=g0_xy_ads0+gb_xy0
        g0_ll(2,4)=g0_xz_ads0+gb_xz0
        g0_ll(3,3)=g0_yy_ads0+gb_yy0
        g0_ll(3,4)=g0_yz_ads0+gb_yz0
        g0_ll(4,4)=g0_zz_ads0+gb_zz0

!CHECKED WITH MATHEMATICA UP TO HERE

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0=',L,i,j,k,x0,y0,z0,rho0
!       write (*,*) 'g0_ll(1,1)=',g0_ll(1,1)
!       write (*,*) 'g0_ll(1,2)=',g0_ll(1,2)
!       write (*,*) 'g0_ll(1,3)=',g0_ll(1,3)
!       write (*,*) 'g0_ll(1,4)=',g0_ll(1,4)
!       write (*,*) 'g0_ll(2,2)=',g0_ll(2,2) 
!       write (*,*) 'g0_ll(2,3)=',g0_ll(2,3) 
!       write (*,*) 'g0_ll(2,4)=',g0_ll(2,4) 
!       write (*,*) 'g0_ll(3,3)=',g0_ll(3,3) 
!       write (*,*) 'g0_ll(3,4)=',g0_ll(3,4) 
!       write (*,*) 'g0_ll(4,4)=',g0_ll(4,4) 

        g0_ll_x(1,1,1)   =0
     &                   +gb_tt_t
        g0_ll_x(1,1,2)   =g0_tt_ads_x
     &                   +gb_tt_x
        g0_ll_x(1,1,3)   =g0_tt_ads_y
     &                   +gb_tt_y
        g0_ll_x(1,1,4)   =g0_tt_ads_z
     &                   +gb_tt_z
        g0_ll_xx(1,1,1,1)=0
     &                   +gb_tt_tt
        g0_ll_xx(1,1,1,2)=0
     &                   +gb_tt_tx
        g0_ll_xx(1,1,1,3)=0
     &                   +gb_tt_ty
        g0_ll_xx(1,1,1,4)=0
     &                   +gb_tt_tz
        g0_ll_xx(1,1,2,2)=g0_tt_ads_xx
     &                   +gb_tt_xx
        g0_ll_xx(1,1,2,3)=g0_tt_ads_xy
     &                   +gb_tt_xy
        g0_ll_xx(1,1,2,4)=g0_tt_ads_xz
     &                   +gb_tt_xz
        g0_ll_xx(1,1,3,3)=g0_tt_ads_yy
     &                   +gb_tt_yy
        g0_ll_xx(1,1,3,4)=g0_tt_ads_yz
     &                   +gb_tt_yz
        g0_ll_xx(1,1,4,4)=g0_tt_ads_zz
     &                   +gb_tt_zz

!
!       write (*,*) 'gb_tt_nm1(i,j,k),gb_tt_n(i,j,k),gb_tt_np1(i,j,k)='
!     &             ,gb_tt_nm1(i,j,k),gb_tt_n(i,j,k),gb_tt_np1(i,j,k)


        g0_ll_x(1,2,1)   =0
     &                   +gb_tx_t
        g0_ll_x(1,2,2)   =0
     &                   +gb_tx_x
        g0_ll_x(1,2,3)   =0
     &                   +gb_tx_y
        g0_ll_x(1,2,4)   =0
     &                   +gb_tx_z
        g0_ll_xx(1,2,1,1)=0
     &                   +gb_tx_tt
        g0_ll_xx(1,2,1,2)=0
     &                   +gb_tx_tx
        g0_ll_xx(1,2,1,3)=0
     &                   +gb_tx_ty
        g0_ll_xx(1,2,1,4)=0
     &                   +gb_tx_tz
        g0_ll_xx(1,2,2,2)=0
     &                   +gb_tx_xx
        g0_ll_xx(1,2,2,3)=0
     &                   +gb_tx_xy
        g0_ll_xx(1,2,2,4)=0
     &                   +gb_tx_xz
        g0_ll_xx(1,2,3,3)=0
     &                   +gb_tx_yy
        g0_ll_xx(1,2,3,4)=0
     &                   +gb_tx_yz
        g0_ll_xx(1,2,4,4)=0
     &                   +gb_tx_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'g0_ll(1,2)=',g0_ll(1,2)
!       write (*,*) 'g0_ll_x(1,2,1)=',g0_ll_x(1,2,1)
!       write (*,*) 'g0_ll_x(1,2,2)=',g0_ll_x(1,2,2)
!       write (*,*) 'g0_ll_x(1,2,3)=',g0_ll_x(1,2,3)
!       write (*,*) 'g0_ll_x(1,2,4)=',g0_ll_x(1,2,4)
!       write (*,*) 'g0_ll_xx(1,2,1,1)=',g0_ll_xx(1,2,1,1)
!       write (*,*) 'g0_ll_xx(1,2,1,2)=',g0_ll_xx(1,2,1,2)
!       write (*,*) 'g0_ll_xx(1,2,1,3)=',g0_ll_xx(1,2,1,3)
!       write (*,*) 'g0_ll_xx(1,2,1,4)=',g0_ll_xx(1,2,1,4)
!       write (*,*) 'g0_ll_xx(1,2,2,2)=',g0_ll_xx(1,2,2,2)
!       write (*,*) 'g0_ll_xx(1,2,2,3)=',g0_ll_xx(1,2,2,3)
!       write (*,*) 'g0_ll_xx(1,2,2,4)=',g0_ll_xx(1,2,2,4)
!       write (*,*) 'g0_ll_xx(1,2,3,3)=',g0_ll_xx(1,2,3,3)
!       write (*,*) 'g0_ll_xx(1,2,3,4)=',g0_ll_xx(1,2,3,4)
!       write (*,*) 'g0_ll_xx(1,2,4,4)=',g0_ll_xx(1,2,4,4)

        g0_ll_x(1,3,1)   =0
     &                   +gb_ty_t
        g0_ll_x(1,3,2)   =0
     &                   +gb_ty_x
        g0_ll_x(1,3,3)   =0
     &                   +gb_ty_y
        g0_ll_x(1,3,4)   =0
     &                   +gb_ty_z
        g0_ll_xx(1,3,1,1)=0
     &                   +gb_ty_tt
        g0_ll_xx(1,3,1,2)=0
     &                   +gb_ty_tx
        g0_ll_xx(1,3,1,3)=0
     &                   +gb_ty_ty
        g0_ll_xx(1,3,1,4)=0
     &                   +gb_ty_tz
        g0_ll_xx(1,3,2,2)=0
     &                   +gb_ty_xx
        g0_ll_xx(1,3,2,3)=0
     &                   +gb_ty_xy
        g0_ll_xx(1,3,2,4)=0
     &                   +gb_ty_xz
        g0_ll_xx(1,3,3,3)=0
     &                   +gb_ty_yy
        g0_ll_xx(1,3,3,4)=0
     &                   +gb_ty_yz
        g0_ll_xx(1,3,4,4)=0
     &                   +gb_ty_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'g0_ll(1,3)=',g0_ll(1,3)
!       write (*,*) 'g0_ll_x(1,3,1)=',g0_ll_x(1,3,1)
!       write (*,*) 'g0_ll_x(1,3,2)=',g0_ll_x(1,3,2)
!       write (*,*) 'g0_ll_x(1,3,3)=',g0_ll_x(1,3,3)
!       write (*,*) 'g0_ll_x(1,3,4)=',g0_ll_x(1,3,4)
!       write (*,*) 'g0_ll_xx(1,3,1,1)=',g0_ll_xx(1,3,1,1)
!       write (*,*) 'g0_ll_xx(1,3,1,2)=',g0_ll_xx(1,3,1,2)
!       write (*,*) 'g0_ll_xx(1,3,1,3)=',g0_ll_xx(1,3,1,3)
!       write (*,*) 'g0_ll_xx(1,3,1,4)=',g0_ll_xx(1,3,1,4)
!       write (*,*) 'g0_ll_xx(1,3,2,2)=',g0_ll_xx(1,3,2,2)
!       write (*,*) 'g0_ll_xx(1,3,2,3)=',g0_ll_xx(1,3,2,3)
!       write (*,*) 'g0_ll_xx(1,3,2,4)=',g0_ll_xx(1,3,2,4)
!       write (*,*) 'g0_ll_xx(1,3,3,3)=',g0_ll_xx(1,3,3,3)
!       write (*,*) 'g0_ll_xx(1,3,3,4)=',g0_ll_xx(1,3,3,4)
!       write (*,*) 'g0_ll_xx(1,3,4,4)=',g0_ll_xx(1,3,4,4)

        g0_ll_x(1,4,1)   =0
     &                   +gb_tz_t
        g0_ll_x(1,4,2)   =0
     &                   +gb_tz_x
        g0_ll_x(1,4,3)   =0
     &                   +gb_tz_y
        g0_ll_x(1,4,4)   =0
     &                   +gb_tz_z
        g0_ll_xx(1,4,1,1)=0
     &                   +gb_tz_tt
        g0_ll_xx(1,4,1,2)=0
     &                   +gb_tz_tx
        g0_ll_xx(1,4,1,3)=0
     &                   +gb_tz_ty
        g0_ll_xx(1,4,1,4)=0
     &                   +gb_tz_tz
        g0_ll_xx(1,4,2,2)=0
     &                   +gb_tz_xx
        g0_ll_xx(1,4,2,3)=0
     &                   +gb_tz_xy
        g0_ll_xx(1,4,2,4)=0
     &                   +gb_tz_xz
        g0_ll_xx(1,4,3,3)=0
     &                   +gb_tz_yy
        g0_ll_xx(1,4,3,4)=0
     &                   +gb_tz_yz
        g0_ll_xx(1,4,4,4)=0
     &                   +gb_tz_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'g0_ll(1,4)=',g0_ll(1,4)
!       write (*,*) 'g0_ll_x(1,4,1)=',g0_ll_x(1,4,1)
!       write (*,*) 'g0_ll_x(1,4,2)=',g0_ll_x(1,4,2)
!       write (*,*) 'g0_ll_x(1,4,3)=',g0_ll_x(1,4,3)
!       write (*,*) 'g0_ll_x(1,4,4)=',g0_ll_x(1,4,4)
!       write (*,*) 'g0_ll_xx(1,4,1,1)=',g0_ll_xx(1,4,1,1)
!       write (*,*) 'g0_ll_xx(1,4,1,2)=',g0_ll_xx(1,4,1,2)
!       write (*,*) 'g0_ll_xx(1,4,1,3)=',g0_ll_xx(1,4,1,3)
!       write (*,*) 'g0_ll_xx(1,4,1,4)=',g0_ll_xx(1,4,1,4)
!       write (*,*) 'g0_ll_xx(1,4,2,2)=',g0_ll_xx(1,4,2,2)
!       write (*,*) 'g0_ll_xx(1,4,2,3)=',g0_ll_xx(1,4,2,3)
!       write (*,*) 'g0_ll_xx(1,4,2,4)=',g0_ll_xx(1,4,2,4)
!       write (*,*) 'g0_ll_xx(1,4,3,3)=',g0_ll_xx(1,4,3,3)
!       write (*,*) 'g0_ll_xx(1,4,3,4)=',g0_ll_xx(1,4,3,4)
!       write (*,*) 'g0_ll_xx(1,4,4,4)=',g0_ll_xx(1,4,4,4)

        g0_ll_x(2,2,1)   =0
     &                   +gb_xx_t
        g0_ll_x(2,2,2)   =g0_xx_ads_x
     &                   +gb_xx_x
        g0_ll_x(2,2,3)   =g0_xx_ads_y
     &                   +gb_xx_y
        g0_ll_x(2,2,4)   =g0_xx_ads_z
     &                   +gb_xx_z
        g0_ll_xx(2,2,1,1)=0
     &                   +gb_xx_tt
        g0_ll_xx(2,2,1,2)=0
     &                   +gb_xx_tx
        g0_ll_xx(2,2,1,3)=0
     &                   +gb_xx_ty
        g0_ll_xx(2,2,1,4)=0
     &                   +gb_xx_tz
        g0_ll_xx(2,2,2,2)=g0_xx_ads_xx
     &                   +gb_xx_xx
        g0_ll_xx(2,2,2,3)=g0_xx_ads_xy
     &                   +gb_xx_xy
        g0_ll_xx(2,2,2,4)=g0_xx_ads_xz
     &                   +gb_xx_xz
        g0_ll_xx(2,2,3,3)=g0_xx_ads_yy
     &                   +gb_xx_yy
        g0_ll_xx(2,2,3,4)=g0_xx_ads_yz
     &                   +gb_xx_yz
        g0_ll_xx(2,2,4,4)=g0_xx_ads_zz
     &                   +gb_xx_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'g0_ll(2,2)=',g0_ll(2,2)
!       write (*,*) 'g0_ll_x(2,2,1)=',g0_ll_x(2,2,1)
!       write (*,*) 'g0_ll_x(2,2,2)=',g0_ll_x(2,2,2)
!       write (*,*) 'g0_ll_x(2,2,3)=',g0_ll_x(2,2,3)
!       write (*,*) 'g0_ll_x(2,2,4)=',g0_ll_x(2,2,4)
!       write (*,*) 'g0_ll_xx(2,2,1,1)=',g0_ll_xx(2,2,1,1)
!       write (*,*) 'g0_ll_xx(2,2,1,2)=',g0_ll_xx(2,2,1,2)
!       write (*,*) 'g0_ll_xx(2,2,1,3)=',g0_ll_xx(2,2,1,3)
!       write (*,*) 'g0_ll_xx(2,2,1,4)=',g0_ll_xx(2,2,1,4)
!       write (*,*) 'g0_ll_xx(2,2,2,2)=',g0_ll_xx(2,2,2,2)
!       write (*,*) 'g0_ll_xx(2,2,2,3)=',g0_ll_xx(2,2,2,3)
!       write (*,*) 'g0_ll_xx(2,2,2,4)=',g0_ll_xx(2,2,2,4)
!       write (*,*) 'g0_ll_xx(2,2,3,3)=',g0_ll_xx(2,2,3,3)
!       write (*,*) 'g0_ll_xx(2,2,3,4)=',g0_ll_xx(2,2,3,4)
!       write (*,*) 'g0_ll_xx(2,2,4,4)=',g0_ll_xx(2,2,4,4)

        g0_ll_x(2,3,1)   =0
     &                   +gb_xy_t
        g0_ll_x(2,3,2)   =g0_xy_ads_x
     &                   +gb_xy_x
        g0_ll_x(2,3,3)   =g0_xy_ads_y
     &                   +gb_xy_y
        g0_ll_x(2,3,4)   =g0_xy_ads_z
     &                   +gb_xy_z
        g0_ll_xx(2,3,1,1)=0
     &                   +gb_xy_tt
        g0_ll_xx(2,3,1,2)=0
     &                   +gb_xy_tx
        g0_ll_xx(2,3,1,3)=0
     &                   +gb_xy_ty
        g0_ll_xx(2,3,1,4)=0
     &                   +gb_xy_tz
        g0_ll_xx(2,3,2,2)=g0_xy_ads_xx
     &                   +gb_xy_xx
        g0_ll_xx(2,3,2,3)=g0_xy_ads_xy
     &                   +gb_xy_xy
        g0_ll_xx(2,3,2,4)=g0_xy_ads_xz
     &                   +gb_xy_xz
        g0_ll_xx(2,3,3,3)=g0_xy_ads_yy
     &                   +gb_xy_yy
        g0_ll_xx(2,3,3,4)=g0_xy_ads_yz
     &                   +gb_xy_yz
        g0_ll_xx(2,3,4,4)=g0_xy_ads_zz
     &                   +gb_xy_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'g0_ll(2,3)=',g0_ll(2,3)
!       write (*,*) 'g0_ll_x(2,3,1)=',g0_ll_x(2,3,1)
!       write (*,*) 'g0_ll_x(2,3,2)=',g0_ll_x(2,3,2)
!       write (*,*) 'g0_ll_x(2,3,3)=',g0_ll_x(2,3,3)
!       write (*,*) 'g0_ll_x(2,3,4)=',g0_ll_x(2,3,4)
!       write (*,*) 'g0_ll_xx(2,3,1,1)=',g0_ll_xx(2,3,1,1)
!       write (*,*) 'g0_ll_xx(2,3,1,2)=',g0_ll_xx(2,3,1,2)
!       write (*,*) 'g0_ll_xx(2,3,1,3)=',g0_ll_xx(2,3,1,3)
!       write (*,*) 'g0_ll_xx(2,3,1,4)=',g0_ll_xx(2,3,1,4)
!       write (*,*) 'g0_ll_xx(2,3,2,2)=',g0_ll_xx(2,3,2,2)
!       write (*,*) 'g0_ll_xx(2,3,2,3)=',g0_ll_xx(2,3,2,3)
!       write (*,*) 'g0_ll_xx(2,3,2,4)=',g0_ll_xx(2,3,2,4)
!       write (*,*) 'g0_ll_xx(2,3,3,3)=',g0_ll_xx(2,3,3,3)
!       write (*,*) 'g0_ll_xx(2,3,3,4)=',g0_ll_xx(2,3,3,4)
!       write (*,*) 'g0_ll_xx(2,3,4,4)=',g0_ll_xx(2,3,4,4)

        g0_ll_x(2,4,1)   =0
     &                   +gb_xz_t
        g0_ll_x(2,4,2)   =g0_xz_ads_x
     &                   +gb_xz_x
        g0_ll_x(2,4,3)   =g0_xz_ads_y
     &                   +gb_xz_y
        g0_ll_x(2,4,4)   =g0_xz_ads_z
     &                   +gb_xz_z
        g0_ll_xx(2,4,1,1)=0
     &                   +gb_xz_tt
        g0_ll_xx(2,4,1,2)=0
     &                   +gb_xz_tx
        g0_ll_xx(2,4,1,3)=0
     &                   +gb_xz_ty
        g0_ll_xx(2,4,1,4)=0
     &                   +gb_xz_tz
        g0_ll_xx(2,4,2,2)=g0_xz_ads_xx
     &                   +gb_xz_xx
        g0_ll_xx(2,4,2,3)=g0_xz_ads_xy
     &                   +gb_xz_xy
        g0_ll_xx(2,4,2,4)=g0_xz_ads_xz
     &                   +gb_xz_xz
        g0_ll_xx(2,4,3,3)=g0_xz_ads_yy
     &                   +gb_xz_yy
        g0_ll_xx(2,4,3,4)=g0_xz_ads_yz
     &                   +gb_xz_yz
        g0_ll_xx(2,4,4,4)=g0_xz_ads_zz
     &                   +gb_xz_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'g0_ll(2,4)=',g0_ll(2,4)
!       write (*,*) 'g0_ll_x(2,4,1)=',g0_ll_x(2,4,1)
!       write (*,*) 'g0_ll_x(2,4,2)=',g0_ll_x(2,4,2)
!       write (*,*) 'g0_ll_x(2,4,3)=',g0_ll_x(2,4,3)
!       write (*,*) 'g0_ll_x(2,4,4)=',g0_ll_x(2,4,4)
!       write (*,*) 'g0_ll_xx(2,4,1,1)=',g0_ll_xx(2,4,1,1)
!       write (*,*) 'g0_ll_xx(2,4,1,2)=',g0_ll_xx(2,4,1,2)
!       write (*,*) 'g0_ll_xx(2,4,1,3)=',g0_ll_xx(2,4,1,3)
!       write (*,*) 'g0_ll_xx(2,4,1,4)=',g0_ll_xx(2,4,1,4)
!       write (*,*) 'g0_ll_xx(2,4,2,2)=',g0_ll_xx(2,4,2,2)
!       write (*,*) 'g0_ll_xx(2,4,2,3)=',g0_ll_xx(2,4,2,3)
!       write (*,*) 'g0_ll_xx(2,4,2,4)=',g0_ll_xx(2,4,2,4)
!       write (*,*) 'g0_ll_xx(2,4,3,3)=',g0_ll_xx(2,4,3,3)
!       write (*,*) 'g0_ll_xx(2,4,3,4)=',g0_ll_xx(2,4,3,4)
!       write (*,*) 'g0_ll_xx(2,4,4,4)=',g0_ll_xx(2,4,4,4)

        g0_ll_x(3,3,1)   =0
     &                   +gb_yy_t
        g0_ll_x(3,3,2)   =g0_yy_ads_x
     &                   +gb_yy_x
        g0_ll_x(3,3,3)   =g0_yy_ads_y
     &                   +gb_yy_y
        g0_ll_x(3,3,4)   =g0_yy_ads_z
     &                   +gb_yy_z
        g0_ll_xx(3,3,1,1)=0
     &                   +gb_yy_tt
        g0_ll_xx(3,3,1,2)=0
     &                   +gb_yy_tx
        g0_ll_xx(3,3,1,3)=0
     &                   +gb_yy_ty
        g0_ll_xx(3,3,1,4)=0
     &                   +gb_yy_tz
        g0_ll_xx(3,3,2,2)=g0_yy_ads_xx
     &                   +gb_yy_xx
        g0_ll_xx(3,3,2,3)=g0_yy_ads_xy
     &                   +gb_yy_xy
        g0_ll_xx(3,3,2,4)=g0_yy_ads_xz
     &                   +gb_yy_xz
        g0_ll_xx(3,3,3,3)=g0_yy_ads_yy
     &                   +gb_yy_yy
        g0_ll_xx(3,3,3,4)=g0_yy_ads_yz
     &                   +gb_yy_yz
        g0_ll_xx(3,3,4,4)=g0_yy_ads_zz
     &                   +gb_yy_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'g0_ll(3,3)=',g0_ll(3,3)
!       write (*,*) 'g0_ll_x(3,3,1)=',g0_ll_x(3,3,1)
!       write (*,*) 'g0_ll_x(3,3,2)=',g0_ll_x(3,3,2)
!       write (*,*) 'g0_ll_x(3,3,3)=',g0_ll_x(3,3,3)
!       write (*,*) 'g0_ll_x(3,3,4)=',g0_ll_x(3,3,4)
!       write (*,*) 'g0_ll_xx(3,3,1,1)=',g0_ll_xx(3,3,1,1)
!       write (*,*) 'g0_ll_xx(3,3,1,2)=',g0_ll_xx(3,3,1,2)
!       write (*,*) 'g0_ll_xx(3,3,1,3)=',g0_ll_xx(3,3,1,3)
!       write (*,*) 'g0_ll_xx(3,3,1,4)=',g0_ll_xx(3,3,1,4)
!       write (*,*) 'g0_ll_xx(3,3,2,2)=',g0_ll_xx(3,3,2,2)
!       write (*,*) 'g0_ll_xx(3,3,2,3)=',g0_ll_xx(3,3,2,3)
!       write (*,*) 'g0_ll_xx(3,3,2,4)=',g0_ll_xx(3,3,2,4)
!       write (*,*) 'g0_ll_xx(3,3,3,3)=',g0_ll_xx(3,3,3,3)
!       write (*,*) 'g0_ll_xx(3,3,3,4)=',g0_ll_xx(3,3,3,4)
!       write (*,*) 'g0_ll_xx(3,3,4,4)=',g0_ll_xx(3,3,4,4)

        g0_ll_x(3,4,1)   =0
     &                   +gb_yz_t
        g0_ll_x(3,4,2)   =g0_yz_ads_x
     &                   +gb_yz_x
        g0_ll_x(3,4,3)   =g0_yz_ads_y
     &                   +gb_yz_y
        g0_ll_x(3,4,4)   =g0_yz_ads_z
     &                   +gb_yz_z
        g0_ll_xx(3,4,1,1)=0
     &                   +gb_yz_tt
        g0_ll_xx(3,4,1,2)=0
     &                   +gb_yz_tx
        g0_ll_xx(3,4,1,3)=0
     &                   +gb_yz_ty
        g0_ll_xx(3,4,1,4)=0
     &                   +gb_yz_tz
        g0_ll_xx(3,4,2,2)=g0_yz_ads_xx
     &                   +gb_yz_xx
        g0_ll_xx(3,4,2,3)=g0_yz_ads_xy
     &                   +gb_yz_xy
        g0_ll_xx(3,4,2,4)=g0_yz_ads_xz
     &                   +gb_yz_xz
        g0_ll_xx(3,4,3,3)=g0_yz_ads_yy
     &                   +gb_yz_yy
        g0_ll_xx(3,4,3,4)=g0_yz_ads_yz
     &                   +gb_yz_yz
        g0_ll_xx(3,4,4,4)=g0_yz_ads_zz
     &                   +gb_yz_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'g0_ll(3,4)=',g0_ll(3,4)
!       write (*,*) 'g0_ll_x(3,4,1)=',g0_ll_x(3,4,1)
!       write (*,*) 'g0_ll_x(3,4,2)=',g0_ll_x(3,4,2)
!       write (*,*) 'g0_ll_x(3,4,3)=',g0_ll_x(3,4,3)
!       write (*,*) 'g0_ll_x(3,4,4)=',g0_ll_x(3,4,4)
!       write (*,*) 'g0_ll_xx(3,4,1,1)=',g0_ll_xx(3,4,1,1)
!       write (*,*) 'g0_ll_xx(3,4,1,2)=',g0_ll_xx(3,4,1,2)
!       write (*,*) 'g0_ll_xx(3,4,1,3)=',g0_ll_xx(3,4,1,3)
!       write (*,*) 'g0_ll_xx(3,4,1,4)=',g0_ll_xx(3,4,1,4)
!       write (*,*) 'g0_ll_xx(3,4,2,2)=',g0_ll_xx(3,4,2,2)
!       write (*,*) 'g0_ll_xx(3,4,2,3)=',g0_ll_xx(3,4,2,3)
!       write (*,*) 'g0_ll_xx(3,4,2,4)=',g0_ll_xx(3,4,2,4)
!       write (*,*) 'g0_ll_xx(3,4,3,3)=',g0_ll_xx(3,4,3,3)
!       write (*,*) 'g0_ll_xx(3,4,3,4)=',g0_ll_xx(3,4,3,4)
!       write (*,*) 'g0_ll_xx(3,4,4,4)=',g0_ll_xx(3,4,4,4)

        g0_ll_x(4,4,1)   =0
     &                   +gb_zz_t
        g0_ll_x(4,4,2)   =g0_zz_ads_x
     &                   +gb_zz_x
        g0_ll_x(4,4,3)   =g0_zz_ads_y
     &                   +gb_zz_y
        g0_ll_x(4,4,4)   =g0_zz_ads_z
     &                   +gb_zz_z
        g0_ll_xx(4,4,1,1)=0
     &                   +gb_zz_tt
        g0_ll_xx(4,4,1,2)=0
     &                   +gb_zz_tx
        g0_ll_xx(4,4,1,3)=0
     &                   +gb_zz_ty
        g0_ll_xx(4,4,1,4)=0
     &                   +gb_zz_tz
        g0_ll_xx(4,4,2,2)=g0_zz_ads_xx
     &                   +gb_zz_xx
        g0_ll_xx(4,4,2,3)=g0_zz_ads_xy
     &                   +gb_zz_xy
        g0_ll_xx(4,4,2,4)=g0_zz_ads_xz
     &                   +gb_zz_xz
        g0_ll_xx(4,4,3,3)=g0_zz_ads_yy
     &                   +gb_zz_yy
        g0_ll_xx(4,4,3,4)=g0_zz_ads_yz
     &                   +gb_zz_yz
        g0_ll_xx(4,4,4,4)=g0_zz_ads_zz
     &                   +gb_zz_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'g0_ll(4,4)=',g0_ll(4,4)
!       write (*,*) 'g0_ll_x(4,4,1)=',g0_ll_x(4,4,1)
!       write (*,*) 'g0_ll_x(4,4,2)=',g0_ll_x(4,4,2)
!       write (*,*) 'g0_ll_x(4,4,3)=',g0_ll_x(4,4,3)
!       write (*,*) 'g0_ll_x(4,4,4)=',g0_ll_x(4,4,4)
!       write (*,*) 'g0_ll_xx(4,4,1,1)=',g0_ll_xx(4,4,1,1)
!       write (*,*) 'g0_ll_xx(4,4,1,2)=',g0_ll_xx(4,4,1,2)
!       write (*,*) 'g0_ll_xx(4,4,1,3)=',g0_ll_xx(4,4,1,3)
!       write (*,*) 'g0_ll_xx(4,4,1,4)=',g0_ll_xx(4,4,1,4)
!       write (*,*) 'g0_ll_xx(4,4,2,2)=',g0_ll_xx(4,4,2,2)
!       write (*,*) 'g0_ll_xx(4,4,2,3)=',g0_ll_xx(4,4,2,3)
!       write (*,*) 'g0_ll_xx(4,4,2,4)=',g0_ll_xx(4,4,2,4)
!       write (*,*) 'g0_ll_xx(4,4,3,3)=',g0_ll_xx(4,4,3,3)
!       write (*,*) 'g0_ll_xx(4,4,3,4)=',g0_ll_xx(4,4,3,4)
!       write (*,*) 'g0_ll_xx(4,4,4,4)=',g0_ll_xx(4,4,4,4)


        ! give values to the metric inverse
        call calc_g0uu(g0_ll(1,1),g0_ll(1,2),g0_ll(1,3),g0_ll(1,4),
     &         g0_ll(2,2),g0_ll(2,3),g0_ll(2,4),
     &         g0_ll(3,3),g0_ll(3,4),g0_ll(4,4),g0_uu)

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'g0_uu(1,1)=',g0_uu(1,1)
!       write (*,*) 'g0_uu(1,2)=',g0_uu(1,2)
!       write (*,*) 'g0_uu(1,3)=',g0_uu(1,3)
!       write (*,*) 'g0_uu(1,4)=',g0_uu(1,4)
!       write (*,*) 'g0_uu(2,2)=',g0_uu(2,2)
!       write (*,*) 'g0_uu(2,3)=',g0_uu(2,3)
!       write (*,*) 'g0_uu(2,4)=',g0_uu(2,4)
!       write (*,*) 'g0_uu(3,3)=',g0_uu(3,3)
!       write (*,*) 'g0_uu(3,4)=',g0_uu(3,4)
!       write (*,*) 'g0_uu(4,4)=',g0_uu(4,4)

!CHECKED WITH MATHEMATICA UP TO HERE

        do a=1,3
          do b=a+1,4
            g0_ll(b,a)=g0_ll(a,b)
            g0_uu(b,a)=g0_uu(a,b) 
            do c=1,4
              g0_ll_x(b,a,c)=g0_ll_x(a,b,c)
            end do
          end do
        end do

        do a=1,4
          do b=1,4
            do c=1,4
              g0_uu_x(a,b,c)=0
              do d=1,4
                do e=1,4
                  g0_uu_x(a,b,c)=g0_uu_x(a,b,c)
     &                          -g0_ll_x(d,e,c)
     &                           *g0_uu(a,d)*g0_uu(b,e)
                end do
     &  
              end do
            end do
          end do
        end do

        do a=1,4
          do b=1,4
            do c=1,4
              do d=1,4
                g0_ll_xx(a,b,c,d)=
     &             g0_ll_xx(min(a,b),max(a,b),min(c,d),max(c,d))
              end do
            end do
          end do
        end do

        ! give values to the Christoffel symbols
        do a=1,4
          do b=1,4
            do c=1,4
              gamma_ull(a,b,c)=0
              do d=1,4
                gamma_ull(a,b,c)=gamma_ull(a,b,c)
     &                          +0.5d0*g0_uu(a,d)
     &                                *(g0_ll_x(c,d,b)
     &                                 -g0_ll_x(b,c,d)
     &                                 +g0_ll_x(d,b,c))
              end do
            end do
          end do
        end do

        ! calculate Christoffel symbol derivatives at point i,j
        !(gamma^a_bc,e = 1/2 g^ad_,e(g_bd,c  + g_cd,b  - g_bc,d)
        !              +   1/2 g^ad(g_bd,ce + g_cd,be - g_bc,de))
        do a=1,4
          do b=1,4
            do c=1,4
              do e=1,4
                gamma_ull_x(a,b,c,e)=0
                do d=1,4
                  gamma_ull_x(a,b,c,e)=gamma_ull_x(a,b,c,e)
     &              +0.5d0*g0_uu_x(a,d,e)*(g0_ll_x(b,d,c)+
     &                     g0_ll_x(c,d,b)-g0_ll_x(b,c,d))
     &              +0.5d0*g0_uu(a,d)*(g0_ll_xx(b,d,c,e)+
     &                     g0_ll_xx(c,d,b,e)-g0_ll_xx(b,c,d,e))
                end do
              end do
            end do
          end do
        end do

        ! give values to the AdS metric
        gads_ll(1,1)=g0_tt_ads0
        gads_ll(1,2)=g0_tx_ads0
        gads_ll(1,3)=g0_ty_ads0
        gads_ll(1,4)=g0_tz_ads0
        gads_ll(2,2)=g0_xx_ads0
        gads_ll(2,3)=g0_xy_ads0
        gads_ll(2,4)=g0_xz_ads0
        gads_ll(3,3)=g0_yy_ads0
        gads_ll(3,4)=g0_yz_ads0
        gads_ll(4,4)=g0_zz_ads0

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'gads_ll(1,1)=',gads_ll(1,1)
!       write (*,*) 'gads_ll(1,2)=',gads_ll(1,2)
!       write (*,*) 'gads_ll(1,3)=',gads_ll(1,3)
!       write (*,*) 'gads_ll(1,4)=',gads_ll(1,4)
!       write (*,*) 'gads_ll(2,2)=',gads_ll(2,2)
!       write (*,*) 'gads_ll(2,3)=',gads_ll(2,3)
!       write (*,*) 'gads_ll(2,4)=',gads_ll(2,4)
!       write (*,*) 'gads_ll(3,3)=',gads_ll(3,3)
!       write (*,*) 'gads_ll(3,4)=',gads_ll(3,4)
!       write (*,*) 'gads_ll(4,4)=',gads_ll(4,4)

        gads_uu(1,1)=g0u_tt_ads0
        gads_uu(1,2)=g0u_tx_ads0
        gads_uu(1,3)=g0u_ty_ads0
        gads_uu(1,4)=g0u_tz_ads0
        gads_uu(2,2)=g0u_xx_ads0
        gads_uu(2,3)=g0u_xy_ads0
        gads_uu(2,4)=g0u_xz_ads0
        gads_uu(3,3)=g0u_yy_ads0
        gads_uu(3,4)=g0u_yz_ads0
        gads_uu(4,4)=g0u_zz_ads0


!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'gads_uu(1,1)=',gads_uu(1,1)
!       write (*,*) 'gads_uu(1,2)=',gads_uu(1,2)
!       write (*,*) 'gads_uu(1,3)=',gads_uu(1,3)
!       write (*,*) 'gads_uu(1,4)=',gads_uu(1,4)
!       write (*,*) 'gads_uu(2,2)=',gads_uu(2,2)
!       write (*,*) 'gads_uu(2,3)=',gads_uu(2,3)
!       write (*,*) 'gads_uu(2,4)=',gads_uu(2,4)
!       write (*,*) 'gads_uu(3,3)=',gads_uu(3,3)
!       write (*,*) 'gads_uu(3,4)=',gads_uu(3,4)
!       write (*,*) 'gads_uu(4,4)=',gads_uu(4,4)

        gads_ll_x(1,1,2)   =g0_tt_ads_x
        gads_ll_x(1,1,3)   =g0_tt_ads_y
        gads_ll_x(1,1,4)   =g0_tt_ads_z
        gads_ll_xx(1,1,2,2)=g0_tt_ads_xx
        gads_ll_xx(1,1,2,3)=g0_tt_ads_xy
        gads_ll_xx(1,1,2,4)=g0_tt_ads_xz
        gads_ll_xx(1,1,3,3)=g0_tt_ads_yy
        gads_ll_xx(1,1,3,4)=g0_tt_ads_yz
        gads_ll_xx(1,1,4,4)=g0_tt_ads_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'gads_ll_x(1,1,2)=',gads_ll_x(1,1,2)
!       write (*,*) 'gads_ll_x(1,1,3)=',gads_ll_x(1,1,3)
!       write (*,*) 'gads_ll_x(1,1,4)=',gads_ll_x(1,1,4)
!       write (*,*) 'gads_ll_xx(1,1,2,2)=',gads_ll_xx(1,1,2,2)
!       write (*,*) 'gads_ll_xx(1,1,2,3)=',gads_ll_xx(1,1,2,3)
!       write (*,*) 'gads_ll_xx(1,1,2,4)=',gads_ll_xx(1,1,2,4)
!       write (*,*) 'gads_ll_xx(1,1,3,3)=',gads_ll_xx(1,1,3,3)
!       write (*,*) 'gads_ll_xx(1,1,3,4)=',gads_ll_xx(1,1,3,4)
!       write (*,*) 'gads_ll_xx(1,1,4,4)=',gads_ll_xx(1,1,4,4)

        gads_ll_x(1,2,2)   =g0_tx_ads_x
        gads_ll_x(1,2,3)   =g0_tx_ads_y
        gads_ll_x(1,2,4)   =g0_tx_ads_z
        gads_ll_xx(1,2,2,2)=g0_tx_ads_xx
        gads_ll_xx(1,2,2,3)=g0_tx_ads_xy
        gads_ll_xx(1,2,2,4)=g0_tx_ads_xz
        gads_ll_xx(1,2,3,3)=g0_tx_ads_yy
        gads_ll_xx(1,2,3,4)=g0_tx_ads_yz
        gads_ll_xx(1,2,4,4)=g0_tx_ads_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'gads_ll_x(1,2,2)=',gads_ll_x(1,2,2)
!       write (*,*) 'gads_ll_x(1,2,3)=',gads_ll_x(1,2,3)
!       write (*,*) 'gads_ll_x(1,2,4)=',gads_ll_x(1,2,4)
!       write (*,*) 'gads_ll_xx(1,2,2,2)=',gads_ll_xx(1,2,2,2)
!       write (*,*) 'gads_ll_xx(1,2,2,3)=',gads_ll_xx(1,2,2,3)
!       write (*,*) 'gads_ll_xx(1,2,2,4)=',gads_ll_xx(1,2,2,4)
!       write (*,*) 'gads_ll_xx(1,2,3,3)=',gads_ll_xx(1,2,3,3)
!       write (*,*) 'gads_ll_xx(1,2,3,4)=',gads_ll_xx(1,2,3,4)
!       write (*,*) 'gads_ll_xx(1,2,4,4)=',gads_ll_xx(1,2,4,4)

        gads_ll_x(1,3,2)   =g0_ty_ads_x
        gads_ll_x(1,3,3)   =g0_ty_ads_y
        gads_ll_x(1,3,4)   =g0_ty_ads_z
        gads_ll_xx(1,3,2,2)=g0_ty_ads_xx
        gads_ll_xx(1,3,2,3)=g0_ty_ads_xy
        gads_ll_xx(1,3,2,4)=g0_ty_ads_xz
        gads_ll_xx(1,3,3,3)=g0_ty_ads_yy
        gads_ll_xx(1,3,3,4)=g0_ty_ads_yz
        gads_ll_xx(1,3,4,4)=g0_ty_ads_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'gads_ll_x(1,3,2)=',gads_ll_x(1,3,2)
!       write (*,*) 'gads_ll_x(1,3,3)=',gads_ll_x(1,3,3)
!       write (*,*) 'gads_ll_x(1,3,4)=',gads_ll_x(1,3,4)
!       write (*,*) 'gads_ll_xx(1,3,2,2)=',gads_ll_xx(1,3,2,2)
!       write (*,*) 'gads_ll_xx(1,3,2,3)=',gads_ll_xx(1,3,2,3)
!       write (*,*) 'gads_ll_xx(1,3,2,4)=',gads_ll_xx(1,3,2,4)
!       write (*,*) 'gads_ll_xx(1,3,3,3)=',gads_ll_xx(1,3,3,3)
!       write (*,*) 'gads_ll_xx(1,3,3,4)=',gads_ll_xx(1,3,3,4)
!       write (*,*) 'gads_ll_xx(1,3,4,4)=',gads_ll_xx(1,3,4,4)


        gads_ll_x(1,4,2)   =g0_tz_ads_x
        gads_ll_x(1,4,3)   =g0_tz_ads_y
        gads_ll_x(1,4,4)   =g0_tz_ads_z
        gads_ll_xx(1,4,2,2)=g0_tz_ads_xx
        gads_ll_xx(1,4,2,3)=g0_tz_ads_xy
        gads_ll_xx(1,4,2,4)=g0_tz_ads_xz
        gads_ll_xx(1,4,3,3)=g0_tz_ads_yy
        gads_ll_xx(1,4,3,4)=g0_tz_ads_yz
        gads_ll_xx(1,4,4,4)=g0_tz_ads_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'gads_ll_x(1,4,2)=',gads_ll_x(1,4,2)
!       write (*,*) 'gads_ll_x(1,4,3)=',gads_ll_x(1,4,3)
!       write (*,*) 'gads_ll_x(1,4,4)=',gads_ll_x(1,4,4)
!       write (*,*) 'gads_ll_xx(1,4,2,2)=',gads_ll_xx(1,4,2,2)
!       write (*,*) 'gads_ll_xx(1,4,2,3)=',gads_ll_xx(1,4,2,3)
!       write (*,*) 'gads_ll_xx(1,4,2,4)=',gads_ll_xx(1,4,2,4)
!       write (*,*) 'gads_ll_xx(1,4,3,3)=',gads_ll_xx(1,4,3,3)
!       write (*,*) 'gads_ll_xx(1,4,3,4)=',gads_ll_xx(1,4,3,4)
!       write (*,*) 'gads_ll_xx(1,4,4,4)=',gads_ll_xx(1,4,4,4)


        gads_ll_x(2,2,2)   =g0_xx_ads_x
        gads_ll_x(2,2,3)   =g0_xx_ads_y
        gads_ll_x(2,2,4)   =g0_xx_ads_z
        gads_ll_xx(2,2,2,2)=g0_xx_ads_xx
        gads_ll_xx(2,2,2,3)=g0_xx_ads_xy
        gads_ll_xx(2,2,2,4)=g0_xx_ads_xz
        gads_ll_xx(2,2,3,3)=g0_xx_ads_yy
        gads_ll_xx(2,2,3,4)=g0_xx_ads_yz
        gads_ll_xx(2,2,4,4)=g0_xx_ads_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'gads_ll_x(2,2,2)=',gads_ll_x(2,2,2)
!       write (*,*) 'gads_ll_x(2,2,3)=',gads_ll_x(2,2,3)
!       write (*,*) 'gads_ll_x(2,2,4)=',gads_ll_x(2,2,4)
!       write (*,*) 'gads_ll_xx(2,2,2,2)=',gads_ll_xx(2,2,2,2)
!       write (*,*) 'gads_ll_xx(2,2,2,3)=',gads_ll_xx(2,2,2,3)
!       write (*,*) 'gads_ll_xx(2,2,2,4)=',gads_ll_xx(2,2,2,4)
!       write (*,*) 'gads_ll_xx(2,2,3,3)=',gads_ll_xx(2,2,3,3)
!       write (*,*) 'gads_ll_xx(2,2,3,4)=',gads_ll_xx(2,2,3,4)
!       write (*,*) 'gads_ll_xx(2,2,4,4)=',gads_ll_xx(2,2,4,4)


        gads_ll_x(2,3,2)   =g0_xy_ads_x
        gads_ll_x(2,3,3)   =g0_xy_ads_y
        gads_ll_x(2,3,4)   =g0_xy_ads_z
        gads_ll_xx(2,3,2,2)=g0_xy_ads_xx
        gads_ll_xx(2,3,2,3)=g0_xy_ads_xy
        gads_ll_xx(2,3,2,4)=g0_xy_ads_xz
        gads_ll_xx(2,3,3,3)=g0_xy_ads_yy
        gads_ll_xx(2,3,3,4)=g0_xy_ads_yz
        gads_ll_xx(2,3,4,4)=g0_xy_ads_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'gads_ll_x(2,3,2)=',gads_ll_x(2,3,2)
!       write (*,*) 'gads_ll_x(2,3,3)=',gads_ll_x(2,3,3)
!       write (*,*) 'gads_ll_x(2,3,4)=',gads_ll_x(2,3,4)
!       write (*,*) 'gads_ll_xx(2,3,2,2)=',gads_ll_xx(2,3,2,2)
!       write (*,*) 'gads_ll_xx(2,3,2,3)=',gads_ll_xx(2,3,2,3)
!       write (*,*) 'gads_ll_xx(2,3,2,4)=',gads_ll_xx(2,3,2,4)
!       write (*,*) 'gads_ll_xx(2,3,3,3)=',gads_ll_xx(2,3,3,3)
!       write (*,*) 'gads_ll_xx(2,3,3,4)=',gads_ll_xx(2,3,3,4)
!       write (*,*) 'gads_ll_xx(2,3,4,4)=',gads_ll_xx(2,3,4,4)

        gads_ll_x(2,4,2)   =g0_xz_ads_x
        gads_ll_x(2,4,3)   =g0_xz_ads_y
        gads_ll_x(2,4,4)   =g0_xz_ads_z
        gads_ll_xx(2,4,2,2)=g0_xz_ads_xx
        gads_ll_xx(2,4,2,3)=g0_xz_ads_xy
        gads_ll_xx(2,4,2,4)=g0_xz_ads_xz
        gads_ll_xx(2,4,3,3)=g0_xz_ads_yy
        gads_ll_xx(2,4,3,4)=g0_xz_ads_yz
        gads_ll_xx(2,4,4,4)=g0_xz_ads_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'gads_ll_x(2,4,2)=',gads_ll_x(2,4,2)
!       write (*,*) 'gads_ll_x(2,4,3)=',gads_ll_x(2,4,3)
!       write (*,*) 'gads_ll_x(2,4,4)=',gads_ll_x(2,4,4)
!       write (*,*) 'gads_ll_xx(2,4,2,2)=',gads_ll_xx(2,4,2,2)
!       write (*,*) 'gads_ll_xx(2,4,2,3)=',gads_ll_xx(2,4,2,3)
!       write (*,*) 'gads_ll_xx(2,4,2,4)=',gads_ll_xx(2,4,2,4)
!       write (*,*) 'gads_ll_xx(2,4,3,3)=',gads_ll_xx(2,4,3,3)
!       write (*,*) 'gads_ll_xx(2,4,3,4)=',gads_ll_xx(2,4,3,4)
!       write (*,*) 'gads_ll_xx(2,4,4,4)=',gads_ll_xx(2,4,4,4)

        gads_ll_x(3,3,2)   =g0_yy_ads_x
        gads_ll_x(3,3,3)   =g0_yy_ads_y
        gads_ll_x(3,3,4)   =g0_yy_ads_z
        gads_ll_xx(3,3,2,2)=g0_yy_ads_xx
        gads_ll_xx(3,3,2,3)=g0_yy_ads_xy
        gads_ll_xx(3,3,2,4)=g0_yy_ads_xz
        gads_ll_xx(3,3,3,3)=g0_yy_ads_yy
        gads_ll_xx(3,3,3,4)=g0_yy_ads_yz
        gads_ll_xx(3,3,4,4)=g0_yy_ads_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'gads_ll_x(3,3,2)=',gads_ll_x(3,3,2)
!       write (*,*) 'gads_ll_x(3,3,3)=',gads_ll_x(3,3,3)
!       write (*,*) 'gads_ll_x(3,3,4)=',gads_ll_x(3,3,4)
!       write (*,*) 'gads_ll_xx(3,3,2,2)=',gads_ll_xx(3,3,2,2)
!       write (*,*) 'gads_ll_xx(3,3,2,3)=',gads_ll_xx(3,3,2,3)
!       write (*,*) 'gads_ll_xx(3,3,2,4)=',gads_ll_xx(3,3,2,4)
!       write (*,*) 'gads_ll_xx(3,3,3,3)=',gads_ll_xx(3,3,3,3)
!       write (*,*) 'gads_ll_xx(3,3,3,4)=',gads_ll_xx(3,3,3,4)
!       write (*,*) 'gads_ll_xx(3,3,4,4)=',gads_ll_xx(3,3,4,4)

        gads_ll_x(3,4,2)   =g0_yz_ads_x
        gads_ll_x(3,4,3)   =g0_yz_ads_y
        gads_ll_x(3,4,4)   =g0_yz_ads_z
        gads_ll_xx(3,4,2,2)=g0_yz_ads_xx
        gads_ll_xx(3,4,2,3)=g0_yz_ads_xy
        gads_ll_xx(3,4,2,4)=g0_yz_ads_xz
        gads_ll_xx(3,4,3,3)=g0_yz_ads_yy
        gads_ll_xx(3,4,3,4)=g0_yz_ads_yz
        gads_ll_xx(3,4,4,4)=g0_yz_ads_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'gads_ll_x(3,4,2)=',gads_ll_x(3,4,2)
!       write (*,*) 'gads_ll_x(3,4,3)=',gads_ll_x(3,4,3)
!       write (*,*) 'gads_ll_x(3,4,4)=',gads_ll_x(3,4,4)
!       write (*,*) 'gads_ll_xx(3,4,2,2)=',gads_ll_xx(3,4,2,2)
!       write (*,*) 'gads_ll_xx(3,4,2,3)=',gads_ll_xx(3,4,2,3)
!       write (*,*) 'gads_ll_xx(3,4,2,4)=',gads_ll_xx(3,4,2,4)
!       write (*,*) 'gads_ll_xx(3,4,3,3)=',gads_ll_xx(3,4,3,3)
!       write (*,*) 'gads_ll_xx(3,4,3,4)=',gads_ll_xx(3,4,3,4)
!       write (*,*) 'gads_ll_xx(3,4,4,4)=',gads_ll_xx(3,4,4,4)

        gads_ll_x(4,4,2)   =g0_zz_ads_x
        gads_ll_x(4,4,3)   =g0_zz_ads_y
        gads_ll_x(4,4,4)   =g0_zz_ads_z
        gads_ll_xx(4,4,2,2)=g0_zz_ads_xx
        gads_ll_xx(4,4,2,3)=g0_zz_ads_xy
        gads_ll_xx(4,4,2,4)=g0_zz_ads_xz
        gads_ll_xx(4,4,3,3)=g0_zz_ads_yy
        gads_ll_xx(4,4,3,4)=g0_zz_ads_yz
        gads_ll_xx(4,4,4,4)=g0_zz_ads_zz

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'gads_ll_x(4,4,2)=',gads_ll_x(4,4,2)
!       write (*,*) 'gads_ll_x(4,4,3)=',gads_ll_x(4,4,3)
!       write (*,*) 'gads_ll_x(4,4,4)=',gads_ll_x(4,4,4)
!       write (*,*) 'gads_ll_xx(4,4,2,2)=',gads_ll_xx(4,4,2,2)
!       write (*,*) 'gads_ll_xx(4,4,2,3)=',gads_ll_xx(4,4,2,3)
!       write (*,*) 'gads_ll_xx(4,4,2,4)=',gads_ll_xx(4,4,2,4)
!       write (*,*) 'gads_ll_xx(4,4,3,3)=',gads_ll_xx(4,4,3,3)
!       write (*,*) 'gads_ll_xx(4,4,3,4)=',gads_ll_xx(4,4,3,4)
!       write (*,*) 'gads_ll_xx(4,4,4,4)=',gads_ll_xx(4,4,4,4)
!
        do a=1,3
          do b=a+1,4
            gads_ll(b,a)=gads_ll(a,b)
            gads_uu(b,a)=gads_uu(a,b)
            do c=1,4
              gads_ll_x(b,a,c)=gads_ll_x(a,b,c)
            end do
          end do
        end do

        do a=1,4
          do b=1,4
            do c=1,4
              gads_uu_x(a,b,c)=
     &              -gads_ll_x(1,1,c)*gads_uu(a,1)*gads_uu(b,1)
     &              -gads_ll_x(1,2,c)*(gads_uu(a,1)*gads_uu(b,2)
     &                               +gads_uu(a,2)*gads_uu(b,1))
     &              -gads_ll_x(1,3,c)*(gads_uu(a,1)*gads_uu(b,3)
     &                               +gads_uu(a,3)*gads_uu(b,1))
     &              -gads_ll_x(1,4,c)*(gads_uu(a,1)*gads_uu(b,4)
     &                               +gads_uu(a,4)*gads_uu(b,1))
     &              -gads_ll_x(2,2,c)*gads_uu(a,2)*gads_uu(b,2)
     &              -gads_ll_x(2,3,c)*(gads_uu(a,2)*gads_uu(b,3)
     &                               +gads_uu(a,3)*gads_uu(b,2))
     &              -gads_ll_x(2,4,c)*(gads_uu(a,2)*gads_uu(b,4)
     &                               +gads_uu(a,4)*gads_uu(b,2))
     &              -gads_ll_x(3,3,c)*gads_uu(a,3)*gads_uu(b,3)
     &              -gads_ll_x(3,4,c)*(gads_uu(a,3)*gads_uu(b,4)
     &                               +gads_uu(a,4)*gads_uu(b,3))
     &              -gads_ll_x(4,4,c)*gads_uu(a,4)*gads_uu(b,4)
            end do
          end do
        end do


        ! give values to the metric deviation
        h0_ll(1,1)=gb_tt0 
        h0_ll(1,2)=gb_tx0
        h0_ll(1,3)=gb_ty0
        h0_ll(1,4)=gb_tz0
        h0_ll(2,2)=gb_xx0
        h0_ll(2,3)=gb_xy0
        h0_ll(2,4)=gb_xz0
        h0_ll(3,3)=gb_yy0
        h0_ll(3,4)=gb_yz0
        h0_ll(4,4)=gb_zz0

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'h0_ll(1,1)=',h0_ll(1,1)
!       write (*,*) 'h0_ll(1,2)=',h0_ll(1,2)
!       write (*,*) 'h0_ll(1,3)=',h0_ll(1,3)
!       write (*,*) 'h0_ll(1,4)=',h0_ll(1,4)
!       write (*,*) 'h0_ll(2,2)=',h0_ll(2,2)
!       write (*,*) 'h0_ll(2,3)=',h0_ll(2,3)
!       write (*,*) 'h0_ll(2,4)=',h0_ll(2,4)
!       write (*,*) 'h0_ll(3,3)=',h0_ll(3,3)
!       write (*,*) 'h0_ll(3,4)=',h0_ll(3,4)
!       write (*,*) 'h0_ll(4,4)=',h0_ll(4,4)
  
        h0_uu(1,1)=g0_uu(1,1)-gads_uu(1,1)
        h0_uu(1,2)=g0_uu(1,2)
        h0_uu(1,3)=g0_uu(1,3)
        h0_uu(1,4)=g0_uu(1,4)
        h0_uu(2,2)=g0_uu(2,2)-gads_uu(2,2)
        h0_uu(2,3)=g0_uu(2,3)-gads_uu(2,3)
        h0_uu(2,4)=g0_uu(2,4)-gads_uu(2,4)
        h0_uu(3,3)=g0_uu(3,3)-gads_uu(3,3)
        h0_uu(3,4)=g0_uu(3,4)-gads_uu(3,4)
        h0_uu(4,4)=g0_uu(4,4)-gads_uu(4,4)

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'h0_uu(1,1)=',h0_uu(1,1)
!       write (*,*) 'h0_uu(1,2)=',h0_uu(1,2)
!       write (*,*) 'h0_uu(1,3)=',h0_uu(1,3)
!       write (*,*) 'h0_uu(1,4)=',h0_uu(1,4)
!       write (*,*) 'h0_uu(2,2)=',h0_uu(2,2)
!       write (*,*) 'h0_uu(2,3)=',h0_uu(2,3)
!       write (*,*) 'h0_uu(2,4)=',h0_uu(2,4)
!       write (*,*) 'h0_uu(3,3)=',h0_uu(3,3)
!       write (*,*) 'h0_uu(3,4)=',h0_uu(3,4)
!       write (*,*) 'h0_uu(4,4)=',h0_uu(4,4)

        h0_ll_x(1,1,1)   =g0_ll_x(1,1,1)
        h0_ll_x(1,1,2)   =g0_ll_x(1,1,2)-gads_ll_x(1,1,2)
        h0_ll_x(1,1,3)   =g0_ll_x(1,1,3)-gads_ll_x(1,1,3)
        h0_ll_x(1,1,4)   =g0_ll_x(1,1,4)-gads_ll_x(1,1,4)
        h0_ll_xx(1,1,1,1)=g0_ll_xx(1,1,1,1)
        h0_ll_xx(1,1,1,2)=g0_ll_xx(1,1,1,2)
        h0_ll_xx(1,1,1,3)=g0_ll_xx(1,1,1,3)
        h0_ll_xx(1,1,1,4)=g0_ll_xx(1,1,1,4)
        h0_ll_xx(1,1,2,2)=g0_ll_xx(1,1,2,2)-gads_ll_xx(1,1,2,2)
        h0_ll_xx(1,1,2,3)=g0_ll_xx(1,1,2,3)-gads_ll_xx(1,1,2,3)
        h0_ll_xx(1,1,2,4)=g0_ll_xx(1,1,2,4)-gads_ll_xx(1,1,2,4)
        h0_ll_xx(1,1,3,3)=g0_ll_xx(1,1,3,3)-gads_ll_xx(1,1,3,3)
        h0_ll_xx(1,1,3,4)=g0_ll_xx(1,1,3,4)-gads_ll_xx(1,1,3,4)
        h0_ll_xx(1,1,4,4)=g0_ll_xx(1,1,4,4)-gads_ll_xx(1,1,4,4)

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'h0_ll_x(1,1,1)=',h0_ll_x(1,1,1)
!       write (*,*) 'h0_ll_x(1,1,2)=',h0_ll_x(1,1,2)
!       write (*,*) 'h0_ll_x(1,1,3)=',h0_ll_x(1,1,3)
!       write (*,*) 'h0_ll_x(1,1,4)=',h0_ll_x(1,1,4)
!       write (*,*) 'h0_ll_xx(1,1,1,1)=',h0_ll_xx(1,1,1,1)
!       write (*,*) 'h0_ll_xx(1,1,1,2)=',h0_ll_xx(1,1,1,2)
!       write (*,*) 'h0_ll_xx(1,1,1,3)=',h0_ll_xx(1,1,1,3)
!       write (*,*) 'h0_ll_xx(1,1,1,4)=',h0_ll_xx(1,1,1,4)
!       write (*,*) 'h0_ll_xx(1,1,2,2)=',h0_ll_xx(1,1,2,2)
!       write (*,*) 'h0_ll_xx(1,1,2,3)=',h0_ll_xx(1,1,2,3)
!       write (*,*) 'h0_ll_xx(1,1,2,4)=',h0_ll_xx(1,1,2,4)
!       write (*,*) 'h0_ll_xx(1,1,3,3)=',h0_ll_xx(1,1,3,3)
!       write (*,*) 'h0_ll_xx(1,1,3,4)=',h0_ll_xx(1,1,3,4)
!       write (*,*) 'h0_ll_xx(1,1,4,4)=',h0_ll_xx(1,1,4,4)

        h0_ll_x(1,2,1)   =g0_ll_x(1,2,1)
        h0_ll_x(1,2,2)   =g0_ll_x(1,2,2)
        h0_ll_x(1,2,3)   =g0_ll_x(1,2,3)
        h0_ll_x(1,2,4)   =g0_ll_x(1,2,4)
        h0_ll_xx(1,2,1,1)=g0_ll_xx(1,2,1,1)
        h0_ll_xx(1,2,1,2)=g0_ll_xx(1,2,1,2)
        h0_ll_xx(1,2,1,3)=g0_ll_xx(1,2,1,3)
        h0_ll_xx(1,2,1,4)=g0_ll_xx(1,2,1,4)
        h0_ll_xx(1,2,2,2)=g0_ll_xx(1,2,2,2)
        h0_ll_xx(1,2,2,3)=g0_ll_xx(1,2,2,3)
        h0_ll_xx(1,2,2,4)=g0_ll_xx(1,2,2,4)
        h0_ll_xx(1,2,3,3)=g0_ll_xx(1,2,3,3)
        h0_ll_xx(1,2,3,4)=g0_ll_xx(1,2,3,4)
        h0_ll_xx(1,2,4,4)=g0_ll_xx(1,2,4,4)

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'h0_ll_x(1,2,1)=',h0_ll_x(1,2,1)
!       write (*,*) 'h0_ll_x(1,2,2)=',h0_ll_x(1,2,2)
!       write (*,*) 'h0_ll_x(1,2,3)=',h0_ll_x(1,2,3)
!       write (*,*) 'h0_ll_x(1,2,4)=',h0_ll_x(1,2,4)
!       write (*,*) 'h0_ll_xx(1,2,1,1)=',h0_ll_xx(1,2,1,1)
!       write (*,*) 'h0_ll_xx(1,2,1,2)=',h0_ll_xx(1,2,1,2)
!       write (*,*) 'h0_ll_xx(1,2,1,3)=',h0_ll_xx(1,2,1,3)
!       write (*,*) 'h0_ll_xx(1,2,1,4)=',h0_ll_xx(1,2,1,4)
!       write (*,*) 'h0_ll_xx(1,2,2,2)=',h0_ll_xx(1,2,2,2)
!       write (*,*) 'h0_ll_xx(1,2,2,3)=',h0_ll_xx(1,2,2,3)
!       write (*,*) 'h0_ll_xx(1,2,2,4)=',h0_ll_xx(1,2,2,4)
!       write (*,*) 'h0_ll_xx(1,2,3,3)=',h0_ll_xx(1,2,3,3)
!       write (*,*) 'h0_ll_xx(1,2,3,4)=',h0_ll_xx(1,2,3,4)
!       write (*,*) 'h0_ll_xx(1,2,4,4)=',h0_ll_xx(1,2,4,4)

        h0_ll_x(1,3,1)   =g0_ll_x(1,3,1)
        h0_ll_x(1,3,2)   =g0_ll_x(1,3,2)
        h0_ll_x(1,3,3)   =g0_ll_x(1,3,3)
        h0_ll_x(1,3,4)   =g0_ll_x(1,3,4)
        h0_ll_xx(1,3,1,1)=g0_ll_xx(1,3,1,1)
        h0_ll_xx(1,3,1,2)=g0_ll_xx(1,3,1,2)
        h0_ll_xx(1,3,1,3)=g0_ll_xx(1,3,1,3)
        h0_ll_xx(1,3,1,4)=g0_ll_xx(1,3,1,4)
        h0_ll_xx(1,3,2,2)=g0_ll_xx(1,3,2,2)
        h0_ll_xx(1,3,2,3)=g0_ll_xx(1,3,2,3)
        h0_ll_xx(1,3,2,4)=g0_ll_xx(1,3,2,4)
        h0_ll_xx(1,3,3,3)=g0_ll_xx(1,3,3,3)
        h0_ll_xx(1,3,3,4)=g0_ll_xx(1,3,3,4)
        h0_ll_xx(1,3,4,4)=g0_ll_xx(1,3,4,4)

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'h0_ll_x(1,3,1)=',h0_ll_x(1,3,1)
!       write (*,*) 'h0_ll_x(1,3,2)=',h0_ll_x(1,3,2)
!       write (*,*) 'h0_ll_x(1,3,3)=',h0_ll_x(1,3,3)
!       write (*,*) 'h0_ll_x(1,3,4)=',h0_ll_x(1,3,4)
!       write (*,*) 'h0_ll_xx(1,3,1,1)=',h0_ll_xx(1,3,1,1)
!       write (*,*) 'h0_ll_xx(1,3,1,2)=',h0_ll_xx(1,3,1,2)
!       write (*,*) 'h0_ll_xx(1,3,1,3)=',h0_ll_xx(1,3,1,3)
!       write (*,*) 'h0_ll_xx(1,3,1,4)=',h0_ll_xx(1,3,1,4)
!       write (*,*) 'h0_ll_xx(1,3,2,2)=',h0_ll_xx(1,3,2,2)
!       write (*,*) 'h0_ll_xx(1,3,2,3)=',h0_ll_xx(1,3,2,3)
!       write (*,*) 'h0_ll_xx(1,3,2,4)=',h0_ll_xx(1,3,2,4)
!       write (*,*) 'h0_ll_xx(1,3,3,3)=',h0_ll_xx(1,3,3,3)
!       write (*,*) 'h0_ll_xx(1,3,3,4)=',h0_ll_xx(1,3,3,4)
!       write (*,*) 'h0_ll_xx(1,3,4,4)=',h0_ll_xx(1,3,4,4)

        h0_ll_x(1,4,1)   =g0_ll_x(1,4,1)
        h0_ll_x(1,4,2)   =g0_ll_x(1,4,2)
        h0_ll_x(1,4,3)   =g0_ll_x(1,4,3)
        h0_ll_x(1,4,4)   =g0_ll_x(1,4,4)
        h0_ll_xx(1,4,1,1)=g0_ll_xx(1,4,1,1)
        h0_ll_xx(1,4,1,2)=g0_ll_xx(1,4,1,2)
        h0_ll_xx(1,4,1,3)=g0_ll_xx(1,4,1,3)
        h0_ll_xx(1,4,1,4)=g0_ll_xx(1,4,1,4)
        h0_ll_xx(1,4,2,2)=g0_ll_xx(1,4,2,2)
        h0_ll_xx(1,4,2,3)=g0_ll_xx(1,4,2,3)
        h0_ll_xx(1,4,2,4)=g0_ll_xx(1,4,2,4)
        h0_ll_xx(1,4,3,3)=g0_ll_xx(1,4,3,3)
        h0_ll_xx(1,4,3,4)=g0_ll_xx(1,4,3,4)
        h0_ll_xx(1,4,4,4)=g0_ll_xx(1,4,4,4)

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'h0_ll_x(1,4,1)=',h0_ll_x(1,4,1)
!       write (*,*) 'h0_ll_x(1,4,2)=',h0_ll_x(1,4,2)
!       write (*,*) 'h0_ll_x(1,4,3)=',h0_ll_x(1,4,3)
!       write (*,*) 'h0_ll_x(1,4,4)=',h0_ll_x(1,4,4)
!       write (*,*) 'h0_ll_xx(1,4,1,1)=',h0_ll_xx(1,4,1,1)
!       write (*,*) 'h0_ll_xx(1,4,1,2)=',h0_ll_xx(1,4,1,2)
!       write (*,*) 'h0_ll_xx(1,4,1,3)=',h0_ll_xx(1,4,1,3)
!       write (*,*) 'h0_ll_xx(1,4,1,4)=',h0_ll_xx(1,4,1,4)
!       write (*,*) 'h0_ll_xx(1,4,2,2)=',h0_ll_xx(1,4,2,2)
!       write (*,*) 'h0_ll_xx(1,4,2,3)=',h0_ll_xx(1,4,2,3)
!       write (*,*) 'h0_ll_xx(1,4,2,4)=',h0_ll_xx(1,4,2,4)
!       write (*,*) 'h0_ll_xx(1,4,3,3)=',h0_ll_xx(1,4,3,3)
!       write (*,*) 'h0_ll_xx(1,4,3,4)=',h0_ll_xx(1,4,3,4)
!       write (*,*) 'h0_ll_xx(1,4,4,4)=',h0_ll_xx(1,4,4,4)


        h0_ll_x(2,2,1)   =g0_ll_x(2,2,1)
        h0_ll_x(2,2,2)   =g0_ll_x(2,2,2)-gads_ll_x(2,2,2)
        h0_ll_x(2,2,3)   =g0_ll_x(2,2,3)-gads_ll_x(2,2,3)
        h0_ll_x(2,2,4)   =g0_ll_x(2,2,4)-gads_ll_x(2,2,4)
        h0_ll_xx(2,2,1,1)=g0_ll_xx(2,2,1,1)
        h0_ll_xx(2,2,1,2)=g0_ll_xx(2,2,1,2)
        h0_ll_xx(2,2,1,3)=g0_ll_xx(2,2,1,3)
        h0_ll_xx(2,2,1,4)=g0_ll_xx(2,2,1,4)
        h0_ll_xx(2,2,2,2)=g0_ll_xx(2,2,2,2)-gads_ll_xx(2,2,2,2)
        h0_ll_xx(2,2,2,3)=g0_ll_xx(2,2,2,3)-gads_ll_xx(2,2,2,3)
        h0_ll_xx(2,2,2,4)=g0_ll_xx(2,2,2,4)-gads_ll_xx(2,2,2,4)
        h0_ll_xx(2,2,3,3)=g0_ll_xx(2,2,3,3)-gads_ll_xx(2,2,3,3)
        h0_ll_xx(2,2,3,4)=g0_ll_xx(2,2,3,4)-gads_ll_xx(2,2,3,4)
        h0_ll_xx(2,2,4,4)=g0_ll_xx(2,2,4,4)-gads_ll_xx(2,2,4,4)

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'h0_ll_x(2,2,1)=',h0_ll_x(2,2,1)
!       write (*,*) 'h0_ll_x(2,2,2)=',h0_ll_x(2,2,2)
!       write (*,*) 'h0_ll_x(2,2,3)=',h0_ll_x(2,2,3)
!       write (*,*) 'h0_ll_x(2,2,4)=',h0_ll_x(2,2,4)
!       write (*,*) 'h0_ll_xx(2,2,1,1)=',h0_ll_xx(2,2,1,1)
!       write (*,*) 'h0_ll_xx(2,2,1,2)=',h0_ll_xx(2,2,1,2)
!       write (*,*) 'h0_ll_xx(2,2,1,3)=',h0_ll_xx(2,2,1,3)
!       write (*,*) 'h0_ll_xx(2,2,1,4)=',h0_ll_xx(2,2,1,4)
!       write (*,*) 'h0_ll_xx(2,2,2,2)=',h0_ll_xx(2,2,2,2)
!       write (*,*) 'h0_ll_xx(2,2,2,3)=',h0_ll_xx(2,2,2,3)
!       write (*,*) 'h0_ll_xx(2,2,2,4)=',h0_ll_xx(2,2,2,4)
!       write (*,*) 'h0_ll_xx(2,2,3,3)=',h0_ll_xx(2,2,3,3)
!       write (*,*) 'h0_ll_xx(2,2,3,4)=',h0_ll_xx(2,2,3,4)
!       write (*,*) 'h0_ll_xx(2,2,4,4)=',h0_ll_xx(2,2,4,4)

        h0_ll_x(2,3,1)   =g0_ll_x(2,3,1)
        h0_ll_x(2,3,2)   =g0_ll_x(2,3,2)-gads_ll_x(2,3,2)
        h0_ll_x(2,3,3)   =g0_ll_x(2,3,3)-gads_ll_x(2,3,3)
        h0_ll_x(2,3,4)   =g0_ll_x(2,3,4)-gads_ll_x(2,3,4)
        h0_ll_xx(2,3,1,1)=g0_ll_xx(2,3,1,1)
        h0_ll_xx(2,3,1,2)=g0_ll_xx(2,3,1,2)
        h0_ll_xx(2,3,1,3)=g0_ll_xx(2,3,1,3)
        h0_ll_xx(2,3,1,4)=g0_ll_xx(2,3,1,4)
        h0_ll_xx(2,3,2,2)=g0_ll_xx(2,3,2,2)-gads_ll_xx(2,3,2,2)
        h0_ll_xx(2,3,2,3)=g0_ll_xx(2,3,2,3)-gads_ll_xx(2,3,2,3)
        h0_ll_xx(2,3,2,4)=g0_ll_xx(2,3,2,4)-gads_ll_xx(2,3,2,4)
        h0_ll_xx(2,3,3,3)=g0_ll_xx(2,3,3,3)-gads_ll_xx(2,3,3,3)
        h0_ll_xx(2,3,3,4)=g0_ll_xx(2,3,3,4)-gads_ll_xx(2,3,3,4)
        h0_ll_xx(2,3,4,4)=g0_ll_xx(2,3,4,4)-gads_ll_xx(2,3,4,4)

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'h0_ll_x(2,3,1)=',h0_ll_x(2,3,1)
!       write (*,*) 'h0_ll_x(2,3,2)=',h0_ll_x(2,3,2)
!       write (*,*) 'h0_ll_x(2,3,3)=',h0_ll_x(2,3,3)
!       write (*,*) 'h0_ll_x(2,3,4)=',h0_ll_x(2,3,4)
!       write (*,*) 'h0_ll_xx(2,3,1,1)=',h0_ll_xx(2,3,1,1)
!       write (*,*) 'h0_ll_xx(2,3,1,2)=',h0_ll_xx(2,3,1,2)
!       write (*,*) 'h0_ll_xx(2,3,1,3)=',h0_ll_xx(2,3,1,3)
!       write (*,*) 'h0_ll_xx(2,3,1,4)=',h0_ll_xx(2,3,1,4)
!       write (*,*) 'h0_ll_xx(2,3,2,2)=',h0_ll_xx(2,3,2,2)
!       write (*,*) 'h0_ll_xx(2,3,2,3)=',h0_ll_xx(2,3,2,3)
!       write (*,*) 'h0_ll_xx(2,3,2,4)=',h0_ll_xx(2,3,2,4)
!       write (*,*) 'h0_ll_xx(2,3,3,3)=',h0_ll_xx(2,3,3,3)
!       write (*,*) 'h0_ll_xx(2,3,3,4)=',h0_ll_xx(2,3,3,4)
!       write (*,*) 'h0_ll_xx(2,3,4,4)=',h0_ll_xx(2,3,4,4)

        h0_ll_x(2,4,1)   =g0_ll_x(2,4,1)
        h0_ll_x(2,4,2)   =g0_ll_x(2,4,2)-gads_ll_x(2,4,2)
        h0_ll_x(2,4,3)   =g0_ll_x(2,4,3)-gads_ll_x(2,4,3)
        h0_ll_x(2,4,4)   =g0_ll_x(2,4,4)-gads_ll_x(2,4,4)
        h0_ll_xx(2,4,1,1)=g0_ll_xx(2,4,1,1)
        h0_ll_xx(2,4,1,2)=g0_ll_xx(2,4,1,2)
        h0_ll_xx(2,4,1,3)=g0_ll_xx(2,4,1,3)
        h0_ll_xx(2,4,1,4)=g0_ll_xx(2,4,1,4)
        h0_ll_xx(2,4,2,2)=g0_ll_xx(2,4,2,2)-gads_ll_xx(2,4,2,2)
        h0_ll_xx(2,4,2,3)=g0_ll_xx(2,4,2,3)-gads_ll_xx(2,4,2,3)
        h0_ll_xx(2,4,2,4)=g0_ll_xx(2,4,2,4)-gads_ll_xx(2,4,2,4)
        h0_ll_xx(2,4,3,3)=g0_ll_xx(2,4,3,3)-gads_ll_xx(2,4,3,3)
        h0_ll_xx(2,4,3,4)=g0_ll_xx(2,4,3,4)-gads_ll_xx(2,4,3,4)
        h0_ll_xx(2,4,4,4)=g0_ll_xx(2,4,4,4)-gads_ll_xx(2,4,4,4)

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'h0_ll_x(2,4,1)=',h0_ll_x(2,4,1)
!       write (*,*) 'h0_ll_x(2,4,2)=',h0_ll_x(2,4,2)
!       write (*,*) 'h0_ll_x(2,4,3)=',h0_ll_x(2,4,3)
!       write (*,*) 'h0_ll_x(2,4,4)=',h0_ll_x(2,4,4)
!       write (*,*) 'h0_ll_xx(2,4,1,1)=',h0_ll_xx(2,4,1,1)
!       write (*,*) 'h0_ll_xx(2,4,1,2)=',h0_ll_xx(2,4,1,2)
!       write (*,*) 'h0_ll_xx(2,4,1,3)=',h0_ll_xx(2,4,1,3)
!       write (*,*) 'h0_ll_xx(2,4,1,4)=',h0_ll_xx(2,4,1,4)
!       write (*,*) 'h0_ll_xx(2,4,2,2)=',h0_ll_xx(2,4,2,2)
!       write (*,*) 'h0_ll_xx(2,4,2,3)=',h0_ll_xx(2,4,2,3)
!       write (*,*) 'h0_ll_xx(2,4,2,4)=',h0_ll_xx(2,4,2,4)
!       write (*,*) 'h0_ll_xx(2,4,3,3)=',h0_ll_xx(2,4,3,3)
!       write (*,*) 'h0_ll_xx(2,4,3,4)=',h0_ll_xx(2,4,3,4)
!       write (*,*) 'h0_ll_xx(2,4,4,4)=',h0_ll_xx(2,4,4,4)

        h0_ll_x(3,3,1)   =g0_ll_x(3,3,1)
        h0_ll_x(3,3,2)   =g0_ll_x(3,3,2)-gads_ll_x(3,3,2)
        h0_ll_x(3,3,3)   =g0_ll_x(3,3,3)-gads_ll_x(3,3,3)
        h0_ll_x(3,3,4)   =g0_ll_x(3,3,4)-gads_ll_x(3,3,4)
        h0_ll_xx(3,3,1,1)=g0_ll_xx(3,3,1,1)
        h0_ll_xx(3,3,1,2)=g0_ll_xx(3,3,1,2)
        h0_ll_xx(3,3,1,3)=g0_ll_xx(3,3,1,3)
        h0_ll_xx(3,3,1,4)=g0_ll_xx(3,3,1,4)
        h0_ll_xx(3,3,2,2)=g0_ll_xx(3,3,2,2)-gads_ll_xx(3,3,2,2)
        h0_ll_xx(3,3,2,3)=g0_ll_xx(3,3,2,3)-gads_ll_xx(3,3,2,3)
        h0_ll_xx(3,3,2,4)=g0_ll_xx(3,3,2,4)-gads_ll_xx(3,3,2,4)
        h0_ll_xx(3,3,3,3)=g0_ll_xx(3,3,3,3)-gads_ll_xx(3,3,3,3)
        h0_ll_xx(3,3,3,4)=g0_ll_xx(3,3,3,4)-gads_ll_xx(3,3,3,4)
        h0_ll_xx(3,3,4,4)=g0_ll_xx(3,3,4,4)-gads_ll_xx(3,3,4,4)

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'h0_ll_x(3,3,1)=',h0_ll_x(3,3,1)
!       write (*,*) 'h0_ll_x(3,3,2)=',h0_ll_x(3,3,2)
!       write (*,*) 'h0_ll_x(3,3,3)=',h0_ll_x(3,3,3)
!       write (*,*) 'h0_ll_x(3,3,4)=',h0_ll_x(3,3,4)
!       write (*,*) 'h0_ll_xx(3,3,1,1)=',h0_ll_xx(3,3,1,1)
!       write (*,*) 'h0_ll_xx(3,3,1,2)=',h0_ll_xx(3,3,1,2)
!       write (*,*) 'h0_ll_xx(3,3,1,3)=',h0_ll_xx(3,3,1,3)
!       write (*,*) 'h0_ll_xx(3,3,1,4)=',h0_ll_xx(3,3,1,4)
!       write (*,*) 'h0_ll_xx(3,3,2,2)=',h0_ll_xx(3,3,2,2)
!       write (*,*) 'h0_ll_xx(3,3,2,3)=',h0_ll_xx(3,3,2,3)
!       write (*,*) 'h0_ll_xx(3,3,2,4)=',h0_ll_xx(3,3,2,4)
!       write (*,*) 'h0_ll_xx(3,3,3,3)=',h0_ll_xx(3,3,3,3)
!       write (*,*) 'h0_ll_xx(3,3,3,4)=',h0_ll_xx(3,3,3,4)
!       write (*,*) 'h0_ll_xx(3,3,4,4)=',h0_ll_xx(3,3,4,4)

        h0_ll_x(3,4,1)   =g0_ll_x(3,4,1)
        h0_ll_x(3,4,2)   =g0_ll_x(3,4,2)-gads_ll_x(3,4,2)
        h0_ll_x(3,4,3)   =g0_ll_x(3,4,3)-gads_ll_x(3,4,3)
        h0_ll_x(3,4,4)   =g0_ll_x(3,4,4)-gads_ll_x(3,4,4)
        h0_ll_xx(3,4,1,1)=g0_ll_xx(3,4,1,1)
        h0_ll_xx(3,4,1,2)=g0_ll_xx(3,4,1,2)
        h0_ll_xx(3,4,1,3)=g0_ll_xx(3,4,1,3)
        h0_ll_xx(3,4,1,4)=g0_ll_xx(3,4,1,4)
        h0_ll_xx(3,4,2,2)=g0_ll_xx(3,4,2,2)-gads_ll_xx(3,4,2,2)
        h0_ll_xx(3,4,2,3)=g0_ll_xx(3,4,2,3)-gads_ll_xx(3,4,2,3)
        h0_ll_xx(3,4,2,4)=g0_ll_xx(3,4,2,4)-gads_ll_xx(3,4,2,4)
        h0_ll_xx(3,4,3,3)=g0_ll_xx(3,4,3,3)-gads_ll_xx(3,4,3,3)
        h0_ll_xx(3,4,3,4)=g0_ll_xx(3,4,3,4)-gads_ll_xx(3,4,3,4)
        h0_ll_xx(3,4,4,4)=g0_ll_xx(3,4,4,4)-gads_ll_xx(3,4,4,4)

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'h0_ll_x(3,4,1)=',h0_ll_x(3,4,1)
!       write (*,*) 'h0_ll_x(3,4,2)=',h0_ll_x(3,4,2)
!       write (*,*) 'h0_ll_x(3,4,3)=',h0_ll_x(3,4,3)
!       write (*,*) 'h0_ll_x(3,4,4)=',h0_ll_x(3,4,4)
!       write (*,*) 'h0_ll_xx(3,4,1,1)=',h0_ll_xx(3,4,1,1)
!       write (*,*) 'h0_ll_xx(3,4,1,2)=',h0_ll_xx(3,4,1,2)
!       write (*,*) 'h0_ll_xx(3,4,1,3)=',h0_ll_xx(3,4,1,3)
!       write (*,*) 'h0_ll_xx(3,4,1,4)=',h0_ll_xx(3,4,1,4)
!       write (*,*) 'h0_ll_xx(3,4,2,2)=',h0_ll_xx(3,4,2,2)
!       write (*,*) 'h0_ll_xx(3,4,2,3)=',h0_ll_xx(3,4,2,3)
!       write (*,*) 'h0_ll_xx(3,4,2,4)=',h0_ll_xx(3,4,2,4)
!       write (*,*) 'h0_ll_xx(3,4,3,3)=',h0_ll_xx(3,4,3,3)
!       write (*,*) 'h0_ll_xx(3,4,3,4)=',h0_ll_xx(3,4,3,4)
!       write (*,*) 'h0_ll_xx(3,4,4,4)=',h0_ll_xx(3,4,4,4)

        h0_ll_x(4,4,1)   =g0_ll_x(4,4,1)
        h0_ll_x(4,4,2)   =g0_ll_x(4,4,2)-gads_ll_x(4,4,2)
        h0_ll_x(4,4,3)   =g0_ll_x(4,4,3)-gads_ll_x(4,4,3)
        h0_ll_x(4,4,4)   =g0_ll_x(4,4,4)-gads_ll_x(4,4,4)
        h0_ll_xx(4,4,1,1)=g0_ll_xx(4,4,1,1)
        h0_ll_xx(4,4,1,2)=g0_ll_xx(4,4,1,2)
        h0_ll_xx(4,4,1,3)=g0_ll_xx(4,4,1,3)
        h0_ll_xx(4,4,1,4)=g0_ll_xx(4,4,1,4)
        h0_ll_xx(4,4,2,2)=g0_ll_xx(4,4,2,2)-gads_ll_xx(4,4,2,2)
        h0_ll_xx(4,4,2,3)=g0_ll_xx(4,4,2,3)-gads_ll_xx(4,4,2,3)
        h0_ll_xx(4,4,2,4)=g0_ll_xx(4,4,2,4)-gads_ll_xx(4,4,2,4)
        h0_ll_xx(4,4,3,3)=g0_ll_xx(4,4,3,3)-gads_ll_xx(4,4,3,3)
        h0_ll_xx(4,4,3,4)=g0_ll_xx(4,4,3,4)-gads_ll_xx(4,4,3,4)
        h0_ll_xx(4,4,4,4)=g0_ll_xx(4,4,4,4)-gads_ll_xx(4,4,4,4)

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0,dx=',L,i,j,k,x0,y0,z0,rho0,dx
!       write (*,*) 'h0_ll_x(4,4,1)=',h0_ll_x(4,4,1)
!       write (*,*) 'h0_ll_x(4,4,2)=',h0_ll_x(4,4,2)
!       write (*,*) 'h0_ll_x(4,4,3)=',h0_ll_x(4,4,3)
!       write (*,*) 'h0_ll_x(4,4,4)=',h0_ll_x(4,4,4)
!       write (*,*) 'h0_ll_xx(4,4,1,1)=',h0_ll_xx(4,4,1,1)
!       write (*,*) 'h0_ll_xx(4,4,1,2)=',h0_ll_xx(4,4,1,2)
!       write (*,*) 'h0_ll_xx(4,4,1,3)=',h0_ll_xx(4,4,1,3)
!       write (*,*) 'h0_ll_xx(4,4,1,4)=',h0_ll_xx(4,4,1,4)
!       write (*,*) 'h0_ll_xx(4,4,2,2)=',h0_ll_xx(4,4,2,2)
!       write (*,*) 'h0_ll_xx(4,4,2,3)=',h0_ll_xx(4,4,2,3)
!       write (*,*) 'h0_ll_xx(4,4,2,4)=',h0_ll_xx(4,4,2,4)
!       write (*,*) 'h0_ll_xx(4,4,3,3)=',h0_ll_xx(4,4,3,3)
!       write (*,*) 'h0_ll_xx(4,4,3,4)=',h0_ll_xx(4,4,3,4)
!       write (*,*) 'h0_ll_xx(4,4,4,4)=',h0_ll_xx(4,4,4,4)

        h0_uu_x(1,1,1)=g0_uu_x(1,1,1)
        h0_uu_x(1,1,2)=g0_uu_x(1,1,2)-gads_uu_x(1,1,2)
        h0_uu_x(1,1,3)=g0_uu_x(1,1,3)-gads_uu_x(1,1,3)
        h0_uu_x(1,1,4)=g0_uu_x(1,1,4)-gads_uu_x(1,1,4)

        h0_uu_x(1,2,1)=g0_uu_x(1,2,1)
        h0_uu_x(1,2,2)=g0_uu_x(1,2,2)
        h0_uu_x(1,2,3)=g0_uu_x(1,2,3)
        h0_uu_x(1,2,4)=g0_uu_x(1,2,4)

        h0_uu_x(1,3,1)=g0_uu_x(1,3,1) 
        h0_uu_x(1,3,2)=g0_uu_x(1,3,2) 
        h0_uu_x(1,3,3)=g0_uu_x(1,3,3) 
        h0_uu_x(1,3,4)=g0_uu_x(1,3,4)

        h0_uu_x(1,4,1)=g0_uu_x(1,4,1)
        h0_uu_x(1,4,2)=g0_uu_x(1,4,2)
        h0_uu_x(1,4,3)=g0_uu_x(1,4,3)
        h0_uu_x(1,4,4)=g0_uu_x(1,4,4)

        h0_uu_x(2,2,1)=g0_uu_x(2,2,1)
        h0_uu_x(2,2,2)=g0_uu_x(2,2,2)-gads_uu_x(2,2,2)
        h0_uu_x(2,2,3)=g0_uu_x(2,2,3)-gads_uu_x(2,2,3)
        h0_uu_x(2,2,4)=g0_uu_x(2,2,4)-gads_uu_x(2,2,4)

        h0_uu_x(2,3,1)=g0_uu_x(2,3,1)
        h0_uu_x(2,3,2)=g0_uu_x(2,3,2)-gads_uu_x(2,3,2)
        h0_uu_x(2,3,3)=g0_uu_x(2,3,3)-gads_uu_x(2,3,3)
        h0_uu_x(2,3,4)=g0_uu_x(2,3,4)-gads_uu_x(2,3,4)

        h0_uu_x(2,4,1)=g0_uu_x(2,4,1)
        h0_uu_x(2,4,2)=g0_uu_x(2,4,2)-gads_uu_x(2,4,2)
        h0_uu_x(2,4,3)=g0_uu_x(2,4,3)-gads_uu_x(2,4,3)
        h0_uu_x(2,4,4)=g0_uu_x(2,4,4)-gads_uu_x(2,4,4)

        h0_uu_x(3,3,1)=g0_uu_x(3,3,1)
        h0_uu_x(3,3,2)=g0_uu_x(3,3,2)-gads_uu_x(3,3,2)
        h0_uu_x(3,3,3)=g0_uu_x(3,3,3)-gads_uu_x(3,3,3)
        h0_uu_x(3,3,4)=g0_uu_x(3,3,4)-gads_uu_x(3,3,4)

        h0_uu_x(3,4,1)=g0_uu_x(3,4,1)
        h0_uu_x(3,4,2)=g0_uu_x(3,4,2)-gads_uu_x(3,4,2)
        h0_uu_x(3,4,3)=g0_uu_x(3,4,3)-gads_uu_x(3,4,3)
        h0_uu_x(3,4,4)=g0_uu_x(3,4,4)-gads_uu_x(3,4,4)

        h0_uu_x(4,4,1)=g0_uu_x(4,4,1)
        h0_uu_x(4,4,2)=g0_uu_x(4,4,2)-gads_uu_x(4,4,2)
        h0_uu_x(4,4,3)=g0_uu_x(4,4,3)-gads_uu_x(4,4,3)
        h0_uu_x(4,4,4)=g0_uu_x(4,4,4)-gads_uu_x(4,4,4)

        do a=1,3
          do b=a+1,4
            h0_ll(b,a)=h0_ll(a,b)
            h0_uu(b,a)=h0_uu(a,b)
            do c=1,4
              h0_ll_x(b,a,c)=h0_ll_x(a,b,c)
              h0_uu_x(b,a,c)=h0_uu_x(a,b,c)
            end do
          end do
        end do

        ! give values to the source functions
        A_l(1)=Hb_t0*(1-rho0**2)
        A_l(2)=Hb_x0*(1-rho0**2)
        A_l(3)=Hb_y0*(1-rho0**2)
        A_l(4)=Hb_z0*(1-rho0**2)

!      write(*,*) 'DEBUG from misc.f'
!      write(*,*) 'L,x0,y0,z0,rho0,dx=',L,x0,y0,z0,rho0,dx
!      write(*,*) 'A_l(1),A_l(2),A_l(3),A_l(4)='
!     &           ,A_l(1),A_l(2),A_l(3),A_l(4)

        A_l_x(1,1)=Hb_t_t*(1-rho0**2)
        A_l_x(1,2)=Hb_t_x*(1-rho0**2)
     &            -Hb_t0*2*x0
        A_l_x(1,3)=Hb_t_y*(1-rho0**2)
     &            -Hb_t0*2*y0
        A_l_x(1,4)=Hb_t_z*(1-rho0**2)
     &            -Hb_t0*2*z0

        A_l_x(2,1)=Hb_x_t*(1-rho0**2)
        A_l_x(2,2)=Hb_x_x*(1-rho0**2)
     &            -Hb_x0*2*x0
        A_l_x(2,3)=Hb_x_y*(1-rho0**2)
     &            -Hb_x0*2*y0
        A_l_x(2,4)=Hb_x_z*(1-rho0**2)
     &            -Hb_x0*2*z0

        A_l_x(3,1)=Hb_y_t*(1-rho0**2)
        A_l_x(3,2)=Hb_y_x*(1-rho0**2)
     &            -Hb_y0*2*x0
        A_l_x(3,3)=Hb_y_y*(1-rho0**2)
     &            -Hb_y0*2*y0
        A_l_x(3,4)=Hb_y_z*(1-rho0**2)
     &            -Hb_y0*2*z0

        A_l_x(4,1)=Hb_z_t*(1-rho0**2)
        A_l_x(4,2)=Hb_z_x*(1-rho0**2)
     &            -Hb_z0*2*x0
        A_l_x(4,3)=Hb_z_y*(1-rho0**2)
     &            -Hb_z0*2*y0
        A_l_x(4,4)=Hb_z_z*(1-rho0**2)
     &            -Hb_z0*2*z0

!      write (*,*) 'Hb_t0,Hb_t_t,Hb_t_x,Hb_t_y,Hb_t_z='
!     &            ,Hb_t0,Hb_t_t,Hb_t_x,Hb_t_y,Hb_t_z
!      write(*,*)  'Hb_t_np1(i,j,k),Hb_t_nm1(i,j,k)='
!     &            ,Hb_t_np1(i,j,k),Hb_t_nm1(i,j,k)
!      write(*,*) 'A_l_x(1,1),A_l_x(1,2),A_l_x(1,3),A_l_x(1,4)='
!     &           ,A_l_x(1,1),A_l_x(1,2),A_l_x(1,3),A_l_x(1,4)
!      write(*,*) 'A_l_x(2,1),A_l_x(2,2),A_l_x(2,3),A_l_x(2,4)='
!     &           ,A_l_x(2,1),A_l_x(2,2),A_l_x(2,3),A_l_x(2,4)
!      write(*,*) 'A_l_x(3,1),A_l_x(3,2),A_l_x(3,3),A_l_x(3,4)='
!     &           ,A_l_x(3,1),A_l_x(3,2),A_l_x(3,3),A_l_x(3,4)
!      write(*,*) 'A_l_x(4,1),A_l_x(4,2),A_l_x(4,3),A_l_x(4,4)='
!     &           ,A_l_x(4,1),A_l_x(4,2),A_l_x(4,3),A_l_x(4,4)

        ! give values to the AdS source functions
!!!!2=1 version!!!!!!
!        Hads_l(1)=0
!        Hads_l(2)=(-2*x0*(3 + rho0**2))/(-1 + (rho0**2)**2)
!        Hads_l(3)=1/y0 - (2*y0*(3 + rho0**2))/(-1 + (rho0**2)**2)
!        Hads_l(4)=0
!!!!!!!!!!!!!!!!!!!!!!

        Hads_l(1)=0
        Hads_l(2)=(2*x0*(-L**2*(5+x0**4+2*y0**2+2*z0**2+(y0**2+z0**2)**2
     &            +2*x0**2*(1+y0**2+z0**2))*(-1+rho0**2)
     &            -4*(2+rho0**2+rho0**4)))
     &            /((-1+rho0**2)*(1+rho0**2)
     &            *(4*rho0**2+L**2*(-1+rho0**2)**2))
        Hads_l(3)=(2*y0*(-L**2*(5+x0**4+2*y0**2+2*z0**2+(y0**2+z0**2)**2
     &            +2*x0**2*(1+y0**2+z0**2))*(-1+rho0**2)
     &            -4*(2+rho0**2+rho0**4)))
     &            /((-1+rho0**2)*(1+rho0**2)
     &            *(4*rho0**2+L**2*(-1+rho0**2)**2))
        Hads_l(4)=-((4*z0)/(-1+rho0**2))-(2*z0)/(1+rho0**2)
     &            +(4*z0*(4+L**2*(-3+rho0**2)))
     &            /(4*rho0**2+L**2*(-1+rho0**2)**2)

!       write (*,*) 'L,i,j,k,x0,y0,z0,rho0=',L,i,j,k,x0,y0,z0,rho0
!       write (*,*) ' Hads_l(2)=' ,Hads_l(2)
!       write (*,*) ' Hads_l(3)=' ,Hads_l(3)
!       write (*,*) ' Hads_l(4)=' ,Hads_l(4)


        ! give values to the scalar field
        phi10_x(1)=phi1_t*(1-rho0**2)**2
        phi10_x(2)=phi1_x*(1-rho0**2)**2
     &            +phi10*(-4*x0)*(1-rho0**2)
        phi10_x(3)=phi1_y*(1-rho0**2)**2
     &            +phi10*(-4*y0)*(1-rho0**2)
        phi10_x(4)=phi1_z*(1-rho0**2)**2
     &            +phi10*(-4*z0)*(1-rho0**2)

        phi10_xx(1,1)=phi1_tt*(1-rho0**2)**2
        phi10_xx(1,2)=phi1_tx*(1-rho0**2)**2
     &               +phi1_t*(-4*x0)*(1-rho0**2)
        phi10_xx(1,3)=phi1_ty*(1-rho0**2)**2
     &               +phi1_t*(-4*y0)*(1-rho0**2)
        phi10_xx(1,4)=phi1_tz*(1-rho0**2)**2
     &               +phi1_t*(-4*z0)*(1-rho0**2)

        phi10_xx(2,2)=phi1_xx*(1-rho0**2)**2
     &               +phi1_x*(2)*(-4*x0)*(1-rho0**2)
     &               +phi10*(-4*(1-rho0**2)+8*x0**2)
        phi10_xx(2,3)=phi1_xy*(1-rho0**2)**2
     &               +phi1_x*(-4*y0)*(1-rho0**2)
     &               +phi1_y*(-4*x0)*(1-rho0**2)
     &               +phi10*(-4*x0)*(-2*y0)
        phi10_xx(2,4)=phi1_xz*(1-rho0**2)**2
     &               +phi1_x*(-4*z0)*(1-rho0**2)
     &               +phi1_z*(-4*x0)*(1-rho0**2)
     &               +phi10*(-4*x0)*(-2*z0)

        phi10_xx(3,3)=phi1_yy*(1-rho0**2)**2
     &               +phi1_y*(2)*(-4*y0)*(1-rho0**2)
     &               +phi10*(-4*(1-rho0**2)+8*y0**2)
        phi10_xx(3,4)=phi1_yz*(1-rho0**2)**2
     &               +phi1_y*(-4*z0)*(1-rho0**2)
     &               +phi1_z*(-4*y0)*(1-rho0**2)
     &               +phi10*(-4*y0)*(-2*z0)

        phi10_xx(4,4)=phi1_zz*(1-rho0**2)**2
     &               +phi1_z*(2)*(-4*z0)*(1-rho0**2)
     &               +phi10*(-4*(1-rho0**2)+8*z0**2)

!!!CHECKED WITH MATHEMATICA UP TO HERE!!

        do a=1,3
          do b=a+1,4
            phi10_xx(b,a)=phi10_xx(a,b)
          end do
        end do

        ! calculate Riemann tensor at point i,j
        !(R^a_bcd =gamma^a_bd,c - gamma^a_bc,d
        !          +gamma^a_ce gamma^e_bd - gamma^a_de gamma^e_bc)
        do a=1,4
          do b=1,4
            do c=1,4
              do d=1,4
                riemann_ulll(a,b,c,d)=
     &                gamma_ull_x(a,b,d,c)-gamma_ull_x(a,b,c,d)
                do e=1,4
                   riemann_ulll(a,b,c,d)=riemann_ulll(a,b,c,d)
     &               +gamma_ull(a,c,e)*gamma_ull(e,b,d)
     &               -gamma_ull(a,d,e)*gamma_ull(e,b,c)
                end do
              end do
            end do
          end do
        end do

        ! calculate Ricci tensor at point i,j
        !(R_bd = R^a_bad)
        do b=1,4
          do d=1,4
            ricci_ll(b,d)=0
            do a=1,4
              ricci_ll(b,d)=ricci_ll(b,d)+riemann_ulll(a,b,a,d)
            end do
          end do
        end do

        ! calculate raised Ricci tensor at point i,j
        !(R_a^b = R_ad g^db)
        do a=1,4
          do b=1,4
            ricci_lu(a,b)=0
            do d=1,4
              ricci_lu(a,b)=ricci_lu(a,b)+ricci_ll(a,d)*g0_uu(d,b)
            end do
          end do
        end do

        ! calculate Ricci scalar
        !(R = R_a^a)
        ricci=0
        do a=1,4
          ricci=ricci+ricci_lu(a,a)
        end do
  
        ! calculates Einstein tensor at point i,j
        !(G_ab = R_ab - 1/2 R g_ab)
        do a=1,4
          do b=1,4
            einstein_ll(a,b)=ricci_ll(a,b)-0.5d0*ricci*g0_ll(a,b)
          end do
        end do

        ! calculates stress-energy tensor at point i,j 
        !(T_ab = 2*phi1,a phi1,b - (phi1,c phi1,d) g^cd g_ab + ...)
        grad_phi1_sq=0
        do a=1,4
          do b=1,4
            grad_phi1_sq=grad_phi1_sq
     &                  +phi10_x(a)*phi10_x(b)*g0_uu(a,b)
          end do
        end do

        do a=1,4
          do b=1,4
            set_ll(a,b)=
     &            phi10_x(a)*phi10_x(b)
     &           -g0_ll(a,b)*(grad_phi1_sq/2)
          end do
        end do

        return
        end

