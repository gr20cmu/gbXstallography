c=====================================================================
      subroutine VLen (v, temp)
c=====================================================================
      include 'common.fi'

       temp=((v(1)**2)+(v(2)**2)+(v(3)**2))**0.5

      return
      end

c=====================================================================
      subroutine normalize (x, v)
c=====================================================================
      include 'common.fi'
       call Vlen(x,temp)
       if (temp.eq.0.0) then
        v(1)=2.0
        v(2)=0.0
        v(3)=0.0
       else
        v(1)=x(1)/temp
        v(2)=x(2)/temp
        v(3)=x(3)/temp
       endif
       return
      end

c=====================================================================
      subroutine GToAA (g, v, temp1)
c=====================================================================
      include 'common.fi'

       temp1=acos2((g(1,1)+g(2,2)+g(3,3)-1.)/2.) ! modified to acos2 220222
       v(1)=(g(2,3)-g(3,2))/(2.*sin(temp1))
       v(2)=(g(3,1)-g(1,3))/(2.*sin(temp1))
       v(3)=(g(1,2)-g(2,1))/(2.*sin(temp1))

      return
      end

c=====================================================================
      subroutine DToRad3 (v, euler)
c=====================================================================
      include 'common.fi'

       euler(1)=v(1)*(pi/180.0)
       euler(2)=v(2)*(pi/180.0)
       euler(3)=v(3)*(pi/180.0)

      return
      end

c=====================================================================
      subroutine RadToD3 (v, euler)
c=====================================================================
      include 'common.fi'

       euler(1)=v(1)*(180.0/pi)
       euler(2)=v(2)*(180.0/pi)
       euler(3)=v(3)*(180.0/pi)

      return
      end

c=====================================================================
      subroutine Cross (n1, n2, vector_out)
c=====================================================================
      include 'common.fi'

       vector_out(1)=(n1(2)*n2(3))-(n1(3)*n2(2))
       vector_out(2)=(n1(3)*n2(1))-(n1(1)*n2(3))
       vector_out(3)=(n1(1)*n2(2))-(n1(2)*n2(1))

      return
      end

c=====================================================================
      subroutine AngTween (n1, n2, theta)
c=====================================================================
      include 'common.fi'

       call VLen(n1, l1)
       call VLen(n2, l2)
       temp = (n1(1)*n2(1))+(n1(2)*n2(2))+(n1(3)*n2(3))
       theta = acos2(temp/(l1*l2))

      return
      end

c=====================================================================
      subroutine DisG (O, nsy, n1, n2, g3nw, x, theta)
c=====================================================================
      include 'common.fi'
       !nsy is nsymm, n1, n2, are eulers, x is axis, theta is angle

       call EToG(n1,g1nw)
       call EToG(n2,g2nw)

       theta = 3.1415
       do i=1,nsy
        do j=1,nsy
         call symop(O,i,so1)
         call MToM (so1,g1nw,gg1)
         call symop(O,j,so2)
         call MToM(so2,g2nw,gg2)
         call trans(gg2,gtnw)
         call MToM(gg1,gtnw,delgnw)
         call GToAA(delgnw, v, mo)
         if (mo.lt.theta) then
          x=v
          g3nw=delgnw
          theta=mo
         endif
         call trans (gg1, gtnw)
         call MToM (gg2, gtnw, delgnw)
         call GToAA(delgnw, v, mo)
         if (mo.lt.theta) then
          x=v
          g3nw=delgnw
          theta=mo
         endif
        enddo
       enddo

      return
      end

c=====================================================================
      subroutine DisGFZ (O, ii, n1, n2, g3, x, theta)
c=====================================================================
      include 'common.fi'
       !this is same as DisG, but will produce an axis in the 
       !fundamental zone for cubic.

       call EToG(n1,g1)
       call EToG(n2,g2)

       theta = pi/2
       do i=1,ii
        do j=1,ii
         call symop(O,i,so1)
         call MToM (so1,g1,gg1)
         call symop(O,j,so2)
         call MToM(so2,g2,gg2)
         call trans(gg2,gt)
         call MToM(gg1,gt,delg)
         call GToAA(delg, v, temp)
         if (v(1).gt.0.0.and.v(2).gt.0.0.and.v(3).gt.0.0) then
          if (temp.lt.theta) then
           x=v
           g3=delg
           theta=temp
          endif
         endif
         call trans (gg1, gt)
         call MToM (gg2, gt, delg)
         call GToAA(delg, v, temp)
         if (v(1).gt.0.0.and.v(2).gt.0.0.and.v(3).gt.0.0) then
          if (temp.lt.theta) then
           x=v
           g3=delg
           theta=temp
          endif
         endif
        enddo
       enddo

      return
      end


!DisQ

