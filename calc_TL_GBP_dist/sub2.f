c version 03/26/2024
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

c=====================================================================
      subroutine GetFileSize (fname0,header0,nnline0)
c=====================================================================
      include 'common.fi'
      nnline0 = 0 ! Initialize counter
      open(50, file=trim(fname0), status='old')
 60   continue
      read(50,*,end=61) ! read line, if end of file, goto 61
      nnline0=nnline0+1  ! Increment counter
      goto 60
 61   continue
      close(50)  ! we are finished, close the file, report result.
      ! Check each line for the string '#'. If found, we know this line
      ! is part of the header.  If not, we have reached the data section.
      open(50, file=trim(fname0), status='old')! open the file
      header0 = 0 ! Initialize counter
      hash0 = 0   ! Initialize flag
      do i=1,nnline0  ! loop over each line in the file
       read(50, "(a)") inline0    ! reads line into a string
       hash0 = index(inline0, '#')! if the character "#" in the string, hash>1.
       if (hash0.ne.0) then       ! check the flag
	    header0=header0+1         ! when "#" found, it is a header line
        hash0=0                   ! reset the flag
       else
        goto 80 ! when hash is = 0, we have reached the end of the header
       endif
      enddo     ! closes i1=1,nnline
 80   close(50) ! close the file, report the results

      return
      end

c=====================================================================
      subroutine GetNumberedFileName (i,keyword0,extension0,fname0)
c=====================================================================
      include 'common.fi'

      msf=int(i/100) !determine the most significant digit in the three digit file number
      nsf=int((i-(msf*100))/10) !determine the next (second) most significant digit
      lsf=i-(100*msf)-(10*nsf)  !determine the least (third) most significant digit
      one=char(48+msf)              !make it a string
      two=char(48+nsf)              !make it a string
      three=char(48+lsf)            !make it a string
      !create the filename by concatenating the parts
      fname0=trim(keyword0)//'_'//one//two//three//trim(extension0)

      return
      end

c=====================================================================
      subroutine DetectRadians (fname0,header0,nnline0,rad0)
c=====================================================================
      include 'common.fi'
      rad0=1
      open(50, file=fname0, status='old')! open the file
      do i=1,header0
       read(50, "(a)") inline0    ! reads line into a string
      enddo
      do i=1,nnline0-header0
       read(50, *) euler(1),euler(2),euler(3)
       if (euler(1)/3.1415.gt.2.1) rad0=0
       if (euler(2)/3.1415.gt.2.1) rad0=0
       if (euler(3)/3.1415.gt.2.1) rad0=0
      enddo
      close(50)

      return
      end

c=====================================================================
      subroutine GetAngDimensions (fname0,header0,nnline0,step0,xmax0,ymax0)
c=====================================================================
      include 'common.fi'
      !This first part determines the step size in the x-direction
      open (50, file=fname0, status='old')
      x01=0
      x02=0
      do i=1,header0+2
       if (i.gt.header0) then
        read(50,*) r01,r02,r03,r04,r05
        if(i.eq.header0+1) x01=r04
        if(i.eq.header0+2) x02=r04
       else
        read(50, "(a)") inline0
       endif
      enddo
      close(50)
      step0 = x02 - x01
      !The next part determines the number of pixels in x and y
      open (50, file=fname0, status='old')
      xmax0=0
      ymax0=0
      do i=1,nnline0
       if (i.gt.header0) then
        read(50,*) r01,r02,r03,r04,r05
        if(r04.GT.xmax0) xmax0=r04   !find x and y max
        if(r05.GT.ymax0) ymax0=r05
       else
        read(50, "(a)") inline0
       endif
      enddo
      close(50)
      xmax0=int(xmax0/step0)+1
      ymax0=int(ymax0/step0)+1

      return
      end

c=====================================================================
      subroutine GetMaxDist (omapA,x1,y1,xmax0,ymax0,l2)
c=====================================================================
      include 'common.fi'
      !(orientation map data, x coord corner, y coord corner, max x, max y, MaxDist)
      l2=0.0
      do i=1,xmax0
       do j=1,ymax0
        if (omapA(i,j,6).eq.0.0) goto 100
        v(1)=omapA(i,j,4)-x1
        v(2)=omapA(i,j,5)-y1
        v(3)=0.0
        call VLen(v,l1)
        if (l1.gt.l2) l2=l1
 100    continue
       enddo
      enddo

      return
      end

c=====================================================================
      subroutine DistortImage (omapA,x1,y1,xmax0,ymax0,l1,TX0,TY0)
c=====================================================================
      include 'common.fi'
      !(orientation map data, max x, max y, x coord corner, y coord corner, MaxDist, trans x, trans y)
      omapB=omapA
      do i=1,ymax0
       do j=1,xmax0
        !if (omapA(j,i,6).eq.0.0) then
         !omapB(j,i,4)=omapA(j,i,4)
         !omapB(j,i,5)=omapA(j,i,5)
         !goto 110
        !endif
        v(1)=abs(omapA(j,i,4)-x1)
        v(2)=abs(omapA(j,i,5)-y1)
        v(3)=0.0
        call VLen(v,l2)
        omapB(j,i,1)=omapA(j,i,1)
        omapB(j,i,2)=omapA(j,i,2)
        omapB(j,i,3)=omapA(j,i,3)
        omapB(j,i,4)=omapA(j,i,4)+(TX0*(1.0-(l2/l1)))
        omapB(j,i,5)=omapA(j,i,5)+(TY0*(1.0-(l2/l1)))
        omapB(j,i,6)=omapA(j,i,6)
 110    continue
       enddo
      enddo
      omapA=omapB

      return
      end

c=====================================================================
      subroutine FixedGrid (omapA,xmax0,ymax0,step0)
c=====================================================================
      include 'common.fi'
      !(orientation map data, max x, max y, x coord corner, y coord corner, MaxDist, trans x, trans y)
      omapB=0.0
      do i=1,ymax0
       do j=1,xmax0
        omapB(j,i,4)=step0*(j-1)
        omapB(j,i,5)=step0*(i-1)
       enddo
      enddo

      do i=1,ymax0
       do j=1,xmax0
        ii=nint(omapA(j,i,4)/step0)+1
        jj=nint(omapA(j,i,5)/step0)+1
        !omapB(ii,jj,1)=omapA(j,i,1)
        !omapB(ii,jj,2)=omapA(j,i,2)
        !omapB(ii,jj,3)=omapA(j,i,3)
        !omapB(ii,jj,4)=((ii-1)*step0)
        !omapB(ii,jj,5)=((jj-1)*step0)
        !omapB(ii,jj,6)=omapA(j,i,6)
        omapB(ii,jj,1)=omapA(j,i,1)
        omapB(ii,jj,2)=omapA(j,i,2)
        omapB(ii,jj,3)=omapA(j,i,3)
        omapB(ii,jj,4)=((ii-1)*step0)
        omapB(ii,jj,5)=((jj-1)*step0)
        omapB(ii,jj,6)=omapA(j,i,6)
       enddo
      enddo
      omapA=omapB
      return
      end

c=====================================================================
      subroutine FillMissingPixels (omapA,xmax0,ymax0)
c=====================================================================
      include 'common.fi'
      omapB=omapA
      do i=2,ymax0-1 ! we don't need to check border pixels
         do j=2,xmax0-1
          if(omapA(j,i,6).eq.0.0) then
           !find number of non-zero
           NNNot0=0.0
           if(omapA(j-1,i,6).ne.0.0)NNNot0(1)=1
           if(omapA(j+1,i,6).ne.0.0)NNNot0(2)=1
           if(omapA(j,i-1,6).ne.0.0)NNNot0(3)=1
           if(omapA(j,i+1,6).ne.0.0)NNNot0(4)=1
           NNNot0Sum=NNNot0(1)+NNNot0(2)+NNNot0(3)+NNNot0(4)
           if(NNNot0Sum.ge.2) then
            if (NNNot0(2).eq.1.and.omapA(j+1,i,6).ne.0.0) then
             omapB(j,i,1)=omapA(j+1,i,1)
             omapB(j,i,2)=omapA(j+1,i,2)
             omapB(j,i,3)=omapA(j+1,i,3)
             omapB(j,i,6)=omapA(j+1,i,6)
             goto 120
            endif
            if (NNNot0(3).eq.1.and.omapA(j,i-1,6).ne.0.0) then
             omapB(j,i,1)=omapA(j,i-1,1)
             omapB(j,i,2)=omapA(j,i-1,2)
             omapB(j,i,3)=omapA(j,i-1,3)
             omapB(j,i,6)=omapA(j,i-1,6)
             goto 120
            endif
            if (NNNot0(1).eq.1.and.omapA(j-1,i,6).ne.0.0) then
             omapB(j,i,1)=omapA(j-1,i,1)
             omapB(j,i,2)=omapA(j-1,i,2)
             omapB(j,i,3)=omapA(j-1,i,3)
             omapB(j,i,6)=omapA(j-1,i,6)
             goto 120
            endif
            if (NNNot0(4).eq.1.and.omapA(j,i+1,6).ne.0.0) then
             omapB(j,i,1)=omapA(j,i+1,1)
             omapB(j,i,2)=omapA(j,i+1,2)
             omapB(j,i,3)=omapA(j,i+1,3)
             omapB(j,i,6)=omapA(j,i+1,6)
             goto 120
            endif

            endif
           else
            continue
           endif
 120       continue
          enddo ! closes j
         enddo ! closes i
         omapA=omapB
      return
      end

c=====================================================================
      subroutine WriteAng (omapA,xmax0,ymax0,fname0)
c=====================================================================
      include 'common.fi'

      open (50, file=trim(fname0), status='unknown')
      do j=1,ymax0
       do i=1,xmax0
        write(50,"(5f13.5,I8)") omapA(i,j,1),omapA(i,j,2),omapA(i,j,3),omapA(i,j,4),omapA(i,j,5),int(omapA(i,j,6))
       enddo
      enddo
      close(50)

      return
      end

c=====================================================================
      subroutine DisQ (n1, n2, theta)
c=====================================================================
      include 'common.fi'
      !n1, n2, are eulers, theta is disorientation angle

      call EToG(n1,g1nw)
      call EToG(n2,g2nw)
      call trans(g2nw,gtnw)        ! take the transpose of the second orientation
      call MToM(g1nw,gtnw,delgnw)      ! calculate the misorientation (dg)
      call GToE(delgnw,euler)            ! return the Euler angles for this misorientation
      !convert euler angles to a quaternion
      sigma=0.5*(euler(1)+euler(3))
      delta=0.5*(euler(1)-euler(3))
      Co=cos(euler(2)/2.0)
      Si=sin(euler(2)/2.0)
      q(1)=Co*cos(sigma)
      q(2)=Si*cos(delta)
      q(3)=Si*sin(delta)
      q(4)=Co*sin(sigma)

      !These operations follow the instructions in section 1.3.3.4 of Sutton and Balluffi
      !Make all four components positive.
      q=ABS(q)
      !new sort routine: Order so that q4≥q3≥q2≥q1
      do i=1,4
       k=0
       do j=1,4
        if(q(i).gt.q(j)) k=k+1
       enddo
       q0(k+1)=q(i)
      enddo

      !compare three choices and take the maximum:
      !e4
      !(e4+e3)/sqrt(2)
      !(e4+e1+e2+e3)/2
      choice(1)=q0(4)
      choice(2)=(q0(3)+q0(4))/1.414213
      choice(3)=(q0(1)+q0(2)+q0(3)+q0(4))/2.0
      maxchoice=MAXVAL(choice)

      !2cos^-1(largest of the three) = disorientation
      theta=2.0*acos2(maxchoice)
      return
      end
c=====================================================================
      subroutine GToAASing (g, v, temp1)
c=====================================================================
      include 'common.fi'
       !special G to axis-angle for case of 180° rotations.
       !based on information at:
       !https://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/index.htm
       !https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
       !by Martin John Baker
       if (g(1,1).gt.g(2,2).and.g(1,1).gt.g(3,3)) then
        Si=2.0*sqrt(1.0+g(1,1)-g(2,2)-g(3,3))
        q(1)=(g(3,2)-g(2,3))/Si
        q(2)=Si*0.25
        q(3)=(g(1,2)+g(2,1))/Si
        q(4)=(g(1,3)+g(3,1))/Si
       else
        if (g(2,2).gt.g(3,3)) then
         Si=2.0*sqrt(1.0+g(2,2)-g(1,1)-g(3,3))
         q(1)=(g(1,3)-g(3,1))/Si
         q(2)=(g(1,2)+g(2,1))/Si
         q(3)=Si*0.25
         q(4)=(g(2,3)+g(3,2))/Si
        else
         Si=2.0*sqrt(1.0+g(3,3)-g(1,1)-g(2,2))
         q(1)=(g(2,1)-g(1,2))/Si
         q(2)=(g(1,3)+g(3,1))/Si
         q(3)=(g(2,3)+g(3,2))/Si
         q(4)=Si*0.25
        endif
       endif
       
       temp1=2.0*acos2(q(1))
       v(1)=(q(2))/(sqrt(1-(q(1)*q(1))))
       v(2)=(q(3))/(sqrt(1-(q(1)*q(1))))
       v(3)=(q(4))/(sqrt(1-(q(1)*q(1))))

      return
      end


