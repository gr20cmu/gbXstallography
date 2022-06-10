C*****************************************
       program display_tjs
C*****************************************
       implicit none
	   include 'common.fi'

       write(6,*) '================================================='
	   write(6,*) ' PROGRAM TO DISPLAY TJ LINE SEGMENTS ON EBSD MAPS '
	   write(6,*) ' '
	   write(6,*) ' rohrer@cmu.edu'
	   write(6,*) ' version 06/10/2022'
       write(6,*) '================================================='
	   
	   version = 'version 06/10/2022'
       !This version is set for a maximum of 10,000 grains. If the image has more
       !than this number, increase the dimensions of 'grain' in common.
       !The pixels must be assigned grain average orientations
       !version 04/15/2022: fixed scaling of labels
       !version 06/10/2022: added variable mode capability

       !some constants that might be useful
       pi = 4.0*atan(1.0)

       !some statements used to format typed output
  1    format(A)
  14   format(6f7.3,4f9.3)
  15   format(2f8.2,A)
  17   format(A,2f9.3,A)
  18   format(3f7.3,A)
  19   format(A,A,A,I10,A)
  25   format(A,A,A)

       !The program requires a file, 'input.txt' that specifies all the parameters
       !Read the lower part of the input file for the definition of these parameters
       open(21, file='input.txt',status='old')
       read(21,*)keyword1    !.ang file base name
       read(21,*)mode        !0 for unsorted junctions, 1 for sorted.
       read(21,*)keyword2    !segment file base name of the output file
       read(21,*)first       !read the first file number
	   read(21,*)last        !read the last file number
       read(21,*)grid        !square or hex grid
       read(21,*)XSTEP       !length of step in x direction
       read(21,*)YSTEP       !length of step in y direction
       read(21,*)NROWS       !number of rows in the map
       read(21,*)NCOLS_ODD   !number of columns in odd numbered rows
       read(21,*)NCOLS_EVEN  !number of columns in even numbered rows
       read(21,*)XOff        !X offset for map
       read(21,*)YOff        !Y offset for map
       read(21,*)XSc         !X scale factor for the map
       read(21,*)YSc         !X scale factor for the map
       read(21,*)radian      !flag for units of Euler angles
       close(21)             !closes the file with input parameters
       !The next section just checks the input to make sure it is in bounds
       If (first.lt.0.OR.first.gt.999) then
        write(6,1)'line 3 of input is out of bounds, must be integer between 000 and 999'
        goto 9000
       endif
       If (last.lt.0.OR.first.gt.999) then
        write(6,1)'line 4 of input is out of bounds, must be integer between 000 and 999'
        goto 9000
       endif
       If (grid.ne.0.AND.grid.ne.1) then
        write(6,1)'line 5 of input is not a valid grid: use only 0 = HEX, 1 = SQR'
        goto 9000
       endif
       If (NROWS.gt.3000) then
        write(6,1)'Number of rows exceeds size limit of 3000: redimension map, edit this check,  and recompile'
        goto 9000
       endif
       If (NCOLS_ODD.gt.3000) then
        write(6,1)'Number of columns exceeds size limit of 3000: redimension map, edit this check, and recompile'
        goto 9000
       endif


       !------------------------------------------------------------------
       ! Begin the loop that sequentially analyzes all files in the list
       !------------------------------------------------------------------

       do nfile=first,last ! This is the loop over the different files (orientation maps)
        !Begin with the first file; first we have to create a string for the file name
        msf=int(nfile/100) !determine the most significant digit in the three digit file number
        nsf=int((nfile-(msf*100))/10) !determine the next (second) most significant digit
        lsf=nfile-(100*msf)-(10*nsf)  !determine the least (third) most significant digit
        one=char(48+msf)              !make it a string
        two=char(48+nsf)              !make it a string
        three=char(48+lsf)            !make it a string
        !create the filename by concatenating the base name, _, and the three digits of the file number.
        fname1=trim(keyword1)//'_'//one//two//three//'.ang'
        point = index(fname1,'.') - 1
        fname2=trim(keyword2)//'_'//one//two//three//'.txt'
        open (30, file=fname1,status='unknown') !This opens the first data file
        open (31, file=fname2,status='unknown')
        !The next section of code determines the length of the .ang files (nnline)
        !and the length of the header (header)
        nnline1 = 0  !counter
 100    continue
        read(30,*,end=101)
        nnline1=nnline1+1
        goto 100
 101    continue
        close(30)
        nnline2 = 0  !counter
 105    continue
        read(31,*,end=110)
        nnline2=nnline2+1
        goto 105
 110    continue
        close(31)
        !Check each line for the string '#'. If found, we know this line is part
        !of the header.  If not, we have reached the data section.
        open(30, file=fname1,status='unknown')
        header1 = 0
        hash1 = 0
        do i1=1,nnline1
         read(30, "(a)") inline
         hash1 = index(inline, '#')
         if (hash1.ne.0) then
          header1=header1+1
          hash1=0
         else
          goto 120
         endif
        enddo ! closes the i1=1,nnline loop
 120    close(30)
        open(31, file=fname2,status='unknown')
        header2 = 0
        hash2 = 0
        do i1=1,nnline2
         read(31, "(a)") inline
         hash2 = index(inline, '#')
         if (hash2.ne.0) then
          header2=header2+1
          hash2=0
         else
          goto 125
         endif
        enddo ! closes the i1=1,nnline loop
 125    close(31)

        !Here, inform the user of what we found and which file we are working on
        write(6,19)'File ',trim(fname1),' has ',nnline1,' lines.'
        write(6,19)'File ',trim(fname2),' has ',nnline2,' lines.'

        !initialize the map matrix
        do i1=1,3000
         do i2=1,3000
          do i3=1,3
           map(i1,i2,i3)=0.0
          enddo
         enddo
        enddo

  
       !------------------------------------------------------------------
       ! The next block of code finds grains, assigns the grain ID numbers,
       ! and determined the location of their 2D centroids
       !------------------------------------------------------------------

        write(6,1)'finding grains ...'
        open(30, file=fname1,status='unknown')
        do i1=1,header1  !Read the header lines
         read(30, "(a)") inline
        enddo
        !here, use the first pixel to define the first grain in the grain list
        grains = 1
        read(30,*) e1(1),e1(2),e1(3),xc,yc
        if (radian.eq.0) then
         e2=e1
         call DToRad3 (e2, e1)
        endif
        map(1,1,1)=float(grains)
        map(1,1,2)=xc
        map(1,1,3)=yc
        grain(1,1)=float(grains)
        grain(2,1)=e1(1)
        grain(3,1)=e1(2)
        grain(4,1)=e1(3)
        close(30)
        !now we need to start at the beginning of the file
        open(30, file=fname1,status='unknown')
        do i1=1,header1  !Read the header lines
         read(30, "(a)") inline
        enddo
        !now we go through all of the other pixels and load to a 2d matrix
        do i1=1,NROWS
         do i2=1,NCOLS_EVEN
          read(30,*) e1(1),e1(2),e1(3),xc,yc
          if (radian.eq.0) then
           e2=e1
           call DToRad3 (e2, e1)
          endif
          match=0
          do i3=1,grains !check the Euler angles from this line against those in the list
           !to handle data without grain averaged Euler angles, we would need to set a
           !threshold, exclude zero euler angles, and low CI data that is not orientation averaged.
           if (e1(1).eq.grain(2,i3).and.e1(2).eq.grain(3,i3).and.e1(3).eq.grain(4,i3)) then
            match=match+1 !the new Euler angles match one in the list, so we can stop
            ind=grain(1,i3)!this is the index of the grain it matched
            goto 500
           else
            continue
           endif
          enddo !ends the i3=1,grains loop
 500      continue
          if (match.eq.0) then !if we get here and match=0, it is a new grain
           if (grains.ge.9999) then ! conditional to catch errors
            write(6,1)'there are too many grains (limit is 10000): redimension grain matrix, edit this check,  and recompile'
            goto 9000
           endif
           grains=grains+1 !found a new grain
           grain(1,grains)=float(grains) !add it to the grain list
           grain(2,grains)=e1(1)
           grain(3,grains)=e1(2)
           grain(4,grains)=e1(3)
           ind=grains
          endif
          map(i1,i2,1)=float(ind) !all points have to be added to the map
          map(i1,i2,2)=xc
          map(i1,i2,3)=yc
         enddo ! this closes the i2=1,NCOLS_EVEN loop
         !if it is an odd column, there is one extra value for the hex grid
         if (mod(i1,2).ne.0.AND.grid.eq.0) then !odd row, read one more value!
          read(30,*) e1(1),e1(2),e1(3),xc,yc
         endif ! This value is not added to the map (ignored)
        enddo  ! i1=1,NROWS
        close(30)


       !------------------------------------------------------------------
       ! The next block of code reads the segment data
       !------------------------------------------------------------------
        !Read the segment data
        open(31, file=fname2,status='unknown')
        do i1=1,header2  !Read the header lines
         read(31, "(a)") inline
        enddo

        segments=nnline2-header2
        do i1=1,segments
         read(31,*)segment(i1,1),segment(i1,2),segment(i1,3),segment(i1,4),segment(i1,5),segment(i1,6),segment(i1,7),segment(i1,8),segment(i1,9),segment(i1,10)
        enddo
        close(31)

     








        !------------------------------------------------------------------
        ! The rest of this program writes out the data that was requested
        ! by the input file.
        !------------------------------------------------------------------


  
        !write IPF map with nodes and lines to a postscript file

        open (59, file='map_'//fname1(1:point)//'.ps',status='unknown') ! open the postscript file
        write (59,1)'%!PS' ! this is the first line of every postscript file
        n_up(1)=0.0 ! reference direction for the IPF map
        n_up(2)=0.0
        n_up(3)=1.0
        do i1=1,NROWS
         do i2=1,NCOLS_EVEN
          ind=int(map(i1,i2,1))
          e1(1)=grain(2,ind)
          e1(2)=grain(3,ind)
          e1(3)=grain(4,ind)
          if(e1(1).eq.0.0.and.e1(2).eq.0.0) then
           goto 3000
          endif

          call EToG (e1, g_or)
          call MToV (g_or,n_up,n_or) !transform it to the crystal reference frame
          do i3=1,3                  !make the components positive definite
           n_or(i3)=abs(n_or(i3))
          enddo
          if (n_or(2).gt.n_or(3)) then !order components so 3>2>1
           tem=n_or(2)
           n_or(2)=n_or(3)
           n_or(3)=tem
          endif
          if (n_or(1).gt.n_or(2)) then
           tem=n_or(1)
           n_or(1)=n_or(2)
           n_or(2)=tem
          endif
          if (n_or(2).gt.n_or(3)) then
           tem=n_or(2)
           n_or(2)=n_or(3)
           n_or(3)=tem
          endif
          eta = (180.0/pi)*atan2(n_or(2),n_or(1)) !projected angle from [101]
          alpha=(180.0/pi)*acos2(n_or(3)) !angle from [001]
          redr=alpha/54.73561 !ratio of angle from max value
          if (redr.gt.1.0) then
           redr=1.0
          endif
          red=1.0-redr !red channel level
          green=abs(eta-45.0)/45.0 ! green channel
          blue=(1.0-green)*redr  !blue channel
          green=green*redr
          red=red**2
          green=green**2
          blue=blue**2
          max_color=red !normalize
          if (green.gt.max_color) then
           max_color=green
          endif
          if (blue.gt.max_color) then
           max_color=blue
          endif
          red=red/max_color
          green=green/max_color
          blue=blue/max_color
          write (59,1)'1 setlinecap'
          write (59,1)'1 setlinewidth'
          write (59,18)red,green,blue,' setrgbcolor'
          write(59,17)'newpath ',XOff+(XSc*map(i1,i2,2)),YOff+(YSc*map(i1,i2,3)),' moveto 0 0 rlineto stroke'
 3000     continue
         enddo ! closes the i2=1,NCOLS_EVEN loop
        enddo ! closes the i1=1,NROWS loop


        !Draw all the segments on the map
        if (mode.eq.0) then
        !for a list of triples
        do i1=1,segments
         write(59,15)XOff+(XSc*segment(i1,7)),YOff+(YSc*segment(i1,8)),' newpath moveto'
         write(59,15)XOff+(XSc*segment(i1,9)),YOff+(YSc*segment(i1,10)),' lineto'
         write(59,1)'0 setgray'
         write(59,1)'0 1 0 setrgbcolor'
         write(59,1)'stroke'
        enddo

        else

        !This part is for sorted matches
        write (59,1)'1.0 setlinewidth '
        do i1=1,segments/6
         do i2=1,6
          ind=6*(i1-1)+i2
          if (i2.eq.1) then !change for first (1) or second (4) three
           write(59,15)XOff+(XSc*segment(ind,7)),YOff+(YSc*segment(ind,8)),' newpath moveto'
           write(59,15)XOff+(XSc*segment(ind,9)),YOff+(YSc*segment(ind,10)),' lineto'
           write(59,1)'1 0 0 setrgbcolor'
           write(59,1)'stroke'
          endif
          if (i2.eq.2) then!change for first (2) or second (5) three
           write(59,15)XOff+(XSc*segment(ind,7)),YOff+(YSc*segment(ind,8)),' newpath moveto'
           write(59,15)XOff+(XSc*segment(ind,9)),YOff+(YSc*segment(ind,10)),' lineto'
           write(59,1)'0 1 0 setrgbcolor'
           write(59,1)'stroke'
          endif
          if (i2.eq.3) then!change for first (3) or second (6) three
           write(59,15)XOff+(XSc*segment(ind,7)),YOff+(YSc*segment(ind,8)),' newpath moveto'
           write(59,15)XOff+(XSc*segment(ind,9)),YOff+(YSc*segment(ind,10)),' lineto'
           write(59,1)'0 0 1 setrgbcolor'
           write(59,1)'stroke'
          endif

          if (i2.eq.4) then !change for first (1) or second (4) three
           write(59,15)XOff+(XSc*segment(ind,7)),YOff+(YSc*segment(ind,8)),' newpath moveto'
           write(59,15)XOff+(XSc*segment(ind,9)),YOff+(YSc*segment(ind,10)),' lineto'
           write(59,1)'1 0 0 setrgbcolor'
           write(59,1)'stroke'
          endif
          if (i2.eq.5) then!change for first (2) or second (5) three
           write(59,15)XOff+(XSc*segment(ind,7)),YOff+(YSc*segment(ind,8)),' newpath moveto'
           write(59,15)XOff+(XSc*segment(ind,9)),YOff+(YSc*segment(ind,10)),' lineto'
           write(59,1)'0 1 0 setrgbcolor'
           write(59,1)'stroke'
          endif
          if (i2.eq.6) then!change for first (3) or second (6) three
           write(59,15)XOff+(XSc*segment(ind,7)),YOff+(YSc*segment(ind,8)),' newpath moveto'
           write(59,15)XOff+(XSc*segment(ind,9)),YOff+(YSc*segment(ind,10)),' lineto'
           write(59,1)'0 0 1 setrgbcolor'
           write(59,1)'stroke'
          endif
         enddo
        enddo

        write (59,1)'0.1 setlinewidth '
        write (59,1)'/Helvetica 12 selectfont'
        write(59,1)'1 0 0 setrgbcolor'
        write(59,15)XOff+(10.),YOff+(YSc*(NROWS*YSTEP))+32,' moveto'
        mystring='first'
        write(59,25)'(',mystring,') show'
        write(59,1)'0 1 0 setrgbcolor'
        write(59,15)XOff+(10.),YOff+(YSc*(NROWS*YSTEP))+17,' moveto'
        mystring='second'
        write(59,25)'(',mystring,') show'
        write(59,1)'0 0 1 setrgbcolor'
        write(59,15)XOff+(10.),YOff+(YSc*(NROWS*YSTEP))+2,' moveto'
        mystring='third'
        write(59,25)'(',mystring,') show'

        endif

        write(59,1)'showpage'


        close(59)


       enddo !closes the nfile=first,last loop

 9000  continue
  
       write(6,1) 'Program complete'

       end






	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
