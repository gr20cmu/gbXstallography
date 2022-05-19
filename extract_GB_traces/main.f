C*****************************************
       program extract_GB_traces
C*****************************************
       implicit none
	   include 'common.fi'

       write(6,1) '================================================='
	   write(6,1) ' PROGRAM TO EXTRACT GB TRACES FROM EBSD MAPS '
	   write(6,1) ' '
	   write(6,1) ' rohrer@cmu.edu'
	   write(6,1) ' version 04/12/2022'
       write(6,1) '================================================='
	   
	   version = 'version 04/12/2022'
       !This version is set for a maximum of 10,000 grains. If the image has more
       !than this number, increase the dimensions of 'grain' in common.
       !Similarly, there must be fewer than 30,000 trinodes and 100,000 boundary nodes
       !If the map is more than 3000 X 3000, increase the dimensions of 'map' in common
       !If a single pair of grain IDs is associated with more the 10 nodes, the boundary is
       !ignored.
       !In this version, the pixels must be assigned grain average orientations
       !Pixels with 0 0 0 Euler angles are assigned grain ID zero, and boundaries with
       !these fields are ignored.
       !04/12/22: updated symmetry.f (to correct trigonal) sub2.f (to improved DisG)

       !This version works, but will ignore a few relatively rare cases.
       ! - It ignores some very short segments on complex boundaries.  Such boundaries
       !   are irrelevant and probably erroneous because they approach the data resolution.
       ! - It misses some peninsular grains when the nodes at the base of the peninsula are
       !   separated by a distance a few pixels.
       ! - It misrepresented the shape of peninsular grains that have thicknesses
       !   that are less than or equal to two pixel spacings.
       ! - if a boundary has an odd number of nodes greater than three, it is skipped.  There
       !   are very few such boundaries.
       ! - If a boundary has a single node, it is skipped. There are very few such boundaries.
       !   To account for this case, find the gb node furthest from the single node, and make
       !   this a new node.  Then find the connecting lines.
       ! - A segment is only broken once based on a distance threshold, so if the maximum
       !   number of steps per segment is too large, the segment might not track the
       !   boundary as closely as desired. 

       ! Future work:
       ! Need to provide opportunity to shift reference frame of segments

       !some constants that might be useful
       pi = 4.0*atan(1.0)
       eps = 1.e-6
       t30 = 0.57735

       !some statements used to format typed output
  1    format(A)
  2    format(A,A)
  3    format(A,I8,A)
  4    format(5f8.5)
  5    format(I10,3f8.5,2f8.2)
  6    format(I10,2f10.5)
  7    format(2f10.4,3I8)
  8    format(4x,2I8)
  9    format(3I8)
  10   format(8x,f6.3,A,f6.3,4x,f6.3,A,f6.3)
  11   format(f6.3,A,f6.3,4x,f6.3,A,f6.3,4x,f6.3,A,f6.3)
  12   format(2f10.5,3I8)
  13   format(25I6)
  14   format(6f7.3,4f9.3)
  15   format(2f8.2,A)
  16   format(2I8,4f7.3,7I8)
  17   format(A,2f9.3,A)
  18   format(3f7.3,A)
  19   format(A,A,A,I10,A)
  20   format(A,I6,A,I6,A)
  21   format(A,I4)
  22   format(A,f11.6)
  23   format(2f10.5,2I8)
  24   format(I4)
  25   format(A,A,A)

       out=0

       !The program requires a file, 'input.txt' that specifies all the parameters
       !Read the lower part of the input file for the definition of these parameters
       open(21, file='input.txt',status='old')
       read(21,*)keyword     !file base name
       read(21,*)out_fname   !base name of the output file
       read(21,*)first       !read the first file number
	   read(21,*)last        !read the last file number
       read(21,*)source      !source of the orientation map
       read(21,*)grid        !square or hex grid
       read(21,*)XSTEP       !length of step in x direction
       read(21,*)YSTEP       !length of step in y direction
       read(21,*)NROWS       !number of rows in the map
       read(21,*)NCOLS_ODD   !number of columns in odd numbered rows
       read(21,*)NCOLS_EVEN  !number of columns in even numbered rows
       read(21,*)SegMax      !maximal line segment
       read(21,*)DistTol     !tolerance distance between the line segment and the GB
       read(21,*)radian      !flag for units of Euler angles
       read(21,1)inline2     !user added comment
       read(21,*)XOff        !X offset for map
       read(21,*)YOff        !Y offset for map
       read(21,*)XSc         !X scale factor for the map
       read(21,*)YSc         !X scale factor for the map
       read(21,*)out(1)      !output flag for the images
       read(21,*)out(2)      !output flag for the grain list
       read(21,*)out(3)      !output flag for the grain ID map
       read(21,*)out(4)      !output flag for the tj node list
       read(21,*)out(5)      !output flag for the gb node list
       read(21,*)out(6)      !output flag for the gb data list
       read(21,*)out(7)      !output flag for writing grain ID numbers on the map
       close(21)             !closes the file with input parameters
       !The next section just checks the input to make sure parameters are in bounds
       If (first.lt.0.OR.first.gt.999) then
        write(6,1)'line 3 of input is out of bounds, must be integer between 000 and 999'
        goto 9000
       endif
       If (last.lt.0.OR.first.gt.999) then
        write(6,1)'line 4 of input is out of bounds, must be integer between 000 and 999'
        goto 9000
       endif
       If (source.ne.0) then
        write(6,1)'line 5 of input is not a valid source: use only 0 = TSL'
        goto 9000
       endif
       If (grid.ne.0.AND.grid.ne.1) then
        write(6,1)'line 6 of input is not a valid grid: use only 0 = HEX, 1 = SQR'
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
       If (SegMax.lt.3) then
        write(6,1)'The number of segments on line 12 of the input file must be and integer greater than or equal to 3'
        goto 9000
       endif
       If (SegMax.gt.NROWS) then
        write(6,1)'The number of segments on line 12 of the input file must be less than the number of rows'
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
        fname=trim(keyword)//'_'//one//two//three//'.ang'
        point = index(fname,'.') - 1
        base2=trim(keyword)//'_'//one//two//three
        open (30, file=fname,status='unknown') !This opens the first data file
        !The next section of code determines the length of the .ang files (nnline)
        !and the length of the header (header)
        nnline = 0  !counter
 100    continue
        read(30,*,end=101)
        nnline=nnline+1
        goto 100
 101    continue
        close(30)
        !Check each line for the string '#'. If found, we know this line is part
        !of the header.  If not, we have reached the data section.
        open(30, file=fname,status='unknown')
        header = 0
        hash = 0
        do i1=1,nnline
         read(30, "(a)") inline
         hash = index(inline, '#')
         if (hash.ne.0) then
          header=header+1
          hash=0
         else
          goto 120
         endif
        enddo ! closes the i1=1,nnline loop
 120    close(30)
        search_fail=0
        !Here, inform the user of what we found and which file we are working on
        write(6,19)'File ',trim(fname),' has ',nnline,' lines.'
        !write(6,3)'There are',header,' lines in the header'
        !The next 44 lines simply initilize a bunch of the matrices used to store information
        !initialize the grain matrix
        do i1=1,10000
         do i2=1,6
          grain(i2,i1)=0.0
         enddo
        enddo
        !initialize the map matrix
        do i1=1,3000
         do i2=1,3000
          do i3=1,3
           map(i1,i2,i3)=0.0
          enddo
         enddo
        enddo
        !initialize the node matrix
        do i1=1,30000
         do i2=1,6
          node(i1,i2)=0.0
         enddo
        enddo
        !initialize the grain boundary matrix
        do i1=1,300000
         do i2=1,4
          gb(i1,i2)=0.0
         enddo
        enddo
        !initialize the grain boundary matrix
        do i1=1,10000
         do i2=1,300
          GBDat(i1,i2)=0.0
         enddo
        enddo
        !initialize the grain boundary trace matrix
        do i1=1,200
         do i2=1,2
          trace(i1,i2)=0.0
         enddo
        enddo
        !initialize the segment matrix
        do i1=1,10000
         do i2=1,10
          segment(i1,i2)=0.0
         enddo
        enddo
        c_seg=0.0
        dist=0.0
        NewNode=0.0

       !------------------------------------------------------------------
       ! The next block of code finds grains, assigns the grain ID numbers,
       ! and determined the location of their 2D centroids
       !------------------------------------------------------------------

        write(6,1)'finding grains ...'
        open(30, file=fname,status='unknown')
        do i1=1,header  !Read the header lines
         read(30, "(a)") inline
        enddo
        !here, use the first pixel to define the first grain in the grain list
        !modify so the grain list does not contain 0 0 0 grains
        !also need to not identify boundaries with 0 0 0 points
        !on map, label with ID zero


        !Begin by finding grain number 1 to seed the grain list
        do i1=1,NROWS
         do i2=1,NCOLS_EVEN
          read(30,*) e1(1),e1(2),e1(3),xc,yc
          if (radian.eq.0) then
           e2=e1
           call DToRad3 (e2, e1)
          endif
          if (e1(1).eq.0.0.AND.e1(2).eq.0.0.AND.e1(3).eq.0.0) then
           continue !part of grain 0, goto next line
          else
           grains = 1   !first non-zero orientation assigned to grain 1
           grain(1,1)=float(grains)
           grain(2,1)=e1(1)
           grain(3,1)=e1(2)
           grain(4,1)=e1(3)
           goto 200 !exit search loop
          endif
         enddo ! this closes the i2=1,NCOLS_EVEN loop
        enddo  ! i1=1,NROWS
 200    continue
        close(30)

        !now we need to start at the beginning of the file
        open(30, file=fname,status='unknown')
        do i1=1,header  !Read the header lines
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
          !first deal with special case of grain 0
          if (e1(1).eq.0.0.AND.e1(2).eq.0.0.AND.e1(3).eq.0.0) then
           map(i1,i2,1)=0.0 !assign grain ID zero
           map(i1,i2,2)=xc
           map(i1,i2,3)=yc
           goto 510 !skip to next line, nothing to see here
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
 510      continue
         enddo ! this closes the i2=1,NCOLS_EVEN loop
         !if it is an odd column, there is one extra value for the hex grid
         if (mod(i1,2).ne.0.AND.grid.eq.0) then !odd row, read one more value!
          read(30,*) e1(1),e1(2),e1(3),xc,yc
         endif ! This value is not added to the map (ignored)
        enddo  ! i1=1,NROWS
        close(30)

        !Here, find the grain centroid locations.  This is just for diagnostics,
        !making it possible to label the grains with their IDs on the map
        do i1=1,grains
         CenTotX=0
         CenTotY=0
         ct=0
         do i2=1,NROWS
          do i3=1,NCOLS_EVEN
           if (grain(1,i1).eq.map(i2,i3,1)) then
            !part of grain in consideration, so add coordinates, increment a counter
            ct=ct+1
            CenTotX=CenTotX+map(i2,i3,2)
            CenTotY=CenTotY+map(i2,i3,3)
           else
            !keep looking
            continue
           endif
          enddo
         enddo
         CentX=CenTotX/float(ct)!divide coordinates by counter to get centroid
         CentY=CenTotY/float(ct)
         grain(5,i1)=CentX!save centroid
         grain(6,i1)=CentY
        enddo ! end the i1=1,grains loop and goto next grain ID

        write(6,3)'There are :',grains,' grains'

        !------------------------------------------------------------------
        ! The next block of code finds all multi nodes (these are trinodes
        ! for hex grids and tri- and quadnodes for square grids), and then
        ! finds all grain boundary nodes (binodes).  This is different for
        ! the hex and square grids, so the code divides at this point, using
        ! a if, else, endif structure.
        !------------------------------------------------------------------

        if (grid.eq.0) then  !The is case for the hex grid
         !next, we have to identify all trinodes
         write(6,1)'finding tri-nodes ...'
         nodes=0
         !there is a special case for the first row, so we account for that first
         i1=1
         do i2=1,NCOLS_EVEN
          if (map(i1,i2,1).ne.map(i1,i2+1,1)) then
           if (map(i1,i2+1,1).eq.0.0) then
            goto 580
           endif
           nodes=nodes+1
           node(nodes,1)=map(i1,i2,2)+(XSTEP/2.0)
           node(nodes,2)=map(i1,i2,3)-(t30*(XSTEP/2.0))
           node(nodes,3)=map(i1,i2,1)
           node(nodes,4)=map(i1,i2+1,1)
           node(nodes,5)=0.0
          endif
 580      continue
         enddo ! closes the i2=1,NCOLS_EVEN loop

         !there is another special case for column one, so we account for that here
         i2=1
         do i1=2,NROWS
          if (map(i1,i2,1).ne.map(i1+1,i2,1)) then !possible trijuction
           if (map(i1+1,i2,1).eq.0.0) then
            goto 590
           endif
           !this is an edge junction - record and move on
           nodes=nodes+1
           node(nodes,1)=map(i1,i2,2)
           node(nodes,2)=map(i1,i2,3)+(t30*XSTEP)
           node(nodes,3)=map(i1,i2,1)
           node(nodes,4)=map(i1+1,i2,1)
           node(nodes,5)=0
           goto 590
          endif
 590      continue
         enddo ! closes the i1=2,NROWS loop

         !This next section handles all the other cases
         do i1=2,NROWS
          do i2=2,NCOLS_EVEN
           !Check the trinode at x+xstep/2,y+tan30*Xstep/2
           !surrounded by grain IDs i1,i2 to i1,i2+1 and i1+1,i2
           if (map(i1,i2,1).eq.map(i1,i2+1,1)) then !not a trinode
            goto 600
           else
            if (map(i1,i2,1).eq.map(i1+1,i2,1)) then !not a trinode
             goto 600
            else
             if (map(i1,i2+1,1).eq.map(i1+1,i2,1)) then !not a trinode
              goto 600
             else ! it is a trinode
              if (nodes.ge.29999) then
               write(6,1) 'Too many nodes: redimension node matrix, edit this check, recompile'
               goto 9000
              endif
              nodes=nodes+1
              node(nodes,1)=map(i1,i2,2)+(XSTEP/2.0)
              node(nodes,2)=map(i1,i2,3)+(t30*(XSTEP/2.0))
              node(nodes,3)=map(i1,i2,1)
              node(nodes,4)=map(i1,i2+1,1)
              node(nodes,5)=map(i1+1,i2,1)
             endif
            endif
           endif
 600       continue
           !Check the trinode at x,y+tan30*Xstep
           !surrounded by grain IDs i1,i2 to i1+1,i2-1 (check zero), i1+1,i2
           if (map(i1,i2,1).eq.map(i1+1,i2,1)) then !not a trinode
            goto 610
           else
            if (map(i1,i2,1).eq.map(i1+1,i2-1,1)) then !not a trinode
             goto 610
            else
             if (map(i1+1,i2,1).eq.map(i1+1,i2-1,1)) then !not a trinode
              goto 610
             else ! it is a trinode
              if (nodes.ge.29999) then
               write(6,1) 'Too many nodes: redimension node matrix, edit this check, recompile'
               goto 9000
              endif
              nodes=nodes+1
              node(nodes,1)=map(i1,i2,2)
              node(nodes,2)=map(i1,i2,3)+(t30*XSTEP)
              node(nodes,3)=map(i1,i2,1)
              node(nodes,4)=map(i1+1,i2,1)
              node(nodes,5)=map(i1+1,i2-1,1)
             endif
            endif
           endif
 610       continue
          enddo ! This ends the i2=1,NCOLS_EVEN
         enddo ! This ends the i1=1,NROWS loop

         write(6,3)'There are :',nodes,' nodes'

         !Identify all boundaries: at each pixel (i1,i2), check two gbs
         !with indices i1-1,i2 and i1,i2+1
         !exclude any with GID = 0
         write(6,1)'finding grain boundary nodes ...'
         gbs=0
         !First, take care of the special case in the first row
         i1=1
         do i2=1,NCOLS_EVEN
          if (map(i1,i2,1).eq.map(i1,i2+1,1)) then !check the neighbor
           goto 630 !no boundary
          else
           if (map(i1,i2+1,1).eq.0.0) then !check for limit of domain
            goto 630 ! no boundary
           else
            continue
           endif
           gbs=gbs+1
           gb(gbs,1)=map(i1,i2,2)+(XSTEP/2.0)
           gb(gbs,2)=map(i1,i2,3)
           gb(gbs,3)=map(i1,i2,1)
           gb(gbs,4)=map(i1,i2+1,1)
          endif
 630      continue
         enddo ! closes the i2=1,NCOLS_EVEN loop
         ! now look at all the rest
         do i1=2,NROWS
          do i2=1,NCOLS_EVEN
           if (map(i1,i2,1).eq.map(i1-1,i2,1)) then !check first one
            goto 640 !no boundary
           else !boundary
            if (map(i1-1,i2,1).eq.0.0) then !check for limit of domain
             goto 640
            else
             continue
            endif
            if (gbs.ge.299999) then
             write(6,1) 'Too many gb nodes: redimension gb matrix, edit this check, recompile'
             goto 9000
            endif
            gbs=gbs+1
            gb(gbs,1)=map(i1,i2,2)-(XSTEP/4.0)
            gb(gbs,2)=map(i1,i2,3)-(t30*0.75*XSTEP)
            gb(gbs,3)=map(i1,i2,1)
            gb(gbs,4)=map(i1-1,i2,1)
           endif
 640       continue

           if (map(i1,i2,1).eq.map(i1,i2+1,1)) then
            goto 700 !no boundary
           else !boundary
            if (map(i1,i2+1,1).eq.0.0) then !check for limit of domain
             goto 700
            else
             continue
            endif
            if (gbs.ge.299999) then
             write(6,1) 'Too many gb nodes: redimension gb matrix, edit this check, recompile'
             goto 9000
            endif
            gbs=gbs+1
            gb(gbs,1)=map(i1,i2,2)+(XSTEP/2.0)
            gb(gbs,2)=map(i1,i2,3)
            gb(gbs,3)=map(i1,i2,1)
            gb(gbs,4)=map(i1,i2+1,1)
           endif
 700      continue
          enddo ! This ends the i2=1,NCOLS_EVEN
         enddo ! This ends the i1=1,NROWS loop

         write(6,3)'There are :',gbs,' grain boundary nodes'

         !This completes the node search for data on a hex grid
         !The next section of code does the same thing for the square grid

        else          !assuming the data is on a square grid
         write(6,1)'finding tri-nodes ...'
         nodes=0
         !there is a special case for the first row, so we account for that first
         i1=1
         do i2=1,NCOLS_EVEN
          if (map(i1,i2,1).ne.map(i1,i2+1,1)) then !they have different GIDs, so trinode
           if (map(i1,i2+1,1).eq.0.0) then
            goto 705 !the means we are at the end of the domain, so not a trinode
           endif
           nodes=nodes+1
           node(nodes,1)=map(i1,i2,2)+(XSTEP) !x-coord of node
           node(nodes,2)=map(i1,i2,3) !y-coord of node
           node(nodes,3)=map(i1,i2,1) ! first grain ID
           node(nodes,4)=map(i1,i2+1,1) ! second grain ID
           node(nodes,5)=0.0 ! third grain ID
          endif
 705      continue
         enddo ! closes the i2=1,NCOLS_EVEN loop

         !there is another special case for column one, so we account for that here
         i2=1
         do i1=2,NROWS
          if (map(i1,i2,1).ne.map(i1+1,i2,1)) then !possible trijuction
           if (map(i1+1,i2,1).eq.0.0) then
            goto 710
           endif
           !this is an edge junction - record and move on
           nodes=nodes+1
           node(nodes,1)=map(i1,i2,2) !x-coord of node
           node(nodes,2)=map(i1,i2,3)+(YSTEP) !y-coord of node
           node(nodes,3)=map(i1,i2,1)
           node(nodes,4)=map(i1+1,i2,1)
           node(nodes,5)=0
          endif
 710      continue
         enddo ! closes the i1=2,NROWS loop

         do i1=2,NROWS
          do i2=2,NCOLS_EVEN
          position=0
          new_or=0
          if (map(i1,i2,1).ne.map(i1+1,i2,1)) then
           new_or=new_or+1
           position(1)=1
          endif
          if (map(i1,i2,1).ne.map(i1,i2+1,1)) then
           new_or=new_or+1
           position(2)=1
          endif
          if (map(i1,i2,1).ne.map(i1+1,i2+1,1)) then
           new_or=new_or+1
          position(3)=1
          endif
          if (map(i1+1,i2,1).eq.map(i1,i2+1,1)) then
           new_or=new_or-1
           position(2)=0
          endif
          if (map(i1+1,i2,1).eq.map(i1+1,i2+1,1)) then
           new_or=new_or-1
           position(3)=0
          endif
          if (map(i1+1,i2+1,1).eq.map(i1,i2+1,1)) then
           new_or=new_or-1
           position(3)=0
          endif
          if (new_or.ge.2) then ! this is a trinode or a quad node
           if (nodes.ge.29999) then
            write(6,1) 'Too many nodes: redimension node matrix, edit this check, recompile'
            goto 9000
           endif
           if (new_or.eq.3) then ! this is for a quad node
            nodes=nodes+1
            node(nodes,1)=map(i1,i2,2)+(XSTEP/2.0)
            node(nodes,2)=map(i1,i2,3)+(YSTEP/2.0)
            node(nodes,3)=map(i1,i2,1)
            node(nodes,4)=map(i1+1,i2,1)
            node(nodes,5)=map(i1,i2+1,1)
            node(nodes,6)=map(i1+1,i2+1,1)
           else                  ! this is for a tri node
            nodes=nodes+1
            node(nodes,1)=map(i1,i2,2)+(XSTEP/2.0)
            node(nodes,2)=map(i1,i2,3)+(YSTEP/2.0)
            node(nodes,3)=map(i1,i2,1)
            if (position(1).eq.1.AND.position(2).eq.1) then
             node(nodes,4)=map(i1+1,i2,1)
             node(nodes,5)=map(i1,i2+1,1)
             goto 712
            endif
            if (position(1).eq.1.AND.position(3).eq.1) then
             node(nodes,4)=map(i1+1,i2,1)
             node(nodes,5)=map(i1+1,i2+1,1)
             goto 712
            endif
            if (position(2).eq.1.AND.position(3).eq.1) then
             node(nodes,4)=map(i1,i2+1,1)
             node(nodes,5)=map(i1+1,i2+1,1)
             goto 712
            endif
            write(6,*)'no node assigned at ',i1,i2
           endif
          endif
 712      continue
         enddo ! This ends the i2=1,NCOLS_EVEN
        enddo ! This ends the i1=1,NROWS loop

        write(6,3)'There are :',nodes,' nodes'
         

        !Identify all boundaries: at each pixel (i1,i2), check two gbs
        !with indices i1+1,i2 and i1,i2+1
        !exclude any with GID = 0
        write(6,1)'finding grain boundary nodes ...'
        gbs=0

        do i1=1,NROWS
         do i2=1,NCOLS_EVEN
          if (map(i1,i2,1).eq.map(i1+1,i2,1)) then !check first one
           goto 715 !no boundary
          else !boundary
           if (map(i1+1,i2,1).eq.0.0) then !check for limit of domain
            goto 715
           else
            continue
           endif
           if (gbs.ge.299999) then
            write(6,1) 'Too many gb nodes: redimension gb matrix, edit this check, recompile'
            goto 9000
           endif
           gbs=gbs+1
           gb(gbs,1)=map(i1,i2,2)
           gb(gbs,2)=map(i1,i2,3)+(YSTEP/2.0)
           gb(gbs,3)=map(i1,i2,1)
           gb(gbs,4)=map(i1+1,i2,1)
          endif
 715      continue

          if (map(i1,i2,1).eq.map(i1,i2+1,1)) then
           goto 720 !no boundary
          else !boundary
           if (map(i1,i2+1,1).eq.0.0) then !check for limit of domain
            goto 720
           else
            continue
           endif
           if (gbs.ge.299999) then
            write(6,1) 'Too many gb nodes: redimension gb matrix, edit this check, recompile'
            goto 9000
           endif
           gbs=gbs+1
           gb(gbs,1)=map(i1,i2,2)+(XSTEP/2)
           gb(gbs,2)=map(i1,i2,3)
           gb(gbs,3)=map(i1,i2,1)
           gb(gbs,4)=map(i1,i2+1,1)
          endif
 720     continue
         enddo ! This ends the i2=1,NCOLS_EVEN
        enddo ! This ends the i1=1,NROWS loop

        write(6,3)'There are :',gbs,' grain boundary nodes'

        endif  ! This closes the section of code for finding the nodes

        !------------------------------------------------------------------
        ! The next block of code gathers together all of the node information,
        ! building a table for all grain boundaries that have a binode.  The
        ! table includes the two grain ids, the number of multinodes and
        ! binodes, then a list of up to 10 multinodes and 200 binodes.
        ! This information is needed to draw in the segments.
        !------------------------------------------------------------------

        !Create a table, grain ID 1, grain ID2, number of GB nodes, number of tri nodes,
        !up to 10 trinode IDs, all GB nodes with same IDs
        write(6,1)'finding related boundary nodes ...'
        !This section of code orders the grain IDs so the smaller one is first
        do i1=1,gbs
         if (int(gb(i1,3)).gt.int(gb(i1,4))) then
          tem=gb(i1,4)
          gb(i1,4)=gb(i1,3)
          gb(i1,3)=tem
         else
          continue
         endif
        enddo ! closes the i1=1,gbs loop
        !put an initial non-zero entry into the GBDat array
        do i1=1,gbs
         if (gb(i1,3).eq.0.0) then
          continue
         else
          GBDat(1,1)=gb(i1,3)
          GBDat(1,2)=gb(i1,4)
          GBDat(1,3)=1
          GBDat(1,15)=1
          boundaries=1
          goto 730
         endif
        enddo
 730    continue

        do i1=2,gbs
         if (gb(i1,3).eq.0.0) then
          goto 740 !ignore boudnaries with grain 0
         endif
         match=0
         do i2=1,boundaries
          if (GBDat(i2,1).eq.0.0) then
           goto 735 !ignore boudnaries with grain 0
          endif
          if ((gb(i1,3).eq.GBDat(i2,1)).and.(gb(i1,4).eq.GBDat(i2,2))) then
           match=match+1
           !put this node in table row where it matched
           GBDat(i2,3)=GBDat(i2,3)+1
           GBDat(i2,14+int(GBDat(i2,3)))=i1
          else
           continue !keep looking
          endif
 735      continue
         enddo ! this closes the i2=1,boundaries loop
         !if you get here and match = 0, then it is a new GI pair, so make a new entry in the table
         if (match.eq.0) then
          boundaries=boundaries+1
          GBDat(boundaries,1) = gb(i1,3)
          GBDat(boundaries,2) = gb(i1,4)
          GBDat(boundaries,3) = 1
          GBDat(boundaries,15) = i1
         else
          continue
         endif
 740     continue
        enddo ! this ends the i1=2,gbs loop
        write(6,3)'There are :',boundaries,' grain boundaries'

        !Next, find the nodes at the end of each boundary
        NodeMatch=0
        do i1=1,boundaries
         do i2=1,nodes
          match=0
          if (int(node(i2,3)).eq.int(GBDat(i1,1))) then
           match=match+1
           goto 800
          endif
          if (int(node(i2,3)).eq.int(GBDat(i1,2))) then
           match=match+1
          endif
 800      continue
          if (int(node(i2,4)).eq.int(GBDat(i1,1))) then
           match=match+1
           goto 810
          endif
          if (int(node(i2,4)).eq.int(GBDat(i1,2))) then
           match=match+1
          endif
 810      continue
          if (int(node(i2,5)).eq.int(GBDat(i1,1))) then
           match=match+1
           goto 815
          endif
          if (int(node(i2,5)).eq.int(GBDat(i1,2))) then
           match=match+1
          endif
          !if (int(node(i2,6)).ne.0) then
 815      continue
           if (int(node(i2,6)).eq.int(GBDat(i1,1))) then
            match=match+1
            goto 820
           endif
           if (int(node(i2,6)).eq.int(GBDat(i1,2))) then
            match=match+1
           endif
          !endif

 820      continue
          if (match.lt.2) then
           goto 825
          else
           NodeMatch=NodeMatch+1
           GBDat(i1,4)=GBDat(i1,4)+1
           GBDat(i1,4+NodeMatch)=i2
          endif
 825      continue
         enddo  ! close the i2=1,nodes loop
 830     continue
         if(NodeMatch.gt.10) then
          write(6,*)'WARNING: boundary ',i1,' associated with more than 10 nodes!'
         endif
         NodeMatch=0
        enddo ! close the i1=1,boundaries loop

        !------------------------------------------------------------------
        ! The next block of code is for writing the grain boundary line
        ! segments.  It creates segments from multinode to multinode,
        ! along the path of binodes.  Several different geometries are
        ! accounted for in different sections.
        !------------------------------------------------------------------

        write(6,1)'finding grain boundary trace segments ...'
        segments=0
        NewNodes=0
        !This is the most common case, when a single boundary type has an even number of nodes.
        do i1=1,boundaries ! loop over all identified boundaries
         trace=0.0
         c_seg=0.0
         NodeList=0.0
         search=0 ! keep track of attempts to find new boundaries
         if (int(GBDat(i1,4)).lt.2) then !ignore boundaries that do not have 2 or more TJ nodes
          goto 1000
         endif
         if (mod(int(GBDat(i1,4)),2).ne.0) then !ignore boundaries that have an odd number of TJ nodes
          goto 1000
         endif
         if (int(GBDat(i1,4)).gt.10) then !ignore boundaries with more than 10 TJ nodes
          goto 1000
         endif
         !condition to eliminate the very short 1 to 2 pixel long gbs on GBs with ≥ 4 nodes
         if ((GBDat(i1,3)/GBDat(i1,4)).lt.1.5.and.GBDat(i1,4).gt.2.0) then
          goto 1000
         endif
         !Start at a node, find the closest GB node, make a line.  At every step, check
         !if you are close to one of the other nodes of the same boundary.  If you are close
         !to a node, draw a line and start search at a new node.  If not, find the
         !next cloest GB node and continue until SegMax.
         !start from the final node and proceed as before.
         TJNum=int(GBDat(i1,4))  !Number of TJs associated with this GB
         do i2=1,TJNum  ! load the node IDs to the node list
          NodeList(i2,1)=GBDat(i1,4+i2)
          NodeList(i2,2)=1.0
         enddo ! this ends the i2=1,TJNum loop
         TraceLength=int(GBDat(i1,3))  !Number of gb nodes
         do i2=1,TraceLength  ! load the segment ids to the vector trace
          trace(i2,1)=GBDat(i1,14+i2)
          trace(i2,2)=1.0
         enddo ! this ends the i2=1,TraceLength loop
         do i2=1,TJNum
          if (NodeList(i2,2).eq.0.0) then  !if 0, this node has already been used, skip it
           goto 910
          endif
          NodeList(i2,2)=0 ! We are using this node, so mark it as used.
          StartPointX=node(int(GBDat(i1,4+i2)),1) !the start coordinate is the coordinate of the first node
          StartPointY=node(int(GBDat(i1,4+i2)),2)
          LastPointX=StartPointX
          LastPointY=StartPointY
 901      continue !after writing a line of SegMax, redefine the start point
          StartPointX=LastPointX
          StartPointY=LastPointY
          search=search+1
          if (search.gt.199) then
           search_fail=search_fail+1
           goto 1000
          endif
          segs=0
 902      continue
          dist=100
          do i3=1,TraceLength ! we find the closest GB node to this tridnode
           if (trace(i3,2).eq.0.0) then  !if 0, this GB node has already been used
            goto 905  !goto next gb node
           endif
           !create a vector from the start coordinate to the gb node of interest
           vec(1)=LastPointX-gb(int(trace(i3,1)),1)
           vec(2)=LastPointY-gb(int(trace(i3,1)),2)
           vec(3)=0
           call VLen(vec,length) ! get length of vector
           if (length.lt.dist) then ! if the vector is the shortest one, remember it.
            dist=length
            ind=i3
           endif
 905       continue
          enddo  !ends the i3=1,TraceLength loop
          segs=segs+1
          !Now we know the gb node closest to the start/last point
          trace(ind,2)=0.0 ! mark as used
          LastPointX=gb(int(trace(ind,1)),1)
          LastPointY=gb(int(trace(ind,1)),2)
          c_seg(segs,1)=LastPointX
          c_seg(segs,2)=LastPointY
          !Now we check if last point is near one of the unused nodes
          do i4=1,TJNum
           if (NodeList(i4,2).eq.0.0) then  !if 0, this node has already been used, skip it
            goto 906
           endif
           !Find distance from last point to this node
           vec(1)=LastPointX-node(int(GBDat(i1,4+i4)),1)
           vec(2)=LastPointY-node(int(GBDat(i1,4+i4)),2)
           vec(3)=0.0
           call VLen(vec,length) ! get length of vector
           if (length.gt.2.0*XSTEP) then
            goto 906 !check another node
           else
            !Draw a line from the start point to this node
            ep1x=StartPointX
            ep1y=StartPointY
            ep2x=node(int(GBDat(i1,4+i4)),1)
            ep2y=node(int(GBDat(i1,4+i4)),2)
            !Check the distance to boundary and break segent if it is too large
            do i5=1,segs
             px=c_seg(i5,1)
             py=c_seg(i5,2)
             dist=abs((ep2x-ep1x)*(ep1y-py)-(ep1x-px)*(ep2y-ep1y))/(((ep2x-ep1x)**2)+((ep2y-ep1y)**2))**0.5
             if (dist.gt.DistTol*XSTEP) then
              !write(6,*)'distance too large',dist
              !write(6,*)px,py
              !If dist too large, make two segs from ep1x,ep1y to px,py and from px,py to ep2x,ep2y
              NewNodes=NewNodes+1
              NewNode(NewNodes,1)=px
              NewNode(NewNodes,2)=py
              segments=segments+1
              segment(segments,1)=grain(2,int(GBDat(i1,1)))
              segment(segments,2)=grain(3,int(GBDat(i1,1)))
              segment(segments,3)=grain(4,int(GBDat(i1,1)))
              segment(segments,4)=grain(2,int(GBDat(i1,2)))
              segment(segments,5)=grain(3,int(GBDat(i1,2)))
              segment(segments,6)=grain(4,int(GBDat(i1,2)))
              segment(segments,7)=ep1x
              segment(segments,8)=ep1y
              segment(segments,9)=px
              segment(segments,10)=py
              segments=segments+1
              segment(segments,1)=grain(2,int(GBDat(i1,1)))
              segment(segments,2)=grain(3,int(GBDat(i1,1)))
              segment(segments,3)=grain(4,int(GBDat(i1,1)))
              segment(segments,4)=grain(2,int(GBDat(i1,2)))
              segment(segments,5)=grain(3,int(GBDat(i1,2)))
              segment(segments,6)=grain(4,int(GBDat(i1,2)))
              segment(segments,7)=px
              segment(segments,8)=py
              segment(segments,9)=ep2x
              segment(segments,10)=ep2y
              segs=0
              NodeList(i4,2)=0.0   !mark the TJ as used
              goto 910
             endif
            enddo  ! ends the i5=1,SegMax loop
            !Here if it passes the distance check, write the segment
            vec(1)=abs(ep1x-ep2x)
            vec(2)=abs(ep1y-ep2y)
            vec(3)=0.0
            call VLen(vec,length) ! get length of vector
            if (length.eq.0.0) then !if this happens, something is wrong, skip it.
             search_fail=search_fail+1
             goto 1000
            endif
            segments=segments+1
            segment(segments,1)=grain(2,int(GBDat(i1,1)))
            segment(segments,2)=grain(3,int(GBDat(i1,1)))
            segment(segments,3)=grain(4,int(GBDat(i1,1)))
            segment(segments,4)=grain(2,int(GBDat(i1,2)))
            segment(segments,5)=grain(3,int(GBDat(i1,2)))
            segment(segments,6)=grain(4,int(GBDat(i1,2)))
            segment(segments,7)=ep1x
            segment(segments,8)=ep1y
            segment(segments,9)=ep2x
            segment(segments,10)=ep2y
            NodeList(i4,2)=0.0   !mark the TJ as used
            goto 910
           endif
 906       continue
          enddo  ! This ends the i4=1,TJNum loop
          ! Has the segment reached the max length so that we need to define another node?
          if (segs.eq.SegMax) then !write a segment
          ep1x=StartPointX ! begin by defining end point coordinates
          ep1y=StartPointY
          ep2x=LastPointX
          ep2y=LastPointY
          !this is where we check the distance from line
          do i5=1,segs
           !test to make sure gb nodes are not more than 2 times XSTEP away
           px=c_seg(i5,1)
           py=c_seg(i5,2)
           dist=abs((ep2x-ep1x)*(ep1y-py)-(ep1x-px)*(ep2y-ep1y))/(((ep2x-ep1x)**2)+((ep2y-ep1y)**2))**0.5
           if (dist.gt.DistTol*XSTEP) then
            !write(6,*)'distance too large',dist
            !write(6,*)px,py
            !If dist too large, make two segs from ep1x,ep1y to px,py and from px,py to ep2x,ep2y
            NewNodes=NewNodes+1
            NewNode(NewNodes,1)=px
            NewNode(NewNodes,2)=py
            segments=segments+1
            segment(segments,1)=grain(2,int(GBDat(i1,1)))
            segment(segments,2)=grain(3,int(GBDat(i1,1)))
            segment(segments,3)=grain(4,int(GBDat(i1,1)))
            segment(segments,4)=grain(2,int(GBDat(i1,2)))
            segment(segments,5)=grain(3,int(GBDat(i1,2)))
            segment(segments,6)=grain(4,int(GBDat(i1,2)))
            segment(segments,7)=ep1x
            segment(segments,8)=ep1y
            segment(segments,9)=px
            segment(segments,10)=py
            segments=segments+1
            segment(segments,1)=grain(2,int(GBDat(i1,1)))
            segment(segments,2)=grain(3,int(GBDat(i1,1)))
            segment(segments,3)=grain(4,int(GBDat(i1,1)))
            segment(segments,4)=grain(2,int(GBDat(i1,2)))
            segment(segments,5)=grain(3,int(GBDat(i1,2)))
            segment(segments,6)=grain(4,int(GBDat(i1,2)))
            segment(segments,7)=px
            segment(segments,8)=py
            segment(segments,9)=ep2x
            segment(segments,10)=ep2y
            segs=0
            goto 901
           endif
          enddo  ! ends the i5=1,SegMax loop
          !end distance check, save the segment
          vec(1)=abs(ep1x-ep2x)
          vec(2)=abs(ep1y-ep2y)
          vec(3)=0.0
          call VLen(vec,length) ! get length of vector
          if (length.eq.0.0) then !if this happens, something is wrong, skip it.
           search_fail=search_fail+1
           goto 1000
          endif
          segments=segments+1
          segment(segments,1)=grain(2,int(GBDat(i1,1)))
          segment(segments,2)=grain(3,int(GBDat(i1,1)))
          segment(segments,3)=grain(4,int(GBDat(i1,1)))
          segment(segments,4)=grain(2,int(GBDat(i1,2)))
          segment(segments,5)=grain(3,int(GBDat(i1,2)))
          segment(segments,6)=grain(4,int(GBDat(i1,2)))
          segment(segments,7)=ep1x
          segment(segments,8)=ep1y
          segment(segments,9)=ep2x
          segment(segments,10)=ep2y
          segs=0
          goto 901
         else
          goto 902
         endif
 910     continue
        enddo  ! ends the i2=1,TJNum loop
 1000   continue
       enddo  ! ends the i1=1,boundaries loop

       !This next section considers only cases where a grain pair has three nodes.
       !This seems to happen on square grids.  Find the two closest nodes - make
       !a segment, then connect the second to the most distant.
    
       do i1=1,boundaries ! loop over all identified boundaries
        trace=0.0
        dis=0.0
        c_seg=0.0
        search=0 ! keep track of attempts to find new boundaries
        if (int(GBDat(i1,4)).ne.3) then !ignore boundaries that do not have 3 TJ nodes
         goto 1050
        endif
        n1x=node(int(GBDat(i1,4+1)),1)
        n1y=node(int(GBDat(i1,4+1)),2)
        n2x=node(int(GBDat(i1,4+2)),1)
        n2y=node(int(GBDat(i1,4+2)),2)
        n3x=node(int(GBDat(i1,4+3)),1)
        n3y=node(int(GBDat(i1,4+3)),2)
        vec(1)=abs(n1x-n2x)
        vec(2)=abs(n1y-n2y)
        vec(3)=0
        call VLen(vec,d12) ! get length of vector
        vec(1)=abs(n1x-n3x)
        vec(2)=abs(n1y-n3y)
        vec(3)=0
        call VLen(vec,d13) ! get length of vector
        vec(1)=abs(n2x-n3x)
        vec(2)=abs(n2y-n3y)
        vec(3)=0
        call VLen(vec,d23) ! get length of vector
        !save the distances and the node IDs
        dis(1,1)=d12
        dis(1,2)=GBDat(i1,4+1)
        dis(1,3)=GBDat(i1,4+2)
        dis(2,1)=d13
        dis(2,2)=GBDat(i1,4+1)
        dis(2,3)=GBDat(i1,4+3)
        dis(3,1)=d23
        dis(3,2)=GBDat(i1,4+2)
        dis(3,3)=GBDat(i1,4+3)
        !next, we have to order these distances
        if (dis(2,1).gt.dis(3,1)) then !order components so 3>2>1
         tem=dis(2,1)
         dis(2,1)=dis(3,1)
         dis(3,1)=tem
         tem=dis(2,2)
         dis(2,2)=dis(3,2)
         dis(3,2)=tem
         tem=dis(2,3)
         dis(2,3)=dis(3,3)
         dis(3,3)=tem
        endif
        if (dis(1,1).gt.dis(2,1)) then
         tem=dis(1,1)
         dis(1,1)=dis(2,1)
         dis(2,1)=tem
         tem=dis(1,2)
         dis(1,2)=dis(2,2)
         dis(2,2)=tem
         tem=dis(1,3)
         dis(1,3)=dis(2,3)
         dis(2,3)=tem
        endif
        if (dis(2,1).gt.dis(3,1)) then
         tem=dis(2,1)
         dis(2,1)=dis(3,1)
         dis(3,1)=tem
         tem=dis(2,2)
         dis(2,2)=dis(3,2)
         dis(3,2)=tem
         tem=dis(2,3)
         dis(2,3)=dis(3,3)
         dis(3,3)=tem
        endif
        !It is important to start the search at the node furthest from the most distant node
        if (dis(1,3).eq.dis(3,2).OR.dis(1,3).eq.dis(3,3)) then
         tem=dis(1,2)
         dis(1,2)=dis(1,3)
         dis(1,3)=tem
        endif
        !now, dis(1,*) is the smallest and dis(2,*) is the next
        TraceLength=int(GBDat(i1,3))
        do i2=1,TraceLength  ! load the segment ids to the vector trace
         trace(i2,1)=GBDat(i1,14+i2)
         trace(i2,2)=1.0
        enddo ! this ends the i2=1,TraceLength loop
        StartPointX=node(int(dis(1,2)),1) !the start coordinate is the coordinate of the first node
        StartPointY=node(int(dis(1,2)),2)
        LastPointX=StartPointX
        LastPointY=StartPointY
 1001   continue !after writing a line of SegMax, redefine the start point
        StartPointX=LastPointX
        StartPointY=LastPointY
        search=search+1
        if (search.gt.199) then
         search_fail=search_fail+1
         goto 1050
        endif
        segs=0
 1002   continue
        dist=100
        do i3=1,TraceLength ! we find the closest GB node to this tridnode
         if (trace(i3,2).eq.0.0) then  !if 0, this GB node has already been used
          goto 1005  !goto next gb node
         endif
         !create a vector from the start coordinate to the gb node of interest
         vec(1)=LastPointX-gb(int(trace(i3,1)),1)
         vec(2)=LastPointY-gb(int(trace(i3,1)),2)
         vec(3)=0
         call VLen(vec,length) ! get length of vector
         if (length.lt.dist) then ! if the vector is the shortest one, remember it.
          dist=length
          ind=i3
         endif
 1005    continue
        enddo  !ends the i3=1,TraceLength loop
        segs=segs+1
        !Now we know the gb node closest to the start/last point
        trace(ind,2)=0.0 ! mark as used
        LastPointX=gb(int(trace(ind,1)),1)
        LastPointY=gb(int(trace(ind,1)),2)
        c_seg(segs,1)=LastPointX
        c_seg(segs,2)=LastPointY
        !Now we check if last point is near the end point
        !Find distance from last point to this node
        vec(1)=LastPointX-node(int(dis(1,3)),1)
        vec(2)=LastPointY-node(int(dis(1,3)),2)
        vec(3)=0.0
        call VLen(vec,length) ! get length of vector
        if (length.gt.2.0*XSTEP) then
         goto 1006 !look for the next gb node
        else
         !Draw a line from the start point to this node
         ep1x=StartPointX
         ep1y=StartPointY
         ep2x=node(int(dis(1,3)),1)
         ep2y=node(int(dis(1,3)),2)
         !Check the distance to boundary and break segent if it is too large
         do i5=1,segs
          px=c_seg(i5,1)
          py=c_seg(i5,2)
          dist=abs((ep2x-ep1x)*(ep1y-py)-(ep1x-px)*(ep2y-ep1y))/(((ep2x-ep1x)**2)+((ep2y-ep1y)**2))**0.5
          if (dist.gt.DistTol*XSTEP) then
           !write(6,*)'distance too large',dist
           !write(6,*)px,py
           !If dist too large, make two segs from ep1x,ep1y to px,py and from px,py to ep2x,ep2y
           NewNodes=NewNodes+1
           NewNode(NewNodes,1)=px
           NewNode(NewNodes,2)=py
           segments=segments+1
           segment(segments,1)=grain(2,int(GBDat(i1,1)))
           segment(segments,2)=grain(3,int(GBDat(i1,1)))
           segment(segments,3)=grain(4,int(GBDat(i1,1)))
           segment(segments,4)=grain(2,int(GBDat(i1,2)))
           segment(segments,5)=grain(3,int(GBDat(i1,2)))
           segment(segments,6)=grain(4,int(GBDat(i1,2)))
           segment(segments,7)=ep1x
           segment(segments,8)=ep1y
           segment(segments,9)=px
           segment(segments,10)=py
           segments=segments+1
           segment(segments,1)=grain(2,int(GBDat(i1,1)))
           segment(segments,2)=grain(3,int(GBDat(i1,1)))
           segment(segments,3)=grain(4,int(GBDat(i1,1)))
           segment(segments,4)=grain(2,int(GBDat(i1,2)))
           segment(segments,5)=grain(3,int(GBDat(i1,2)))
           segment(segments,6)=grain(4,int(GBDat(i1,2)))
           segment(segments,7)=px
           segment(segments,8)=py
           segment(segments,9)=ep2x
           segment(segments,10)=ep2y
           segs=0
           !NodeList(i4,2)=0.0   !mark the TJ as used
           goto 1010
          endif
         enddo  ! ends the i5=1,SegMax loop
         vec(1)=abs(ep1x-ep2x)
         vec(2)=abs(ep1y-ep2y)
         vec(3)=0.0
         call VLen(vec,length) ! get length of vector
         if (length.eq.0.0) then !if this happens, something is wrong, skip it.
          search_fail=search_fail+1
          goto 1050
         endif
         segments=segments+1
         segment(segments,1)=grain(2,int(GBDat(i1,1)))
         segment(segments,2)=grain(3,int(GBDat(i1,1)))
         segment(segments,3)=grain(4,int(GBDat(i1,1)))
         segment(segments,4)=grain(2,int(GBDat(i1,2)))
         segment(segments,5)=grain(3,int(GBDat(i1,2)))
         segment(segments,6)=grain(4,int(GBDat(i1,2)))
         segment(segments,7)=ep1x
         segment(segments,8)=ep1y
         segment(segments,9)=ep2x
         segment(segments,10)=ep2y
         goto 1010
        endif
 1006   continue

        ! Has the segment reached the max length so that we need to define another node?
        if (segs.eq.SegMax) then !write a segment
        ep1x=StartPointX ! begin by defining end point coordinates
        ep1y=StartPointY
        ep2x=LastPointX
        ep2y=LastPointY
        !this is where we should check the distance from line
        do i5=1,segs
         px=c_seg(i5,1)
         py=c_seg(i5,2)
         dist=abs((ep2x-ep1x)*(ep1y-py)-(ep1x-px)*(ep2y-ep1y))/(((ep2x-ep1x)**2)+((ep2y-ep1y)**2))**0.5
         if (dist.gt.DistTol*XSTEP) then
          !write(6,*)'distance too large',dist
          !write(6,*)px,py
          !If dist too large, make two segs from ep1x,ep1y to px,py and from px,py to ep2x,ep2y
          NewNodes=NewNodes+1
          NewNode(NewNodes,1)=px
          NewNode(NewNodes,2)=py
          segments=segments+1
          segment(segments,1)=grain(2,int(GBDat(i1,1)))
          segment(segments,2)=grain(3,int(GBDat(i1,1)))
          segment(segments,3)=grain(4,int(GBDat(i1,1)))
          segment(segments,4)=grain(2,int(GBDat(i1,2)))
          segment(segments,5)=grain(3,int(GBDat(i1,2)))
          segment(segments,6)=grain(4,int(GBDat(i1,2)))
          segment(segments,7)=ep1x
          segment(segments,8)=ep1y
          segment(segments,9)=px
          segment(segments,10)=py
          segments=segments+1
          segment(segments,1)=grain(2,int(GBDat(i1,1)))
          segment(segments,2)=grain(3,int(GBDat(i1,1)))
          segment(segments,3)=grain(4,int(GBDat(i1,1)))
          segment(segments,4)=grain(2,int(GBDat(i1,2)))
          segment(segments,5)=grain(3,int(GBDat(i1,2)))
          segment(segments,6)=grain(4,int(GBDat(i1,2)))
          segment(segments,7)=px
          segment(segments,8)=py
          segment(segments,9)=ep2x
          segment(segments,10)=ep2y
          segs=0
          !NodeList(i4,2)=0.0   !mark the TJ as used
          goto 1001
         endif
        enddo  ! ends the i5=1,SegMax loop
        vec(1)=abs(ep1x-ep2x)
        vec(2)=abs(ep1y-ep2y)
        vec(3)=0.0
        call VLen(vec,length) ! get length of vector
        if (length.eq.0.0) then !if this happens, something is wrong, skip it.
         search_fail=search_fail+1
         goto 1050
        endif
        segments=segments+1
        segment(segments,1)=grain(2,int(GBDat(i1,1)))
        segment(segments,2)=grain(3,int(GBDat(i1,1)))
        segment(segments,3)=grain(4,int(GBDat(i1,1)))
        segment(segments,4)=grain(2,int(GBDat(i1,2)))
        segment(segments,5)=grain(3,int(GBDat(i1,2)))
        segment(segments,6)=grain(4,int(GBDat(i1,2)))
        segment(segments,7)=ep1x
        segment(segments,8)=ep1y
        segment(segments,9)=ep2x
        segment(segments,10)=ep2y
        segs=0
        goto 1001
       else
        goto 1002
       endif
 1010  continue

        StartPointX=node(int(dis(2,2)),1) !the start coordinate is the coordinate of the first node
        StartPointY=node(int(dis(2,2)),2)
        LastPointX=StartPointX
        LastPointY=StartPointY
 1011   continue !after writing a line of SegMax, redefine the start point
        StartPointX=LastPointX
        StartPointY=LastPointY
        search=search+1
        if (search.gt.199) then
         search_fail=search_fail+1
         goto 1050
        endif
        segs=0
 1012   continue
        dist=100
        do i3=1,TraceLength ! we find the closest GB node to this tridnode
         if (trace(i3,2).eq.0.0) then  !if 0, this GB node has already been used
          goto 1015  !goto next gb node
         endif
         !create a vector from the start coordinate to the gb node of interest
         vec(1)=LastPointX-gb(int(trace(i3,1)),1)
         vec(2)=LastPointY-gb(int(trace(i3,1)),2)
         vec(3)=0
         call VLen(vec,length) ! get length of vector
         if (length.lt.dist) then ! if the vector is the shortest one, remember it.
          dist=length
          ind=i3
         endif
 1015    continue
        enddo  !ends the i3=1,TraceLength loop
        segs=segs+1
        !Now we know the gb node closest to the start/last point
        trace(ind,2)=0.0 ! mark as used
        LastPointX=gb(int(trace(ind,1)),1)
        LastPointY=gb(int(trace(ind,1)),2)
        c_seg(segs,1)=LastPointX
        c_seg(segs,2)=LastPointY
        !Now we check if last point is near the end point
        !Find distance from last point to this node
        vec(1)=LastPointX-node(int(dis(2,3)),1)
        vec(2)=LastPointY-node(int(dis(2,3)),2)
        vec(3)=0.0
        call VLen(vec,length) ! get length of vector
        if (length.gt.2.0*XSTEP) then
         goto 1016 !look for the next gb node
        else
         !Draw a line from the start point to this node
         ep1x=StartPointX
         ep1y=StartPointY
         ep2x=node(int(dis(2,3)),1)
         ep2y=node(int(dis(2,3)),2)
         !Check the distance to boundary and break segent if it is too large
         do i5=1,segs
          px=c_seg(i5,1)
          py=c_seg(i5,2)
          dist=abs((ep2x-ep1x)*(ep1y-py)-(ep1x-px)*(ep2y-ep1y))/(((ep2x-ep1x)**2)+((ep2y-ep1y)**2))**0.5
          if (dist.gt.DistTol*XSTEP) then
           !write(6,*)'distance too large',dist
           !write(6,*)px,py
           !If dist too large, make two segs from ep1x,ep1y to px,py and from px,py to ep2x,ep2y
           NewNodes=NewNodes+1
           NewNode(NewNodes,1)=px
           NewNode(NewNodes,2)=py
           segments=segments+1
           segment(segments,1)=grain(2,int(GBDat(i1,1)))
           segment(segments,2)=grain(3,int(GBDat(i1,1)))
           segment(segments,3)=grain(4,int(GBDat(i1,1)))
           segment(segments,4)=grain(2,int(GBDat(i1,2)))
           segment(segments,5)=grain(3,int(GBDat(i1,2)))
           segment(segments,6)=grain(4,int(GBDat(i1,2)))
           segment(segments,7)=ep1x
           segment(segments,8)=ep1y
           segment(segments,9)=px
           segment(segments,10)=py
           segments=segments+1
           segment(segments,1)=grain(2,int(GBDat(i1,1)))
           segment(segments,2)=grain(3,int(GBDat(i1,1)))
           segment(segments,3)=grain(4,int(GBDat(i1,1)))
           segment(segments,4)=grain(2,int(GBDat(i1,2)))
           segment(segments,5)=grain(3,int(GBDat(i1,2)))
           segment(segments,6)=grain(4,int(GBDat(i1,2)))
           segment(segments,7)=px
           segment(segments,8)=py
           segment(segments,9)=ep2x
           segment(segments,10)=ep2y
           segs=0
           !NodeList(i4,2)=0.0   !mark the TJ as used
           goto 1020
          endif
         enddo  ! ends the i5=1,SegMax loop
         vec(1)=abs(ep1x-ep2x)
         vec(2)=abs(ep1y-ep2y)
         vec(3)=0.0
         call VLen(vec,length) ! get length of vector
         if (length.eq.0.0) then !if this happens, something is wrong, skip it.
          search_fail=search_fail+1
          goto 1050
         endif
         segments=segments+1
         segment(segments,1)=grain(2,int(GBDat(i1,1)))
         segment(segments,2)=grain(3,int(GBDat(i1,1)))
         segment(segments,3)=grain(4,int(GBDat(i1,1)))
         segment(segments,4)=grain(2,int(GBDat(i1,2)))
         segment(segments,5)=grain(3,int(GBDat(i1,2)))
         segment(segments,6)=grain(4,int(GBDat(i1,2)))
         segment(segments,7)=ep1x
         segment(segments,8)=ep1y
         segment(segments,9)=ep2x
         segment(segments,10)=ep2y
         goto 1020
        endif
 1016   continue

        ! Has the segment reached the max length so that we need to define another node?
        if (segs.eq.SegMax) then !write a segment
        ep1x=StartPointX ! begin by defining end point coordinates
        ep1y=StartPointY
        ep2x=LastPointX
        ep2y=LastPointY
        !this is where we should check the distance from line
        do i5=1,segs
         px=c_seg(i5,1)
         py=c_seg(i5,2)
         dist=abs((ep2x-ep1x)*(ep1y-py)-(ep1x-px)*(ep2y-ep1y))/(((ep2x-ep1x)**2)+((ep2y-ep1y)**2))**0.5
         if (dist.gt.DistTol*XSTEP) then
          !write(6,*)'distance too large',dist
          !write(6,*)px,py
          !If dist too large, make two segs from ep1x,ep1y to px,py and from px,py to ep2x,ep2y
          NewNodes=NewNodes+1
          NewNode(NewNodes,1)=px
          NewNode(NewNodes,2)=py
          segments=segments+1
          segment(segments,1)=grain(2,int(GBDat(i1,1)))
          segment(segments,2)=grain(3,int(GBDat(i1,1)))
          segment(segments,3)=grain(4,int(GBDat(i1,1)))
          segment(segments,4)=grain(2,int(GBDat(i1,2)))
          segment(segments,5)=grain(3,int(GBDat(i1,2)))
          segment(segments,6)=grain(4,int(GBDat(i1,2)))
          segment(segments,7)=ep1x
          segment(segments,8)=ep1y
          segment(segments,9)=px
          segment(segments,10)=py
          segments=segments+1
          segment(segments,1)=grain(2,int(GBDat(i1,1)))
          segment(segments,2)=grain(3,int(GBDat(i1,1)))
          segment(segments,3)=grain(4,int(GBDat(i1,1)))
          segment(segments,4)=grain(2,int(GBDat(i1,2)))
          segment(segments,5)=grain(3,int(GBDat(i1,2)))
          segment(segments,6)=grain(4,int(GBDat(i1,2)))
          segment(segments,7)=px
          segment(segments,8)=py
          segment(segments,9)=ep2x
          segment(segments,10)=ep2y
          segs=0
          !NodeList(i4,2)=0.0   !mark the TJ as used
          goto 1011
         endif
        enddo  ! ends the i5=1,SegMax loop
        !end distance check, save the segment
        vec(1)=abs(ep1x-ep2x)
        vec(2)=abs(ep1y-ep2y)
        vec(3)=0.0
        call VLen(vec,length) ! get length of vector
        if (length.eq.0.0) then !if this happens, something is wrong, skip it.
         search_fail=search_fail+1
         goto 1050
        endif
        segments=segments+1
        segment(segments,1)=grain(2,int(GBDat(i1,1)))
        segment(segments,2)=grain(3,int(GBDat(i1,1)))
        segment(segments,3)=grain(4,int(GBDat(i1,1)))
        segment(segments,4)=grain(2,int(GBDat(i1,2)))
        segment(segments,5)=grain(3,int(GBDat(i1,2)))
        segment(segments,6)=grain(4,int(GBDat(i1,2)))
        segment(segments,7)=ep1x
        segment(segments,8)=ep1y
        segment(segments,9)=ep2x
        segment(segments,10)=ep2y
        segs=0
        goto 1011
       else
        goto 1012
       endif

 1020  continue

 1050  continue

      enddo  ! ends the i1=1,boundaries loop

       !The next case considers only island grains - those grains with zero multinodes

       do i1=1,boundaries ! loop over all identified boundaries
        if (int(GBDat(i1,4)).ne.0) then !ignore boundaries that do not have 0 TJ nodes
         goto 1200
        endif
        !Here, we start at the first grain boundary node and progress until
        !we get back to it.
        TraceLength=int(GBDat(i1,3))  !Number of gb nodes
        do i2=1,TraceLength  ! load the segment ids to the vector trace
         trace(i2,1)=GBDat(i1,14+i2)
         trace(i2,2)=1.0
        enddo ! this ends the i2=1,TraceLength loop
        FromOrigin=0
        OriginPointX=gb(int(trace(1,1)),1) !start at the the first gb node in the list
        OriginPointY=gb(int(trace(1,1)),2)
        trace(1,2)=0.0 ! mark as used
        LastPointX=OriginPointX
        LastPointY=OriginPointY
 1110   continue !after writing a line of SegMax, redefine the start point
        StartPointX=LastPointX
        StartPointY=LastPointY
        segs=0
 1120   continue
        dist=100
        do i3=2,TraceLength ! we find the closest GB node to the origin
         if (trace(i3,2).eq.0.0) then  !if 0, this GB node has already been used
          goto 1130  !goto next gb node
         endif
         !create a vector from the start coordinate to the gb node of interest
         vec(1)=LastPointX-gb(int(trace(i3,1)),1)
         vec(2)=LastPointY-gb(int(trace(i3,1)),2)
         vec(3)=0
         call VLen(vec,length) ! get length of vector
         if (length.lt.dist) then ! if the vector is the shortest one, remember it.
          dist=length
          ind=i3
         endif
 1130     continue
         enddo  !ends the i3=2,TraceLength loop
         segs=segs+1
         FromOrigin=FromOrigin+1
         !Now we know the gb node closest to the start/last point
         trace(ind,2)=0.0 ! mark as used
         LastPointX=gb(int(trace(ind,1)),1)
         LastPointY=gb(int(trace(ind,1)),2)
         c_seg(segs,1)=LastPointX
         c_seg(segs,2)=LastPointY
         !Now we check if last point is near the origin point
         if (FromOrigin.gt.(TraceLength-4)) then ! close the loop
          ep1x=StartPointX
          ep1y=StartPointY
          ep2x=OriginPointX
          ep2y=OriginPointY
          !Check the distance to boundary and break segent if it is too large
          do i5=1,segs
           px=c_seg(i5,1)
           py=c_seg(i5,2)
           dist=abs((ep2x-ep1x)*(ep1y-py)-(ep1x-px)*(ep2y-ep1y))/(((ep2x-ep1x)**2)+((ep2y-ep1y)**2))**0.5
           if (dist.gt.DistTol*XSTEP) then
            !write(6,*)'distance too large',dist
            !write(6,*)px,py
            !If dist too large, make two segs from ep1x,ep1y to px,py and from px,py to ep2x,ep2y
            NewNodes=NewNodes+1
            NewNode(NewNodes,1)=px
            NewNode(NewNodes,2)=py
            segments=segments+1
            segment(segments,1)=grain(2,int(GBDat(i1,1)))
            segment(segments,2)=grain(3,int(GBDat(i1,1)))
            segment(segments,3)=grain(4,int(GBDat(i1,1)))
            segment(segments,4)=grain(2,int(GBDat(i1,2)))
            segment(segments,5)=grain(3,int(GBDat(i1,2)))
            segment(segments,6)=grain(4,int(GBDat(i1,2)))
            segment(segments,7)=ep1x
            segment(segments,8)=ep1y
            segment(segments,9)=px
            segment(segments,10)=py
            segments=segments+1
            segment(segments,1)=grain(2,int(GBDat(i1,1)))
            segment(segments,2)=grain(3,int(GBDat(i1,1)))
            segment(segments,3)=grain(4,int(GBDat(i1,1)))
            segment(segments,4)=grain(2,int(GBDat(i1,2)))
            segment(segments,5)=grain(3,int(GBDat(i1,2)))
            segment(segments,6)=grain(4,int(GBDat(i1,2)))
            segment(segments,7)=px
            segment(segments,8)=py
            segment(segments,9)=ep2x
            segment(segments,10)=ep2y
            segs=0
            !NodeList(i4,2)=0.0   !mark the TJ as used
            goto 1200
           endif
          enddo  ! ends the i5=1,SegMax loop
          vec(1)=abs(ep1x-ep2x)
          vec(2)=abs(ep1y-ep2y)
          vec(3)=0.0
          call VLen(vec,length) ! get length of vector
          if (length.eq.0.0) then !if this happens, something is wrong, skip it.
           search_fail=search_fail+1
           goto 1200
          endif
          segments=segments+1
          segment(segments,1)=grain(2,int(GBDat(i1,1)))
          segment(segments,2)=grain(3,int(GBDat(i1,1)))
          segment(segments,3)=grain(4,int(GBDat(i1,1)))
          segment(segments,4)=grain(2,int(GBDat(i1,2)))
          segment(segments,5)=grain(3,int(GBDat(i1,2)))
          segment(segments,6)=grain(4,int(GBDat(i1,2)))
          segment(segments,7)=ep1x
          segment(segments,8)=ep1y
          segment(segments,9)=ep2x
          segment(segments,10)=ep2y
          goto 1200 !end path, check next boundary
         endif
         ! Has the segment reached the max length so that we need to define another node?
         if (segs.eq.SegMax) then !write a segment
          !define the endpoint coordinates
          ep1x=StartPointX
          ep1y=StartPointY
          ep2x=LastPointX
          ep2y=LastPointY
          !Check the distance to boundary and break segent if it is too large
          do i5=1,segs
           px=c_seg(i5,1)
           py=c_seg(i5,2)
           dist=abs((ep2x-ep1x)*(ep1y-py)-(ep1x-px)*(ep2y-ep1y))/(((ep2x-ep1x)**2)+((ep2y-ep1y)**2))**0.5
           if (dist.gt.DistTol*XSTEP) then
            !write(6,*)'distance too large',dist
            !write(6,*)px,py
            !If dist too large, make two segs from ep1x,ep1y to px,py and from px,py to ep2x,ep2y
            NewNodes=NewNodes+1
            NewNode(NewNodes,1)=px
            NewNode(NewNodes,2)=py
            segments=segments+1
            segment(segments,1)=grain(2,int(GBDat(i1,1)))
            segment(segments,2)=grain(3,int(GBDat(i1,1)))
            segment(segments,3)=grain(4,int(GBDat(i1,1)))
            segment(segments,4)=grain(2,int(GBDat(i1,2)))
            segment(segments,5)=grain(3,int(GBDat(i1,2)))
            segment(segments,6)=grain(4,int(GBDat(i1,2)))
            segment(segments,7)=ep1x
            segment(segments,8)=ep1y
            segment(segments,9)=px
            segment(segments,10)=py
            segments=segments+1
            segment(segments,1)=grain(2,int(GBDat(i1,1)))
            segment(segments,2)=grain(3,int(GBDat(i1,1)))
            segment(segments,3)=grain(4,int(GBDat(i1,1)))
            segment(segments,4)=grain(2,int(GBDat(i1,2)))
            segment(segments,5)=grain(3,int(GBDat(i1,2)))
            segment(segments,6)=grain(4,int(GBDat(i1,2)))
            segment(segments,7)=px
            segment(segments,8)=py
            segment(segments,9)=ep2x
            segment(segments,10)=ep2y
            segs=0
            !NodeList(i4,2)=0.0   !mark the TJ as used
            goto 1110
           endif
          enddo  ! ends the i5=1,SegMax loop
          ! save this line segment
          vec(1)=abs(ep1x-ep2x)
          vec(2)=abs(ep1y-ep2y)
          vec(3)=0.0
          call VLen(vec,length) ! get length of vector
          if (length.eq.0.0) then !if this happens, something is wrong, skip it.
           search_fail=search_fail+1
           goto 1200
          endif
          segments=segments+1
          segment(segments,1)=grain(2,int(GBDat(i1,1)))
          segment(segments,2)=grain(3,int(GBDat(i1,1)))
          segment(segments,3)=grain(4,int(GBDat(i1,1)))
          segment(segments,4)=grain(2,int(GBDat(i1,2)))
          segment(segments,5)=grain(3,int(GBDat(i1,2)))
          segment(segments,6)=grain(4,int(GBDat(i1,2)))
          segment(segments,7)=ep1x
          segment(segments,8)=ep1y
          segment(segments,9)=ep2x
          segment(segments,10)=ep2y
          segs=0
          goto 1110
         else
          goto 1120
         endif

 1200   continue
       enddo  ! ends the i1=1,boundaries loop

       write(6,3)'There are :',segments,' grain boundary trace segments'
       write(6,*)search_fail,' failed searches'
       write(6,21)'The number of traces subdivided = ',NewNodes

        !------------------------------------------------------------------
        ! The rest of this program writes out the data that was requested
        ! by the input file.
        !------------------------------------------------------------------

        !Write out the segment data
        !open(43, file=trim(out_fname)//'_'//fname(1:point)//'.txt',status='unknown')
        open(43, file=trim(out_fname)//'_'//one//two//three//'.txt',status='unknown')
        write(43,2)'# File written by extract_gb_traces: ',version
        write(43,2)'# Using data from file: ',fname
        if (source.eq.0) then
         write(43,1)'# data from TSL'
        endif
        if (grid.eq.0) then
         write(43,1)'# data is on a HEX grid'
        else
         write(43,1)'# data is on a square grid'
        endif
        write(43,22)'# length of step in x direction: ',XSTEP
        write(43,22)'# length of step in y direction: ',YSTEP
        write(43,21)'# number of rows: ',NROWS
        write(43,21)'# number of columns in odd numbered rows: ',NCOLS_ODD
        write(43,21)'# number of columns in even numbered rows: ',NCOLS_EVEN
        write(43,21)'# maximal number of pixels in each line segments: ',SegMax
        write(43,2)'# ',inline2
        write(43,3)'# There are :',grains,' grains'
        write(43,3)'# There are :',nodes,' trinodes'
        write(43,3)'# There are :',gbs,' grain boundary nodes'
        write(43,3)'# There are :',boundaries,' grain pairs with boundaries'
        write(43,21)'# The number of traces subdivided: ',NewNodes
        !Write some information about the types of nodes.  This is mainly diagnostic
        ct=0
        do i1=1,boundaries
         if (int(GBDat(i1,4)).ne.2) then
          ct=ct+1
         endif
        enddo
        write(43,20)'# ',ct,' of ',boundaries,' boundaries do not have 2 trinodes'
        ct=0
        do i1=1,boundaries
         if (int(GBDat(i1,4)).eq.0) then
          ct=ct+1
         endif
        enddo
        write(43,20)'# ',ct,' of ',boundaries,' boundaries have 0 trinode'
        ct=0
        do i1=1,boundaries
         if (int(GBDat(i1,4)).eq.1) then
          ct=ct+1
         endif
        enddo
        write(43,20)'# ',ct,' of ',boundaries,' boundaries have 1 trinode'
        ct=0
        do i1=1,boundaries
         if (int(GBDat(i1,4)).eq.2) then
          ct=ct+1
         endif
        enddo
        write(43,20)'# ',ct,' of ',boundaries,' boundaries have 2 trinodes'
        ct=0
        do i1=1,boundaries
         if (int(GBDat(i1,4)).eq.3) then
          ct=ct+1
         endif
        enddo
        write(43,20)'# ',ct,' of ',boundaries,' boundaries have 3 trinodes'
        ct=0
        do i1=1,boundaries
          if (int(GBDat(i1,4)).eq.4) then
           ct=ct+1
          endif
         enddo
        write(43,20)'# ',ct,' of ',boundaries,' boundaries have 4 trinodes'
        ct=0
        do i1=1,boundaries
         if (int(GBDat(i1,4)).eq.5) then
          ct=ct+1
         endif
        enddo
        write(43,20)'# ',ct,' of ',boundaries,' boundaries have 5 trinodes'
        ct=0
        do i1=1,boundaries
         if (int(GBDat(i1,4)).eq.6) then
          ct=ct+1
         endif
        enddo
        write(43,20)'# ',ct,' of ',boundaries,' boundaries have 6 trinodes'
        ct=0
        do i1=1,boundaries
         if (int(GBDat(i1,4)).eq.7) then
          ct=ct+1
         endif
        enddo
        write(43,20)'# ',ct,' of ',boundaries,' boundaries have 7 trinodes'
        ct=0
        do i1=1,boundaries
         if (int(GBDat(i1,4)).eq.8) then
          ct=ct+1
         endif
        enddo
        write(43,20)'# ',ct,' of ',boundaries,' boundaries have 8 trinodes'
        !should write all the parameters to the header here
        write (43,1)'#  phi1    PHI   phi2   phi1    PHI   phi2       x1       y1       x2       y2'
        do i1=1,segments
         write(43,14)segment(i1,1),segment(i1,2),segment(i1,3),segment(i1,4),segment(i1,5),segment(i1,6),segment(i1,7),segment(i1,8),segment(i1,9),segment(i1,10)
        enddo
        close(43)
  
        !write IPF map with nodes and lines to a postscript file
        if (out(1).eq.1) then
        open (59, file='map_'//fname(1:point)//'.ps',status='unknown') ! open the postscript file
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
        write (59,1)'1.0 setlinewidth '
        do i1=1,segments
         write(59,15)XOff+(XSc*segment(i1,7)),YOff+(YSc*segment(i1,8)),' newpath moveto'
         write(59,15)XOff+(XSc*segment(i1,9)),YOff+(YSc*segment(i1,10)),' lineto'
         write(59,1)'0 setgray'
         write(59,1)'stroke'
        enddo

        !draw all the gb nodes on the map
        write (59,1)'0.1 setlinewidth '
        do i1=1,gbs
         write (59,1)'0 setgray'
         write(59,15)XOff+(XSc*gb(i1,1)),YOff+(YSc*gb(i1,2)),' 0.6 0 360 arc closepath'
         write(59,1)'gsave'
         !write (59,1)'1 0 0 setrgbcolor fill'
         write (59,1)'grestore'
         write (59,1)'stroke'
        enddo

        !draw all the tri nodes on the map
        write (59,1)'0.1 setlinewidth '
        do i1=1,nodes
         write (59,1)'1 0 0 setrgbcolor'
         write(59,15)XOff+(XSc*node(i1,1)),YOff+(YSc*node(i1,2)),' 0.6 0 360 arc closepath'
         write(59,1)'gsave'
         write (59,1)'1 0 0 setrgbcolor fill'
         write (59,1)'grestore'
         write (59,1)'stroke'
        enddo

        !draw all the new nodes on the map
        write (59,1)'0.1 setlinewidth '
        do i1=1,NewNodes
         write (59,1)'0 0 1 setrgbcolor'
         write(59,15)XOff+(XSc*NewNode(i1,1)),YOff+(YSc*NewNode(i1,2)),' 1.2 0 360 arc closepath'
         write(59,1)'gsave'
         write (59,1)'0 0 1 setrgbcolor fill'
         write (59,1)'grestore'
         write (59,1)'stroke'
        enddo

        !write(59,1)'showpage'
        !close(59) ! close the IPF map
        !endif
     
        !Mark all grain IDs on the map
        if (out(7).eq.1) then
         write (59,1)'0.1 setlinewidth '
         write (59,1)'/Helvetica 8 selectfont'
         write (59,1)'0 setgray'
         do i1=1,grains
          !write (59,1)'0 1 0 setrgbcolor'
          write(59,15)XOff+(XSc*grain(5,i1)),YOff+(YSc*grain(6,i1)),' moveto'
          write(IDLabel,24)int(grain(1,i1))
          mystring=trim(adjustl(IDLabel))
          write(59,25)'(',mystring,') show'
         enddo
        endif

        !here we put some labels on the page

        write (59,1)'/Helvetica 12 selectfont'
        write (59,1)'0 setgray'
        write(59,15)XOff+(10.),YOff+(YSc*(NROWS*YSTEP))+62,' moveto'
        mystring2=trim(fname)
        write(59,25)'(',mystring2,') show'

        write (59,1)'1.0 setlinewidth '
        write(59,15)XOff+(10.),YOff+(YSc*(NROWS*YSTEP))+51,' newpath moveto'
        write(59,15)XOff+(25.),YOff+(YSc*(NROWS*YSTEP))+51,' lineto'
        write(59,1)'0 setgray'
        write(59,1)'stroke'
        write (59,1)'/Helvetica 12 selectfont'
        write (59,1)'0 setgray'
        write(59,15)XOff+(30.),YOff+(YSc*(NROWS*YSTEP))+47,' moveto'
        mystring2='GB line segment'
        mystring2=trim(mystring2)
        write(59,25)'(',mystring2,') show'

        write (59,1)'0.1 setlinewidth '
        write (59,1)'1 0 0 setrgbcolor'
        write(59,15)XOff+(18.),YOff+(YSc*(NROWS*YSTEP))+36,' moveto'
        write(59,15)XOff+(18.),YOff+(YSc*(NROWS*YSTEP))+36,' 0.6 0 360 arc closepath'
        write(59,1)'gsave'
        write (59,1)'1 0 0 setrgbcolor fill'
        write (59,1)'grestore'
        write (59,1)'stroke'
        write(59,15)XOff+(30.),YOff+(YSc*(NROWS*YSTEP))+32,' moveto'
        mystring2='GB trinode'
        mystring2=trim(mystring2)
        write(59,25)'(',mystring2,') show'

        write(59,15)XOff+(18.6),YOff+(YSc*(NROWS*YSTEP))+21,' newpath moveto'
        write (59,1)'0.1 setlinewidth '
        write (59,1)'0 setgray'
        write(59,15)XOff+(18.),YOff+(YSc*(NROWS*YSTEP))+21,' 0.6 0 360 arc closepath'
        write(59,1)'gsave'
        write (59,1)'grestore'
        write (59,1)'stroke'
        write(59,15)XOff+(30.),YOff+(YSc*(NROWS*YSTEP))+17,' moveto'
        mystring2='GB binode'
        mystring2=trim(mystring2)
        write(59,25)'(',mystring2,') show'

  
        write (59,1)'0.1 setlinewidth '
        write (59,1)'0 0 1 setrgbcolor'
        write(59,15)XOff+(18.),YOff+(YSc*(NROWS*YSTEP))+6,' moveto'
        write(59,15)XOff+(18.),YOff+(YSc*(NROWS*YSTEP))+6,' 1.2 0 360 arc closepath'
        write(59,1)'gsave'
        write (59,1)'0 0 1 setrgbcolor fill'
        write (59,1)'grestore'
        write (59,1)'stroke'
        write(59,15)XOff+(30.),YOff+(YSc*(NROWS*YSTEP))+2,' moveto'
        mystring2='seg break node'
        mystring2=trim(mystring2)
        write(59,25)'(',mystring2,') show'


        write(59,1)'showpage'
        close(59) ! close the IPF map

       endif



        !Write out the grain list, if requested
        if (out(2).eq.1) then
         open(43, file=trim(out_fname)//'_'//one//two//three//'.txt',status='unknown')
         open (40, file='grain_list_'//fname(1:point)//'.txt',status='unknown')
         do i1=1,25  !Read the header lines from the segment file and write to output
          read(43, "(a)") inline
          write(40, "(a)") inline
         enddo
         close(43)
         write(40,2)'# list of grains in: ',trim(fname)
         write(40,2)'# written by extract_GB_traces, ',version
         write(40,1)'# grain ID    phi1     PHI    phi2       X       Y'
         do i1=1,grains
          write(40,5)int(grain(1,i1)),grain(2,i1),grain(3,i1),grain(4,i1),grain(5,i1),grain(6,i1)
         enddo
         close(40)
        endif

        !Write out the grain ID map file, if requested
         if (out(3).eq.1) then
          open(43, file=trim(out_fname)//'_'//one//two//three//'.txt',status='unknown')
          open (40, file='grain_ID_map_'//fname(1:point)//'.txt',status='unknown')
          do i1=1,25  !Read the header lines from the segment file and write to output
           read(43, "(a)") inline
           write(40, "(a)") inline
          enddo
          close(43)
          write(40,2)'# Map of grain IDs in: ',trim(fname)
          write(40,2)'# written by extract_GB_traces, ',version
          write(40,1)'# grain ID   x cood    y coord'
          do i1=1,NROWS
           do i2=1,NCOLS_EVEN
            write(40,6)int(map(i1,i2,1)),map(i1,i2,2),map(i1,i2,3)
           enddo
          enddo
          close(40)
         endif

         !Write out the triple junction node list, if requested
         if (out(4).eq.1) then
          open(43, file=trim(out_fname)//'_'//one//two//three//'.txt',status='unknown')
          open (41, file='grain_tj_node_list_'//fname(1:point)//'.txt',status='unknown')
          do i1=1,25  !Read the header lines from the segment file and write to output
           read(43, "(a)") inline
           write(41, "(a)") inline
          enddo
          close(43)
          write(41,2)'# list of grains in: ',trim(fname)
          write(41,2)'# written by extract_GB_traces, ',version
          write(41,1)'#        X         Y    GID1    GID2    GID3'
          do i1=1,nodes
          w rite(41,12)node(i1,1),node(i1,2),int(node(i1,3)),int(node(i1,4)),int(node(i1,5))
          enddo
          close(41)
         endif

         !Write out the list of grain boundary nodes, if requested
         if (out(5).eq.1) then
          open(43, file=trim(out_fname)//'_'//one//two//three//'.txt',status='unknown')
          open (42, file='gb_node_list_'//fname(1:point)//'.txt',status='unknown')
          do i1=1,25  !Read the header lines from the segment file and write to output
           read(43, "(a)") inline
           write(42, "(a)") inline
          enddo
          close(43)
          write(42,2)'# list of grain boundary nodes in: ',trim(fname)
          write(42,2)'# written by extract_GB_traces, ',version
          write(42,1)'#        X         Y    GID1    GID2'
          do i1=1,gbs
           write(42,23)gb(i1,1),gb(i1,2),int(gb(i1,3)),int(gb(i1,4))
          enddo
          close(42)
         endif

         !Write out the list of the GB_dat file, if requested
         if (out(6).eq.1) then
          open(43, file=trim(out_fname)//'_'//one//two//three//'.txt',status='unknown')
          open (44, file='gb_data_'//fname(1:point)//'.txt',status='unknown')
          do i1=1,25  !Read the header lines from the segment file and write to output
           read(43, "(a)") inline
           write(44, "(a)") inline
          enddo
          close(43)
          write(44,2)'# list of grain boundary information: ',trim(fname)
          write(44,2)'# written by extract_GB_traces, ',version
          write(44,1)'# GID1  GID2  #gbs  #tjs   tj1   tj2   tj3   tj4   tj5   tj6   tj7   tj8   tj9  tj10   gb1   gb2   gb3   gb4   gb5   gb6   gb7   gb8   gb9  gb10  gb11'
          do i1=1,boundaries
           write(44,13)int(GBDat(i1,1)),int(GBDat(i1,2)),int(GBDat(i1,3)),int(GBDat(i1,4)),int(GBDat(i1,5)),int(GBDat(i1,6)),int(GBDat(i1,7)),int(GBDat(i1,8)),int(GBDat(i1,9)),int(GBDat(i1,10)),int(GBDat(i1,11)),int(GBDat(i1,12)),int(GBDat(i1,13)),int(GBDat(i1,14)),int(GBDat(i1,15)),int(GBDat(i1,16)),int(GBDat(i1,17)),int(GBDat(i1,18)),int(GBDat(i1,19)),int(GBDat(i1,20)),int(GBDat(i1,21)),int(GBDat(i1,22)),int(GBDat(i1,23)),int(GBDat(i1,24)),int(GBDat(i1,25))
          enddo
          close(44)
         endif


       !make a map of the nodes -
       !open (52, file='node_map.ps',status='unknown')
       !write (52,1)'%!PS'
       !do i1=1,nodes
        !write (52,1)'1 0 0 setrgbcolor'
        !write (52,1)'0.1 setlinewidth '
        !write(52,15)50+(6*node(i1,1)),50+(6*node(i1,2)),' 0.6 0 360 arc closepath'
        !write(52,1)'gsave'
        !write (52,1)'1 0 0 setrgbcolor fill'
        !write (52,1)'grestore'
        !write (52,1)'stroke'
       !enddo
       !write(52,1)'showpage'
       !close(52)






       enddo !closes the nfile=first,last loop

 9000  continue

       write(6,1) 'Program complete'

       end






	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
