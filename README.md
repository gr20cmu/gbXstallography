# gbXstallography
programs used to compute and plot distributions of grain boundary properties

1. condition_segments.  The program writes lists of grain boundary line segments in a format that is compatible with all subsequent programs
in the gbXstallography repository.  If the line segments were written by TSL, then they will have 12, 14, or 21 columns.  All programs in 
this repsository expect files with a 10 or 12 column format.  Each grain boundary line segment has 10 essential pieces of information:  
The three Euler angles for the crystal on one side of the boundary, the three Euler angles for the crystal on the other side of the boundary, 
the x and y coordinates of the start of the segment, and the x and y coordinates of the end of the segment.  
The program also eliminates lines if the Euler angles are exactly zero or if the line segment has zero length.
There are two main modes of operation:
- concatenate mode.  All line segments from multiple files are written to a single file.
- do not concatenate mode.  Each file is rewritten in the desired formate, but with the string '_cnd_' added to the name.


