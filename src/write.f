!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!* writetofile writes the positions to file every savenum              *
!* It rewrites the positions in the same file, deleting the            *
!* old position, thus not taking too much space.  The output is written*
!* to two alternating files so if the program crashes while writing    *
!* the data will not be lost                                           *
!* Modified 10/6 to write failures also                                *
!***********************************************************************

      subroutine writetofile(filei)
      include 'MD.com'
      integer filei
      character(LEN=20) filename
      character(LEN=15) traj_fmt
      character(LEN=15) traj_fmt_min
      parameter(traj_fmt="(3F8.3)")
      parameter(traj_fmt_min="(3F12.3)")
      real res
      dimension Res(AN)

      if(startt .eq. 1 .or. startt .eq. 2)then
         filename = initval
      elseif(startt .eq. 3)then
         filename = conf
      endif

      if(filei .eq. 0)then
      open(61, FILE= finalpx1, access = 'SEQUENTIAL', 
     Q   status = 'unknown')
      elseif(filei .eq. 1)then
      open(61, FILE= 'temp', access = 'SEQUENTIAL', 
     Q   status = 'unknown')
      endif


	call findtemp(KE)
	write(61,*) ' # There has been', writecount, 'steps'
	write(61,*) '# The current temp is', temp
        write(61,*) '# Number of atoms'
	write(61,*) AN
        write(61,*) '# X(i)', 'Y(i)', 'Z(i)'

        if(Trajectory .ne. 'NO')then
            if(minTrajOut .ne. 'YES')then
               write(6,*) writecount
            endif
	      endif

      do  i=1, AN
        write(61,FMTX) BeadIndex(i), GroupIndex(i),AtType(i),
     Q                 ResID(i), X(i), Y(i), Z(i), ms(i)
      end do


      if (ExcludeStatsFromTraj .eq. 'YES') then
      ! this is intend to avoid printing static atoms coordinates
      ! and to be more economic with the size of the files
      ! AMIR , 28-3-2010
         do  i=1, DynamicAtomRange(DynLength)
            if(Trajectory .ne. 'NO')then
               write(6,list) X(i), Y(i), Z(i)
            endif
         end do
      else
         do  i=1, AN
            if(Trajectory .ne. 'NO')then
               if(minTrajOut .ne. 'YES')then
                  write(6,traj_fmt) X(i), Y(i), Z(i)
               else
                  write(6,traj_fmt_min) X(i), Y(i), Z(i)
               endif
            endif
         end do
      endif

        do i=1, AN
        write(61,list2) i, Vx(i), Vy(i), Vz(i)
        end do

         write(61,*) MDT
         do i=1,MDT
            write(61,*) ChainLength(i)
         enddo

      close(61, status = 'KEEP')

      open(55, file = output, access = 'SEQUENTIAL',	
     Q   status = 'unknown')	

	if (startt .eq. 1)then
	write(55,*) 'The values used for momenta for this run were 
     Q generated randomly.  The position was read from file', filename  
	else 

	write(55,*) 'The position and momenta values used for this 
     Q run were read from file ' , filename
	endif

	write(55,*) 'It has been', writecount, 'of ', stepstop, 'steps'
	write(55,*) 'The leap-frog method was used in this simulation'
        write(55,*) 'The current temperature is', temp

	if (thermint .ne. 0)then
	write(55,*) 'The desired temperature was', T
	endif
	write(55,*) 'The tau being used is', tau
        write(55,*) 'Conformation is ', conf

	write(55,*) 'The total number of possible 2-body 
     Q interactions is', NC
        write(55,*) 'The total number of possible 3-body
     Q interactions is', tripi


	close (55, status = 'KEEP')	


	END

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^End of writetofile^^^^^^^^^^^^^^^^^^^^^^
