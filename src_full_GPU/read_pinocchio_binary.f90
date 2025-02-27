! *****************************************************************
! *                        PINOCCHIO  V4.1                        *
! *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
! *****************************************************************
! 
! This code was written by
! Pierluigi Monaco
! Copyright (C) 2016
! 
! web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html
! 
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!


   
! This code reads a catalog written by pinocchio in binary format 
! and checks its integrity

! to compile the code:
!    ifort -O -o read_pinocchio_binary read_pinocchio_binary.f90


program read_pinocchio

  implicit none
  real(kind=8)::z,q(3),x(3),v(3),M,theta,phi,v_los,obsz,zme, zpe, zap
  integer(kind=8)::id,n,name,ntrees_tot,nbranch_tot
  logical::First
  integer(kind=4)::nproc,ngood,igood,Npart,NSlices,nblocks,ThisSlice,ntrees,nbranch,nbranch_tree,itree,ibranch,thistree,nick, ll, mw, mass, mam
  character(len=100)::runname,filename,redshift

  write(*,'(a)')'This is an example code to read pinocchio catalog in binary format'
  write(*,'(a)')'Pinocchio generates catalogs at fixed redshift in files'
  write(*,'(a)')'like pinocchio.0.0000.SecondEuclidMocks0001.catalog.out,'
  write(*,'(a)')'where SecondEuclidMocks0001 is the run name'
  write(*,'(a)')'and files like pinocchio.SecondEuclidMocks0001.plc.out'
  write(*,'(a)')'for the catalog on the past-light-cone'
  write(*,'(a)')'or pinocchio.SecondEuclidMocks0001.histories.out'
  write(*,'(a)')'for the halo merger trees'
  write(*,'(a)')
  write(*,'(a)')'  Usage: read_pinocchio_binary runname redshift'
  write(*,'(a)')'     runname is the name of the run for all the files'
  write(*,'(a)')'     redshift must be equal to the redshift given in the file name'
  write(*,'(a)')'     to read the past-light-cone enter the string "plc" in place of redshift'
  write(*,'(a)')'     to read the histories file enter string "hist" in place of redshift'
  write(*,'(a)')
 
  call getarg(1,runname)
  call getarg(2,redshift)

  write(*,'("runname: ",a)')runname
  write(*,'("redshift: ",a)')redshift

  if (runname=='' .or. redshift=='') then
     stop
  endif

  if (trim(redshift) == "plc") then

     ! CATALOG ON THE PAST LIGHT CONE
     filename="pinocchio."//trim(runname)//".plc.out"
     write(*,'("Opening file ",a)')trim(filename)
     write(*,'("I will write out the first ten halos")')

     open(10,file=filename,status='old',form='unformatted',err=20)
     n=0
     do while(.true.)     
        read(10,end=2,err=30)id,z,x,v,M,theta,phi,v_los,obsz
        if (n<10) write(*,'(1x,i12,7(1x,f16.6),1x,e15.8,4(1x,f16.6))')id,z,x,v,M,theta,phi,v_los,obsz
        n=n+1
     enddo     
     
2    close(10)
     write(*,'("I found ",i," halos in the file")')n
     write(*,'("This is the last halo in the catalog:")')
     write(*,'(1x,i12,7(1x,f16.6),1x,e15.8,4(1x,f16.6))')id,z,x,v,M,theta,phi,v_los,obsz

  else if (trim(redshift) == "hist") then

     ! MERGER HISTORIES
     filename="pinocchio."//trim(runname)//".histories.out"
     write(*,'("Opening file ",a)')trim(filename)
     write(*,'("I will write out the first tree")')

     open(10,file=filename,status='old',form='unformatted',err=15)
     read(10,err=99)NSlices
     write(*,'("The run used ",i," slices")')NSlices

     ntrees_tot=0
     nbranch_tot=0
     First=.true.
     do ThisSlice=1,NSlices
        read(10,err=99)ntrees,nbranch
        ntrees_tot=ntrees_tot+ntrees
        nbranch_tot=nbranch_tot+nbranch

        do itree=1,ntrees
           read(10,err=99)thistree,nbranch_tree

           do ibranch=1,nbranch_tree
              read(10,err=99)name,nick,ll,mw,mass,mam,zme,zpe,zap

              if (First) then
                 write(*,'(1x,i12,1x,i6,1x,i6,1x,i6,1x,i9,1x,i9,1x,f9.4,1x,f9.4,1x,f9.4)')name,nick,ll,mw,mass,mam,zme,zpe,zap
              endif

           enddo
           First=.false.

        enddo
     enddo

     write(*,'("I found ",i," trees and ",i," branches in the file")')ntrees_tot,nbranch_tot
     write(*,'("This is the last branch in the catalog:")')
     write(*,'(1x,i12,1x,i6,1x,i6,1x,i6,1x,i9,1x,i9,1x,f9.4,1x,f9.4,1x,f9.4)')name,nick,ll,mw,mass,mam,zme,zpe,zap

     close(10)

   else

     ! CATALOG AT FIXED REDSHIFT
     filename="pinocchio."//trim(redshift)//"."//trim(runname)//".catalog.out"
     write(*,'("Opening file ",a)')trim(filename)
     write(*,'("I will write out the first ten halos")')

     open(10,file=filename,status='old',form='unformatted',err=10)
     n=0

     read(10,err=44)nproc,NSlices
     write(*,'("The file has been written by ",i," tasks")')nproc
     write(*,'("The box has been fragmented in ",i," slices")')NSlices
     go to 45
44   rewind(10)
     read(10,err=99)nproc
     write(*,'("The file has been written by ",i," tasks")')nproc
     write(*,'("This old version catalog does not give the number of slices used")')
     NSlices=-1
45   continue

     nblocks=0
     do while (.true.)
        read(10,err=99,end=21)ngood
        n=n+ngood
        nblocks=nblocks+1
        write(*,'("  found ",i," halos in block N. ",i)')ngood,nblocks
        do igood=1,ngood
           read(10)id, M, q, x, v, Npart
           if (nblocks==1 .and. igood<=10) &
                write(*,'(1x,i12,1x,e13.6,9(1x,f10.2),1x,i12)') id, M, q, x, v, Npart 
        enddo
     enddo

21   close(10)
     write(*,'("I found ",i," halos in the file")')n
     write(*,'("This is the last halo in the catalog:")')
     write(*,'(1x,i12,1x,e13.6,9(1x,f10.2),1x,i12)') id, M, q, x, v, Npart 

     if (NSlices<0) then
        write(*,'("   nproc=",i,",  nblocks=",i,",  reconstructed nslices=",i)')nproc,nblocks,nblocks/nproc
        if (mod(nblocks,nproc) /= 0) then
           write(*,'("WARNING: the numbers of blocks and tasks are not consistent!")')
        endif
     else
        write(*,'("   nproc=",i,",  nblocks=",i,",  read nslices=",i)')nproc,nblocks,NSlices
	if (nblocks /= nproc*NSlices) then
           write(*,'("WARNING: the number of blocks is not as expected, nproc * nslices = ",i)')nproc*NSlices
        endif

     endif
  endif

  stop

10 write(*,'("Error: catalog file ",a," not found")')trim(filename)
  stop

99 write(*,'("Error in reading catalog file ",a)')trim(filename)
  stop

20 write(*,'("Error: PLC file ",a," not found")')trim(filename)
  stop

30 write(*,'("Error in reading PLC file",a)')trim(filename)
  stop

15 write(*,'("Error: history file ",a," not found")')trim(filename)
  stop

end program read_pinocchio
