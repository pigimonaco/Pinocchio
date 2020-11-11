 # *****************************************************************
 # *                        PINOCCHIO  V4.1                        *
 # *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 # *****************************************************************
 
 # This code was written by
 # Pierluigi Monaco
 # Dipartimento di Fisica, Universita` di Trieste
 # Copyright (C) 2016
 
 # web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 # The original code was developed with:
 # Tom Theuns           Institute for Computational Cosmology, University of Durham 
 # Giuliano Taffoni     INAF - Osservatorio Astronomico di Trieste

 # This program is free software; you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation; either version 2 of the License, or
 # (at your option) any later version.
 
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 
 # You should have received a copy of the GNU General Public License
 # along with this program; if not, write to the Free Software
 # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


"""
Routine for reading Pinocchio's binary catalogs and PLCs.
Usage:

# PLC
import ReadPinocchio as rp
myplc = rp.plc("pinocchio.example.plc.out")
print myplc.redshift

# CATALOGS
import ReadPinocchio as rp
mycat = rp.catalog("pinocchio.0.0000.example.catalog.out")
print mycat.Mass

# HISTORIES
import ReadPinocchio as rp
myhist = rp.histories("pinocchio.example.histories.out")
print myhist.name

Written by Emiliano Munari

LAST MODIFIED:
Pierlugi Monaco 15/04/2016

"""


import numpy as np
import os
import sys
from struct import unpack



class catalog:
    def __init__(self,filename):
        if (not os.path.exists(filename)):
            print "file not found:", filename
            sys.exit()

        print 'opening file '+filename
        f = open(filename,'rb')

        Ngroups = 0L

        dummy = f.read(4)
        NTasksPerFile = np.fromfile(f,dtype=np.int32,count=1)[0]
        NSlices = np.fromfile(f,dtype=np.int32,count=1)[0]
        dummy = f.read(4)

        print 'This file has been written by ',NTasksPerFile,' tasks'
        if NSlices>10:
            lite=True
            record_length=NSlices
            NSlices=1
            print 'This is light output format, record length: %d'%record_length
        else:
            lite=False
            record_length=-1
            if NSlices==1:
                print 'The box has been fragmented in 1 slice'
            else:
                print 'The box has been fragmented in %d slices'%NSlices

        if lite:

            # Count the number of halos
            for iproc in xrange(NTasksPerFile):
                dummy = f.read(4)
                ngood = np.fromfile(f,dtype=np.int32,count=1)[0]
                print ' +++ found %d halos...'%ngood
                dummy = f.read(4)

                Ngroups += ngood
                f.seek(ngood*record_length+8,1)

            # Go back to the starting point (NTasksPerFile already read)
            f.seek(16)

            print 'Total number of halos: ',Ngroups
            print '+++++++++++++++++++'

            self.name = np.empty(Ngroups,dtype=np.uint64)
            self.Mass = np.empty(Ngroups,dtype=np.float32)
            self.pos = np.empty((Ngroups,3),dtype=np.float32)
            self.vel = np.empty((Ngroups,3),dtype=np.float32)

            if record_length == 40:
                record_dtype = np.dtype( [ ( 'name' , np.uint64 ), \
                                           ( 'Mass' , np.float32 ), \
                                           ( 'pos' , np.float32,3 ), \
                                           ( 'vel' , np.float32,3 ), \
                                           ( 'pad' , np.int32 ) ])
            elif record_length == 36:
                record_dtype = np.dtype( [ ( 'name' , np.uint64 ), \
                                           ( 'Mass' , np.float32 ), \
                                           ( 'pos' , np.float32,3 ), \
                                           ( 'vel' , np.float32,3 ) ])
            else:
                print "I do not recognise the record length!"
                sys.exit(1)


            startid = 0
            stopid = 0
            for iproc in xrange(NTasksPerFile):
                dummy = f.read(4)
                ngood = np.fromfile(f,dtype=np.int32,count=1)[0]
                dummy = f.read(4)
                print 'Reading %d halos...'%ngood
                stopid += ngood
                dummy = f.read(4)
                catalog = np.fromfile(f,dtype=record_dtype,count=ngood)
                dummy = f.read(4)

                self.name[startid:stopid] = catalog['name']            
                self.Mass[startid:stopid] = catalog['Mass']            
                self.pos[startid:stopid] = catalog['pos']            
                self.vel[startid:stopid] = catalog['vel']            
                del catalog
                startid = stopid


        else:

            # Count the number of halos
            record_length=None
            for islice in xrange(NSlices):
                for iproc in xrange(NTasksPerFile):
                    dummy = f.read(4)
                    ngood = np.fromfile(f,dtype=np.int32,count=1)[0]
                    print ' +++ found %d halos...'%ngood
                    dummy = f.read(4)
                    if record_length is None:
                        record_length = np.fromfile(f,dtype=np.int32,count=1)[0]
                        f.seek(-4,1)
                    Ngroups += ngood
                    f.seek(ngood*(record_length+8),1)

            # Go back to the starting point (NTasksPerFile already read)
            f.seek(16) 

            print 'Total number of halos: ',Ngroups
            print '+++++++++++++++++++'

            self.name = np.empty(Ngroups,dtype=np.uint64)
            self.Mass = np.empty(Ngroups,dtype=np.float64)
            self.posin = np.empty((Ngroups,3),dtype=np.float64)
            self.pos = np.empty((Ngroups,3),dtype=np.float64)
            self.vel = np.empty((Ngroups,3),dtype=np.float64)
            self.Npart = np.empty(Ngroups,dtype=np.int64)

            record_dtype = np.dtype( [     ( 'dummy' , np.int32 ), \
                                           ( 'name' , np.uint64 ), \
                                           ( 'Mass' , np.float64 ), \
                                           ( 'posin' , np.float64,3 ), \
                                           ( 'pos' , np.float64,3 ), \
                                           ( 'vel' , np.float64,3 ), \
                                           ( 'Npart' , np.int32 ),\
                                           ( 'pad' , np.int32 ),\
                                           ( 'dummy2' , np.int32 )])

            startid = 0
            stopid = 0
            for islice in xrange(NSlices):
                for iproc in xrange(NTasksPerFile):
                    dummy = f.read(4)
                    ngood = np.fromfile(f,dtype=np.int32,count=1)[0]
                    dummy = f.read(4)
                    print 'Reading %d halos...'%ngood
                    stopid += ngood
                    catalog = np.fromfile(f,dtype=record_dtype,count=ngood)

                    self.name[startid:stopid] = catalog['name']            
                    self.Mass[startid:stopid] = catalog['Mass']            
                    self.posin[startid:stopid] = catalog['posin']            
                    self.pos[startid:stopid] = catalog['pos']            
                    self.vel[startid:stopid] = catalog['vel']            
                    self.Npart[startid:stopid] = catalog['Npart']  
                    del catalog
                    startid = stopid

        f.close()
    
    

class plc:
    def __init__(self,filename,memsplit=10):
        if (not os.path.exists(filename)):
            print "file not found:", filename
            sys.exit()

        print 'opening file '+filename
        f=open(filename,'rb')

        header = np.fromfile(f,dtype=np.int32,count=3)

        lite=(header[0]==4)

        if lite:

            end=os.path.getsize(filename)

            record_length = header[1]

            print 'This is light output format'
            if record_length != 40:
                print "I do not recognise the record length!"
                sys.exit(1)

            record_dtype = np.dtype( [ ( 'name' , np.uint64 ), \
                                       ( 'Mass' , np.float32 ), \
                                       ( 'theta' , np.float32 ), \
                                       ( 'phi' , np.float32 ), \
                                       ( 'obs_redshift' , np.float32 ),  \
                                       ( 'true_redshift' , np.float32 ),  \
                                       ( 'vel' , np.float32,3 ) ] ) 

            # Count the number of halos
            Ngroups = 0L

            while f.tell()<end:
                nmore = np.fromfile(f,dtype=np.int32,count=3)[1]
                print ' +++ found %d halos...'%nmore
                Ngroups += nmore

                f.seek(nmore * record_length + 8, 1)

            print 'Total number of halos: ',Ngroups
            print '+++++++++++++++++++'

            self.name = np.empty(Ngroups,dtype=np.uint32)
            self.Mass = np.empty(Ngroups,dtype=np.float32)
            self.theta = np.empty(Ngroups,dtype=np.float32)
            self.phi = np.empty(Ngroups,dtype=np.float32)
            self.obs_redshift = np.empty(Ngroups,dtype=np.float32)
            self.true_redshift = np.empty(Ngroups,dtype=np.float32)
            self.vel = np.empty((Ngroups,3),dtype=np.float32)

            f.seek(12)
            startid = 0
            while f.tell()<end:
                nmore = np.fromfile(f,dtype=np.int32,count=3)[1]
                stopid = startid+nmore
                print 'Reading %d halos...'%nmore

                dummy = f.read(4)
                plc = np.fromfile(f,dtype=record_dtype,count=stopid-startid)
                dummy = f.read(4)

                self.name[startid:stopid] = plc['name']
                self.Mass[startid:stopid] = plc['Mass']
                self.theta[startid:stopid] = plc['theta']
                self.phi[startid:stopid] = plc['phi']
                self.obs_redshift[startid:stopid] = plc['obs_redshift']
                self.true_redshift[startid:stopid] = plc['true_redshift']
                self.vel[startid:stopid] = plc['vel']

                del plc
                startid = stopid

        else:
            
            f.seek(0)

            Ngroups = os.path.getsize(filename)/(13*8+4+4)
            print "reading %d halos"%Ngroups

            self.name = np.empty(Ngroups,dtype=np.uint64)
            self.redshift = np.empty(Ngroups,dtype=np.float64)
            self.pos = np.empty((Ngroups,3),dtype=np.float64)
            self.vel = np.empty((Ngroups,3),dtype=np.float64)
            self.Mass = np.empty(Ngroups,dtype=np.float64)
            self.theta = np.empty(Ngroups,dtype=np.float64)
            self.phi = np.empty(Ngroups,dtype=np.float64)
            self.vlos = np.empty(Ngroups,dtype=np.float64)
            self.obsz = np.empty(Ngroups,dtype=np.float64)

            record_dtype = np.dtype( [ ( 'dummy' , np.int32 ), \
                                       ( 'name' , np.uint64 ), \
                                       ( 'redshift' , np.float64 ),  \
                                       ( 'pos' , np.float64,3 ), \
                                       ( 'vel' , np.float64,3 ), \
                                       ( 'Mass' , np.float64 ), \
                                       ( 'theta' , np.float64 ), \
                                       ( 'phi' , np.float64 ), \
                                       ( 'vlos' , np.float64 ), \
                                       ( 'obsz' , np.float64 ), \
                                       ( 'dummy2' , np.int32 ) ] )        

            Npartial = Ngroups/memsplit
            startid = 0
            stopid = Npartial
            for i in xrange(memsplit):
                if (i == memsplit-1):
                    stopid = stopid + Ngroups % memsplit
            
                plc = np.fromfile(f,dtype=record_dtype,count=stopid-startid)            
                self.name[startid:stopid] = plc['name']
                self.redshift[startid:stopid] = plc['redshift']
                self.pos[startid:stopid] = plc['pos']
                self.vel[startid:stopid] = plc['vel']
                self.Mass[startid:stopid] = plc['Mass']
                self.theta[startid:stopid] = plc['theta']
                self.phi[startid:stopid] = plc['phi']
                self.vlos[startid:stopid] = plc['vlos']
                self.obsz[startid:stopid] = plc['obsz']

                del plc
                startid = stopid
                stopid += Npartial

        f.close()


        
class histories:
    def __init__(self,filename):
        if (not os.path.exists(filename)):
            print "file not found:", filename
            sys.exit()

        f = open(filename,'rb')

        Ngroups = 0L

        header = np.fromfile(f,dtype=np.int32,count=3)
        NSlices = header[1]

        if NSlices>10:
            lite=True
            record_length=NSlices
            NSlices=1
            print 'This is light output format, record length: %d'%record_length

            record_dtype = np.dtype( [  ( 'dummy'          , np.int32 ), \
                                        ( 'name'           , np.uint64 ), \
                                        ( 'nickname'       , np.int32 ), \
                                        ( 'link'           , np.int32 ), \
                                        ( 'merged_with'    , np.int32 ), \
                                        ( 'mass_at_merger' , np.int32 ), \
                                        ( 'mass_of_main'   , np.int32 ), \
                                        ( 'z_merging'      , np.float32 ), \
                                        ( 'z_peak'         , np.float32 ), \
                                        ( 'z_appear'       , np.float32 ), \
                                        ( 'dummy2'         , np.int32 ) ] )        
        else:
            lite=False
            record_length=-1
            if NSlices==1:
                print 'The box has been fragmented in 1 slice'
            else:
                print 'The box has been fragmented in %d slices'%NSlices

            record_dtype = np.dtype( [  ( 'dummy'          , np.int32 ), \
                                        ( 'name'           , np.uint64 ), \
                                        ( 'nickname'       , np.int32 ), \
                                        ( 'link'           , np.int32 ), \
                                        ( 'merged_with'    , np.int32 ), \
                                        ( 'mass_at_merger' , np.int32 ), \
                                        ( 'mass_of_main'   , np.int32 ), \
                                        ( 'z_merging'      , np.float64 ), \
                                        ( 'z_peak'         , np.float64 ), \
                                        ( 'z_appear'       , np.float64 ), \
                                        ( 'pad'            , np.int32 ) , \
                                        ( 'dummy2'         , np.int32 ) ] )        

       # Count the number of halos
        Total=0
        for islice in xrange(NSlices):

            header = np.fromfile(f,dtype=np.int32,count=4)
            Ntrees = header[1]
            Nbranches = header[2]
            Total+=Nbranches
            
            print 'Slice N. ',islice,': ',Ntrees,' trees and ',Nbranches,' branches'

            for itree in xrange(Ntrees):
                header = np.fromfile(f,dtype=np.int32,count=4)
                mytree = header[1]
                mynbranch = header[2]

                f.seek(mynbranch*record_dtype.itemsize,1)

        # Go back to the starting point (NTasksPerFile already read)
        print '+++++++++++++++++++'
        f.seek(12)

        print 'Total number of halos: ',Total

        self.name           = np.empty(Total,dtype=np.uint64)
        self.nickname       = np.empty(Total,dtype=np.int32 )
        self.link           = np.empty(Total,dtype=np.int32 )
        self.merged_with    = np.empty(Total,dtype=np.int32 )
        self.mass_at_merger = np.empty(Total,dtype=np.int32 )
        self.mass_of_main   = np.empty(Total,dtype=np.int32 )
        self.z_merging      = np.empty(Total,dtype=np.float64)
        self.z_appear       = np.empty(Total,dtype=np.float64)
        self.z_peak         = np.empty(Total,dtype=np.float64)

        startid = 0
        stopid = 0
        for islice in xrange(NSlices):
            header = np.fromfile(f,dtype=np.int32,count=4)
            Ntrees = header[1]
            Nbranches = header[2]

            for itree in xrange(Ntrees):
                header = np.fromfile(f,dtype=np.int32,count=4)
                mytree = header[1]
                mynbranch = header[2]

                stopid += mynbranch
                catalog = np.fromfile(f,dtype=record_dtype,count=mynbranch)

                self.name          [startid:stopid] = catalog['name']
                self.nickname      [startid:stopid] = catalog['nickname']
                self.link          [startid:stopid] = catalog['link']
                self.merged_with   [startid:stopid] = catalog['merged_with']
                self.mass_at_merger[startid:stopid] = catalog['mass_at_merger']
                self.mass_of_main  [startid:stopid] = catalog['mass_of_main']
                self.z_merging     [startid:stopid] = catalog['z_merging']
                self.z_appear      [startid:stopid] = catalog['z_appear']
                self.z_peak        [startid:stopid] = catalog['z_peak']
                del catalog
                startid = stopid

        f.close()

