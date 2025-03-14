# HEADER

"""
Routine for reading Pinocchio's binary catalogs and PLCs.
Usage:

# CATALOGS
import ReadPinocchio5 as rp
mycat = rp.catalog("pinocchio.0.0000.example.catalog.out")
print(mycat.Mass)

# PLC
import ReadPinocchio5 as rp
myplc = rp.plc("pinocchio.example.plc.out")
print(myplc.redshift)

# HISTORIES
import ReadPinocchio5 as rp
myhist = rp.histories("pinocchio.example.histories.out")
print(myhist.name)

Written by Pierluigi Monaco, Matteo Biagetti and Emiliano Munari

LAST MODIFIED:
Pierlugi Monaco 30/12/22

"""


import numpy as np
import os
import sys
import copy
import struct


VERBOSE=False

class catalog:

    def __init__(self,filename,silent=False,first_file=None,last_file=None):

        if VERBOSE:
            silent=False

        # checks that the filename contains 'catalog'
        if not 'catalog' in filename:
            print("Are you sure you are providing the right file name?")
            if 'plc' in filename:
                print("...this looks like a plc file, please read it with rp.plc(filename)")
            elif 'histories' in filename:
                print("...this looks like a histories file, please read it with rp.histories(filename)")
            return None

        # checks that the input file ends by ".out"
        last_ext=filename.rfind('.')
        if filename[last_ext:]!='.out':
            print("The catalog file should end with .out, the file number extension will be checked by the code")
            return None

        # checks that the file exists, or that there are multiple files, and in case count them
        if not os.path.exists(filename):

            if not os.path.exists(filename+'.0'):

                print("file {} or {} not found:".format(filename,filename+'.0'))
                return None

            else:

                Nfiles=1
                while os.path.exists(filename+'.{}'.format(Nfiles)):
                    Nfiles+=1
                if not silent:
                    print("The catalog is written in {} files".format(Nfiles))

        else:

            Nfiles=1
            if not silent:
                print("The catalog is written in 1 file")

        # opens the (first) file and reads the record length
        if Nfiles==1:
            if not silent:
                print('reading header of file '+filename)
            reading = np.fromfile(filename,dtype=np.int32,count=10)
        else:
            if not silent:
                print('reading header of file '+filename+'.0')
            reading = np.fromfile(filename+'.0',dtype=np.int32,count=10)

        # number of tasks that write into a single file
        NTasksPerFile = reading[1]
        if not silent:
            print('This file has been written by {} tasks'.format(NTasksPerFile))

        # the header gives either the number of slices (always < 10) or the record length
        # in the first case the record length can be read from the 8th integer
        if reading[2]>10:
            newRun=True
            record_length=reading[2]
            Nslices=1
            if not silent:
                print('This is new output format, record length: {}'.format(record_length))
        else:
            newRun=False
            record_length=reading[7]
            Nslices=reading[2]
            if not silent:
                print('This is classic output format, record length: {}'.format(record_length))
                if Nslices==1:
                    print('The box has been fragmented in 1 slice')
                else:
                    print('The box has been fragmented in {} slices'.format(Nslices))


        # sets the record
        if record_length==96:          # this is the classic format in double precision

            self.cat_dtype=[ ('name',  np.int64),
                             ('Mass',  np.float64),
                             ('posin', np.float64, 3),
                             ('pos',   np.float64, 3),
                             ('vel',   np.float64, 3),
                             ('npart', np.int32) ]

            stored_dtype  =[ ('fort',  np.int32),
                             ('name',  np.int64),
                             ('Mass',  np.float64),
                             ('posin', np.float64, 3),
                             ('pos',   np.float64, 3),
                             ('vel',   np.float64, 3),
                             ('npart', np.int32),
                             ('pad'  , np.int32),
                             ('trof',  np.int32) ]


        elif record_length==56:

            if newRun:                     # the new format has posin in a different position

                self.cat_dtype=[ ('name',  np.int64),
                                 ('Mass',  np.float32),
                                 ('pos',   np.float32, 3),
                                 ('vel',   np.float32, 3),
                                 ('posin', np.float32, 3),
                                 ('npart', np.int32) ]

                stored_dtype  =[ ('name',  np.int64),
                                 ('Mass',  np.float32),
                                 ('pos',   np.float32, 3),
                                 ('vel',   np.float32, 3),
                                 ('posin', np.float32, 3),
                                 ('npart', np.int32),
                                 ('pad'  , np.int32) ]

            else:                          # this is the classic format in single precision

                self.cat_dtype=[ ('name',  np.int64),
                                 ('Mass',  np.float32),
                                 ('posin', np.float32, 3),
                                 ('pos',   np.float32, 3),
                                 ('vel',   np.float32, 3),
                                 ('npart', np.int32) ]

                stored_dtype  =[ ('fort',  np.int32),
                                 ('name',  np.int64),
                                 ('Mass',  np.float32),
                                 ('posin', np.float32, 3),
                                 ('pos',   np.float32, 3),
                                 ('vel',   np.float32, 3),
                                 ('npart', np.int32),
                                 ('pad'  , np.int32),
                                 ('trof',  np.int32) ]

        elif record_length==48:        # this is the new light format

            self.cat_dtype=[ ('name',  np.int64),
                             ('Mass',  np.float32),
                             ('pos',   np.float32, 3),
                             ('vel',   np.float32, 3),
                             ('posin', np.float32, 3) ]

            stored_dtype = self.cat_dtype


        elif record_length==40:        # this was used in NewClusterMocks

            self.cat_dtype=[ ('name',  np.int64),
                             ('Mass',  np.float32),
                             ('pos',   np.float32, 3),
                             ('vel',   np.float32, 3),
                             ('npart', np.int32) ]

            stored_dtype  =[ ('fort',  np.int32),
                             ('name',  np.int64),
                             ('Mass',  np.float32),
                             ('pos',   np.float32, 3),
                             ('vel',   np.float32, 3),
                             ('npart', np.int32),
                             ('trof',  np.int32) ]

        else:
            print("sorry, I do not recognize this record length")
            return None

        # decides what files to read
        if Nfiles>1:
            if first_file is None:
                first_file=0
            elif first_file<0:
                first_file=0
            elif first_file>Nfiles:
                first_file=Nfiles
            if last_file is None:
                last_file=Nfiles
            else:
                last_file += 1      
                if last_file<first_file:
                    last_file=first_file+1
                elif last_file>Nfiles:
                    last_file=Nfiles
            # this is to be used to define a pythonic range
            if not silent:
                print("I will read files in the python range from {} to {}".format(first_file,last_file))
        else:
            first_file=0
            last_file=Nfiles

        # prepares to read the file(s)
        self.data=None
        NhalosPerFile=np.zeros(Nfiles,dtype=np.int64)

        # loops on the files to be read
        for myfile in range(first_file,last_file):

            if Nfiles==1:
                myfname = filename
            else:
                myfname = filename+'.{}'.format(myfile)

            if VERBOSE:
                print("reading file {}".format(myfname))

            # reads the file as a binary object
            with open(myfname,'rb') as f:
                bindata=f.read()
                f.close()
            filesize = len(bindata)

            # reads the number of halos contained in each header
            Nblocks=NTasksPerFile * Nslices     # this is the number of blocks to read
            Nwritten=0                          # some tasks may have no halos to write
            pos=16
            vechalo=[]
            while pos < filesize:
                # reads the number of halos to read
                vec=struct.Struct('iii').unpack(bindata[pos:pos+12])[1]
                pos+=12
                vechalo.append(vec)
                if vec>0:
                    if newRun:
                        pos += 8+vec*record_length
                    else:
                        pos += vec*(record_length+8)

            vechalo = np.asarray(vechalo)
            #print(vechalo.size, Nblocks, vechalo)
            NhalosPerFile[myfile]=vechalo.sum()
            Nwritten=(vechalo>0).sum()

            if VERBOSE:
                print(f"Found {Nwritten} non-void blocks in this file over {Nblocks}")

            # checks that the lenght of data is as expected
            if newRun:
                file_length = 16 + Nblocks*12 + Nwritten*8 + NhalosPerFile[myfile]*record_length
            else:
                file_length = 16 + Nblocks*12 + NhalosPerFile[myfile]*(record_length+8)

            if file_length != filesize:
                print(f'ERROR: inconsistency in the file size, should be {file_length} but I found {filesize} bytes')
                return None
            elif VERBOSE:
                print(f'predicted file length {file_length} matches with file size {filesize}')

            # format to select information from the byte object
            if newRun:
                cleanForm='16x '
                for block in range(Nblocks):
                    if vechalo[block]>0:
                        cleanForm+='16x {}s 4x '.format(vechalo[block]*record_length)
                    else:
                        cleanForm+='12x '

            else:
                cleanForm='16x '
                for block in range(Nblocks):
                    if vechalo[block]>0:
                        cleanForm+='12x {}s '.format(vechalo[block]*(record_length+8))
                    else:
                        cleanForm+='12x '

            # removes all unwanted information from the binary structure
            try:
                cleaned = b''.join(struct.unpack(cleanForm, bindata))
            except:
                print("ERROR: I do not recognise the data structure!")
                return None
            del bindata

            # reads the catalog from the cleaned bynary structure
            thiscat = np.frombuffer(cleaned, dtype=stored_dtype)
            del cleaned

            # removes unwanted columns from the catalog
            if self.data is None:
                self.data = np.zeros(NhalosPerFile[myfile], dtype=self.cat_dtype)
            else:
                self.data.resize(self.data.shape[0]+NhalosPerFile[myfile])

            for name in self.data.dtype.names:
                self.data[name][-NhalosPerFile[myfile]:]=thiscat[name]
            del thiscat

            if not silent:
                print("done with file {}".format(myfname))

        if not silent:
            print("Reading catalog done, {} groups found".format(NhalosPerFile.sum()))

        # Create few pointers to make it compatible
        self.Mass  = self.data['Mass']
        self.pos   = self.data['pos']
        self.Npart = self.data['npart']
        self.vel   = self.data['vel']

class plc:

    def __init__(self,filename,silent=False,first_file=None,last_file=None,onlyNfiles=False):

        if VERBOSE:
            silent=False

        # checks that the filename contains 'plc'
        if not 'plc' in filename:
            print("Are you sure you are providing the right file name?")
            if 'catalog' in filename:
                print("...this looks like a catalog file, please read it with rp.catalog(filename)")
            elif 'histories' in filename:
                print("...this looks like a histories file, please read it with rp.histories(filename)")
            return None

        # checks that the input file ends by ".out"
        last_ext=filename.rfind('.')
        if filename[last_ext:]!='.out':

            print("The catalog file should end with .out, the file number extension will be checked by the code")
            return None

        # checks that the file exists, of that there are multiple files
        if not os.path.exists(filename):

            if not os.path.exists(filename+'.0'):

                print("file {} or {} not found:".format(filename,filename+'.0'))
                return None

            else:

                Nfiles=1
                while os.path.exists(filename+'.{}'.format(Nfiles)):
                    Nfiles+=1
                if not silent:
                    print("The catalog is written in {} files".format(Nfiles))

        else:

            Nfiles=1
            if not silent:
                print("The catalog is written in 1 file")

        self.Nfiles=Nfiles
        if onlyNfiles:
            return

        # opens the (first) file and reads the record length
        if Nfiles==1:
            if not silent:
                print('reading header of file '+filename)
            reading = np.fromfile(filename,dtype=np.int32,count=3)
        else:
            if not silent:
                print('reading header of file '+filename+'.0')
            reading = np.fromfile(filename+'.0',dtype=np.int32,count=3)

        # reads the record length
        if reading[0]==4:
            record_length=reading[1]
            newRun=True
            if not silent:
                if record_length==32:
                    print('This is new light output format, record length: {}'.format(record_length))
                else:
                    print('This is new full output format, record length: {}'.format(record_length))
        else:
            record_length=reading[0]
            newRun=False
            if not silent:
                print('This is classic output format, record length: {}'.format(record_length))
                  
        # sets the record
        if record_length==104:

            self.cat_dtype = [ ( 'name'  , np.uint64 ),
                               ( 'truez' , np.float64 ),
                               ( 'pos'   , np.float64,3 ),
                               ( 'vel' , np.float64,3 ),
                               ( 'Mass' , np.float64 ),
                               ( 'theta' , np.float64 ),
                               ( 'phi' , np.float64 ),
                               ( 'vlos' , np.float64 ),
                               ( 'obsz' , np.float64 ) ]

            stored_dtype   = [ ('fort',  np.int32),
                               ( 'name'  , np.uint64 ),
                               ( 'truez' , np.float64 ),
                               ( 'pos'   , np.float64,3 ),
                               ( 'vel' , np.float64,3 ),
                               ( 'Mass' , np.float64 ),
                               ( 'theta' , np.float64 ),
                               ( 'phi' , np.float64 ),
                               ( 'vlos' , np.float64 ),
                               ( 'obsz' , np.float64 ) ,
                               ('trof',  np.int32) ]

        elif record_length==56:

            self.cat_dtype = [ ( 'name'  , np.uint64 ),
                               ( 'truez' , np.float32 ),
                               ( 'pos'   , np.float32,3 ),
                               ( 'vel' , np.float32,3 ),
                               ( 'Mass' , np.float32 ),
                               ( 'theta' , np.float32 ),
                               ( 'phi' , np.float32 ),
                               ( 'vlos' , np.float32 ),
                               ( 'obsz' , np.float32 ) ]
                  
            stored_dtype   = [ ('fort',  np.int32),
                               ( 'name'  , np.uint64 ),
                               ( 'truez' , np.float32 ),
                               ( 'pos'   , np.float32,3 ),
                               ( 'vel' , np.float32,3 ),
                               ( 'Mass' , np.float32 ),
                               ( 'theta' , np.float32 ),
                               ( 'phi' , np.float32 ),
                               ( 'vlos' , np.float32 ),
                               ( 'obsz' , np.float32 ) ,
                               ('trof',  np.int32) ]

        elif record_length==32:

            self.cat_dtype = [ ( 'name'  , np.uint64 ),
                               ( 'truez' , np.float32 ),
                               ( 'Mass'  , np.float32 ),
                               ( 'theta' , np.float32 ),
                               ( 'phi'   , np.float32 ),
                               ( 'obsz'  , np.float32 ) ]

            stored_dtype   = [ ( 'name'  , np.uint64 ),
                               ( 'truez' , np.float32 ),
                               ( 'Mass'  , np.float32 ),
                               ( 'theta' , np.float32 ),
                               ( 'phi'   , np.float32 ),
                               ( 'obsz'  , np.float32 ),
                               ( 'pad'   , np.float32 ) ]

        else:
            print("sorry, I do not recognize this record length")
            return None

        # decides what files to read
        if Nfiles>1:
            if first_file is None:
                first_file=0
            elif first_file<0:
                first_file=0
            elif first_file>Nfiles:
                first_file=Nfiles
            if last_file is None:
                last_file=Nfiles
            else:
                last_file += 1      
                if last_file<first_file:
                    last_file=first_file+1
                elif last_file>Nfiles:
                    last_file=Nfiles
            # this is to be used to define a pythonic range
            if not silent:
                print("I will read files in the python range from {} to {}".format(first_file,last_file))
        else:
            first_file=0
            last_file=Nfiles

        # prepares to read the file(s)
        self.data=None
        NhalosPerFile=np.zeros(Nfiles,dtype=np.int64)

        # loops on the files to be read
        for myfile in range(first_file,last_file):

            if Nfiles==1:
                myfname = filename
            else:
                myfname = filename+'.{}'.format(myfile)

            if VERBOSE:
                print("reading file {}".format(myfname))

            # reads the file as a binary object
            with open(myfname,'rb') as f:
                bindata=f.read()
                f.close()
            filesize = len(bindata)
            
            if newRun:

                # scrolls the file to map the blocks

                if VERBOSE:
                    print("scrolling file to reconstruct the blocks:")
                cleanForm='12x '
                pos=12
                vechalo=[]
                Nblocks=0
                while pos < filesize:
                    # reads the number of halos to read
                    vec=struct.Struct('iii').unpack(bindata[pos:pos+12])[1]
                    pos+=12
                    vechalo.append(vec)
                    if vec==0:
                        print("THIS SHOULD NOT HAPPEN!")
                        return None
                    Nblocks+=1
                    pos += 8+vec*record_length
                    cleanForm+='16x {}s 4x '.format(vec*record_length)

                vechalo = np.asarray(vechalo)
                NhalosPerFile[myfile]=vechalo.sum()

                if VERBOSE:
                    print(f"data are written in {Nblocks} blocks")

                # removes all unwanted information from the binary structure
                try:
                    cleaned = b''.join(struct.unpack(cleanForm, bindata))
                    if VERBOSE:
                        print("cleaning of binary object done")
                except:
                    print("ERROR: I do not recognise the data structure!")
                    return None
                del bindata

                # reads the catalog from the cleaned bynary structure
                thiscat = np.frombuffer(cleaned, dtype=stored_dtype)
                if VERBOSE:
                    print("catalog extracted")
                del cleaned

            else:

                # reads the catalog directly from the binary object
                thiscat = np.frombuffer(bindata, dtype=stored_dtype)
                if VERBOSE:
                    print("catalog extracted")
                del bindata
                NhalosPerFile[myfile]=len(thiscat)

            if newRun:
                file_length = 12 + Nblocks*20 + NhalosPerFile[myfile]*record_length
            else:
                file_length = NhalosPerFile[myfile]*(record_length+8)
            if file_length != filesize:
                print(f'ERROR: inconsistency in the file size, should be {file_length} but I found {filesize} bytes')
                return None
            elif VERBOSE:
                print(f'predicted file length {file_length} matches with file size {filesize}')


            # removes unwanted columns from the catalog
            if self.data is None:
                self.data = np.zeros(NhalosPerFile[myfile], dtype=self.cat_dtype)
            else:
                self.data.resize(self.data.shape[0]+NhalosPerFile[myfile])

            for name in self.data.dtype.names:
                self.data[name][-NhalosPerFile[myfile]:]=thiscat[name]
            del thiscat


            if not silent:
                print("done with file {}".format(myfname))

        if not silent:
            print("Reading plc done, {} groups found".format(NhalosPerFile.sum()))



        
class histories:

    def __init__(self,filename,silent=False,first_file=None,last_file=None):

        if VERBOSE:
            silent=False

        # checks that the filename contains 'catalog'
        if not 'histories' in filename:
            print("Are you sure you are providing the right file name?")
            if 'plc' in filename:
                print("...this looks like a plc file, please read it with rp.plc(filename)")
            elif 'catalog' in filename:
                print("...this looks like a catalog file, please read it with rp.catalog(filename)")
            return None

        # checks that the input file ends by ".out"
        last_ext=filename.rfind('.')
        if filename[last_ext:]!='.out':

            print("The history file should end with .out, the file number extension will be checked by the code")
            return None

        # checks that the file exists, of that there are multiple files
        if not os.path.exists(filename):

            if not os.path.exists(filename+'.0'):

                print("file {} or {} not found:".format(filename,filename+'.0'))
                return None

            else:

                Nfiles=1
                while os.path.exists(filename+'.{}'.format(Nfiles)):
                    Nfiles+=1
                if not silent:
                    print("The catalog is written in {} files".format(Nfiles))

        else:

            Nfiles=1
            if not silent:
                print("The catalog is written in 1 file")

        # opens the (first) file and reads the record length
        if Nfiles==1:
            if not silent:
                print('opening file '+filename)
            reading = np.fromfile(filename,dtype=np.int32,count=12)
            FileLength = os.path.getsize(filename)
        else:
            if not silent:
                print('opening file '+filename+'.0')
            reading = np.fromfile(filename+'.0',dtype=np.int32,count=12)
            FileLength = os.path.getsize(filename+'.0')

        # reads the header and the record length
        if reading[1]>10:
            newRun=True
            record_length=reading[1]
            Nslices=1
            Light = FileLength==np.int64(record_length)*np.int64(reading[5])+np.int64(4*7)
            if not silent:
                if not Light:
                    print('This is new output format, record length: {}'.format(record_length))
                else:
                    print('This is new light output format, record length: {}'.format(record_length))
        else:

            if not silent:
                print("WARNING: reading of V4 histories is very slow")

            newRun=False
            record_length=reading[11]
            Nslices=reading[1]
            if not silent:
                print('This is classic output format, record length: {}'.format(record_length))
                if Nslices==1:
                    print('The box has been fragmented in 1 slice')
                else:
                    print('The box has been fragmented in {} slices'.format(Nslices))


        if record_length==40:

            self.cat_dtype=[ ( 'name'           , np.uint64 ),
                             ( 'nickname'       , np.int32 ),
                             ( 'link'           , np.int32 ),
                             ( 'merged_with'    , np.int32 ),
                             ( 'mass_at_merger' , np.int32 ),
                             ( 'mass_of_main'   , np.int32 ),
                             ( 'z_merging'      , np.float32 ),
                             ( 'z_peak'         , np.float32 ),
                             ( 'z_appear'       , np.float32 ) ]

        elif record_length==56:

            self.cat_dtype=[ ( 'name'           , np.uint64 ),
                             ( 'nickname'       , np.int32 ),
                             ( 'link'           , np.int32 ),
                             ( 'merged_with'    , np.int32 ),
                             ( 'mass_at_merger' , np.int32 ),
                             ( 'mass_of_main'   , np.int32 ),
                             ( 'z_merging'      , np.float64 ),
                             ( 'z_peak'         , np.float64 ),
                             ( 'z_appear'       , np.float64 ) ]

            stored_dtype  =[ ( 'fort'           , np.int32 ),
                             ( 'name'           , np.uint64 ),
                             ( 'nickname'       , np.int32 ),
                             ( 'link'           , np.int32 ),
                             ( 'merged_with'    , np.int32 ),
                             ( 'mass_at_merger' , np.int32 ),
                             ( 'mass_of_main'   , np.int32 ),
                             ( 'pad'            , np.int32 ),
                             ( 'z_merging'      , np.float64 ),
                             ( 'z_peak'         , np.float64 ),
                             ( 'z_appear'       , np.float64 ),
                             ( 'trof'           , np.int32 ) ]

        else:
            print("sorry, I do not recognize this record length")
            return None

        # decides what files to read
        if Nfiles>1:
            if first_file is None:
                first_file=0
            elif first_file<0:
                first_file=0
            elif first_file>Nfiles:
                first_file=Nfiles
            if last_file is None:
                last_file=Nfiles
            else:
                last_file += 1      
                if last_file<first_file:
                    last_file=first_file+1
                elif last_file>Nfiles:
                    last_file=Nfiles
            # this is to be used to define a pythonic range
            if not silent:
                print("I will read files in the python range from {} to {}".format(first_file,last_file))
        else:
            first_file=0
            last_file=Nfiles

        # prepares to read the file(s)
        self.data=None
        self.Nbranches=None
        TTotal=np.int64(0)
        BTotal=np.int64(0)

        for myfile in range(first_file,last_file):

            if Nfiles==1:
                myfname = filename
            else:
                myfname = filename+'.{}'.format(myfile)

            if VERBOSE:
                print("reading file {}".format(myfname))

            # reads the file as a binary object
            with open(myfname,'rb') as f:
                bindata=f.read()
                f.close()
            FileLength = len(bindata)
            Tthisfile, Bthisfile = struct.Struct('ii').unpack(bindata[16:24])
            TTotal += Tthisfile
            BTotal += Bthisfile

            if newRun:

                if Light:
                    
                    # light format is straightforward to read
                    if VERBOSE:
                        print("reading in data")

                    if self.data is None:
                        self.data = np.copy(np.frombuffer(bindata[28:], dtype=self.cat_dtype))
                    else:
                        thiscat = np.frombuffer(bindata[28:], dtype=self.cat_dtype)
                        self.data.resize(self.data.shape[0]+len(thiscat))
                        self.data[-len(thiscat):]=np.copy(thiscat)
                        del thiscat

                else:

                    # reading standard V5 format
                    cleanForm='28x '
                    pos=28
                    Nblocks=0
                    while pos < FileLength:
                        # reads the number of branches
                        Nthisblock = struct.Struct('iii').unpack(bindata[pos:pos+12])[1]
                        pos+=16
                        Nbranches = np.frombuffer(bindata[pos:pos+Nthisblock*4],
                                                  dtype=np.int32,count=Nthisblock)
                        Bthisblock=Nbranches.sum()
                        pos+=Nthisblock*4 + 8 + Bthisblock*record_length + 4
                        cleanForm+='{}x {}s 4x'.format(12+Nthisblock*4+12,Bthisblock*record_length)
                        Nblocks+=1

                        if self.Nbranches is None:
                            self.Nbranches = np.copy(Nbranches)
                        else:
                            self.Nbranches.resize(self.Nbranches.shape[0]+Nthisblock)
                            self.Nbranches[-Nthisblock:]=np.copy(Nbranches)
                        del Nbranches

                    if VERBOSE:
                        print(f"data are written in {Nblocks} blocks")

                    # removes all unwanted information from the binary structure
                    try:
                        cleaned = b''.join(struct.unpack(cleanForm, bindata))
                        if VERBOSE:
                            print("cleaning of binary object done")
                    except:
                        print("ERROR: I do not recognise the data structure!")
                        return None
                    del bindata

                    thiscat = np.frombuffer(cleaned, dtype=self.cat_dtype)
                    if self.data is None:
                        self.data = np.copy(thiscat)
                    else:
                        self.data.resize(self.data.shape[0]+len(thiscat))
                        self.data[-len(thiscat):]=np.copy(thiscat)
                    del thiscat

            else:

                # first reads the total number of branches
                pos=12
                Bthisfile=0
                Tthisfile=0
                while pos < FileLength:
                    Nt, Nb = struct.Struct('iiii').unpack(bindata[pos:pos+16])[1:3]
                    pos+=16 + Nt*16 + Nb*(record_length+8)
                    Tthisfile+=Nt
                    Bthisfile+=Nb

                if self.data is None:
                    cNb=0
                    cTr=0
                    self.data      = np.empty(Bthisfile, dtype=self.cat_dtype)
                    self.Nbranches = np.empty(Tthisfile, dtype=np.int32)
                else:
                    cNb=self.Nbranches.shape[0]
                    cTr=self.data.shape[0]
                    self.data.resize(self.data.shape[0]+Bthisfile)
                    self.Nbranches.resize(self.Nbranches.shape[0]+Tthisfile)

                pos=12
                while pos < FileLength:
                    Nt = struct.Struct('iiii').unpack(bindata[pos:pos+16])[1]
                    pos+=16
                    for t in range(Nt):
                        Nb = struct.Struct('iiii').unpack(bindata[pos:pos+16])[2]
                        pos+=16
                        self.Nbranches[cNb]=Nb
                        cNb+=1
                        thiscat=np.frombuffer(bindata[pos:pos+Nb*(record_length+8)],dtype=stored_dtype,count=Nb)
                        for name in self.data.dtype.names:
                            self.data[name][cTr:cTr+Nb]=thiscat[name]
                        del thiscat
                        pos+=Nb*(record_length+8)
                        cTr+=Nb


            if newRun:
                if Light:
                    file_length = 28 + Bthisfile*record_length
                else:
                    file_length = 28 + Nblocks*28 + Tthisfile*4 + Bthisfile*record_length
            else:
                file_length = 12 + Nslices*16 + Tthisfile*16 + Bthisfile*(record_length+8)
            if FileLength != file_length:
                print("Error: inconsistent length for the file {}, I expected {} and found {}".format(myfname,FileLength,file_length))
                return None
            elif VERBOSE:
                print("File size of {} is as expected".format(myfname))


        if newRun and not Light:
            # cumulate pointers
            self.pointers = np.cumsum(np.insert(self.Nbranches,0,0))[:TTotal]
        else:
            if not silent:
                print("Building pointers...")
            po = np.empty(TTotal+1, dtype=np.int32)
            self.Nbranches = np.empty(TTotal, dtype=np.int32)

            po[0]=0
            for i in range(TTotal):
                self.Nbranches[i]=self.data[po[i]]['nickname']
                po[i+1]=po[i]+self.Nbranches[i]

            self.pointers = po[0:TTotal]
            del po
            
        self.Ntrees=TTotal
        self.Nbranches_tot=BTotal

        if not silent:
            print("Reading catalog done")
