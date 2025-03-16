import numpy as np
from astropy.io import fits
import ReadPinocchio5 as rp
import sys

'''
This script converts Pinocchio binary output to fits

Usage: python3 Pinocchio2fits.py [pinocchio parameter file]
'''

VERSION='v5.0'
CATALOGS=True
PLC=True
HISTORIES=True

params={
    'RunFlag'                        :['     ','string  ','name of the run',''],
    'OutputList'                     :['     ','string  ','filename with required output redshifts',''],
    'BoxSize'                        :['     ','int     ','physical size of the box [Mpc or Mpc/h]',''],
    'BoxInH100'                      :['False','bool    ','specify that the box is in Mpc/h',''],
    'GridSize'                       :['     ','int     ','number of grid points per side',''],
    'RandomSeed'                     :['     ','int     ','random seed for initial conditions',''],
    'FixedIC'                        :['False','bool    ','modulus in ICs is fixed to the average',''],
    'PairedIC'                       :['False','bool    ','phase in ICs is shifted by PI',''],
    'Omega0'                         :['     ','float   ','Omega_0 (total matter)',''],
    'OmegaLambda'                    :['     ','float   ','Omega_Lambda',''],
    'OmegaBaryon'                    :['     ','float   ','Omega_b (baryonic matter)',''],
    'Hubble100'                      :['     ','float   ','little h',''],
    'Sigma8'                         :['     ','float   ','sigma8; if 0, it is computed from P(k)',''],
    'PrimordialIndex'                :['     ','float   ','spectral index n_s',''],
    'DEw0'                           :['     ','float   ','w0 of dark energy EoS',''],
    'DEwa'                           :['     ','float   ','wa of dark energy EoS',''],
    'TabulatedEoSfile'               :['     ','string  ','dark energy EoS tabulated in file',''],
    'FileWithInputSpectrum'          :['     ','string  ','P(k) tabulated in file',''],
    'InputSpectrum_UnitLength_in_cm' :['     ','float   ','units of tabulated P(k) [cm]',''],
    'WDM_PartMass_in_kev'            :['     ','float   ','WDM cut [keV]',''],
    'BoundaryLayerFactor'            :['     ','float   ','width of boundary layer for fragmentation',''],
    'MaxMem'                         :['     ','int     ','max available memory to an MPI task [Mbyte]',''],
    'MaxMemPerParticle'              :['     ','int     ','max available memory [bytes per particle]',''],
    'PredPeakFactor'                 :['     ','float   ','guess for the number of peaks in subvolume',''],
    'CatalogInAscii'                 :['False','bool    ','catalogs are written in ascii, not binary',''],
    'OutputInH100'                   :['False','bool    ','units are in H=100 instead of true H value',''],
    'NumFiles'                       :['     ','int     ','number of files for each catalog',''],
    'MinHaloMass'                    :['     ','int     ','smallest halo in output [particles]',''],
    'AnalyticMassFunction'           :['     ','int     ','analytic mass function given in *.mf.out',''],
    'WriteTimelessSnapshot'          :['False','bool    ','write timeless snapshot',''],
    'DoNotWriteCatalogs'             :['False','bool    ','skip writing catalogs (including PLC)',''],
    'DoNotWriteHistories'            :['False','bool    ','skip writing merger histories',''],
    'StartingzForPLC'                :['     ','float   ','starting (highest) redshift for PLC',''],
    'LastzForPLC'                    :['     ','float   ','final (lowest) redshift for PLC',''],
    'PLCAperture'                    :['     ','float   ','cone aperture for PLC [deg]',''],
    'PLCProvideConeData'             :['False','bool    ','read vertex and direction of cone',''],
    'PLCCenter'                      :['     ','3 floats','cone vertex [same coordinates as BoxSize]',''],
    'PLCAxis'                        :['     ','3 floats','un-normalized direction of cone axis',''],
    'CTtableFile'                    :['     ','string  ','filename with collapse time table to read in',''],
    'CAMBMatterFileTag'              :['     ','string  ','label for CAMB matter power spectrum files',''],
    'CAMBTransferFileTag'            :['     ','string  ','label for CAMB transfer function files',''],
    'CAMBRunName'                    :['     ','string  ','name of CAMB run',''],
    'CAMBRedshiftsFile'              :['     ','string  ','list of redshifts of CAMB power spectra','']}


def read_param_file(param_fname):


    try:
        f=open(param_fname)
    except:
        print(f'Error: parameter file {param_fname} not found')
        return 1

    print(f'Parameters in file: {param_fname}')

    for line in f:

        if line[0]=='#' or line[0]=='%':
            continue
        
        keys=line.split()

        if len(keys)==0 or keys[0]=='%':
            continue

        if keys[0] not in params:
            print(f"WARNING: unrecognized parameter, {keys[0]}")

        else:

            if params[keys[0]][0]=='False':
                params[keys[0]][0]='True'
            elif keys[0]=='PLCCenter' or keys[0]=='PLCAxis':
                params[keys[0]][0]=keys[1]+' '+keys[2]+' '+keys[3]
            else:
                params[keys[0]][0]=keys[1]

    f.close()    

    return 0


def create_fits_file(output_file, column_data, column_names):

    # Create a FITS binary table extension
    columns = [fits.Column(name=col, format='D', array=np.array(data)) for col, data in zip(column_names, column_data)]
    table = fits.BinTableHDU.from_columns(columns)

    # Save the FITS file
    table.writeto(output_file, overwrite=True)


if len(sys.argv)<2:

    print('Usage: python3 Pinocchio2fits.py [pinocchio parameter file]')
    sys.exit(0)

param_fname=sys.argv[1]

if read_param_file(param_fname):
    sys.exit(0)

runflag=params['RunFlag'][0]
print(f'RunFlag: {runflag}')

output_fname=params['OutputList'][0]

print(f'Output list in file: {output_fname}')

try:
    outputs=np.loadtxt(output_fname,unpack=True)
except:
    print(f'Error: outputs file {output_fname} not found')
    sys.exit(0)

print(f'Output list: {outputs}')

def write_fits(fname,cat_type,cat):

    print(f'Writing fits catalog {fname}')

    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header.append(('CODE',f'Pinocchio {VERSION}','https://github.com/pigimonaco/Pinocchio'))
    hdul = fits.HDUList([primary_hdu])
    hdu = fits.BinTableHDU(cat.data)
    hdu.name=cat_type
    counter=1
    hdu.header.append(('NHALOS',cat.Nhalos,'Number of halos in catalog'))
    for key in params:
        hdu.header.append((f'PAR{counter}',key,params[key][2]))
        hdu.header.append((f'VAL{counter}',params[key][0],params[key][3]))
        counter+=1
    hdul.append(hdu)
    hdul.writeto(fname, overwrite=True)

    print(f'Writing done')


def write_fits_histories(fname,cat):

    print(f'Writing fits catalog {fname}')

    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header.append(('CODE',f'Pinocchio {VERSION}','https://github.com/pigimonaco/Pinocchio'))
    hdul = fits.HDUList([primary_hdu])
    data = fits.BinTableHDU(cat.data)
    myarr=np.empty((cat.Ntrees,),dtype=[('Nbranches', np.uint32), ('pointers', np.uint32)])
    myarr['Nbranches']=cat.Nbranches
    myarr['pointers']=cat.pointers
    pointers = fits.BinTableHDU(myarr)
    data.name='histories'
    pointers.name='pointers'
    counter=1
    data.header.append(('NTREES',cat.Ntrees,'number of trees'))
    data.header.append(('NBRANCH',cat.Nbranches_tot,'number of branches'))
    for key in params:
        data.header.append((f'PAR{counter}',key,params[key][2]))
        data.header.append((f'VAL{counter}',params[key][0],params[key][3]))
        counter+=1
    pointers.header.append(('NTREES',cat.Ntrees,'number of trees'))
    pointers.header.append(('NBRANCH',cat.Nbranches_tot,'number of branches'))
    hdul.append(data)
    hdul.append(pointers)
    hdul.writeto(fname, overwrite=True)

    print(f'Writing done')

#
# catalogs
#

if CATALOGS:
    for z in outputs:

        pin_fname=f'pinocchio.{z:6.4f}.{runflag}.catalog.out'
        print(f'Reading catalog {pin_fname}')

        cat=rp.catalog(pin_fname)

        print(f'The catalog contains the following fields: {cat.data.dtype.names}')

        fits_fname=pin_fname[:-3]+'fits'

        write_fits(fits_fname,'catalog',cat)


if PLC:
    pin_fname=f'pinocchio.{runflag}.plc.out'
    print(f'Reading catalog {pin_fname}')

    cat=rp.plc(pin_fname)

    print(f'The catalog contains the following fields: {cat.data.dtype.names}')

    fits_fname=pin_fname[:-3]+'fits'

    write_fits(fits_fname,'plc',cat)


if HISTORIES:
    pin_fname=f'pinocchio.{runflag}.histories.out'
    print(f'Reading catalog {pin_fname}')

    cat=rp.histories(pin_fname)

    print(f'The catalog contains the following fields: {cat.data.dtype.names}')

    fits_fname=pin_fname[:-3]+'fits'

    write_fits_histories(fits_fname,cat)




