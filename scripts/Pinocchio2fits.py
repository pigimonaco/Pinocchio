import numpy as np
from astropy.io import fits

# Prototype for converting Pinocchio PLC output to a fits catalog 
# This should be extended to final catalogs and merger histories
# Should we keep the .npz format?
# Should we integrated also the ReadPinocchio5 in case (as normally done) the final output is in binary format?

def read_columns_from_text(input_file, column_indices, skip_lines):
    
    # Create empty arrays for each column
    column_data = [[] for _ in column_indices]

    with open(input_file, 'r') as file:
        for _ in range(skip_lines):
            next(file)  # Skip the specified number of lines
            
        for line in file:
            parts = line.split()  # Split the line by whitespace
            
            # Check if the line has enough columns
            if max(column_indices) < len(parts):
                for i, col_idx in enumerate(column_indices):
                    column_data[i].append(float(parts[col_idx]))  # Store data

    return column_data


def read_columns_from_npz(input_file, column_names):

    # Load the NumPy .npz file
    cat = np.load(input_file)

    # Create empty arrays for each column
    column_data = [cat[column_name] for column_name in column_names]

    return column_data

def create_fits_file(output_file, column_data, column_names):

    # Create a FITS binary table extension
    columns = [fits.Column(name=col, format='D', array=np.array(data)) for col, data in zip(column_names, column_data)]
    table = fits.BinTableHDU.from_columns(columns)

    # Save the FITS file
    table.writeto(output_file, overwrite=True)

def read_and_save_to_fits(input_file, output_fits_file, column_indices, column_names, skip_lines=0):
    
    # Determine the file type based on extension
    if input_file.lower().endswith('.txt'):

        # Read specific columns from the text file
        data = read_columns_from_text(input_file, column_indices, skip_lines)

    elif input_file.lower().endswith('.npz'):

        # Read specific columns from the npz file
        column_names = ['mass', 'redshift', 'theta', 'phi']
        data = read_columns_from_npz(input_file, column_names)
    else:
        raise ValueError("Unsupported file type. Only .txt and .npz files are supported.")

    # Create and save FITS file
    create_fits_file(output_fits_file, data, column_names)


###########################################################################################################################

# Specify your input and output file paths

input_text_file = "converted_plc_mock0001.npz"
output_fits_file = f'{input_text_file}.fits'

skip_lines = 11
column_indices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
column_names = ['group_ID', 'True_Z', 'X_POS', 'Y_POS', 'Z_POS', 'VEL_X', 'VEL_Y', 'VEL_Z', 'MASS', 'THETA', 'PHI', 'PEC_VEL', 'Observed_Z']

# Skip lines, column indices, and column names are not used for .npz files
read_and_save_to_fits(input_text_file, output_fits_file, None, None)
