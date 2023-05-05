import numpy as np
from PIL import Image
from astropy.io import fits

from .config import RB
from .config import __version__


def chsuf(i_string, new_suffix, strip_prefix=False, trim_prefix=False,
          basename=False):
    """
    Change string suffixes. This is useful for generating file names.

    Inputs:
        i_string   - The input string containing a suffix of the form '.suffix'
        new_suffix - The new suffix that is desired
      strip_prefix - Optionally strip off the prefix. This is useful when
                     you do not want results written into the data directory.
       trim_prefix - Same as strip_prefix
          basename - Same as strip_prefix
        
    Returns:
        A string having the specified suffix
    """
    # Strip off the prefix if requested
    if strip_prefix or trim_prefix or basename:
        i_string = os.path.basename(i_string)

    # Get everything before the prefix
    here = i_string.rfind('.')
    result = i_string[:here]
    result += new_suffix

    # Done
    return(result)



def make_lowpass_filter(f_half_power, w_cutoff, n, d=1.0):
    """
    Make a lowpass Fourier filter
    
    Parameters: f_half_power, number
                  Half power frequency
                w_cutoff, number
                  Width of cosine cutoff. The response transitions
                  from 1x to 0x over this range of frequencies
                n, int
                  Number of samples in the timeseries
                d, float
                  Sample spacing (inverse of the sampling rate). Defaults to 1.
    """
    
    # Make frequencies vector
    freq = np.fft.rfftfreq(n, d=d)
    
    # Build a cosine wave that is approriately shifted
    cos = (1 + np.cos(np.pi*(freq - f_half_power)/w_cutoff+np.pi/2))/2
    
    # Construct low-pass filter with cosine rolloff
    filt = np.where(freq<=f_half_power-w_cutoff/2,1,cos)
    filt = np.where(freq<=f_half_power+w_cutoff/2,filt,0)
    
    # Done!
    return(filt)



def png2fits(input_filename, output_filename, overwrite=False):
    """
    Convert PNG format background masks to FITS
    
    Parameters: input_filename, string
                  Input file name
                output_filename, string
                  Output file name
                overwrite, bool
                  Overwrite the output file
    Notes:
    1) The input files should be NIRSpec images in the DMS format
       with spectral dispersion running horizontally. This is what
       the STScI pipeline provides.
    2) Pixel values equal to pure red, = ff0000x, are taken to be background.
       Background pixels should be known to be either blanked off within the
       instrument, or unilluminated insofar as practical
    3) I typically create the input file by using SAOImage DS9 to export
       a 2048x2048 pixel grayscale image (one per NIRSpec channel). I then
       use the GNU Image Manipulation Program (GIMP) to color selected
       pixels perfectly red. Something automated would be better...
    """
    # Get data
    B = np.array(Image.open(input_filename))
    
    # Set pure red = True and everything else false
    B = (B[:,:,0] == 255) & (B[:,:,1] == 0) & (B[:,:,2] == 0)

    # Always ignore the reference pixels. The reference pixel border is always
    # 4 pixels wide in NIRSpec
    B[:RB,:] = False
    B[-RB:,:] = False
    B[:,:RB] = False
    B[:,-RB:] = False

    # Correct orientation to match STScI DMS
    B = B[::-1,:]
    
    # Convert to a data type available in FITS
    B = np.array(B, dtype=np.uint8)
    
    # Write results
    hdul = fits.PrimaryHDU(B)
    hdul.header['comment'] = 'Pixels = 1 are background'
    hdul.header['history'] = 'Created by png2fits Rev. ' + __version__
    hdul.writeto(output_filename, overwrite=overwrite)