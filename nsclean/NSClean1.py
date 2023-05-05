import os

# If you want to use a GPU, set the NSCLEAN_USE_GPU
# environment variable == 'YES'
if os.getenv('NSCLEAN_USE_CUPY') == 'YES':
    import cupy as cp
else:
    import numpy as cp

class NSClean1:
    """
    RECOMMENATION - If you are able to specify which pixels should be treated as background for
    your observations, consider using NSClean, the enhanced second-generation
    NIRSpec cleaner. It does a much better job of removing correlated read noise
    than this tool.

    NSClean1 is a first generation NIRSpec cleaner. It uses default areas of background
    pixels to fit and subtract a highly smoothed background model. Compared to the
    NIRSpec Team's "rolling median" technique, it marginally more uniform and more
    tunable.
    """
    
    # Class variables. These are the same for all instances.
    ny     = 2048     # Number of lines in detector space
    nx     = 2048     # Number of columns in detector space
    n      = ny*nx    # Total number of pixels
    stride = 128      # Fit this many Fourier vectors at a time. This is
                      #   sized to allow the software to run on mainstream GPUs.
                      #   For development, I used a GPU with 8 GB of RAM.
    sigrej = 3        # Standard deviation threshold for flagging statistical
                      #   outliers.
    sigma_okern = 64  # Statistical outliers are filled using a Gauss-weighted
                      #   average of neighbors. This is the standard deviation
                      #   of that kernel. Bigger is more robust vs large groups of outliers.
    
    def __init__(self, opmode, detector, buffer_sigma=2.5):
        """
        RECOMMENDATION: If you are able to specify which pixels should be gtreated
        as background, use NSClean. It produces much cleaner results than this
        first generation tool.

            Parameters: opmode, string
                          Operating mode selected from {'MOS','IFU'}
                        detector, string
                          A string selected from {'NRS1','NRS2'}
                        buffer_sigma, float
                          Standard deviation of the buffer's Gaussian smooth. This is an
                          optional 1-dimensional smooth in the spectral dispersion direction
                          only. It is done after fitting the background.
        
        """
        # Pick off kwargs
        if (detector != 'NRS1') & (detector != 'NRS2'):
            print('ERROR: Invalid detector. Detector must be NRS1 or NRS2')
            return
        self.detector = detector
        self.buffer_sigma = buffer_sigma
        
        # Set mode
        if opmode=='MOS':
            
            # Define number of Fourier frequencies
            self.nvec = 1024+128 # Use only blanked-off rows behind the MSA cruciform
            self.kw   = 256      # The apodizer rolls the resoponse off over this many frequencies
            
            # Make a mask of the pixels to use for background modeling.
            # Use only pixels that are known to be blanked off by the MSA cruiciform.
            # Exclude reference pixels because these have already been used by IRS^2.
            # This is in detector coordinates with the IRS^2 "zipper" running along
            # the bottom as displayed in ds9.
            self.bgpx = cp.zeros((self.ny,self.nx), dtype=bool)
            self.bgpx[4:-4,830:930] = True    # From Stephan Birkmann
            self.bgpx[4:-4,1110:1200] = True  # From Stephan Birkmann
             
        elif opmode=='IFU':
            
            # Define number of Fourier frequencies to fit
            self.nvec = 1536+128 # Assumes blanked off strips near
                                 #   bottom, middle, and top in DMS format.
                                 #   DMS format has the dispersion running
                                 #   horizontally as displayed in SAOImage ds9.
                                 #   2048 frequencies are enough to resolve the
                                 #   change in the "picture frame" noise pattern
                                 #   between the middle of an array and the edges.
            self.kw   = 256      # "Kill width". Apodizer rolls to zero over this many frequencies.
                        
            # Make a mask of the dark pixels to fit to model the background.
            # The ones near the middle are from Stephan Birkmann. It looks like
            # the first and last few rows in DMS format also stay pretty dark. These
            # definitions are in detector coordinates with the IRS^2 zipper running
            # along the bottom.
            self.bgpx = cp.zeros((self.ny,self.nx), dtype=bool)
            self.bgpx[4:-4,830:916] = True    # From Stephan, but reduced to 86 cols.
            self.bgpx[4:-4,1110:1196] = True  # From Stephan, but reduced to 86 cols.
            self.bgpx[4:-4,4:90] = True       # Near bottom, 86 cols.
            self.bgpx[4:-4,-90:-4] = True     # Near top, 86 cols.
            
        else:
            print('NSClean Error: Mode must be either MOS or IFU')
            return
            
        # Make a vector of the pixel indices to fit. This must be a row vector
        # to broadcast correctly.
        self.m = cp.arange(self.n, dtype=cp.float32)[self.bgpx.flatten()].reshape((1,-1))
            
        # Make a vector of the Fourier frequencies to fit. This must
        # be a column vector to broadcast correctly.
        self.k = cp.arange(self.nvec, dtype=cp.float32).reshape((-1,1))
                        
        # Build the apodizer. The data are not fully sampled. To not push things, this
        # rolls the response to near zero for the highest fitted Fourier vector.
        _k = 2*cp.pi/(2*self.kw)  # Wave number for a cosine shaped roll off
        _x = cp.arange(self.kw)   # X-values for generating the cosine rolloff
        self.apodizer = cp.ones(self.nvec, dtype=cp.float32) # Generate apodizer
        self.apodizer[-self.kw:] = 0.5*(cp.cos(_k*_x)+1) # Roll off at high frequency
        
        # Build a 1-dimensional Gaussian kernel for "buffing". Buffing is in the
        # dispersion direction only. In detector coordinates, this is axis zero. Even though
        # the kernel is 1-dimensional, we must still use a 2-dimensional array to 
        # represent it. I tried broadcasting a vector, but that made a kernel 2048 
        # columns wide (in detector space).
        _y = cp.arange(self.ny) # Row indices
        _mu = self.ny//2+1 # Center of kernel
        _sigma = self.buffer_sigma # Standard deviation of kernel
        _gkern = cp.exp(-((_y-_mu)/_sigma)**2/2)/_sigma/cp.sqrt(2*cp.pi) # Centered kernel as a vector
        gkern = cp.zeros((self.ny,self.nx), dtype=cp.float32) # 2D kernel template
        gkern[:,_mu] = _gkern # Copy in the kernel. Normalization is already correct.
        gkern = cp.fft.ifftshift(gkern) # Shift for Numpy/CuPy
        self.fgkern = cp.array(cp.fft.rfft2(gkern), dtype=cp.complex64) # FFT for fast convolution
        
        # Build the 1-dimensional kernel for filling statistical outliers
        _x = cp.arange(self.m.shape[1])
        _mu = self.m.shape[1]//2+1
        _sigma = self.sigma_okern
        okern = cp.exp(-((_x-_mu)/_sigma)**2/2)/_sigma/cp.sqrt(2*cp.pi)
        okern = cp.fft.ifftshift(okern) # Orient for Numpy/CuPy
        self.fokern = cp.array(cp.fft.rfft(okern), dtype=cp.complex64)
        
    def fit(self, D):
        """
            fit(D::CuArray{Float32,2})
        
        Fit a background model to the supplied frame of data.
        
            Parameters: D::CuArray{Float32,2}
                          A NIRSpec image. The image must be in the detector-space
                          orientation with the IRS^2 zipper running along the bottom as
                          displayed in SAOImage DS9.
               Returns: B::CuArray{Float32,2}
                          The fitted background model.
        """

        # ***** Fill statistical outliers with a locally smoothed values
        
        # Extract just the background pixels and robustly estimate the mean and
        # standard deviation using median absolute deviation.
        d = D[self.bgpx]
        mu = cp.median(d)
        sigma = 1.4826*cp.median(cp.abs(d-mu))
        
        # Set statistical outliers and nearest neighbors identically =0.
        # This will occasionally generate false positives, but it is computationally efficient.
        d = cp.where(cp.logical_and(mu-self.sigrej*sigma<=d,d<=mu+self.sigrej*sigma), d, 0.)
        
        # Fill in the zeros with locally smoothed values
        g = cp.array(cp.where(d==0, 0., 1.), dtype=cp.float32) # Good sample mask
        f = cp.fft.irfft(cp.fft.rfft(d) * self.fokern, len(d)) / \
                    cp.fft.irfft(cp.fft.rfft(g) * self.fokern, len(g)) # Fill values
        d = cp.array(cp.where(d==0., f, d), dtype=cp.float32) # Backfill
        # ***** End filling statistical outliers *****
        
        # Project out Fourier vectors
        rfft = cp.zeros((self.n//2+1), dtype=cp.complex64)
        for i in cp.arange(0,self.nvec,self.stride):
            # Build one stride of the Fourier matrix
            CuF = cp.float32(self.n)*cp.exp(cp.complex64(-2*cp.pi*1J)*
                      self.m*self.k[i:i+self.stride,:]/cp.float32(self.n))/\
                            cp.float32(len(self.m.flatten()))
            # Compute Fourier vectors for this stride
            rfft[i:i+self.stride] = cp.matmul(CuF, d)
            
        # Apodize
        rfft[:self.nvec] *= self.apodizer
        
        # Invert the FFT to get the background model. This is in
        # detector coordiantes
        B = cp.fft.irfft(rfft, self.n).reshape((self.ny,self.nx))
                                               
        # Done!
        return(B)
    

    def clean(self, D, buff=True):
        """
            clean(D)
            
        "Clean" NIRspec images by fitting and subtracting the instrumental background. 
        This is intended to improve the residual vertical banding that is sometimes seen
        in NRS2. Because the NRS2 banding is not picked up by the reference pixels,
        normal IRS^2 processing does not completely remove it.
        
        Unlike IRS^2, this is an ad-hoc correction. We model the background by fitting it
        in Fourier space. There is an option to "improve" the result by "buffing" in the
        spectral dispersion direction. "Buffing" is just smoothing using a 1-dimensional
        Gaussian.
        
            Parameters: D::array_like
                          The input data. This should be the normal end result of
                          Stage 1 processing.
               Returns: D::array_like
                          The data, but with less striping and the background subtracted.
                        buff, bool
                          "Buff" the fitted spectrum by applying a slight Gaussian blur
                          in the spectral dispersion direction.
        """
        
        # Transform the data to detector space with the IRS^2
        # zipper running along the bottom.
        if self.detector=='NRS2':
            D = D.transpose()[::-1] # Go to detector space
        else:
            D = D.transpose()       # No flip required for NRS1
            
        # Fit, optionally buff, and subtract the background model
        B = self.fit(D) # Background model
        if buff is True:
            B = cp.fft.irfft2(cp.fft.rfft2(B) * self.fgkern, s=B.shape) # Buff
        D -= B # Subtract background model
        
        # Transform back to DMS space
        if self.detector=='NRS2':
            D = D[::-1].transpose() # Go back to DMS space
        else:
            D = D.transpose()       # No flip required for NRS1
            
        # Done
        return(D)
