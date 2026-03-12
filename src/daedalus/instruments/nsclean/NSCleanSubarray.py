import os
if os.getenv('NSCLEAN_USE_CUPY') == 'YES':
    import cupy as cp
    import numpy as np
else:
    import numpy as np
    import numpy as cp

def mad(d, median=True):
    """
        Median absolute deviation
        
    Computes the median and the median absolute deviation (MAD). For normally
    distributed data, multiply the MAD by 1.4826 to approximate standard deviation. 
    If median=True (the default), this returns a tuple containing (median, MAD).
    Otherwise, only the MAD is returned.
    
        Parameters: d, ndarray
                      The input data
                    median, bool
                      Return both the median and MAD
    """
    d = d[cp.isfinite(d)] # Exclude NaNs
    m = cp.median(d)
    mad = cp.median(cp.abs(d-m))
    if median is True:
        return(m, mad)
    else:
        return(mad)
    
    

class NSCleanSubarray:
    """
    NSCleanSubarray is the base class for removing residual correlated
    read noise from generic JWST near-IR Subarray images.  It is
    intended for use on Level 1 pipeline products, i.e. slope images.
    """
    
    # Class variables. These are the same for all instances.
    nloh = cp.int32(12)       # New line overhead in pixels
    tpix = cp.float32(10.e-6) # Pixel dwell time in seconds
    sigrej = cp.float32(4.0)  # Standard deviation threshold for flagging
                              #   statistical outliers.
        
    def __init__(self, D, M, fc=(1061, 1211, 49943, 49957),
                exclude_outliers=True, weights_kernel_sigma=None):
        """
        Background modeling and subtraction for generic JWST near-IR
        subarrays.

        Parameters: D, ndarray::float32
                      Two-dimensional input data.
                    M, ndarray::bool
                      Two dimensional background pixels mask. Pixels =True are modeled as background.
                      Pixels=False are excluded from the background model.
                    fc, tuple
                      Apodizing filter definition. These parameters are tunable. They
                      happen to work well for NIRSpec BOTS exposures.
                        1) Unity gain for f < fc[0]
                        2) Cosine roll-off from fc[0] to fc[1]
                        3) Zero gain from fc[1] to fc[2]
                        4) Cosine roll-on from fc[2] to fc[3]
                    exclude_outliers, bool
                      Exclude statistical outliers and their nearest neighbors
                      from the background pixels mask.
                    weights_kernel_sigma, float
                      Standard deviation of the 1-dimensional Gaussian kernel
                      that is used to approximate background sample density. This
                      is ad-hoc. See the NSClean journal article for more information. The
                      default for subarrays results in nearly equal weighting of all background
                      samples.
        Notes:
        1) NSCleanSubarray works in detector coordinates. Both the data and mask
           need to be transposed and flipped so that slow-scan runs from bottom
           to top as displayed in SAOImage DS9. The fast scan direction is 
           required to run from left to right.
        """
        # Definitions
        self.D = cp.array(D, dtype=cp.float32)
        self.M = cp.array(M, dtype=cp.bool_)
        self.ny = cp.int32(D.shape[0]) # Number of pixels in slow scan direction
        self.nx = cp.int32(D.shape[1]) # Number of pixels in fast scan direction
        self.fc = cp.array(fc, dtype=cp.float32)
        self.n = cp.int32(self.ny * (self.nx + self.nloh)) # Number of ticks in clocking pattern
        self.rfftfreq = cp.array(cp.fft.rfftfreq(self.n, self.tpix),
                                    dtype=cp.float32)  # Fourier frequencies in clocking pattern
        
        # We will weight be the inverse of the local sample density in time. We compute the local
        # sample density by convolution using a Gaussian. Define the standard deviation of the
        # Gaussian here. This is ad-hoc.
        if weights_kernel_sigma == None:
            self.weights_kernel_sigma = 1/((fc[0]+fc[1])/2)/self.tpix/2/4
        else:
            self.weights_kernel_sigma = weights_kernel_sigma
        
        # The mask potentially contains NaNs. Exclude them.
        self.M[cp.isnan(self.D)] = False
        
        # The mask potentially contains statistical outliers.
        # Optionally exclude them.
        if exclude_outliers is True:
            m,s = mad(self.D[self.M]) # Compute median and median absolute deviation
            s *= 1.4826 # Convert MAD to std
            vmin = m - self.sigrej*s # Minimum value to keep
            vmax = m + self.sigrej*s # Maximum value to keep
            # NaNs were causing problems in the following line. Briefly change them
            # to inf.
            self.D[cp.isnan(self.D)] = cp.inf # Get rid of NaNs
            bdpx = cp.array(cp.where(cp.logical_or(self.D<vmin,self.D>vmax),
                                          1, 0), dtype=cp.float32) # Flag statistical outliers
            self.D[cp.isinf(self.D)] = cp.nan # Restore NaNs
            bdpx[cp.logical_not(self.M)] = 0 # We don't need to worry about non-background pixels
            # Also flag 4 nearest neighbors
            bdpx = bdpx +\
                        cp.roll(bdpx, (+1,0), axis=(0,1)) +\
                        cp.roll(bdpx, (-1,0), axis=(0,1)) +\
                        cp.roll(bdpx, (0,+1), axis=(0,1)) +\
                        cp.roll(bdpx, (0,-1), axis=(0,1))
            # bdpx now contains the pixels to exclude from the background pixels
            # mask. Exclude them.
            self.M[bdpx != 0] = False
            self.D -= m # STUB - Median subtract
        
        # Build the apodizing filter. This has unity gain at low frequency to capture 1/f. It 
        # also has unity gain at Nyquist to capture alternating column noise. At mid frequencies, the gain is
        # zero.
        # Unity gain at low frequencies
        self.apodizer = cp.zeros(len(self.rfftfreq), dtype=cp.float32)
        self.apodizer[self.rfftfreq<self.fc[0]] = 1.0
        # Cosine roll-off
        here = cp.logical_and(self.fc[0]<=self.rfftfreq, self.rfftfreq<self.fc[1])
        self.apodizer[here] = 0.5*(cp.cos((cp.pi/(self.fc[1]-self.fc[0])) *
                                          (self.rfftfreq[here] - self.fc[0]))+1)
        # Cosine roll-on
        here = cp.logical_and(self.fc[2]<=self.rfftfreq, self.rfftfreq<self.fc[3])
        self.apodizer[here] = 0.5*(cp.cos((cp.pi/(self.fc[3]-self.fc[2])) *
                                          (self.rfftfreq[here] - self.fc[2]) + cp.pi)+1)
        # Unity gain between f[3] and end
        self.apodizer[self.rfftfreq>=self.fc[3]] = 1.0
        

        
    def fit(self, return_fit=False, weight_fit=False):
        """
        Fit a background model to the data.
        
        Parameters: return_fit:bool
                      Return the Fourier transform.
                    weight_fit:bool
                      Use weighted least squares as described in the NSClean paper.
                      Turn off by default. For subarrays it is TBD if this is necessary.
        """
        
        # To build the incomplete Fourier matrix, we require the index of each
        # clock tick of each valid pixel in the background samples. For consistency with
        # numpy's notation, we call this 'm' and require it to be a column vector.
        _x = cp.arange(self.nx).reshape((1,-1))
        _y = cp.arange(self.ny).reshape((-1,1))
        m = (_y*(self.nx+self.nloh) + _x)[self.M].reshape((-1,1))
        
        # Define which Fourier vectors to fit. For consistency with numpy, call this k.
        k = cp.arange(len(self.rfftfreq))[self.apodizer>0.].reshape((1,-1))
        
        # Build the incomplete Fourier matrix
        B = cp.array(cp.exp(2*cp.pi*1J*m*k/self.n)/m.shape[0], dtype=cp.complex64)
        
        # Weighted NSClean fitting
        if weight_fit is True:
            
            # Build the weight matrix. Weight by the reciprocal of the local background
            # sample density in time. Roughly approximate the local density, P, using
            # the reciprocal of convolution by a Gaussian kernel.
            _x = cp.arange(self.n, dtype=cp.float32) # x-values for building a 1-D Gaussian
            _mu = cp.float32(self.n//2+1) # Center point of Gaussian
            _sigma = cp.float32(self.weights_kernel_sigma) # Standard deviation of Gaussian
            W = cp.exp(-(_x-_mu)**2/_sigma**2/2)/_sigma/cp.sqrt(2*cp.pi) # Build centered Gaussian
            FW = cp.fft.rfft(cp.fft.ifftshift(W)) # Forward FFT
            _M = cp.hstack((self.M, cp.zeros((self.ny,self.nloh), dtype=cp.bool_))).flatten() # Add new line overhead to mask
            with np.errstate(divide='ignore'):
                P = 1/cp.fft.irfft(cp.fft.rfft(cp.array(_M, dtype=cp.float32)) * FW, self.n) # Compute weights
            P = P[_M==True] # Keep only background samples
        
            # NSClean's weighting requires the Moore-Penrose invers of A = P*B.
            #     $A^+ = (A^H A)^{-1} A^H$
            A = P.reshape((-1,1)) * B # P is diagonal. Hadamard product is most RAM efficient
            AH = cp.conjugate(A.transpose()) # Hermitian transpose of A
            pinv_PB = cp.matmul(cp.linalg.inv(cp.matmul(AH, A)), AH)
            
        else:
            # Unweighted fit
            pinvB = cp.linalg.pinv(B)
        
        # Solve for the (approximate) Fourier transform of the background samples.
        rfft = cp.zeros(len(self.rfftfreq), dtype=cp.complex64)
        if weight_fit is True:
            rfft[self.apodizer>0.] = cp.matmul(pinv_PB * P.reshape((1,-1)), self.D[self.M])
        else:
            rfft[self.apodizer>0.] = cp.matmul(pinvB, self.D[self.M])
        
        
        # Numpy requires that the forward transform multiply
        # the data by n. Correct normalization.
        rfft *= self.n / m.shape[0]
        
        # Invert the apodized Fourier transform to build the background model for this integration
        self.model = cp.fft.irfft(rfft*self.apodizer, self.n).reshape((self.ny,-1))[:,:self.nx]
        
        # Done
        if return_fit is True:
            return(rfft)
        
        
        
    def clean(self, weight_fit=True):
        """
        Clean the data
        
        Parameters: weight_fit:bool
                      Use weighted least squares as described in the NSClean paper.
                      Otherwise, it is a simple unweighted fit.

        """ 
        self.fit(weight_fit=weight_fit)           # Fit the background model
        self.D -= self.model # Overwrite data with cleaned data
        return(self.D)
