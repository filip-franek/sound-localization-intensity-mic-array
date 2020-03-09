This is the time-domain blind audio source separation method called T-ABCD by Zbynìk Koldovský,
Petr Tichavský and Jiøí Málek.

release: July 20, 2011


Instalation and execution
-------------------------

unzip the package into a directory, run Matlab, load microphone recordings into rows of a matrix x, 
and type (in Matlab)

>> tddeconv(x)


Usage of GUI
------------

Click signal to listen to either in the "Mixed signals" plot or in the "Separated signals" plot. Note 
that the correct sampling frequency should be selected.

"Mixed signals" plot > Click left and right button to select, respectively, the offset of 
the block of data for ICA and the end of the block (the green are in the plot).

"SEPARATION" button > Start the separation.

"Save results" button > Save last separated signals into variable "shat".

Filter length > Corresponds to the length of the separation filter in case that mu=1. This number times
number of input channels equals the dimension of the observation space to that the ICA algorithm 
is applied [0]. The higher the length, the higher the computational burden. Select between 3-40.

mu parameter > Parameter of Laguerre separating filters. mu=1 => FIR, mu~=1 => IIR. The values should
be from 0 to 2 [0,9].

ICA method > Method used for the ICA decomposition. EFICA exploits nonGaussianity of sources while 
BGSEP (formerly reffered to as BGL) utilizes their nonstationarity. BGSEP is faster then 
EFICA since only second-order statistics are used; see [4],[6]. Block EFICA combines both approaches [5], 
thus, is computationally more expensive. Analogously, BARBI combines block-stationarity
and block-spectral-diversity of sources [7]. The applicability of BARBI within T-ABCD is questionable: it might
fail as well as perform best; see [8].

Similarity measure > The criterion of similarity of independent components. Projections were described 
in [0-3]. GCC-PHAT coefficients were first used in [11].

Clustering method > The Hierarchical clst. algorithm [1,2] or the Fuzzy C-Means algorithm [3].

Weighting type > Normal [1], Toeplitz-structure-based [12], binary (hard) [1,2], fuzzy [3],
"constrained" (adjustable by user through AdjustWeights.m).

Level of weighting > The alpha parameter in (10) in [1,3] that controls "hardness" of weighting in the
fuzzy reconstruction mechanism. The higher the value, the higher interference suppression but at
the cost of higher spectral distortion. Default value is 2. The best Signal-to-Distortion ratio (SDR) 
is usually achieved for alpha close to 1.

# Subbands > The number of subbands. One means fullband separation (default). The subband processing 
approach was described in [10].

"Record 5 secs" button > Record 5 seconds of stereo signal from a standard audio device (works 
well with MS Windows XP; not tested elsewhere), and restart tddeconv with the recorded data.

Save results > Store the separated signals into 'shat' and the estimated microphone responses (images) into
'est_resp'.

Sampling frequency > sampling frequency for playback and recording 

Length of plot > range of x-axis in the "Mixed signals" plot

Offset and Length of block for ICA > The beginning and the length of block of data where ICA
decomposition is computed.

References
----------

[0] Z. Koldovský and P. Tichavský, "Time-Domain Blind Separation of Audio Sources on the basis 
of a Complete ICA Decomposition of an Observation Space", accepted for publication in 
IEEE Trans. on Speech, Audio and Language Processing, April 2010.

[1] Z. Koldovský and P. Tichavský, "Time-domain Blind Audio Source Separation Using 
Advanced Component Clustering and Reconstruction", The Joint 
Workshop on Hands-free Speech Communication and Microphone Arrays (HSCMA 2008), 
May 6-8, Trento, Italy, 2008.

[2] Z. Koldovský and P. Tichavský, "Time-Domain Blind Audio Source Separation Using 
Advanced ICA Methods", Proceedings of 8th Annual Conference of the International 
Speech Communication Association (Interspeech 2007), pp. 846-849, August 2007.

[3] J. Málek, Z. Koldovský, J. Žïánský and J. Nouza, "Enhancement of Noisy Speech 
Recordings via Blind Source Separation", Proceedings of the 9th Annual Conference 
of the International Speech Communication Association (Interspeech 2008), pp. 159-162, 
ISSN: 1990-9772, September 22-26, Brisbane, Australia, 2008.

EFICA

[4] Z. Koldovský, P. Tichavský and E. Oja, "Efficient Variant of Algorithm FastICA 
for Independent Component Analysis Attaining the Cramér-Rao Lower Bound", IEEE Trans. 
on Neural Networks, Vol. 17, No. 5, Sept 2006.

Block EFICA

[5] Z. Koldovský, J. Málek, P. Tichavský, Y. Deville, and S. Hosseini, "Blind 
Separation of Piecewise Stationary NonGaussian Sources", accepted for publication 
in Signal Processing, 2009.

BGSEP

[6] P. Tichavský and A. Yeredor,  Fast Approximate Joint Diagonalization Incorporating 
Weight Matrices, IEEE Tr. on Signal Processing, Vol. 57, No. 3, pp. 878-891, March 2009.

BARBI

[7] P. Tichavský, A. Yeredor, and Z. Koldovský, "A Fast Asymptotically Efficient Algorithm 
for Blind Separation of a Linear Mixture of Block-Wise Stationary Autoregressive Processes, 
ICASSP 2009, pp. 3133-3136, Taipei, Taiwan, April 2009. 

Further extensions and comparisons

[8] Z. Koldovský and P. Tichavský, "A Comparison of Independent Component and Independent 
Subspace Analysis Algorithms," EUSIPCO 2009 , pp. 1447-1451, Glasgow, Scotland, 
August 24-28, 2009.

[9] Z. Koldovský, P. Tichavský, and J. Málek, "Time-domain Blind Audio Source Separation
Method Producing Separating Filters of Generalized Feedforward Structure," Proc. of LVA/ICA 2010, 
St. Malo, France, Sept. 2010.

[10] Z. Koldovský, P. Tichavský, and J. Málek, "Subband Blind Audio Source Separation
Using a Time-Domain Algorithm and Tree-Structured QMF Filter Bank," Proc. of LVA/ICA 2010, 
St. Malo, France, Sept. 2010.

[11] J. Málek, Z. Koldovský, and P. Tichavský, "Adaptive Time-Domain Blind Separation of
Speech Signals, " Proc. of LVA/ICA 2010, St. Malo, France, Sept. 2010.

[12] Z. Koldovský, J. Málek, and P. Tichavský, "Blind Speech Separation in Time-Domain 
Using Block-Toeplitz Structure of Reconstructed Signal Matrices," Interspeech 2011, 
Florence, Italy, Aug. 2011. 


Links
-----

http://itakura.ite.tul.cz/zbynek/downloads.htm

http://si.utia.cas.cz/downloadPT.htm

Contacts
--------
zbynek.koldovsky_at_tul.cz
tichavsk_at_utia.cas.cz    
jiri.malek_at_tul.cz
(_at_ = @)

Enjoy it! :)
