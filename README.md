# Shuffle3DGRE
An open (via Pulseq) 3D gradient echo sequence supporting multiple Cartesian trajectory options for reducing physiological fluctuation artifacts 
====================================================================
% 
% General

% -------

%
% This scripts can be used to generate any of a family of 3D gradient echo
% sequences with different RF phase increments, TRs, and (multiple) TEs.
% The main purpose is to make different trajectories (on the 2D phase
% encoding grid) accessible.

%
% Standard user parameters are set in the 'Sequence params' section of the
% code. Most should be self explanatory (together with the comments).
%

% Shuffling/Scrambling of acquisition order

% -----------------------------------------
%

% A unique subset of parameters controls "Shuffling" the order in which the
% 2D phase encodes (denoted as PE and 3D) are sampled. Several predefined
% order schemes are available via the Setup.fOrdering parameter. The
% options include Local and Segmented Scrambling of the “ordered” case,
% i.e., line-by-line, (as in the reference [1]). The scrambling “window” is
% set by a positive cutoff frequency, where an empty cutoff frequency
% parameter means no scrambling, while a zero cutoff frequency means full
% randomization. The options also include a Generalized Hilbert curve order
% (reference [2]) and three spiral acquisition variants; all three of which
% can be combined with Local Scrambling.
%

% For more schemes, new functions have to be written. These function should
% be easily integrate here, as long as their input and output adhere to the
% "standard" (see Setup.fOrdering and the functions it can point to, for
% more information).
%

% For testing and comparing order schemes there is also an option to
% "interleave" two trajectories, i.e., alternate between sampling a readout
% using the PE & 3D from one order scheme and one using PE & 3D from a
% second (optional) order. This doubles the total acquisition time, but
% allows to "simultaneously" acquire two trajectories. This should reduce
% the uncertainty of different motion or physiological fluctuation between
% two separate acquisitions when trying to compare different trajectories.
% (However, interleaving two trajectories might affect the eddy current
% introduced, compared to separate scanning.)
% To enable this option just define a second ordering, i.e., use a
% non-empty Setup.fOrdering2 (and associated parameters).
%

% Automatic reconstruction on the scanner (Siemens interpreter v1.4.3)

% --------------------------------------------------------------------
%

% The generated sequences can be used with Siemens’ automatic image
% reconstruction. As of 2024.08.19, we tested the scripts using the Siemens
% v1.4.3 interpreter. Siemens v1.4.3 interpreter works for scans without
% acceleration by simply choosing 'ICE 3D' for 'Data handling' in the
% Special Card of the scanner UI.  However, for GRAPPA acceleration this
% interpreter does not always fill all parameters correctly. If fixed
% manually (using xbuilder), a retrospective reconstruction can be applied.

% The following parameter should be patched:

%   MEAS.ucEnableNoiseAdjust = 1
%   MEAS.sPat.lAccelFactPE   = value of Setup.AccelerationPE
%   MEAS.sPat.lAccelFact3D   = value of Setup.Acceleration3D
%   MEAS.sPat.lRefLinesPE    = value of Setup.RefLinesPE
%   MEAS.sPat.lRefLines3D    = value of Setup.RefLines3D
%   MEAS.sPat.ucPATMode      = 2
%   MEAS.sPat.ucRefScanMode  = 2
%   
%   YAPS.lFirstFourierLine      = minimal line acquired (Excluding noise scan)
%   YAPS.lFirstRefLine          = minimal reference line acquired
%   YAPS.lFirstFourierPartition = minimal partition (3D) acquired 
%                                 (Excluding noise scan)
%   YAPS.lFirstRefPartition     = minimal reference partition (3D) acquired
%

% For GRAPPA acceleration (along PE only) to be reconstructed automatically
% on the scanner, a small patch of the interpreter source code is required
% (as well as matching the GRAPPA parameters in the scanner UI and in the
% script). 

% For retro-reconstruction of elliptical scanning, several additional
% parameter updates are required, as well as patching the raw data itself
% (works on VE12U and probably some other IDEA versions). The steps below
% should be followed to reconstruct the images on the scanner:

% 1) Add to the last acquired line the MDH flag: MDH_LASTSCANINMEAS.

% 2) Patch the following parameters (in addition to any GRAPPA related
%    updates):

%      MEAS.sKSpace.ucEnableEllipticalScanning = 1
%      YAPS.alICEProgramPara[25]  = 8
%      DICOM.alICEProgramPara[25] = 8

% 3) Run retro-recon.
%

%
% References/Credit

% -----------------

% For background on Local Scrambling and some other scrambling see:
% [1]  https://doi.org/10.1002/mrm.29790
%

% The "Gilbert" trajectory (Generalized Hilbert curve) used here is a
% Matlab implementation of the Python code from:
% [2] https://github.com/jakubcerveny/gilbert/blob/master/gilbert2d.py
%

% Written by: Amir Seginer, Roni Avishay, Rita Schmidt

	
% NOTE: 'Shuffling' and 'Scrambling' may be used here indistinguishably.
 
% NOTE: sequence events (created by mr.makeXXX() functions) are denoted
%       here by variables with an initial small letter 'e'.

