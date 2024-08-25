% =========================================================================
% A Pulseq 3D gradient echo sequence with flexible Cartesian trajectories
% =========================================================================
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


%% Sequence params (user controlled)

% Mark start of preplimiary setups (and initialization and ...)
tStart = tic ;
fprintf(1, 'Preliminary setup ... ') ;

% NOTE: Relation between X/Y/Z and RO/PE/3D is defined below (look for
%       AxisRO/AxisPE/Axis3D)

% NOTE: [DEFINITIONS] to be passed to the interpreter are placed in
%       the structure Defs4Interpreter. (Defined towards the end of the
%       script.)
%       

clear Setup ; % clear any previous settings.

% ----------------------
% Standard 3D GRE params
% ----------------------

Setup.FOV_RoPe3d = [220e-3, 220e-3, 220e-3] ; % [m] RO x PE x 3D
% Setup.FOV_RoPe3d = [200, 220e-3, 220e-3] ; % [m] RO x PE x 3D
Setup.SlabThickness = Setup.FOV_RoPe3d(3) ; % [m] empty mean no slab selection, i.e. excite all
Setup.FlipAngleDeg = 15 ; % [Deg]
Setup.N_RO = 128 ; 
Setup.N_PE = 128 ; % <--- Siemens: must be rounded value of an integer percent of N_RO!!!
Setup.N_3D = 128 ; 
Setup.DimFast = '3D' ; % 'PE' or '3D'; which dimension is stepped through first
Setup.AccelerationPE = 3 ; % integer. > 1 - acceleration.
Setup.Acceleration3D = 1 ; % integer. > 1 - acceleration.
Setup.RefLinesPE = 24 ; % Number of fully sampled lines at center, along PE
Setup.RefLines3D = 24 ; % Number of fully sampled lines at center, along 3D
% Array(!) of TEs. Non-positive values will be interperted as using minimum
% TE possible (for that echo). Length of array is number of echos.
Setup.TE = -1 ; % [s] - non-positive values mean use shortest TE.
Setup.TR = 20e-3 ; % [s]
Setup.BWPerPixel = 260 ; % [Hz/pixel] updated later to obey ADC raster

% Define slab-selective excitation (in case it is requested)
% Original Definitions are based on Siemens a_gre (VE12U)
Setup.RFExciteSlabDuration = 2e-3 ; % [s]
Setup.RFExciteSlabBWTimeProd = 12.5 ; % [-] should update with B0 and flip angle?
Setup.RFExciteSlabNumSamples = 500 ; 

% Define non-selective excitation 
% Original Definitions are based on Siemens a_gre (VE12U)
Setup.RFExciteAllDuration = 0.1e-3 ; % [s]


% When acquiring with multiple TEs, should the RO gradient alternate in
% sign (faster) or always have the same sign (slower). Switching sign is
% refered to as bi-polar
Setup.bBipolarROGrads = true ;

% RF spoiling control
Setup.RFSpoilIncDeg = 117 ; % [deg]

% Spoiler control - one of two ways
%
% Spoiler area can be defined either explicitly in mT*us/m or relative to
% the net area of gradients required to achieve desired resolution. Two
% sets of variables are given at least one must be empty.
Setup.SpoilerArea_RO = [] ; % [mT*us/m]
Setup.SpoilerArea_PE = [] ; % [mT*us/m]
Setup.SpoilerArea_3D = [] ; % [mT*us/m]
Setup.SpoilerAreaFactor_RO = 2 ;
Setup.SpoilerAreaFactor_PE = 2 ;
Setup.SpoilerAreaFactor_3D = 2 ;
% % DEBUG/TESTING:
% warning('Zero spoiling along Y and Z for testing!') ;
% Setup.SpoilerAreaFactor_RO = 0 ;
% Setup.SpoilerAreaFactor_PE = 0 ;
% Setup.SpoilerAreaFactor_3D = 0 ;


% Dummy scans - minimum time(!) of dummy scans
Setup.MinDurDummyScans = 5 ; % [s]
% % DEBUG/TESTING:
% warning('MinDurDummyScans is set to zero/small for testing!') ;
% Setup.MinDurDummyScans = 0*TR ; % [s]


% ---------
% Shuffling
% ---------

% Cutoff frequency (if relevant): 
% -------------------------------

% - [] (empty) - no shuffling.
% - 0 - full randomization.
% - > 0 - Cutoff to use for scrambling (if relevant)
Setup.CutoffFreq = [] ; 1/20 ; % [Hz] 


% Set trajectory
% --------------

% set function handle to function used to order the trajectory.
% All functions accept the same inputs, although might not use them all.
Setup.fOrdering = @Ordering.Ordered_LocalShuffle ;
% Setup.fOrdering = @Ordering.SpiralRS_LocalShuffle ;
% Setup.fOrdering = @Ordering.Ordered_SegmentedShuffle ;
% Setup.fOrdering = @Ordering.SpiralSquare_LocalShuffle ;
% Setup.fOrdering = @Ordering.Gilbert ;
% Setup.fOrdering = @Ordering.Spiral_LocalShuffle ;
% Optional structure of extra parameters which are specific to the ordering
% function used.
Setup.OrderingExtraParamsStruct = [] ;
% Setup.OrderingExtraParamsStruct.bElliptic = true ;

% Random seed for shuffling
Setup.RandomSeed = 0 ;

% optional second defintion of a trajectory, for interleaving two
% trajectories. This will double the acquisition time, but allows to
% measure two trajectories/ordering practically simulatenously for a fair
% comparison. This is just for testing/comparing. Typically should not be
% used and fOrdering2 should be empty (or not defined at all).
% Interleaving means that the PE & 3D encoding of the acquisitions will
% alternate between those set by he first trajectory and those of the
% second trajectory.
Setup.fOrdering2 = [] ; % Empty - no 2nd trajectory to interleave.
% warning('A second trajectory (to be interleaved) is defined.')
% Setup.fOrdering2 = @Ordering.Ordered_LocalShuffle ;
Setup.OrderingExtraParamsStruct2 = Setup.OrderingExtraParamsStruct ;
Setup.RandomSeed2 = Setup.RandomSeed ;
Setup.CutoffFreq2 = [] ; % [] (empty) - no shuffling.
% warn user if interleaving is in effect.
if (~isempty(Setup.fOrdering2))
  warning('NOTE: InTeRlEaViNg of two trajectories is in place.')
end



% -------------
% Set system HW
% -------------

% Select HW system
[SystemNormal, ...
 SystemsStructArray, ...
 SystemNormalIdx] = Systems.SiemensTerraXR() ; % Terra (7T) defintions

% Set gradient limits for different elements of the sequence

% Default system to use
SystemDefault = SystemNormal ;
% Slab selection excitation including refocus grad
SystemSlabSelect = SystemNormal ;
% RO gradient during ADC
SystemRO = SystemNormal ;
% RO/PE/3D prephasers, rephase, and rewinder for monopolar-multiechos
SystemRePrePhase = SystemNormal ;
% Spoilers
SystemSpoilers = SystemNormal ;

% ----------------------------
% Set orienation (non-oblique)
% ----------------------------

% Set axes (X/Y/Z vs. RO/PE/3D - oblique not supported) 
% In Siemens interpreter the defintions here must agree with the
% 'Orientation mapping' setting.

% mapping of RO/PE/3D to X/Y/Z
Setup.AxisRO = 'x' ;
Setup.AxisPE = 'y' ;
Setup.Axis3D = 'z' ; % Slab/Paritions

% Flip or not X/Y/Z to match patient positive/negative directions. Usefull
% if reconstruction is done by system to get correct orientation of images.
Setup.SignCorr.x = -1 ;
Setup.SignCorr.y = -1 ;
Setup.SignCorr.z = -1 ;


% -------------
% Store results (not part of Setup struct)
% -------------

% Full path of .seq file to save results into. If SeqName is not defined,
% or empty, will not save a file.
% NOTE: SeqName is also used to store the name in the sequence file
SeqFolder = fullfile(pwd, 'SeqFiles') ; 
% Name of function we call to set the trajectory. Removes 'Ordering.' at 
% its start, if included.
fOrderingStr = regexprep(func2str(Setup.fOrdering), '^Ordering\.', '') ;
CutoffFreqStr = '' ;
if (~isempty(Setup.CutoffFreq))
  CutoffFreqStr = sprintf('_CutFreq%g', Setup.CutoffFreq) ;
end
% Sequence name without(!) '.seq' extension.
if (isfield(Setup, 'fOrdering2') && ~isempty(Setup.fOrdering2))
  fOrdering2Str = regexprep(func2str(Setup.fOrdering2), '^Ordering\.', '') ;
  CutoffFreqStr2 = '' ;
  if (~isempty(Setup.CutoffFreq2))
    CutoffFreqStr2 = sprintf('_CutFreq%g', Setup.CutoffFreq2) ;
  end
  SeqName = sprintf('3D_GRE_%dx%dx%d_PEx%d_3Dx%d_BW%d-%s%s-%s%s', ...
                    Setup.N_RO, Setup.N_PE, Setup.N_3D, ...
                    Setup.AccelerationPE, Setup.Acceleration3D, ...
                    Setup.BWPerPixel, ...
                    fOrderingStr, CutoffFreqStr, ...
                    fOrdering2Str, CutoffFreqStr2) ; 
else
  SeqName = sprintf('3D_GRE_%dx%dx%d_PEx%d_3Dx%d_BW%d-%s%s', ...
                    Setup.N_RO, Setup.N_PE, Setup.N_3D, ...
                    Setup.AccelerationPE, Setup.Acceleration3D, ...
                    Setup.BWPerPixel, fOrderingStr, CutoffFreqStr) ; 
end
% warning('SeqName is overwritten with a temp name!!!!')
% SeqName = 'test' ;
% warning('SeqName is overwritten to be empty!!!!')
% SeqName = [] ;

% -------------
% Misc.
% -------------

% Supress the plots at the end. (Saves time, but no feedback is given).
bSupressPulseqPlots = false ;

%% Initialize Actual Params

% Some parameters used will be different than those requested in the Setup
% structure. These actual parameters will be stored in the structure
% Actual.

% Initialize Actual to be the same as Setup
Actual = Setup ;

%%

% For Siemens Recon - ensure N_PE is consistent with UI
Actual.N_PE = round((round(Setup.N_PE/Setup.N_RO*100)/100) * Setup.N_RO) ;
if (Actual.N_PE ~= Setup.N_PE)
  warning(['N_PE updated to be consistent with Siemens UI (integer ' ...
           'percent of N_RO): %d --> %d'], ...
          Setup.N_PE, Actual.N_PE) ;
end



%% Initialize Sequence

% Moved here in case we need it early.
% (For example set IDs to some events that are used many times, without
%  change. This speeds up the adding of blocks to the sequence.)
Seq = mr.Sequence(SystemDefault) ;


%% Helper (anonymous) functions

% Checks if a number is practically an integer:
IsIntValue = @(n) abs(rem(n,1)) <= eps(1) ;


%% Round times to rasters (if necessary)

% -------------------------------------------------------------------------
% Convinient definitions
% -------------------------------------------------------------------------

% short names for quantities we'll need later.
ADCRaster = SystemDefault.adcRasterTime ; % [s] 
GradRaster = SystemDefault.gradRasterTime ; % [s]
RFRaster = SystemDefault.rfRasterTime ; % [s] 
BlockRaster = SystemDefault.blockDurationRaster ; % [s]

% Number of TEs defined
NumTEs = numel(Actual.TE) ;


% -------------------------------------------------------------------------
% Ensure TE & TR are on gradient raster time
% -------------------------------------------------------------------------

% Round TEs (supports multi TEs for multiple echos)
bTEPos = (Actual.TE > 0) ; % non-positive TEs mean use shortest TE
Actual.TE(bTEPos) = round(Setup.TE(bTEPos) / GradRaster) * GradRaster ;
% Warn user if this has an effect.
for TECounter = 1:NumTEs
  if (bTEPos(TECounter) && ...
      abs(Actual.TE(TECounter) - Setup.TE(TECounter)) > eps(GradRaster))
    warning('TE(%d) updated to be on the gradient raster: %f ms --> %f ms', ...
            TECounter, Setup.TE(TECounter)*1e3, Actual.TE(TECounter)*1e3) ;
  end
end

% Round user defined TR
Actual.TR = round(Setup.TR / GradRaster) * GradRaster ;
% Warn user if this has an effect.
if (abs(Actual.TR - Setup.TR) > eps(GradRaster))
  warning('TR updated to be on the gradient raster: %f ms --> %f ms', ...
          Setup.TR*1e3, Actual.TR*1e3) ;
end


% Mark end of preplimiary setups (and initialization and ...)
tEnd = toc(tStart) ;
fprintf(1, 'Done. (%g s)\n', tEnd) ;

%% Event definitions

% Mark start of event defintions
tStart = tic ;
fprintf(1, 'Event definitions ... ') ;


% ------------------------------------------
% Short(er) named variables, for convenience
% ------------------------------------------

% mapping of RO/PE/3D to X/Y/Z
AxisRO = Actual.AxisRO ;
AxisPE = Actual.AxisPE ;
Axis3D = Actual.Axis3D ; % Slab/Paritions

% Flip or not X/Y/Z to match patient positive/negative directions. Usefull
% if reconstruction is done by system to get correct orientation of images.
SignCorr = Actual.SignCorr ; % a structure


% -------------------------------------------------------------------------
% Excitation events - Based on Siemens SEQ::RF_NORMAL for a_gre (VE12U)
% -------------------------------------------------------------------------


% Derived parameters
FlipAngleRad = Actual.FlipAngleDeg * pi/180 ;
ExciteDwellTime = Actual.RFExciteSlabDuration / ...
                  Actual.RFExciteSlabNumSamples ;

% generate "events"
% We use the Hann filter for apodization (0.5) instead of the Hamming
% filter (0.46). The Hann filtered slice profile decays faster away from
% the slice, although the first side lobe of the Hamming case is smaller.
% (This also seems to be what Siemens does.)
if (isempty(Actual.SlabThickness)) % non-selective excitation
  [eExciteRF] = mr.makeBlockPulse(FlipAngleRad, SystemDefault, ...
                                 'duration', Actual.RFExciteAllDuration, ...
                                 'use', 'excitation') ;
  % Update Actual parameters
  Actual.RFExciteAllDuration = eExciteRF.shape_dur ;
else % slab selective excitation
  [eExciteRF, ...
   eExciteGrad, ...
   eExciteRefocGrad] = ...
                mr.makeSincPulse(FlipAngleRad, SystemSlabSelect, ...
                                'duration', Actual.RFExciteSlabDuration, ...
                                'timeBwProduct', Actual.RFExciteSlabBWTimeProd, ...
                                'dwell', ExciteDwellTime, ...
                                'apodization', 0.5, ... % 0.5 - Hann
                                'sliceThickness', Actual.SlabThickness, ...
                                'use', 'excitation') ;
  % Fix direction & sign (of two gradients)
  eExciteGrad.channel = Axis3D ;
  eExciteGrad.amplitude = SignCorr.(Axis3D) * eExciteGrad.amplitude ;
  eExciteGrad.flatArea = SignCorr.(Axis3D) * eExciteGrad.flatArea ;
  eExciteRefocGrad.channel = Axis3D ;
  eExciteRefocGrad.amplitude = SignCorr.(Axis3D) * eExciteRefocGrad.amplitude ;
  eExciteRefocGrad.flatArea = SignCorr.(Axis3D) * eExciteRefocGrad.flatArea ;

  % Update Actual parameters
  Actual.RFExciteSlabDuration = eExciteRF.shape_dur ;
  Actual.RFExciteSlabNumSamples = numel(eExciteRF.t) ;
end




% -------------------------------------------------------------------------
% RO: ADC and gradient
% -------------------------------------------------------------------------

% Start by defining the dwell (rounded to adcRasterTime)
ADCDwell = round(1/(Actual.N_RO*Actual.BWPerPixel)/ADCRaster)*ADCRaster ; % [s] 

% Get ADC duration and the duration rounded up to the gradient raster
ADCDuration = Actual.N_RO * ADCDwell ;
% Update actual BW per pixel used
Actual.BWPerPixel = 1/ADCDuration ; % [Hz/pixel]

% Find the ADC duration rounded to the twice(!) gradient raster. (We want
% to ensure the center of the flat area is on the gradient raster, so TEs
% can fall on the gradient raster
% NOTE: We implicitly assume below (with the ADCDelay) that the RF raster
%       is finer than the gradient raster. (That we do not have to worry
ADCDurGradRaster = ceil(ADCDuration / GradRaster/2) * 2*GradRaster ; % [s]

% define RO gradient event + fix sign?
ROAmp = SignCorr.(AxisRO) / (ADCDwell*Actual.FOV_RoPe3d(1)) ; % [Hz/m]
eRO = mr.makeTrapezoid(AxisRO, SystemRO, ...
                       'flatTime', ADCDurGradRaster, ...
                       'amplitude', ROAmp) ;
% define ADC event inc. delay relative to start of RO gradient.
% Note: delay is in RF raster time, not ADC raster time.
ADCDelay = eRO.riseTime + (ADCDurGradRaster - ADCDuration)/2 ; % [s]
ADCDelay = ceil(ADCDelay/RFRaster) * RFRaster ; % Ensure on RF(!) raster.
eADC = mr.makeAdc(Actual.N_RO, SystemDefault, ...
                  'dwell', ADCDwell , ...
                  'delay', ADCDelay) ; 
% Round up total time of eADC to the block raster time. In rare(?) cases
% the total time may be longer than the RO gradient and not be on the block
% raster time. This will be used to slightly extend the block duration to
% be on the block raster.
ROADCGrossDurBlockRaster = ceil(mr.calcDuration(eADC) / BlockRaster) * ...
                           BlockRaster ; % [s]
% Create a matching "delay" (played in parallel to the other events in the
% block, so will just ensure its duration.
eROADCDelay = mr.makeDelay(ROADCGrossDurBlockRaster) ;

% Set ADC labels (Siemens MDH flags)
% For fast reuse also define an ID for each.

% Marks GRAPPA line which is used as reference (center of k-space).
eLabelRef = mr.makeLabel('SET','REF', true) ;
  eLabelRef.id = Seq.registerLabelEvent(eLabelRef) ;
% Marks GRAPPA line which is used as reference (center of k-space) AND(!)
% also would be acquired if no reference lines were acquired (e.g., SENSE),
% i.e., would also be acquired due to the subsampling only.
eLabelRefAndImage = mr.makeLabel('SET','IMA', true) ;
  eLabelRefAndImage.id = Seq.registerLabelEvent(eLabelRefAndImage) ;
% Negates eLabelRef.
eLabelNoRef = mr.makeLabel('SET','REF', false) ;
  eLabelNoRef.id = Seq.registerLabelEvent(eLabelNoRef) ;
% Negates eLabelRefAndImage.
eLabelNoRefAndImage = mr.makeLabel('SET','IMA', false) ;
  eLabelNoRefAndImage.id = Seq.registerLabelEvent(eLabelNoRefAndImage) ;

% Labels to mark if ADC should be time reversed. Used in bipolar multi-echo
% acquisitions.
eLabelRevOff = mr.makeLabel('SET', 'REV', false) ;
  eLabelRevOff.id = Seq.registerLabelEvent(eLabelRevOff) ;
eLabelRevOn = mr.makeLabel('SET', 'REV', true) ;
  eLabelRevOn.id = Seq.registerLabelEvent(eLabelRevOn) ;


% Labels to mark which trajectory was sampled, in case two were sampled
% (interleaved).
eLabelTraj1 = mr.makeLabel('SET', 'SET', 0) ;
  eLabelTraj1.id = Seq.registerLabelEvent(eLabelTraj1) ;
eLabelTraj2 = mr.makeLabel('SET', 'SET', 1) ;
  eLabelTraj2.id = Seq.registerLabelEvent(eLabelTraj2) ;

% Labels to set echo, in case of multi echo acquisition.

% Allocate memory first by creating the last echo label (probably not very
% important).
clear eLabelEchos ; % clear from previous (possible) run
% set last echo in array (to initialize structure array)
eLabelEchos(NumTEs) = mr.makeLabel('SET', 'ECO', NumTEs-1) ; 
% set all remaining echo labels
for TECounter = 1:(NumTEs-1)
  eLabelEchos(TECounter) = mr.makeLabel('SET', 'ECO', TECounter-1) ;
end
% Add IDs to all labels (done separately after previous for-loop because
% once one ID is set the output of mr.makeLabel no longer has the same
% structure as eLabelEchos(n) (does not contain the 'id' field).
for TECounter = 1:NumTEs
  eLabelEchos(TECounter).id = Seq.registerLabelEvent(eLabelEchos(TECounter)) ;
end


% -------------------------------------------------------------------------
% Max. RO prephase
% -------------------------------------------------------------------------

% Max absolute area of the RO prephase(!) gradient
GradPrephaseMaxAbsAreaEO = abs(eRO.area/2) ; % [1/m] ;


% -------------------------------------------------------------------------
% Max. PE prehase
% -------------------------------------------------------------------------

% % DEBUG/TESTING:
% warning('Overriding N_PE, N_3D to a small value (for testing)') ;
% Actual.N_PE = 2 ; Actual.N_3D = 2 ;

% Center when counting from 1
KPECenterIdx = Ordering.Utils.FFTCenterIndex(Actual.N_PE) ; 

% Find start and end indices of y phase encoding.
KPEIdxMax = Actual.N_PE - KPECenterIdx ;
KPEIdxMin = KPEIdxMax - Actual.N_PE + 1 ;
% Find largest index of the two, in absolute value. (The negative one,
% if they are not equal.)
KPEAbsIdxMax = max(abs([KPEIdxMax, KPEIdxMin])) ;
% Max absolute area of the prephase gradient
GradMaxAbsAreaPE = KPEAbsIdxMax/Actual.FOV_RoPe3d(2) ; % [1/m] ;

% -------------------------------------------------------------------------
% Max. 3D prehase
% -------------------------------------------------------------------------

% Center when counting from 1
K3DCenterIdx = Ordering.Utils.FFTCenterIndex(Actual.N_3D) ; 

% Find start and end indices of z phase encoding.
K3DIdxMax = Actual.N_3D - K3DCenterIdx ;
K3DIdxMin = K3DIdxMax - Actual.N_3D + 1 ;
% Find largest index of the two, in absolute value. (The negative one,
% if they are not equal.)
K3DAbsIdxMax = max(abs([K3DIdxMax, K3DIdxMin])) ;
% Max absolute area of the prephase gradient
GradMaxAbsArea3D = K3DAbsIdxMax/Actual.FOV_RoPe3d(3) ; % [1/m] ;


% -------------------------------------------------------------------------
% Single total prephaser event (supports any oblique)
% -------------------------------------------------------------------------

% Our prephasers (RO, PE, and 3D) are going to have the same timings (ramp
% up, flat duration, and ramp down) and will also supprt any oblique
% direction. For this we'll first define a trapazoid with an area that
% combines all the prephasers (combines areas via a root of the sum of
% squares). Once we have that we'll defined three separate gradients with
% the same timing with the respective area (or max area).

% We start by defining a gradient that supports the combined area
% (root-sum-of-squares) of all directions. The AxisRO direction used is
% arbitrary.
PrephaseMaxAbsArea = sqrt(GradPrephaseMaxAbsAreaEO^2 + GradMaxAbsAreaPE^2 + ...
                           GradMaxAbsArea3D^2) ;
ePrephase = mr.makeTrapezoid(AxisRO, SystemRePrePhase, ...
                             'area', PrephaseMaxAbsArea) ;

% -------------------------------------------------------------------------
% PE prephaser/rephaser
% -------------------------------------------------------------------------

% Define PE prephase (Amplitude Will be updated later)
ePrephasePE = ePrephase ;
ePrephasePE.channel = AxisPE ;

% Define PE rephaser. (Amplitude Will be updated later according to actual
% ePrephasePE used)
eRephasePE = ePrephasePE ;

% The amplitude of the prephase gradient for 1 k-space step, for the
% gradient just defined. (We will use an amplitude which is an integer
% multiple of this.)
PrePhaseAmpStepPE = SignCorr.(AxisPE) * ...
                    (GradMaxAbsAreaPE/PrephaseMaxAbsArea * ...
                    ePrephasePE.amplitude) / KPEAbsIdxMax ; % [Hz/m]

% -------------------------------------------------------------------------
% 3D prephaser/rephaser
% -------------------------------------------------------------------------

% Define 3D prephase (Amplitude Will be updated later)
ePrephase3D = ePrephase ;
ePrephase3D.channel = Axis3D ;

% Define 3D rephaser. (Amplitude Will be updated later according to actual
% ePrephase3D used)
eRephase3D = ePrephase3D ;

% The amplitude of the prephase gradient for 1 k-space step, for the
% gradient just defined. (We will use an amplitude which is an integer
% multiple of this.)
PrePhaseAmpStep3D = SignCorr.(Axis3D) * ...
                    (GradMaxAbsArea3D/PrephaseMaxAbsArea * ...
                    ePrephase3D.amplitude) / K3DAbsIdxMax ; % [Hz/m]

% -------------------------------------------------------------------------
% RO prephaser/rephaser (and TE filler delays)
% -------------------------------------------------------------------------

% Define RO prephaser event
% Its area should be half the area of the RO gradient (we shifted the ADC
% to be centered at the center of the RO gradient).
ePrephaseRO = ePrephase ;
ePrephaseRO.channel = AxisRO ;
% Set ampilitude. Because it depends on eRO, whose sign has already been
% corrected, we do not need to correct sign here.
ePrephaseRO.amplitude = (-eRO.area/2) * ...
                       ePrephase.amplitude/PrephaseMaxAbsArea ;

% Define RO rephaser
eRephaseRO = ePrephaseRO ;
eRephaseRO.amplitude = -ePrephaseRO.amplitude ;

% Define RO prephaser for multiple mono-polar TEs (bBipolarROGrads is
% false). In this case after the acquisition RO gradient we have to negate
% it before we can run it again for the next TE.
eROUndo = mr.makeTrapezoid(AxisRO, SystemRePrePhase, ...
                            'area', -eRO.area) ;

% % DEBUG/TESTING:
% warning('Nulling RO gradients prephase, rephase, and RO itself')
% eRO.amplitude = 0 ;
% ePrephaseRO.amplitude = 0 ;
% eRephaseRO.amplitude = 0 ;


% -------------------------------------------------------------------------
% Spoiler(s)
% -------------------------------------------------------------------------

% Similar to prephaser/rephasers we will first find the area in each
% direction, then combine as root of sum of squares, then define a single
% gradient to cover that and finally make copies of per direction with
% matching amplitudes.

% RO spoiler area
% -------------------

% For conveniece later on
if (~isfield(Actual, 'SpoilerArea_RO'))
  Actual.SpoilerArea_RO = [] ; 
end
if (~isfield(Actual, 'SpoilerAreaFactor_RO'))
  Actual.SpoilerAreaFactor_RO = [] ; 
end

% How is spoiler defined? Absolute area or relative area
bSpoilerAreaDefinedRO = ~isempty(Actual.SpoilerArea_RO) ;
bSpoilerFactorDefinedRO = ~isempty(Actual.SpoilerAreaFactor_RO) ;


% Set area of spoiler
if (bSpoilerAreaDefinedRO && ~bSpoilerFactorDefinedRO) 
  % Area explcitly defined in mT*us/m. Translate to Hz*s/m
  SpoilAreaRO = Actual.SpoilerArea_RO * 1e-3 * SystemDefault.gamma ; % [Hz*s/m]
  % Ensure 
elseif (bSpoilerFactorDefinedRO && ~bSpoilerAreaDefinedRO )
  SpoilAreaRO = Actual.SpoilerAreaFactor_RO * Actual.N_RO/Actual.FOV_RoPe3d(1) ; % [Hz*s/m = 1/m]
elseif (~bSpoilerAreaDefinedRO && ~bSpoilerFactorDefinedRO)
%   error(['No RO spoiler is defined. Either ''SpoilerArea_RO'' or ' ...
%          '''SpoilerAreaFactor_RO'' must be defined'])
  warning(['No RO spoiler is defined. Neither ''SpoilerArea_RO'' nor ' ...
         '''SpoilerAreaFactor_RO''. Assuming zero.'])
  SpoilAreaRO = 0 ;
else % both are defined
  error(['RO spoiler is defined twice, both in ''SpoilerArea_RO'' and ' ...
         'in ''SpoilerAreaFactor_RO''. Set one to be empty.'])
end


% PE spoiler area
% -------------------

% For conveniece later on
if (~isfield(Actual, 'SpoilerArea_PE'))
  Actual.SpoilerArea_PE = [] ; 
end
if (~isfield(Actual, 'SpoilerAreaFactor_PE'))
  Actual.SpoilerAreaFactor_PE = [] ; 
end

% How is spoiler defined? Absolute area or relative area
bSpoilerAreaDefinedPE = ~isempty(Actual.SpoilerArea_PE) ;
bSpoilerFactorDefinedPE = ~isempty(Actual.SpoilerAreaFactor_PE) ;



% Set area of spoiler
if (bSpoilerAreaDefinedPE && ~bSpoilerFactorDefinedPE) 
  % Area explcitly defined in mT*us/m. Translate to Hz*s/m
  SpoilAreaPE = Actual.SpoilerArea_PE * 1e-3 * SystemDefault.gamma ; % [Hz*s/m]
elseif (bSpoilerFactorDefinedPE && ~bSpoilerAreaDefinedPE )
  SpoilAreaPE = Actual.SpoilerAreaFactor_PE * Actual.N_PE/Actual.FOV_RoPe3d(2) ; % [Hz*s/m = 1/m]
elseif (~bSpoilerAreaDefinedPE && ~bSpoilerFactorDefinedPE)
%   error(['No PE spoiler is defined. Either ''SpoilerArea_PE'' or ' ...
%          '''SpoilerAreaFactor_PE'' must be defined'])
  warning(['No PE spoiler is defined. Neither ''SpoilerArea_PE'' nor ' ...
         '''SpoilerAreaFactor_PE''. Assuming zero.'])
  SpoilAreaPE = 0 ;
else % both are defined
  error(['PE spoiler is defined twice, both in ''SpoilerArea_PE'' and ' ...
         'in ''SpoilerAreaFactor_PE''. Set one to be empty.'])
end


% 3D spoiler area
% -------------------

% For conveniece later on
if (~isfield(Actual, 'SpoilerArea_3D'))
  Actual.SpoilerArea_3D = [] ; 
end
if (~isfield(Actual, 'SpoilerAreaFactor_3D'))
  Actual.SpoilerAreaFactor_3D = [] ; 
end

% How is spoiler defined? Absolute area or relative area
bSpoilerAreaDefined3D = ~isempty(Actual.SpoilerArea_3D) ;
bSpoilerFactorDefined3D = ~isempty(Actual.SpoilerAreaFactor_3D) ;


if (bSpoilerAreaDefined3D && ~bSpoilerFactorDefined3D) 
  % Area explcitly defined in mT*us/m. Translate to Hz*s/m
  SpoilArea3D = Actual.SpoilerArea_3D * 1e-3 * SystemDefault.gamma ; % [Hz*s/m]
elseif (bSpoilerFactorDefined3D && ~bSpoilerAreaDefined3D )
  SpoilArea3D = Actual.SpoilerAreaFactor_3D * Actual.N_3D/Actual.FOV_RoPe3d(3) ; % [Hz*s/m = 1/m]
elseif (~bSpoilerAreaDefined3D && ~bSpoilerFactorDefined3D)
%   error(['No 3D spoiler is defined. Either ''SpoilerArea_3D'' or ' ...
%          '''SpoilerAreaFactor_3D'' must be defined'])
  warning(['No 3D spoiler is defined. Neither ''SpoilerArea_3D'' nor ' ...
         '''SpoilerAreaFactor_3D''. Assuming zero.'])
  SpoilArea3D = 0 ;
else % both are defined
  error(['3D spoiler is defined twice, both in ''SpoilerArea_3D'' and ' ...
         'in ''SpoilerAreaFactor_3D''. Set one to be empty.'])
end


% Single total spoiler event (supports any oblique)
% -------------------------------------------------

% We start by defining a gradient that supports the combined area
% (root-sum-of-squares) of all directions. The AxisRO direction used is
% arbitrary.
SpoilMaxAbsArea = sqrt(SpoilAreaRO^2 + SpoilAreaPE^2 + ...
                       SpoilArea3D^2) ;
if (SpoilMaxAbsArea == 0 )
  % Mark that no spoiler is used
  bSpoilers = false ;
else
  % Mark that spoilers are used
  bSpoilers = true ;

  % Define dummy "RO" spoiler event
  eSpoil = mr.makeTrapezoid(AxisRO, SystemSpoilers, ...
                            'area', SpoilMaxAbsArea) ;
  SpoilDuration = eSpoil.riseTime + eSpoil.flatTime + eSpoil.fallTime ;
  
  % RO spoiler event
  % --------------------
  eSpoilRO =  mr.makeTrapezoid(AxisRO, SystemSpoilers, ...
                               'duration', SpoilDuration, ...
                               'riseTime', eSpoil.riseTime, ...
                               'fallTime', eSpoil.fallTime, ...
                               'area', SignCorr.(AxisRO) * SpoilAreaRO) ;
  
  
  % PE spoiler event
  % --------------------
  eSpoilPE =  mr.makeTrapezoid(AxisPE, SystemSpoilers, ...
                               'duration', SpoilDuration, ...
                               'riseTime', eSpoil.riseTime, ...
                               'fallTime', eSpoil.fallTime, ...
                               'area', SignCorr.(AxisPE) * SpoilAreaPE) ;
  
  % 3D spoiler event
  % --------------------
  eSpoil3D =  mr.makeTrapezoid(Axis3D, SystemSpoilers, ...
                               'duration', SpoilDuration, ...
                               'riseTime', eSpoil.riseTime, ...
                               'fallTime', eSpoil.fallTime, ...
                               'area', SignCorr.(Axis3D) * SpoilArea3D) ;

end

  
% -------------------------------------------------------------------------
% TR filler delay
% -------------------------------------------------------------------------

% Initialize a zero delay that will be modified later on
% to ensure we get the desired TR (if possible).
eTRFill = mr.makeDelay(1) ; % dummy time. (Delay of zero is not allowed)
eTRFill.delay = 0 ; % [s] force delay of zero.



% Mark end of event defintions
tEnd = toc(tStart) ;
fprintf(1, 'Done. (%g s)\n', tEnd) ;


%% Set PE & 3D ordering

% Mark start of generating ordering scheme
tStart = tic ;
fprintf(1, 'Generating ordering scheme ... ') ;


% Is there (valid) acceleration along PE or 3D
bAccelerationPE = false ; % intialize no acceleration along PE
bAcceleration3D = false ; % intialize no acceleration along 3D
if (isfield(Actual, 'AccelerationPE') && ~isempty(Actual.AccelerationPE))
  if (~IsIntValue(Actual.AccelerationPE) || Actual.AccelerationPE < 1)
    error('AccelerationPE is expected to be a positive integer.')
  end
  if (Actual.AccelerationPE > 1)
    bAccelerationPE = true ;
  end
end
if (isfield(Actual, 'Acceleration3D') && ~isempty(Actual.Acceleration3D))
  if (~IsIntValue(Actual.Acceleration3D) || Actual.Acceleration3D < 1)
    error('Acceleration3D is expected to be a positive integer.')
  end
  if (Actual.Acceleration3D > 1)
    bAcceleration3D = true ;
  end
end

% Switch from PE and 3D to fast (#1) and slow (#2) dimensions.
% (Which is fast/slow depends on 'DimFast'.)

% initialize
Acceleration1 = 1 ;
RefLines1 = [] ;
Acceleration2 = 1 ;
RefLines2 = [] ;


if (strcmpi(Actual.DimFast, 'PE')) % Step through PE first (and then 3D)

  % stepped through first
  N1 = Actual.N_PE ; 
  K1IdxMin = KPEIdxMin ;
  K1IdxMax = KPEIdxMax ;
  K1CenterIdx = KPECenterIdx ;
  bAcceleration1 = bAccelerationPE ;
  if (bAcceleration1)
    Acceleration1 = Actual.AccelerationPE ;
    RefLines1 = Actual.RefLinesPE ;
  end

  % stepped through second
  N2 = Actual.N_3D ; 
  K2IdxMin = K3DIdxMin ;
  K2IdxMax = K3DIdxMax ;
  K2CenterIdx = K3DCenterIdx ;
  bAcceleration2 = bAcceleration3D ;
  if (bAcceleration2)
    Acceleration2 = Actual.Acceleration3D ;
    RefLines2 = Actual.RefLines3D ;
  end

elseif (strcmpi(Actual.DimFast, '3D'))  % Step through 3D first (and then PE)

  % stepped through first
  N1 = Actual.N_3D ; 
  K1IdxMin = K3DIdxMin ;
  K1IdxMax = K3DIdxMax ;
  K1CenterIdx = K3DCenterIdx ;
  bAcceleration1 = bAcceleration3D ;
  if (bAcceleration1)
    Acceleration1 = Actual.Acceleration3D ;
    RefLines1 = Actual.RefLines3D ;
  end

  % stepped through second
  N2 = Actual.N_PE ; 
  K2IdxMin = KPEIdxMin ;
  K2IdxMax = KPEIdxMax ;
  K2CenterIdx = KPECenterIdx ;
  bAcceleration2 = bAccelerationPE ;
  if (bAcceleration2)
    Acceleration2 = Actual.AccelerationPE ;
    RefLines2 = Actual.RefLinesPE ;
  end

else
  error(['DimFast is ''%s'' but should be either ''PE'' or ''3D'' ' ...
         '(case insensitive). Aborting'], Actual.DimFast) ;
end

% Set GRAPPA masks (even if no acceleration is used)
if (bAcceleration1 || bAcceleration2)
  [bSample, ...
   bRefFull, ...
   bImagingAndRefFull] = ...
                   Ordering.Utils.GrappaMasks(N1, N2, ...
                                              Acceleration1, Acceleration2, ...
                                              RefLines1, RefLines2) ;
else
  [bSample, ...
   bRefFull, ...
   bImagingAndRefFull] = Ordering.Utils.GrappaMasks(N1, N2, 1, 1) ;
end

% Set order of sampling
% Set handle fOrdering to get the desired order.
% NOTE: IOut and JOut are zero for k=0.

% product of cutoff frequency and TR
TRtimesCutoffFreq = Actual.CutoffFreq * Actual.TR ;
if (isfield(Actual, 'fOrdering2') && ...
    isa(Actual.fOrdering2, 'function_handle'))
  % fOrdering2 is defined so we are interleaving two trajectories. Thus
  % effectively the TR of each trajectory is doubled (we sample from each
  % trajectory only every second TR).
  TRtimesCutoffFreq = Actual.CutoffFreq * (2*Actual.TR) ;
end
[orderMat, ...
 IOut, JOut, ...
 SampledReorder] = Actual.fOrdering(N1, N2, ...
                                    TRtimesCutoffFreq, ...
                                    bSample, ...
                                    Actual.OrderingExtraParamsStruct, ...
                                    Actual.RandomSeed) ;

% Extract subset of Ref/ImagingAndRef from actually sampled samples to
% match the length of IOut and Jout. Recall that IOut and JOut are possibly
% reorderd relative to original bRefFull and bImagingAndRefFull, so we have
% to reorder them and cut them short. This is what we have orderOut for.
bRef = bRefFull(bSample(:)) ;
bRef = bRef(SampledReorder(:)) ;
bImagingAndRef = bImagingAndRefFull(bSample(:)) ;
bImagingAndRef = bImagingAndRef(SampledReorder(:)) ;


% Should we interleave two trajectories?
bInterleaveTrajectories = false ;
if (isfield(Actual, 'fOrdering2') && ~isempty(Actual.fOrdering2))
  % If we reached here we are interleaving two trajectories, so the
  % effective TR of each is doubled. (Same as for first trajectory above.)
  TRtimesCutoffFreq2 = Actual.CutoffFreq2 * (2*Actual.TR) ;

  % Generate trajectory #2
  [orderMat2, ...
   IOut2,JOut2, ...
   SampledReorder2] = Actual.fOrdering2(N1, N2, ...
                                        TRtimesCutoffFreq2, ...
                                        bSample, ...
                                        Actual.OrderingExtraParamsStruct2, ...
                                        Actual.RandomSeed2) ;
  % Ensure trajectory is of the same length as previous one
  if (numel(IOut) ~= numel(IOut2))
    error(['The lengthd of the two trajectories to be interleaved are ' ...
           'different (#1: %d and #2: %d). Aborting.'], ...
           numel(IOut), numel(IOut2)) ;
  end
  
  % Now let us interleave the two
  VecTemp = zeros(numel(IOut) + numel(IOut2), 1) ;
  % Update IOut
  VecTemp(1:2:end) = IOut(:) ;
  VecTemp(2:2:end) = IOut2(:) ;
  IOut = VecTemp ;
  % Update JOut
  VecTemp(1:2:end) = JOut(:) ;
  VecTemp(2:2:end) = JOut2(:) ;
  JOut = VecTemp ;

  % Set matching Ref/ImagingAndRef from actually sampled samples to
  bRef2 = bRefFull(bSample(:)) ;
  bRef2 = bRef2(SampledReorder2(:)) ;
  bImagingAndRef2 = bImagingAndRefFull(bSample(:)) ;
  bImagingAndRef2 = bImagingAndRef2(SampledReorder2(:)) ;

  % Inteleave Ref/ImagingAndRef of the two trajectories
  % Update bRef
  VecTemp(1:2:end) = bRef(:) ;
  VecTemp(2:2:end) = bRef2(:) ;
  bRef = VecTemp ;
  % Update bImagingAndRef
  VecTemp(1:2:end) = bImagingAndRef(:) ;
  VecTemp(2:2:end) = bImagingAndRef2(:) ;
  bImagingAndRef = VecTemp ;

  % done with VecTemp
  clear VecTemp ;

  % Mark that we are interleaving two trajectories (set different counter
  % in the header for them).
  bInterleaveTrajectories = true ;
end


% conversion to k-space representation format (indices are zero for k=0).
IOut = IOut-Ordering.Utils.FFTCenterIndex(N1);
JOut = JOut-Ordering.Utils.FFTCenterIndex(N2);



% Get list of PE and 3D indices to sample: return from using fast dimension
% (1) and slow dimension (2).
if (strcmpi(Actual.DimFast, AxisPE)) % Step through Y (#1) first and then Z (#2)
  PE3DOrder = [IOut(:), JOut(:)] ;
else
  PE3DOrder = [JOut(:), IOut(:)] ;
end

% % DEBUG: Show bRef and bImagingAndRef are marked correctly.
% figure ; 
%   plot(PE3DOrder(:,1), PE3DOrder(:,2), '.') ;
%   hold on
%   plot(PE3DOrder(bRef(:),1), PE3DOrder(bRef(:),2), 'o') ;
%   plot(PE3DOrder(bImagingAndRef(:),1), PE3DOrder(bImagingAndRef(:),2), '^', ...
%        'MarkerSize', 12) ;


% Mark end of generating ordering scheme
tEnd = toc(tStart) ;
fprintf(1, 'Done. (%g s)\n', tEnd) ;

%% Create sequence cell/structre

% Mark start of building sequence blocks
tStart = tic ;
fprintf(1, 'Building sequence ... \n') ;


% Initialize variables
RFSpoilInc = Actual.RFSpoilIncDeg * pi/180 ; % [rad]
% What is the actual minimum TR possible for current sequence (found when
% generating the sequence)
MinTRActual = 0 ;

% translate FOV from RO/PE/3D to X/Y/Z:
[~, FOVXYZOrder] = sort([AxisRO, AxisPE, Axis3D], 'ascend') ;
Actual.FOV = Actual.FOV_RoPe3d(FOVXYZOrder) ;


% How many dummy-scans do we want
NumDummyScans = ceil(Actual.MinDurDummyScans/Actual.TR) ;

% Initialize sequence
% Seq = mr.Sequence(SystemDefault) ;  % moved to ~start
RFSpoilPhase = eExciteRF.phaseOffset ; % [rad]
RFSpoilPhaseShift = 0 ; % [rad] Next step of RF phase (updates)

% Start with a noise scan - only on the first repetition.
% The line and partition of the noise scan must be on one of the "imaging"
% scans. For simplicity we set it as the center of k-space because that
% is always(?) sampled.

% Round up eADC total duration (including delays and dead time) to be on
% the block duration raster time. Hopefully in the future we can remove the
% extra dummy delay below of duration NoiseBlockDur. (This delay is in
% parallel to the rest of the events in the block, not before them.)
NoiseBlockDur = BlockRaster * ...
               ceil(mr.calcDuration(eADC)/BlockRaster - 1e-6) ;
Seq.addBlock(eADC, ...
             mr.makeDelay(NoiseBlockDur), ... % Temp bug fix!!!
             mr.makeLabel('SET', 'ONCE', 1), ... % only on 1st repetition.
             mr.makeLabel('SET', 'LIN', KPECenterIdx - 1), ...
             mr.makeLabel('SET', 'PAR', K3DCenterIdx - 1), ...
             mr.makeLabel('SET', 'NOISE', true), ...
             eLabelNoRef, ...
             eLabelNoRefAndImage) ;
Seq.addBlock(mr.makeLabel('SET', 'NOISE', false), ...
             mr.makeLabel('SET', 'ONCE', 0)) ; % reset 'ONCE'


% Label used to mark which trajectory we use, in case we interleave two.
% We initialize it here for the first one (in case there is no
% interleaving).
eLabelTraj = eLabelTraj1 ;

% Label used to mark if we should reverse the ADC (in bi-polar multi-echo,
% each time we flip the sign of the RO gradient).
eLabelReverse = eLabelRevOff ; % Don't reverse the first time

% Loop over TR: non-positive counter values mean dummy scan (ky=kz=0). 
% Otherwise, advance in the PE3DOrder table of phase encodes (PE & 3D) 
% NOTE: REFERENCE SCANS NOT SUPPORTED YET!!!!
for TRCounter = (-NumDummyScans+1):size(PE3DOrder, 1)


  % -----------------------------------------------------------------------
  % Preliminaries
  % -----------------------------------------------------------------------
  
  % Reset duration of current TR
  TimeInTR = 0 ; % [s]

  % Set y (PE) and z (3D) indices
  if (TRCounter <= 0) % dummy scans)
    KPEIdx = 0 ; % k = 0 for PE
    K3DIdx = 0 ; % k = 0  for 3D
    bADCOn = false ; % should signal be acquired?
  else
    KPEIdx = PE3DOrder(TRCounter, 1) ;
    K3DIdx = PE3DOrder(TRCounter, 2) ;
    bADCOn = true ; % should signal be acquired?
  end

  % -----------------------------------------------------------------------
  % Update RF spoiling phase to use in this TR
  % -----------------------------------------------------------------------

  % update RF phase to use this round
  RFSpoilPhase = RFSpoilPhase + RFSpoilPhaseShift ; % [rad]

  % upadte RFSpoilPhaseShift for next round ;
  RFSpoilPhaseShift = RFSpoilPhaseShift + RFSpoilInc ; % [rad]
  
  % -----------------------------------------------------------------------
  % Excitation block (without refocusing grad.)
  % -----------------------------------------------------------------------
  
  % update RF spoiling phase of RF pulse
  eExciteRF.phaseOffset = RFSpoilPhase ; % [rad]

  % Add excitation block
  if (isempty(Actual.SlabThickness)) % non-selective excitation
    Seq.addBlock(eExciteRF) ;
  else % slab selective excitation
    Seq.addBlock(eExciteRF,eExciteGrad) ;
  end

  %

  % Update duration within TR
  TimeInTR = TimeInTR + Seq.blockDurations(end) ;

  % Initialize time from TE: time between RF center (from start of block)
  % to end of the block.
  TimeFromExcite = Seq.blockDurations(end) - ...
                   (eExciteRF.delay + eExciteRF.shape_dur/2) ; % [s]


  % -----------------------------------------------------------------------
  % Slab refocusing(?) + prephasing block
  % -----------------------------------------------------------------------
  
  % Update y prephaser amplitude
  ePrephasePE.amplitude =  KPEIdx * PrePhaseAmpStepPE ; % [Hz/m]
  % Update z prephaser amplitude
  ePrephase3D.amplitude =  K3DIdx * PrePhaseAmpStep3D ; % [Hz/m]

  % Add block for slab refocusing?
  if (~isempty(Actual.SlabThickness)) % slab selective excitation
    Seq.addBlock(eExciteRefocGrad) ;

    % Update duration within TR
    TimeInTR = TimeInTR + Seq.blockDurations(end) ;
    % Update time from "excitation"
    TimeFromExcite = TimeFromExcite + Seq.blockDurations(end) ;
  end

  % Add block of prephasers
  Seq.addBlock(ePrephaseRO, ePrephasePE, ePrephase3D) ;


  % Update duration within TR
  TimeInTR = TimeInTR + Seq.blockDurations(end) ;
  % Update time from "excitation"
  TimeFromExcite = TimeFromExcite + Seq.blockDurations(end) ;

  % -----------------------------------------------------------------------
  % multi(?)-echo RO acquisition
  % -----------------------------------------------------------------------

  % Set phase of ADC to match RF spoiling phase. This way the RF spoiling
  % phase will not affect the acquired signal
  eADC.phaseOffset = eExciteRF.phaseOffset ;

  if (bInterleaveTrajectories)
    if (mod(TRCounter, 2) == 1)
      eLabelTraj = eLabelTraj1 ;
    else
      eLabelTraj = eLabelTraj2 ;
    end
  end

  % Store RO amplitude as we go in to echo loop (in case of a multi-echo
  % bi-polar acquisition, the sign will alternate).
  ROAmp0 = eRO.amplitude ;

  for TECounter = 1:NumTEs

    % Prepare for bi-polar or mono-polar RO grads in multi TE
    % -------------------------------------------------------

    % For TEs beyond the first we have to either switch the RO direction
    % alternatingly, or insert a rephaser before we can re-use our RO
    % gradient
    if (TECounter > 1)
      if (Actual.bBipolarROGrads)
        % We are in bi-polar multi TE mode, so switch sign of RO gradient
        eRO.amplitude = (-1)^(TECounter-1)*ROAmp0 ;
        % Switch between reversing ADC or not.
        if (mod(TECounter, 2)) % even
          eLabelReverse = eLabelRevOn ;
        else % Odd echo
          eLabelReverse = eLabelRevOff ;
        end

      else % monopolar case
  
        % We have to fully undo the last RO gradient before we can run the
        % next.
        Seq.addBlock(eROUndo) ;
  
        % Update duration within TR
        TimeInTR = TimeInTR + Seq.blockDurations(end) ;
        % Update time from "excitation"
        TimeFromExcite = TimeFromExcite + Seq.blockDurations(end) ;
      end
    end

    % Insert TE filler delay and the acqusition 
    % ------------------------------------------

    % Set filler delay to achieve requested TE (rounded up later)
    if (Actual.TE(TECounter) < 0)
      % Negative TE is interperted as using the minimal TE possible.
      TEFill = 0 ;
    else % try and achieve requested TE
      TEFill = Actual.TE(TECounter) - ...
               (TimeFromExcite + ADCDelay + ADCDuration/2) ;
    end

    % Sanity check
    if (TEFill < 0)
      error(['Cannot achieve desired TE[%d] = %f ms. ' ...
             'Minimum possible is %f ms.'], ...
            TECounter, 1e3*TE(TECounter), 1e3*(TE(TECounter) - TEFill))
    end

    % Round up to gradient raster
    TEFill = ceil(TEFill/GradRaster)*GradRaster ;

    % Set ADC labels (PE and partition). Note that the first index of each
    % label is zero (instead of marking k=0 position as zero index).
    % NOTE: We set the labels explicitly (using 'SET'), because the order
    %       may be arbitrary (depending on the order within PE3DOrder).
    eLabelPE = mr.makeLabel('SET', 'LIN', KPEIdx - KPEIdxMin) ; % PE
    eLabel3D = mr.makeLabel('SET', 'PAR', K3DIdx - K3DIdxMin) ; % 3D

    eLabelRefUse = eLabelNoRef ;
    eLabelRefAndImageUse = eLabelNoRefAndImage ;
    % TRCounter > 0 only when bADCOn, so we add that to our test
    if (bADCOn  && bRef(TRCounter)) % parallel imaging reference line
      % Reference line, so mark it as such.
      eLabelRefUse = eLabelRef ; 
      if (bImagingAndRef(TRCounter)) % Also a regular imaging line
        % reference and(!) imaging line, so mark it as such.
        eLabelRefAndImageUse = eLabelRefAndImage ;
      end
    end

    % add block:

    % Add delay to events (and remove it after adding to block)
    eRO.delay = eRO.delay + TEFill ;
    eADC.delay = eADC.delay + TEFill ;
    eROADCDelay.delay = eROADCDelay.delay + TEFill ;
    

    if (bADCOn)
      % ADC used
      Seq.addBlock(eRO, eADC, ...
                   eROADCDelay, ... ensure we are on the block raster
                   eLabelPE, eLabel3D, eLabelEchos(TECounter), ...
                   eLabelRefUse, eLabelRefAndImageUse, ...
                   eLabelTraj, ...
                   eLabelReverse) ; 
    else
       % no ADC
      Seq.addBlock(eRO, ...
                   eROADCDelay) ; ... ensure consistancy + the block raster
    end

    % remove extra delay (for next round)
    eRO.delay = eRO.delay - TEFill ;
    eADC.delay = eADC.delay - TEFill ;
    eROADCDelay.delay = eROADCDelay.delay - TEFill ;


    % update actual TE
    Actual.TE(TECounter) = TimeFromExcite + TEFill + ...
                          ADCDelay + ADCDuration/2 ;

    
    % Update duration within TR
    TimeInTR = TimeInTR + Seq.blockDurations(end) ;
    % Update time from "excitation"
    TimeFromExcite = TimeFromExcite + Seq.blockDurations(end) ;
  end
  
  % -----------------------------------------------------------------------
  % Rephasers
  % -----------------------------------------------------------------------

  % Update PE rephaser amplitude according to prephaser
  eRephasePE.amplitude =  -ePrephasePE.amplitude ; % [Hz/m]
  % Update 3D rephaser amplitude according to prephaser
  eRephase3D.amplitude =  -ePrephase3D.amplitude ; % [Hz/m]
  % Update RO rephaser amplitude according to prephaser
  if (Actual.bBipolarROGrads) % bi-polar RO gradients (alternating sign)
    eRephaseRO.amplitude =  (-1)^NumTEs * ePrephaseRO.amplitude ; % [Hz/m]
  else % mono-polar
    eRephaseRO.amplitude =  -ePrephaseRO.amplitude ; % [Hz/m]
  end

  % Add block of rephasers. 
  Seq.addBlock(eRephaseRO, eRephasePE, eRephase3D) ;

  % Update duration within TR
  TimeInTR = TimeInTR + Seq.blockDurations(end) ;

  % -----------------------------------------------------------------------
  % Spoilers
  % -----------------------------------------------------------------------

  if (bSpoilers)
    % Add block of spoilers. 
    Seq.addBlock(eSpoilRO, eSpoilPE, eSpoil3D) ;
  
    % Update duration within TR
    TimeInTR = TimeInTR + Seq.blockDurations(end) ;
  end

  % -----------------------------------------------------------------------
  % Update minimum possible TR (for information only)
  % -----------------------------------------------------------------------

  % update minimum possible TR (before adding the TR fill time)
  MinTRActual = max(MinTRActual, TimeInTR) ;

  % -----------------------------------------------------------------------
  % TR Fill block
  % -----------------------------------------------------------------------

  % Set filler delay to achieve requested TR (rounded up later)
  TRFill = Actual.TR - TimeInTR ;

 
  % Sanity check
  if (TRFill < -eps(0))
     error(['Total time (%f ms) of blocks within current TR (#%d) is ' ...
            'longer than desired TR (%f ms)!'], ...
           1e3*TimeInTR, TRCounter, 1e3*Actual.TR) ;
  end

  % update delay of eTRFill
  eTRFill.delay = TRFill ;
  % Add delay to the sequence
  Seq.addBlock(eTRFill) ;

  % Update duration within TR
  TimeInTR = TimeInTR + Seq.blockDurations(end) ;
  
end

% print TE used in last TR (hopefully the same as the rest)
fprintf(1, '  Actual TEs (in last TR): %f', 1e3*Actual.TE(1)) ;
  if (numel(Actual.TE) > 1)
    fprintf(1, ', %f', 1e3*Actual.TE(2:end))
  end
  fprintf(1, ' ms\n') ;
% Print minimum TR for current parameters.
fprintf(1, '  Min. possible TR for current sequence is: %f s.\n', MinTRActual) ;

% Mark end of filling sequence blocks
tEnd = toc(tStart) ;
fprintf(1, 'Done. (%g s)\n', tEnd) ;



%% Update [DEFINITIONS] of sequence including dummy comments (the setup)

% Mark start updating defintions
tStart = tic ;
fprintf(1, 'Updating [DEFINITIONS] ... ') ;


% Set definitions to be used by interpreter
clear Defs4Interpreter ; % clear from previous runs.
Defs4Interpreter.Name = SeqName ;
Defs4Interpreter.FOV = Actual.FOV ;
Defs4Interpreter.kSpaceCenterLine = KPECenterIdx - 1 ; % count from 0
Defs4Interpreter.kSpaceCenterPartition = K3DCenterIdx - 1 ; % count from 0
Defs4Interpreter.TE = Actual.TE ;
Defs4Interpreter.TR = Actual.TR ;

% Add Setup and Actual structs as comments in Seq (and then in .seq file),
% by exploiting the defintions mechanism.
% as well as add definition (Defs4Interpreter above) for use by
% interpreter, if it can read them.
Seq = Services.UpdateSeqDefs(Seq, Setup, Actual, Defs4Interpreter) ;


% Mark end of updating defintions
tEnd = toc(tStart) ;
fprintf(1, 'Done. (%g s)\n', tEnd) ;


%% Save sequence as a .seq file, if a file name/path is defined

if (exist('SeqName', 'var') && ~isempty(SeqName))

  % Mark start writing sequence to file
  tStart = tic ;
  fprintf(1, 'Writing .seq file ... ') ;

  if (exist('SeqFolder', 'var') )
    % ensure SeqFolder exists
    if (~isfolder(SeqFolder))
      mkdir(SeqFolder) ;
    end
    Seq.write(fullfile(SeqFolder, [SeqName, '.seq'])) ;
  else
    Seq.write([SeqName, '.seq']) ;
  end

  % Mark end of writing sequence to file
  tEnd = toc(tStart) ;
  fprintf(1, 'Done. (%g s)\n', tEnd) ;
  
end

%% check whether the timing of the sequence is correct

% Mark start checking sequence timing
tStart = tic ;
fprintf(1, 'Checking sequence timing ... ');

[ok, error_report]=Seq.checkTiming;

% Mark end of writing sequence to file
tEnd = toc(tStart) ;
fprintf(1, 'Done. (%g s)\n', tEnd) ;

if (ok)
  fprintf(1, 'Timing check passed successfully.\n');
else
  fprintf(1, '\n');
  fprintf(1, 'Timing check failed! Error listing follows:\n');
  fprintf(1, '\n');
  fprintf(1, [error_report{:}]);
  % When  the list is long the heading (that these are errors) is lost so
  % we repeat it at the end
  fprintf(1, '\n');
  fprintf(1, 'Timing check failed! Error listing above this line.\n');
end


%% plot sequence?

if (~bSupressPulseqPlots)

  % Mark start plotting
  tStart = tic ;
  fprintf(1, 'Plotting sequence ... ');
  
  % Plot waveforms/diagram of sequence
  % ----------------------------------
  
  % 1 sec before and after dummy scan end
  Seq.plot('timeRange', Actual.MinDurDummyScans + [-1,1]) ; 
% warning('Plot all!!!')
% Seq.plot() ;   
  
  % Plot k-space trajectory
  % -----------------------
  [ktraj_adc, t_adc, ktraj, t_ktraj, ...
   t_excitation, t_refocusing] = Seq.calculateKspacePP() ;
  figure; 
    plot(t_ktraj.', ktraj.') ;
    hold on ;
  
    xlabel('t [s]') ;
    ylabel('k [1/m]') ;
    title('k-space components as functions of time') ;
  
    % calculateKspacePP() should return physical x, y, z and not logical RO
    % PE and 3D.
    legend('k_x', 'k_y', 'k_z') ;
  
  
  % Plot k-space trajectory (3D)
  % ---------------------------
  figure; 
    plot3(ktraj(1,:), ktraj(2,:), ktraj(3,:),'b') ;
    hold on ; 
    plot3(ktraj_adc(1,:), ktraj_adc(2,:), ktraj_adc(3,:), 'r.') ; 
  
    % calculateKspacePP() should return physical x, y, z and not logical RO
    % PE and 3D.
    xlabel('k_x [1/m]') ;
    ylabel('k_y [1/m]') ;
    zlabel('k_z [1/m]') ;
    grid on ;
    title('3D k-space') ;
    legend('trajectory', 'ADC')
  
  % Mark end of plotting
  tEnd = toc(tStart) ;
  fprintf(1, 'Done. (%g s)\n', tEnd) ;

end




