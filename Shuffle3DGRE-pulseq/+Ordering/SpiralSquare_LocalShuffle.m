function [OrderMat,IOut,JOut,MaskOrderOut] = ...
                             SpiralSquare_LocalShuffle(Length1, Length2, ...
                                                       TRtimesCutoffFreq, ...
                                                       bSampleMask, ...
                                                       ExtraParamsStruct, ...
                                                       RandomSeed, ...
                                                       bPlot)
% [OrderMat,IOut,JOut,MaskOrderOut] = ...
%                            SpiralSquare_LocalShuffle(Length1, Length2, ...
%                                                      TRtimesCutoffFreq, ...
%                                                      bSampleMask, ...
%                                                      ExtraParamsStruct, ...
%                                                      RandomSeed, ...
%                                                      bPlot)
% Generates a "square" spiral trajectory within a Cartesian matrix of
% dimensions Length1 x Length2. Optionally applies mask(s) on that
% trajectory and returned a locally shuffled version of the trajectory.
%
% Inputs:
% -------
% - Length1 - length of 1st dimension of the output OrderMat.
% - Length2 - length of 2nd dimension of the output OrderMat.
% - TRtimesCutoffFreq - Product of the scanning TR (in seconds) and the
%   desired cutoff frequency. Here TR (repetition time) is the time between
%   samples. Possible values:
%   * [] (empty) - No shuffling. Keep original order.
%   * 0 - Random scrambling of the order.
%   * > 0 - Apply Local Shuffling according to the value.
%   [Default: [] (empty)]
% - bSampleMask - An optional mask to apply to the base spiral trajectory. 
%   False elements in the mask mark unsampled elements. If empty ([]), no
%   additional masking is done.
%   [Default: [] (empty)]
% - ExtraParamsStruct - Structure containing extra controls of pre-shuffled
%   trajectory. Fields are:
%   * ExtraParamsStruct.bElliptic - If true an elliptic mask is applied on the
%     trajectory (corners of k-space are not sampled).
%     [Default: false]
%   * ExtraParamsStruct.bInsideOut - If true the trajectory starts at the
%     center of the matrix and spirals our. Otherwise, the trajectory is
%     reversed.
%     [Default: true]
% - RandomSeed - A seed for the random number generator (enables
%   repeatability of output). Value is passed to rng(), so any valid seed
%   of rng() is allowed.
%   [Default: 0]
% - bPlot - If true, displays an image of OrderMat with a plot of the
%   trajectory on top of it.
%   [Default: false]
%
% Outputs:
% --------
% - OrderMat - A 2D matrix of dimensions Length1 x Length2 whos values
%   give the new order of elements to sample (starting from 1). Unsampled
%   elemets are marked by -1 (controlled internally by NonSampledValue).
% - IOut - The ordered index along the columns(!) of the sampled elements
%   (ordered according to OrderMat). IOut, JOut together give the ordered
%   position inside OrderMat. (As in Matlab, indices start from 1, not
%   zero).
% - JOut - The ordered index along the rows(!) of the sampled elements
%   (ordered according to OrderMat). IOut, JOut together give the ordered
%   position inside OrderMat. (As in Matlab, indices start from 1, not
%   zero).
% - MaskOrderOut - Gives the reordering of bSampleMask due to the shuffling:
%   bSampleMask --> bSampleMask(MaskOrderOut(:)). It accounts for the fact
%   that not all points may have been sampled in OrderMatIn, so
%   MaskOrderOut may be shorter(!) than samples inside bSampleMask.
%   It is usefull if we have other masks or flags that are related to
%   bSampleMask and we should know how to reorder them.
%
% NOTE: Local Shuffling currently treats all samples as if they were
%       sampled on an equi-spaced Cartesian grid. The samples are
%       trasnformed to a 1D array which is shuffled. 
%       Due to unsampled points in OrderMatIn (e.g., elliptic scanning) or
%       due in bSampleMask (GRAPPA with only the center densely sampled)
%       the actual sampling may have regions which are more densely or less
%       densely sampled. This is not taken into account and a smarter
%       shuffling scheme is necessary.


%% Initial setup - set defaults for missing inputs

  % Set Some defaults
  TRtimesCutoffFreqDefault = [] ; % empty means no shuffling.
  ExtraParamsStructDefault.bElliptic = false ;
  ExtraParamsStructDefault.bInsideOut = true ;
  RandomSeedDefault = 0 ;
  bPlotDefault = false ;

  % Handle partial input (use defaults instead)
  switch nargin
    case 1
      Length2 = Length1 ;
      TRtimesCutoffFreq = TRtimesCutoffFreqDefault  ;
      bSampleMask = [] ;
      ExtraParamsStruct = ExtraParamsStructDefault ;
      RandomSeed = RandomSeedDefault ;
      bPlot = bPlotDefault ;
    case 2
      TRtimesCutoffFreq = TRtimesCutoffFreqDefault  ;
      bSampleMask = [] ;
      ExtraParamsStruct = ExtraParamsStructDefault ;
      RandomSeed = RandomSeedDefault ;
      bPlot = bPlotDefault ;
    case 3
      bSampleMask = [] ;
      ExtraParamsStruct = ExtraParamsStructDefault ;
      RandomSeed = RandomSeedDefault ;
      bPlot = bPlotDefault ;
    case 4
      ExtraParamsStruct = ExtraParamsStructDefault ;
      RandomSeed = RandomSeedDefault ;
      bPlot = bPlotDefault ;
    case 5
      RandomSeed = RandomSeedDefault ;
      bPlot = bPlotDefault ;
    case 6
      bPlot = bPlotDefault ;
    case 7
      % Do nothing
    otherwise
      error('Too many or too few input arguments.') ;
  end

  % Fill ExtraParamsStruct with defaults, if necesary
  ExtraParamNamesCell = fieldnames(ExtraParamsStructDefault) ;
  for FieldCounter = 1:numel(ExtraParamNamesCell)
    Field = ExtraParamNamesCell{FieldCounter} ;
    if(~isfield(ExtraParamsStruct, Field))
      ExtraParamsStruct.(Field) = ExtraParamsStructDefault.(Field) ;
    end
  end

  
  % Generate initial spiral on whole grid
  bCountFromZero = false ;
  bPlot0 = false ;
  [OrderMat0] = ...
            Ordering.BaseTraj.SpiralSquare(Length1, Length2, ...
                                           ExtraParamsStruct.bInsideOut, ...
                                           ExtraParamsStruct.bElliptic, ...
                                           bCountFromZero, bPlot0) ;

  % Shuffle the order
  bPlot = (bPlot || nargout < 1) ; % plot if no outputs
  [OrderMat,IOut,JOut,MaskOrderOut] = ...
            Ordering.Utils.LocalShuffleInputTrajectory(OrderMat0, ...
                                                       TRtimesCutoffFreq, ...
                                                       bSampleMask, ...
                                                       RandomSeed, ...
                                                       bPlot) ;

end
