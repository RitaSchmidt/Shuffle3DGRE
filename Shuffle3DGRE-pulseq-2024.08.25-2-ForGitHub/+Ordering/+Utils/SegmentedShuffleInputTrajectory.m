function [OrderMat,IOut,JOut,MaskOrderOut] = ...
                       SegmentedShuffleInputTrajectory(OrderMatIn, ...
                                                       TRtimesCutoffFreq, ...
                                                       bMask, RandomSeed, ...
                                                       bPlot)
% [OrderMat,IOut,JOut,MaskOrderOut] = ...
%                      SegmentedShuffleInputTrajectory(OrderMatIn, ...
%                                                      TRtimesCutoffFreq, ...
%                                                      bMask, RandomSeed, ...
%                                                      bPlot)
% Accepts an input matrix OrderMatIn whos values mark the order of sampling
% the matrix (values below 1 are assumed to not be sampled). An (optional)
% subsampling mask is first applied and then Segmented Shuffling of the 
% order is performed.
%
% Inputs:
% -------
% - OrderMatIn - a 2D matrix whos values are the order of sampling the
%   elements of the matrix. Values below(!) 1 mark elements not so be
%   sampled. Values are not required to be continuous, the order is
%   determined by sort().
% - TRtimesCutoffFreq - Product of the scanning TR (in seconds) and the
%   desired cutoff frequency. Here TR (repetition time) is the time between
%   samples. Possible values:
%   * [] (empty) - No shuffling. Keep original order.
%   * 0 - Random scrambling of the order.
%   * > 0 - Apply Segmented Shuffling according to the value.
%   [Default: [] (empty)]
% - bMask - An optional mask to apply to OrderMatIn. False elements in the
%   mask mark unsampled elements. If empty ([]), no additional masking is
%   done (OrderMatIn might already contain unsampled elements).
%   [Default: [] (empty)]
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
% - OrderMat - A 2D matrix of the same dimensions as OrderMatIn whos value
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
% - MaskOrderOut - Gives the reordering of bMask due to the shuffling:
%   bMask --> bMask(MaskOrderOut(:)). It accounts for the fact that not all
%   points may have been sampled in OrderMatIn, so MaskOrderOut may be
%   shorter(!) than samples inside bMask.
%   It is usefull if we have other masks or flags that are related to bMask
%   and we should know how to reorder them.
%
% NOTE: Segmented Shuffling currently treats all samples as if they were
%       sampled on an equi-spaced Cartesian grid. The samples are
%       trasnformed to a 1D array which is shuffled. 
%       Due to unsampled points in OrderMatIn (e.g., elliptic scanning) or
%       due in bMask (GRAPPA with only the center densely sampled) the
%       actual sampling may have regions which are more densely or less
%       densely sampled. This is not taken into account and a smarter
%       shuffling scheme is necessary.


%% Initial setup - set defaults for missing inputs

  % Set Some defaults
  TRtimesCutoffFreqDefault = [] ; % empty means no shuffling.
  bMaskDefault = [] ;
  RandomSeedDefault = 0 ;
  bPlotDefault = false ;

  % Value used to marked unsmapled element in OrderMat/OrderMatIn.
  % Should be smaller than any "sampled" index/value in OrderMatIn (indices
  % mark the order of sampling)
  NonSampledValue = -1 ; 

  % Handle partial input (use defaults instead)
  switch nargin
    case 1
      TRtimesCutoffFreq = TRtimesCutoffFreqDefault  ;
      bMask = bMaskDefault ;
      RandomSeed = RandomSeedDefault ;
      bPlot = bPlotDefault ;
    case 2
      bMask = bMaskDefault ;
      RandomSeed = RandomSeedDefault ;
      bPlot = bPlotDefault ;
    case 3
      RandomSeed = RandomSeedDefault ;
      bPlot = bPlotDefault ;
    case 4
      bPlot = bPlotDefault ;
    case 5
      % Do nothing
    otherwise
      error('Too many or too few input arguments.') ;
  end

  % For convenience
  [Length1, Length2] = size(OrderMatIn) ;

  
  % Apply subsampling given by bMask on OrderMatIn
  if (~isempty(bMask))
    OrderMatIn(~bMask(:)) = NonSampledValue ;
  end

  % get initial IOut and JOut (pre-shufling)
  % First find how to reorder OrderMatIn, so that values will be sorted
  [~, SortOrder] = sort(OrderMatIn(:));
  % When sorting, the first elemets will be those that are not sampled.
  NumNonSampled = nnz(OrderMatIn <= NonSampledValue) ;
  % Now we can define the non-shuffled(!) IOut and JOut
  [IOut, JOut] = ind2sub([Length1, Length2], ...
                         SortOrder((NumNonSampled+1):end)) ;

  % Determine ShuffleDiameter (or full scrambling)
  if (isempty(TRtimesCutoffFreq))
    IJReorder = 1:numel(IOut) ; % no reordering
    ShuffleDiameter = 1 ; % for PlotResults()
  elseif (TRtimesCutoffFreq == 0) % Peform full scrambling
    % Set the seed for repeatability.
    rng(RandomSeed) ;
    % Full ramdomization of order
    IJReorder = randperm(numel(IOut)) ;
    ShuffleDiameter = inf ; % for PlotResults()
  elseif(TRtimesCutoffFreq < 0)
    error('TRtimesCutoffFreq < 0. Can be empty ([]), 0, or positive') ;
  else
    % Internally, ShuffleSegmented round ShuffleDiameter, so we round it
    % here as well.
    ShuffleDiameter = round(1/TRtimesCutoffFreq) ;

    % Shuffle a 1D vector of the same length as (the new) IOut and Jout.
    [~, IJReorder] = ...
                Ordering.Utils.ShuffleSegmented(numel(IOut), 1, ...
                                            ShuffleDiameter, ...
                                            RandomSeed,  ...
                                            false, ... % not meandering
                                            false, ... % Count from 1
                                            false) ;   % no plotting
  end
  

  % Reorder IOut and JOut
  IOut = IOut(IJReorder(:)) ;
  JOut = JOut(IJReorder(:)) ;

  
  % Fill output OrderMat (after reordering and applying bMask)
  OrderMat = NonSampledValue * ones(Length1, Length2) ;
  OrderMat(sub2ind([Length1, Length2], IOut(:), JOut(:)))  = 1:numel(IOut) ;

  % Fill MaskOrderOut
  % MaskOrderOut is the reordring that bMask(!) underwent, but we should
  % account for the fact that not all points may have been sampled in
  % OrderMatIn, so MaskOrderOut may be shorter than samples inside bMask.
  if (~isempty(bMask))
    Vec1 = OrderMat(bMask(:)) ;
  else
    Vec1 = OrderMat(:) ;
  end
  % move to the end "samples" marked as not sampled (order index < 1).
  Vec1(Vec1(:) < 1) = Length1*Length2 + 1 ;
  % Get MaskOrderOut (inclding unsampled elements)
  [~, MaskOrderOut] = sort(Vec1) ;
  % Remove elements at the end that are not sampled
  MaskOrderOut = MaskOrderOut(1:numel(IOut)) ;


  % if bPlot is true, we to do it explicitly here.
  if (bPlot)
    bCountFromZero = false ;
    
    TitleStr = sprintf('Shuffle Segmented %dx%d ("shuffle diameter" = %g)', ...
                       Length1, Length2, ShuffleDiameter) ;
    Ordering.Utils.PlotIdxPerPosAndTraj(OrderMat, IOut, JOut, ...
                                        bCountFromZero, TitleStr)
  end


end