function [OrderMat,IOut,JOut,MaskOrderOut] = Gilbert(Length1, Length2, ...
                                                     ~, ... TRtimesCutoffFreq
                                                     bSampleMask, ...
                                                     ~, ... ExtraParamsStruct
                                                     ~, ... RandomSeed
                                                     bPlot)
% [OrderMat,IOut,JOut,MaskOrderOut] = Gilbert(Length1, Length2, ...
%                                             TRtimesCutoffFreq, ...
%                                             bSampleMask, ...
%                                             ExtraParamsStruct, ...
%                                             RandomSeed, ...
%                                             bPlot)
% Generates a Generalized Hilbert curve (a Gilbert curve) within a matrix
% of dimensions Length1 x Length2. Optionally applies mask on the matrix
% before(!) generating the trajectory. The mask must be such that a curve
% can be defined on it. Here this means that remaining elements can be
% compacted into a rectangular matrix. (Removing all blank rows and
% columns, leaves a matrix without any holes. Maybe a better generalization
% is also possible (or at least necessary).
%
% Inputs:
% -------
% - Length1 - length of 1st dimension of the output OrderMat.
% - Length2 - length of 2nd dimension of the output OrderMat.
% - TRtimesCutoffFreq - Dummy variable, not used here.
% - bSampleMask - An optional mask to apply before(!) generating the
%   Gilbert curve. The mask must generate a compact matrix that the Gilbert
%   curve can be defined on. (See "introduction" above.)
%   False elements in the mask mark unsampled elements. If empty ([]), no
%   additional masking is done.
%   [Default: [] (empty)]
% - ExtraParamsStruct - Dummy variable, not used here.
% - RandomSeed - Dummy variable, not used here.
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
% NOTE: The "Gilbert" trajectory used here is a Matlab implementatin of the
%       Python code from:  
%       https://github.com/jakubcerveny/gilbert/blob/master/gilbert2d.py




%% Initial setup - set defaults for missing inputs

  % Set Some defaults

  % TRtimesCutoffFreqDefault = [] ; % empty means no shuffling.
  % RandomSeedDefault = [] ;
  % ExtraParamsStructDefault = [] ;
  bPlotDefault = false ;

  % Handle partial input (use defaults instead)
  switch nargin
    case 1
      Length2 = Length1 ;
      % TRtimesCutoffFreq = TRtimesCutoffFreqDefault  ;
      bSampleMask = [] ;
      % ExtraParamsStruct = ExtraParamsStructDefault ;
      % RandomSeed = RandomSeedDefault ;
      bPlot = bPlotDefault ;
    case 2
      % TRtimesCutoffFreq = TRtimesCutoffFreqDefault  ;
      bSampleMask = [] ;
      % ExtraParamsStruct = ExtraParamsStructDefault ;
      % RandomSeed = RandomSeedDefault ;
      bPlot = bPlotDefault ;
    case 3
      bSampleMask = [] ;
      % ExtraParamsStruct = ExtraParamsStructDefault ;
      % RandomSeed = RandomSeedDefault ;
      bPlot = bPlotDefault ;
    case 4
      % ExtraParamsStruct = ExtraParamsStructDefault ;
      % RandomSeed = RandomSeedDefault ;
      bPlot = bPlotDefault ;
    case 5
      % RandomSeed = RandomSeedDefault ;
      bPlot = bPlotDefault ;
    case 6
      bPlot = bPlotDefault ;
    case 7
      % Do nothing
    otherwise
      error('Too many or too few input arguments.') ;
  end


  % We assume here that we do not count from zero. 
  % (It is probably not enough to change this to true if you do want to
  %  count from zero.) 
  bCountFromZero = false ; % True requires more changes!!

  % Generate trajectory/outputs
  % ---------------------------

  % We handle the cases of a sample mask or none separately

  if (~isempty(bSampleMask)) % A sampling mask has been given

    % Is bSampleMask of the exected size (Length1 x Length2)
    if (~ismatrix(bSampleMask) || ...
        size(bSampleMask, 1) ~= Length1 || ...
        size(bSampleMask, 2) ~= Length2)
      error('bSampleMask dimensions are not Length1 x Length2 (not %dx%d)', ...
            Length1, Length2) ;
    end
  
    % Can bSampleMask be "compacted" to a smaller "fully sampled" matrix?
    % Find unsampled rows
    bUnsampledRows = all(~bSampleMask, 2) ;
    % Find unsampled columns
    bUnsampledCols = all(~bSampleMask, 1) ;
    % Remove unsampled rows and columns
    bCompactMask = bSampleMask(~bUnsampledRows, :) ;
    bCompactMask = bCompactMask(:, ~bUnsampledCols) ;
    % Are there still unsampled elements?
    if (any(~bCompactMask))
      error(['bSampleMask cannot be compacted into a full matrix by ' ...
             'discarding unsampled rows and unsampled columns']) ;
    end

    % Get size of compact sampled matrix
    [LCompact1, LCompact2] = size(bCompactMask) ;

    % Generate a Gilbert (generalized Hilbert curve) trajectory inside the
    % compat matrix
    bPlot0 = false ;
    [OrderMatCompact] = Ordering.BaseTraj.Gilbert2D(LCompact1, LCompact2, ...
                                                    bCountFromZero, bPlot0) ;

    % Move from compat matrix to full matrix

    % Initialize
    NonSampledValue = -1 ;
    OrderMat = NonSampledValue * ones(Length1, Length2) ;
    % Fill OrderMat
    OrderMat(bSampleMask(:)) = OrderMatCompact(:) ;

    % Get IOut and JOut

    % Give an runing index to each element in the sample (will be used to
    % generate subscripts I,J of each sampled element).
    PositionIdxs = reshape(1:(Length1*Length2), Length1, Length2) ;
    % Get subset of PositionIdxs, should be the same size as
    % OrderMatCompact 
    PositionIdxs = PositionIdxs(bSampleMask) ; 
    % Find order of elements inside the matrix
    [~, SortOrder] = sort(OrderMatCompact(:)) ;
    % Use the order within the matrix to re-order PositionIdxs
    PositionIdxs = PositionIdxs(SortOrder(:)) ;
    % Get the I,J subscripts of the trajectory from the re-ordered
    % PositionIdxs.
    [IOut, JOut] = ind2sub([Length1, Length2], PositionIdxs) ;

    % Set MaskOrderOut
    MaskOrderOut = SortOrder(:) ;

  else % No sampling mask, fill the whole matrix.

    bPlot0 = false ;
    [OrderMat, ...
     IOut, JOut] = Ordering.BaseTraj.Gilbert2D(Length1, Length2, ...
                                               bCountFromZero, bPlot0) ;
    MaskOrderOut = [] ;
  end

  % Plot results, if requested
  % --------------------------

  if (bPlot || nargout < 1)
    figure ;
      % Note: in imagesc(JEdges, IEdges, IdxPerPosOut), 'JEdges' describes
      % the x-axis, not the first dimension ...
%       IEdges = [1, Length1] ;
%       JEdges = [1, Length2] ;
%       imagesc(JEdges, IEdges, IdxPerPosOut) ;
      imagesc(OrderMat) ;
      hold on ;
      % Note: plot uses plot(x,y), but we want matrix dimension I,J, where
      % I is along the y-axis, so we use plot(JOut, IOut)
      plot(JOut, IOut, 'Color', [0.85, 0.33, 0.1]) ;
      title(sprintf('Gilbert %dx%d', Length1, Length2)) ;
      xlabel('J (''Length2'')') ;
      ylabel('I (''Length1'')') ;
      axis ij

  end


end
