function [IdxPerPosOut, IOut, JOut] = Ordered(Length1, Length2, ...
                                              bMeandering, ...
                                              bElliptic, ...
                                              bCountFromZero, ...
                                              bPlot)
% [IdxPerPosOut, IOut, JOut] = Ordered(Length1, Length2, ...
%                                      bMeandering, ...
%                                      bElliptic
%                                      bCountFromZero, ...
%                                      bPlot)
% 
% A column-by-column filling of a matrix of size: Length1 x Length2 
%  
% Inputs:
% -------
% - Length1 - Size of desired matrix along first dimension. (For a matrix,
%   this is the length of each column.)
% - Length2 - Size of desired matrix along second dimension. (For a matrix,
%   this is the length of each row.)
%   [Default: Length1]
% - bMeandering - If true, will have the columns of the matrix alternating 
%   between going down and up. This way the position in the matrix will be
%   more continuous.
%   [Default: false]
%                 bMeandering = false        bMeandering = true
%                               
%                   1   5    9   13            1   8    9   16
%                   2   6   10   14    --->    2   7   10   15
%                   3   7   11   15            3   6   11   14
%                   4   8   12   16            4   5   12   13
% - bElliptic - If true, only an elliptical region in k-space/matrix is
%   covered. 
%   [Default: false]
% - bCountFromZero - Whether outputs count from zero (C++) or one (Matlab). 
%   Indices include the coordinates in final matrix (IOut, JOut) AND the
%   step counter within IdxPerPosOut matrix. 
%   [Default: false]
% - bPlot - If true, a line plot of the curve will be drawn. However, if no
%   output is requested (nargout < 1), then a plot will be shown.
%   [Default: false]
%
% Outputs:
% --------
% - IdxPerPosOut - A Length1 x Length2 matrix, with element holding the 
%   step in which the curve passes through the element (counts from 1).
% - IOut, JOut - the subscript order in which the curve passes through,
%   i.e., the coordinates of the curve. IOut gives the Length1 coordinates
%   and JOut is the Length2 coordnates (both counting from 1).


  % Defaults
  bMeanderingDefault = false ;
  bEllipticDefault = false ;
  bCountFromZeroDefault = false ;
  bPlotDefault = false ;

  switch nargin
    case 1
      Length2 = Length1 ;
      bElliptic = bEllipticDefault ;
      bMeandering = bMeanderingDefault ;
      bCountFromZero = bCountFromZeroDefault ; 
      bPlot = bPlotDefault ;
    case 2
      bElliptic = bEllipticDefault ;
      bMeandering = bMeanderingDefault ;
      bCountFromZero = bCountFromZeroDefault ; 
      bPlot = bPlotDefault ;
    case 3
      bMeandering = bMeanderingDefault ;
      bCountFromZero = bCountFromZeroDefault ; 
      bPlot = bPlotDefault ;
    case 4
      bCountFromZero = bCountFromZeroDefault ; 
      bPlot = bPlotDefault ;
    case 5
      bPlot = bPlotDefault ;
    case 6
      % Do nothing
    otherwise
      error('Too many or too few input arguments.') ;
  end

  % In case of a non-integer, round it.
  Length1 = round(Length1) ;
  Length2 = round(Length2) ;

  %% Set up eliiptical mask (or full k-space)
  
  if(bElliptic) % setup elliptical mask
    Mask = Ordering.Utils.EllipticMask(Length1, Length2) ;
  else
    Mask = true(Length1, Length2) ;
  end


  %% Fill IdxPerPosOut

  % total number of points in the curve.
  NumSamples = nnz(Mask) ;

  % Initialize IdxPerPosOut
  IdxPerPosOut = -1 * ones(Length1, Length2) ; % -1 marks unsampled

  
  if (bMeandering)
    % because the elliptic mask may not be symmetric to flipping up-down we
    % we need to flip it as well
    MaskFlipped = Mask ;
    MaskFlipped(:, 2:2:end) = flipud(Mask(:, 2:2:end)) ;

    IdxPerPosOut(MaskFlipped(:)) = ...
                           (1-bCountFromZero):(NumSamples-bCountFromZero) ;
    % Flip everything
    IdxPerPosOut(:, 2:2:end) = flipud(IdxPerPosOut(:, 2:2:end)) ;
  else
    IdxPerPosOut(Mask(:)) = (1-bCountFromZero):(NumSamples-bCountFromZero) ;
  end

  %% Determine IOut, JOut

  PositionIdxs = (1:(Length1*Length2)).' ;
  PositionIdxs = PositionIdxs(Mask(:)) ;
  [~, SortOrder] = sort(IdxPerPosOut(Mask(:))) ;
  PositionIdxs = PositionIdxs(SortOrder) ;
  [IOut, JOut] = ind2sub([Length1, Length2], PositionIdxs) ;
  if (bCountFromZero)
    IOut = IOut - 1 ;
    JOut = JOut - 1 ;
  end
  

  %% Plot?
  if (bPlot || nargout < 1)
    TitleStr = sprintf('%s (%dx%d)', mfilename, Length1, Length2) ;
    Ordering.Utils.PlotIdxPerPosAndTraj(IdxPerPosOut, IOut, JOut, ...
                                        bCountFromZero, TitleStr)
  end
  return ;

end
