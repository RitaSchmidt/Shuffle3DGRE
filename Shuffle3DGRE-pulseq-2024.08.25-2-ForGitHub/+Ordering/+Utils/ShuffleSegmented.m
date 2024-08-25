function [IdxPerPosOut, IOut, JOut] = ShuffleSegmented(Length1, Length2, ...
                                                       ShuffleDiameter, ...
                                                       RandomSeed,  ...
                                                       bMeandering, ...
                                                       bCountFromZero, ...
                                                       bPlot)
% [IdxPerPosOut, IOut, JOut] = ShuffleSegmented(Length1, Length2, ...
%                                               ShuffleDiameter, ...
%                                               RandomSeed,  ...
%                                               bCountFromZero, ...
%                                               bPlot)
% 
% Shuffle by segments the order of elements in a matrix of size:
%   Length1 x Length2 
% Treats matrix as a 1D long ordered series, splits it to similarly sized
% segments and randomly permutes the order within each such segment. The
% Segments are of approximate length of ~ShuffleDiameter.
% The lengths of the segment being shuffled are not exactly
% ShuffleDiameter, because the length may not be divisible by
% ShuffleDiameter, but tries to be approximately so.
%
% Inputs:
% -------
% - Length1 - Size of desired matrix along first dimension. (For a matrix,
%   this is the length of each column.)
% - Length2 - Size of desired matrix along second dimension. (For a matrix,
%   this is the length of each row.)
%   [Default: Length1]
% - ShuffleDiameter - The (approximate) length of segments in which the
%   order of steps is shuffled randomly.
%   No shuffling is done if ShuffleDiameter < 2. If ShuffleDiameter ==
%   Length1*Length2, function will shuffle all elemets randomly.
%   [Default: Length1 * Length2]
% - RandomSeed - A seed for the random generator. Passed to an rng()
%   command. For consistancy between runs, use the same seed always.
%   [Default: 0]
% - bMeandering - If true, pre-shuffled version (i.e., "ordered" version)
%   will have the columns of the matrix alternating between going down and
%   up. This way the position in the matrix will be more continuous after
%   shuffling.
%   [Default: false]
%                 bMeandering = false        bMeandering = true
%                               
%                   1   5    9   13            1   8    9   16
%                   2   6   10   14    --->    2   7   10   15
%                   3   7   11   15            3   6   11   14
%                   4   8   12   16            4   5   12   13
% - bCountFromZero - Whether indices count from zero or one. Indices
%   include the coordinates in final matrix (IOut, JOut) AND the step 
%   counter within IdxPerPosOut matrix.
%   [Default: false]
% - bPlot - If true, a line plot of the curve will be drawn. However, if no
%   output is requested (nargout < 1), then a plot will be shown.
%   [Default: false]
%
% Outputs:
% --------
% - IdxPerPosOut - A Length1 x Length2 matrix, with element holding the 
%   step in which the curve passes through the element (counts from zero or
%   one, depending on bCountFromZero).
% - IOut, JOut - the matrix subscripts in the order in which the curve
%   passes through, i.e., the coordinates of the curve. IOut gives the
%   Length1 coordinates and JOut is the Length2 coordnates (both counting
%   from zero or one, depending on bCountFromZero).
%
% NOTE: When counting from 1, the relation between IdxPerPosOut and IOut,
%       JOut is: 
%         [~, SortOrder] = sort(IdxPerPosOut(:)) ;
%         [IOut, JOut] = ind2sub(size(IdxPerPosOut), SortOrder) ;
%


  switch nargin
    case 1
      Length2 = Length1 ;
      ShuffleDiameter = Length1 * Length2 ;
      RandomSeed = 0 ;
      bMeandering = false ;
      bCountFromZero = false ; 
      bPlot = false ;
    case 2
      ShuffleDiameter = Length1 * Length2 ;
      RandomSeed = 0 ;
      bMeandering = false ;
      bCountFromZero = false ; 
      bPlot = false ;
    case 3
      RandomSeed = 0 ;
      bMeandering = false ;
      bCountFromZero = false ; 
      bPlot = false ;
    case 4
      bMeandering = false ;
      bCountFromZero = false ; 
      bPlot = false ;
    case 5
      bCountFromZero = false ; 
      bPlot = false ;
    case 6
      bPlot = false ;
    case 7
      % Do nothing
    otherwise
      error('Too many or too few input arguments.') ;
  end

  % In case of a non-integer, round it.
  ShuffleDiameter = round(ShuffleDiameter) ;
  Length1 = round(Length1) ;
  Length2 = round(Length2) ;

  % total number of points in the curve.
  NIn = Length1 * Length2 ;
  
  % Initialize
  IdxPerPosOut = 1:NIn ;
  
  if (ShuffleDiameter < 2)
    % ShuffleDiameter = 1, means nothing to do. Anything smaller we will
    % also disregard.
    if (bCountFromZero)
      IdxPerPosOut = reshape(0:(NIn-1), Length1, Length2) ;
      [IOut, JOut] = ndgrid(0:(Length1-1), 0:(Length2-1)) ;
    else
      IdxPerPosOut = reshape(1:NIn, Length1, Length2) ;
      [IOut, JOut] = ndgrid(1:Length1, 1:Length2) ;
    end

    % Handle "meandering" case
    if (bMeandering)
      IOut(:,2:2:end) = flipud(IOut(:,2:2:end)) ;
      IdxPerPosOut(:,2:2:end) = flipud(IdxPerPosOut(:,2:2:end)) ;
    end
    
    % Tranform into column vectors
    IOut = IOut(:) ;
    JOut = JOut(:) ;
    
    
    
    if (bPlot || nargout < 1)
      PlotResults(IdxPerPosOut, IOut, JOut, Length1, Length2, ...
                  bCountFromZero, ShuffleDiameter) ;
    end
    return ;
  end
  
  % Reset the seed for repeatability. (Should add as a possible input)
  rng(RandomSeed) ;
  
  % shuffle each of consecutive "blocks" of size ~ShuffleDiameter
  Idx = 0 ;
  NSegs1 = max(1, round(NIn/ShuffleDiameter)) ;
  for Counter1 = 1:NSegs1
    Remaining = NIn - Idx ;
    RangeLength = round(Remaining / (NSegs1 - Counter1 + 1)) ;
    Range = Idx + (1:RangeLength) ;
    IdxPerPosOut(Range) = Range(randperm(RangeLength)) ;
    Idx = Range(end) ;
  end
  
  
  % Reshape in to a matrix
  IdxPerPosOut = reshape(IdxPerPosOut, Length1, Length2) ;
  % Handle meandering case
  if (bMeandering)
    IdxPerPosOut(:,2:2:end) = flipud(IdxPerPosOut(:,2:2:end)) ;
  end
  
  % Determine IOut, JOut, if needed
  [~, SortOrder] = sort(IdxPerPosOut(:)) ;
  [IOut, JOut] = ind2sub([Length1, Length2], SortOrder) ;
  if (bCountFromZero)
    IOut = IOut - 1 ;
    JOut = JOut - 1 ;
  end
  
  % Plot results (if desired)
  if (bPlot || nargout < 1)
    PlotResults(IdxPerPosOut, IOut, JOut, Length1, Length2, ...
                bCountFromZero, ShuffleDiameter) ;
  end

end


function PlotResults(IdxPerPosOut, IOut, JOut, Length1, Length2, ...
                     bCountFromZero, ShuffleDiameter)

  IEdges = [1, Length1] ;
  JEdges = [1, Length2] ;
  if (bCountFromZero)
    IEdges = IEdges - 1 ;
    JEdges = JEdges - 1 ;
  end

  figure ;
    % Note: in imagesc(JEdges, IEdges, IdxPerPosOut), 'JEdges' describes
    % the x-axis, not the first dimension ...
    imagesc(JEdges, IEdges, IdxPerPosOut) ;
    hold on ;
    % Note: plot uses plot(x,y), but we want matrix dimension I,J, where
    % I is along the y-axis, so we use plot(JOut, IOut)
    plot(JOut, IOut, 'Color', [0.85, 0.33, 0.1]) ;
    title(sprintf('Shuffle Local %dx%d ("shuffle diameter" = %d)', ...
                  Length1, Length2, ShuffleDiameter)) ;
    xlabel('J (''Length2'')') ;
    ylabel('I (''Length1'')') ;
    axis ij

end
