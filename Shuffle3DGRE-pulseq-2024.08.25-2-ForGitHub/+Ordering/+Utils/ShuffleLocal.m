function [IdxPerPosOut, IOut, JOut] = ShuffleLocal(Length1, Length2, ...
                                                   ShuffleDiameter, ...
                                                   RandomSeed,  ...
                                                   bMeandering, ...
                                                   bCountFromZero, ...
                                                   bPlot, ...
                                                   SwapDistsWeightsIn)
% [IdxPerPosOut, IOut, JOut] = ShuffleLocal(Length1, Length2, ...
%                                           ShuffleDiameter, ...
%                                           RandomSeed,  ...
%                                           bCountFromZero, ...
%                                           bPlot, ...
%                                           WeightsIn)
% 
% Locally and randomly shuffle the order of elements in a matrix of size:
%   Length1 x Length2 
% Treats matrix as a 1D long ordered series. Attempts to shuffle the
% elements along the 1D vector so their shifts will obey a given
% distribution. The standard distribution is a symmetric triangle around a
% zero shift (most likely) with a "base" size ShuffleDiameter. Use
% SwapDistsWeightsIn to modify the distributions shape (see help below).
%  
% Inputs:
% -------
% - Length1 - Size of desired matrix along first dimension. (For a matrix,
%   this is the length of each column.)
% - Length2 - Size of desired matrix along second dimension. (For a matrix,
%   this is the length of each row.)
%   [Default: Length1]
% - ShuffleDiameter - The (approximate) support distance for the shuffling.
%   Beyond this (approximate) a sample cannot be moved (shuffled).
%   No shuffling is done if ShuffleDiameter < 2.
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
% - SwapDistsWeightsIn - An option to override the internally used
%   distribution of swaps performed. Used, only if non-empty.
%   If a scalar, instead of the default distribution (triangular) uses the
%   distribution to the power of 'SwapDistsWeightsIn'. The zero swap
%   distancce (no swap) is given a factor of one half after the power is
%   applied. (This is also applied when SwapDistsWeightsIn is not given, or
%   empty.)
%   Otherwise (and not empty), assumed a vector of length 
%   2*ShuffleDiameter + 1. This vector (after internal normaliztion, to a
%   sum of 1) is treated as the desired distribution of swaps from
%   -ShuffleDiameter to +ShuffleDiameter. (No extra "correction" is
%   performed internally, i.e., not division by 2 of probability for a zero
%   swap.)
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
      SwapDistsWeightsIn = [] ; 
    case 2
      ShuffleDiameter = Length1 * Length2 ;
      RandomSeed = 0 ;
      bMeandering = false ;
      bCountFromZero = false ; 
      bPlot = false ;
      SwapDistsWeightsIn = [] ; 
    case 3
      RandomSeed = 0 ;
      bMeandering = false ;
      bCountFromZero = false ; 
      bPlot = false ;
      SwapDistsWeightsIn = [] ; 
    case 4
      bMeandering = false ;
      bCountFromZero = false ; 
      bPlot = false ;
      SwapDistsWeightsIn = [] ; 
    case 5
      bCountFromZero = false ; 
      bPlot = false ;
      SwapDistsWeightsIn = [] ; 
    case 6
      bPlot = false ;
      SwapDistsWeightsIn = [] ; 
    case 7
      SwapDistsWeightsIn = [] ; 
    case 8
      % Do nothing
    otherwise
      error('Too many or too few input arguments.') ;
  end
  
  % For debugging
  bPlotAndStats = false ;
%   bPlotAndStats = true ;
  % Show plots every PlotStep swaps (if bPlotAndStats is true).
  PlotStep = 5000 ;

  % In case of a non-integer, round it.
  ShuffleDiameter = round(ShuffleDiameter) ;
  Length1 = round(Length1) ;
  Length2 = round(Length2) ;

  % total number of points in the curve.
  NIn = Length1 * Length2 ;
  
  % Initialize
  IdxPerPosOut = 1:NIn ;
  
  if (ShuffleDiameter < 1)
    % ShuffleDiameter = 0, means nothing to do. Anything smaller we will
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
  
  bIdxsUsed = false(size(IdxPerPosOut)) ;

  % Range of if indices around "current" center (chanegs all the time) on
  % which to consider swapping.
  IdxRangeAroundZero = -ShuffleDiameter:ShuffleDiameter ;
  % The set of allowed swap distances (typically the same as
  % IdxRangeAroundZero, but not necessarily).
  SwapDists = -ShuffleDiameter:ShuffleDiameter ;
  % Weights to give to each swapping distance. (Translated later to
  % probablilities.)
  if (isempty(SwapDistsWeightsIn))
    SwapDistsWeights = ShuffleDiameter - abs(SwapDists) ;
 
    % In general, for every sample moved +d in time, another is
    % moved -d in time. However, a zero distance swap affects only a single
    % sample. Thus, if we count how many shifts of each kind occurs:
    % Num Samples = 2*sum(#swaps != 0) + (#swaps == 0).
    % We want this to look triangular, so 
    SwapDistsWeights(SwapDists == 0) = SwapDistsWeights(SwapDists == 0) / 2 ;
  else
    if isscalar(SwapDistsWeightsIn)
      % Assume it is a power of the "default"
      SwapDistsWeights = ...
                 (ShuffleDiameter - abs(SwapDists)).^(SwapDistsWeightsIn) ;

      % In general, for every sample moved +d in time, another is
      % moved -d in time. However, a zero distance swap affects only a single
      % sample. Thus, if we count how many shifts of each kind occurs:
      % Num Samples = 2*sum(#swaps != 0) + (#swaps == 0).
      % We want this to look triangular, so 
      SwapDistsWeights(SwapDists == 0) = SwapDistsWeights(SwapDists == 0) / 2 ;
    else
      if (length(SwapDistsWeightsIn) ~= 2*ShuffleDiameter + 1)
        error('length of SwapDistsWeightsIn should be 2*ShuffleDiameter + 1') ;
      else
        SwapDistsWeights = Row(SwapDistsWeightsIn) ;
      end
      
    end
  end
% SwapDistsWeights = SwapDistsWeights.^2 ; 
% SwapDistsWeights = SwapDistsWeights.^3.0 ; 
% SwapDistsWeights = cos(0.5*pi/ShuffleDiameter * SwapDists) ;
% SwapDistsWeights = 1 - cos(pi/ShuffleDiameter * SwapDists) ;
% SwapDistsWeights = zeros(size(SwapDists)) + 2 ;
  
  % SwapDistsRemainCount - Count of swap distances we should generate.
  % Initialyy based on SwapDistsWeights and the number of smaples, this is
  % updated during the processing, so we know how many swaps of each
  % distance remain.
  % CDFSwapDistsInit - Initial cumulative distribution fucntion (CDF)
  % according to the desired probabilty density (from SwapDistsWeights).
  CDFSwapDistsInit = cumsum(SwapDistsWeights) / sum(SwapDistsWeights) ;
  % The actual SwapDistsRemainCount. Calculated so we will get intiger
  % values and cover all samples.
  SwapDistsRemainCount = diff([0, round(NIn*CDFSwapDistsInit)]) ;

  % A changing (evolving) probability for each swap distance to occur.
  % Based on SwapDistsRemainCount, i.e., the remaining (desired) number of
  % swaps of each distance. Initially based on SwapDistsWeights
  SwapDistsProbabilitiesUse = zeros(size(SwapDistsWeights)) ;
  
  % "Statistics" for debugging (and open a new figure)
  if (bPlotAndStats)

    % How many swaps of each distance have occured so far
    Stats_SwapDistsCounter = zeros(size(SwapDists)) ;
    % How many times, had to perform a swap of distance zero, i.e., could
    % not swap with any other sample, so did not move/swap.
    Stats_ForcedDistZero = 0 ;
    % Records the swap distance in the order they were performed. 
    % Initialized to a low impossible value, so can easily discard
    % unused/unfilled elements later.
    Stats_SwapDistsChronological = zeros(NIn, 1) - 2*ShuffleDiameter ;
    % Counter of how many swaps have so far been performed
    Stats_SwapsCounter = 0 ; 

    % open figure for progress, if desired.
    hDebug = figure ;
  end

  % loop in a random order of all samples. (probably not efficient)
  for CurrIdx = randperm(NIn)
    
    % Treat current sample only if it has not previously been swapped.
    if (~bIdxsUsed(CurrIdx))
      % Range of samples around current one to possibly swap with.
      Range = CurrIdx + IdxRangeAroundZero ;
      % Exclude points beyons the whole set of samples
      bInRange = ((Range > 0) & (Range <= NIn)) ;
      % Exclude samples (within range) previously swapped
      bInRange(bInRange) = ~(bIdxsUsed(Range(bInRange))) ;
      % Exclude swap distances of which we have passed our quota.
      % Zero distance swapping is always allowed, as a last resort.
      bInRange = bInRange & (SwapDistsRemainCount > 0) ;
      
      % random number used to decide which swap distance to choose. (See
      % below.)
      p = rand(1) ;

      % Choose swap distance, only if there is an option to something other
      % than zero.
      if (any(bInRange) && nnz(bInRange) > 1)
        % Update probability of possible swap distances. Based on currently
        % available distances and on the remaining quota of each swap
        % distance. (Doesn't work very well.)
        SwapDistsProbabilitiesUse(bInRange) = ...
                                     SwapDistsRemainCount((bInRange)) / ...
                                     sum(SwapDistsRemainCount(bInRange)) ;


        % find cumulative probability (CDF) for each swap distance
        pDist = cumsum(SwapDistsProbabilitiesUse(bInRange)) / ...
                sum(SwapDistsProbabilitiesUse(bInRange)) ;
        % Swap distance available, i.e., to samples not previously swapped.
        DistAvailable = SwapDists(bInRange) ;
        % Of the available swaps, which matches the random value generated.
        SwapDistSubIdx = find( pDist > p, 1, 'first') ;
        % Get actual swap distance.
        Dist = DistAvailable(SwapDistSubIdx) ;
            
      else % Only possible swap is no swap.
        Dist = 0 ;
        
        if (bPlotAndStats)
          Stats_ForcedDistZero = Stats_ForcedDistZero + 1 ;  
        end
      end

      % Perform swapping
      
      PosNew = CurrIdx + Dist ;

      IdxPerPosOut(PosNew) = CurrIdx ;
      IdxPerPosOut(CurrIdx) = PosNew ;
      
      % Update book keeping

      bIdxsUsed(PosNew) = true ;
      bIdxsUsed(CurrIdx) = true ;
      
      % update SwapDistsRemainCount, but do not allow negative values.
      DistIdx = find(SwapDists == Dist, 1) ;
      SwapDistsRemainCount(DistIdx) = ...
                                max(SwapDistsRemainCount(DistIdx) - 1, 0) ;



      % Plotting and "statistics" update for debug
      if (bPlotAndStats)

        Stats_SwapDistsCounter(DistIdx) = Stats_SwapDistsCounter(DistIdx) + 1 ;
        Stats_SwapsCounter = Stats_SwapsCounter + 1 ;
        Stats_SwapDistsChronological(Stats_SwapsCounter) = Dist ;

        figure(hDebug) ;

        if (mod(Stats_SwapsCounter, PlotStep) == 0)
          yyaxis left ;
            plot(SwapDists, Stats_SwapDistsCounter/sum(Stats_SwapDistsCounter)) ;
            hold on
            plot(SwapDists, SwapDistsWeights/sum(SwapDistsWeights), 'k-') ;
            hold off ;
          yyaxis right ;
            % % Distribution of swap distances.
            % plot(SwapDists, DistsRemain /sum(DistsRemain)) ;
          title(sprintf('SwapsCounter = %d, ForcedDistZero = %d', Stats_SwapsCounter, Stats_ForcedDistZero)) ;
        %   legend('p curr', 'p desired', 'p "future"') ;
          legend('p curr', 'p desired') ;
          pause(0.2) ;
        end
      end

      
    end
    
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

  if (bPlotAndStats)
    figure(hDebug) ;
    yyaxis left ;
      plot(SwapDists, Stats_SwapDistsCounter/sum(Stats_SwapDistsCounter)) ;
      hold on
      plot(SwapDists, SwapDistsWeights/sum(SwapDistsWeights), 'k-') ;
      hold off ;
    yyaxis right ;
      % % Distribution of swap distances.
      % plot(SwapDists, DistsRemain /sum(DistsRemain)) ;
    title(sprintf('SwapsCounter = %d, ForcedDistZero = %d', Stats_SwapsCounter, Stats_ForcedDistZero)) ;
  %   legend('p curr', 'p desired', 'p "future"') ;
    legend('p curr', 'p desired') ;
    pause(0.2) ;
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
