function [ReorderOut] = PeanoMeanderingRandmizedReorder(nOrder, pV1, pH1, ...
                                                        Seed, ...
                                                        bPlot, MarkerSize, ...
                                                        ColorOrder)

% [ReorderOut] = PeanoRandmizedReorder(nOrder, pV1, pH1, Seed, bPlot, ...
%                                      MarkerSize, ColorOrder)
%
% Generates a 2D space filling curve, based on the Peano curve.
% Fills a grid of points of size 3^nOrder x 3^nOrder, returns its order 
% and optionally plots it.
%
% Inputs
% ------
% 
% - nOrder - Determines size of ooutput matrix: 3^nOrder x 3^nOrder
% - pV1 - probability of using a "vertical" building block, of "type 1"
%   instead of "vertical" building block, of "type 2".
%   [Default: 1]
% - pH1 - probability of using a "horizontal" building block, of "type 1"
%   instead of "horizontal" building block, of "type 2".
%   [Default: 1]
% - Seed - Input for rng() command called. Used at start to set the seed.
%   [Default: 0]
% - bPlot - Whether to plot the curve or not. The color of the grid points
%   gradually chages to help visualize the curve's path. (But see also
%   ColorOrder.)
%   [Default: false]
% - MarkerSize - Size of marker used in plot (to mark the grid points).
% - ColorOrder - To help visualize the path of the curve, the color of the 
%   marker (marking the grid point) gradually changes. To further help, the
%   user can control when the color is reset. If ColorOrder is set to 
%   nOrder, the color changes gradually from start to end. If, however, it
%   is set to nOrder-1, the color changs gradually within each sub-block 
%   (3x3 of those) of the whole data, and so on. ColorOrder, means that all
%   markers have the same color. (Non-integer numbers will produce strange
%   results.)
%
% Outputs
% -------
%
% - ReorderOut - A 3^nOrder x 3^nOrder matrix, such that each element
%   stores its own order within the curve. Starting with the element
%   holding zero (0), then the element holding one (1) and so on.

  switch nargin
    case 1
      pV1 = 1 ;
      pH1 = 1 ;
      Seed = 0 ;
      bPlot = false ;
      MarkerSize = 12 ;
      ColorOrder = [] ;
    case 2
      pH1 = 1 ;
      Seed = 0 ;
      bPlot = false ;
      MarkerSize = 12 ;
      ColorOrder = [] ;
    case 3
      Seed = 0 ;
      bPlot = false ;
      MarkerSize = 12 ;
      ColorOrder = [] ;
    case 4
      bPlot = false ;
      MarkerSize = 12 ;
      ColorOrder = [] ;
    case 5
      MarkerSize = 12 ;
      ColorOrder = [] ;
    case 6
      ColorOrder = [] ;
    case 7
      % do nothing ; (we have all inputs)
    otherwise
      error('Too many or too few input arguments.') ;
  end

  % Set random generator seed (for repeatability)
  if isempty(Seed)
    Seed = 0 ;
  end
  rng(Seed) ;

  % Generated recursively
  bVertical = true ;
  bForward = true ; 
  ReorderOut = PeanoRecurse(nOrder, pV1, pH1, bVertical, bForward) ;
  
  % plot?
  if (bPlot)
    [X, Y] = meshgrid(1:(3^nOrder), 1:(3^nOrder)) ;
    [~, SortedOrder] = sort(ReorderOut(:)) ;
    X = X(SortedOrder(:)) ;
    Y = Y(SortedOrder(:)) ;
    Length = numel(X) ;
    if (isempty(ColorOrder))
      ColorIdx = mod(linspace(0, Length-1, Length), Length) ;
    else
      ColorIdx = mod(linspace(0, Length-1, Length),9^ColorOrder) ;
    end
    ColorIdx = ColorIdx(:) ;
    
    figure ; 
      plot(X, Y, '-')
      hold on
      scatter(X, Y, MarkerSize, ColorIdx, 'filled') ;
      axis ij ;
      axis tight ;
  end
  
end

function [ReorderOut] = PeanoRecurse(nOrder, pV1, pH1, bVertical, bForward)

  % nOrder = 0 is a trivial case of a single element with counter zero
  % along the peano curve.
  if (nOrder < 1)
    ReorderOut = 0 ;
    return ;
  end

  p = pV1 ;
  q = pH1 ;
  n = nOrder - 1 ;

  if (bVertical)
    bV1 = (rand()< pV1) ;
    
    if (bV1)
      M11 = PeanoRecurse(n, p, q,  bVertical,  bForward) ;
      M21 = PeanoRecurse(n, p, q, ~bVertical,  bForward) ;
      M31 = PeanoRecurse(n, p, q, ~bVertical, ~bForward) ;
      M32 = PeanoRecurse(n, p, q, ~bVertical, ~bForward) ;
      M22 = PeanoRecurse(n, p, q, ~bVertical,  bForward) ;
      M12 = PeanoRecurse(n, p, q,  bVertical, ~bForward) ;
      M13 = PeanoRecurse(n, p, q,  bVertical,  bForward) ;
      M23 = PeanoRecurse(n, p, q,  bVertical,  bForward) ;
      M33 = PeanoRecurse(n, p, q,  bVertical,  bForward) ;
      d = 9^n * [0 3 4
                 1 2 5
                 8 7 6] ;
    else
      M11 = PeanoRecurse(n, p, q, ~bVertical,  bForward) ;
      M21 = PeanoRecurse(n, p, q, ~bVertical, ~bForward) ;
      M31 = PeanoRecurse(n, p, q,  bVertical,  bForward) ;
      M32 = PeanoRecurse(n, p, q,  bVertical, ~bForward) ;
      M22 = PeanoRecurse(n, p, q, ~bVertical, ~bForward) ;
      M12 = PeanoRecurse(n, p, q, ~bVertical,  bForward) ;
      M13 = PeanoRecurse(n, p, q,  bVertical,  bForward) ;
      M23 = PeanoRecurse(n, p, q,  bVertical,  bForward) ;
      M33 = PeanoRecurse(n, p, q,  bVertical,  bForward) ;
      d = 9^n * [0 1 2
                 7 6 3
                 8 5 4] ;
    end
  else
    bH1 = (rand()< pH1) ;
    
    if (bH1)
      M11 = PeanoRecurse(n, p, q,  bVertical,  bForward) ;
      M21 = PeanoRecurse(n, p, q,  bVertical, ~bForward) ;
      M31 = PeanoRecurse(n, p, q,  bVertical,  bForward) ;
      M32 = PeanoRecurse(n, p, q,  bVertical,  bForward) ;
      M22 = PeanoRecurse(n, p, q, ~bVertical,  bForward) ;
      M12 = PeanoRecurse(n, p, q, ~bVertical,  bForward) ;
      M13 = PeanoRecurse(n, p, q, ~bVertical, ~bForward) ;
      M23 = PeanoRecurse(n, p, q, ~bVertical, ~bForward) ;
      M33 = PeanoRecurse(n, p, q,  bVertical,  bForward) ;
      
      d = 9^n * [0 1 8
                 3 2 7
                 4 5 6] ;
    else
      M11 = PeanoRecurse(n, p, q, ~bVertical,  bForward) ;
      M21 = PeanoRecurse(n, p, q, ~bVertical,  bForward) ;
      M31 = PeanoRecurse(n, p, q,  bVertical,  bForward) ;
      M32 = PeanoRecurse(n, p, q,  bVertical,  bForward) ;
      M22 = PeanoRecurse(n, p, q, ~bVertical, ~bForward) ;
      M12 = PeanoRecurse(n, p, q, ~bVertical, ~bForward) ;
      M13 = PeanoRecurse(n, p, q,  bVertical,  bForward) ;
      M23 = PeanoRecurse(n, p, q,  bVertical, ~bForward) ;
      M33 = PeanoRecurse(n, p, q,  bVertical,  bForward) ;
      
      d = 9^n * [0 7 8
                 1 6 5
                 2 3 4] ;
    end
  end
  
  
  if (bForward)
    ReorderOut = [d(1,1) + M11, d(1,2) + M12, d(1,3) + M13
                  d(2,1) + M21, d(2,2) + M22, d(2,3) + M23
                  d(3,1) + M31, d(3,2) + M32, d(3,3) + M33] ;
  else
    ReorderOut = [d(3,3) + M33, d(3,2) + M32, d(3,1) + M31
                  d(2,3) + M23, d(2,2) + M22, d(2,1) + M21
                  d(1,3) + M13, d(1,2) + M12, d(1,1) + M11] ;
  end
  
  return ;
end