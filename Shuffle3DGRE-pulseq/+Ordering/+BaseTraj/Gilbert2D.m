function [IdxPerPosOut, IOut, JOut] = Gilbert2D(Length1, Length2, ...
                                                bCountFromZero, bPlot)

% [IdxPerPosOut, IOut, JOut] = Gilbert2D(Length1, Length2, bPlot)
% 
% Generalized Hilbert ('gilbert') space-filling curve for arbitrary-sized
% 2D rectangular grids. Generates discrete 2D coordinates to fill a
% rectangle of size (width x height).
%
% A Matlab implementatin of the Python code from: 
% https://github.com/jakubcerveny/gilbert/blob/master/gilbert2d.py
%
% Inputs:
% -------
% - Length1 - Size of desired matrix along first dimension. (For a matrix,
%   this is the length of each column.)
% - Length2 - Size of desired matrix along second dimension. (For a matrix,
%   this is the length of each row.)
%   [Default: Length1]
% - bCountFromZero - Whether indices count from zero or one. Indices
%   include the coordinates in final matrix (IOut, JOut) and step counter
%   within IdxPerPosOut matrix.
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
%
% NOTE: The relation between IdxPerPosOut and IOut, JOut is:
%   [~, SortOrder] = sort(IdxPerPosOut(:)) ;
%   [IOut, JOut] = ind2sub(size(IdxPerPosOut), SortOrder) ; 

  switch nargin
    case 1
      Length2 = Length1 ;
      bCountFromZero = false ;
      bPlot = false ;
    case 2
      bCountFromZero = false ;
      bPlot = false ;
    case 3
      bPlot = false ;
    case 4
      % Do nothing.
    otherwise
      error('Too many or too few input arguments.') ;
  end

  bInit = true ; 
  
  if (Length2 >= Length1)
    [IdxPerPosOut, IOut, JOut] = generate2d(0, 0, Length2, 0, 0, Length1, ...
                                            bInit) ;
  else
    [IdxPerPosOut, IOut, JOut] = generate2d(0, 0, 0, Length1, Length2, 0, ...
                                            bInit) ;
  end
  
  % The results of generate2d start counting from 1 (Matlab style), if we
  % want to count from zero, have to remove 1 from results:
  if (bCountFromZero)
    IdxPerPosOut = IdxPerPosOut - 1 ;
    IOut = IOut - 1 ;
    JOut = JOut - 1 ;
  end
  
  % Plot results (if desired)
  if (bPlot || nargout < 1)
    TitleStr = sprintf('%s (%dx%d)', mfilename, Length1, Length2) ;
    Ordering.Utils.PlotIdxPerPosAndTraj(IdxPerPosOut, IOut, JOut, ...
                                        bCountFromZero, TitleStr)
  end
end

function [IdxPerPosOut, IOut, JOut] = generate2d(x, y, ax, ay, bx, by, bInit)

  persistent Reorder ;
  persistent Counter ;
  persistent I ;
  persistent J ;
  
  w = abs(ax + ay) ;
  h = abs(bx + by) ;
  
  if (bInit) % initialize static variables
    Reorder = zeros(ay+by, ax+bx) ;
    Counter = 0 ; 
    I = zeros((ay+by)*(ax+bx), 1) ;
    J = zeros((ay+by)*(ax+bx), 1) ;
      
  end
  
  dax = sign(ax) ;
  day = sign(ay) ;
  dbx = sign(bx) ;
  dby = sign(by) ;
  
  if (h == 1)
    % trivial row fill
    for ii = 0:(w-1)
      Counter = Counter + 1 ;
      Reorder(y+1, x+1) = Counter ;
      I(Counter) = y+1 ;
      J(Counter) = x+1 ;
      x = x + dax ;
      y = y + day ;
    end
    return ;
  end
  
  if (w == 1)
    % trivial column fill
    for ii = 0:(h-1)
      Counter = Counter + 1 ;
      Reorder(y+1, x+1) = Counter ;
      I(Counter) = y+1 ;
      J(Counter) = x+1 ;
      x = x + dbx ;
      y = y + dby ;
    end
    return ;
  end
  
  ax2 = floor(ax/2) ;
  ay2 = floor(ay/2) ;
  bx2 = floor(bx/2) ;
  by2 = floor(by/2) ;

  w2 = abs(ax2 + ay2) ;
  h2 = abs(bx2 + by2) ;
  
  if (2*w > 3*h)
      if (mod(w2, 2)&& (w > 2))
          % prefer even steps
          ax2 = ax2 + dax ;
          ay2 = ay2 + day ;
      end

      % long case: split in two parts only
      generate2d(x, y, ax2, ay2, bx, by, false) ;
      generate2d(x+ax2, y+ay2, ax-ax2, ay-ay2, bx, by, false) ;

  else
    if (mod(h2, 2) && (h > 2))
      % prefer even steps
      bx2 = bx2 + dbx ;
      by2 = by2 + dby ;
    end

    % standard case: one step up, one long horizontal, one step down
    generate2d(x, y, bx2, by2, ax2, ay2, false) ;
    generate2d(x+bx2, y+by2, ax, ay, bx-bx2, by-by2, false) ;
    generate2d(x+(ax-dax)+(bx2-dbx), y+(ay-day)+(by2-dby), ...      
               -bx2, -by2, -(ax-ax2), -(ay-ay2), false) ;
  end
  
  if (nargout > 0)
    IdxPerPosOut = Reorder ;
    IOut = I ;
    JOut = J ;
  end
  
end


