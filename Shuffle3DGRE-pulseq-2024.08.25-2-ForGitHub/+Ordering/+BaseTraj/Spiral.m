function [IdxPerPosOut, IOut, JOut] = Spiral(Length1, Length2, ...
                                             bInsideOut, ...
                                             bElliptic, ...
                                             bEllipticByRadius, ...
                                             bCountFromZero, ...
                                             bPlot)
                                                    
% [IdxPerPosOut, IOut, JOut] = Spiral(Length1, Length2, ...
%                                     bInsideOut, ...
%                                     bElliptic, bEllipticByRadius, ...
%                                     bCountFromZero, bPlot)
% 
% Shuffle the order of elements in a matrix to mimic as an inside-out
% "circular" spiral.
%
% Inputs:
% -------
% - Length1 - Size of desired matrix along first dimension. (For a matrix,
%   this is the length of each column.)
% - Length2 - Size of desired matrix along second dimension. (For a matrix,
%   this is the length of each row.)
%   [Default: Length1]
% - bInsideOut - If true starts from center of matrix and wind out.
%   Otherwise, starts outside and winds inward.
%   If true, will not (necessarily) start from the center, but a bit off,
%   so can cover the whole matrix by unwinding. (Just as in the reverese
%   order will not end at the center.)
% - bElliptic - If true, only an elliptical region in k-space is covered. 
%   (If all of k-space is covered, large jumps will be required when the
%   curve reaches beyond the k-space rectangle, to the next point when it
%   comes back in. 
%   [Default: true]
% - bEllipticByRadius - Only relevant when bElliptic is true. If
%   bEllipticByRadius is false, the spiral stops the first time it goes
%   out of the matrix. Otherwise, when true, "continues" out of the matrix
%   and returns back in, until the threshold radius is reached. The true
%   case can have more samples and cover the edges of k-space, however, it 
%   may include large jumps between going outside the matrix and back in.
%   [Default: false]
% - bEllipticByRadius - If true (and bElliptic is true), the mask used will 
%   be defined by a limiting "radius". (If the matrix is not square, an
%   effective radius is used, scaling one dimension.)
%   If false, the mask is based on the spiral generated (as long as it does
%   not step outside of the rectangular k-space. This avoids curve jumps
%   larger than sqrt(2). Currently(?) number samples in this case is the
%   minimum between the radius limited number and the spiral limited
%   number. (Minimal testng showed that the radius based mask limits the
%   number.)
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
  bInsideOutDefault = true ;
  bEllipticDefault = true ;
  bEllipticByRadiusDefault = false ;
  bCountFromZeroDefault = false ;
  bPlotDefault = false ;

  switch nargin
    case 1
      Length2 = Length1 ;
      bInsideOut = bInsideOutDefault ;
      bElliptic = bEllipticDefault ;
      bEllipticByRadius = bEllipticByRadiusDefault ;
      bCountFromZero = bCountFromZeroDefault ;
      bPlot = bPlotDefault ;
    case 2
      bInsideOut = bInsideOutDefault ;
      bElliptic = bEllipticDefault ;
      bEllipticByRadius = bEllipticByRadiusDefault ;
      bCountFromZero = bCountFromZeroDefault ;
      bPlot = bPlotDefault ;
    case 3
      bElliptic = bEllipticDefault ;
      bEllipticByRadius = bEllipticByRadiusDefault ;
      bCountFromZero = bCountFromZeroDefault ;
      bPlot = bPlotDefault ;
    case 4
      bEllipticByRadius = bEllipticByRadiusDefault ;
      bCountFromZero = bCountFromZeroDefault ;
      bPlot = bPlotDefault ;
    case 5
      bCountFromZero = bCountFromZeroDefault ;
      bPlot = bPlotDefault ;
    case 6
      bPlot = bPlotDefault ;
    case 7
      % Do nothing
    otherwise
      error('Too many or too few input arguments.') ;
  end
  
  % Ensure bCountFromZero is zero or one, only
  bCountFromZero = logical(bCountFromZero) ;
  
  % Round integer values
  Length1 = round(Length1) ;
  Length2 = round(Length2) ;
  LengthOrig1 = Length1 ;
  LengthOrig2 = Length2 ;
  
  % Add "buffer" rows and columns (for use with masking when spiral goes
  % out of matrix for the first/last time).
  if (bElliptic && ~bEllipticByRadius) % "Cheat"
    Length1 = LengthOrig1 + 2 ;
    Length2 = LengthOrig2 + 2 ;
    % Because we expand the first I, J we actually use will be 2, so to
    % shift the final result back to 1, need to remove 1 (done at the end).
    ExpansionFix1 = -1 ;
    ExpansionFix2 = -1 ;
  else
    ExpansionFix1 = 0 ;
    ExpansionFix2 = 0 ;
  end

  % Total number of voxels in the matrix.
  NumVoxelsTot = Length1*Length2 ;

  
  % get k-space "indices" centered around zero.
  [ii, jj]=ndgrid((1:Length1) - Ordering.Utils.FFTCenterIndex(Length1), ...
                  (1:Length2) - Ordering.Utils.FFTCenterIndex(Length2)) ;

  % Scaling so we can handle rectangular matrices. Similar to 'a', 'b' in
  % defintions of ellipses, only their reciprocals (e.g., a --> 1/a), I
  % think. 
  a = LengthOrig1/min(LengthOrig1, LengthOrig2) ;
  b = LengthOrig2/min(LengthOrig1, LengthOrig2) ;

  % Find angle and distance to each voxel from center
  theta0 = atan2(  b *  ii, a * jj) ;
  r = sqrt(ii.^2 + jj.^2) ;

  % We assume a spiral so that:
  %   ii = a*theta/(2*pi)*sin(theta)
  %   jj = b*theta/(2*pi)*cos(theta)
  % We found theta0 for all points, but it is constrained to -pi to pi. We
  % now find n, so that out theta will be theta0 + 2*pi * n
  n = r ./ ...
         sqrt(a^2 * sin(theta0).^2 + b^2 * cos(theta0).^2) - theta0/(2*pi) ;
  % round n (because it must be an integer)
  n = round(n) ;
  
  % Generate the "spiraling" theta
  theta = theta0 + 2*pi*n ;

  
  if(bElliptic)
    
    % Generate a mask from an effective radius. 

    % We need this mask, even if we do not really use it, to find the
    % number of points it includes. (We use this number to possibly limit
    % the number of points in the other masking method.)

    % Effective radius of all points ("elliptical"), acounting for the
    % rectangular shape of the matrix.
    rEff = sqrt(b^2*ii.^2 + a^2*jj.^2) ;

    % Actually Define the mask
    % We would like to use floor below (for odd lengths), but slight
    % rounding errors might throw us off. Since the results should be
    % either integer or half integer, subtractinbg 0.1 and rounding, should
    % work as well. (Assuming rounding errors are much smaller than 0.1)
    % Mask = (rEff <= floor(max(Length1/2, Length2/2))) ;
    RLim = round(max(LengthOrig1/2, LengthOrig2/2) -0.1) ;
    Mask = (rEff <= RLim + eps(RLim)) ;
    if (LengthOrig1 ~= Length1)
      Mask(1,:) = false ;
      Mask(end,:) = false ;
      Mask(:,1) = false ;
      Mask(:, end) = false ;
    end



    % Give points outside of mask an extreme enough value, so they will be
    % at the end after sorting (ascending or descending). Needed only for
    % radius based mask. (The other methods just keeps the first/last set
    % of point after sorting.)
    if (bEllipticByRadius)
      if (bInsideOut)
        % Will sort ascending, so points outside should be "too" large (and
        % appear at the end)
        theta(~Mask) = 2*pi*(Length1+Length2) ;
      else
        % Will sort descending, so points outside should be "too" small
        % (smaller than -pi) and appear at the end.
        theta(~Mask) = -4 ;
      end
    end
  end
  

  % Sort points according to the spiralling theta (inwards or outwards).
  if (bInsideOut)
    [~, SortOrder] = sort(theta(:), 'ascend') ;
  else
    [~, SortOrder] = sort(theta(:), 'descend') ;
  end
  
  % Find order of voxels ("coordinates" IOut, JOut), includes points which
  % may be masked later.
  [IOut, JOut] = ind2sub([Length1, Length2], SortOrder(:)) ;

  
  % Determine which voxels to keep after masking.
  % Depends on masking method and spiraling direction (in or out). 
  if (bElliptic)
    % number of points in radius-based mask (used to limit the number of
    % points even when the radius mask is not used explicitly used ).
    NumVoxelsInMask = nnz(Mask) ;
    
    if (bEllipticByRadius)
      % Sorting already put the extra (masked) points at the end, so just
      % use the first NumVoxelsInMask points.
      SortRangeUse = 1:NumVoxelsInMask ;
    else
      % External edge of spiral is at most where it goes out of matrix
      % rectangle (could be before, if there are "too many" points).
      % If going inside-out, it is the first point the next step is greater
      % than one voxel (straight or diagonal), i.e. > sqrt(2). Since the
      % next largest step is 2, then to avoid rounding errors, we just
      % check comapred to 1.5 (instead of sqrt(2)). If going in the
      % opposite direction (outside-in) it is the last(!) time a step is
      % larger than allowed + 1 (need the next point).
      
      MinIOut = 2 ;
      MaxIOut = Length1 - 1  ;
      MinJOut = 2 ;
      MaxJOut = Length2 - 1 ;
      
      % Determine start and end indices along the curve.
      if (bInsideOut)
        % We start at the center, which is always included, and wind out
        % to the edge.
        SortStartIdx = 1 ;
        SortEndIdx = min(find(IOut(:) < MinIOut | ...
                              JOut(:) < MinJOut | ...
                              IOut(:) > MaxIOut | ...
                              JOut(:) > MaxJOut , ...
                              1, 'first' ), ...
                         NumVoxelsInMask) ;
      else
        % We start at the edge and wind in to the center, which is always
        % included
        SortStartIdx = max(find(IOut(:) < MinIOut | ...
                                JOut(:) < MinJOut | ...
                                IOut(:) > MaxIOut | ...
                                JOut(:) > MaxJOut , ...
                                1, 'last' ), ...
                           NumVoxelsTot - NumVoxelsInMask + 1) ;
        SortEndIdx = NumVoxelsTot ;
        
      end

      % The Range of points to use, once we found SortStartIdx, SortEndIdx.
      SortRangeUse = SortStartIdx:SortEndIdx ;

    end
  else % No mask (use everything).
    NumVoxelsIn = LengthOrig1 * LengthOrig2 ;
    SortRangeUse = (NumVoxelsTot - NumVoxelsIn + 1):NumVoxelsTot ;
  end

  % double version of binary variables, for numerics ('Colon operands
  % cannot be logical')
  bdCountFromZero = double(bCountFromZero) ;
  bdCountFromOne  = double(~bCountFromZero) ;

  % Limit indices to the first or last NumVoxelsIn (the rest are outside
  % our "mask"). + Fix for extra + 2 along columns and rows
  IOut = IOut(SortRangeUse) - bdCountFromZero + ExpansionFix1 ;
  JOut = JOut(SortRangeUse) - bdCountFromZero + ExpansionFix2 ;

  % Fill IdxPerPosOut, with the step number at each position
  
  % initialize to value outside of our mask
  IdxPerPosOut = zeros(LengthOrig1, LengthOrig2) -1 ;
  
  % fill (depends on whether we count from 1 or zero)
  NumVoxelsIn = SortRangeUse(end) - (SortRangeUse(1)) + 1 ;
  % SortOrder might not be adequate for the size of the final output, so we
  % fill based on IOut & JOut, but need to account for bdCountFromZero.
  FillOrder = sub2ind([LengthOrig1, LengthOrig2], ...
                      IOut(:) + bdCountFromZero, ...
                      JOut(:) + bdCountFromZero) ;
  IdxPerPosOut(FillOrder) = bdCountFromOne:(NumVoxelsIn - bdCountFromZero) ;
  
  
  % Plot results (if desired)
  if (bPlot || nargout < 1)
    TitleStr = sprintf('%s (%dx%d)', mfilename, LengthOrig1, LengthOrig2) ;
    Ordering.Utils.PlotIdxPerPosAndTraj(IdxPerPosOut, IOut, JOut, ...
                                        bCountFromZero, TitleStr)
  end
  
end
