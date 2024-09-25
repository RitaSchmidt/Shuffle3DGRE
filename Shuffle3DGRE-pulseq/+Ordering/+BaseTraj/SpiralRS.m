function [IdxPerPosOut, IOut, JOut] = SpiralRS(Length1, Length2, ...
                                               bInsideOut, ...
                                               bElliptic, ...
                                               bCountFromZero, ...
                                               bPlot)
                                                    
% [IdxPerPosOut, IOut, JOut] = SpiralRS(Length1, Length2, ...
%                                       bInsideOut, bCountFromZero, bPlot)
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
%   [Default: true]
% - bElliptic - If true, only an elliptical region in k-space is covered. 
%   [Default: true]
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
  bCountFromZeroDefault = false ;
  bPlotDefault = false ;

  switch nargin
    case 1
      Length2 = Length1 ;
      bInsideOut = bInsideOutDefault ;
      bElliptic = bEllipticDefault ;
      bCountFromZero = bCountFromZeroDefault ; 
      bPlot = bPlotDefault ;
    case 2
      bInsideOut = bInsideOutDefault ;
      bElliptic = bEllipticDefault ;
      bCountFromZero = bCountFromZeroDefault ; 
      bPlot = bPlotDefault ;
    case 3
      bElliptic = bEllipticDefault ;
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
  
  % Ensure bCountFromZero is zero or one, only
  bCountFromZero = logical(bCountFromZero) ;
  % double version of binary variables, for numerics ('Colon operands
  % cannot be logical')
  bdCountFromZero = double(bCountFromZero) ;
  bdCountFromOne  = double(~bCountFromZero) ;
  
  % Round integer values
  Length1 = round(Length1) ;
  Length2 = round(Length2) ;
  
  % Total number of voxels in the matrix.
  NumVoxelsTot = Length1*Length2 ;
  

%   % Sort samples from the inside out. Returns the 1D index of samples in a
%   % matrix, but reordered according to the chornlogical order (of inside
%   % out spiraling).
%   [xx yy]=ndgrid(linspace(-1,1,NFast),linspace(-1,1,NSlow));
%   [thmat rmat]=cart2pol(xx,yy);
%   radstep=sqrt(2)/sqrt(NSlow*NFast);
%   rsteps=0:radstep:sqrt(2);
%   SortOrder=[];
%   for ii=1:numel(rsteps)-1
%       Ir=find(rmat(:)<rsteps(ii+1)&rmat(:)>=rsteps(ii));
%       [~, is] = sort(thmat(Ir)) ;
%       SortOrder=[SortOrder; Ir(is)]; 
%   end
% %   Ir=find(setdiff(1:NSamples,SortOrder));
% %   [~, is] = sort(thmat(Ir)) ;


  
  % Sort samples from the inside out. Returns the 1D index of samples in a
  % matrix, but reordered according to the chornlogical order (of inside
  % out spiraling).

  % get k-space "indices" centered around zero.
  Range1 = (1:Length1) - Ordering.Utils.FFTCenterIndex(Length1) ;
  Range2 = (1:Length2) - Ordering.Utils.FFTCenterIndex(Length2) ;
  MaxAbs1 = max(abs(Range1)) ;
  MaxAbs2 = max(abs(Range2)) ;
  % Changed original code, because typically in k-space (even number of
  % samples), the origin is slightly shifted in the acquisition matrix,
  % i.e. of the form -2 -1 0 1 (for 4 samples) instead of -2 -0.7 0.7 2. 
  % Original code was:
  % [xx, yy]=ndgrid(linspace(-1,1,Length1),linspace(-1,1,Length2)) ;
  % This changes the results
  [xx, yy]=ndgrid(Range1/MaxAbs1, Range2/MaxAbs2) ;
                
  [thmat, rmat]=cart2pol(xx, yy);
  radstep=sqrt(2)/sqrt(NumVoxelsTot);
  if (bElliptic)
    % with a mask we want to limit the "radius" to 1. Expanded slightly
    % beyond 1 (1 + 2*eps(1)) to also include cases of 1. Not a perfect
    % solution.
    rsteps = linspace(0, 1+2*eps(1), 1/radstep + 1) ;
  else % No mask
    % Added an extra value to rsteps which is much beyond sqrt(2), so we
    % won't miss any samples by accident.
    rsteps=[0:radstep:sqrt(2), 2];
  end
  SortOrder = zeros(NumVoxelsTot, 1) ;
  SamplesCounter = 0 ; 
  for ii=1:(numel(rsteps)-1)
      Ir=find(rmat(:)<rsteps(ii+1)&rmat(:)>=rsteps(ii));
      [~, is] = sort(thmat(Ir)) ;
      NumFound = numel(Ir) ;
      if (NumFound > 0)
        SortOrder(SamplesCounter + (1:NumFound)) = Ir(is); 
        SamplesCounter = SamplesCounter + NumFound ;
      end
  end
  
  % Sorting above was from the inside out. If we want to reverse the
  % direction, just flip the resulting vector.
  if (~bInsideOut)
    SortOrder(1:SamplesCounter) = flipud(SortOrder(1:SamplesCounter)) ;
  end
  
  % Fill IdxPerPosOut (only position not masked). Value of -1 are values
  % masked. 
  IdxPerPosOut = zeros(Length1, Length2) - 1 ;
  IdxPerPosOut(SortOrder(1:SamplesCounter)) = ...
                        bdCountFromOne:(SamplesCounter - bdCountFromZero) ;


  
  % Find order of voxels ("coordinates" IOut, JOut), includes points which
  % may be masked later.
  [IOut, JOut] = ind2sub([Length1, Length2], SortOrder(1:SamplesCounter)) ;

  if (bCountFromZero)
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
