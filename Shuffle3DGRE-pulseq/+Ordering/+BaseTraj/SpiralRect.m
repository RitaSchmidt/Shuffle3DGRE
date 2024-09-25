function [IdxPerPosOut, IOut, JOut] = SpiralRect(Length1, Length2, ...
                                                 bInsideOut, ...
                                                 bElliptic, ...
                                                 bCountFromZero, ...
                                                 bPlot)
                                                    
% [IdxPerPosOut, IOut, JOut] = RectSpiral(Length1, Length2, ...
%                                         bInsideOut, ...
%                                         bElliptic, ...
%                                         bCountFromZero, ...
%                                         bPlot)
% 
% Order the elements in a matrix to mimic as an outside-in rectangular
% spiral. Shuffle the order of elements in a matrix to mimic as an
% inside-out rectangular spiral. The aspect ratio of the spiral is set to
% fit exactly inside the matrix without any jumps (no need to 'lift the
% pen' drawing the spiral). If an ellipric mask is applied only the
% elements inside are kept and then jumps allowed (skipping the removed
% elements).

%
% Inputs:
% -------
% - Length1 - Size of desired matrix along first dimension. (For a matrix,
%   this is the length of each column.)
% - Length2 - Size of desired matrix along second dimension. (For a matrix,
%   this is the length of each row.)
%   [Default: Length1]
% - bInsideOut - If true starts from center of matrix and winds out.
%   Otherwise, starts outside and winds inward.
%   If true, will not (necessarily) start from the center, but a bit off,
%   so can cover the whole matrix by unwinding. (Just as in the reverese
%   order will not end at the center.)
%   [Default: true]
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
%   step in which the curve passes through the element (counts from 1
%   unless bCountFromZero is true). 
% - IOut, JOut - the subscript order in which the curve passes through,
%   i.e., the coordinates of the curve. IOut gives the Length1 coordinates
%   and JOut is the Length2 coordnates (both counting from 1 unless
%   bCountFromZero is true). 

warning('Work in progress!')

%% Initialize/Defaults

  % Defaults
  bInsideOutDefault = true ;
  bEllipticDefault = false ;
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
  
  %% Set up eliiptical mask (or full k-space)
  
  if(bElliptic) % setup elliptical mask
   Mask = Ordering.Utils.EllipticMask(Length1, Length2) ;
  else
    Mask = true(Length1, Length2) ;
  end

  
  %% Fill 

  % Total number of elements in our matrix
  NIn = nnz(Mask) ;
  
  if (bCountFromZero)
    CountFromZeroCorr = 0 ;
  else
    CountFromZeroCorr = 1 ;
  end
  
  % Initialize our counter
  Counter = 0 ;
  % We always fill the matrix from the outside inwards
  if (bInsideOut)
    StepSign = - 1 ;
    StepCurr = NIn + CountFromZeroCorr ;
  else
    StepSign = + 1 ;
    StepCurr = CountFromZeroCorr - 1 ;
  end
  % Initialize IdxPerPosOut
  IdxPerPosOut = -1*ones(Length1, Length2) ;
  
  MMin = 1 ;
  MMax = Length1 ;
  NMin = 1 ;
  NMax = Length2 ;
  
  % Direction of advancemtn index: 0:+I ; 1:+J ; 2:-I ; 3:-J
  DirIdx = 0 ; 
  
%   [I0, J0] = ind2sub([Length1, Length2], find(Mask(:), 1, 'first')) ; 
%   PosCurr = [I0, J0] ;
  PosCurr = [1,1] ;
  
  while (Counter < NIn)
    switch DirIdx
      case 0 % +I
        FillMin = max(PosCurr(1), find(Mask(:, PosCurr(2)), 1, 'first')) ;
        FillMax = min(MMax, find(Mask(:, PosCurr(2)), 1, 'last')) ;
%         Delta = MMax - PosCurr(1) + 1 ;
%         IdxPerPosOut(PosCurr(1):MMax, PosCurr(2)) = StepCurr + ...
%                                                     StepSign*(1:Delta) ;
        Delta = FillMax - FillMin + 1 ;
        IdxPerPosOut(FillMin:FillMax, PosCurr(2)) = StepCurr + ...
                                                       StepSign*(1:Delta) ;
        PosCurr = [MMax, PosCurr(2) + 1] ;
        NMin = NMin + 1 ;
      case 1 % +J
        FillMin = max(PosCurr(2), find(Mask(PosCurr(1), :), 1, 'first')) ;
        FillMax = min(NMax, find(Mask(PosCurr(1), :), 1, 'last')) ;
%         Delta = NMax - PosCurr(2) + 1 ;
%         IdxPerPosOut(PosCurr(1), PosCurr(2):NMax) = StepCurr + ...
%                                                     StepSign*(1:Delta) ;
        Delta = FillMax - FillMin + 1 ;
        IdxPerPosOut(PosCurr(1), FillMin:FillMax) = StepCurr + ...
                                                       StepSign*(1:Delta) ;
        PosCurr = [PosCurr(1) - 1, NMax] ;
        MMax = MMax - 1 ;
      case 2 % -I
        FillMin = max(MMin, find(Mask(:, PosCurr(2)), 1, 'first')) ;
        FillMax = min(MMax, find(Mask(:, PosCurr(2)), 1, 'last')) ;
%         Delta = PosCurr(1) - MMin + 1 ;
%         IdxPerPosOut(PosCurr(1):-1:MMin, PosCurr(2)) = StepCurr + ...
%                                                        StepSign*(1:Delta) ;
        Delta = FillMax - FillMin + 1 ;
        IdxPerPosOut(FillMax:-1:FillMin, PosCurr(2)) = StepCurr + ...
                                                          StepSign*(1:Delta) ;
        PosCurr = [MMin, PosCurr(2) - 1] ;
        NMax = NMax - 1 ;
      case 3 % -J
        FillMin = max(NMin, find(Mask(PosCurr(1), :), 1, 'first')) ;
        FillMax = min(NMax, find(Mask(PosCurr(1), :), 1, 'last')) ;
%         Delta = PosCurr(2) - NMin + 1 ;
%         IdxPerPosOut(PosCurr(1), PosCurr(2):-1:NMin) = StepCurr + ...
%                                                        StepSign*(1:Delta) ;
        Delta = FillMax - FillMin + 1 ;
        IdxPerPosOut(PosCurr(1), FillMax:-1:FillMin) = StepCurr + ...
                                                          StepSign*(1:Delta) ;
        PosCurr = [PosCurr(1) + 1, NMin] ;
        MMin = MMin + 1 ;
      otherwise
        error('DirIdx is expected to be 0, 1, 2, or 3. Something is wrong') ;
    end
    
    Counter = Counter + Delta ;
    StepCurr = StepCurr + StepSign*Delta ;
    DirIdx = mod(DirIdx + 1, 4) ;
  end
  
  % Generate IOut, JOut.
  [~, SortOrder] = sort(IdxPerPosOut(:)) ;
  % Only keep indices within the mask (with total of NIn elements).
  SortOrder = SortOrder((end-NIn+1):end) ;
  [IOut, JOut] = ind2sub([Length1, Length2], SortOrder) ;
  
  % Corrections due to bCountFromZero (true or false).
  if (bCountFromZero)
    IOut = IOut - 1 ;
    JOut = JOut - 1 ;
  else
    % Do nothing
  end
  
  
  %% Plot results (if desired)
  if (bPlot || nargout < 1)
    TitleStr = sprintf('%s (%dx%d)', mfilename, Length1, Length2) ;
    Ordering.Utils.PlotIdxPerPosAndTraj(IdxPerPosOut, IOut, JOut, ...
                                        bCountFromZero, TitleStr)
  end
  
end
