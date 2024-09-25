function [IdxPerPosOut, IOut, JOut] = SpiralSquare(Length1, Length2, ...
                                                   bInsideOut, ...
                                                   bElliptic, ...
                                                   bCountFromZero, ...
                                                   bPlot)

% [IdxPerPosOut, IOut, JOut] = Spiral(Length1, Length2, ...
%                                     bInsideOut, ...
%                                     bElliptic, ...
%                                     bCountFromZero, ...
%                                     bPlot)
% 
% Order the elements in a matrix to mimic as an inside-out square spiral.
% If the spiral reaches the end of the matrix , continues as if the matrix
% is a window in to a larger matrix where the spiral fits (adds either
% columns alternatingly on either side or rows alternatingly at the bottom
% top and bottom. Below is an illustration:
%
%               ˄ ˄  ˄ ˃-˃-˃-˃  ˅ ˅ 
%               | |  | |     |  | |    
%               ˄ ˄  ˄ ˄ ˃-˃ ˅  ˅ ˅ 
%               | |  | | | | |  | |     
%               ˄ ˄  ˄ ˄ . ˅ ˅  ˅ ˅
%               | |  | |   | |  | |                      
%               ˄ ˄  ˄ ˄-˂-˂ ˂  ˅ ˅
%               | |  |       |  | |   
%               ˄ ˄  ˂-˂-˂-˂-˂  ˂ ˅
%
% Inputs:
% -------
% - Length1 - Size of desired matrix along first dimension. (For a matrix,
%   this is the length of each column.)
% - Length2 - Size of desired matrix along second dimension. (For a matrix,
%   this is the length of each row.)
%   [Default: Length1]
% - bInsideOut - If true starts from center of matrix and winds out.
%   Otherwise, starts outside and winds inward (reverse order).
%   [Default: true]
% - bElliptic - If true, only an elliptical region in k-space/matrix is
%   covered. 
%   [Default: false]
% - bCountFromZero - Whether outputs count from zero (C++) or one (Matlab). 
%   Indices include the coordinates in final matrix (IOut, JOut) AND the
%   step counter within IdxPerPosOut matrix. 
%   [Default: false]
% - bPlot - If true, will plot the a line plot of the curve on top of the
%   IdxPerPosOut. (Accounts for bCountFromZero)
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
  
  % Ensure bCountFromZero is true or false, only
  bCountFromZero = logical(bCountFromZero) ;
  
  % Round integer values
  Length1 = round(Length1) ;
  Length2 = round(Length2) ;

  % Initialize starting position for inside-out trajectory (counting from
  % 1).
  I0 = Ordering.Utils.FFTCenterIndex(Length1) ;
  J0 = Ordering.Utils.FFTCenterIndex(Length2) ;
  PosIdx0 = sub2ind([Length1, Length2], I0, J0) ;


  %% Set up eliiptical mask (or full k-space)
  
  if(bElliptic) % setup elliptical mask
   Mask = Ordering.Utils.EllipticMask(Length1, Length2) ;
  else
    Mask = true(Length1, Length2) ;
  end



  %% Fill IdxPerPosOut (inside  out)

  % Initialize
  % ----------

  IdxPerPosOut = -1*ones(Length1, Length2) ;
  
  if bInsideOut
    SampleCounter =  -bCountFromZero ;
    TrajSign = 1 ;
  else
    SampleCounter =  nnz(Mask) + 1 - bCountFromZero ;
    TrajSign = -1 ;
  end
  
  % Chronological list of positions in the IdxPerPosOut matrix.
  Positions = zeros(nnz(Mask),1) ;
  if (Mask(PosIdx0))
    SampleCounter = SampleCounter + TrajSign ;
    IdxPerPosOut(PosIdx0) = SampleCounter ;
    % Positions count (for now, from 1). 
    % Used to set IOut and JOut using ind2sub(), and matalb expects that
    % indices are counted from 0. We will shift IOut and JOut after we get
    % them from ind2sub().
    Positions(SampleCounter + bCountFromZero) = PosIdx0 ;
  end

      
  
    
  StepSize = 1 ;
  StepDirList = [-1,  0
              0, -1
              1,  0
              0,  1];

  PosIdxCurr = PosIdx0 ;
  DirCounter = 1 ;

  StepLimits = [Length1, Length2] ;
  
  while (StepSize < StepLimits(logical(StepDirList(DirCounter,:))))
  
      PositionsFill = PosIdxCurr + ...
                  StepDirList(DirCounter, 1)*(1:StepSize) + ...
                  StepDirList(DirCounter, 2)*Length1*(1:StepSize) ;
      PosIdxNext = PositionsFill(end) ;
  
      % Apply Mask on Postitions
      PositionsFill = PositionsFill(Mask(PositionsFill(:))) ;
      % Update OrderMat
      NumPositions = numel(PositionsFill) ;
      IdxPerPosOut(PositionsFill(:)) = SampleCounter + TrajSign*(1:NumPositions) ;
      Positions(SampleCounter + ...
                TrajSign*(1:NumPositions) + ...
                bCountFromZero) = PositionsFill(:) ;
      SampleCounter = SampleCounter + TrajSign*NumPositions ;
      
      LastStepSize = StepSize ;
      if (mod(DirCounter, 2) == 0)
        StepSize = StepSize + 1 ;
      end
      DirCounter = mod(DirCounter, 4) + 1 ;
  
      PosIdxCurr = PosIdxNext ; % for next round
  
  end
  
  % We now switch to a different filling 
  
  
  SkipStep = LastStepSize + 1 ;
  MaxLength = max(Length1, Length2) ;
  
  StepDir = StepDirList(DirCounter, :) ;
  SkipDirCounter = mod(DirCounter, 4) + 1 ;
  SkipDir = StepDirList(SkipDirCounter, :) ;
  
  StepSizeNew = StepSize ;
  StepSize = StepSize -1 ;
  
  while (SkipStep <= MaxLength)
  
      PositionsFill = PosIdxCurr + ...
                  StepDir(1)*(1:StepSize) + ...
                  StepDir(2)*Length1*(1:StepSize) ;
      PosIdxNext = PositionsFill(end) + ...
                   SkipDir(1)*SkipStep + ...
                   SkipDir(2)*Length1*SkipStep + ...
                   (StepDir(1) + StepDir(2)*Length1) ;
  
  
      % Apply Mask on Postitions
      PositionsFill = PositionsFill(Mask(PositionsFill(:))) ;
      % Update OrderMat
      NumPositions = numel(PositionsFill) ;
      IdxPerPosOut(PositionsFill(:)) = SampleCounter + TrajSign*(1:NumPositions) ;
      Positions(SampleCounter + ...
                TrajSign*(1:NumPositions) + ...
                bCountFromZero) = PositionsFill(:) ;
      SampleCounter = SampleCounter + TrajSign*NumPositions ; 
      
      PosIdxCurr = PosIdxNext ; % for next round
      SkipDir = -SkipDir ;
      SkipStep = SkipStep + 1 ;
      % reverese direction ;
      StepDir = -StepDir ;
      StepSize = StepSizeNew ;
  
  
  end
  
  [IOut, JOut] = ind2sub([Length1, Length2], Positions(:)) ;
  % "Fix" IOut and JOut if we count from zero
  if (bCountFromZero)
    IOut = IOut - bCountFromZero ;
    JOut = JOut - bCountFromZero ;
  end

  %% Plot results (if desired)
  if (bPlot || nargout < 1)
    TitleStr = sprintf('%s (%dx%d)', mfilename, Length1, Length2) ;
    Ordering.Utils.PlotIdxPerPosAndTraj(IdxPerPosOut, IOut, JOut, ...
                                        bCountFromZero, TitleStr)
  end
  
end
