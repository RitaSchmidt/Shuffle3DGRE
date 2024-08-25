function [bSample, bRef, bImagingAndRef] = GrappaMasks(N1, N2, ...
                                                       Acceleration1, ...
                                                       Acceleration2, ...
                                                       RefLines1, ...
                                                       RefLines2)

% [bSample, bRef, bImagingAndRef] = GrappaMasks(N1, N2, ...
%                                               Acceleration1, ...
%                                               Acceleration2, ...
%                                               RefLines1, ...
%                                               RefLines2)
%
% Returns logical masks/matrices of dimension N1xN2 for implementing GRAPPA
% subsampling of a k-space grid.
%
% Inputs:
% -------
% - N1 - The size of the 1st dimension of the k-space grid/matrix. (Number
%   of rows in the grid).
% - N2 - The size of the 2nd dimension of the k-space grid/matrix. (Number
%   of columns in the grid).
% - Acceleration1 - A positive integer(!) defining the acceleration along
%   the 1st dimension. The step from one sampled row to the next will be
%   Acceleration1.
% - Acceleration2 - A positive integer(!) defining the acceleration along
%   the 2nd dimension. The step from one sampled column to the next will be
%   Acceleration2.
% - RefLines1 - Number of rows(!) centered around the "center row" that will
%   be fully sampled (used as reference to solve the missing/skipped
%   samples). The "center row" is defined via FFTCenterIndex().
% - RefLines2 - Number of columns(!) centered around the "center column"
%   that will be fully sampled (used as reference to solve the
%   missing/skipped samples). The "center column" is defined via
%   FFTCenterIndex().
%
% Outputs:
% --------
% - bSample - Logcial masks of which samples are actually sampled and which
%   not.
% - bRef - Logical masks of all sampled point that shall serve as reference
%   for solving the missed/skipped samples.
% - bImagingAndRef - Logical mask of samples that are both(!) part of the
%   reference smaples and would also be sampled if no reference scan.

  %% Helper (anonymous) functions
  
  % Checks if a number is practically an integer:
  IsIntValue = @(n) abs(rem(n,1)) <= eps(1) ;

  %% Define parameters

  % Do we have an acceletation in each direction?
  bAcceleration1 = (Acceleration1 > 1) ;
  bAcceleration2 = (Acceleration2 > 1) ;

  % define parameters, only if there is acceleration
  if (bAcceleration1)
    % Center when counting from 1
    Dir1CenterIdx = Ordering.Utils.FFTCenterIndex(N1) ; 

    % Limits of indices of GRAPPA reference region (center of k-space),
    % where k=0 matches index of zero. 

    % Center when counting from 1
    Dir1RefCenterIdx = Ordering.Utils.FFTCenterIndex(RefLines1) ; 
    
    Dir1RefIdxMax = min(N1, (Dir1CenterIdx + (RefLines1 - Dir1RefCenterIdx)) ) ;
    Dir1RefIdxMin = max(1, (Dir1CenterIdx + (1 - Dir1RefCenterIdx)) ) ;

  end

  if (bAcceleration2)
    % Center when counting from 1
    Dir2CenterIdx = Ordering.Utils.FFTCenterIndex(N2) ; 

    % Limits of indices of GRAPPA reference region (center of k-space),
    % where k=0 matches index of zero. 

    % Center when counting from 1
    Dir2RefCenterIdx = Ordering.Utils.FFTCenterIndex(RefLines2) ; 

    Dir2RefIdxMax = min(N2, (Dir2CenterIdx + (RefLines2 - Dir2RefCenterIdx)) ) ;
    Dir2RefIdxMin = max(1, (Dir2CenterIdx + (1 - Dir2RefCenterIdx)) ) ;

  end
  
  %% Define masks along each dimension separately

  if (bAcceleration1)
    % Subsample with a step of Acceleration1, but ensure we include k=0.
    bImaging1 = IsIntValue(((1:N1)-Dir1CenterIdx)/Acceleration1).' ; % column
  
    % Center points within
    bRef1 = false(N1, 1) ; % column
    bRef1(Dir1RefIdxMin:Dir1RefIdxMax) = true ;
  else
    bImaging1 = true(N1, 1) ; % column
    bRef1 = false(N1, 1) ; % column
  end

  if (bAcceleration2)
    % Subsample with a step of Acceleration2, but ensure we include k=0.
    bImaging2 = IsIntValue(((1:N2)-Dir2CenterIdx)/Acceleration2) ; % row
  
    % Center points within
    bRef2 = false(1, N2) ; % row
    bRef2(Dir2RefIdxMin:Dir2RefIdxMax) = true ;
  else
    bImaging2 = true(1, N2) ; % row
    bRef2 = false(1, N2) ; % row
  end

  %% Join 1D masks to 2D masks
  
  % Reference/Imaging/Both booleans of dimension N1xN2 (modified below)

  % Mark range at center of k-space to use as reference
  if (bAcceleration1 && bAcceleration2)
    % Only keep the center of both dimensions
    bRef = (bRef1 & bRef2) ;
  else
      % Only one of the dimensions is used, so use it.
      bRef = (bRef1 | bRef2) ;
  end

  % Which samples to sample if no reference is scanned (or scanned
  % separately)
  bImaging = (bImaging1 & bImaging2) ;

  % All sampled positions (reference and/or imaging)
  bSample = (bImaging | bRef) ;

  % Which samples have a dual use, both reference and "imaging". Some
  % reconstruction algorithems need this.
  bImagingAndRef = (bImaging & bRef) ;


end

