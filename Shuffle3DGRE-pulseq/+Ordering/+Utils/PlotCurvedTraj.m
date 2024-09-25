function PlotCurvedTraj(I, J, StepSizeOfMaxAmp, MaxAmp, Marker, ExtraControls)

% PlotCurvedTraj(I, J, StepSizeOfMaxAmp, MaxAmp, Marker, ExtraControls)
%
% Plots a trajectory given by Nx2 set of coordinates [I(:), J(:)]. To
% better visualize trajectories which go back and forth on a grid, the
% positions are not connected by straight lines, but by curves. Up to a
% limit the greater the distance between the points, the further away is
% the curve from the stright line connecting the two points.
%
% Inputs:
% -------
% - I - An array of corrdinates along the 1st dimension of a matrix.
%   Internally turned to a 1D column vector. Must have the same number of
%   elements as J.
% - J - An array of corrdinates along the 2nd dimension of a matrix.
%   Internally turned to a 1D column vector. Must have the same number of
%   elements as I.
% - StepSizeOfMaxAmp - In general, the further the distance between 
%   two coordinates ([I(m), J(m)] and [I(m+1), J(m+1)], the further away
%   from a straight line is the curve connecting the two. However, this is
%   true only to a maximum distance between the points. Above this
%   (StepSizeOfMaxAmp) the deviation from the striaght line is fixed.
%   (See also MaxAmp.)
%   [Default: std() of all step sizes - arbitrary]
% - MaxAmp - The maximum distance the curve connecting to points can have
%   from straight line connecting the two. Units are the same as those of
%   the inputs I and J (assumed indices in a matrix, but can be any).
%   [Default: 0.4]
% - Marker - If not empty, the points at the coordinates I,J will also be
%   marked with this marker (just as in regular plotting). Note: all points
%   in I,J will be marked (even if not all connections are plotted).
%   [Default: '.']
% - ExtraControls - Structure containing extra controls of the plotting:
%   > MaxConnects - If not empty, plots only the first MaxConnects + 1
%     point in I and J.
%     [Default: [] (empty) => plot all]
%   > ExtraPointsPerConnection - How many points in each curve connecting
%     two I,J points, not counting the edge coordinates.
%   > CurvePower - The curve connecting two point I,J points is a
%     polynomial of order CurvePower (single maxima/minima at center
%     between the two points).
%     [Default: 3]
%   > ColorCycleLength - Cycles through the full colormap every this many
%     path segments.
%     [Default: 64]
%   > Colormap - The colormap (an Nx3 RGB array) to be used.
%     [Default: jet(ColorCycleLength)]
%   > BGColor - Bacground color of plot
%     [Default: [0.8, 0.8, 0.8] - light gray]
%   > LineStyle - The line style to be used to plot the curves.
%     [Default: '-' (solid line)]
%   > LineWidth - The linewidth to be used to plot the curves.
%     [Default: 2]
%   > bPlotConnectionWithPatch - Two options to plot "connections" between 
%     samples, i.e., the trajectory. Either use the patch command, easy to
%     give each connection a different value, and therefore color,
%     according to the colormap). Alternatively, just plot a matrix of
%     x-values and y-values. In this case, each column will be colored
%     differently according to the current "color order".
%     [Default: true]
%   > MarkerSize - Size of markers to plot (if plotted).
%     [Default: 12]
%   > MarkerEdgeColor - Edge color of markers to plot (if plotted).
%     [Default: 'k' - black]
%   > MarkerFaceColor - Face color (fill) markers to plot (if plotted).
%     [Default: 'k' - black]

%% Initial setup - set defaults for missing inputs

  % Set Some defaults
  StepSizeOfMaxAmpDefault = std (sqrt(diff(I).^2 + diff(J).^2)) ;
  MaxAmpDefault = 0.4 ;
  MarkerDefault = '.' ; % empty means no marker
  ExtraControlsDefault.MaxConnects = [] ;
  ExtraControlsDefault.ExtraPointsPerConnection = 8 ;
  ExtraControlsDefault.CurvePower = 3 ;
  ExtraControlsDefault.ColorCycleLength = 64 ;
  if isempty(ExtraControlsDefault.ColorCycleLength)
    ExtraControlsDefault.Colormap = jet ;
  else
    ColormapLengthDefault = ExtraControlsDefault.ColorCycleLength * ...
                            ExtraControlsDefault.ExtraPointsPerConnection ;
    ExtraControlsDefault.Colormap = jet(ColormapLengthDefault) ;
  end
  ExtraControlsDefault.BGColor = 0.8 * [1,1,1] ;
  ExtraControlsDefault.LineStyle = '-' ;
  ExtraControlsDefault.LineWidth = 2 ;
  ExtraControlsDefault.bPlotConnectionWithPatch = true ;
  ExtraControlsDefault.MarkerSize = 12 ;
  ExtraControlsDefault.MarkerEdgeColor = 'k' ;
  ExtraControlsDefault.MarkerFaceColor = 'k' ;
  
  % Handle partial input (use defaults instead)
  switch nargin
    case 2
      StepSizeOfMaxAmp = StepSizeOfMaxAmpDefault ;
      MaxAmp = MaxAmpDefault ;
      Marker = MarkerDefault ;
      ExtraControls = ExtraControlsDefault ;
    case 3
      MaxAmp = MaxAmpDefault ;
      Marker = MarkerDefault ;
      ExtraControls = ExtraControlsDefault ;
    case 4
      Marker = MarkerDefault ;
      ExtraControls = ExtraControlsDefault ;
    case 5
      ExtraControls = ExtraControlsDefault ;
    case 6
      % Do nothing.
    otherwise
      error('Too many or too few input arguments.') ;
  end
      
  % Fill ExtraControls with defaults, if necesary
  if (isempty(ExtraControls))
    ExtraControls = ExtraControlsDefault ;
  else
    ExtraParamNamesCell = fieldnames(ExtraControlsDefault) ;
    for FieldCounter = 1:numel(ExtraParamNamesCell)
      Field = ExtraParamNamesCell{FieldCounter} ;
      if(~isfield(ExtraControls, Field))
        ExtraControls.(Field) = ExtraControlsDefault.(Field) ;
      end
    end
  end
  


%% render

  % -----------------------------------------------------------------------
  % General controls
  % -----------------------------------------------------------------------
  
  % Renderer = 'painters' ; % 'opengl', 'painters'
  bPlotConnectionWithPatch = ExtraControls.bPlotConnectionWithPatch ;
  MaxConnects = ExtraControls.MaxConnects ; % empty plot all connections
  
  
  % Path plot settings
  ColorCycleLength = ExtraControls.ColorCycleLength ;
  bReversePathColor = false ;
  % Controls for curved connections between samples (2D indices)
  % number of extra points plotted between two sample positions
  ExtraPointsPerConnection = ExtraControls.ExtraPointsPerConnection ;
  % power used to define curve (2 means parabola).
  CurvePower = ExtraControls.CurvePower ;
  % Style of line (use 'none' to hide display)
  % LineStyle = 'none' ;
  LineStyle = ExtraControls.LineStyle ;
  LineWidth = ExtraControls.LineWidth ;
  BGColor = ExtraControls.BGColor ;

  % Marker settings
  % Marker = ExtraControls.Marker ; % empty means no marker.
  MarkerSize = ExtraControls.MarkerSize ;
  MarkerEdgeColor = ExtraControls.MarkerEdgeColor ;
  MarkerFaceColor = ExtraControls.MarkerFaceColor ;

  
  
  
  % ------------------------------------
  % Turn jumps in matrix to curved paths
  % ------------------------------------
    
  PathSimple = [I(:), J(:)] ;
  % If we want to limit the number of connections shown, then discard the
  % initial points in the path.
  if (~isempty(MaxConnects) && (MaxConnects+1) < numel(I) )
    PathSimple = PathSimple((end-MaxConnects):end) ;
  end

  % Generate curved path
  PathCurved = PolyLine2InterSpline2D(PathSimple, ExtraPointsPerConnection, ...
                                      MaxAmp, CurvePower, StepSizeOfMaxAmp);
  NSamples = size(PathCurved, 1) ;
      
  % ---------------------                                    
  % path curve
  % --------------------- 
  
  figure ; 


  % Two options to plot "connections" between samples, i.e., the
  % trajectory. Either use the patch command, easy to give each
  % connection a different value, and therefore color, according to the
  % colormap). Alternatively, just plot a matrix of x-values and
  % y-values. In this case, each column will be colored differently
  % according to the current "color order".

  if (bPlotConnectionWithPatch)
    % Switched to patch, which allows the plot color to gradually change.
    % Set color in reverse order, so there will be a good contrast with
    % the background image
%       PathColors = [(reshape(repmat(1:(numel(I)-1), ...
%                                    ExtraPointsPerConnection +1, 1), ...
%                             [], 1))
%                     numel(I)
%                     NaN] ;
    if (bReversePathColor)
%       PathColors(1:(end-1)) = flipud(PathColors(1:(end-1))) ;
      PathColors = [linspace(NSamples, 1, size(PathCurved,1)), ...
                    NaN].' ;
    else
      PathColors = [linspace(1, NSamples, size(PathCurved,1)), ...
                    NaN].' ;
    end

    % modify so we cycle through full colormap every ColorCycleLength
    % path segments
    if (~isempty(ColorCycleLength))
%       PathColors = mod((PathColors-1)*NSamples/ColorCycleLength, ...
%                        NSamples) + 1;
      PathColors = mod((PathColors-1), ColorCycleLength) + 1;
    end


    % Actual plot
    patch([PathCurved(:,2) ; NaN], ...
          [PathCurved(:,1) ; NaN], ...
          PathColors, ...
          'EdgeColor','interp', 'LineWidth', LineWidth, ...
          'LineStyle', LineStyle) ;
  else
    Connections = zeros((ExtraPointsPerConnection+2), (NSamples-1), 2) ;

    Connections(1:(end-1), :, :) = reshape(PathCurved(1:(end-1),:), ...
                                           (1+ExtraPointsPerConnection), ...
                                           (NSamples-1), 2) ;

    Connections(end, :, :) = reshape(PathSimple(2:end, :), ...
                                     1, (NSamples-1), 2) ;


    % Plot 
    plot(Connections(:,:,2), Connections(:,:,1), ...
         'LineWidth', LineWidth, 'LineStyle', LineStyle) ;
  end

  colormap(ExtraControls.Colormap) ;

  axis ij ;
  set(gca, 'Color', BGColor) ; % set backgournd color

  % Add markers?
  if (~isempty(Marker))
    hold on ;
    plot(J(:), I(:), 'LineStyle', 'none', 'Marker', Marker, ...
         'MarkerSize', MarkerSize, ...
         'MarkerEdgeColor', MarkerEdgeColor, ...
         'MarkerFaceColor', MarkerFaceColor) ;
  end
   
  


end


function [PathNew] = PolyLine2InterSpline2D(PathIn, ExtraPointsPerSegment, ...
                                            MaxAmp, CurvePower, StepOfMaxAmp)
% [PathNew] = PolyLine2InterSpline2D(PathIn, ExtrsPointsPerSegment, ...
%                                    MaxAmp, CurvePower, StepOfMaxAmp)
%
% "Interpolates" a polyline, inserting within each segment (between two
% point) ExtraPointsPerSegment which are on a curve between the two. The
% distance of the curve from the straight line between the points depends
% on the distance between the points. The farther away the points are, the
% farther the curve is from the straight line connecting the points.
%
% Inputs:
% -------
% - ??????????
% 
% Outputs:
% --------
% No outputs 

  % -----------------------------------------------------------------------
  % Sanity checks
  % -----------------------------------------------------------------------

  % Check if input is 2D or not.
  if(isempty(PathIn) || ~ismatrix(PathIn))
    error('PathIn is expected to be a non-empty 2D array') ;
  end
  % Check at least one dimension is of length 2
  if (~any(size(PathIn) == 2))
    error('At least one dimension of PathIn is expected to be of length 2') ;
  end
  
  if (any(size(PathIn) == 1))
    % nothing to "interpolate", just a single point.
    PathNew = PathIn ;
    return ;
  end
  
  % -----------------------------------------------------------------------
  % set defaults
  % -----------------------------------------------------------------------
  
  switch nargin
    case 3
      CurvePower = 3 ;
      StepOfMaxAmp = 0 ;
    case 4
      StepOfMaxAmp = 0 ;
    case 5
      % Do nothing
    otherwise
      error('Too many or too few input arguments') ;
  end

  % -----------------------------------------------------------------------
  % Ensure Nx2 array
  % -----------------------------------------------------------------------
  
  % For convenience work with arrays of size Nx2:
  PermuteDims = [1, 2] ;
  if (size(PathIn,2) > 2)
    PermuteDims = [2, 1] ;
  end
  Path = permute(PathIn, PermuteDims);
  clear PathIn;
  
  % -----------------------------------------------------------------------
  % Simple linear interpolation of segments (starting point)
  % -----------------------------------------------------------------------
  
  % number of input points
  NIn = size(Path, 1) ;
  % number of output points
  NOut = (NIn - 1) * (1 + ExtraPointsPerSegment) + 1;
  
  InterpIndicesIn = (1:NIn).' ;
  InterpIndicesOut = (1:(1/(ExtraPointsPerSegment+1)):NIn).' ;
  
  PathNew = zeros(length(InterpIndicesOut), 2) ;
  
  % interpolate, before any "shifting" of interpolated points
  PathNew(:, 1) = interp1(InterpIndicesIn, Path(:, 1), InterpIndicesOut) ;
  PathNew(:, 2) = interp1(InterpIndicesIn, Path(:, 2), InterpIndicesOut) ;
  
  % -----------------------------------------------------------------------
  % Find "shift vectors" perpendiular to each segment
  % -----------------------------------------------------------------------
  
  % Start with normalized shifts per segment - unit normal per segment
  
  % Matrix to rotate 2D row vector by 90 degrees
  Rotate90Right = [ 0 1
                   -1 0] ;
  % Not normalizes yet
  SegNormals = diff(Path,1, 1) * Rotate90Right ;
  % Normalize
  SegLength =  sqrt(sum((SegNormals).^2, 2)) ;
  SegNormals = SegNormals ./ SegLength ;
  % Scale
  SegNormalScale = SegLength / StepOfMaxAmp ;
  SegNormalScale(SegNormalScale>1) = 1 ;
  SegNormals = SegNormals .* SegNormalScale ;
  
  % Set shifts per interpolated point (normalized for now) 
  
%   ShiftNormals = zeros(NOut-1, 2) ; % all points except last
%   
%   ShiftNormals = ...
%            reshape(ShiftNormals, (ExtrsPointsPerSegment+1), (NIn - 1), 2) ;
  
  % Fill shift per interpolated point, with normalized shift per segment
  ShiftNormals = repmat(reshape(SegNormals, 1, (NIn - 1), 2), ...
                        ExtraPointsPerSegment+1, 1, 1) ;
  
  % Scale shift along each segment
  d = (ExtraPointsPerSegment+1)/2 ;
  % define ShiftScaleFactors in steps
  % set to: d, d - 1, ... 0, 1, ..., d -1 (final d belongs to next segment)
  ShiftScaleFactors = abs((0:ExtraPointsPerSegment).' - d) ;
  % "Normalize": (max is 1)
  ShiftScaleFactors = ShiftScaleFactors/d ;
  % Apply power (to curve it)
  ShiftScaleFactors = ShiftScaleFactors.^CurvePower ;
  % Shift so center is 1, instead of edges
  ShiftScaleFactors = 1 - ShiftScaleFactors ;
  % Scale, to Height ;
  ShiftScaleFactors = MaxAmp * ShiftScaleFactors ;
  
  
%   % set to: 0, 1, ..., d, d-1, ..., 1
%   ShiftScaleFactors = d -abs(ShiftScaleFactors) ;
%   % Normalize and apply power factor so peak is Height
%   ShiftScaleFactors = (ShiftScaleFactors/d).^CurvePower ;
%   % Set peak to Height (it is now 1)
%   ShiftScaleFactors = Height * ShiftScaleFactors;
  
  % Now we can scale each shift as a function of position along segment.
  ShiftNormals = ShiftNormals .* ShiftScaleFactors ;
  
  % Reshape ShiftNormals:
  ShiftNormals = reshape(ShiftNormals, NOut-1, 2) ;
  
  % Add shift to PathNew to generate (approx.) curve per segment
  PathNew(1:(NOut-1),:) = PathNew(1:(NOut-1),:) + ShiftNormals ;
  
  
  % -----------------------------------------------------------------------
  % Undo dimension permutation
  % -----------------------------------------------------------------------
  PathNew = ipermute(PathNew, PermuteDims);
  
  
  % -----------------------------------------------------------------------
  % The End
  % -----------------------------------------------------------------------
  
  return ;

end