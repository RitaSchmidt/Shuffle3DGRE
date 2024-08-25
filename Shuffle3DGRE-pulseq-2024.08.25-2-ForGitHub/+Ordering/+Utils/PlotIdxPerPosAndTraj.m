function PlotIdxPerPosAndTraj(IdxPerPosOut, IOut, JOut, ...
                              bCountFromZero, TitleStr)
% PlotIdxPerPosAndTraj(IdxPerPosOut, IOut, JOut, Length1, Length2, ...
%                      bCountFromZero, TitleStr)
%
% Helper function to plot a Length1 x Length1 array containig the sampling
% order (IdxPerPosOut) as well as the trajectory as a line plot.
%
% Inputs:
% -------
% - IdxPerPosOut - A 2D array. The values of the elements determine their
%   sampling order. Values are assumed to count (without gaps) from either
%   1 or 0.
% - IOut - [IOut(n), JOut(n)] is the n-th position in the matrix to be
%   sampled (i.e., subscripts of elements inside the matrix). IOut gives
%   the subscript along the 1st dimension. Can either start counting from
%   zero, or from one.
% - JOut - [IOut(n), JOut(n)] is the n-th position in the matrix to be
%   sampled (i.e., subscripts of elements inside the matrix). JOut gives
%   the subscript along the 2nd dimension. Can either start counting from
%   zero, or from one.
% - bCountFromZero - If true counting in IdxPerPosOut, IOut, and JOut is
%   done from zero, otherwise, from 1 (Matlab style).
% - TitleStr - a short string to appear at the title of the resulting
%   figure (in addition to other information).

  [Length1, Length2] = size(IdxPerPosOut) ;
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
    title(TitleStr) ;
    xlabel('J (''dimension 2'')') ;
    ylabel('I (''dimension 1'')') ;
    % Orient axes so that 1st dimension is going from top to bottom (and
    % 2nd dimension from left to right).
    axis ij

end
