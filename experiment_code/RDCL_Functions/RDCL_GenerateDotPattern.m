function dotPattern = RDCL_GenerateDotPattern(gridSize, nDots)

%
% Function: RDCL_GenerateDotPattern
%
% Generates a random 2D dot pattern.  Dots uniformly distributed over the
%   grid.
%
% Arguments:
%   gridSize = 1x2 vector that gives the horizontal and vertical size of
%     the grid on which to draw the pattern, in pixels.
%   nDots = the number of dots to draw.
%
% Return:
%   dotPattern = a structure containing
%     dotPattern{1} = a matrix of the dot positions in pixels. The lower-
%       left corner of the grid is (0, 0).
%     dotPattern{2} = the gridSize.
%
% Example:
%   dotPattern = RDCL_GenerateDotPattern([100 100], 9);
%
% Andrew Cohen
% 6/24/08
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Change history:
%

% Constants

% Calculate the dot positions
dotPattern{1} = round(rand(nDots, 2).*repmat(gridSize, [nDots, 1]));

% Store the grid size
dotPattern{2} = gridSize;
