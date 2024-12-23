function RDCL_DrawSquareWithLine(SquareLen, LineOrientation, varargin)

% Function: RDCL_DrawSquareWithLine
%
% Draws a square with a vertical line starting from its base.
%
% Arguments:
%    SquareLen = length of square sides, in pixels. No default.
%    LineOrientation = how is line oriented in square? No default.
%       - options - 'vtop' (starts from top, goes downward), 'vbot' (starts %        from base, goes up), 'rhoz' (starts on right hand side of square, goes left), 'lhoz' (starts on left hand size of square, goes right)
%    varargin:
%       'SquareLocs' = where the square is centered on screen, 1x2 vector.
%           * Default = Center of the screen
%       'SquareCol' = color of the square border in rgb, 1x3 vector
%           * Default = white ([255 255 255])
%       'SquareThickness' = thickness of the square border
%           * Default = 2 pixels
%       'SquareFillCol' = fill color of square
%           * Default = color of the screen
%       'LineLen' = length of the line to be drawn, in pixels.
%           * Default = 1/4 of the square side length
%       'LinePos' = where is line positioned relative to a particular corner?
%           if vtop - distance from top left corner (on x axis)
%           if vbot - distance from bottom left corner (on x axis)
%           if rhoz - distance from bottom right corner (on y axis)
%           if lhoz - distance from top left corner (on y axis)
%       'LineCol' = color of the line in rgb, 1x3 vector
%           * Default = white
%       'LineThickness' = thicknes of the line in pixels.
%           * Default = 2 pixels, Max = 7
%
% Return:
%   none
%
% Example:
%   square = RDCL_DrawSquareWithLine(250, 'rhoz');
%
% Notes:
%
% Sean Conway
% August 2020
%
% Change history:
%   September 2020 - changed function to draw a line projecting inward from any part of square - now called 'RDCL_DrawSquareWithLine'

% Constants -------------------------------------------
global screen_ptr;

[screen_middle_x screen_middle_y] = RDCL_GetScreenMiddle();

% Default values---------------------------------------
    % square_len
    square_len = SquareLen;

    % line_orientation
    line_orientation = LineOrientation;

    % Square location - Default in center of the screen.
    square_locs = [screen_middle_x screen_middle_y];

    % Square color - Default is white.
    square_col =  [255 255 255];

    % Square thickness - Default is 2 pixels
    square_thickness = 2;

    % Square fill color - Have to draw a separate, smaller square inside it.
    % So the default for this parameter is empty until specified to fill.
    square_fill_col = [];

    % Line length - Default is 1/4 length of square side.
    line_len = square_len*.25;

    % Line Position - X Axis. Default is in the center of square.
    line_pos = .5*square_len;

    % Line color - Default is white
    line_col = [255 255 255];

    % Line thickness - Default is 2, Thickest it can be is 8
    line_thickness = 2;

% Varargin---------------------------------------

if nargin > 1
  for k = 1:2:length(varargin)

    user_opts = varargin{k};
    if strcmp(user_opts, 'SquareLocs')
      square_locs = varargin{k+1};
    elseif strcmp(user_opts, 'SquareCol')
      square_col = varargin{k+1};
    elseif strcmp(user_opts, 'SquareThickness')
      square_thickness = varargin{k+1};
    elseif strcmp(user_opts, 'SquareFillCol')
      square_fill_col = varargin{k+1};
    elseif strcmp(user_opts,'LineLen')
      line_len = varargin{k+1};
    elseif strcmp(user_opts, 'LinePos')
      line_pos = varargin{k+1};
    elseif strcmp(user_opts, 'LineCol')
      line_col = varargin{k+1};
    elseif strcmp(user_opts, 'LineThickness')
      line_thickness = varargin{k+1};
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OUTER SQUARE

    % Get Square into a 1x4 vector
    square_vec = [0 0 square_len square_len];

    % Centers the square where the user wants it on the screen
    square = CenterRectOnPointd(square_vec, square_locs(1),square_locs(2));
    Screen('FrameRect', screen_ptr, square_col, square, square_thickness);

% INNER SQUARE (Only applicable if fill color is specified)
    % Draw inner square if fill color is specified
    if (length(square_fill_col)==3)
        inner_square_len = square_len-(2*square_thickness);
        inner_square_vec = [0 0 inner_square_len inner_square_len];
        inner_square = CenterRectOnPointd(inner_square_vec, square_locs(1),square_locs(2));
        Screen('FillRect', screen_ptr, square_fill_col, inner_square);
    end

% Depends on what orientation you want the line at
switch(line_orientation)

% VERTICAL LINE STARTING FROM BOTTOM OF SQUARE
    case 'vbot'
        % Picks out the vertical line's x axis location on this particular square
        line_locs_x = square(1) + line_pos;

        % Find starting and ending point for the line on the y axis
        line_locs_y_bottom = square(4);
        line_locs_y_top = line_locs_y_bottom - line_len;

        % Get starting & ending points (both x & y) in a 1x4 vector
        line_locs_full_vec = [line_locs_x line_locs_y_bottom line_locs_x line_locs_y_top];

        % Draw vertical line
        Screen('DrawLine', screen_ptr, line_col, line_locs_full_vec(1),line_locs_full_vec(2),line_locs_full_vec(3),line_locs_full_vec(4), line_thickness);

% VERTICAL LINE STARTING FROM TOP OF SQUARE
    case 'vtop'
        % Picks out the vertical line's x axis location on this particular square
        line_locs_x = square(1) + line_pos;

        % Find starting and ending point for the line on the y axis
        line_locs_y_top = square(2);
        line_locs_y_bottom = line_locs_y_top - line_len;

        % Get starting & ending points (both x & y) in a 1x4 vector
        line_locs_full_vec = [line_locs_x line_locs_y_bottom line_locs_x line_locs_y_top];

        % Draw vertical line
        Screen('DrawLine', screen_ptr, line_col, line_locs_full_vec(1),line_locs_full_vec(2),line_locs_full_vec(3),line_locs_full_vec(4), line_thickness);

% HORIZONTAL LINE STARTING FROM RIGHT SIDE OF SQUARE
    case 'rhoz'
        % Pick out line x axis location on this particular square
        line_locs_x_right = square(3);

        % Find line end point
        line_locs_x_left = square(3) - line_len;

        % Where does line belong on y axis (i.e., distance from bottom right corner)
        line_locs_y = square(4) - line_pos;

        % Get starting & ending points (both x & y) in a 1x4 vector
        line_locs_full_vec = [line_locs_x_right line_locs_y line_locs_x_left line_locs_y];

        % Draw horizontal line
        Screen('DrawLine', screen_ptr, line_col, line_locs_full_vec(1),line_locs_full_vec(2),line_locs_full_vec(3),line_locs_full_vec(4), line_thickness);

% HORIZONTAL LINE STARTING FROM LEFT SIDE OF SQUARE
    case 'lhoz'

        % Pick out line x axis location on this particular square
        line_locs_x_right = square(1);

        % Find line end point
        line_locs_x_left = square(1) + line_len;

        % Where does line belong on y axis (i.e., distance from bottom right corner)
        line_locs_y = square(4) - line_pos;

        % Get starting & ending points (both x & y) in a 1x4 vector
        line_locs_full_vec = [line_locs_x_right line_locs_y line_locs_x_left line_locs_y];

        % Draw horizontal line
        Screen('DrawLine', screen_ptr, line_col, line_locs_full_vec(1),line_locs_full_vec(2),line_locs_full_vec(3),line_locs_full_vec(4), line_thickness);

% error (line_orientation has no default)
    otherwise
      error('Error: Specify a line orientation.');

end



% end of function