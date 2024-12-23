function RDCL_DrawDotPattern(dotPattern, varargin)

%
% Function: RDCL_DrawDotPattern
%
% Draws a dot pattern
%
% Arguments:
%   dotPattern = the dot pattern to draw.
%   varargin:
%     'DrawDots' = whether to draw the dots.
%     'DrawLines' = whether to draw the lines.
%     'DotColor' = the dot color.
%     'LineColor' = the line color.
%     'DotSize' = the dot size.
%     'LineWidth' = the line width.
%     'Scale' = scaling factor.
%     'Loc' = screen location of the center of the ORIGINAL grid, i.e.,
%       after perturbation, some points may be outside this grid.
%     'Connect' = connect first and last dots with line?
%
% Return:
%   none
%
% Example:
%   dotPattern = RDCL_GenerateDotPattern([500 500], 9);
%   RDCL_DrawDotPattern(dotPattern);
%
% Notes:
%
% Andrew Cohen
% 6/25/08
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Change history:
%

% Constants
global screen_ptr

[screen_middle_x screen_middle_y] = RDCL_GetScreenMiddle();
gray = RDCL_GetGray();

nDots = size(dotPattern{1}, 1);

% Defaults
drawDots = 1;
drawLines = 1;
dotColor = [0 0 0];
lineColor = [0 0 0];
dotSize = 20;
lineWidth = 1;
scaleSize = 1;
loc = [screen_middle_x screen_middle_y];
connect = 1;

% Get user options
if nargin > 2
    for k = 1:2:length(varargin)

        option_text = varargin{k};
        if strcmp(option_text, 'DrawDots')
            drawDots = varargin{k+1};
        elseif strcmp(option_text,'DrawLines')
            drawLines = varargin{k+1};
        elseif strcmp(option_text, 'DotColor')
            dotColor = varargin{k+1};
        elseif strcmp(option_text, 'LineColor')
            lineColor = varargin{k+1};
        elseif strcmp(option_text, 'DotSize')
            dotSize = varargin{k+1};
        elseif strcmp(option_text, 'LineWidth')
            lineWidth = varargin{k+1}(1);
        elseif strcmp(option_text, 'Scale')
            scaleSize = varargin{k+1};
        elseif strcmp(option_text, 'Loc')
            loc = varargin{k+1};
        elseif strcmp(option_text, 'Connect')
            connect = varargin{k+1};
        else
            error('Incorrect RDCL_DrawDotPattern option');
        end

    end
end

% Scale pattern
x_cent = round(dotPattern{2}(1)/2);
y_cent = round(dotPattern{2}(2)/2);

newDotPattern = dotPattern{1} - repmat([x_cent y_cent], [nDots 1]);
newDotPattern = newDotPattern*scaleSize;
newDotPattern = newDotPattern + repmat([x_cent y_cent], [nDots 1]);

% Shift locations
newDotPattern = newDotPattern + repmat([loc(1)-x_cent loc(2)-y_cent], [nDots 1]);

% Draw the lines
if drawLines
    for p = 1:nDots-1
        Screen('DrawLine', screen_ptr, lineColor, newDotPattern(p, 1), ...
            newDotPattern(p, 2), newDotPattern(p+1, 1), ...
            newDotPattern(p+1, 2), lineWidth);
    end

    if connect == 1
        Screen('DrawLine', screen_ptr, lineColor, newDotPattern(1, 1), ...
            newDotPattern(1, 2), newDotPattern(nDots, 1), ...
            newDotPattern(nDots, 2), lineWidth);
    end
end

% Draw the dots
if drawDots
    for p = 1:nDots
        Screen('FillOval', screen_ptr, dotColor, ...
            [newDotPattern(p, 1) - round(dotSize/2) ...
            newDotPattern(p, 2) - round(dotSize/2) ...
            newDotPattern(p, 1) + round(dotSize/2) ...
            newDotPattern(p, 2) + round(dotSize/2)]);
    end
end
