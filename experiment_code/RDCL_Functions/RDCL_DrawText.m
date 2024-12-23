function RDCL_DrawText(x, y, txt, varargin)

% function RDCL_DrawText
%
% Displays text.
%
% Arguments:
%    x, y = x & y coords for center of the text
%    txt = the text to draw
%    varargin:
%      'Col' = color [r g b]. Default is white.
%
% Returns:
%    none
%
% Example:
%   RDCL_DrawText(100, 100, "Hello, world", scr_ptr);
%
% Andrew Cohen
% 8/6/18
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Notes:
%
% Known Bugs:
%   none
%
% Change History:

% Draw_Text %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
global screen_ptr

% Defaults
col = [255 255 255];

% Get user options
if nargin >= 2
    for k = 1:2:length(varargin)

        option_text = varargin{k};
%        disp(option_text);
        if strcmp(option_text, 'Col')
            col = varargin{k+1}(1);
        else
            error('RDCL_DrawText: Incorrect option');
        end

    end
end

[norm_bounds offset_bounds] = Screen('TextBounds', screen_ptr, txt);
Screen('DrawText', screen_ptr, txt, ...
x - round(norm_bounds(3)/2), y - (round(norm_bounds(4)/2)),
col);
