function [text_x text_y] = RDCL_CenterText(text)

%
% Function: RDCL_CenterText
%
% Determines upper left corner to center text on the screen.
%
% Arguments:
%   text = the text to display.
%
% Return:
%   text_x, text_y = the upper left corner of the rect to center the text
%     on the screen.
%
% Example:
%   [text_x text_y] = RDCL_CenterText('hi');
%
% Notes:
%   none.
%
% Andrew Cohen
% 5/17/06
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Change history:
%

% Constants
global screen_ptr

% Center the text
[screen_middle_x screen_middle_y] = RDCL_GetScreenMiddle();
[norm_bounds_rect offset_bounds_rect] = Screen('TextBounds', screen_ptr, text);

text_x = screen_middle_x - round(norm_bounds_rect(3)/2);
text_y = screen_middle_y - round(norm_bounds_rect(4)/2);
