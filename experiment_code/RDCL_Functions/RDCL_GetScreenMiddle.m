function [screen_middle_x screen_middle_y] = RDCL_GetScreenMiddle()

% function RDCL_GetScreenMiddle
%
% Coordinates for the middle of the screen.
%
% Arguments:
%    none.
%
% Returns:
%    screen_middle_x, screen_middle_y: x and y middle coordinates.
%
% Example:
%   [screen_middle_x screen_middle_y] = RDCL_GetScreenMiddle();
%
% Andrew Cohen
% 5/15/06
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Known Bugs:
%   none
%
% Change History:
%

% Constants
global screen_ptr;

screen_rect = Screen('Rect', screen_ptr);

screen_middle_x = round(screen_rect(3)/2);
screen_middle_y = round(screen_rect(4)/2);
