function gray = RDCL_GetGray()

% function RDCL_GetGray
%
% Defines gray as midway between white and black.
%
% Arguments:
%   none.
%
% Returns:
%    gray = gray in [r g b];
%
% Example:
%   gray = RDCL_GetGray();
%
% Andrew Cohen
% 5/9/06
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Known Bugs:
%   none
%
% Change History:
%

% Constants
global screen_ptr;

white = WhiteIndex(screen_ptr);
black = BlackIndex(screen_ptr);

gray = (white+black)/2;
if round(gray) == white
    gray = black;
end
