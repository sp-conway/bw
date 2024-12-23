function new_val = RDCL_Round(old_val, precision)

% function RDCL_Round
%
% Rounds to a specified precision.
%
% Arguments:
%    old_val = the original value.
%    precision = how many decimal places to use.
%
% Returns:
%    new_val = the rounded value.
%
% Example:
%   RDCL_Round(10.1234, 2);
%
% Andrew Cohen
% 5/25/19
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Known Bugs:
%   none
%
% Change History:
%

new_val = old_val*10^precision;
new_val = floor(new_val);
new_val = new_val*10^-precision;
