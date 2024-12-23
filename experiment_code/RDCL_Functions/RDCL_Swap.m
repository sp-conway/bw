function [new_x1 new_x2] = RDCL_Swap(x1, x2)

%
% Function: RDCL_Swap
%
% Swaps the values of two variables.
%
% Arguments:
%   x1, x2 = the two variables to swap.
%
% Return:
%   new_x1, new_x2 = the swapped variables.
%
% Example:
%   [new_x1 new_x2] = RDCL_Swap(x1, x2);
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

new_x1 = x2;
new_x2 = x1;
