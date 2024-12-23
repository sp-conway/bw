function RDCL_WaitForKeyPress(varargin)

% function RDCL_WaitForKeyPress
%
% Waits for any keypress.
%
% Arguments:
%   varargin:
%     input_device = the number of the device that will get button press.
%
% Returns:
%    none
%
% Example:
%   RDCL_WaitForKeyPress(1);
%
% Andrew Cohen
% 5/5/06
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Known Bugs:
%   none
%
% Change History:
% 12/6/06 - Michael Ross - modified to directly call low-level KbWait
% Psychtoolbox function.
% 7/27/11 - alc - Removed RDCL_Constants

if nargin > 0
    KbWait(varargin{1});
else
    KbWait();
end

return;
