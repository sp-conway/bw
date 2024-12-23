function RDCL_Finish(varargin)

% function RDCL_Finish
%
% End the experiment.
%
% Arguments:
%   varargin = minimum priority level (usually 0);
%
% Returns:
%   none.
%
% Example:
%   RDCL_Finish(0);
%
% Andrew Cohen
% 5/8/06
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Known Bugs:
%   none
%
% Change History:
%   7/26/11 - alc - Removed call to RDCL_Constants.
%     hardcoded default min_priority_level to 0.
%

if nargin > 0
    min_priority_level = varargin{1};
else
    min_priority_level = 0;
end

ShowCursor;  % Show the cursor
fclose('all');  % Close the output files
Screen('CloseAll'); % Close the psychtoolbox screens
Priority(min_priority_level);
Snd('Close'); % Close the sound channel
