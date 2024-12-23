function num_frames = RDCL_Sec2Frames(time, ifi)

% function RDCL_Sec2Frames
%
% Converts seconds to number of frames.
%
% Arguments:
%    time = the stimulus time in seconds.
%    ifi = interflip interval
%
% Returns:
%    num_frames = the number of frames.
%
% Example:
%   num_frames = RDCL_Sec2Frames(1.00, ifi);
%
% Andrew Cohen
% 5/15/06
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Known Bugs:
%   none
%
% Change History:

num_frames = round(time ./ ifi);
