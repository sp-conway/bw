function RDCL_DummyCalls()

% function RDCL_DummyCalls
%
% Do dummy calls to GetSecs, WaitSecs, KbCheck to make sure
% they are loaded and ready when we need them - without delays
% in the wrong moment:
%
% Arguments:
%    none
%
% Returns:
%    none.
%
% Example:
%   RDCL_DummyCalls();
%
% Andrew Cohen
% 7/25/11
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Known Bugs:
%   none
%
% Change History:
%

KbCheck;
WaitSecs(0.1);
GetSecs();
