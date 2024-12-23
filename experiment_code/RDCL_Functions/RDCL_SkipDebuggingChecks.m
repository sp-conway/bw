function RDCL_SkipDebuggingChecks(skip)

% function RDCL_SkipDebuggingChecks
%
% Skips sync tests and sets visual debug level to low to speed up
%   start of experiment.
%
% Do NOT skip the checks when testing your experiment.
%
% Arguments:
%    skip = whether to skip ('TRUE' is yes).
%
% Returns:
%    none.
%
% Example:
%   RDCL_SkipDebuggingChecks('TRUE');
%
% Andrew Cohen
% 5/15/06
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Known Bugs:
%   none
%
% Change History:
%   7/25/11: Set default behavior if not skipped.
%            Added Verbosity and SuppressAllWarnings
%

if strcmp('TRUE', skip)
	% For experiment running

    Screen('Preference', 'SkipSyncTests', 1);
    Screen('Preference', 'VisualDebugLevel', 1);
	Screen('Preference', 'Verbosity', 1);
	Screen('Preference', 'SuppressAllWarnings', 1);

    disp('PACLab-WARNING:  YOU ARE SKIPPING THE SYNC TEST, VISUAL DEBUG LEVEL IS LOW, VERBOSITY IS LOW, AND ALL WARNING SUPPRESSED!!!');

else
	% For experiment testing

    Screen('Preference', 'SkipSyncTests', 0);
	Screen('Preference', 'VisualDebugLevel', 4);
	Screen('Preference', 'Verbosity', 3);
	Screen('Preference', 'SuppressAllWarnings', 0);

end
