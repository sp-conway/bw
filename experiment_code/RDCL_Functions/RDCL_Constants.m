%
% RDCL_Constants
%
% Holds constants for PACLab functions
%
% Note:
%   Assumes screen_ptr has been set as a global.
%
% Andrew Cohen
% 5/6/06
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Change History:
% 12/6/06 - Michael Ross - removed def_input_device, was incorrect for
%   version 3.0.x of Psychtoolbox
% 7/25/11 - alc - updated for use in linux/octave. Changed data_file_postfix
%   to .mat. Removed data_file_label, header_file_label, data_folder, data_file_postfix;
%   removed colors; removed some text constants.
% 7/27/11 - alc - moved constants for DisplayTextFile to that m-file.

global screen_ptr

% Escape keys (to quit experiment)
if ismac
	esc_keys = {'ESCAPE' 'LeftShift' 'RightShift'};
elseif isunix
	esc_keys = {'Escape' 'Shift_L' 'Shift_R'};
end

% Screen capture
def_screen_capture_type = 'bmp';
