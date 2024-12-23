function system_info = RDCL_GetSystemInfo()

% function RDCL_GetSystemInfo
%
% Gets system information
%
% Arguments:
% 	 none
%
% Returns:
%    system_info.
%		os = operating system
%		PTB_ver = PTB version
%  		language = e.g., octave
%		language_ver = langauge version
%		RDCL_ver = PACLab version
%		RDCL_date = date for PACLab version
%		comp_name = computer name
%
% Example:
%   system_info = RDCL_GetSystemInfo();
%
% Andrew Cohen
% 7/28/11
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Notes:
%
% Known Bugs:
%   none
%
% Change History:

v = Screen('Version');

system_info.os = v.os;
system_info.PTB_ver = v.version;
system_info.language = v.language;
system_info.language_ver = version;
system_info.RDCL_ver = RDCL_Version().version;
system_info.RDCL_date = RDCL_Version().date;
system_info.comp_name = RDCL_GetComputerName();
