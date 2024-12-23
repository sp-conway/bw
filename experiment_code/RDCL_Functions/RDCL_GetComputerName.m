function computer_name = RDCL_GetComputerName()

% function RDCL_GetComputerName
%
% Get the name of the computer. Mac OSX dependent.
%
% Arguments:
%   none.
%
% Returns:
%    computer_name = the name of the computer
%
% Example:
%   computer_name = RDCL_GetComputerName();
%
% Andrew Cohen
% 5/5/06
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Known Bugs:
%   none
%
% Change History:
% 5/9/06: Changed method for getting comptuer name to Screen function
% 7/25/11: Added method for unix
% 4/20/13 alc Get rid of newline at end of computer_name

if ismac

  info = Screen('Computer');
  computer_name = info.machineName;

elseif isunix

  [status computer_name] = system('hostname');
	computer_name = computer_name(1:(end-1));

end
