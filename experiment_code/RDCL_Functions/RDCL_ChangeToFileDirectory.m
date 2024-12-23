function RDCL_ChangeToFileDirectory(file_name)

% function RDCL_ChangeToFileDirectory
%
% Changed to the directory of the named file.
%
% Arguments:
%   file_name = the name of the file.
%
% Returns:
%    none
%
% Example:
%   RDCL_ChangeToFileDirectory(mfilename);
%
% Andrew Cohen
% 5/5/06
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Known Bugs:
%   none
%
% Change History:
%

main_file = sprintf("%s.m", file_name);
eval(sprintf("main_path=which(\'%s\');", main_file));
main_path = main_path(1:end-length(main_file));
cd(main_path);
