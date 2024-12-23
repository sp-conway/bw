function out_file_ptr = RDCL_OpenOutputFiles(sub_file_name, exp_name, header_text);

% function RDCL_OpenOutputFiles
%
% Opens output files: data and header, if needed. Writes to header file
%
% Arguments:
%   sub_file_name = name of subject data file.
%   exp_name = name of the experiment.
%   varargin:
%     'HeaderText' = text for header, a cell with each entry a line of
%       text.
%
% Returns:
%    out_file_ptr = data output file pointer
%
% Example:
%
%
% Andrew Cohen
% 5/5/06
%
% Known Bugs:
%   none
%
% Change History:
%

% Constants
RDCL_Constants

% Open output file
out_file_ptr = fopen(sub_file_name, 'w');

% Create prefixes
out_header_prefix = [data_folder exp_name header_file_label];   % What to call the header file (to store info about data format, etc)

% Write to header file. One per experiment
tmp_file_name = [out_header_prefix data_file_postfix];

if exist(tmp_file_name) == 0

    out_header_ptr = fopen(tmp_file_name, 'w');

    for k = 1:length(header_text)

        fprintf(out_header_ptr, '%s\n', header_text{k});

    end

    fclose(out_header_ptr);

end

