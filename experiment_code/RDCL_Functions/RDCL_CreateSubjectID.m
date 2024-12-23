
function [sub_ID sub_num sub_file_name] = RDCL_CreateSubjectID(exp_name, varargin)

% function RDCL_CreateSubjectID
%
% Creates a unique subject ID.
%
% Arguments:
%   exp_name = the name of the experiment.
%   varargin:
%     'ComputerName' = if want to include computer name in subject ID set
%       to 'TRUE'
%			'DataFolder' = folder for data
%			'DataFileLabel' = label to mark data file, e.g., _data_
%			'DataFilePostfix' = postfix for data file, e.g., .dat
%
% Returns:
%    sub_ID = ID of subject
%    sub_num = subject number (for that computer!)
%    sub_file_name = where is data going?
%
% Example:
%   [sub_ID sub_num sub_file_name] = RDCL_CreateSubjectID(exp_name, 'ComputerName', 'TRUE');
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
% 2/9/07 alc Added sub_num to return
% 4/21/13 alc added defaults/options for data prefix/folder/postfix

% Constants
RDCL_Constants

% Defaults
computer_name = '';
data_folder = 'Data/';
data_file_label = '_data_';
data_file_postfix = '.dat';

% Create prefixes
out_file_prefix = [data_folder exp_name data_file_label];   % What to call the output file

k = 1;
while k < length(varargin)

    if strcmp('ComputerName', varargin{k})

      computer_name = RDCL_GetComputerName;

    elseif strcmp('DataFolder', varargin{k})

			tmp_data_folder = varargin{k + 1};

			if (exist(tmp_data_folder, 'dir') == 7)
				data_folder = tmp_data_folder
			elseif (!isstr(tmp_data_folder))
				error('Data folder argument must be a string in RDCL_CreateSubjectID');
			elseif (!mkdir(tmp_data_folder))
					error('Cannot make data folder in RDCL_CreateSubjectID');
			end

			k = k + 1;

   elseif strcmp('DataFileLabel', varargin{k})

			tmp_data_file_label = varargin{k + 1};

			if (!isstr(tmp_data_file_label))
				error('Data file label argument must be a string in RDCL_CreateSubjectID');
			end

			k = k + 1;

   elseif strcmp('DataFilePostfix', varargin{k})

			tmp_data_file_postfix = varargin{k + 1};

			if (!isstr(tmp_data_file_postfix))
				error('Data file postfix argument must be a string in RDCL_CreateSubjectID');
			end

			k = k + 1;

   else

		error('Incorrect arguments for RDCL_CreateSubjectID');

   end

	k = k + 1;

end

done = 0;
sub_num = 1;
sub_ID = [computer_name '_' num2str(sub_num)];
while ~done
    sub_file_name = [out_file_prefix sub_ID data_file_postfix];
    if exist(sub_file_name) == 0
        done = 1;
    else
        sub_num = sub_num + 1;
        sub_ID = [computer_name '_' num2str(sub_num)];
    end
end
