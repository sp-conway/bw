function [sub_info file_info] = RDCL_GetSubInfo(sub_info, file_info, varargin)

% function RDCL_GetSubInfo
%
% Uses a dialog box to input subject information.
%
% Arguments:
%   varargin:
%     'IDPrefix' = if want to add a prefix to the subject ID (str)
%     'OtherInputs' = prompts for information beyond sub num (cell of strings)
%       Takes 2 arguments, the first is the prompt, the 2nd is the type of
%       input, 'str'=string or 'int'=int. Prompts must be valid
%       variable names.
%			'UseJava' = Use a java pop-up, requires java package. (default = 0).
%
% Returns:
% 		sub_info.
%			sub_num (as input)
%			sub_ID ([prefix_S_#])
%			other prompt info by names given in OtherInputs
%       file_info.
%			data_file_name_mat (data_folder/sub_ID.mat)
%			data_file_name_dat (data_folder/sub_ID.dat)
%
% Example:
%   sub_info = RDCL_GetSubInfo('IDPrefix', exp_name, 'OtherInputs', {'cond', 'gender'});
%
% Andrew Cohen
% 7/25/11
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Notes:
%   The java package from http://octave.sourceforge.net/ must be installed.
%		<In synaptic package manager install octave3.2-headers, default-jdk, openjdk-6-jdk>
%		export JAVA_HOME=/usr/lib/jvm/java-6-openjdk
%		export PATH=$PATH:$JAVA_HOME
%		<assuming amd64 architecture>
%		export LD_LIBRARY_PATH=/usr/lib/jvm/java-1.6.0-openjdk/jre/lib/amd64
%		cd /usr/lib/jvm/java-6-openjdk/jre/lib/amd64/
%		<if needed>
%		mkdir client
%		sudo ln -s ../server/libjvm.so .
%		<Download java pkg from sourceforge (e.g., java-1.2.8.tar.gz)>
%		<start octave as root>
%		sudo octave
%		<in octave>
%		SETENV("JAVA_HOME", "/usr/lib/jvm/java-6-openjdk");
%		<install java pkg, uninstall first if already there (check with pkg list) using pkg uninstall java>
%		pkg install -verbose -auto java-1.2.8.tar.gz
%		(You may need to update the versions/architecture used.)
%		(Instructions for Ubuntu 11.10, octave3.2.4, java 1.2.8, jdk 6)
%   Assumes output file is a .mat file.
%   Only accepts numeric subject numbers.
%
%		If you use java, put this call before you open the screen.
%   If you don't use java, you must open the screen first.
%
%   Sub_ID 0 is only for debugging purposes.
%
%   The above may be out of date.
%   !!!java version not tested!!!
%
% Known Bugs:
%   none
%
% Change History:
% 5/17/19 alc - Fixed other inputs non-java version.
% 2/5/18 alc - if Sub_ID is 0, overwrites old file
% 5/10/13 alc - added use java option
% 4/15/14 alc - added .mat and .dat output file options
%

global screen_ptr;

% Initialize variables
add_prefix = 0;
other_inputs = {};
other_inputs_types = {};
use_java = 0;

% Handle the optional arguments
in_ind = 1;
k = 1;
while k < length(varargin)

  if strcmp('IDPrefix', varargin{k})
		add_prefix = 1;
		prefix = varargin{k + 1};
		k = k + 1;
	elseif strcmp('OtherInputs', varargin{k})
		other_inputs{in_ind} = char(varargin{k + 1});
		other_inputs_types{in_ind} = varargin{k + 2};
    in_ind = in_ind + 1;
		k = k + 2;
	elseif strcmp('UseJava', varargin{k})
		use_java = varargin{k + 1};
		k = k + 1;
  else
    error('RDCL_GetSubInfo: Incorrect option');
  end

	k = k + 1;

end

% Java startup and information
if (use_java)

  % Concatenate the input names to get the prompts
  prompts = [{'Subject Number'}, other_inputs];

    % The dialog box title
  title = 'Subject Information';

  % No default information
  defaults = cell(length(prompts), 1);

	pkg load java

end

% Get the information
% Input a subject number and see if it is unused
done = 0;
while ~done

	if (use_java) % JAVA

		% Input dialog box
		input_vals = inputdlg(prompts, title, 1, defaults);

		% Exit if hit cancel
		if isempty(input_vals)
			error('RDCL_GetSubInfo: Hit cancel');
		end

		% Subject number
		sub_num_str = input_vals{1};
		sub_info.sub_num = str2num(sub_num_str);

		% Input a bad number, retry
		if isempty(sub_info.sub_num)
			defaults = input_vals;
			defaults{1} = [];
		  continue;
		end

	else % NOT JAVA

		% Get subject number
		sub_info.sub_num = RDCL_GetInteger('Prompt', 'Subject Number: ');

		sub_num_str = num2str(sub_info.sub_num);

	end

	% Create the subject ID
	sub_info.sub_ID = ['S_' sub_num_str];

	if add_prefix
		sub_info.sub_ID = [prefix '_' sub_info.sub_ID];
	end

	file_info.data_file_name_mat = [file_info.data_folder sub_info.sub_ID '.mat'];
	file_info.data_file_name_dat = [file_info.data_folder sub_info.sub_ID '.dat'];

	% For eyetracking (has to be short)
	file_info.edf_file_name = [file_info.data_folder 'S_' sub_num_str '.edf'];

  % Only exit if the subject number is not being used
  if ( (exist(file_info.data_file_name_mat) == 0  && exist(file_info.data_file_name_dat) == 0) || ...
       sub_num_str == '0' )
		done = 1;
  end

end

% Get other information if not using java
if (!use_java)

  for i = 1:length(other_inputs)

    WaitSecs(.25);

    if strcmp(other_inputs_types{i}, 'int')

      input_val = RDCL_GetInteger('Prompt', [other_inputs{i} ': ']);
    	eval(sprintf("sub_info.%s = %d;", other_inputs{i}, input_val));

    elseif strcmp(other_inputs_types{i}, 'str')

      input_val = RDCL_GetString('Prompt', [other_inputs{i} ': ']);
    	eval(sprintf("sub_info.%s = '%s';", other_inputs{i}, input_val));

    endif

  end

end

% Close up java
if (use_java)

  % Other information
  for i = 1:length(other_inputs)

    eval(sprintf("sub_info.%s = '%s'", other_inputs{i}, input_vals{i + 1}));

  end

  % Unload the java package
	pkg unload java

end
