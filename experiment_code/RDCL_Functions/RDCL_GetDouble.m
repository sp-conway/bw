function the_integer = RDCL_GetDouble(varargin)

% function RDCL_GetDouble
%
% Double input.
%
% Arguments:
%   varargin:
%     'Location' = [x y] here the input is displayed. Default is center of the screen. Left edge of text.
%     'Font' = the font name.
%			'FontSize' = the font size in pixels.
%			'Color' = the text color in RGB.
%			'InputDevice' = the number of the device that will get button press.
%			'MaxLength' = the maximum number of input characters.
%
% Returns:
% 		the_double
%
% Example:
%   the_double = RDCL_GetDouble();
%
% Andrew Cohen
% 5/10/13
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Notes:
%
% Known Bugs:
%   none
%
% Change History:
%

global screen_ptr

ints_n_per = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '.'};

done = 0;
while(~done)

	the_string = RDCL_GetString(varargin{:}, 'LegalChars', ints_n_per);

	the_double = str2num(the_string);

	if(~isempty(the_double))
		done = 1;
	end

end
