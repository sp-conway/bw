function the_string = RDCL_GetString(varargin)

	% function RDCL_GetString
	%
	% String input.
	%
	% Arguments:
	%   varargin:
	%     'Location' = [x y] here the input is displayed. Default is center of the screen. Left edge of text.
	%     'Font' = the font name.
	%			'FontSize' = the font size in pixels.
	%			'Color' = the text color in RGB.
	%			'InputDevice' = the number of the device that will get button press.
	%			'MaxLength' = the maximum number of input characters.
	% 		'LegalChars' = a cell array of the legal keys.
	%			'Prompt' = Text to show before the input.
	%			'MaxTime' = Max allowed time (in s). Default = no limit. Rough est of time.
	%
	% Returns:
	% 		the_string
	%
	% Example:
	%   the_string = RDCL_GetString();
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
	%  5/17/19 alc major change to make text input faster.
	%  6/13/18 Moved default x position so right edge of prompt is at x middle.
	%  6/13/18 Added a maximum time option
	%  6/13/18 Reduced keypress time

	global screen_ptr

	KbName('UnifyKeyNames')

	% Defaults
	[text_loc_x text_loc_y] = RDCL_GetScreenMiddle();
	text_font = 'Arial';
	text_font_size = 32;
	text_color = BlackIndex(screen_ptr);
	max_length = 50;
	prompt = '';
	input_device = -1;
	resps = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', ...
	'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', ...
	'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'space', '_'};
	max_time = Inf;

	% Handle the optional arguments
	k = 1;
	while k < length(varargin)

		if strcmp('Location', varargin{k})
			text_loc = varargin{k + 1};
			text_loc_x = text_loc(1);
			text_loc_y = text_loc(2);
			k = k + 1;
		elseif strcmp('Font', varargin{k})
			text_font = varargin{k + 1};
			k = k + 1;
		elseif strcmp('FontSize', varargin{k})
			text_font_size = varargin{k + 1};
			k = k + 1;
		elseif strcmp('Color', varargin{k})
			text_color = varargin{k + 1};
			k = k + 1;
		elseif strcmp('MaxLength', varargin{k})
			max_length = varargin{k + 1};
			k = k + 1;
		elseif strcmp('InputDevice', varargin{k})
			input_device = varargin{k + 1};
			k = k + 1;
		elseif strcmp('LegalChars', varargin{k})
			resps = varargin{k + 1};
			k = k + 1;
		elseif strcmp('Prompt', varargin{k})
			prompt = varargin{k + 1};
			k = k + 1;
		elseif strcmp('MaxTime', varargin{k})
			max_time = varargin{k + 1};
			k = k + 1;
		else
			error('RDCL_GetString: Incorrect option');
		end

		k = k + 1;

	end

	% Prompt to left of middle
	if (!isnull(prompt))
		[norm_bounds offset_bounds] = Screen('TextBounds', screen_ptr, prompt);
		text_loc_x = text_loc_x -   norm_bounds(3);
	end

	% Set options
	[old_font_name old_font_number] = Screen('TextFont', screen_ptr, text_font);
	old_text_size = Screen('TextSize', screen_ptr, text_font_size);
	old_text_style = Screen('TextStyle', screen_ptr, 0);

	% Constants
	cursor = '|';

	% Initialize the text
	the_text = '';

	% Clear the screen
	gray = RDCL_GetGray();
	Screen('FillRect', screen_ptr, gray);
	Screen('DrawText', screen_ptr, [prompt the_text cursor], text_loc_x, text_loc_y, text_color);
	Screen('Flip', screen_ptr);

	% Get and display the input
	done = 0;
	% prev_key_down = 0;
	% key_down = 0;
	% only_shift_down = 0;
	prev_char = '';
	accept_repeats = true;
	tic_id = tic;
	while ~done

		if input_device > 0
			[key_down time_secs key_code] = KbCheck(input_device);
		else
			[key_down time_secs key_code] = KbCheck();
		end

		the_char = KbName(key_code);

		if key_down == 0 || strcmp(the_char, 'LeftShift') || strcmp(the_char, 'RightShift')
			accept_repeats = true;
		endif

		if sum(key_code) == 2 && ismember('LeftShift', the_char)
			use_caps = true;
			the_char = char(setdiff(the_char, 'LeftShift'));
		elseif sum(key_code) == 2 && ismember('RightShift', the_char)
			use_caps = true;
			the_char = char(setdiff(the_char, 'RightShift'));
		elseif sum(key_code) > 1
			the_char = '';
		else
			use_caps = false;
		endif

		switch the_char
			case '1!' the_char = '1';
			case '2@' the_char = '2';
			case '3#' the_char = '3';
			case '4$' the_char = '4';
			case '5%' the_char = '5';
			case '6^' the_char = '6';
			case '7&' the_char = '7';
			case '8*' the_char = '8';
			case '9(' the_char = '9';
			case '0)' the_char = '0';
		end

		if strcmp(prev_char, the_char) && ~accept_repeats
			the_char = '';
		endif

		%if strcmp(the_char, 'BackSpace') && strcmp(prev_char, the_char)
	%		WaitSecs(.1);
%		endif

		if strcmp(the_char, 'Return') || strcmp(the_char, 'Enter')
			done = true;
		elseif strcmp(the_char, 'BackSpace')
			if ~isempty(the_text)
				the_text = the_text(1:length(the_text)-1);
				WaitSecs(.1);
			endif
			prev_char = the_char;
		else
			if length(the_text) < max_length && ismember(tolower(the_char), resps)
				if strcmp(the_char, 'space')
					the_text = [the_text ' '];
				else
					if use_caps
						the_text = [the_text toupper(the_char)];
					else
						the_text = [the_text the_char];
					endif
				endif
				prev_char = the_char;
				accept_repeats = false;
			endif
		end

		text_output = [prompt the_text cursor];
		Screen('DrawText', screen_ptr, text_output, text_loc_x, text_loc_y, text_color);
		Screen('Flip', screen_ptr);

		if (toc(tic_id) > max_time)
			done = true;
		 	the_text = '';
		end

	end

	% Return value
	the_string = the_text;

	% Reset fonts
	Screen('TextFont', screen_ptr, old_font_name);
	Screen('TextSize', screen_ptr, old_text_size);
	Screen('TextStyle', screen_ptr, old_text_style);

end
