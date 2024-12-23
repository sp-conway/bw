function [resp_key resp_num time_secs] = RDCL_GetResponse(legal_inputs, input_device, wait_time)

% function RDCL_GetResponse
%
% Waits for one of a particular set of keypresses.
%
% Arguments:
%   legal_inputs = a cell structure of legal keys and their associated
%     response numbers (for return value)
%   input_device = the number of the device that will get button press (optional).
%   wait_time = maximum length of time in seconds the computer will wait for a
%     response (optional).
%
%
% Returns:
%    resp_key = the key that was pressed. -1 if timed out.
%    resp_num = the user defined number associated with the pressed key. -1 if timed out.
%    time_secs = when the key was pressed. -1 if timed out.
%
% Example:
%   [resp_key resp_num time_secs] = RDCL_GetResponse({{'a', 0}, {'b', 1}, {'c'}});
%
% Andrew Cohen
% 5/8/06
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Notes:
%   You should call GetSecs() before it is used for timing, so the function is preloaded.
%     That is, put a call to GetSecs() into your program, before this function is called.
%
% Known Bugs:
%   none
%
% Change History:
%
% 5/17/19 - alc - added escape squence. only response when 1 key is hit (not > 1)
% 12/6/06 - Michael Ross - Swapped order of arguments, made input_device optional.
% 9/11/07 - Matt Zivot - Added an optional wait time.
% 9/12/07 - Andrew Cohen - Changed the way wait_time worked. Made timing more precise.
% 7/9/08 - Andrew Cohen - Fixed bug in wait time, now checks if nargin <= 2, not < 2.
% 7/27/11 - alc - Removed RDCL_Constants.

% Get the time the function started
start_time = GetSecs();

% Constants
RDCL_Constants;
timed_out = -1;

esc_keys = {'LeftControl', 'RightShift', 'LeftAlt'}; % Lctrl Rshift Lalt

num_legal_inputs = length(legal_inputs);

for i = 1:num_legal_inputs
    resps(i) = legal_inputs{i}(1);

    if length(legal_inputs{i}) == 2
        resp_num(i) = legal_inputs{i}{2};
    else
        resp_num(i) = 0; % default
    end
end

if nargin <= 2
    wait_time = Inf;
end

key_down = 0;
while ~key_down && GetSecs() - start_time < wait_time

    WaitSecs(0.001);
%    WaitSecs(0.050);

    if nargin > 1
        [key_down time_secs key_code] = KbCheck(input_device);

    else
        [key_down time_secs key_code] = KbCheck();
    end

    if key_down
	      KbReleaseWait();
        resp_key = KbName(key_code);

        % Special stop sequence
        if sum(key_code) == 3
          if strcmp(resp_key(1), esc_keys{1}) ...
                && strcmp(resp_key(2), esc_keys{2}) ...
                && strcmp(resp_key(3), esc_keys{3})
            error('Experimenter Stop!'); % Goes to catch to close things up
          else
            key_down = 0;
          endif
        elseif sum(key_code) == 1 % Make sure only 1 key being pressed
          [tf loc] = ismember(resp_key, resps);
          if tf == 1
            resp_num = resp_num(loc);
          else
            key_down = 0;
          endif
        else
          key_down = 0;
        endif

        % [tf loc] = ismember(resp_key, resps);
        % if tf == 1
        %     resp_num = resp_num(loc);
        % elseif numel(resp_key) >= 3
        %     % Check for experimenter quit
        %     if strcmp(resp_key(1), esc_keys(1)) ...
        %             && strcmp(resp_key(2), esc_keys(2)) ...
        %             && strcmp(resp_key(3), esc_keys(3))
        %         error('Experimenter Stop!'); % Goes to catch to close things up
        %     else
        %         key_down = 0;
        %     end
        % else
        %     key_down = 0;
        % end

    end
end

if key_down == 0
    resp_key = timed_out;
    resp_num = timed_out;
    time_secs = timed_out;
end
