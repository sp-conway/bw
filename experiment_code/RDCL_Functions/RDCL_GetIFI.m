function ifi_info = RDCL_GetIFI(varargin)

% function RDCL_GetIFI
%
% Measures the flip interval to be used for timing.

% Perform extra calibration pass to estimate monitor refresh
% interval. We want at least 100 valid samples, requiring standard
% deviation of the measurements below 50 microseconds, but timing out
% after 20 seconds if we can't get that level of accuracy:
%
% Arguments:
%    varargin:
%      'debug' = 'T' or 'F'. F skips tests.
%
% Returns:
%   ifi_info.
%				duration = time between flips (s).
%       nvalid = number of valid samples.
%       stddev = standard deviation of samples (s).
%
% Example:
%   ifi_info = RDCL_GetIFI();
%
% Andrew Cohen
% 5/9/06
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Known Bugs:
%   none
%
% Change History:
%
%   5/15/06: Added skip_sampling option; calculates max_priority, not sent
%     in as argument.
%   5/23/06: To prevent update from being called all the time,
%     changed priority back to old_priority at end, not 0.
%   12/6/06 - Michael Ross - corrected error message
%   7/25/11 - alc - simplified to one call to GetFlipInterval
%   7/26/11 - alc - Added wait message.
%		4/21/13 - alc Changed text characteristics

% Constants
global screen_ptr

% Defaults
debug = false;

% Get user options
if nargin >= 2
    for k = 1:2:length(varargin)

        option_text = varargin{k};

        if strcmp(option_text, 'debug')
            debug = varargin{k+1}(1);
        else
            error('RDCL_GetIFI: Incorrect option');
        end

    end
end

if debug == 'T'

  % Change text style
  RDCL_FontSetup('TextFont', 'Times', 'TextSize', 32, 'TextStyle', 0);

  % Waiting message
  break_message = 'SKIPPING measuring monitor refresh interval ...';
  DrawFormattedText(screen_ptr, break_message, 'center', 'center', [255 255 255]);
  Screen('Flip', screen_ptr);

  WaitSecs(.1);

  [ifi_info.duration ifi_info.nvalid ifi_info.stdev ] = ...
	  Screen('GetFlipInterval', screen_ptr, 1, 0.0005, 2);

  % Clear the screen
  Screen('Flip', screen_ptr);

else

  % Store the original priority
  old_priority = Priority(MaxPriority(screen_ptr));

  % Change text style
  RDCL_FontSetup('TextFont', 'Times', 'TextSize', 32, 'TextStyle', 0);

  % Waiting message
  break_message = 'Measuring monitor refresh interval (~20s)...';
  DrawFormattedText(screen_ptr, break_message, 'center', 'center', [255 255 255]);
  Screen('Flip', screen_ptr);

  WaitSecs(2);

  [ifi_info.duration ifi_info.nvalid ifi_info.stdev ] = ...
	  Screen('GetFlipInterval', screen_ptr, 10, 0.0005, 2);
  %[ifi_info.duration ifi_info.nvalid ifi_info.stdev ] = ...
  %	Screen('GetFlipInterval', screen_ptr, 100, 0.00005, 20);

  % Clear the screen
  Screen('Flip', screen_ptr);

  % Go back to the original priority
  Priority(old_priority);

end
