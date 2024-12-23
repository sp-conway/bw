function RDCL_Beep(varargin)

% function RDCL_Beep
%
% Displays a text file.
%
% Arguments:
%    varargin:
%      'Vol' = Volume.
%      'Freq' = Frequency (Hz).
%      'Dur' = Duration (ms).
%      'Chan' = Which single channel? All others silent. Default = all, -1.
%      'Dev' = Which device? Defult = first with output channels, -1.
%			 'NumChan' = Hardcode the number of channels for the device, -1.
%
% Returns:
%    none
%
% Example:
%   RDCL_Beep('Vol', .2, 'Freq', 500, 'Dur', 1000);
%
% Andrew Cohen
% 7/27/11
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Notes:
% 	Uses the default sampling rate.
%   Closes all audio devices at start, none at end! Do it manually with PsychPortAudio('Close');
%
% Known Bugs:
%   none
%
% Change History:
% 4/22/13 - alc - Added NumChan option

% Defaults
vol = .2;
freq = 200;
dur = 500;
chan = -1; % -1 = All, which output channel to use
dev = -1; % -1 = First
num_chan = -1; % -1 = Use all

% Get user options
if nargin >= 2
    for k = 1:2:length(varargin)

        option_text = varargin{k};

        if strcmp(option_text, 'Vol')
            vol = varargin{k+1}(1);
        elseif strcmp(option_text,'Freq')
            freq = varargin{k+1}(1);
        elseif strcmp(option_text, 'Dur')
            dur = varargin{k+1}(1);
        elseif strcmp(option_text, 'Chan')
            chan = varargin{k+1}(1);
        elseif strcmp(option_text, 'Dev')
            dev = varargin{k+1}(1);
        elseif strcmp(option_text, 'NumChan')
            num_chan = varargin{k+1}(1);
        else
            error('RDCL_Beep: Incorrect option');
        end

    end
end

% Get device information
devices = PsychPortAudio('GetDevices');

% Find the first output device, if needed
if dev == -1

	% For later indexing
	dev_vector_index = min(find([devices.NrOutputChannels] > 0));

	dev_index = [devices.DeviceIndex](dev_vector_index);

else

	if dev > max([devices.DeviceIndex]) || dev < 0

		error('RDCL_Beep: Device index out of bounds.');

	else

		% For later indexing
		dev_vector_index = dev + 1;

		dev_index = dev;

	end

end

% Find How many output channels the device has
if num_chan == -1

	num_chan = [devices.NrOutputChannels](dev_vector_index);

end

if chan ~= -1 && (chan < 1 || chan > num_chan)

	error('RDCL_Beep: Channel out of bounds.');

end

% Determine the beep
beep_single_chan = vol*MakeBeep(freq, dur/1000);

% Fill in the channels
if chan == -1

	beep = repmat(beep_single_chan, [num_chan 1]);

else

	for i = 1:num_chan

		if i == chan

			beep(i, :) = beep_single_chan;

		else

			beep(i, :) = zeros(size(beep_single_chan));

		end

	end

end

% Closes all audio devices
PsychPortAudio('Close');

% Open the audo device
pa_handle = PsychPortAudio('Open', dev_index);

% Fill the audio playback buffer:
PsychPortAudio('FillBuffer', pa_handle, beep);

% Start audio playback
PsychPortAudio('Start', pa_handle);
