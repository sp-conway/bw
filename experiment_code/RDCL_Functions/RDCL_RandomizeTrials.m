function trial_order = RDCL_RandomizeTrials(num_trials)

% function RDCL_RandomizeTrails
%
% Returns a vector of integers 1:num_trials with the numbers in random
%   order.
%
% Arguments:
%   num_trials = the number of trials.
%
% Returns:
%    trial_order = the random vector or trials.
%
% Example:
%   trial_order = RDCL_RandomizeTrials(20);
%
% Andrew Cohen
% 5/8/06
% /* Copyright (c) 2019 Andrew L. Cohen 
%
% Known Bugs:
%   none
%
% Change History:
%

trial_order = [1:num_trials]';
tmp = rand(num_trials, 1);
trial_order = [tmp trial_order];
trial_order = sortrows(trial_order);
trial_order = trial_order(:, 2);
