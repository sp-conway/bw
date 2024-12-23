function [seed_rand seed_randn] = RDCL_SeedRandomNumbers()

% function RDCL_SeedRandomNumbers
%
% Seeds the random numbers
%
% Arguments:
%    none.
%
% Returns:
%    The seed for the rand and randn functions.
%
% Example:
%   [seed_rand seed_randn] = RDCL_SeedRandomNumbers();
%
% Andrew Cohen
% 5/15/06
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Known Bugs:
%   none
%
% Change History:
%   7/25/11 - alc - Added updated code for Matlab.
%   7/25/11 - alc - Separated Matlab and Octave code.
%

if exist('OCTAVE_VERSION')
	% Ocatve

	% Random number constants
    seed_rand = sum(round(sum(37*clock)));
	seed_randn = sum(round(sum(183*clock)));

	% Seed random number generators
	rand('state', seed_rand);
	randn('state', seed_randn);

else
	% Matlab

	seed_rand = sum(round(sum(37*clock)));
	seed_randn = -1;

	seed = RandStream.create('mt19937ar', 'seed', seed_rand);
	RandStream.setDefaultStream(seed);

end
