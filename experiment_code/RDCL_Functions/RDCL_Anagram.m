
function [rt timed_out] = RDCL_Anagram(varargin)

% function RDCL_Anagram
%
% Complete an anagram.
%
% Arguments:
%    varargin:
%      'Zipf' = [lo hi] the frequency range in Zipf's measure, Default = [3.5 4.5]
%      'WPM' = [lo hi] the frequency range in words per million, Default = [3.5 4.5]
%      'WordLength' = [lo hi] range of word length in letters, Default = [4 5]
%      'BGCol' = [r g b] bg color,  Default = gray
%      'StrCol' = [r g b] string color,  Default = white
%      'MaxTime' = scalar (s). Max amount of time to finish, Default = no limit. Rough estimate.
%      'TimedOutPrompt' = str text to show if timed out
%      'TimedOutPromptTime' = scalar (s). How long to show the prompt
%      'DB' = string. The database to use. Default = 'RDCL_Anagram_DB_NR.csv';
%
% Uses database in RDCL_Anagram_DB_NR.csv (NR stands for no-repeats)
%  From https://www.ugent.be/pp/experimentele-psychologie/en/research/documents/subtlexus
%  Excel file with 60,384 words that have a frequency higher than 1 per million.
%  Zipf = LOG(SUBTLwf + 1) + 3, SUBTLwf = per million (Adrian Staub's suggestion)
%  Zipf frequency:  low < 3.5, 3.5 <= med < 4.5, 4.5 > high
%  anagram_change.m has been used to remove all words from the database which:
%   1) Start with a capital letter.
%   2) Can be unscrambled to make a different word (no-repeats).
%
% If both WPM and Zipf set, WPM is used.
%
% Returns:
%    rt = amount of time to solve (not precise, uses tic toc)
%    timed_out = whether they finished in time
%
% Example:
%   RDCL_Anagram('Frequency', [3.5 4.5], 'Word_Length', [4 5]);
%
% Andrew Cohen
% 6/13/18
% /* Copyright (c) 2019 Andrew L. Cohen */
%
% Known Bugs:
%   none
%
% Change History:

global screen_ptr

% Default options
zipf_range = [3.5 4.5]; % Low < 3.5; 3.5 <= Med < 4.5; High >= 4.5 (Zipf)
wpm_range = [exp(zipf_range(1)-3)-1 exp(zipf_range(1)-3)-1];
wpm_set = 0;
word_len = [4 4];
bg_col = RDCL_GetGray();
str_col = [255 255 255];
max_time = Inf;
timed_out_prompt = '';
timed_out_prompt_time = 2;
db = 'RDCL_Resources/RDCL_Anagram_DB_NR.csv';

% Get user options
if nargin >= 2
    for k = 1:2:length(varargin)

        option_text = varargin{k};
        if strcmp(option_text, 'Zipf')
            zipf_range = varargin{k+1};
        elseif strcmp(option_text,'WPM')
            wpm_range = varargin{k+1};
            wpm_set = 1;
        elseif strcmp(option_text,'WordLength')
            word_len = varargin{k+1};
        elseif strcmp(option_text,'BGCol')
            bg_col = varargin{k+1};
        elseif strcmp(option_text,'StrCol')
            str_col = varargin{k+1};
        elseif strcmp(option_text,'MaxTime')
            max_time = varargin{k+1}(1);
        elseif strcmp(option_text,'TimedOutPrompt')
            timed_out_prompt = varargin{k+1};
        elseif strcmp(option_text,'TimedOutPromptTime')
            timed_out_prompt_time = varargin{k+1}(1);
        elseif strcmp(option_text,'DB')
            db = varargin{k+1};
        else
            error('RDCL_Anagram: Incorrect option');
        end

    end
end

% WPM to Zipf
% Exclusion based on Zipf
if (wpm_set)
  zipf_range = [log(wpm_range(1)+1)+3 log(wpm_range(2)+1)+3];
end

% Read in the database
[wd tmp tmp] = fileparts(mfilename('fullpath'));
db_name = [wd '/' db];
[words,wpm,zipf,lens] = textread(db_name, '%s %f %f %d', 'delimiter', ' ');

% Select word of right length and frequency
keep_len = find(lens>=word_len(1) & lens<=word_len(2) == 1);
keep_zipf = find(zipf>=zipf_range(1) & zipf<=zipf_range(2) == 1);
keep = intersect(keep_len, keep_zipf);
if(length(keep)==0) error('RDCL_Anagram: No words selected.'); end
selected_word_num = keep(randi(size(keep,1),1));
the_word = tolower(words{selected_word_num});

% Scramble it
ps = perms(1:length(the_word));
r = randi(length(ps)-1); % Don't allow unscrambled words
p = ps(r,:);
word_scram = the_word(p);

% Show it and get response
done = false;
tic_id = tic;
timed_out = false;
while (~done)

  WaitSecs(.2);

  % Get their response
  the_prompt = ['Unscramble the letters "' word_scram '":  '];
  input_str = RDCL_GetString('Prompt', the_prompt, 'Color', str_col, 'MaxTime', max_time);

  % Is  the input string in the dictionary?
  IndexC = strfind(words, input_str);
  Index = find(not(cellfun('isempty', IndexC)));

  if (length(input_str)==length(the_word) && ...
    ~isempty(Index) && ...
    strcmp(sort(the_word),sort(input_str)))
    done = true;
  end

  if (toc(tic_id) > max_time)
    done = true;
    timed_out = true;
  end

end

rt = toc(tic_id);

if (timed_out)
  [mid_x mid_y] = RDCL_GetScreenMiddle();
  tmp_text = timed_out_prompt;
  [norm_bounds offset_bounds] = Screen('TextBounds', screen_ptr, tmp_text);
  Screen('DrawText', screen_ptr, tmp_text, ...
    mid_x - round(norm_bounds(3)/2), mid_y - (round(norm_bounds(4)/2)),
    str_col);
  Screen('Flip', screen_ptr);
  WaitSecs(timed_out_prompt_time);
end

clear IndexC Index words wpm zipf lens keep keep_len keep_zipf ps
