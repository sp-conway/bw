function [mood_info] = RDCL_MoodCheck(screen_rect, data_folder, sub_ID, sub_cond, exp_name)

% function RDCL_MoodCheck
%
% Administers a scale to check the effectiveness of a mood induction.
%
% Arguments:
%   screen_rect = screen dimensions for setting text location
%   data_folder = where to store the data (same as main experiment)
%   sub_ID = subject ID (same as main experiment)
%   sub_cond = subject condition (same as main experiment)
%   exp_name = name of main experiment
%
% Returns:
%   mood_info = struct including the time, trial number, mood item, and response for each trial
%
% Example:
%   mood_info = RDCL_MoodCheck(S.rect, file_info.data_folder, sub_info.sub_ID, sub_info.sub_cond, exp_info.name);
%
% Andrea Cataldo
% 3/19/14
% /* Copyright (c) 2019 Andrew L. Cohen and Andrea Cataldo */
%
% Known Bugs:
%
% Change History:

global screen_ptr

% Define basic colors
white = WhiteIndex(screen_ptr);
black = BlackIndex(screen_ptr);
gray = RDCL_GetGray();

% Constants

mood_check_txt = {'HAPPY', 'ANGRY', 'ANXIOUS', 'CALM', ...
      'IRRITATED', 'SAD', 'CONTENT', 'WORRIED', 'AFRAID', ...
      'GOOD', 'ANNOYED', 'NERVOUS', 'FEARFUL', 'RELAXED', ...
      'GLOOMY', 'MISERABLE', 'NEUTRAL', 'SCARED', 'EAGER', 'EXCITED'};

text_start_y = 100;
text_incr = 35;

% Open output files

mood_info.data_file_name = [data_folder 'mood_' sub_ID '.mat']; % Stores to same data file as main experiment

out_file_ptr = fopen(mood_info.data_file_name, 'w');

header_text{1} = 'RDCL_MoodCheck: checks the effectiveness of a mood induction';
header_text{2} = 'exp_name: name of main experiment';
header_text{3} = 'cond: determined by main experiment';
header_text{4} = 'mood_check_text: mood word presented in each trial (1-20)';
header_text{5} = 'resp_num: 1=not at all, 2, 3=a little, 4, 5=very much';
header_text{6} = 'sub_ID exp_name cond date time trial mood_check_txt resp_num';

% Clear the screen to start
Screen('FillRect', screen_ptr, gray);
[feedback_vbl feedback_sos feedback_timestamp] = Screen('Flip', screen_ptr);

% Display instructions
for i = 1:9

  switch(i)
    case 1,
       tmp_text = ' ';
    case 2,
      tmp_text = 'Please indicate how you felt while you were recalling and writing about your past';
    case 3,
      tmp_text = 'memory. For each item, please choose a number that best corresponds to your';
    case 4,
      tmp_text = 'feelings at the moment that you were thinking about your autobiographical memory.';
    case 5,
      tmp_text = ' ';
    case 6,
      tmp_text = 'After each response, you will automatically move on to the next item.';
    case 7,
      tmp_text = 'There will be 20 items total.';
    case 8,
      tmp_text = ' ';
    case 9,
      tmp_text = 'When you are ready, press SPACEBAR to begin.';
  end

  [norm_bounds offset_bounds] = Screen('TextBounds', ...
    screen_ptr, tmp_text);
  Screen('TextFont', screen_ptr, 'Arial');
  Screen('TextSize', screen_ptr, 32);
  Screen('DrawText', screen_ptr, tmp_text, ...
    round(screen_rect(3)/2 - norm_bounds(3)/2), ...
    text_start_y + text_incr*i, black);
  [stim_vbl stim_sos stim_timestamp] = Screen('Flip', screen_ptr, feedback_vbl + .750, 1);

end

RDCL_WaitForKeyPress(1);

Screen('FillRect', screen_ptr, gray);
[feedback_vbl feedback_sos feedback_timestamp] = Screen('Flip', screen_ptr);

% Run trials

for mood_question = 1:length(mood_check_txt);

  % Display scale items
  for i = 1:9

    switch(i)
      case 1,
        tmp_text = ' ';
      case 2,
        tmp_text = ['# ' num2str(mood_question)];
      case 3,
        tmp_text = ' ';
      case 4,
        tmp_text = ' ';
      case 5,
        tmp_text = mood_check_txt{mood_question};
      case 6,
        tmp_text = ' ';
      case 7,
        tmp_text = ' ';
      case 8,
        tmp_text = 'Not at all            A little            Very much';
      case 9,
        tmp_text = '1.......2.......3.......4.......5';
    end

    [norm_bounds offset_bounds] = Screen('TextBounds', ...
      screen_ptr, tmp_text);
    Screen('TextFont', screen_ptr, 'Arial');
    Screen('TextSize', screen_ptr, 32);
    Screen('DrawText', screen_ptr, tmp_text, ...
      round(screen_rect(3)/2 - norm_bounds(3)/2), ...
      round(screen_rect(4)/4) + text_incr*i, black);
    [stim_vbl stim_sos stim_timestamp] = Screen('Flip', screen_ptr, feedback_vbl + .750, 1);

  end

  % Get the response
  [resp_key resp_num time_secs] = ...
    RDCL_GetResponse({{'1!', 1}, {'2@', 2}, {'3#', 3}, ...
    {'4$', 4}, {'5%', 5}});

  [feedback_vbl feedback_sos feedback_timestamp] = Screen('Flip', screen_ptr);

  % Store data in a struct
  mood_info.m_time(mood_question, :) = clock;
  mood_info.mood_question(mood_question) = mood_question;
  mood_info.mood_check_txt = mood_check_txt;
  mood_info.m_resp_num(mood_question) = resp_num;

  % Store data to a file
  time = clock;
  fprintf(out_file_ptr, '%s %s %d %s %s %d %s %d \n', ...
      sub_ID, exp_name, sub_cond, date, [num2str(time(4)) ':' num2str(time(5)) ':' num2str(round(time(6)))], ...
      mood_question, mood_check_txt{mood_question}, resp_num);

  % Increment mood trial
  mood_question = mood_question + 1;

end
