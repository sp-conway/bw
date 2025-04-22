% Fal  2024
% bw choice CE exp
% Sean Conway
function bw_exp
   
    % PTB SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    more off;
    global screen_ptr;
    debug=false;
    if debug
      PsychDebugWindowConfiguration;
    endif

    % ADD NECESSARY FILE PATHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % FUNCTIONS
    addpath('./RDCL_Functions');    
    addpath('./utility_functions');
    
    % TRIAL PARAMS
    addpath('./choice_trial_params');

    % DATA PATHS
    choice_destination_folder = './data';
    addpath(choice_destination_folder);
    
    computer_num = 2;

    % MORE SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %KB NAME
    KbName('UnifyKeyNames');
    
    if ~debug
      % Dummy calls
      RDCL_DummyCalls();
    endif
    
    % seed rng
    RDCL_SeedRandomNumbers();
      
    % EXPERIMENT SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % NUMBER OF BLOCKS 
    n_choice_blocks = 8;

    % number of practice trials 
    n_choice_prac = 3;

    %% max vertical jitter
    vjit = 15; 

    % COLORS
    white = [255 255 255];
    black = [0 0 0];
    red = [255 179 179];
  
    % rectangle distance
    rect_dist = .05;
    
    % rectangle fills
    rect_fill_avail=[178 190 181];
    rect_fill_unavail=red;
    
    % minimum and maximum stimulus value (for filler trials)
    min_stim = 56;
    max_stim = 195;
    
    % mouse lag
    lag_time = .008;

    %%%%%%% GET SUBJECT NUMBER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Prompt text
    prompt = ('Please enter the participant number:');
    
    % Get Experimenter input
    answer = inputdlg(prompt);
    sub_n = deal(answer{:});
    sub_n = num2str(sub_n);

    % OUTPUT FILE NAMES
    choice_output_name = ['bw_' sub_n '.csv'];      

    % Check to see if subject number has been used
    while exist(choice_output_name)
      sub_n_error = ('That participant number has already been used. Please double check the runsheet and enter a new number.');
      error_message = inputdlg(sub_n_error);
      sub_n = deal(error_message{:});
      sub_n = num2str(sub_n);
      choice_output_name = ['bw_' sub_n '.csv'];      
    end  

    % IMPORTANT %
    % % ASSIGN SUBJECT TO CONDITION
    % odd subject # means best/worst, even subject # worst/best
    if rem(str2double(sub_n),2)>0
      bw_cond = "bw";
    else
      bw_cond = "wb";
    endif 

    % figure out text prompt order based on experimental condition
    if(bw_cond=="bw")
      text_prompt1 = "Click on the rectangle with the LARGEST area";
      text_prompt2 = "Click on the rectangle with the SMALLEST area";
    elseif(bw_cond=="wb")
      text_prompt1 = "Click on the rectangle with the SMALLEST area";
      text_prompt2 = "Click on the rectangle with the LARGEST area";
    endif 
    
    % START EXPERIMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Open Screen  
    [screen_ptr s_rect] = Screen('OpenWindow', 0, 0, [], 32, 2);
    
    % FIND SCREEN CENTER
    s_middle_x = s_rect(3)*.5;
    s_middle_y = s_rect(4)*.5;

    % rectangle distance in pixels
    rect_dist_px = round(s_middle_x * rect_dist);
    % GET CONSENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rect border thickness
    submit_rect_thick = 1;
    
    if ~debug
      % SHOW CONSENT FORM
      cons_info = imfinfo('consent.png');
      cons_width = 2000;
      cons_height = 1600;
      %Get the image for instructions
      cons= imread('consent.png');
      cons_scale = .4;
      cons_rect = [s_middle_x-round(cons_scale*cons_width/1.5) 
                    s_middle_y-round(cons_scale*cons_height/1.5) 
                    s_middle_x+round(cons_scale*cons_width/1.5) 
                    s_middle_y+round(cons_scale*cons_height/1.5)];
      
      % Put the image into the buffer      
      agree_rect = CenterRectOnPoint([0 0 100 60], s_middle_x-300, cons_rect(4)+75);
      not_agree_rect = CenterRectOnPoint([0 0 300 60], s_middle_x+300, cons_rect(4)+75);      
        
      %Put the image into the buffer
      Screen('FillRect', screen_ptr, white);
      Screen('PutImage', screen_ptr, cons, cons_rect);
      Screen('FrameRect', screen_ptr, black, agree_rect, submit_rect_thick);
      Screen('FrameRect', screen_ptr, black, not_agree_rect, submit_rect_thick);
      RDCL_DrawText(s_middle_x-300, cons_rect(4)+65, "I agree", 'Col',black);
      RDCL_DrawText(s_middle_x+300, cons_rect(4)+65, "I do not agree", 'Col',black);
      Screen('Flip', screen_ptr); 

      % CHECK IF THEY ACCEPTED OR DECLINED TO PARTICIPATE
      deciding=true;
      while deciding
      [clicks click_x click_y] = GetClicks(screen_ptr);
        if click_x > agree_rect(1) && click_x < agree_rect(3) && click_y > agree_rect(2) && click_y < agree_rect(4)
          deciding=false;
          declined=false;
        elseif  click_x > not_agree_rect(1) && click_x < not_agree_rect(3) && click_y > not_agree_rect(2) && click_y < not_agree_rect(4)
          deciding=false;
          declined=true;
        endif
      endwhile  
     if declined
       sca
     endif     
    
     % Clear Screen
     Screen('FillRect', screen_ptr, white);
     Screen('Flip', screen_ptr); 
     WaitSecs(1);
    endif
     
    % Data file
    fp_choice = fopen(fullfile(choice_destination_folder,choice_output_name),'w');     

    % Write to output file
    fprintf(fp_choice, 'sub_n, computer_n, bw_cond, effect, set, distance, diag, block_n, trial_n, h1, w1, h2, w2, h3, w3, choice, rt \n');

    % READ IN CHOICE TRIAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % NUMERIC PARAMS
    h1 = dlmread("choice_trial_params/h1.txt");
    h2 = dlmread("choice_trial_params/h2.txt");
    h3 = dlmread("choice_trial_params/h3.txt");
    w1 = dlmread("choice_trial_params/w1.txt");
    w2 = dlmread("choice_trial_params/w2.txt");
    w3 = dlmread("choice_trial_params/w3.txt");
    distance = dlmread("choice_trial_params/distance.txt");
    diag = dlmread("choice_trial_params/diag.txt");

    % STRING PARAMS
    sets = read_params("choice_trial_params/set.txt","%s");
    effect = read_params("choice_trial_params/effect.txt","%s");

    % diagonal ranges
    d3_r = 120:195;
    d2_r = 90:165;
    d1_r = 60:135;

    % diagonal intercepts
    d1_int = 195;
    d2_int = 255;
    d3_int = 315;

    % minimum and maximum stimulus value (for filler trials)
    min_stim = 56;
    max_stim = 195;

    % DETERMINE N TRIALS
    n_choice_trials=numel(h1);

    % quadrant locs for drawing rectangles
    rect_quadrant = [s_middle_x-(s_middle_x/2) s_middle_y-(s_middle_y/2) s_middle_x+(s_middle_x/2) s_middle_y+(s_middle_y/2)];

    % actual "rectangle" where rectangles go
    rectangle_box_c = [s_middle_x s_middle_y];
    rectangle_box = CenterRectOnPoint(rect_quadrant, s_middle_x, s_middle_y);

    % width & height of the rectangle box
    rectangle_box_w = rectangle_box(3)-rectangle_box(1);
    rectangle_box_h = rectangle_box(4)-rectangle_box(2);

    % CHOICE INSTRUCITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Screen('FillRect', screen_ptr, white);
    RDCL_DrawText(s_middle_x, s_middle_y-150, "Welcome to the experiment!", 'Col', black);
    RDCL_DrawText(s_middle_x, s_middle_y, "Press 'space' to continue the instructions.", 'Col', black);
    Screen('Flip', screen_ptr); 
    WaitSecs(1);
    RDCL_GetResponse({{'space', 0}});
    
    
    Screen('FillRect', screen_ptr, white);
    RDCL_DrawText(s_middle_x, s_middle_y-300, "In this experiment, you will be making decisions about rectangles.", 'Col', black);
    RDCL_DrawText(s_middle_x, s_middle_y-150, "Press 'space' to continue the instructions.", 'Col', black);
    Screen('Flip', screen_ptr); 
    WaitSecs(1);
    RDCL_GetResponse({{'space', 0}});

    Screen('FillRect', screen_ptr, white);
    RDCL_DrawText(s_middle_x, s_middle_y-300, "On each trial, you will see 3 rectangles on the screen, labeled ‘1’, ‘2’, and ‘3’.  ", 'Col', black);
    RDCL_DrawText(s_middle_x, s_middle_y-250, "You will select both the largest and the smallest rectangle on each trial.", 'Col', black);
    RDCL_DrawText(s_middle_x, s_middle_y-200, "You select a rectangle using the left mouse button", 'Col', black);
    RDCL_DrawText(s_middle_x, s_middle_y-100, "Press ‘space’ to continue the instructions.", 'Col', black);
    Screen('Flip', screen_ptr); 
    WaitSecs(1);
    RDCL_GetResponse({{'space', 0}});
    
    if(bw_cond=="bw")
      Screen('FillRect', screen_ptr, white);
      RDCL_DrawText(s_middle_x, s_middle_y-300, "First, you will select the largest rectangle, after which it will turn pink.", 'Col', black);
      RDCL_DrawText(s_middle_x, s_middle_y-250, "After you make your first selection, you cannot go back.", 'Col', black);
      RDCL_DrawText(s_middle_x, s_middle_y-200, "From the remaining two rectangles, you will select the smallest rectangle.", 'Col', black);
      RDCL_DrawText(s_middle_x, s_middle_y-100, "Press ‘space’ to continue the instructions.", 'Col', black);
      Screen('Flip', screen_ptr); 
      WaitSecs(1);
      RDCL_GetResponse({{'space', 0}});
    elseif(bw_cond=="wb")
      Screen('FillRect', screen_ptr, white);
      RDCL_DrawText(s_middle_x, s_middle_y-300, "First, you will select the smallest rectangle, after which it will turn pink.", 'Col', black);
      RDCL_DrawText(s_middle_x, s_middle_y-250, "After you make your first selection, you cannot go back.", 'Col', black);
      RDCL_DrawText(s_middle_x, s_middle_y-200, "From the remaining two rectangles, you will select the largest rectangle.", 'Col', black);
      RDCL_DrawText(s_middle_x, s_middle_y-100, "Press ‘space’ to continue the instructions.", 'Col', black);
      Screen('Flip', screen_ptr); 
      WaitSecs(1);
      RDCL_GetResponse({{'space', 0}});
    endif  
    

    Screen('FillRect', screen_ptr, white);
    RDCL_DrawText(s_middle_x, s_middle_y-300, "Some trials will be harder than others. We only ask that you try your best.", 'Col', black);
    RDCL_DrawText(s_middle_x, s_middle_y-250, "We will keep track of how you did and let you know at the end of the experiment.", 'Col', black);
    RDCL_DrawText(s_middle_x, s_middle_y-150, "Press ‘space’ to continue the instructions.", 'Col', black);
    Screen('Flip', screen_ptr); 
    WaitSecs(1);
    RDCL_GetResponse({{'space', 0}});

    Screen('FillRect', screen_ptr, white);
    RDCL_DrawText(s_middle_x, s_middle_y-300, "We will start with a few practice trials.", 'Col', black);
    RDCL_DrawText(s_middle_x, s_middle_y-200, "Press ‘space’ to start the practice trials.", 'Col', black);
    Screen('Flip', screen_ptr); 
    WaitSecs(1);
    RDCL_GetResponse({{'space', 0}});

    % PRACTICE TRIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for trial_index = 1:n_choice_prac
      set_tmp = "NA";
      distance_tmp = 999;
      effect_tmp = "practice"
      diag_tmp = 999;
      block_num=0;

      r1_tmp = sample_filler(min_stim,max_stim);
      r2_tmp = sample_filler(min_stim,max_stim);
      r3_tmp = sample_filler(min_stim,max_stim);
      w_all_tmp_unshuffled = [r1_tmp(1) r2_tmp(1) r3_tmp(1)];
      h_all_tmp_unshuffled = [r1_tmp(2) r2_tmp(2) r3_tmp(2)];

      % shuffle rects 
      rects_order_tmp = RDCL_RandomizeTrials(3);

      % assigned shuffled rects to h1, h2, h3
      h1_tmp = round(h_all_tmp_unshuffled(rects_order_tmp(1)));
      h2_tmp = round(h_all_tmp_unshuffled(rects_order_tmp(2)));
      h3_tmp = round(h_all_tmp_unshuffled(rects_order_tmp(3)));

      % assigned shuffled rects to w1, w2, w3
      w1_tmp = round(w_all_tmp_unshuffled(rects_order_tmp(1)));
      w2_tmp = round(w_all_tmp_unshuffled(rects_order_tmp(2)));
      w3_tmp = round(w_all_tmp_unshuffled(rects_order_tmp(3)));

      % array of all heights
      h_all_tmp_shuffled = [h1_tmp h2_tmp h3_tmp];

      % array of all widths
      w_all_tmp_shuffled = [w1_tmp w2_tmp w3_tmp];

      % Rectangle 1 center x pos
      r1_x = rectangle_box_c(1)-(.5*w1_tmp)-(rect_dist_px)-(.5*w2_tmp);
      
      % Rectangle 2 center x pos
      r2_x = rectangle_box_c(1);

      % rectangle 3 center x pos
      r3_x = rectangle_box_c(1)+(.5*w2_tmp)+(rect_dist_px)+(.5*w3_tmp);
      
      % jitter
      r1_jit = add_noise(vjit);
      r2_jit = add_noise(vjit);
      r3_jit = add_noise(vjit);
      
      % Rectangle 1 center y pos
      r1_y = round(rectangle_box_c(2)+(.5*rect_dist_px)+.5*h1_tmp+r1_jit);

      % Rectangle 2 center y pos
      r2_y = round(rectangle_box_c(2)-(.5*rect_dist_px)-.5*h2_tmp+r2_jit);

      % Rectangle 3 center y pos
      r3_y = round(rectangle_box_c(2)+(.5*rect_dist_px)+.5*h3_tmp+r3_jit); 
                     
      % RECRTANGLE LABELS - y
      rtxt_y = max([round(r1_y+h1_tmp)
                    round(r2_y+h2_tmp)
                    round(r3_y+h3_tmp)])+round(.05*rectangle_box_h);
      
      % Rectangle 1 text center
      r1_txt_x = r1_x;
      r1_txt_y = rtxt_y;

      % Rectangle 2 text center
      r2_txt_x = r2_x;
      r2_txt_y = rtxt_y;

      % Rectangle 3 text center
      r3_txt_x = r3_x;
      r3_txt_y = rtxt_y;
      
      % Center rectangles 
      r1_rect = CenterRectOnPoint([0 0 w1_tmp h1_tmp], r1_x, r1_y);
      r2_rect = CenterRectOnPoint([0 0 w2_tmp h2_tmp], r2_x, r2_y);
      r3_rect = CenterRectOnPoint([0 0 w3_tmp h3_tmp], r3_x, r3_y);

      % flip screen to white
      Screen('FillRect', screen_ptr, white);
      Screen('Flip', screen_ptr); 
      WaitSecs(.1);

      % draw rectangles
      Screen('FillRect', screen_ptr, rect_fill_avail, r1_rect);
      Screen('FillRect', screen_ptr, rect_fill_avail, r2_rect);
      Screen('FillRect', screen_ptr, rect_fill_avail, r3_rect);

      % draw labels
      RDCL_DrawText(r1_txt_x, r1_txt_y, "1", 'Col', black);
      RDCL_DrawText(r2_txt_x, r2_txt_y, "2", 'Col', black);
      RDCL_DrawText(r3_txt_x, r3_txt_y, "3", 'Col', black);

      % Draw prompt
      RDCL_DrawText(s_middle_x, s_rect(4)*.25, text_prompt1, 'Col', black);

      % Position mouse
      SetMouse(s_middle_x, rectangle_box_c(2));

      % FLIP SCREEN SHOW EVERYTHING
      % SHOW EVERYTHING
      [trial_onset_time] = Screen('Flip', screen_ptr);

      % GET FIRST CHOICE - BEST if BW, WORST if WB
      choosing1=true;
      while choosing1
        WaitSecs(lag_time);
      [clicks click_x click_y] = GetClicks(screen_ptr);
        if click_x > r1_rect(1) && click_x <= r1_rect(3) && click_y > r1_rect(2) && click_y<=r1_rect(4)
          choice_tmp1 = 1;
          rt_tmp = GetSecs()-trial_onset_time;        
          choosing1=false;
        elseif click_x > r2_rect(1) && click_x <= r2_rect(3) && click_y > r2_rect(2) && click_y<=r2_rect(4)
          choice_tmp1 = 2;
          rt_tmp = GetSecs()-trial_onset_time;
          choosing1=false;
        elseif click_x > r3_rect(1) && click_x <= r3_rect(3) && click_y > r3_rect(2) && click_y<=r3_rect(4)
          choice_tmp1 = 3;
          rt_tmp = GetSecs()-trial_onset_time;
          choosing1=false;
        endif
      endwhile

      % WRITE DATA
      fprintf(fp_choice, 
      '%s, %d, %s, %s, %s, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',
      sub_n, computer_num, bw_cond, effect_tmp, set_tmp, distance_tmp, diag_tmp, block_num, trial_index, 
      h1_tmp, w1_tmp, h2_tmp, w2_tmp, h3_tmp, w3_tmp, choice_tmp1, rt_tmp);
      fflush(fp_choice);
      
      % FIGURE OUT NEW COLORS
      if(choice_tmp1==1)
        rect1_fill=rect_fill_unavail;
        rect2_fill=rect_fill_avail;
        rect3_fill=rect_fill_avail;
      elseif(choice_tmp1==2)
        rect1_fill=rect_fill_avail;
        rect2_fill=rect_fill_unavail;
        rect3_fill=rect_fill_avail;
      elseif(choice_tmp1==3)
        rect1_fill=rect_fill_avail;
        rect2_fill=rect_fill_avail;
        rect3_fill=rect_fill_unavail;
      endif
      
      % NOW GET SECOND CHOICE
      % draw rectangles
      Screen('FillRect', screen_ptr, rect1_fill, r1_rect);
      Screen('FillRect', screen_ptr, rect2_fill, r2_rect);
      Screen('FillRect', screen_ptr, rect3_fill, r3_rect);

      % draw labels
      RDCL_DrawText(r1_txt_x, r1_txt_y, "1", 'Col', black);
      RDCL_DrawText(r2_txt_x, r2_txt_y, "2", 'Col', black);
      RDCL_DrawText(r3_txt_x, r3_txt_y, "3", 'Col', black);

      % PROMPT - WORST IF BW, BEST IF WB
      RDCL_DrawText(s_middle_x, s_rect(4)*.25, text_prompt2, 'Col', black);
      [trial_onset_time] = Screen('Flip', screen_ptr);

      % GET SECOND CHOICE
      choosing2=true;
      while choosing2
          WaitSecs(lag_time);
      [clicks click_x click_y] = GetClicks(screen_ptr);
        if click_x > r1_rect(1) && click_x <= r1_rect(3) && click_y > r1_rect(2) && click_y<=r1_rect(4) && choice_tmp1 != 1
          choice_tmp2 = 1;
          rt_tmp = GetSecs()-trial_onset_time;        
          choosing2=false;
        elseif click_x > r2_rect(1) && click_x <= r2_rect(3) && click_y > r2_rect(2) && click_y<=r2_rect(4) && choice_tmp1 != 2
          choice_tmp2 = 2;
          rt_tmp = GetSecs()-trial_onset_time;
          choosing2=false;
        elseif click_x > r3_rect(1) && click_x <= r3_rect(3) && click_y > r3_rect(2) && click_y<=r3_rect(4) && choice_tmp1 != 3
          choice_tmp2 = 3;
          rt_tmp = GetSecs()-trial_onset_time;
          choosing2=false;
        endif
      endwhile

      % WRITE DATA
      fprintf(fp_choice, 
      '%s, %d, %s, %s, %s, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',
      sub_n, computer_num, bw_cond, effect_tmp, set_tmp, distance_tmp, diag_tmp, block_num, trial_index, 
      h1_tmp, w1_tmp, h2_tmp, w2_tmp, h3_tmp, w3_tmp, choice_tmp2, rt_tmp);
      fflush(fp_choice);
    endfor

    % FINISHED PRACTICE TRIALS
    % FOLLOW-UP INSTRUCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Screen('FillRect', screen_ptr, white);
    RDCL_DrawText(s_middle_x, s_middle_y-300, "You finished the practice trials.", 'Col', black);
    RDCL_DrawText(s_middle_x, s_middle_y-200, "Press space to continue.", 'Col', black);
    Screen('Flip', screen_ptr); 
    WaitSecs(1);
    RDCL_GetResponse({{'space', 0}});
    
    Screen('FillRect', screen_ptr, white);
    RDCL_DrawText(s_middle_x, s_middle_y-300, "The experimental trials will start next.", 'Col', black);
    RDCL_DrawText(s_middle_x, s_middle_y-200, "There will be a break halfway through.", 'Col', black);
    RDCL_DrawText(s_middle_x, s_middle_y-100, "Press space to continue to the experimental trials.", 'Col', black);
    Screen('Flip', screen_ptr); 
    WaitSecs(1);
    RDCL_GetResponse({{'space', 0}});

    % INIT VARIABLE TO KEEP TRACK OF HOW MANY FIRST CHOICES CORRECT
    choice1_correct = 0;

    % INIT VARIABLE TO KEEP TRACK OF HOW MANY SECOND CHOICES CORRECT
    choice2_correct = 0;

    % FIGURE OUT HOW MANY TOTAL CHOICE TRIALS in EXP
    choice_total = n_choice_trials*n_choice_blocks;
    

    % Loop through blocks
    for block_num = 1:n_choice_blocks

      % Give break halfway through
      if  ( (block_num -1) / n_choice_blocks ) == .5
        Screen('FillRect', screen_ptr, white);
        RDCL_DrawText(s_middle_x, s_middle_y-300, "You are halfway through the experiment!", 'Col', black);
        RDCL_DrawText(s_middle_x, s_middle_y-200, "Take a short break.", 'Col', black);
        RDCL_DrawText(s_middle_x, s_middle_y-100, "The experiment will resume in 30 seconds.", 'Col', black);
        Screen('Flip', screen_ptr); 
        WaitSecs(30);
      endif

      % Randomize trial order in this block
      trial_order = RDCL_RandomizeTrials(n_choice_trials);

      % Loop through trials
      for trial_index = 1:n_choice_trials 
        
        % Set screen to blank first
        Screen('FillRect', screen_ptr, white);
        Screen('Flip', screen_ptr); 
        WaitSecs(.1);

        % Get shuffled trial number
        trial_num = trial_order(trial_index);

        % get set
        set_tmp = sets{trial_num};

        % get distance
        distance_tmp = distance(trial_num);

        % get effect (trial type, i.e. attraction / filler / catch)
        effect_tmp = effect{trial_num};

        % get diagonal
        diag_tmp = diag(trial_num);

        % GET RECT HEIGHT & WIDTH
        if strcmp(effect_tmp,"catch") 
          r1_tmp = sample_diag(d1_r, d1_int);
          r2_tmp = sample_diag(d2_r, d2_int);
          r3_tmp = sample_diag(d3_r, d3_int);
          w_all_tmp_unshuffled = [r1_tmp(1) r2_tmp(1) r3_tmp(1)];
          h_all_tmp_unshuffled = [r1_tmp(2) r2_tmp(2) r3_tmp(2)];
        elseif strcmp(effect_tmp,"filler")
          r1_tmp = sample_filler(min_stim,max_stim);
          r2_tmp = sample_filler(min_stim,max_stim);
          r3_tmp = sample_filler(min_stim,max_stim);
          w_all_tmp_unshuffled = [r1_tmp(1) r2_tmp(1) r3_tmp(1)];
          h_all_tmp_unshuffled = [r1_tmp(2) r2_tmp(2) r3_tmp(2)];
        elseif strcmp(effect_tmp,"attraction")
            h_all_tmp_unshuffled = [h1(trial_num) h2(trial_num) h3(trial_num)];
            w_all_tmp_unshuffled = [w1(trial_num) w2(trial_num) w3(trial_num)];
        endif

        % shuffle rects 
        rects_order_tmp = RDCL_RandomizeTrials(3);

        % assigned shuffled rects to h1, h2, h3
        h1_tmp = round(h_all_tmp_unshuffled(rects_order_tmp(1)));
        h2_tmp = round(h_all_tmp_unshuffled(rects_order_tmp(2)));
        h3_tmp = round(h_all_tmp_unshuffled(rects_order_tmp(3)));

        % assigned shuffled rects to w1, w2, w3
        w1_tmp = round(w_all_tmp_unshuffled(rects_order_tmp(1)));
        w2_tmp = round(w_all_tmp_unshuffled(rects_order_tmp(2)));
        w3_tmp = round(w_all_tmp_unshuffled(rects_order_tmp(3)));

        % array of all heights
        h_all_tmp_shuffled = [h1_tmp h2_tmp h3_tmp];

        % array of all widths
        w_all_tmp_shuffled = [w1_tmp w2_tmp w3_tmp];

        % Rectangle 1 center x pos
        r1_x = rectangle_box_c(1)-(.5*w1_tmp)-(rect_dist_px)-(.5*w2_tmp);

        % Rectangle 2 center x pos
        r2_x = rectangle_box_c(1);

        % rectangle 3 center x pos
        r3_x = rectangle_box_c(1)+(.5*w2_tmp)+(rect_dist_px)+(.5*w3_tmp);
        
        % sample noise for vertical position
        r1_jit = add_noise(vjit);
        r2_jit = add_noise(vjit);
        r3_jit = add_noise(vjit);

        % Rectangle 1 center y pos
        r1_y = round(rectangle_box_c(2)+(.5*rect_dist_px)+.5*h1_tmp+r1_jit);

        % Rectangle 2 center y pos
        r2_y = round(rectangle_box_c(2)-(.5*rect_dist_px)-.5*h2_tmp+r2_jit);

        % Rectangle 3 center y pos
        r3_y = round(rectangle_box_c(2)+(.5*rect_dist_px)+.5*h3_tmp+r3_jit);           
	

        % RECTANGLE LABELS
        rtxt_y = max([round(r1_y+h1_tmp)
                      round(r2_y+h2_tmp)
                      round(r3_y+h3_tmp)])+round(.05*rectangle_box_h);

        % Rectangle 1 text center
        r1_txt_x = r1_x;
        r1_txt_y = rtxt_y;

        % Rectangle 2 text center
        r2_txt_x = r2_x;
        r2_txt_y = rtxt_y;

        % Rectangle 3 text center
        r3_txt_x = r3_x;
        r3_txt_y = rtxt_y;
            
        % Center rectangles 
        r1_rect = CenterRectOnPoint([0 0 w1_tmp h1_tmp], r1_x, r1_y);
        r2_rect = CenterRectOnPoint([0 0 w2_tmp h2_tmp], r2_x, r2_y);
        r3_rect = CenterRectOnPoint([0 0 w3_tmp h3_tmp], r3_x, r3_y);
        
        Screen('FillRect', screen_ptr, white);
        Screen('Flip', screen_ptr); 
        WaitSecs(.1);

        % draw rectangles
        Screen('FillRect', screen_ptr, rect_fill_avail, r1_rect);
        Screen('FillRect', screen_ptr, rect_fill_avail, r2_rect);
        Screen('FillRect', screen_ptr, rect_fill_avail, r3_rect);

        % draw labels
        RDCL_DrawText(r1_txt_x, r1_txt_y, "1", 'Col', black);
        RDCL_DrawText(r2_txt_x, r2_txt_y, "2", 'Col', black);
        RDCL_DrawText(r3_txt_x, r3_txt_y, "3", 'Col', black);

        % DRAW PROMPT
        % BEST IF BW, WORST IF WB
        RDCL_DrawText(s_middle_x, s_rect(4)*.25, text_prompt1, 'Col', black);

	% Position mouse
        SetMouse(s_middle_x, round(((r1_y + r2_y)/2 + r3_y)/2) );

        % SHOW EVERYTHING
        [trial_onset_time] = Screen('Flip', screen_ptr);

        % GET CHOICE1
        choosing1=true;
        while choosing1

        [click_x click_y button] = GetMouse(screen_ptr);

        if button(1,1)
        
          if click_x > r1_rect(1) && click_x <= r1_rect(3) && click_y > r1_rect(2) && click_y<=r1_rect(4)
            rt_tmp = GetSecs()-trial_onset_time;    
	    choice_tmp1 = 1;         
            choice1_area_tmp = h1_tmp*w1_tmp;
	    choosing1=false;
          elseif click_x > r2_rect(1) && click_x <= r2_rect(3) && click_y > r2_rect(2) && click_y<=r2_rect(4)
            rt_tmp = GetSecs()-trial_onset_time;

            choice_tmp1 = 2;
            
            
            choice1_area_tmp = h2_tmp*w2_tmp;
	    choosing1=false;
          elseif click_x > r3_rect(1) && click_x <= r3_rect(3) && click_y > r3_rect(2) && click_y<=r3_rect(4)
            rt_tmp = GetSecs()-trial_onset_time;
            choice_tmp1 = 3;
            choice1_area_tmp = h3_tmp*w3_tmp;
            choosing1=false;
          endif
          WaitSecs(lag_time);
        endif
        endwhile

        % WRITE DATA
        fprintf(fp_choice, 
        '%s, %d, %s, %s, %s, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',
        sub_n, computer_num, bw_cond, effect_tmp, set_tmp, distance_tmp, diag_tmp, block_num, trial_index, 
        h1_tmp, w1_tmp, h2_tmp, w2_tmp, h3_tmp, w3_tmp, choice_tmp1, rt_tmp);
        fflush(fp_choice);
        
        % FIGURE OUT COLORS FOR CHOICE 2
        if(choice_tmp1==1)
          rect1_fill=rect_fill_unavail;
          rect2_fill=rect_fill_avail;
          rect3_fill=rect_fill_avail;
        elseif(choice_tmp1==2)
          rect1_fill=rect_fill_avail;
          rect2_fill=rect_fill_unavail;
          rect3_fill=rect_fill_avail;
        elseif(choice_tmp1==3)
          rect1_fill=rect_fill_avail;
          rect2_fill=rect_fill_avail;
          rect3_fill=rect_fill_unavail;
        endif
        
        % draw rectangles
        Screen('FillRect', screen_ptr, rect1_fill, r1_rect);
        Screen('FillRect', screen_ptr, rect2_fill, r2_rect);
        Screen('FillRect', screen_ptr, rect3_fill, r3_rect);

        % draw labels
        RDCL_DrawText(r1_txt_x, r1_txt_y, "1", 'Col', black);
        RDCL_DrawText(r2_txt_x, r2_txt_y, "2", 'Col', black);
        RDCL_DrawText(r3_txt_x, r3_txt_y, "3", 'Col', black);

        % DRAW PROMPT
        % WORST IF BW, BEST IF WB
        RDCL_DrawText(s_middle_x, s_rect(4)*.25, text_prompt2, 'Col', black);

        % SHOW EVERYTHING
        [trial_onset_time] = Screen('Flip', screen_ptr);

        % GET CHOICE 2
        choosing2=true;
        while choosing2
          [click_x click_y button] = GetMouse(screen_ptr);

          if button(1,1)
          
            if click_x > r1_rect(1) && click_x <= r1_rect(3) && click_y > r1_rect(2) && click_y<=r1_rect(4) && choice_tmp1 != 1
              rt_tmp = GetSecs()-trial_onset_time;    
              choice_tmp2 = 1;
              choice2_area_tmp = h1_tmp*w1_tmp;
              choosing2=false;    
            elseif click_x > r2_rect(1) && click_x <= r2_rect(3) && click_y > r2_rect(2) && click_y<=r2_rect(4) && choice_tmp1 != 2
              rt_tmp = GetSecs()-trial_onset_time;
        choice_tmp2 = 2;
              choice2_area_tmp = h2_tmp*w2_tmp;      	
              choosing2=false;
                  
            elseif click_x > r3_rect(1) && click_x <= r3_rect(3) && click_y > r3_rect(2) && click_y<=r3_rect(4) && choice_tmp1 != 3
              rt_tmp = GetSecs()-trial_onset_time;
        choice_tmp2 = 3;
              choice2_area_tmp = h3_tmp*w3_tmp;
              choosing2=false;
            endif
            WaitSecs(lag_time);
          endif 
        endwhile

        fprintf(fp_choice, 
        '%s, %d, %s, %s, %s, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n',
        sub_n, computer_num, bw_cond, effect_tmp, set_tmp, distance_tmp, diag_tmp, block_num, trial_index, 
        h1_tmp, w1_tmp, h2_tmp, w2_tmp, h3_tmp, w3_tmp, choice_tmp2, rt_tmp);
        fflush(fp_choice);
      
        % FIGURE  OUT WHETHER WORST AND BEST CHOICES WERE CORRECT
        if bw_cond=="bw"
          if choice1_area_tmp == max([h1_tmp*w1_tmp h2_tmp*w2_tmp h3_tmp*w3_tmp])
            choice1_correct=choice1_correct+1
          endif
          if choice2_area_tmp == min([h1_tmp*w1_tmp h2_tmp*w2_tmp h3_tmp*w3_tmp])
            choice2_correct=choice2_correct+1  
          endif   
        elseif bw_cond=="wb"
          if choice1_area_tmp == min([h1_tmp*w1_tmp h2_tmp*w2_tmp h3_tmp*w3_tmp])
            choice1_correct=choice1_correct+1
          endif
          if choice2_area_tmp == max([h1_tmp*w1_tmp h2_tmp*w2_tmp h3_tmp*w3_tmp])
            choice2_correct=choice2_correct+1  
          endif
        endif         
    endfor
    endfor
    disp(choice_total);
    
    % FIGURE OUT TOTAL PERCENTAGES CORRECT
    if bw_cond=="bw"
      bchoice_correct_perc = (choice1_correct/choice_total)*100;
      disp("best");
disp(bchoice_correct_perc);
      wchoice_correct_perc = (choice2_correct/choice_total)*100;
	disp("worst");
disp(wchoice_correct_perc);
      %disp(bchoice_correct_perc);
      %disp(wchoice_correct_perc);
    elseif bw_cond=="wb"
      bchoice_correct_perc = (choice2_correct/choice_total)*100;
      disp("best");
disp(bchoice_correct_perc);
      wchoice_correct_perc = (choice1_correct/choice_total)*100;
	disp("worst");
disp(wchoice_correct_perc);
      %disp(bchoice_correct_perc);
      %disp(wchoice_correct_perc);
    endif  

    % END OF EXPERIMENT
    Screen('FillRect', screen_ptr, white);
disp("flip");
    RDCL_DrawText(s_middle_x, s_middle_y-300, "Congratulations! You reached the end of the experiment", 'Col', black);
disp("1");
    RDCL_DrawText(s_middle_x, s_middle_y-250, "Your percentage correct for the largest rectangles was:", 'Col', black);
disp("2");
    RDCL_DrawText(s_middle_x, s_middle_y-200, [num2str(bchoice_correct_perc) "%"], 'Col', black);
disp("3");
    RDCL_DrawText(s_middle_x, s_middle_y-150, "Your percentage correct for the smallest rectangles was:", 'Col', black);
disp("4");
    RDCL_DrawText(s_middle_x, s_middle_y-100, [num2str(wchoice_correct_perc) "%"], 'Col', black);
disp("5");
    %catch 
     % RDCL_DrawText(s_middle_x, s_middle_y-200, "91% of largest choices correct", 'Col', black);
     % RDCL_DrawText(s_middle_x, s_middle_y-100, "87% of largest choices correct", 'Col', black);
   % end
    
    RDCL_DrawText(s_middle_x, s_middle_y, "Press ‘space’ to see the debriefing form.", 'Col', black);
disp("6");
    Screen('Flip', screen_ptr); 
disp("7");
    RDCL_GetResponse({{'space', 0}});
disp("8");

    % SHOW DEBRIEF
    debrief_info = imfinfo('Judg_choice_debrief.png');
    debrief_width = debrief_info.Width;
    debrief_height = debrief_info.Height;

    %Get the image for instructions
    debrief = imread('Judg_choice_debrief.png');

    %Where and how big should the image be?
    debrief_scale = .55;
    debrief_rect = [s_middle_x-round(debrief_scale*debrief_width/1.5) 
          s_middle_y-round(debrief_scale*debrief_height/1.5) 
          s_middle_x+round(debrief_scale*debrief_width/1.5) 
          s_middle_y+round(debrief_scale*debrief_height/1.5)];

    % Put the image into the buffer
    Screen('PutImage', screen_ptr, debrief, debrief_rect);

    %Open Screen with instructions displayed. When they hit space the next page shows up
    Screen('Flip', screen_ptr);
    RDCL_GetResponse({{'space', 0}}); 
    sca
end
