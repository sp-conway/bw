more off;
global screen_ptr;
addpath('..')
KbName('UnifyKeyNames')

[screen_ptr screen_rect] = Screen('OpenWindow', 0, 0, [0 0 1024 768], 32, 2);

  done = false;
  bar_vals = [25 50 25];
  prev_bar_vals = bar_vals;
  total = sum(bar_vals);
  while (~done)

    bar_vals(1) = RDCL_SliderBar(bar_vals(1), 'CenterY', 200);
    bar_vals(2) = RDCL_SliderBar(bar_vals(2), 'CenterY', 400);
    bar_vals(3) = RDCL_SliderBar(bar_vals(3), 'CenterY', 600);

    changed = find (bar_vals != prev_bar_vals);
    not_changed = find (bar_vals == prev_bar_vals);

    if (!isempty(changed))
      if (sum(bar_vals(not_changed))> 0)
        bar_vals(not_changed) = ...
          (total - bar_vals(changed))*bar_vals(not_changed)/sum(bar_vals(not_changed));
      else
        bar_vals(not_changed) = ...
          (total - bar_vals(changed))/length(not_changed);
      endif
    endif

    prev_bar_vals = bar_vals;

    Screen('Flip', screen_ptr);

  end
