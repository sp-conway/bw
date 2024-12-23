function stim = sample_filler(min_x, max_x)
    x = min_x:max_x;
    y = min_x:max_x;
    x = shuffle(x)(1);
    y = shuffle(y)(2);
    stim = [x y];
end
