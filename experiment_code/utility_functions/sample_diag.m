function stim = sample_diag(x, intercept)
    w = shuffle(x)(1);
    h = intercept - w;
    stim = [w h];
end
