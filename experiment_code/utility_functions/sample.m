function samp = sample(v,n)
    s = shuffle(v);
    samp = s(1:n);
end 
