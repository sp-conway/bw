function noise = add_noise(max_noise)
    dir = sample([-1 1],1);
    noise = dir*(rand()*max_noise);
end 
