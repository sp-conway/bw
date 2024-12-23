function par = read_params(file,format)
    f=fopen(file,"r");
    d=textscan(f,"%s",'Delimiter','\n');
    par=d{:}
    fclose(f);
end 