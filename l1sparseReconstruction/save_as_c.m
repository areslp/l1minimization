function [ ] = save_as_c(file)
data=load(file);

[pathstr,name,ext] = fileparts(file);

ff=fopen([name '_c' ext],'w');

for i=1:size(data,1)
    for j=1:size(data,2)
        fprintf(ff,'%f ',data(i,j));
    end
    fprintf(ff,'\n');
end

fclose(ff);
end

