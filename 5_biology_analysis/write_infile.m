function write_infile(filename,source_matrix )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
fid=fopen(filename,'w');
[x,y]=size(source_matrix);

for i=1:x
    for j=1:y-1
        fprintf(fid,'%f\t',source_matrix(i,j));
    end
    fprintf(fid,'%f\n',source_matrix(i,y));
end
fclose(fid);


end

