function [] = write_mesh(points,normals,file)
f=fopen(file,'w');
n=size(points,1);
for i=1:n
    fprintf(f,'%f %f %f %f %f %f\n',points(i,1),points(i,2),points(i,3),normals(i,1),normals(i,2),normals(i,3));
end
fclose(f);
