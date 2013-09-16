function [ ] = draw_points_and_normals( points,normals,t )
figure;
p_num = size(points,1);
title(t);
axis equal;
hold on;
% plot3(points(:,1),points(:,2),points(:,3));
scatter3(points(:,1),points(:,2),points(:,3),3,[1,0,0]);
for i=1:p_num
    point = points(i,:);
    p_end = point + 0.1*normals(i,:);
    mt = [point;p_end];
    arrow3(point,p_end,'r-1',0);
end
hold off;
axis equal;
end
