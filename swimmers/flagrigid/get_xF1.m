function [xF1_] = get_xF1(vtx,elt, alpha, beta) 

% Gypsilab path
run('addpathGypsilab.m')

% Meshes
mesh_head = msh(vtx, elt)
mean_vtx = mean(mesh_head.vtx);
x_S = [mean_vtx(1);0;0;]; %because xOy and xOZ sym.


%intersection point
if alpha ~= 0. | beta ~=0. 
  intersection_points = intersection_surf_line(x_S, mesh_head,alpha, beta);
    [maxi,argmax] = max(intersection_points(:,1));
    xF1_ = intersection_points(argmax,:); 
    
    xF2_ = [xF1_(1), xF1_(2), -xF1_(3)];
else 
    xF1_ = zeros(1,3);
    xF1_(1) = max(mesh_head.vtx(:,1));
    
end

