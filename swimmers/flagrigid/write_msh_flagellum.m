function [] = write_msh_flagellum(N, L, rad,step_length,height,rad_section,section_type,ke_type, filename) 

disp("Parameters for the flagellum mesh : ")
% Gypsilab path
run('addpathGypsilab.m')

mesh_tail=mshHelix_head_tail_separated(N, L, rad,step_length,height,rad_section,section_type,ke_type);
disp("Mesh written in file : " + filename+"mesh_tail.ply");
mshWritePly(filename+"mesh_tail.ply",mesh_tail);



end