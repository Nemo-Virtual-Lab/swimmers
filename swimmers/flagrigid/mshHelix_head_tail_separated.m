function [mesh,length_centerline,V] = mshHelix_head_tail_separated(N,L,rad,step_length,height,rad_section,section_type,ke_type,varargin)
%+========================================================================+
%|                                                                        |
%|                 OPENMSH - LIBRARY FOR MESH MANAGEMENT                  |
%|           openMsh is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : mshHelix_head_tail_separated.m                |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Luca Berti                                    |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Build mesh for a special helix                |
%|  `---'  |                                                              |
%+========================================================================+


tail_length_is_ok = true;
accuracy = 1e-3;
desired_length = L;
%desired_length = 3;
%desired_length = 8; % Phan-Thien


while tail_length_is_ok
    t = linspace(0,height,N);
    %k_E=0.2*2*pi/step_length;% Thick tail simulations
    if strcmp(ke_type, 'schum')
        k_E=0.333*2*pi/step_length;% Schum
    elseif strcmp(ke_type, 'pt')
        k_E=2*pi/step_length; % Phan-Thien
    end

    x=t;
    y=rad*(1-exp(-k_E.^2*t.^2)).*cos(2*pi/step_length*t);
    z=rad*(1-exp(-k_E.^2*t.^2)).*sin(2*pi/step_length*t);
    length_centerline=0;
    %disp(size(t,2))
    %disp(x(2)-x(1))
    
    for i=1:size(t,2)-1
        length_centerline=length_centerline+norm([x(i+1)-x(i);
                                                  y(i+1)-y(i);
                                                  z(i+1)-z(i)]);
                                                                                           
                                              %if length_centerline <=desired_length + accuracy && length_centerline >=desired_length-accuracy
                                                  %length_centerline
                                              %   tail_length_is_ok = false;
                                              %    disp('ok length')
                                              %   break
                                              %elseif length_centerline >=desired_length + accuracy
                                               %   length_centerline
                                               %   N = N+5;
                                                  %disp('too long')
                                               %   break
                                              %end
                                              %disp(length_centerline)
                                              %disp(height)
                                            
                                              %disp(N)
                                             
    end

    if length_centerline <=desired_length-accuracy
       %   length_centerline;
          height = height + 0.0001; 
          %disp('too short')
    elseif length_centerline >= desired_length + accuracy
	  height = height - 0.0001;
    else
	  tail_length_is_ok = false;
    end

end

% t=t(1:i-1);
% x=x(1:i-1);
% y=y(1:i-1);
% z=z(1:i-1);
%t=t(1:i+1);
%x=x(1:i+1);
%y=y(1:i+1);
%z=z(1:i+1);

length_centerline=0;
for i=1:size(t,2)-1
    length_centerline=length_centerline+norm([x(i+1)-x(i);
                                              y(i+1)-y(i);
                                              z(i+1)-z(i)]);
end
length_centerline
N=length(t)
xm=x(end)
Nlamb = xm/step_length
% Points describing the section -> they should be given as an input, and
% the first point should be the one attached to the center of the helix

Section_area = pi*rad_section^2; % Compute the area section to calculate the iso-area tails with increasing perimeter
if section_type == "circular"
    section = mshDisk(10,rad_section);
    boundary=mshBoundary(section);
elseif section_type == "square"
    side_square = sqrt(Section_area);
    section = mshSquare(20,[side_square; side_square]);
    boundary = mshBoundary(section);
elseif section_type == "rectangle"
    side_rectangle1 = sqrt(Section_area)*3/2 ;
    side_rectangle2 = sqrt(Section_area)*2/3 ;
    section = mshSquare(20,[side_rectangle1; side_rectangle2]);
    boundary = mshBoundary(section);
elseif section_type == "triangle"
    side_triangle = sqrt(Section_area/sqrt(3))*2;
    top_vertex = [0 rad_section 0];
    left_vertex = [-cos(-pi/6)*rad_section sin(-pi/6)*rad_section 0];
    right_vertex = [cos(-pi/6)*rad_section sin(-pi/6)*rad_section 0];
    origin = [0,0,0];
    X = [top_vertex;left_vertex;right_vertex;origin];
    DT = delaunayTriangulation(X(:,1),X(:,2));
    Points = [DT.Points zeros(size(X(:,1),1),1)];
    section = msh(Points,DT.ConnectivityList);
    boundary = mshBoundary(section);
else
    disp('Error on the type of section');
end

% Construction of the local reference frame - The parallel transported
% frame was chosen with respect to the Frenet-Serret 
% Ref. for the algorithm https://janakiev.com/blog/framing-parametric-curves/
tangent=[ones(1,N);
    -2*pi/step_length*rad*sin(2*pi/step_length*t).*(1-exp(-k_E.^2.*t.^2))+cos(2.*pi.*t/step_length).*rad.*2.*t.*k_E.^2.*exp(-k_E.^2.*t.^2);
    2.*pi/step_length.*rad.*cos(2.*pi/step_length.*t).*(1-exp(-k_E.^2.*t.^2))+sin(2.*pi.*t/step_length).*rad.*2.*t.*k_E.^2.*exp(-k_E.^2.*t.^2)];
%tangent=tangent./norm(tangent);
normal=[zeros(1,N);
    -2.*pi/step_length.*rad.*(2.*pi/step_length.*cos(2.*pi.*t/step_length).*(1-exp(-k_E.^2.*t.^2))+sin(2.*pi/step_length.*t).*2.*t.*k_E.^2.*exp(-k_E.^2.*t.^2))+...
    rad.*(2.*k_E.^2.*exp(-k_E.^2.*t.^2).*cos(2.*pi.*t/step_length)-2.*pi/step_length.*sin(2.*pi/step_length).*2.*t.*k_E.^2.*exp(-k_E.^2.*t.^2)-cos(2.*pi.*t/step_length).*(2.*t.*k_E.^2).^2.*exp(-k_E.^2.*t.^2));
    2.*pi/step_length.*rad.*(-2.*pi/step_length.*sin(2.*pi.*t/step_length).*(1-exp(-k_E.^2.*t.^2))+cos(2.*pi.*t/step_length).*2.*t.*k_E.^2.*exp(-k_E.^2.*t.^2))+...
    rad.*(2.*pi/step_length.*cos(2.*pi/step_length.*t).*2.*t.*k_E.^2.*exp(-k_E.^2.*t.^2)+sin(2.*pi.*t/step_length).*2.*t.*k_E.^2.*exp(-k_E.^2.*t.^2)-sin(2.*pi.*t/step_length).*(2.*t.*k_E.^2).^2.*exp(-k_E.^2.*t.^2))
    ];
%normal=normal./norm(normal);
for i=1:N
    normal(:,i)=normal(:,i)/norm(normal(:,i));
    tangent(:,i)=tangent(:,i)/norm(tangent(:,i));
end
R = zeros(3,3);
for i =1:N-1
    bin = cross(tangent(:,i),tangent(:,i+1));
    if norm(bin) < 1e-8
        normal(:,i+1) = normal(:,i);
    else
        bin = bin/norm(bin);
        theta = acos(dot(tangent(:,i),tangent(:,i+1)));
        % ROtation matrix of angle theta around vector bin; it allows
        % parallel transport of the normal vector normal
        R(1,1) = cos(theta) + bin(1)^2*(1-cos(theta));
        R(2,2) = cos(theta) + bin(2)^2*(1-cos(theta));
        R(3,3) = cos(theta) + bin(3)^2*(1-cos(theta));
        R(1,2) = bin(1)*bin(2)*(1-cos(theta))-bin(3)*sin(theta);
        R(2,1) = bin(1)*bin(2)*(1-cos(theta))+bin(3)*sin(theta);
        R(1,3) = bin(1)*bin(3)*(1-cos(theta))+bin(2)*sin(theta);
        R(3,1) = bin(1)*bin(3)*(1-cos(theta))-bin(2)*sin(theta);
        R(2,3) = bin(3)*bin(2)*(1-cos(theta))-bin(1)*sin(theta);
        R(3,2) = bin(3)*bin(2)*(1-cos(theta))+bin(1)*sin(theta);
        normal(:,i+1) = R*normal(:,i);
    end
end

binormal=cross(tangent,normal);
boundary_helix=([x(1) y(1) z(1)]'*ones(size(section.vtx(:,1),1),1)'+(section.vtx(:,1)*normal(:,1)')'+(section.vtx(:,2)*binormal(:,1)')')';
%boundary.vtx
for i=2:size(t,2)-1
    boundary_helix=[boundary_helix;
                    ([x(i) y(i) z(i)]'*ones(size(boundary.vtx(:,1),1),1)'+(boundary.vtx(:,1)*normal(:,i)')'+(boundary.vtx(:,2)*binormal(:,i)')')'];
end
 boundary_helix=[boundary_helix;
     ([x(end) y(end) z(end)]'*ones(size(section.vtx(:,1),1),1)'+(section.vtx(:,1)*normal(:,end)')'+(section.vtx(:,2)*binormal(:,end)')')'];
%plot3(boundary_helix(:,1),boundary_helix(:,2),boundary_helix(:,3),'-')

% Assemble the vertices matrix + point in the middle
X=[x' y' z';
    boundary_helix];
t=MyRobustCrust(X);
if section_type == "circular"
    if isempty(varargin) == 0 
        [V,S] = alphavol(X,cell2mat(varargin(1)),1);
        %S = alphaShape(X,cell2math(varargin(1)));
        disp('varargin is non empty');
    else
        [V,S] = alphavol(X,rad_section*2.5,1);
    end
    [V,S] = alphavol(X,rad_section*log(L),1); 
    %S = alphaShape(X, rad_section*10);
    %[V,S] = alphavol(X,rad_section*0.75,1);
elseif section_type == "square"
    %[V,S] = alphavol(X,side_square*0.5,1);
     %S = alphaShape(X, side_square*10);
elseif section_type == "rectangle"
    %[V,S] = alphavol(X,side_rectangle1,1);
     %S = alphaShape(X, side_rectangle1*10);
elseif section_type == "triangle"
    %[V,S] = alphavol(X,side_triangle,1);
    %[V,S] = alphavol(X,side_triangle*0.4,1);
    %S = alphaShape(X, side_triangle*10);
end
% Delaunay triangulation
DT        = delaunay(X);
%[elt,vtx] = freeBoundary(DT);

% Mesh
mesh = msh(X,S.tri);
%mesh = msh(X,alphaTriangulation(S));
mesh=mesh.bnd;
end














% 
% binormal=cross(tangent,normal);
% boundary_helix=([]')';%[x(1) y(1) z(1)]'*ones(size(section.vtx(:,1),1),1)'+(section.vtx(:,1)*normal(:,1)')'+(section.vtx(:,2)*binormal(:,1)')')';
% %boundary.vtx
% 
% 
% nb_test = 100;
% rad_ratio = 1.01;
% 
% occ = rad_ratio^(nb_test);
% for i=2:2+nb_test
%     section_start = mshDisk(10,rad_section/occ);
%     boundary_start=mshBoundary(section_start);   
%     boundary_helix=[boundary_helix;
%                     ([x(i) y(i) z(i)]'*ones(size(boundary_start.vtx(:,1),1),1)'+(boundary_start.vtx(:,1)*normal(:,i)')'+(boundary_start.vtx(:,2)*binormal(:,i)')')'];
%     occ = occ/rad_ratio;
% 
% end
% 
% 
% for i=2+nb_test+1:size(t,2)-1-nb_test+1
%     boundary_helix=[boundary_helix;
%                     ([x(i) y(i) z(i)]'*ones(size(boundary.vtx(:,1),1),1)'+(boundary.vtx(:,1)*normal(:,i)')'+(boundary.vtx(:,2)*binormal(:,i)')')'];
% end
% 
% 
% occ = rad_ratio;
% for i=size(t,2)-nb_test+1:size(t,2)-1
% 
% 
%     section_end = mshDisk(10,rad_section/occ);
%     boundary_end=mshBoundary(section_end);   
%     boundary_helix=[boundary_helix;
%                     ([x(i) y(i) z(i)]'*ones(size(boundary_end.vtx(:,1),1),1)'+(boundary_end.vtx(:,1)*normal(:,i)')'+(boundary_end.vtx(:,2)*binormal(:,i)')')'];
%    occ = occ*rad_ratio;
% 
% end


 %boundary_helix=[boundary_helix;
 %    ([x(end) y(end) z(end)]'*ones(size(section.vtx(:,1),1),1)'+(section.vtx(:,1)*normal(:,end)')'+(section.vtx(:,2)*binormal(:,end)')')'];
%plot3(boundary_helix(:,1),boundary_helix(:,2),boundary_helix(:,3),'-')

% Assemble the vertices matrix + point in the middle
% X=[x' y' z';
%     boundary_helix];
% t=MyRobustCrust(X);
% if section_type == "circular"
%     if isempty(varargin) == 0 
%         [V,S] = alphavol(X,cell2mat(varargin(1)),1);
%         %S = alphaShape(X,cell2math(varargin(1)));
%         disp('varargin is non empty');
%     else
%         [V,S] = alphavol(X,rad_section*2.5,1);
%     end
%     [V,S] = alphavol(X,rad_section*log(L)*2,1); 
%     %S = alphaShape(X, rad_section*10);
%     %[V,S] = alphavol(X,rad_section*0.75,1);
% elseif section_type == "square"
%     %[V,S] = alphavol(X,side_square*0.5,1);
%      %S = alphaShape(X, side_square*10);
% elseif section_type == "rectangle"
%     %[V,S] = alphavol(X,side_rectangle1,1);
%      %S = alphaShape(X, side_rectangle1*10);
% elseif section_type == "triangle"
%     %[V,S] = alphavol(X,side_triangle,1);
%     %[V,S] = alphavol(X,side_triangle*0.4,1);
%     %S = alphaShape(X, side_triangle*10);
% end
% % Delaunay triangulation
% DT        = delaunay(X);
% 
% 
% % Mesh
% mesh = msh(X,S.tri);
% %mesh = msh(X,alphaTriangulation(S));
% mesh=mesh.bnd;
% end
