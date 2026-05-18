function [objective, constraints] = solve_pb(vtx,elt, N, L, rad,step_length,height,rad_section,section_type,ke_type, eps,nb_flag,alpha,gamma, beta, delta, num_pb) 

% BEM for a flagella swimmer
% Input : Mesh path of the head (meshpath_head), (N), radius of the helix (rad), length between 2 waves (step_length), (height) , radius of the section (rad_section), section type of the tail (section_type) 
% Output : cost function


% =========================== Gypsilab path ===========================
run('addpathGypsilab.m')

U0 = -0.0306;
P0 = -24.7512;
numPointsMC = 500000;


% =========================== Meshes generation ===========================
mesh_head = msh(vtx, elt);
mesh_tail=mshHelix_head_tail_separated(N, L, rad,step_length,height,rad_section,section_type, ke_type);
mean_vtx = mean(mesh_head.vtx);
x_S = [mean_vtx(1);0;0;]; %because xOy and xOZ sym.

% Generate convex hull for the flagella
[vtx_hull, elt_hull]= generate_flagella_convex_hull(5*L,L,rad,step_length,rad_section,ke_type, 20);
mmg_mesh = mmg(msh(vtx_hull, elt_hull));
mesh_hull1 = mmg_mesh.Mesh;

% =========================== Insertion point of flagella ===========================
if alpha ~= 0. | beta ~=0. 
  intersection_points = intersection_surf_line(x_S, mesh_head,alpha, beta); % intersection between head and centerline of flagella
    [maxi,argmax] = max(intersection_points(:,1)); % take the point with the largest x coordinate (the one at the back of the head) 
    xF1_ = intersection_points(argmax,:);   
else 
    xF1_ = zeros(1,3);
    xF1_(1) = max(mesh_head.vtx(:,1));  
end

% ============================= Computation of the cost function and the constraints ===========================
A=zeros(6+nb_flag,4); %store the results for the 4 rotation angles of the tail : 6 for the linear and angular velocities, nb_flag for the power dissipated by each flagella
sigma_n = zeros(3*nb_flag*size(mesh_tail.vtx,1) ,4);
viscosity=1;
joint_distance = 2*rad_section;

% Loop over 1 stroke decomposed into 4 rotation angles of the tail
for iii = 2:5

    % Rotate tail and translate it according to the insertion point and joint distance
    tail_rotation=(iii-1)/2; 

    R_rotated_tail=[1 0 0 ;
        0 cos(tail_rotation*pi) -sin(tail_rotation*pi);
        0 sin(tail_rotation*pi) cos(tail_rotation*pi)];

    mesh_tail1=mesh_tail;
    mesh_tail1.vtx = (R_rotated_tail*mesh_tail1.vtx')';

    R_gamma = [cos(gamma) 0 -sin(gamma);
    0 1 0;
    sin(gamma) 0 cos(gamma)];

    R_delta = [cos(delta) -sin(delta) 0;
    sin(delta) cos(delta) 0;
    0 0 1];
    
    mesh_tail1.vtx = ((R_delta*R_gamma*mesh_tail1.vtx') + xF1_'+ joint_distance*[cos(alpha)*cos(beta);sin(beta)*cos(alpha);sin(alpha)])';

    % Only one flagellum 
    if nb_flag == 1

        % Compute finite element spaces and domains for the head and the tail
        Gamma1 = dom(mesh_head,3); % head
        Gamma2 = dom(mesh_tail1,3); % tail
        phi1 = fem(mesh_head,'P1');
        phi2 = fem(mesh_tail1,'P1');

        % Build Green matrix
        tic;
        [G11] = stokesletEQ(Gamma1,phi1,viscosity,0,[]);
        disp('G11 computed');
        [G22] = stokesletEQ(Gamma2,phi2,viscosity,0,[]);
        disp('G22 computed');
        [G12] = crossstokesletEQ(Gamma1,Gamma2,phi1,phi2,viscosity,0,[]);
        disp('G12 computed');
        G21 = G12';
        %[G21] = crossstokesletEQ(Gamma2,Gamma1,phi2,phi1,viscosity,0,[]);
        disp('G21 computed');
      
        G=[G11 G12;
            G21 G22];
        toc
        disp('done G');

        % Build system matrix and right-hand side
        x_fk1 = xF1_'+ joint_distance*[cos(alpha)*cos(beta);sin(beta)*cos(alpha);sin(alpha)];
        omegae_fk1 = -2*pi*R_delta*R_gamma*[1;0;0];

        time=0;
        tic;
        lhs = build_lhs_head_tail_separated(G,Gamma1,Gamma2,phi1,phi2,x_S);
        rhs = build_rhs_head_tail_separated(time,Gamma1,Gamma2,phi1,phi2,omegae_fk1,x_fk1); % the tail's axis is aligned with x axis
        rhs = cell2mat(rhs);
        lhs=cell2mat(lhs);
        toc
        disp('done lhs-rhs');
    
        % Solve the system for surface stresses, linear and angular velocities
        solution = lhs\rhs;
    
        % Obtain the new position, orientation and surface stresses given the old ones
        U=solution(end-5:end-3);
        Omega=solution(end-2:end);
        sigma_n(:,iii-1) = solution(3*size(mesh_head.vtx,1)+1:end-6);

        % Compute the power dissipated by the flagella using the formula : P = -\omega \cdot (\int \sigma_n \wedge (X-XF))
        Ndof1 = size(mesh_tail1.vtx,1);

        sigma_n_1 = sigma_n(1:3*Ndof1, iii-1); %Flagella 1
        sigma_n_11 = sigma_n_1(1:Ndof1);
        sigma_n_21 = sigma_n_1(Ndof1+1:2*Ndof1);
        sigma_n_31 = sigma_n_1(2*Ndof1+1:3*Ndof1);

        r1_1 = @(X) X(:,1)-x_fk1(1); %Flagella 1
        r2_1 = @(X) X(:,2)-x_fk1(2);
        r3_1 = @(X) X(:,3)-x_fk1(3);
       
        % Compute (\int phi \wedge (X-XF)) for x, y and z (tail)
        Ix_1 = integral(Gamma2, phi2, r1_1); 
        Iy_1 = integral(Gamma2, phi2, r2_1);
        Iz_1 = integral(Gamma2, phi2, r3_1);

        % Compute (\int \sigma_n \wedge (X-XF)) for x, y and z (tail)
        I1y_1 = dot(sigma_n_11', Iy_1); 
        I1z_1 = dot(sigma_n_11', Iz_1);
        I2x_1 = dot(sigma_n_21', Ix_1);
        I2z_1 = dot(sigma_n_21', Iz_1);
        I3x_1 = dot(sigma_n_31', Ix_1);
        I3y_1 = dot(sigma_n_31', Iy_1);
        I_lhs_1 = [I2z_1-I3y_1; I3x_1-I1z_1; I1y_1-I2x_1]; 
        
        % Compute power dissipated for the flagella : -\omega \cdot (\int \sigma_n \wedge (X-XF))
        P1 = -dot(omegae_fk1, I_lhs_1); 
       
        % Store the results for the current rotation angle of the tail
        A(1:6,iii-1) = solution(end-5:end);
        A(7,iii-1) = P1;
    
    % Two flagella (symmetric configuration)
    elseif nb_flag == 2 
        
        mesh_tail2=mesh_tail1;

        R_secondtail=[1 0 0 ;
            0 cos(pi) -sin(pi);
            0 sin(pi) cos(pi)];

        mesh_tail2.vtx=(R_secondtail*mesh_tail2.vtx')'; % symmetry with respect to x axis
  
        % Compute finite element spaces and domains for the head and the flagella
        Gamma1 = dom(mesh_head,3); % head
        Gamma2=dom(mesh_tail1,3); % tail up
        Gamma3=dom(mesh_tail2,3); % tail down
        phi1 = fem(mesh_head,'P1');
        phi2 = fem(mesh_tail1,'P1');
        phi3 = fem(mesh_tail2,'P1');

        tic;
        time=0;

        % Build Green matrix
        [G11] = stokesletEQ(Gamma1,phi1,viscosity,0,[]);
        disp('done G11');
        [G22] = stokesletEQ(Gamma2,phi2,viscosity,0,[]);
        disp('done G22');
        [G33] = stokesletEQ(Gamma3,phi3,viscosity,0,[]);
        disp('done G33');
        [G12] = crossstokesletEQ(Gamma1,Gamma2,phi1,phi2,viscosity,0,[]);
        disp('done G12');
        [G13] = crossstokesletEQ(Gamma1,Gamma3,phi1,phi3,viscosity,0,[]);
        disp('done G13');
        [G23] = crossstokesletEQ(Gamma2,Gamma3,phi2,phi3,viscosity,0,[]);
        disp('done G23');
        G=[G11 G12 G13;
            G12' G22 G23;
            G13' G23' G33];
        toc
        disp('done G');

        % Build system matrix and right-hand side
        x_fk1 = xF1_'+ joint_distance*[cos(alpha)*cos(beta);sin(beta)*cos(alpha);sin(alpha)];
        x_fk2 = R_secondtail*x_fk1;

        omegae_fk1 = -2*pi*R_delta*R_gamma*[1;0;0];
        omegae_fk2 = -2*pi*R_secondtail*R_delta*R_gamma*[1;0;0];

        lhs = build_lhs_head_tail_separated_2tails(G,Gamma1,Gamma2,Gamma3,phi1,phi2,phi3,x_S);
        rhs = build_rhs_head_tail_separated_2tails(time,Gamma1,Gamma2,Gamma3,phi1,phi2,phi3,-2*pi*[1;0;0],...
            x_fk1,x_fk2,omegae_fk1,omegae_fk2);
        rhs = cell2mat(rhs);
        lhs=cell2mat(lhs);
     
        % Solve the system for surface stresses, linear and angular velocities
        solution = lhs\rhs;

        % Obtain the new position, orientation and surface stresses given the old ones
        U=solution(end-5:end-3);
        Omega=solution(end-2:end);
        sigma_n(:,iii-1) = solution(3*size(mesh_head.vtx,1)+1:end-6);

        % Compute the power dissipated by the flagellum 1 (because of symmetry, P1=P2) 
        Ndof1 = size(mesh_tail1.vtx,1);
       
        sigma_n_1 = sigma_n(1:3*Ndof1,iii-1); %Flagella 1
        sigma_n_11 = sigma_n_1(1:Ndof1);
        sigma_n_21 = sigma_n_1(Ndof1+1:2*Ndof1);
        sigma_n_31 = sigma_n_1(2*Ndof1+1:3*Ndof1); 

        r1_1 = @(X) X(:,1)-x_fk1(1); %Flagellum 1
        r2_1 = @(X) X(:,2)-x_fk1(2);
        r3_1 = @(X) X(:,3)-x_fk1(3);

        % Compute (\int phi \wedge (X-XF)) for x, y and z (flagellum 1)
        Ix_1 = integral(Gamma2, phi2, r1_1); %F1 : \int phi \wedge (X-XF)
        Iy_1 = integral(Gamma2, phi2, r2_1);
        Iz_1 = integral(Gamma2, phi2, r3_1);
        
        % Compute (\int \sigma_n \wedge (X-XF)) for x, y and z (flagellum 1)
        I1y_1 = dot(sigma_n_11', Iy_1); %F1 : coeff to compute the inf of vector
        I1z_1 = dot(sigma_n_11', Iz_1);
        I2x_1 = dot(sigma_n_21', Ix_1);
        I2z_1 = dot(sigma_n_21', Iz_1);
        I3x_1 = dot(sigma_n_31', Ix_1);
        I3y_1 = dot(sigma_n_31', Iy_1);
        I_lhs_1 = [I2z_1-I3y_1; I3x_1-I1z_1; I1y_1-I2x_1]; 

        % Compute power dissipated for the flagellum 1 : -\omega \cdot (\int \sigma_n \wedge (X-XF))
        P1 = -dot(omegae_fk1, I_lhs_1); 

        % Store the results for the current rotation angle of the flagella
        A(1:6,iii-1) = solution(end-5:end);
        A(7,iii-1) = P1;
        A(8,iii-1) = P1; %P2; because P2=P1
   
    end

    disp(U);
    disp(Omega);
    disp(A(7:end,iii-1));

end

% ============================= Average over the 4 rotation ===========================
mean_vel=mean(A,2);
U=mean_vel(1:3);
Omega=mean_vel(4:6);
P = mean_vel(7:7+nb_flag-1);

disp(U);
disp(Omega);
disp(sum(P));


% ============================= Volume of at least 2 intersection of convex hull =============================
mesh_hull1.vtx = ((R_delta*R_gamma*mesh_hull1.vtx') + xF1_'+ joint_distance*[cos(alpha)*cos(beta);sin(beta)*cos(alpha);sin(alpha)])';%(1+0.134)*[R1*cos(alpha);0;R3*sin(alpha)])' 
hulls = cell(1, nb_flag+1);
hulls{1} = struct('vertices', mesh_head.vtx, 'faces', mesh_head.elt);
hulls{2} = struct('vertices', mesh_hull1.vtx, 'faces', mesh_hull1.elt);

a1 = x_S(1); % inf bound x-component of square domain for the Monte-Carlo start from mass center (useless to take a bigger domain because intersection is in the back of the head)
a2 = min(mesh_head.vtx(:,2));
b2 = max(mesh_head.vtx(:,2));
a3 = min(mesh_head.vtx(:,3));
b3 = max(mesh_head.vtx(:,3));

if nb_flag == 1 
    
    b1 = max([max(mesh_head.vtx(:,1))+0.5, max(mesh_tail1.vtx(:,1))/2.]); % sup bound x-component of square domain (max between the back of the head and the half of the tail because intersection is not further than the half of the tail)

elseif nb_flag == 2

    % Compute the second flagellum's convex hull
    mesh_hull2=mesh_hull1;
    mesh_hull2.vtx=(R_secondtail*mesh_hull2.vtx')';
    hulls{3} = struct('vertices', mesh_hull2.vtx, 'faces', mesh_hull2.elt);
    b1 = max([max(mesh_head.vtx(:,1))+0.5,max(mesh_tail1.vtx(:,1))/2.,max(mesh_tail2.vtx(:,1))/2.]);

end

 
range = [a1, b1; a2, b2; a3, b3];
% Monte-Carlo method  to compute the volume of intersection of convex hulls  
[pointsInside, pointsOutside] = pointsInsideIntersection(hulls, numPointsMC, range);


% ============================= Return cost function and constraints =============================
if num_pb == 1
    objective = U(1);
elseif num_pb == 2
    objective = -(U(1)/sum(P))/(U0/P0); %Add - in order to minimise the objective 
end

constraints = [abs(U(2:3))-eps; size(pointsInside,1)];

