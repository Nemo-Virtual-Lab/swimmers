function rhs = build_rhs_head_tail_separated_2tails(time,Gamma1,Gamma2,Gamma3,phi1,phi2,phi3,omega,x_fk1,x_fk2,omegae_fk1,omegae_fk2)
% BUILDS THE RIGHT-HAND SIDE OF THE SYSTEM WHICH SOLVES FOR SURFACE 
% STRESSES, LINEAR AND ANGULAR VELOCITIES
rhs = cell(3,1);
Ndof1 = size(phi1.unk,1);
Ndof2 = size(phi2.unk,1);
Ndof3 = size(phi3.unk,1);


rhs_aux = zeros((Ndof1+Ndof2+Ndof3)*3,3);


r1_1 = @(X) X(:,1)-x_fk1(1);
r2_1 = @(X) X(:,2)-x_fk1(2);
r3_1 = @(X) X(:,3)-x_fk1(3);
r1_2 = @(X) X(:,1)-x_fk2(1);
r2_2 = @(X) X(:,2)-x_fk2(2);
r3_2 = @(X) X(:,3)-x_fk2(3);


rhs_aux(3*Ndof1+1:3*Ndof1+Ndof2,2)        =  0.5*integral(Gamma2,phi2,r3_1);
rhs_aux(3*Ndof1+1:3*Ndof1+Ndof2,3)        = -0.5*integral(Gamma2,phi2,r2_1);
rhs_aux(3*Ndof1+Ndof2+1:3*Ndof1+2*Ndof2,1) = -0.5*integral(Gamma2,phi2,r3_1);
rhs_aux(3*Ndof1+Ndof2+1:3*Ndof1+2*Ndof2,3) =  0.5*integral(Gamma2,phi2,r1_1);
rhs_aux(3*Ndof1+2*Ndof2+1:3*(Ndof1+Ndof2),1)  =  0.5*integral(Gamma2,phi2,r2_1);
rhs_aux(3*Ndof1+2*Ndof2+1:3*(Ndof1+Ndof2),2)  = -0.5*integral(Gamma2,phi2,r1_1);

rhs_aux(3*(Ndof1+Ndof2)+1:3*(Ndof1+Ndof2)+Ndof3,2)        =  0.5*integral(Gamma3,phi3,r3_2);
rhs_aux(3*(Ndof1+Ndof2)+1:3*(Ndof1+Ndof2)+Ndof3,3)        = -0.5*integral(Gamma3,phi3,r2_2);
rhs_aux(3*(Ndof1+Ndof2)+Ndof3+1:3*(Ndof1+Ndof2)+2*Ndof3,1) = -0.5*integral(Gamma3,phi3,r3_2);
rhs_aux(3*(Ndof1+Ndof2)+Ndof3+1:3*(Ndof1+Ndof2)+2*Ndof3,3) =  0.5*integral(Gamma3,phi3,r1_2);
rhs_aux(3*(Ndof1+Ndof2)+2*Ndof3+1:end,1)  =  0.5*integral(Gamma3,phi3,r2_2);
rhs_aux(3*(Ndof1+Ndof2)+2*Ndof3+1:end,2)  = -0.5*integral(Gamma3,phi3,r1_2);


rhs_1=-[rhs_aux(1:3*Ndof1,:)*omega; % it's made of zeros only
        rhs_aux(3*Ndof1+1:3*(Ndof1+Ndof2),:)*omegae_fk1;
        rhs_aux(3*(Ndof1+Ndof2)+1:end,:)*omegae_fk2]*2; % dec 2019

rhs{1}=rhs_1;


rhs{2} = zeros(3,1);%/(mu*L*U); % 0.5 is due to lhs{2,1} which 
% is multiplied by 0.5 to make it equal to lhs{1,2}; the minus sign is due
% to the fact that the SUM of FORCES = 0 = F_{ext} + \int \sigma n = 0
rhs{3} = zeros(3,1);%/(mu*L^2*U); % net torque over the swimmer is the external one,
% with a minus sign because lhs{3,1} = -lhs{1,3} (cfr. how the lhs is
% constructed); the sum of all torques is zero
end