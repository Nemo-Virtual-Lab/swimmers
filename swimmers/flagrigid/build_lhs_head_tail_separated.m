function lhs = build_lhs_head_tail_separated(G,Gamma1,Gamma2,phi1,phi2,x_S);
% Function building the left-hand side of the global system
% Dimension (3*Ndof+6) x (3*Ndof+6)
% |   G   0.5*lhs21'     lhs13  |   |    0     |
% | 0.5*lhs21       0          0    | = |    0.5*0     |
% | 0.5*lhs31       0          0    |   |-0.5*tau_magn |

Ndof1 = size(phi1.unk,1); %wall
Ndof2 = size(phi2.unk,1); %swimmer

lhs{1,1} = G;

% Dimension 3 x 3*Ndof
%         | 1 1 ... 1 0 0 ... 0 0 .... 0 |
% lhs21 = | 0 0 ... 0 1 1 ... 1 0 .... 0 |
%         | 0 0 ... 0 0 0 ... 0 1 .... 1 |
g=@(X) ones(size(X,1),1);
lhs21 = sparse(3,3*(Ndof1+Ndof2));
lhs21(1,1:Ndof1) = integral(Gamma1,phi1,g);
lhs21(2,Ndof1+1:2*Ndof1) = integral(Gamma1,phi1,g);
lhs21(3,2*Ndof1+1:3*Ndof1) = integral(Gamma1,phi1,g);

lhs21(1,3*Ndof1+1:3*Ndof1+Ndof2) = integral(Gamma2,phi2,g);
lhs21(2,3*Ndof1+Ndof2+1:3*Ndof1+2*Ndof2) = integral(Gamma2,phi2,g);
lhs21(3,3*Ndof1+2*Ndof2+1:end) = integral(Gamma2,phi2,g);

lhs{2,1} = lhs21;
lhs{1,2} = 2*0.5*lhs21';  %  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IL 2 L'HO AGGIUNTO 9/12


% Dimension 3 x 3*Ndof
%                          | (r3-r2)*phi    0         0       |
% lhs31 = int_{\partial S} |      0    (r1-r3)*phi    0       |
%                          |      0         0     (r2-r1)*phi |

r1 = @(X) X(:,1)-x_S(1);
r2 = @(X) X(:,2)-x_S(2);
r3 = @(X) X(:,3)-x_S(3);

% 0.5*Int_{\partial S} \Omega \wedge (x-x_S) \phi(x) \dS(x)
lhs13 =  sparse(3*(Ndof1+Ndof2),3);
lhs13(1:Ndof1,2)        =  0.5*integral(Gamma1,phi1,r3);
lhs13(1:Ndof1,3)        = -0.5*integral(Gamma1,phi1,r2);
lhs13(Ndof1+1:2*Ndof1,1) = -0.5*integral(Gamma1,phi1,r3);
lhs13(Ndof1+1:2*Ndof1,3) =  0.5*integral(Gamma1,phi1,r1);
lhs13(2*Ndof1+1:3*Ndof1,1)  =  0.5*integral(Gamma1,phi1,r2);
lhs13(2*Ndof1+1:3*Ndof1,2)  = -0.5*integral(Gamma1,phi1,r1);

lhs13(3*Ndof1+1:3*Ndof1+Ndof2,2)        =  0.5*integral(Gamma2,phi2,r3);
lhs13(3*Ndof1+1:3*Ndof1+Ndof2,3)        = -0.5*integral(Gamma2,phi2,r2);
lhs13(3*Ndof1+Ndof2+1:3*Ndof1+2*Ndof2,1) = -0.5*integral(Gamma2,phi2,r3);
lhs13(3*Ndof1+Ndof2+1:3*Ndof1+2*Ndof2,3) =  0.5*integral(Gamma2,phi2,r1);
lhs13(3*Ndof1+2*Ndof2+1:end,1)  =  0.5*integral(Gamma2,phi2,r2);
lhs13(3*Ndof1+2*Ndof2+1:end,2)  = -0.5*integral(Gamma2,phi2,r1);

lhs{1,3} = 2*lhs13; %  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IL 2 L'HO AGGIUNTO 9/12
lhs{3,1} = 2*lhs13';

lhs{2,2} = sparse(3,3);
lhs{2,3} = sparse(3,3);
lhs{3,2} = sparse(3,3);
lhs{3,3} = sparse(3,3);
end