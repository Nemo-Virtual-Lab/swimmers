function Mass = Nlink_3D_Mass(N,kLT,kLP,kHT,kHP,kHR,r,l,z,t)
cross_m = @(A) [ 0 -A(3) A(2)
                 A(3) 0  -A(1)
                -A(2) A(1) 0];
Rot_x = @(x) [ 1              , 0 ,          0
               0              ,cos(x), -sin(x)
               0              ,sin(x), cos(x)] ;

Rot_y = @(x) [ cos(x) ,       0,     sin(x)
                0,            1,       0
              -sin(x),        0  , cos(x) ];
          
Rot_z = @(x) [ cos(x) , -sin(x), 0
               sin(x), cos(x),  0
               0  , 0      , 1 ];

theta0 = z(4:6);


R_head =  feval(Rot_x,theta0(1))*feval(Rot_y,theta0(2))*feval(Rot_z,theta0(3));
R_tail  = zeros(3,3,N);
L_tail  = zeros(3,2,N);

    for i = 1:N 
        alpha = z(7 + 2*(i-1) : 6 + 2*i);
        R_tail(:,:,i) = feval(Rot_y,alpha(1))*feval(Rot_z,alpha(2));
        L_tail(:,:,i) = [ 0    sin(alpha(1))
                         1     0  
                            0              cos(alpha(1))]; 
    end
          
 
       

   L_Omega0 = [           1                             0                         sin(theta0(2))  
                         0                         cos(theta0(1))       -cos(theta0(2))*sin(theta0(1))   
                         0                        sin(theta0(1))          cos(theta0(1))*cos(theta0(2))];
              
                     
    F_head = -[R_head*(N*l)*diag([kHT,kHP,kHP])*R_head',zeros(3),zeros(3,2*N)];
    T_head = -[zeros(3),kHR*(N*l)^3*L_Omega0,zeros(3,2*N)];


S_i = zeros(3,3,N);
D = -diag([kLT,kLP,kLP]);
for i = 1:N 
   S_i(:,:,i) = R_head*R_tail(:,:,i)*D*R_tail(:,:,i)'*R_head'; 
end

A_i = zeros(3,2*N+6,N);
ex_cross = feval(cross_m,[1,0,0]');  
for i = 1:N 
      A_i(:,1:3,i) = S_i(:,:,i)*l;
      A_i(:,4:6,i) = R_head*R_tail(:,:,i)*D*R_tail(:,:,i)'*(l*r*ex_cross + l^2/2*feval(cross_m,R_tail(:,:,i)*[1,0,0]')...,
                     +l^2*feval(cross_m,sum(R_tail(:,:,1:i-1),3)*[1,0,0]'))*R_head'*L_Omega0;
      for j = 1:i
          if j<i 
            A_i(:,2*j+5:2*j+6,i) = R_head*R_tail(:,:,i)*D*R_tail(:,:,i)'*R_tail(:,:,j)*l^2*feval(cross_m,[1,0,0]')*L_tail(:,:,j);
          else
            A_i(:,2*j+5:2*j+6,i) = R_head*R_tail(:,:,i)*D*l^2/2*ex_cross*L_tail(:,:,i);  
          end
      end
end

B_i = zeros(3,2*N+6,N);   
       for i = 1:N 
          B_i(:,1:3,i) = S_i(:,:,i)*l^2/2;
          B_i(:,4:6,i) = R_head*R_tail(:,:,i)*D*R_tail(:,:,i)'*(r*l^2/2*ex_cross + l^3/3*feval(cross_m,R_tail(:,:,i)*[1,0,0]')...,
                         +l^3/2*feval(cross_m,sum(R_tail(:,:,1:i-1),3)*[1,0,0]'))*R_head'*L_Omega0;
          for j = 1:i
              if j<i 
                B_i(:,2*j+5:2*j+6,i) = R_head*R_tail(:,:,i)*D*R_tail(:,:,i)'*R_tail(:,:,j)*l^3/2*feval(cross_m,[1,0,0]')*L_tail(:,:,j);
              else
                B_i(:,2*j+5:2*j+6,i) = R_head*R_tail(:,:,i)*D*l^3/3*ex_cross*L_tail(:,:,i);  
              end
          end
       end

M_ij = zeros(3,2*N+6,N,N);
        for i = 1:N
            for j = 1:i
                sumRhRk = zeros(3,1);
                RhRi = R_head*R_tail(:,:,i);
                for k = j:i-1
                    RhRk = R_head*R_tail(:,:,k);
                    sumRhRk = RhRk(:,1) + sumRhRk;
                end
                M_ij(:,:,i,j) = -l*feval(cross_m,sumRhRk)*A_i(:,:,i) -feval(cross_m,RhRi(:,1))*B_i(:,:,i);  
            end 
        end
Mass = zeros(2*N + 6, 2*N + 6); 
Mass(1:3,:) = sum(A_i,3) + F_head;
Mass(4:6,:) = sum(M_ij(:,:,1:N,1),3) -feval(cross_m,r*R_head*[1,0,0]')*sum(A_i,3) +T_head;     
       for i = 1:N
              temp = R_tail(:,:,i)'*R_head'*sum(M_ij(:,:,i:N,i),3);
              Mass(2*i+5:2*i+6,:) = temp(2:3,:);
       end
 

end



