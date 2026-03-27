function RHS = RHS(N,B,M,k_el,z,t)

    cross_m = @(A) [ 0 -A(3) A(2)
                     A(3) 0  -A(1)
                    -A(2) A(1) 0];
    Rot_x = @(x) [ 1              , 0 ,          0
                   0              ,cos(x), -sin(x)
                   0              ,sin(x), cos(x)] ;

    Rot_y = @(x) [ cos(x) ,       0,     sin(x)
                    0,            1,       0
                  - sin(x),        0  , cos(x) ];

    Rot_z = @(x) [ cos(x) , -sin(x), 0
                   sin(x), cos(x),  0
                    0  , 0      , 1 ];
    theta0 = z(4:6);

    R_head =  feval(Rot_x,theta0(1))*feval(Rot_y,theta0(2))*feval(Rot_z,theta0(3));
    R_tail  = zeros(3,3,N);
        for i = 1:N 
            alpha = z(7 + 2*(i-1) : 6 + 2*i);
            R_tail(:,:,i) = feval(Rot_y,alpha(1))*feval(Rot_z,alpha(2));
        end


    Mag = R_head*[M,0,0]';

    B_t = feval(B,t);       

    M_mag = feval(cross_m,Mag)*B_t;

    RHS = zeros(2*N+6,1); 
    RHS(4:6) = -M_mag;
    ex_h = R_head(:,1);
    ex_1 = R_head*R_tail(:,:,1)*[1,0,0]';
    temp = -k_el*(feval(cross_m,ex_1)*ex_h); 
    RHS(7:8) = temp(2:3,:);
    for i = 2:N
        ex_i = R_head*R_tail(:,:,i)*[1,0,0]';
        ex_imoins1 = R_head*R_tail(:,:,i-1)*[1,0,0]';
        temp = -k_el*(feval(cross_m,ex_i)*ex_imoins1); 
        RHS(2*i+5:2*i+6) = temp(2:3);
    end
    
    
    end

