
function [tout,zout] = simu_Nlink3D(N,r,l,kLT,kLP,kHT,kHP,kHR,k_el,M,t0,tf,z0,timeinterp, ux, uy, uz, varargin)



B = @(t,z) [interp1(timeinterp, ux, t),interp1(timeinterp, uy, t), interp1(timeinterp, uz,t)]';



odefun = @(t,z) rhs(N,B,M,k_el,z,t); 
options = odeset('Mass',@(t,z) mass(N,kLT,kLP,kHT,kHP,kHR,r,l,z,t),'MStateDependence','Strong','MassSingular','no','RelTol',1e-4,'AbsTol',1e-6);

if isempty(varargin)
    [tout,zout]= ode15s(odefun,[t0,tf],z0,options);
else
     [tout,zout]= ode15s(odefun,linspace(0,tf,varargin{1}),z0,options);
end

end
