function odeQuad()
global u;

tspan = [0 1];

val = rand(13, 1);
u = rand(4, 1);

opts = odeset('RelTol',1e-5,'AbsTol',1e-4);

val(4:7) = val(4:7)./norm(val(4:7));

y0 = val;
yp0 = f([], y0);

% Compute consistent initial conditions.
[y01,yp01] = decic(@Fimpl,tspan(1),y0,[],yp0,[],opts);

tic
sol = ode15i(@Fimpl,tspan,y01,yp01,opts);
toc

end

function res = Fimpl(t,y,yp)
      % Burgers' equation -- fully implicit formulation
      res = mass(t,y)*yp - f(t,y);
end 

function res = mass(t, y)
     res = diag([ones(3, 1); y(4:7); ones(6, 1)]);
end


function res= f(t, val)

robot = Robot();
[r, q,v,omega,u_,Iges,IM,m,kT,kQ,d,g] = getParams(val, robot);

res =  [v(1) - 2 * v(1) * q(3) ^ 2 - 2 * v(1) * q(4) ^ 2 - 2 * v(2) * q(1) * q(4) + 2 * v(2) * q(2) * q(3) + 2 * v(3) * q(1) * q(3) + 2 * v(3) * q(2) * q(4); 
        2 * v(1) * q(1) * q(4) + 2 * v(1) * q(2) * q(3) + v(2) - 2 * v(2) * q(2) ^ 2 - 2 * v(2) * q(4) ^ 2 - 2 * v(3) * q(1) * q(2) + 2 * v(3) * q(3) * q(4); 
        -2 * v(1) * q(1) * q(3) + 2 * v(1) * q(2) * q(4) + 2 * v(2) * q(1) * q(2) + 2 * v(2) * q(3) * q(4) + v(3) - 2 * v(3) * q(2) ^ 2 - 2 * v(3) * q(3) ^ 2; 
        -(q(1) * q(2) * omega(1)) / 0.2e1 - (q(1) * q(3) * omega(2)) / 0.2e1 - (q(1) * q(4) * omega(3)) / 0.2e1; 
        (q(1) * q(2) * omega(1)) / 0.2e1 - (q(2) * q(4) * omega(2)) / 0.2e1 + (q(2) * q(3) * omega(3)) / 0.2e1; 
        (q(3) * q(4) * omega(1)) / 0.2e1 + (q(1) * q(3) * omega(2)) / 0.2e1 - (q(2) * q(3) * omega(3)) / 0.2e1; 
        (q(1) * q(4) * omega(3)) / 0.2e1 + (q(2) * q(4) * omega(2)) / 0.2e1 - (q(3) * q(4) * omega(1)) / 0.2e1; 
        -omega(2) * v(3) + omega(3) * v(2) + 2 * g * q(1) * q(3) - 2 * g * q(2) * q(4); 
        -omega(3) * v(1) + omega(1) * v(3) - 2 * g * q(1) * q(2) - 2 * g * q(3) * q(4); 
        (kT * u_(1) ^ 2 + kT * u_(2) ^ 2 + kT * u_(3) ^ 2 + kT * u_(4) ^ 2 - m * omega(1) * v(2) + m * omega(2) * v(1) - m * g + 2 * m * g * q(2) ^ 2 + 2 * m * g * q(3) ^ 2) / m; 
        -(-d * kT * u_(2) ^ 2 + d * kT * u_(4) ^ 2 + IM * omega(2) * u_(1) - IM * omega(2) * u_(2) + IM * omega(2) * u_(3) - IM * omega(2) * u_(4) + omega(2) * Iges(3) * omega(3) - omega(3) * Iges(2) * omega(2)) / Iges(1); 
        (-d * kT * u_(1) ^ 2 + d * kT * u_(3) ^ 2 + IM * omega(1) * u_(1) - IM * omega(1) * u_(2) + IM * omega(1) * u_(3) - IM * omega(1) * u_(4) - omega(3) * Iges(1) * omega(1) + omega(1) * Iges(3) * omega(3)) / Iges(2); 
        -(kQ * u_(1) ^ 2 - kQ * u_(2) ^ 2 + kQ * u_(3) ^ 2 - kQ * u_(4) ^ 2 + omega(1) * Iges(2) * omega(2) - omega(2) * Iges(1) * omega(1)) / Iges(3);];

end


function [r, q,v,omega,u_,Iges,IM,m,kT,kQ,d,g] = getParams(val, robot)
    global u;    
    r   = val(1:3);
    q   = val(4:7);
    v   = val(8:10);
    omega   = val(11:13);
    
    u_ = u;
    
    Iges = robot.I;
    IM = robot.I_M;
    m = robot.m;

    kT  = robot.kT;
    kQ  = robot.kQ;
    d   = robot.d;
    g   = 9.81;
end