function odeQuad()
global u;

tspan = [0 0.5];

val = [zeros(3, 1); 1; 0 ;0 ;0; zeros(6, 1)];
%val = rand(13, 1);
u = [0, 10000, 10000+20, 10000];



val(4:7) = val(4:7)./norm(val(4:7));

y0 = val;
yp0 = f([], y0);

%opts = odeset('RelTol',1e-5,'AbsTol',1e-8, 'Mass', @mass, 'MassSingular', 'yes', 'InitialSlope', yp0);
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
% Compute consistent initial conditions.
tic
%[y01,yp01] = decic(@Fimpl,tspan(1),y0,[],yp0,[],opts);

%sol = ode15i(@Fimpl,tspan,y0,yp0, opts);

%tic
sol = ode15i(@Fimpl,tspan,y0,yp0, opts);
sol1 = ode45(@f1, tspan, y0, opts);
%toc

[yaw, pitch, roll] = quat2angle((sol.y(4:7, :))');

% Norm der Quaternionen
 for i = 1:size(sol.y, 1)
      disp('ode15')
      norm(sol.y( 4:7,i)) -1
      disp('ode45')
      norm(sol1.y( 4:7,i)) -1
      disp('Integrationsunterschied')
      norm(abs(sol1.y(1:3) - sol.y(1:3)))
 end
end

% 
function res = Fimpl(t,x,xp)
    
    Fun = f(t, x);
    
    res = zeros(13, 1);
    
    if abs(x(4)) < 1e-7
        x(4) = sign(x(4))* 1e-7;
    end
    
    res(1:3, 1) = xp(1:3, 1) - Fun(1:3, 1);
    res(4:6, 1) = xp(5:7, 1) - Fun(5:7, 1);
    res(7:end-1, 1) = xp(8:end) - Fun(8:end);
    res(end, 1) = x(4:7)'*xp(4:7);
    %res(end, 1) = x(4)^2 + x(5)^2 + x(6)^2 + x(7)^2 - 1;    
end 
% 
% function res = Fimpl2(t,x,xp)
%     
%     y = [x(1:6); 0; x(7:end)];
%     
%     Fun = f(t, y);
%     res = zeros(12, 1);
%     
%     res(1:3, 1) = xp(1:3, 1) - Fun(1:3, 1);
%     res(4:5, 1) = xp(4:5, 1) - Fun(4:5, 1);
%     res(6:end-1, 1) = xp(7:end) - Fun(8:end);
%     
%     res(end, 1) = x(4)^2 + x(5)^2 + x(6)^2 - 1;
%     
% end 

% function res = Fimpl(t, x)
%     res = f(t,x);
% end

function res = mass(t, y)
     res = diag([ones(3, 1); y(4:7); ones(6, 1)]);
end

function [dfdy,dfdyp] = Jac(t,y,yp)
      J = JacY(y);
      Jx = J(1:13, 1:13);
      dfdy = Jx; 
      dfdyp = mass(t,y);
end


function [dfdy] = JacY(y)
robot = Robot();
[r, q,v,omega,u_,Iges,IM,m,kT,kQ,d,g] = getParams(y, robot);
dfdy = [0 0 0 -2 * v(2) * q(4) + 2 * v(3) * q(3) 2 * v(2) * q(3) + 2 * v(3) * q(4) -4 * v(1) * q(3) + 2 * v(2) * q(2) + 2 * v(3) * q(1) -4 * v(1) * q(4) - 2 * v(2) * q(1) + 2 * v(3) * q(2) 1 - 2 * q(3) ^ 2 - 2 * q(4) ^ 2 -2 * q(1) * q(4) + 2 * q(2) * q(3) 2 * q(1) * q(3) + 2 * q(2) * q(4) 0 0 0 0 0 0 0; 0 0 0 2 * v(1) * q(4) - 2 * v(3) * q(2) 2 * v(1) * q(3) - 4 * v(2) * q(2) - 2 * v(3) * q(1) 2 * v(1) * q(2) + 2 * v(3) * q(4) 2 * v(1) * q(1) - 4 * v(2) * q(4) + 2 * v(3) * q(3) 2 * q(1) * q(4) + 2 * q(2) * q(3) 1 - 2 * q(2) ^ 2 - 2 * q(4) ^ 2 -2 * q(1) * q(2) + 2 * q(3) * q(4) 0 0 0 0 0 0 0; 0 0 0 -2 * v(1) * q(3) + 2 * v(2) * q(2) 2 * v(1) * q(4) + 2 * v(2) * q(1) - 4 * v(3) * q(2) -2 * v(1) * q(1) + 2 * v(2) * q(4) - 4 * v(3) * q(3) 2 * v(1) * q(2) + 2 * v(2) * q(3) -2 * q(1) * q(3) + 2 * q(2) * q(4) 2 * q(1) * q(2) + 2 * q(3) * q(4) 1 - 2 * q(2) ^ 2 - 2 * q(3) ^ 2 0 0 0 0 0 0 0; 0 0 0 -(q(2) * omega(1)) / 0.2e1 - (q(3) * omega(2)) / 0.2e1 - (q(4) * omega(3)) / 0.2e1 -(q(1) * omega(1)) / 0.2e1 -(q(1) * omega(2)) / 0.2e1 -(q(1) * omega(3)) / 0.2e1 0 0 0 -(q(1) * q(2)) / 0.2e1 -(q(1) * q(3)) / 0.2e1 -(q(1) * q(4)) / 0.2e1 0 0 0 0; 0 0 0 (q(2) * omega(1)) / 0.2e1 (q(1) * omega(1)) / 0.2e1 - (q(4) * omega(2)) / 0.2e1 + (q(3) * omega(3)) / 0.2e1 (q(2) * omega(3)) / 0.2e1 -(q(2) * omega(2)) / 0.2e1 0 0 0 (q(1) * q(2)) / 0.2e1 -(q(2) * q(4)) / 0.2e1 (q(2) * q(3)) / 0.2e1 0 0 0 0; 0 0 0 (q(3) * omega(2)) / 0.2e1 -(q(3) * omega(3)) / 0.2e1 (q(4) * omega(1)) / 0.2e1 + (q(1) * omega(2)) / 0.2e1 - (q(2) * omega(3)) / 0.2e1 (q(3) * omega(1)) / 0.2e1 0 0 0 (q(3) * q(4)) / 0.2e1 (q(1) * q(3)) / 0.2e1 -(q(2) * q(3)) / 0.2e1 0 0 0 0; 0 0 0 (q(4) * omega(3)) / 0.2e1 (q(4) * omega(2)) / 0.2e1 -(q(4) * omega(1)) / 0.2e1 -(q(3) * omega(1)) / 0.2e1 + (q(2) * omega(2)) / 0.2e1 + (q(1) * omega(3)) / 0.2e1 0 0 0 -(q(3) * q(4)) / 0.2e1 (q(2) * q(4)) / 0.2e1 (q(1) * q(4)) / 0.2e1 0 0 0 0; 0 0 0 2 * g * q(3) -2 * g * q(4) 2 * g * q(1) -2 * g * q(2) 0 omega(3) -omega(2) 0 -v(3) v(2) 0 0 0 0; 0 0 0 -2 * g * q(2) -2 * g * q(1) -2 * g * q(4) -2 * g * q(3) -omega(3) 0 omega(1) v(3) 0 -v(1) 0 0 0 0; 0 0 0 0 4 * g * q(2) 4 * g * q(3) 0 omega(2) -omega(1) 0 -v(2) v(1) 0 2 * kT * u_(1) / m 2 * kT * u_(2) / m 2 * kT * u_(3) / m 2 * kT * u_(4) / m; 0 0 0 0 0 0 0 0 0 0 0 -(IM * u_(1) - IM * u_(2) + IM * u_(3) - IM * u_(4) + Iges(3) * omega(3) - omega(3) * Iges(2)) / Iges(1) omega(2) * (-Iges(3) + Iges(2)) / Iges(1) -IM * omega(2) / Iges(1) (2 * d * kT * u_(2) + IM * omega(2)) / Iges(1) -IM * omega(2) / Iges(1) (-2 * d * kT * u_(4) + IM * omega(2)) / Iges(1); 0 0 0 0 0 0 0 0 0 0 (IM * u_(1) - IM * u_(2) + IM * u_(3) - IM * u_(4) - omega(3) * Iges(1) + Iges(3) * omega(3)) / Iges(2) 0 -omega(1) * (Iges(1) - Iges(3)) / Iges(2) (-2 * d * kT * u_(1) + IM * omega(1)) / Iges(2) -IM * omega(1) / Iges(2) (2 * d * kT * u_(3) + IM * omega(1)) / Iges(2) -IM * omega(1) / Iges(2); 0 0 0 0 0 0 0 0 0 0 omega(2) * (-Iges(2) + Iges(1)) / Iges(3) omega(1) * (-Iges(2) + Iges(1)) / Iges(3) 0 -2 * kQ * u_(1) / Iges(3) 2 * kQ * u_(2) / Iges(3) -2 * kQ * u_(3) / Iges(3) 2 * kQ * u_(4) / Iges(3);];

end


function res= f1(t, val)

robot = Robot();
[r, q,v,omega,u,Iges,IM,m,kT,kQ,d,g] = getParams(val, robot);

res = [(1 - 2 * q(3) ^ 2 - 2 * q(4) ^ 2) * v(1) + (-2 * q(1) * q(4) + 2 * q(2) * q(3)) * v(2) + (2 * q(1) * q(3) + 2 * q(2) * q(4)) * v(3); (2 * q(1) * q(4) + 2 * q(2) * q(3)) * v(1) + (1 - 2 * q(2) ^ 2 - 2 * q(4) ^ 2) * v(2) + (-2 * q(1) * q(2) + 2 * q(3) * q(4)) * v(3); (-2 * q(1) * q(3) + 2 * q(2) * q(4)) * v(1) + (2 * q(1) * q(2) + 2 * q(3) * q(4)) * v(2) + (1 - 2 * q(2) ^ 2 - 2 * q(3) ^ 2) * v(3); -(q(2) * omega(1)) / 0.2e1 - (q(3) * omega(2)) / 0.2e1 - (q(4) * omega(3)) / 0.2e1 + ((1 - q(1) ^ 2 - q(2) ^ 2 - q(3) ^ 2 - q(4) ^ 2) * q(1)); (q(1) * omega(1)) / 0.2e1 - (q(4) * omega(2)) / 0.2e1 + (q(3) * omega(3)) / 0.2e1 + ((1 - q(1) ^ 2 - q(2) ^ 2 - q(3) ^ 2 - q(4) ^ 2) * q(2)); (q(4) * omega(1)) / 0.2e1 + (q(1) * omega(2)) / 0.2e1 - (q(2) * omega(3)) / 0.2e1 + ((1 - q(1) ^ 2 - q(2) ^ 2 - q(3) ^ 2 - q(4) ^ 2) * q(3)); -(q(3) * omega(1)) / 0.2e1 + (q(2) * omega(2)) / 0.2e1 + (q(1) * omega(3)) / 0.2e1 + ((1 - q(1) ^ 2 - q(2) ^ 2 - q(3) ^ 2 - q(4) ^ 2) * q(4)); 1 / m * (-m * (omega(2) * v(3) - omega(3) * v(2)) - m * (-2 * q(1) * q(3) + 2 * q(2) * q(4)) * g); 1 / m * (-m * (omega(3) * v(1) - omega(1) * v(3)) - m * (2 * q(1) * q(2) + 2 * q(3) * q(4)) * g); 1 / m * (kT * (u(1) ^ 2 + u(2) ^ 2 + u(3) ^ 2 + u(4) ^ 2) - m * (omega(1) * v(2) - omega(2) * v(1)) - m * (1 - 2 * q(2) ^ 2 - 2 * q(3) ^ 2) * g); 1 / Iges(1) * (d * kT * u(2) ^ 2 - d * kT * u(4) ^ 2 - (u(1) - u(2) + u(3) - u(4)) * IM * omega(2) - omega(2) * Iges(3) * omega(3) + omega(3) * Iges(2) * omega(2)); 1 / Iges(2) * (-d * kT * u(1) ^ 2 + d * kT * u(3) ^ 2 + (u(1) - u(2) + u(3) - u(4)) * IM * omega(1) - omega(3) * Iges(1) * omega(1) + omega(1) * Iges(3) * omega(3)); 1 / Iges(3) * (-kQ * u(1) ^ 2 + kQ * u(2) ^ 2 - kQ * u(3) ^ 2 + kQ * u(4) ^ 2 - omega(1) * Iges(2) * omega(2) + omega(2) * Iges(1) * omega(1));];
end

function res= f(t, val)

robot = Robot();
[r, q,v,omega,u,Iges,IM,m,kT,kQ,d,g] = getParams(val, robot);

%res = [v(1) - 2 * v(1) * q(3) ^ 2 - 2 * v(1) * q(4) ^ 2 - 2 * v(2) * q(1) * q(4) + 2 * v(2) * q(2) * q(3) + 2 * v(3) * q(1) * q(3) + 2 * v(3) * q(2) * q(4); 2 * v(1) * q(1) * q(4) + 2 * v(1) * q(2) * q(3) + v(2) - 2 * v(2) * q(2) ^ 2 - 2 * v(2) * q(4) ^ 2 - 2 * v(3) * q(1) * q(2) + 2 * v(3) * q(3) * q(4); -2 * v(1) * q(1) * q(3) + 2 * v(1) * q(2) * q(4) + 2 * v(2) * q(1) * q(2) + 2 * v(2) * q(3) * q(4) + v(3) - 2 * v(3) * q(2) ^ 2 - 2 * v(3) * q(3) ^ 2; -(q(2) * omega(1)) / 0.2e1 - (q(3) * omega(2)) / 0.2e1 - (q(4) * omega(3)) / 0.2e1 + q(1) - (q(1) ^ 3) - (q(1) * q(2) ^ 2) - (q(1) * q(3) ^ 2) - (q(1) * q(4) ^ 2); (q(1) * omega(1)) / 0.2e1 - (q(4) * omega(2)) / 0.2e1 + (q(3) * omega(3)) / 0.2e1 + q(2) - (q(2) * q(1) ^ 2) - (q(2) ^ 3) - (q(2) * q(3) ^ 2) - (q(2) * q(4) ^ 2); (q(4) * omega(1)) / 0.2e1 + (q(1) * omega(2)) / 0.2e1 - (q(2) * omega(3)) / 0.2e1 + q(3) - (q(3) * q(1) ^ 2) - (q(3) * q(2) ^ 2) - (q(3) ^ 3) - (q(3) * q(4) ^ 2); -(q(3) * omega(1)) / 0.2e1 + (q(2) * omega(2)) / 0.2e1 + (q(1) * omega(3)) / 0.2e1 + q(4) - (q(4) * q(1) ^ 2) - (q(4) * q(2) ^ 2) - (q(4) * q(3) ^ 2) - (q(4) ^ 3); -omega(2) * v(3) + omega(3) * v(2) + 2 * g * q(1) * q(3) - 2 * g * q(2) * q(4); -omega(3) * v(1) + omega(1) * v(3) - 2 * g * q(1) * q(2) - 2 * g * q(3) * q(4); (kT * u(1) ^ 2 + kT * u(2) ^ 2 + kT * u(3) ^ 2 + kT * u(4) ^ 2 - m * omega(1) * v(2) + m * omega(2) * v(1) - m * g + 2 * m * g * q(2) ^ 2 + 2 * m * g * q(3) ^ 2) / m; -(-d * kT * u(2) ^ 2 + d * kT * u(4) ^ 2 + IM * omega(2) * u(1) - IM * omega(2) * u(2) + IM * omega(2) * u(3) - IM * omega(2) * u(4) + omega(2) * Iges(3) * omega(3) - omega(3) * Iges(2) * omega(2)) / Iges(1); (-d * kT * u(1) ^ 2 + d * kT * u(3) ^ 2 + IM * omega(1) * u(1) - IM * omega(1) * u(2) + IM * omega(1) * u(3) - IM * omega(1) * u(4) - omega(3) * Iges(1) * omega(1) + omega(1) * Iges(3) * omega(3)) / Iges(2); -(kQ * u(1) ^ 2 - kQ * u(2) ^ 2 + kQ * u(3) ^ 2 - kQ * u(4) ^ 2 + omega(1) * Iges(2) * omega(2) - omega(2) * Iges(1) * omega(1)) / Iges(3);];


res  = [(1 - 2 * q(3) ^ 2 - 2 * q(4) ^ 2) * v(1) + (-2 * q(1) * q(4) + 2 * q(2) * q(3)) * v(2) + (2 * q(1) * q(3) + 2 * q(2) * q(4)) * v(3); (2 * q(1) * q(4) + 2 * q(2) * q(3)) * v(1) + (1 - 2 * q(2) ^ 2 - 2 * q(4) ^ 2) * v(2) + (-2 * q(1) * q(2) + 2 * q(3) * q(4)) * v(3); (-2 * q(1) * q(3) + 2 * q(2) * q(4)) * v(1) + (2 * q(1) * q(2) + 2 * q(3) * q(4)) * v(2) + (1 - 2 * q(2) ^ 2 - 2 * q(3) ^ 2) * v(3); -(q(2) * omega(1)) / 0.2e1 - (q(3) * omega(2)) / 0.2e1 - (q(4) * omega(3)) / 0.2e1; (q(1) * omega(1)) / 0.2e1 - (q(4) * omega(2)) / 0.2e1 + (q(3) * omega(3)) / 0.2e1; (q(4) * omega(1)) / 0.2e1 + (q(1) * omega(2)) / 0.2e1 - (q(2) * omega(3)) / 0.2e1; -(q(3) * omega(1)) / 0.2e1 + (q(2) * omega(2)) / 0.2e1 + (q(1) * omega(3)) / 0.2e1; 1 / m * (-m * (omega(2) * v(3) - omega(3) * v(2)) - m * (-2 * q(1) * q(3) + 2 * q(2) * q(4)) * g); 1 / m * (-m * (omega(3) * v(1) - omega(1) * v(3)) - m * (2 * q(1) * q(2) + 2 * q(3) * q(4)) * g); 1 / m * (kT * (u(1) ^ 2 + u(2) ^ 2 + u(3) ^ 2 + u(4) ^ 2) - m * (omega(1) * v(2) - omega(2) * v(1)) - m * (1 - 2 * q(2) ^ 2 - 2 * q(3) ^ 2) * g); 1 / Iges(1) * (d * kT * u(2) ^ 2 - d * kT * u(4) ^ 2 - (u(1) - u(2) + u(3) - u(4)) * IM * omega(2) - omega(2) * Iges(3) * omega(3) + omega(3) * Iges(2) * omega(2)); 1 / Iges(2) * (-d * kT * u(1) ^ 2 + d * kT * u(3) ^ 2 + (u(1) - u(2) + u(3) - u(4)) * IM * omega(1) - omega(3) * Iges(1) * omega(1) + omega(1) * Iges(3) * omega(3)); 1 / Iges(3) * (-kQ * u(1) ^ 2 + kQ * u(2) ^ 2 - kQ * u(3) ^ 2 + kQ * u(4) ^ 2 - omega(1) * Iges(2) * omega(2) + omega(2) * Iges(1) * omega(1));];



% res =  [v(1) - 2 * v(1) * q(3) ^ 2 - 2 * v(1) * q(4) ^ 2 - 2 * v(2) * q(1) * q(4) + 2 * v(2) * q(2) * q(3) + 2 * v(3) * q(1) * q(3) + 2 * v(3) * q(2) * q(4); 
%          2 * v(1) * q(1) * q(4) + 2 * v(1) * q(2) * q(3) + v(2) - 2 * v(2) * q(2) ^ 2 - 2 * v(2) * q(4) ^ 2 - 2 * v(3) * q(1) * q(2) + 2 * v(3) * q(3) * q(4); 
%          -2 * v(1) * q(1) * q(3) + 2 * v(1) * q(2) * q(4) + 2 * v(2) * q(1) * q(2) + 2 * v(2) * q(3) * q(4) + v(3) - 2 * v(3) * q(2) ^ 2 - 2 * v(3) * q(3) ^ 2; 
%          -(q(1) * q(2) * omega(1)) / 0.2e1 - (q(1) * q(3) * omega(2)) / 0.2e1 - (q(1) * q(4) * omega(3)) / 0.2e1; 
%          (q(1) * q(2) * omega(1)) / 0.2e1 - (q(2) * q(4) * omega(2)) / 0.2e1 + (q(2) * q(3) * omega(3)) / 0.2e1; 
%          (q(3) * q(4) * omega(1)) / 0.2e1 + (q(1) * q(3) * omega(2)) / 0.2e1 - (q(2) * q(3) * omega(3)) / 0.2e1; 
%          (q(1) * q(4) * omega(3)) / 0.2e1 + (q(2) * q(4) * omega(2)) / 0.2e1 - (q(3) * q(4) * omega(1)) / 0.2e1; 
%          -omega(2) * v(3) + omega(3) * v(2) + 2 * g * q(1) * q(3) - 2 * g * q(2) * q(4); 
%          -omega(3) * v(1) + omega(1) * v(3) - 2 * g * q(1) * q(2) - 2 * g * q(3) * q(4); 
%          (kT * u_(1) ^ 2 + kT * u_(2) ^ 2 + kT * u_(3) ^ 2 + kT * u_(4) ^ 2 - m * omega(1) * v(2) + m * omega(2) * v(1) - m * g + 2 * m * g * q(2) ^ 2 + 2 * m * g * q(3) ^ 2) / m; 
%          -(-d * kT * u_(2) ^ 2 + d * kT * u_(4) ^ 2 + IM * omega(2) * u_(1) - IM * omega(2) * u_(2) + IM * omega(2) * u_(3) - IM * omega(2) * u_(4) + omega(2) * Iges(3) * omega(3) - omega(3) * Iges(2) * omega(2)) / Iges(1); 
%          (-d * kT * u_(1) ^ 2 + d * kT * u_(3) ^ 2 + IM * omega(1) * u_(1) - IM * omega(1) * u_(2) + IM * omega(1) * u_(3) - IM * omega(1) * u_(4) - omega(3) * Iges(1) * omega(1) + omega(1) * Iges(3) * omega(3)) / Iges(2); 
%          -(kQ * u_(1) ^ 2 - kQ * u_(2) ^ 2 + kQ * u_(3) ^ 2 - kQ * u_(4) ^ 2 + omega(1) * Iges(2) * omega(2) - omega(2) * Iges(1) * omega(1)) / Iges(3);];

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