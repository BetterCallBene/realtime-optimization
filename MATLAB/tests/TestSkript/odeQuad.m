function odeQuad()
global u;

tspan = [0.025 0.05];

val = rand(13, 1);
u = rand(4, 1);

opts = odeset('RelTol',1e-4,'AbsTol',1e-5, 'Jacobian',@Jac);

val(4:7) = val(4:7)./norm(val(4:7));

y0 = val;
yp0 = f([], y0);

% Compute consistent initial conditions.
tic
[y01,yp01] = decic(@Fimpl,tspan(1),y0,[],yp0,[],opts);


sol = ode15i(@Fimpl,tspan,y01,yp01,opts);
toc

yend = sol.y(:, 11);
% Norm der Quaternionen
Q = norm(yend(4:7))
% Differenz zur 1
if abs(Q - 1) < 1e-5 % in der Toleranz von 1e-5
    disp('Yes');
end
end

function res = Fimpl(t,y,yp)
      
      res = mass(t,y)*yp - f(t,y);
end 

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

function res= f(t, val)

robot = Robot();
[r, q,v,omega,u_,Iges,IM,m,kT,kQ,d,g] = getParams(val, robot);


t117 = q(4);
t144 = -2 * t117 ^ 2;
t118 = q(3);
t143 = -2 * t118 ^ 2;
t142 = (d * kT);
t141 = -t117 / 0.2e1;
t140 = -t118 / 0.2e1;
t119 = q(2);
t139 = -t119 / 0.2e1;
t120 = q(1);
t138 = t120 / 0.2e1;
t137 = 0.1e1 - 0.2e1 * t119 ^ 2;
t113 = u_(4);
t97 = (t113 ^ 2);
t115 = u_(2);
t99 = (t115 ^ 2);
t136 = (-t97 - t99);
t116 = u_(1);
t100 = (t116 ^ 2);
t114 = u_(3);
t98 = (t114 ^ 2);
t135 = (t100 + t98);
t134 = t118 * t117;
t133 = t119 * t117;
t132 = t119 * t118;
t131 = t120 * t117;
t130 = t120 * t118;
t129 = t120 * t119;
t128 = 2 * g;
t110 = v(3);
t127 = 2 * t110;
t111 = v(2);
t126 = 2 * t111;
t112 = v(1);
t125 = 2 * t112;
t124 = t143 + t137;
t123 = t130 - t133;
t122 = -t129 - t134;
t121 = ((t114 - t113 + t116 - t115) * IM);
t109 = Iges(1);
t108 = Iges(2);
t107 = Iges(3);
t106 = omega(1);
t105 = omega(2);
t104 = omega(3);
res = [((1 + t143 + t144) * t112) + (-t131 + t132) * t126 + (t130 + t133) * t127; (t131 + t132) * t125 + (t144 + t137) * t111 + (-t129 + t134) * t127; -t123 * t125 - t122 * t126 + t124 * t110; (t106 * t139 + t105 * t140 + t104 * t141) * t120; (t106 * t138 + t105 * t141 + t118 * t104 / 0.2e1) * t119; (t117 * t106 / 0.2e1 + t105 * t138 + t104 * t139) * t118; (t104 * t138 + t119 * t105 / 0.2e1 + t106 * t140) * t117; -t105 * t110 + t104 * t111 + t123 * t128; -t104 * t112 + t106 * t110 + t122 * t128; ((t135 - t136) * kT + (-t106 * t111 + t105 * t112 - t124 * g) * m) / m; -((-t99 + t97) * t142 + ((t107 - t108) * t104 + t121) * t105) / t109; ((-t100 + t98) * t142 + ((-t109 + t107) * t104 + t121) * t106) / t108; -((t108 - t109) * t106 * t105 + (t135 + t136) * kQ) / t107;];


% res =  [v(1) - 2 * v(1) * q(3) ^ 2 - 2 * v(1) * q(4) ^ 2 - 2 * v(2) * q(1) * q(4) + 2 * v(2) * q(2) * q(3) + 2 * v(3) * q(1) * q(3) + 2 * v(3) * q(2) * q(4); 
%         2 * v(1) * q(1) * q(4) + 2 * v(1) * q(2) * q(3) + v(2) - 2 * v(2) * q(2) ^ 2 - 2 * v(2) * q(4) ^ 2 - 2 * v(3) * q(1) * q(2) + 2 * v(3) * q(3) * q(4); 
%         -2 * v(1) * q(1) * q(3) + 2 * v(1) * q(2) * q(4) + 2 * v(2) * q(1) * q(2) + 2 * v(2) * q(3) * q(4) + v(3) - 2 * v(3) * q(2) ^ 2 - 2 * v(3) * q(3) ^ 2; 
%         -(q(1) * q(2) * omega(1)) / 0.2e1 - (q(1) * q(3) * omega(2)) / 0.2e1 - (q(1) * q(4) * omega(3)) / 0.2e1; 
%         (q(1) * q(2) * omega(1)) / 0.2e1 - (q(2) * q(4) * omega(2)) / 0.2e1 + (q(2) * q(3) * omega(3)) / 0.2e1; 
%         (q(3) * q(4) * omega(1)) / 0.2e1 + (q(1) * q(3) * omega(2)) / 0.2e1 - (q(2) * q(3) * omega(3)) / 0.2e1; 
%         (q(1) * q(4) * omega(3)) / 0.2e1 + (q(2) * q(4) * omega(2)) / 0.2e1 - (q(3) * q(4) * omega(1)) / 0.2e1; 
%         -omega(2) * v(3) + omega(3) * v(2) + 2 * g * q(1) * q(3) - 2 * g * q(2) * q(4); 
%         -omega(3) * v(1) + omega(1) * v(3) - 2 * g * q(1) * q(2) - 2 * g * q(3) * q(4); 
%         (kT * u_(1) ^ 2 + kT * u_(2) ^ 2 + kT * u_(3) ^ 2 + kT * u_(4) ^ 2 - m * omega(1) * v(2) + m * omega(2) * v(1) - m * g + 2 * m * g * q(2) ^ 2 + 2 * m * g * q(3) ^ 2) / m; 
%         -(-d * kT * u_(2) ^ 2 + d * kT * u_(4) ^ 2 + IM * omega(2) * u_(1) - IM * omega(2) * u_(2) + IM * omega(2) * u_(3) - IM * omega(2) * u_(4) + omega(2) * Iges(3) * omega(3) - omega(3) * Iges(2) * omega(2)) / Iges(1); 
%         (-d * kT * u_(1) ^ 2 + d * kT * u_(3) ^ 2 + IM * omega(1) * u_(1) - IM * omega(1) * u_(2) + IM * omega(1) * u_(3) - IM * omega(1) * u_(4) - omega(3) * Iges(1) * omega(1) + omega(1) * Iges(3) * omega(3)) / Iges(2); 
%         -(kQ * u_(1) ^ 2 - kQ * u_(2) ^ 2 + kQ * u_(3) ^ 2 - kQ * u_(4) ^ 2 + omega(1) * Iges(2) * omega(2) - omega(2) * Iges(1) * omega(1)) / Iges(3);];

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