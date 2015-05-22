function [ A, B, Q, R ] = fill_Jk( input_args )
% Besetzt die einzelnen Matrizen der approximierten Hessematrix J_k
% Das Minimierungsproblem hat die Form sum_i (0.5*||f_i(s_i,q_i)||^2)

% A    enthält alle A_i = (grad_s_i h_i)^T
% B    enthält alle B_i = (grad_q_i h_i)^T
% Q    enthält alle Q_i = grad_s_i(f_i)*(grad_s_i(f_i))^T
% R    enthält alle R_i = grad_q_i(f_i)*(grad_q_i(f_i))^T

% die Größen der Variablen s und q
n_s      = 13;
n_q      = 4;

% die Größe der zu betrachtenden Zeitschritte
N        = ??;

% Erstelle Platzhalter in der passenden Größe
A = zeros((N-1)*n_s,n_s);
B = zeros((N-1)*n_q,n_q);
R = zeros((N-1)*n_q,n_q);
Q = zeros(N*n_s,n_s);

% Besetze die Matrizen
for i = 1:(N-1)
    A((i-1)*n_s+1:i*n_s;:) = % (grad_s_i h_i)^T
    B((i-1)*n_q+1:i*n_q;:) = % (grad_q_i h_i)^T
    R((i-1)*n_q+1:i*n_q;:) = % grad_q_i(f_i)*(grad_q_i(f_i))^T
    % sollte l_i(q_i) = q_i sein gilt
    % R = ones((N-1)*n_q,n_q))
    Q((i-1)*n_s+1:i*n_s;:) = % grad_s_i(f_i)*(grad_s_i(f_i))^T
end
Q((N-1)*n_s+1:N*n_s;:) = % grad_s_N(f_N)*(grad_s_N(f_N))^T

end

