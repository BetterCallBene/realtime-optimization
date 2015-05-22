function [ DLk ] = fill_DLk( input_args )
% Berechne den Gradienten DLk der Lagrangefunktion L_k

% L_k(y_k) = sum(F_i(s_i,q_i))+F_N(s_N,q_N) + lambda_k^T(x_k-s_k) + sum(
% lamda_k+1^T(h_i(s_i,q_i)-s_i+1)

% die Größen der Variablen lamda, s und q
n_s      = 13;
n_q      = 4;
%n_lambda = n_s;

% die Größe der zu betrachtenden Zeitschritte
N        = ??;

% Erstelle Platzhalter 
DLk = zeros((N-1)*(2*n_s+n_q)+(2*n_s),1);

% Berechne die Werte für DLk
% DLk(1:n_s) bleibt 0, da dieser Wert mit x_k berechnet wird und zurzeit
% noch nicht bekannt ist
DLk(n_s+1:n_s+n_s) = %grad_s_1(F_i(s_1,q_1))-lambda_1+grad_s_1(h_1)*lambda_2
DLk(n_s+n_s+1:n_s+n_s+n_q) = %grad_q_1(F_1(s_1,q_1))+grad_q_1(h_1)*lambda_2
for i = 2:(N-1)
    DLk((i-1)*(2*n_s+n_q)+1:(i-1)*(2*n_s+n_q)+n_s) = h_i-1(s_i-1,q_i-1)-s_i
    DLk((i-1)*(2*n_s+n_q)+n_s+1:(i-1)*(2*n_s+n_q)+n_s+n_s) = %grad_s_i(F_i(s_i,q_i))-lambda_i+grad_s_i(h_i)*lambda_i+1
    DLk((i-1)*(2*n_s+n_q)+n_s+n_s+1:i*(2*n_s+n_q)) = %grad_q_i(F_i(s_i,q_i))+grad_q_i(h_i)*lambda_i+1
end
DLk((N-1)*(2*n_s+n_q)+1:(N-1)*(2*n_s+n_q)+n_s) = %h_N-1(s_N-1,q_N-1)-s_N
DLk((N-1)*(2*n_s+n_q)+n_s+1:(N-1)*(2*n_s+n_q)+n_s+n_s) = %grad_s_N(F_N(s_N))-lambda_N



end

