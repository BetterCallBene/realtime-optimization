function [ P, new_gsi ] = riccati( DLk, A, B, Q, R )
% Berechnet P_i und new_grad_s_i mittels Riccati-Recursion

% die Größen der Variablen lamda, s und q
n_s      = 13;
n_q      = 4;
%n_lambda = n_s;

% die Größe der zu betrachtenden Zeitschritte
N        = ??;

P = zeros(N*n_s,n_s);
new_gsi = zeros(N*n_s,1);

P((N-1)*n_s+1:N*n_s,:) = Q((N-1)*n_s+1:N*n_s;:);
new_gsi((N-1)*n_s+1:N*n_s) = DLk((N-1)*(2*n_s+n_q)+n_s+1:(N-1)*(2*n_s+n_q)+n_s+n_s);

for i = N:2
    %setze zunächst die Matrizen, die in der Formel verwendet werden
    P_i = P((i-1)*n_s+1:i*n_s,:);
    n_gsi = new_gsi((i-1)*n_s+1:i*n_s);
    o_gsim = DLk((i-2)*(2*n_s+n_q)+n_s+1:(i-1)*(2*n_s+n_q)+n_s+n_s);
    gli = DLk((i-1)*(2*n_s+n_q)+1:i*(2*n_s+n_q)+n_s);
    gqim = DLk((i-2)*(2*n_s+n_q)+n_s+n_s+1:(i-1)*(2*n_s+n_q)+n_s+n_s+n_q);
    
    A_im = A((i-2)*n_s+1:(i-1)*n_s,:);
    B_im = B((i-2)*n_q+1:(i-1)*n_q,:);
    Q_im = Q((i-2)*n_s+1:(i-1)*n_s,:);
    R_im = R((i-2)*n_q+1:(i-1)*n_q,:);
    
    % Berechne P_{i-1}=P_im
    P_im = Q_im + A_im'*P_i*A_im - A_im'*P_i*B_im *(R_im + B_im'*P_i*B_im)\(B_im'*P_i*A_im);
    n_gsim = o_gsim + A_im'*P_i*gli + A_im'*n_gsi - (R_im + B_im'*P_i*B_im)\(gqim + B_im'*P_i*gli + B_im'*n_gsi);
    
    % setze das Ergebnis in die Speichermatrix ein
    P((i-2)*n_s+1:(i-1)*n_s,:) = P_im;
    new_gsi((i-2)*n_s+1:(i-1)*n_s) = n_gsim;
end


end

