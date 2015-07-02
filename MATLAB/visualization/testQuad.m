load('Data.mat')

% n = size(v, 1) /17;
% for i= 1:n
%     tmp = v((i -1) * 17 +4:(i-1)*17 + 7);
%     tmp = tmp./norm(tmp);
%     v((i -1) * 17 +4:(i-1)*17 + 7) = tmp;
% end

t = linspace(0, 1, size(v, 1) /17);

visualizeQuad(t, v);