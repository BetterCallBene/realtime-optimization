function  visualizeQuad( t,v)
%VISUALIZEQUAD Summary of this function goes here
%   Detailed explanation goes here

global visualization

Q = zeros(length(t), 12);
array = reshape(v, [17, length(t)])'; 

Q(:, 10:12) = array(:, 1:3);
[Q(:, 4), Q(:, 5), Q(:, 6)] = quat2angle((array(:, 4:7)));
%[Q(:, 7:9)] = array(:, 8:10);
%[Q(:, 10:12)] = array(:, 11:13);

visualization.yout = Q;
visualization.tout = t;

QuadAnim4;

end

