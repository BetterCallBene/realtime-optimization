restoredefaultpath
addpath('helper');
currentpath= cd('..');
addpath(pwd);
addpath(strcat(pwd, '/BasisConfiguration/'));
cd('..');
addpath(strcat(pwd, '/rtopt/'));
cd(currentpath);
addpath('/usr/texbin')
addpath('/usr/local/bin')