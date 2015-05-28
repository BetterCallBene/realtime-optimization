currentpath = cd('..');
addpath(pwd);
cd(currentpath);
copyfile('../BasisConfiguration/script.m', 'script.m')
run generate_scripts

