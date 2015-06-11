
clear
currentpath = cd('..');
cd('..');
rtopt_path = strcat(pwd, '/rtopt/');
cd('..')
python_path = strcat(pwd, '/PYTHON/');
addpath(rtopt_path, '/Library/Frameworks/Maple.framework/Versions/Current/bin');
cd(currentpath)

template_path1 = strcat(rtopt_path, 'templates/template_BasisGenQDyn.m');
script_path1 = strcat(python_path, 'GenerateBasisGenQDyn.py');

template_path2 = strcat(rtopt_path, 'templates/template_GenConstraints.m');
script_path2 = strcat(python_path, 'GenerateGenConstraints.py');

% %prg = strcat('python ', script_path,  ' "', template_path, '"');
prg1 = strcat('python "', script_path1, '"');
prg1 = strcat(prg1, ' "', template_path1, '"');

result1 =  system(prg1);

prg2 = strcat('python "', script_path2, '"');
prg2 = strcat(prg2, ' "', template_path2, '"');
% ToDo: Inklude
result2 =  system(prg2);

if result1 == 0  && result2 == 0
    run runtests
else
    error('Pythonskripte waren nicht erfolgreich')
end