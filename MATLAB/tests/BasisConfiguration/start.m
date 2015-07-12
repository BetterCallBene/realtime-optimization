
currentpath = pwd;
path_maple = strcat(currentpath, '/MapleConfDesktop');
if exist(path_maple, 'dir') == 0
    mkdir('MapleConfDesktop')
end
path_maple_constraints = strcat(currentpath, '/MapleConfDesktop/GenerateConstraints.mw');
if exist(path_maple_constraints, 'file') == 0
    copyfile('../BasisConfiguration/MapleConfDesktop/GenerateConstraints.mw', 'MapleConfDesktop/GenerateConstraints.mw')
end

path_maple_hessejacobi = strcat(currentpath, '/MapleConfDesktop/GenerateJacobiHesse.mw');
if exist(path_maple_hessejacobi, 'file') == 0
    copyfile('../BasisConfiguration/MapleConfDesktop/GenerateJacobiHesse.mw', 'MapleConfDesktop/GenerateJacobiHesse.mw')
end
path_html = strcat(currentpath, '/html');
if exist('html', 'dir') == 0
    mkdir('html')
end
path_rtoptmain = strcat(currentpath, '/runtests.m');
if exist(path_rtoptmain, 'file') == 0
    copyfile('../BasisConfiguration/runtests.m', 'runtests.m')
end
path_init = strcat(currentpath, '/init.m');
if exist(path_init, 'file') == 0
    copyfile('../BasisConfiguration/init.m', 'init.m')
end

path_protokoll = strcat(currentpath, 'html/protokoll.tex');
if exist(path_protokoll, 'file') == 0
    copyfile('../BasisConfiguration/Latex/protokoll.tex', 'html/protokoll.tex')
end

helper_dir = strcat(currentpath, 'helper/');
if exist(helper_dir, 'dir') == 0
    copyfile('../BasisConfiguration/helper/*', 'helper')
end

run init



