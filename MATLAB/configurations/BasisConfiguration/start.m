
currentpath = pwd;

if exist('MapleConfDesktop', 'dir') == 0
    mkdir('MapleConfDesktop')
end

if exist('MapleConfDesktop/GenerateConstraints.mw', 'file') == 0
    copyfile('../BasisConfiguration/MapleConfDesktop/GenerateConstraints.mw', 'MapleConfDesktop/GenerateConstraints.mw')
end

if exist('MapleConfDesktop/GenerateJacobiHesse.mw', 'file') == 0
    copyfile('../BasisConfiguration/MapleConfDesktop/GenerateJacobiHesse.mw', 'MapleConfDesktop/GenerateJacobiHesse.mw')
end

if exist('html', 'dir') == 0
    mkdir('html')
end
path_rtoptmain = strcat(currentpath, '/rtoptmain.m');
if exist(path_rtoptmain, 'file') == 0
    copyfile('../BasisConfiguration/rtoptmain.m', 'rtoptmain.m')
end
path_init = strcat(currentpath, '/init.m');
if exist(path_init, 'file') == 0
    copyfile('../BasisConfiguration/init.m', 'init.m')
end

path_protokoll = strcat(currentpath, 'html/protokoll.tex');
if exist(path_protokoll, 'file') == 0
    copyfile('../BasisConfiguration/Latex/protokoll.tex', 'html/protokoll.tex')
end



