function publishit(scriptname)
    %PUBLISHIT
    publ = fileread(strcat(scriptname, '.m'));
    
    torun = sprintf('run %s', scriptname);
    eval(torun);
    
    publ = strrep(publ, 'PUBLISHABLE = true;', 'PUBLISHABLE = false;');
    

    while(1)
       ind1 = strfind(publ, '$(');
       ind2 = strfind(publ, ')$');
       len = ind2 - ind1;
       if ~isempty(len)
           try
            oldstr = publ(ind1(1): ind1(1) + len(1) +1);
            curstr =     ...
                    num2str(eval(publ(ind1(1) + 2:ind1(1) + len(1) - 1)));

            publ = strrep(publ, oldstr, curstr);
           catch
            publ = strrep(publ, oldstr, '');  
           end
       else
           break;
       end
    end

    %fId = fopen('testpublish.m', 'w+');
    write_to = strcat(scriptname, '_publish.m');
    fId = fopen(write_to, 'w+');
    fwrite(fId, publ);
    fclose(fId);
    
    clear options 
    options.format = 'latex';
    options.stylesheet = 'publish2latex.xsl';
    options.showCode = false;
    options.evalCode = true;
    options.imageFormat = 'png'; %TODO: eps umstellen
    
    publish(write_to, ...
        options ...
    )
    
    currentdir = cd('html');
    
    if exist('../protokoll.pdf', 'file')
        delete('../protokoll.pdf');
    end
    
    if system('/usr/texbin/pdflatex protokoll.tex') == 0
        copyfile('protokoll.pdf', '../protokoll.pdf')
    end
    
    cd(currentdir)
    
    open('protokoll.pdf')
end
%eval(publ(ind1 + 2:ind1 + len - 1))
