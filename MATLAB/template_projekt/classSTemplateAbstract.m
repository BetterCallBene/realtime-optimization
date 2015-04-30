classdef(Abstract) classSTemplateAbstract < handle
    %classSTemplateAbstractSuper Diese Klasse soll ein Hilfsmittel zum schnellen Erstellen einer abstrakten Klasse sein.
    %   Von Abstrakten Klassen k�nnen keine Instanzen erstellt werden. Dies
    %   Klasse schreibt deren Unterklassen vor, welchen properties und
    %   Methoden enthalten sein m�ssen.
    %   Diese Klasse ist von der Klasse handle abgeleitet und wird somit zum "referenzierbaren" Objekt. 
    %   N�here Informationen gibt es hier: http://de.mathworks.com/help/matlab/ref/handle-class.html?searchHighlight=handle
    
    
    properties(Abstract, Constant) % Folgende properties m�ssen in den Unterklassen �berschrieben werden
        %Gravitation
        g;
    end
    
    methods(Abstract) % Folgende methods m�ssen in den Unterklassen �berschrieben werden
        %Gibt die Gravitation zur�ck
        GetGravitation(obj);
    end
    
end

