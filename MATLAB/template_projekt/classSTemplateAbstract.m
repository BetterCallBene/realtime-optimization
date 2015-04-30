classdef(Abstract) classSTemplateAbstract < handle
    %classSTemplateAbstractSuper Diese Klasse soll ein Hilfsmittel zum schnellen Erstellen einer abstrakten Klasse sein.
    %   Von Abstrakten Klassen können keine Instanzen erstellt werden. Dies
    %   Klasse schreibt deren Unterklassen vor, welchen properties und
    %   Methoden enthalten sein müssen.
    %   Diese Klasse ist von der Klasse handle abgeleitet und wird somit zum "referenzierbaren" Objekt. 
    %   Nähere Informationen gibt es hier: http://de.mathworks.com/help/matlab/ref/handle-class.html?searchHighlight=handle
    
    
    properties(Abstract, Constant) % Folgende properties müssen in den Unterklassen überschrieben werden
        %Gravitation
        g;
    end
    
    methods(Abstract) % Folgende methods müssen in den Unterklassen überschrieben werden
        %Gibt die Gravitation zurück
        GetGravitation(obj);
    end
    
end

