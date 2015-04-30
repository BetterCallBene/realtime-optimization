classdef classTemplate < classSTemplate
    %classTemplate Unterklasse von classSTemplate
    %   Beispiel einer abgeleiten Klasse.
    
    
    properties(Constant, Access=private)
        %Mondkonstante
        gMond = 1.622; 
    end
    
    
    methods(Access=public)
        %Constructor der classTemplate
        function obj = classTemplate()
            obj@classSTemplate('classTemplate'); %Aufruf des Superclass Konstruktor  
        end
        %Gibt die Gravitationskonstante des Mondes zurück
        function [neug, altg] = GetGravitation(obj) %Diese Funktion überschreibt die Funktion von classSTemplate
            altg = GetGravitation@classSTemplate(obj); % Ruft die Funktion Gravitation in der Superclass auf.
            neug = obj.gMond;
        end
    end
    
end

