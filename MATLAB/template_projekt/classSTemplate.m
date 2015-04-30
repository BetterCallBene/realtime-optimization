classdef classSTemplate < handle 
    %classSTemplate Diese Klasse soll ein Hilfsmittel zum schnellen Erstellen einer Superclass sein. 
    %   Diese Klasse ist von der Klasse handle abgeleitet und wird somit zum "referenzierbaren" Objekt 
    %   Nähere Informationen gibt es hier: http://de.mathworks.com/help/matlab/ref/handle-class.html?searchHighlight=handle
    
    %Beispiel von Konstanten Attributen
    properties(Constant, Access=private)%In diesem Bereich werden Konstanten definiert
        %Gravitationskonstante
        gErde = 9.81;
    end
    
    properties(GetAccess=private, SetAccess=private) %Auf diese Variablen kann man von außen und allen Unterklassen nicht zugreifen.
        
    end
    
    properties(GetAccess=protected, SetAccess=protected) %Auf diese Variablen kann man von außen nicht zugreifen, aber in
        %allen Unterklassen
        
        %Interner Klassenname
        intName;
    end
    
    properties(GetAccess=public, SetAccess=public)
        %Auf diese Variablen kann man von außen zugreifen
    end
    
    properties(Dependent, GetAccess=public, SetAccess=public) % Die Variablen in diesen Bereich speichern keine Werte und sind deshalb für Get/Set Methoden vorgesehen
        %Öffentlicher Name
        Name;
    end
        
    methods(Access=public) % Der Zugriff in dem Bereich ist public
        function obj = classSTemplate(varargin)% Variable Anzahl von Parameter
            %Dies ist der Konstruktor der Class. Diese Funktion wird immer
            %bei ersten der Klasse aufgerufen.
            
            if (nargin == 1) %
                if isa(varargin{1}, 'char')
                    obj.intName = varargin{1}; % Auf ersten Parameterzugreifen    
                else 
                    error('Name should be a text');
                end
            elseif (nargin == 0)
                obj.intName = 'classSTemplate';
            else
                error('wrong number of inputs'); 
            end
        end
    
        %Gibt die Gravitationskonstante der Erde zurück
        function [altg, neug] = GetGravitation(obj)
            altg = obj.gErde;
            neug = 0;
        end
    end
    
    methods % get/set Methoden erlauben keine Attribute, wie (Access = public)
        %Gibt Name der Klasse zurück
        function value = get.Name(obj)
            value =  obj.intName;
        end
        %Bestimmt den Namen der Klasse
        function set.Name(obj, value)
            obj.intName = value;
        end
    end 
end