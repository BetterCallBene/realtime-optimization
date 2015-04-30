% template_unit_test Dies ist ein kleines Beispiel für einen unit test in
% Klassenform
%   Bei diesem Beispiel schlägt der Test fehl, dabei gibt es
%   unterschiedlichen unter Stufen des fehlschlagens, Verify, Assume,
%   Fatal und Assert
%   Mehr Informationen dazu unter: http://de.mathworks.com/help/matlab/matlab_prog/author-class-based-unit-tests-in-matlab.html
classdef classTemplateUnitTests < matlab.unittest.TestCase

    methods(Test)
        function testclassTemplate(testCase)
            templateClass = classTemplate();
            [g1, g2] = templateClass.GetGravitation();
            testCase.assumeEqual(g1, 5);
        end
        
        function testclassTemplateAbstractAssert(testCase)
            templateClass = classTemplateAbstract();
            g = templateClass.GetGravitation();
            testCase.assertEqual(g, 5);
        end
        
         
         
         
         function testclassTemplateAbstractVerify(testCase)
             template_Abstract_Class = classTemplateAbstract();
             g = template_Abstract_Class.GetGravitation();
             testCase.verfiyEqual(g, 5);
         end
        
       
        function testclassTemplateAbstractFatal(testCase)
           template_Abstract_Class = classTemplateAbstract();
           g = template_Abstract_Class.GetGravitation();
           testCase.fatalEqual(g, 5);
        end
    end
%abstract_superclass = template_abstract_class();
%abstract_superclass.GetGravitation()


%mc = ?template_abstract_superclass; %Greift auf die Metadaten von template_abstract_superclass zu
%if ~mc.Abstract %Test, ob Klasse template_abstract_superclass abstract ist
   % not an abstract class
   
%end

end