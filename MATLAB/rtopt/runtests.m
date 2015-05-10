%RUN_TEST Diese Datei starten die Tests für template.
%   Matlab besitzt sogenanntes unit test framework. Dies ist ein kleines
%   Beispiel dazu
%   Mehr Informationen gibt es unter: http://de.mathworks.com/help/matlab/matlab_prog/create-simple-test-suites.html
import matlab.unittest.TestSuite


%Create Suite from SolverTest Class
%The fromClass method creates a suite from all Test methods in the SolverTest class.

suiteClassQuadrocopter = TestSuite.fromClass(?classQuadrocopter);
result = run(suiteClassQuadrocopter);

suiteClassQuadrocopterDyn = TestSuite.fromClass(?classQuadrocopterDyn);
resultDyn = run(suiteClassQuadrocopterDyn);

%Create Suite from SolverTest Class Definition File
%The fromFile method creates a suite using the name of the file to identify the class.

%suiteFile = TestSuite.fromFile('SolverTest.m');
%result = run(suiteFile);