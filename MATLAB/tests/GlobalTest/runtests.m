%RUN_TEST Diese Datei starten die Tests f�r template.
%   Matlab besitzt sogenanntes unit test framework. Dies ist ein kleines
%   Beispiel dazu
%   Mehr Informationen gibt es unter: http://de.mathworks.com/help/matlab/matlab_prog/create-simple-test-suites.html
import matlab.unittest.TestSuite

global TEST;
TEST = true;


%Create Suite from SolverTest Class
%The fromClass method creates a suite from all Test methods in the SolverTest class.

%suiteBasisQDyn = BasisQDyn();  Check: BB: 05.07.2015
%resultDyn = run(suiteBasisQDyn);


%suiteode15sM = ode15sM();  Check: BB: 05.07.2015
%resultsuiteode15sM = run(suiteode15sM);



%suiteInt = TestInt;  Check: BB: 05.07.2015
%resultInt = run(suiteInt);

%suiteQuadrocopterExt = QuadrocopterExt; %Check: BB: 05.07.2015
%suiteQuadrocopterExt.steadyPoint = []; % SteadyPoint initialisieren.
%resultQuadrocopterExt = run(suiteQuadrocopterExt, 'testWind'); % Nur Windfunktion testen

%suiteMultiShooting = TestSuite.fromClass(?MultiShooting); 
%resultMultiShooting = run(suiteMultiShooting);


suiteConstraints = TestSuite.fromClass(?Constraints);
resultConstraints = run(suiteConstraints);






%suiteCosts = TestSuite.fromClass(?Costs);
%resultCosts = run(suiteCosts);

%Create Suite from SolverTest Class Definition File
%The fromFile method creates a suite using the name of the file to identify the class.

%suiteFile = TestSuite.fromFile('SolverTest.m');
%result = run(suiteFile);