%RUN_TEST Diese Datei starten die Tests f�r template.
%   Matlab besitzt sogenanntes unit test framework. Dies ist ein kleines
%   Beispiel dazu
%   Mehr Informationen gibt es unter: http://de.mathworks.com/help/matlab/matlab_prog/create-simple-test-suites.html
import matlab.unittest.TestSuite

global TEST;
TEST = true;


%Create Suite from SolverTest Class
%The fromClass method creates a suite from all Test methods in the SolverTest class.

%suiteBasisQDyn = BasisQDyn(); %Check
%resultDyn = run(suiteBasisQDyn);

suiteForwEuler = ForwEuler(); %Check
resultForwEuler = run(suiteForwEuler);



%suiteode15sM = ode15sM(); %Check
%resultsuiteode15sM = run(suiteode15sM);

%suiteode45M = ode45M(); %Check
%resultsuiteode45M = run(suiteode45M);


suiteInt = TestInt; %Check
resultInt = run(suiteInt);


%suiteConstraints = TestSuite.fromClass(?Constraints);
%resultConstraints = run(suiteConstraints);

%suiteMultiShooting = TestSuite.fromClass(?MultiShooting); %Check
%resultMultiShooting = run(suiteMultiShooting);




%suiteCosts = TestSuite.fromClass(?Costs);
%resultCosts = run(suiteCosts);

%Create Suite from SolverTest Class Definition File
%The fromFile method creates a suite using the name of the file to identify the class.

%suiteFile = TestSuite.fromFile('SolverTest.m');
%result = run(suiteFile);