%RUN_TEST Diese Datei starten die Tests f�r template.
%   Matlab besitzt sogenanntes unit test framework. Dies ist ein kleines
%   Beispiel dazu
%   Mehr Informationen gibt es unter: http://de.mathworks.com/help/matlab/matlab_prog/create-simple-test-suites.html
import matlab.unittest.TestSuite

global TEST;
TEST = true;


%Create Suite from SolverTest Class
%The fromClass method creates a suite from all Test methods in the SolverTest class.

%suiteBasisQDyn = BasisQDyn();
%resultDyn = run(suiteBasisQDyn);

%suiteForwEuler = ForwEuler();
%resultForwEuler = run(suiteForwEuler);


%suiteRungeKutta = RungeKutta();
%resultRungeKuttaK1 = run(suiteRungeKutta, 'testOde');
%resultRungeKuttaK2 = run(suiteRungeKutta, 'testk2');
%resultRungeKuttaK3 = run(suiteRungeKutta, 'testk3');
%resultRungeKuttaK4 = run(suiteRungeKutta, 'testk4');
%resulttestIntegrate = run(suiteRungeKutta, 'testIntegrate');

%suiteInt = TestInt;
%resultInt = run(suiteInt);

%suiteOdeM = odeM();
%resultOdeM = run(suiteOdeM);

suiteConstraints = TestSuite.fromClass(?Constraints);
resultConstraints = run(suiteConstraints);

%suiteMultiShooting = TestSuite.fromClass(?MultiShooting);
%resultMultiShooting = run(suiteMultiShooting);




%suiteCosts = TestSuite.fromClass(?Costs);
%resultCosts = run(suiteCosts);

%Create Suite from SolverTest Class Definition File
%The fromFile method creates a suite using the name of the file to identify the class.

%suiteFile = TestSuite.fromFile('SolverTest.m');
%result = run(suiteFile);