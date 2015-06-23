%RUN_TEST Diese Datei starten die Tests f�r template.
%   Matlab besitzt sogenanntes unit test framework. Dies ist ein kleines
%   Beispiel dazu
%   Mehr Informationen gibt es unter: http://de.mathworks.com/help/matlab/matlab_prog/create-simple-test-suites.html

import matlab.unittest.TestSuite

global TEST;
TEST = true;

clc;


%Create Suite from SolverTest Class
%The fromClass method creates a suite from all Test methods in the SolverTest class.

% suiteBasisQDyn = TestSuite.fromClass(?BasisQDyn);
% resultDyn = run(suiteBasisQDyn);

% suiteConstraints = TestSuite.fromClass(?Constraints);
% resultConstraints = run(suiteConstraints);

% suiteForwEuler = TestSuite.fromClass(?ForwEuler);
% resultForwEuler = run(suiteForwEuler);
% 
% suiteCosts = TestSuite.fromClass(?Costs);
% resultCosts = run(suiteCosts);
% 
% suiteRiccati = TestSuite.fromClass(?RiccatiManager);
% resultRiccati = run(suiteRiccati);

suiteRTSolver = TestSuite.fromClass(?RealtimeSolver);
resultRTSolver = run(suiteRTSolver);

%Create Suite from SolverTest Class Definition File
%The fromFile method creates a suite using the name of the file to identify the class.

%suiteFile = TestSuite.fromFile('SolverTest.m');
%result = run(suiteFile);