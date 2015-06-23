%RUN_TEST Diese Datei starten die Tests f�r template.
%   Matlab besitzt sogenanntes unit test framework. Dies ist ein kleines
%   Beispiel dazu
%   Mehr Informationen gibt es unter: http://de.mathworks.com/help/matlab/matlab_prog/create-simple-test-suites.html
import matlab.unittest.TestSuite

global TEST;
TEST = true;


%Create Suite from SolverTest Class
%The fromClass method creates a suite from all Test methods in the SolverTest class.

% suiteBasisQDyn = TestSuite.fromClass(?BasisQDyn);
% resultDyn = run(suiteBasisQDyn);                          %Check

% suiteConstraints = TestSuite.fromClass(?Constraints);
% resultConstraints = run(suiteConstraints);                %Check

suiteMultiShooting = TestSuite.fromClass(?MultiShooting);
resultMultiShooting = run(suiteMultiShooting);                %Check
        
% suiteForwEuler = TestSuite.fromClass(?ForwEuler);
% resultForwEuler = run(suiteForwEuler);                    %NOT 

%suiteODE15iM2 = TestSuite.fromClass(?ode15iM2);
%resultODE15iM2 = run(suiteODE15iM2);                      %Check

% suiteCosts = TestSuite.fromClass(?Costs);
% resultCosts = run(suiteCosts);                            %Check

 
% suiteRiccati = TestSuite.fromClass(?RiccatiManager);
% resultRiccati = run(suiteRiccati);                        %Check

%suiteRTSolver = TestSuite.fromClass(?RealtimeSolver);
%resultRTSolver = run(suiteRTSolver);                        %Check

% 
%Create Suite from SolverTest Class Definition File
%The fromFile method creates a suite using the name of the file to identify the class.

%suiteFile = TestSuite.fromFile('SolverTest.m');
%result = run(suiteFile);