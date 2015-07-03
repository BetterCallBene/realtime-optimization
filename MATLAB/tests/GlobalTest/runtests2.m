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

%suiteConstraints = TestSuite.fromClass(?Constraints);
%resultConstraints = run(suiteConstraints);                %Check


% suiteMultiShooting = TestSuite.fromClass(?MultiShooting);
% resultMultiShooting = run(suiteMultiShooting);                %Check
        
% suiteForwEuler = TestSuite.fromClass(?ForwEuler);
% resultForwEuler = run(suiteForwEuler);                    %NOT 

% suiteInt = TestInt;
% resultInt = run(suiteInt, 'testIntegratoren');

% suiteODE15sM = TestSuite.fromClass(?ode15sM);
% resultODE15sM = run(suiteODE15sM);                      %Check

suiteCosts = TestSuite.fromClass(?Costs);
resultCosts = run(suiteCosts);                            %Check
 
suiteCostsXU = TestSuite.fromClass(?CostsXU);
resultCostsXU = run(suiteCostsXU);                        %Check

suiteCostsComplete = TestSuite.fromClass(?CostsComplet);
resultCostsComplete = run(suiteCostsComplete);            %Check

% suiteLagrange = TestSuite.fromClass(?Lagrange);
% resultLagrange = run(suiteLagrange);

suiteRiccati = TestSuite.fromClass(?RiccatiManager);
resultRiccati = run(suiteRiccati);                        %Check

suiteRTSolver = TestSuite.fromClass(?RealtimeSolver);
resultRTSolver = run(suiteRTSolver);                        %Check

% 
%Create Suite from SolverTest Class Definition File
%The fromFile method creates a suite using the name of the file to identify the class.

%suiteFile = TestSuite.fromFile('SolverTest.m');
%result = run(suiteFile);

% If every test ran successfully, notify user
myvars  = who;
someflag = true;
for i = 1:length(myvars)
    tmp = myvars(i);
    tmp = tmp{1};
   if isa(eval(tmp),  'matlab.unittest.TestResult')
     if ~ eval([tmp, '.Passed'])
       someflag = false;
       break;
     end
   end
end

if someflag
    load handel
    sound(y,Fs);
end
