%RUN_TEST Diese Datei starten die Tests f�r template.
%   Matlab besitzt sogenanntes unit test framework. Dies ist ein kleines
%   Beispiel dazu
%   Mehr Informationen gibt es unter: http://de.mathworks.com/help/matlab/matlab_prog/create-simple-test-suites.html
import matlab.unittest.TestSuite

global TEST;
TEST = true;


%Create Suite from SolverTest Class
%The fromClass method creates a suite from all Test methods in the SolverTest class.

suiteBasisQDyn = BasisQDyn();  %Check: BB: 05.07.2015
resultDyn = run(suiteBasisQDyn);


suiteode15sM = ode15sM();  %Check: BB: 05.07.2015
resultsuiteode15sM = run(suiteode15sM);


suiteInt = TestInt;  %Check: BB: 05.07.2015
resultInt = run(suiteInt);


suiteQuadrocopterExt = QuadrocopterExt(); %Check: BB: 05.07.2015
suiteQuadrocopterExt.steadyPoint = []; % SteadyPoint initialisieren.
resultQuadrocopterExt = run(suiteQuadrocopterExt, 'testWind'); % Nur Windfunktion testen

suiteMultiShooting = TestSuite.fromClass(?MultiShooting); %Check: BB: 05.07.2015
resultMultiShooting = run(suiteMultiShooting);

suiteConstraints = TestSuite.fromClass(?Constraints); %Check: BB: 05.07.2015
resultConstraints = run(suiteConstraints);

suiteCosts = TestSuite.fromClass(?Costs);  %Check: BB: 05.07.2015
resultCosts = run(suiteCosts);                           
 
suiteCostsXU = TestSuite.fromClass(?CostsXU); %Check: BB: 05.07.2015
resultCostsXU = run(suiteCostsXU);                        

suiteCostsComplete = TestSuite.fromClass(?CostsComplet); %Check: BB: 05.07.2015
resultCostsComplete = run(suiteCostsComplete);     

suiteRiccati = TestSuite.fromClass(?RiccatiManager); %Check: BB: 05.07.2015
resultRiccati = run(suiteRiccati); 

suiteRTSolver = TestSuite.fromClass(?RealtimeSolver); %Check: BB: 05.07.2015
resultRTSolver = run(suiteRTSolver);                        %Check

suiteLagrange = TestSuite.fromClass(?Lagrange); % Failed: testgetLD_Euler mit 3.080005295425193e-04  und testgetLD_ode15sM mit 1.489192913595154 kann aber auch sein, das die Zeitintervalle zu kurz sind
resultLagrange = run(suiteLagrange);


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
