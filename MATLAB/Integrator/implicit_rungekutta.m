function [t, y] = implicit_rungekutta(k, s, f, y0, t0)

% implizites Runge-Kutter-Verfahren zur numerischen Lösung von Differentialgleichungen
% Eingabewerte:
        % n: Anzahl der auszuführenden Zeitschritte
        % s: Groesse eines Zeitschrittes
        % f: f(t,y), rechte Seite der Differentialgleichung
% Ausgabewerte:
        % t = Vektor aller Zeitschritte
        % y = Matrix der errechneten Funktionswerte zu jedem Zeitschritt


n = length(y0);
t = zeros(1,k+1);
y = zeros(n,k+1);

t(1) = t0;
y(:,1) = y0;

for j = 1:k;
    % k1 und k2 müssen noch über ein Gleichungssystem ausgerechnet werden
    % k1 = feval(f,t(j)+(1/2-sqrt(3)/6)*s,y(:,j)+s*(k1/4+(1/4-sqrt(3)/6)*k2));
    % k2 = feval(f,t(j)+(1/2+sqrt(3)/6)*s,y(:,j)+s*((1/4+sqrt(3)/6)*k1+k2/4));
    k = %Newtonverfahren(f,f_x,A) -> A: Koeffizienten a_ij zur Berechnung der k_i
    y(:,j+1) = y(:,j)+s/2*(k(1)+k(2));
    t(:,j+1) = t(j)+s;
end