function [t,y] = rungekutter(n, s, f1, f2, y0, t0)

% Runge-Kutter-Verfahren zur numerischen Lösung von Differentialgleichungen
% Eigabewerte:
        % n: Anzahl der auszuführenden Zeitschritte
        % s: Groesse eines Zeitschrittes
        % f1: Matrix der Dynamik ausgewertet zu allen Zeitpunkten
        % f2: rechte Seite der Dynamik als Funktion
% Ausgabewerte:
        % t = Vektor aller Zeitschritte
        % y = errechneter Funktionswert zu jedem Zeitschritt


t = zeros(1,n+1);
y = zeros(17,n+1);

t(1) = t0;
y(1) = y0;

for j = 1:n;
    k1 = f1(:,j);
    k2 = feval(f2,y(j)+s/2.*k1);
    k3 = feval(f2,y(j)+s/2.*k2);
    k4 = feval(f2,y(j)+s.*k3);
    y(j+1) = y(j)+s/6.*(k1+2.*k2+2.*k3+k4);
    t(j+1) = t(j)+s;
end