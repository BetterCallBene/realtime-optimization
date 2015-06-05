function [t, y] = rungekutta(k, s, f, y0, t0)

% Runge-Kutter-Verfahren zur numerischen Lösung von Differentialgleichungen
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
    k1 = feval(f,t,y(:,j));
    k2 = feval(f,t,y(:,j)+s/2.*k1);
    k3 = feval(f,t,y(:,j)+s/2.*k2);
    k4 = feval(f,t,y(:,j)+s.*k3);
    y(:,j+1) = y(:,j)+s/6.*(k1+2.*k2+2.*k3+k4);
    t(:,j+1) = t(j)+s;
end