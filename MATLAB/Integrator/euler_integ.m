function [t, y] = euler_integ(k, s, f, y0, t0)

% Einfacher Vorwärts-Euler zur numerischen Loesung von
% Differentialgleichungen
% Eingabewerte:
        % k: Anzahl der Schritte
        % s: Schrittgroesse
        % f: f(t,y), rechte Seite der Differentialgleichung
        % y0: Anfangswert
        % t0: Anfangszeit
% Ausgabewerte:
        % t = Vektor aller Zeitschritte
        % y = Matrix der errechneten Funktionswerte zu jedem Zeitschritt        


n = length(y0);
y = zeros(n,k+1);
t = zeros(1,k+1);

y(:,1) = y0;
t(1) = t0;

for j = 1:k;
    y(:,j+1) = y(:,j) + s .* feval(f,t(j),y(:,j));
    t(j+1) = t(j) + s;
end
