Realtime Optimizaton Ordnerstruktur

MAPLE: ODE's des Quadrocoptermodells und deren Ableitungen
PYTHON: Konvertiert den MAPLE Inhalt (Quadrocopter Modell) in die MATLAB Struktur
MATLAB

	- template_projekt: Template Projektdateien
	- visualization: 
		- 3D Animation des Quadrocopters
		- Beispiele fuer Approximation, Skierfahrer, Horizon und Approximierung
	- tests: Tests fuer den Quadropter  (vor jedem Test das Skript init ausfuehren )
		- Wind: Simulation des Quadrocopter mit starkem Windeinfluss
		- realtimevsfmincon: Vergleich zwischen Realtime Ansatz und fmincon
		- GlobalTest: Test aller Komponenten mit Hilfe von dem Skript runtests
		- BasisKonfiguration: Alle Dateien die benoetigt werden, um einem neuen Test erstellen.
			Erstellen eines Test im tests Ordner -> Neuen Ordner erstellen -> Ordner BasisKonfigurationen in den Pfad MATLAB integrieren -> Befehl start in die CommandBox eingeben -> MAPLE Input Skripte in den neuen Ordner exportieren 
	- rtopt: Ordner der MATLAB Projektdateien
		- Realtime Ordner: Projektdateien speziell fuer den Realtime Ansatz.
			



