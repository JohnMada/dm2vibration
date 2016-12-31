// DM2 vibrations
// Programme principal
// Yedhir MEZACHE
// Jonathan RAVAHIMANANA
// Hasinantenaina RAZAFIMAHALEO
clear;
clc;
exec("dm2_nddl.sci",-1);
//xdel(winsid());

//Donnees du probleme
// Poutre
rho = 1000; // kg/m3
L = .5; // m longueur*
h = .02; // m hauteur de section
b = .01; // m base de section
E = 2e11; // N/m2 module d'Young de la poutre
I = b*h^3; // m4 moment d'inertie de section
EI = E*I;
S = b*h; // m2 surface de section
// Disque
m = 5; // kg masse du disque
R = 0.06; // m rayon du disque
Iz = m*R^2; // moment d'inertie du disqueautour de z
// Approximation
ne = 4; // Nombre d'elements. Doit etre pair pour prendre en compte le fait qu'il y a un disque au milieu
ndisque = ne/2; // Numero de l element qui porte le disque

dx = L/ne;

[K,M] = matrices(EI,dx,rho); // Matrices elementaires

[V, w] = spec(inv(M)*K); // resp. vecteur propre et valeur propre du probleme
