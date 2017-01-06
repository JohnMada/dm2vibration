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
rho = 5000; // kg/m3
L = 1; // m longueur
h = .05; // m hauteur de section
b = .05; // m base de section
E = 2e11; // N/m2 module d'Young de la poutre
I = b*h^3; // m4 moment d'inertie de section
EI = E*I;
S = b*h; // m2 surface de section
// Sollicitation
w = 5000; // 2*pi/s pulsation de l'excitation
// Disque
m = 8; // kg masse du disque
R = 0.15; // m rayon du disque
Iz = m*R^2; // moment d'inertie du disque autour de z
// Approximation
/// 2ddl
//


/// nddl
ne = 4; // Nombre d'elements. Doit etre pair pour prendre en compte le fait qu'il y a un disque au milieu
ndisque = ne/2; // Numero de l element qui porte le disque

dx = L/ne;

[K,M] = assemblage(ne); // Matrices elementaires
C = .5*K + .1*M; // Matrice d'amortissement visqueux,frottement proportionnel

// Analyse modale
[al, be, V] = spec(K,M); // val. propres generalisees

ki = vmodales(V, K); // vecteur des raideurs modales
ci = vmodales(V, C); // vecteur des amortissements modaux
mi = vmodales(V, M); // masses modales

wi = ki./mi; // pulsations modales

B = eigenvscale(mi, V); // vecteurs propres divises par les masses modales associee a leur modes respectifs
// On peut verifier que Bt*M*B = matrice identite

// RVF
t = 0:0.1:3; // duree de l'excitation
F = sollicit(w,1,ne, t);
