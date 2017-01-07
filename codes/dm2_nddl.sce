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

wi = sqrt(ki./mi); // pulsations modales carrees

B = eigenvscale(mi, V); // vecteurs propres divises par les masses modales associee a leur modes respectifs
// On peut verifier que Bt*M*B = matrice identite

// RVF
t = 0:0.1:3; // duree de l'excitation
F = sollicit(w,1,ne, t);

// Deformees
X = linspace(0, L, 50*ne+1); // intervalle [0,L]
// Deformee du mode 1
v1 = defmodale(X, B(:,1));
v2 = defmodale(X, B(:,2));
v3 = defmodale(X, B(:,3));
v4 = defmodale(X, B(:,4));
// Traces de ces modes
figure('figure_name','modes')
plot(X, v1, X, v2, X, v3, X, v4);
txtleg = ['mode 1, w = '+string(wi(1))+' rad/s',
'mode 2, w = '+string(wi(2))+' rad/s',
'mode 3, w = '+string(wi(3))+' rad/s',
'mode 4, w = '+string(wi(4))+' rad/s',
'mode 5, w = '+string(wi(1))+' rad/s'
];
legend(txtleg);
