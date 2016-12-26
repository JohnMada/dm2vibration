// DM2 vibrations
// Programme principal
// Yedhir MEZACHE
// Jonathan RAVAHIMANANA
// Hasinantenaina RAZAFIMAHALEO
clear;
clc;
//exec("dm2.sci",-1);
//xdel(winsid());

M = read('results/matrice_2d_m',2,2); // Matrices obtenue a l'aide de python
K = read('results/matrice_2d_k',2,2); //

[V, w] = spec(inv(M)*K); // resp. vecteur propre et valeur propre du probleme
