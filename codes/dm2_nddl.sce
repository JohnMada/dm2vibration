// DM2 vibrations
// Programme principal
// Yedhir MEZACHE
// Jonathan RAVAHIMANANA
// Hasinantenaina RAZAFIMAHALEO
clear;
clc;
exec("dm2_nddl.sci",-1);
xdel(winsid());

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
m = 20; // kg masse du disque
R = 0.15; // m rayon du disque
Iz = m*R^2; // moment d'inertie du disque autour de z
// Approximation
/// 2ddl
//


/// nddl
ne = 8; // Nombre d'elements. Doit etre pair pour prendre en compte le fait qu'il y a un disque au milieu
ndisque = ne/2; // Numero de l element qui porte le disque

dx = L/ne;

[K,M] = assemblage(ne); // Matrices elementaires
c1 = .1; c2 = .5 
C = c1*K + c2*M; // Matrice d'amortissement visqueux,frottement proportionnel

// Analyse modale
[al, be, V] = spec(K,M); // val. propres generalisees
V = real(V); // Etrangement, scilab considere ces vecteurs comme complexe meme si leur parties imaginaires sont toutes nulles

ki = vmodales(V, K); // vecteur des raideurs modales
ci = vmodales(V, C); // vecteur des amortissements modaux
mi = vmodales(V, M); // masses modales

wi = sqrt(ki./mi); // pulsations modales

B = eigenvscale(mi, V); // vecteurs propres divises par les masses modales associee a leur modes respectifs
// On peut verifier que Bt*M*B = matrice identite

// RVF
dt = .1;
t = 0:dt:5; // duree de l'excitation

// Deformees
X = linspace(0, L, 100*ne+1); // intervalle [0,L]
// Deformee du mode 1
v1 = defmodale(X, B(:,$));
v2 = defmodale(X, B(:,$-1));
v3 = defmodale(X, B(:,$-2));
v4 = defmodale(X, B(:,$-3));
v5 = defmodale(X, B(:,$-4));
// Traces de ces modes
// fonctions de base
b1 = []; b2 = [], b3 = [], b4 = [];
x = [];
for i=1:ne
    b1 = [b1,base(X, ne, i, n1)];
end
for i=1:ne
    b3 = [b3, base(X, ne, i, n2)];
end
for i=1:ne
    b2 = [b2,base(X, ne, i, h1)];
end
for i=1:ne
    b4 = [b4,base(X, ne, i, h2)];
    x = [x, linspace((i-1)*dx,(i-1)*dx + dx,length(base(X, ne, i, h2)))]
end
figure("figure_name","bases",'BackgroundColor',[1,1,1])
plot(x,b1, x,b2, x,b3, x,b4)


figure('figure_name','modes symétriques','BackgroundColor',[1,1,1])
title('Déformées modales, modes symétriques');
plot(X, v1, X, v3, X, v5);
txtlegs = [
'mode 1, w = '+string(wi($))+' rad/s',
'mode 3, w = '+string(wi($-2))+' rad/s',
'mode 5, w = '+string(wi($-4))+' rad/s'
];
legend(txtlegs);

xsave("results/deformees_sym.png");

figure('figure_name','modes anti-symétriques','BackgroundColor',[1,1,1])
title('Déformées modales, modes anti-symétrique');
plot(X, v2, X, v4);
txtlegas = [
'mode 2, w = '+string(wi($-1))+' rad/s',
'mode 4, w = '+string(wi($-3))+' rad/s'
];
legend(txtlegas);
xsave("results/deformees_asym.png");

// Vibration forcee
// Resolution numerique par une methode iterative
w = 1510; // pulsation forcee
f0 = 2000; // N amplitude maximale de la force
F = sollicit(f0, ne);
Fp = reform(F,ne);

q = repmodale(w, wi, ci, Fp, t);
v = B*q;

cnp = diag(B'*C*B); // amortissements modales normalises au sens de la masse
dnp = diag(B'*K*B); // pulsation modales carrees normalisees au sens de la masse

//t = 0:0.1:5; // duree de simulation
//V0 = zeros(length(cnp),1); // La poutre n'est pas deformee au debut de la sollicitation
//Q0 = reform(V0, ne);
//
//A = matiter(cnp, dnp, ne);
//Q = rvf(cnp, dnp, Q0, w, Fp, ne, t);
//V = B*Q(1:2:$,:);

figure("figure_name","reponse forcee",'BackgroundColor',[1,1,1])
plot(t,v(ndisque-1,:), t, v(ndisque,:)); // amplitude de deplacement du disque
title("Déplacement vertical du disque au cours du temps, w = '+string(w)+' rad/s')
xlabel('t [s]')
ylabel('v(x'+string(ndisque - 1) + ',t) [m], theta(x'+string(ndisque-1) + ',t) [rad]')
legend('translation','rotation')

temps = [10, 25, 35, 50]; // *dx
couleur = ['k' 'b' 'r' '--'];
legende = []
figure("figure_name","poutre",'BackgroundColor',[1,1,1])
for i = 1:length(temps)
    ti = temps(i);
    ci = couleur(i)
    plot(X, defmodale(X, v(:,ti)), ci)
    legende = [legende, 't = '+string(ti*dt)+' s'];
end
//plot(X, defmodale(X, v(:,10)), X, defmodale(X, v(:,25)), X, defmodale(X, v(:,50)))
title('Déformation de la poutre au cours du temps, w = '+string(w)+' rad/s')
xlabel('x [m]')
ylabel('v(x, t) [m]')
legend(legende)

// w = w2
Fp = inv(B)*F;
w = wi($-1)*0.85;
q = repmodale(w, wi, ci, Fp, t);
v = B*q; // retour dans l'espace physique

figure("figure_name","reponse forcee",'BackgroundColor',[1,1,1])
plot(t,v(ndisque-1,:), t, v(ndisque,:)); // amplitude de deplacement du disque
title("Réponse du disque au cours du temps, w = '+string(w)+' rad/s')
xlabel('t [s]')
ylabel('v(x'+string(ndisque - 1) + ',t) [m], theta(x'+string(ndisque-1) + ',t) [rad]')
legend('translation','rotation')

figure("figure_name","poutre",'BackgroundColor',[1,1,1])
for i = 1:length(temps)
    ti = temps(i);
    ci = couleur(i)
    plot(X, defmodale(X, v(:,ti)), ci)
    legende = [legende, 't = '+string(ti*dt)+' s'];
end
//plot(X, defmodale(X, v(:,10)), X, defmodale(X, v(:,25)), X, defmodale(X, v(:,50)))
title('Déformation de la poutre au cours du temps, w = '+string(w)+' rad/s')
xlabel('x [m]')
ylabel('v(x, t) [m]')
legend(legende)
