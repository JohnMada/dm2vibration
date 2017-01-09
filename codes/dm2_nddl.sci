function [K,M] = matrices(EI,dx,rho)
    // Construction des matrices de masse et de raideur elementaires 
    K = zeros(4,4);
    M = zeros(4,4);
    // Matrice de raideur
    K(1,2) = 12*EI/dx^2; K(1,3) = -24*EI/dx^3; K(1,4) = 12*EI/dx^2;
    K(2,3) = -12*EI/dx^2; K(2,4) = 4*EI/dx;
    K(3,4) = -12*EI/dx^2;

    K = K+K';
    K(1,1) = 24*EI/dx^3; 
    K(2,2) = 8*EI/dx;
    K(3,3) = 24*EI/dx^3;
    K(4,4) = 8*EI/dx;
    K = K/2;

    // Matrice de masse
    M(1,2) = 11/6*dx; M(1,3) = 9/2 ; M(1,4) = -13/12*dx;
    M(2,3) = 13/12*dx; M(2,4) = -1/4*dx^2;
    M(3,4) = -11/6*dx;

    M = M+M';
    M(1,1) = 13; 
    M(2,2) = 1/3*dx^2;
    M(3,3) = 13;
    M(4,4) = 1/3*dx^2;

    M = M*rho*S*dx/35;  
endfunction

function tab = tbc(Ne)
    // table de connexion
    // Ne : nombre d'elements
    tab = zeros(Ne,4); // 4 colonnes correspondant chacune à un degré de liberté dans l'espace de reference et Ne lignes correspondant chacun a un element
    for i = 1:Ne
        for j = 1:4
            tab(i,j) = j+2*(i-1)
        end
    end
endfunction

function [K, M]=assemblage(Ne)
    // fonction d'assemblage des matrices elementaires en une matrice globale
    // Ne : nombre d'elements
    nddl = 2*Ne + 2; // nombre total de ddl
    dx = L/Ne; // longueur d'un element
    [ke, me] = matrices(EI,dx,rho)
    tab = tbc(Ne); // table de connexion
    // assemblage
    K = zeros(nddl,nddl); M = zeros(nddl,nddl);
    for i=1:Ne/2-1 // elements avant l'element porteur
        for j=1:4
            for l=1:4 // j et l = index dans la matrice elementaire
                idg1 = tab(i,j); idg2 = tab(i,l); // idgi = index dans la matrice globale
                K(idg1,idg2) = K(idg1,idg2) + ke(j,l);
                M(idg1,idg2) = M(idg1,idg2) + me(j,l); 
I = b*h^3; // m4 moment d'inertie de section
EI = E*I;
S = b*h; // m2 surface de section
// Sollicitation
w = 5000; // 2*pi/s pulsation de l'excitation
// Disqueg2) = M(idg1,idg2) + me(j,l);
            end
        end
    end
    // element qui porte le disque
    // i = Ne/2
    med = me; med(3,3) = med(3,3) + m;
    for j=1:4
        for l=1:4 // j et l = index dans la matrice elementaire
            idg1 = tab(Ne/2,j); idg2 = tab(Ne/2,l); // idgi = index dans la matrice globale
            K(idg1,idg2) = K(idg1,idg2) + ke(j,l);
            M(idg1,idg2) = M(idg1,idg2) + med(j,l);
        end
    end
    for i=Ne/2+1:Ne // elements apres l'element porteur
        for j=1:4
            for l=1:4 // j et l = index dans la matrice elementaire
                idg1 = tab(i,j); idg2 = tab(i,l); // idgi = index dans la matrice globale
                K(idg1,idg2) = K(idg1,idg2) + ke(j,l);
                M(idg1,idg2) = M(idg1,idg2) + me(j,l);
            end
        end
    end
    //// Conditions limites
    // Deplacement nul sur le ddl 1 (V1=0)
    // Deplacement nul sur l avant dernier ddl (V_(nddl-1)=0)
    K = K([2:nddl-2,nddl],[2:nddl-2,nddl]);
    M = M([2:nddl-2,nddl],[2:nddl-2,nddl]);
     
endfunction

function W = eigenvscale(c, V)
    // c vecteur de n lignes
    // V matrice des n vecteurs propres du systeme
    n = length(c); // ordre du systeme
    W = zeros(n,n); // Matrice des vecteurs propres normalises
    for j = 1:n
        W(:,j) = V(:,j)/sqrt(c(j));
    end
endfunction

function mp = vmodales(V, M) // Carre scalaire d'un vecteur propre au sens de M
    // mp : vecteur colonne des masses modales du systeme
    // V : matrice des vecteurs propres
    // M : matrice de masse
    // Marche pour K, M et C -> analyse modale
    n = size(V,1); // ordre du systeme
    mp = [];
    for i = 1:n
        mp = [mp;V(:,i)'*M*V(:,i)];
    end
endfunction

// Fonctions de forme
function phixi = n1(xi)
    phixi = (1-xi).^2 .* (0.5 + xi/4);
endfunction
function phixi = n2(xi)
    phixi = (1+xi).^2 .* (0.5 - xi/4);
endfunction
function psyxi = h1(xi)
    psyxi = .25*(1-xi).^2 .* (1+xi);
endfunction
function psyxi = h2(xi)
    psyxi = -.25*(1+xi).^2 .* (1-xi);
endfunction

function phix = base(X, ne, k, Phi)
    // Trace des fonctions de bases dans l'espace physique [0,L]. Translation
    // ne : nombre d'elements
    // dx : longueur d'un element
    // X : vecteur des abscisses [0,L]
    // k : numero de l'element concerne
    dx = L/ne;
    // passage entre la liste des noeuds et la liste des abscisses references dans x
    npoints = length(X); // nombre d'abscisses dans x
    nint = npoints - 1; // nombre d'intervalles dans x. Chaque intervalle separe deux points xk et xk+1
    npint = nint/ne + 1; // nombre de points constituant un intervalle entre xi et xi+1

    xi = linspace(-1, 1, npint); // xi dans l'espace de reference
    xk1 = X((npint-1)*(k+1) + 2 - npint); // abscisse de xk+1
    xk = X((npint-1)*k + 2 - npint); // abscisse de xk
    x = linspace(xk, xk1, npint); //dx*xi/2 + (xk1 + xk)/2; // x dans l'espace physique [xk, xk+1]

    //phix = zeros(1,npoints); // fonction de basse continue par morceau et a support dans [xi,xi+1]
    phix = Phi(xi); // fonction de forme dans l'espace physique
endfunction

function vx = defmodale(X, V)
    // Deformation modale sur [0,L]
    // X : vecteur abscpsse dans [0,L]
    // V : vecteur propres
    // ne : nombre d'elements
    // Sortie :
    // vx : deformation modale de la poutre  
    N = length(V);
    nddl = N+2; // nombre totale de degres de libertes, CL de Dirichlet incluses
    nel = nddl/2 - 1; // nombre d'elements

    npoints = length(X); // nombre d'abscpsses dans x
    nint = npoints - 1; // nombre d'intervalles dans x. Chaque intervalle separe deux points xk et xk+1
    npint = nint/ne + 1; // nombre de points constituant un intervalle entre xi et xi+1

    tab = tbc(nel); // table de connexions
    vx = zeros(1,npoints);

    v = [0; V(1:N-1); 0; V(N)]; // remise des degres de libertes enleves (V1 = 0 et Vnddl = 0)
    for i=1:nel
        idg1 = tab(i,1); // indice globale du V à gauche de l'element i
        idg2 = tab(i,2); // indice globale du theta à gauche de l'element i
        idg3 = tab(i,3); // indice globale du V à droite de l'element i
        idg4 = tab(i,4); // indice globale du theta à droite de l'element i
        vx((npint-1)*i + 2 - npint:(npint-1)*(i+1) + 2 - npint) = v(idg1)*base(X, nel, i, n1) + dx*v(idg2)*base(X, nel, i, h1) + v(idg3)*base(X, nel, i, n2) + dx*v(idg4)*base(X, nel, i, h2);
    end
endfunction

function f = sollicit(f0,ne)
    // Vecteur des forces harmoniques nodales
    // f : vecteur des forces nodales
    // w : pulsation de la sollicitation
    // f0 : amplitude initiale
    // ne : nombre d'elements
    // t : vecteur des temps
    nddl = 2*ne+2;
    f = zeros(nddl,1);
    f(nddl/2-1,:) = f0; // La force ne s'applique que sur la masse
    // Prise en compte des CL
    f = f([2:nddl-2, nddl]);
    //
endfunction

function [syst] = matiter(cp, dp, ne)
    // Retourne la matrice creuse d'iteration pour la resolution numerique de la reponse forcee
    // cp : vecteur des coefficients d amortissement normalises
    // dp : vecteur des pulsations carrees normalises
    nm = length(cp); // nombre de modes
    
    dinf = zeros(2*nm-1,1); // diagonale inferieure de la matrice
    dia = zeros(2*nm,1);
    dsup = repmat([1; 0],nm-1, 1); // diagonale superieure
    dsup = [dsup; 1]; // 1 0 1 0 ...ite 1 0 1 0 ... 1 0 1 0 1

    dinf(1:2:2*nm) = dp; // d1 0 d2 0 ... 0 dp 0 ... 0 dnm
    dia(2:2:2*nm) = cp; // 0 c1 0 ... 0 cp 0 ... 0 cnm

    syst = -diag(dia) - diag(dinf, -1) + diag(dsup, 1);
    //syst = sparse(syst);
endfunction

function Q = reform(q, ne)
    // transforme un vecteur pour l'adapter au systeme d'equation differentielle numerique
    nm = 2*ne + 2 - 2; // nombre de ddl
    q = inv(B)*q; // passage dans l'espace normal
    Q = zeros(2*nm,1);
    Q(2:2:$) = q; 
endfunction

function qdot = iter(t, q)
    // q: vecteur des deplacements modaux dans la base normale
    qdot = A*q + fp*cos(w*t);
endfunction

function Q = rvf(cp, dp, Q0, w, Fp, ne, t)
    // Q0 : conditions initiales dans la base normale
    aa = matiter(cp, dp, ne);
    Q = [];
    for i=1:2:2*ne+1
        fp = Fp(i:i+1);
        A = aa([i,i+1],[i,i+1]);
        qp = ode(Q0(i:i+1), 0, t, iter);
        Q = [Q ; qp];
    end
endfunction

function xp = repmodale(w, wi, ci, fp, t)
    n = length(wi);
    xp = [];
    for i = 1:n
        gammai = (c1 + c2*wi(i)^2)/2/wi(i);
        betai = 1/sqrt((1-w^2/wi(i)^2)^2 + (2*gammai));
        phii = atan(2*gammai*(w/wi(i))/(1 - w^2/wi(i)^2));
        xp = [xp; fp(i)/wi(i)^2*betai*cos(w*t - phii)];
    end   
endfunction
