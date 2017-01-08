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

    //    K(1,:) = 0; M(1,:) = 0;
    //    K(:,1) = 0; M(:,1) = 0;
    //    K(1,1) = 1; M(1,1) = 1;
    //    // Deplacement nul sur l avant dernier ddl (V=0)
    //    K(nddl-1,:) = 0; M(nddl-1,:) = 0;
    //    K(:,nddl-1) = 0; M(:,nddl-1) = 0;
    //    K(nddl-1,nddl-1) = 1; M(nddl-1,nddl-1) = 1;       
endfunction

function f = sollicit(w,f0,ne, t)
    // Vecteur des forces harmoniues nodales
    // f : vecteur des forces nodales
    // w : pulsation de la sollicitation
    // f0 : amplitude initiale
    // ne : nombre d'elements
    // t : vecteur des temps
    nddl = 2*ne+2;
    f = zeros(nddl,length(t));
    f(nddl/2-1,:) = f0*cos(w*t);
    // Prise en compte des CL
    f = f([2:nddl-2, nddl],:);
    //
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
    mp = zeros(n,1);
    for i = 1:n
        mp(i) = V(:,i)'*M*V(:,i);
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

//function phix = base(X, ne, k, H)
//    // Trace des fonctions de bases dans l'espace physique [0,L]. Rotation
//    // ne : nombre d'elements
//    // dx : longueur d'un element
//    // X : vecteur des abscisses [0,L]
//    // k : numero de l'element concerne
//    dx = L/ne;
//    // passage entre la liste des noeuds et la liste des abscisses references dans x
//    npoints = length(X); // nombre d'abscisses dans x
//    nint = npoints - 1; // nombre d'intervalles dans x. Chaque intervalle separe deux points xk et xk+1
//    npint = nint/ne + 1; // nombre de points constituant un intervalle entre xi et xi+1
//
//    xi = linspace(-1, 1, npint); // xi dans l'espace de reference
//    xk1 = X((npint-1)*(k+1) + 2 - npint); // abscisse de xk+1
//    xk = X((npint-1)*k + 2 - npint); // abscisse de xk
//    x = linspace(xk, xk1, npint); //dx*xi/2 + (xk1 + xk)/2; // x dans l'espace physique [xk, xk+1]
//
//    phix = zeros(1,npoints); // fonction de basse continue par morceau et a support dans [xi,xi+1]
//    phix((npint-1)*k + 2 - npint:(npint-1)*(k+1) + 2 - npint) = dx*H(xi); // fonction de forme dans l'espace physique
//endfunction

function vx = defmodale(X, V)
    // Deformation modale sur [0,L]
    // X : vecteur abscisse dans [0,L]
    // V : vecteur propres
    // ne : nombre d'elements
    // Sortie :
    // vx : deformation modale de la poutre  
    N = length(V);
    nddl = N+2; // nombre totale de degres de libertes, CL de Dirichlet incluses
    nel = nddl/2 - 1; // nombre d'elements

    npoints = length(X); // nombre d'abscisses dans x
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
        vx((npint-1)*i + 2 - npint:(npint-1)*(i+1) + 2 - npint) = v(idg1)*base(X, nel, i, n1) + v(idg2)*base(X, nel, i, h1) + v(idg3)*base(X, nel, i, n2) + v(idg4)*base(X, nel, i, h2);
    end
endfunction
