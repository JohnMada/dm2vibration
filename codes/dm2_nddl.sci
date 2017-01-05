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
     
    // Matrice de masse
    M(1,2) = 0.0523809523809524*S*dx^2*rho; M(1,3) = 0.128571429*S*dx*rho ; M(1,4) = -.030952380952381*S*dx^2*rho;
    M(2,3) = .030952380952381*S*dx^2*rho; M(2,4) = -.00714285714285714*S*dx^2*rho;
    M(3,4) = -.0523809523809524*S*dx^2*rho;
        
    M = M+M';
    M(1,1) = 0.37142857142857*S*dx*rho; 
    M(2,2) = .00952380952380952*S*dx^3*rho;
    M(3,3) = .371428571428571*S*dx*rho;
    M(4,4) = .00952380952380952*S*dx^3*rho;  
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
    med = me; med(3,3) = med(3,3) + m
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
