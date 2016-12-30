function u=rvla(y0, t0, t1, Eta) // reponse en vibration libre amortie
    // y0: vecteur des CI
    // t0: temps initial, inutilise
    // t1: intervalle de temps de simulation
    // Eta: facteur d'amortissement visqueux
    // Retourne une matrice avec en premiere ligne le deplacement et en deuxieme ligne la vitesse
    c = Eta*c0 // coefficient de frottement visqueux
    function [Y]=fonction(t, V) // fonction a fournir a ODE pour resolution de l equation diff du second ordre regissant le systeme amorti
        M = [[0, 1];[-k/m, -c/m]] // matrice d iteration
        Y = M*V         
    endfunction
    [u]=ode(y0, t0, t1, fonction);
endfunction

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
    M(1,2) = 12*EI/dx^2; M(1,3) = -24*EI/dx^3; M(1,4) = 12*EI/dx^2;
    M(2,3) = -12*EI/dx^2; M(2,4) = 4*EI/dx;
    M(3,4) = -12*EI/dx^2;
        
    M = M+M';
    M(1,1) = 24*EI/dx^3; 
    M(2,2) = 8*EI/dx;
    M(3,3) = 24*EI/dx^3;
    M(4,4) = 8*EI/dx;  
endfunction

function y=rvlfs(y0, t0, t1, ts, tc, eta) // reponse en vibration avec frottement sec
    // y0: vecteur des CI
    // t0: temps initial, inutilise
    // t1: intervalle de temps de simulation
    // ts: force maximale du frottement statique
    // tc: force de frottement dynamique
    // Utilise la méthode RK4 pour calculer l'oscillation avec frottement sec
    // Retourne une matrice avec en premiere ligne le deplacement et en deuxieme ligne la vitesse
    dt = t1(2)-t1(1); // pas de temps
    n = length(t1); // nombre d'itérations
    y = zeros(2,n);
    y(:,1) = y0; // debut initialisation
    K1 = dt*iterfsd(y0);
    K2 = dt*iterfsd(y0+K1/2);
    K3 = dt*iterfsd(y0+K2/2);
    K4 = dt*iterfsd(y0+K3);
    y(:,2) = y0 + K1/6 + K2/3 + K3/3 + K4/6; // fin initialisation
    for i=2:n-1      
        yi = y(:,i); // contient y_i et dy_i/dt
        K1 = dt*iterfsd(yi);
        K2 = dt*iterfsd(yi+K1/2);
        K3 = dt*iterfsd(yi+K2/2);
        K4 = dt*iterfsd(yi+K3);
        y(:,i+1) = yi + K1/6 + K2/3 + K3/3 + K4/6;  
        if (yi(2)*y(2,i+1)<= 0) // vitesse nulle (ie traverse l'axe des abscisses)      
            if k*abs(yi(1)) < ts // frottement statique
//                K1 = dt*iterfss(yi);
//                K2 = dt*iterfss(yi+K1/2);
//                K3 = dt*iterfss(yi+K2/2);
//                K4 = dt*iterfss(yi+K3);               
                y(:,i+1) = [yi(1);0]; // arret du mouvement
            end
        end    
    end
endfunction
