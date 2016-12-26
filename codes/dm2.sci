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

function fy = iterfsd(y) // F(y) de la methode de RK pour le cas avec frottement sec dynamique
    // y: vecteur colonne 2x1 : 
        // [position
        //vitesse]
    c = eta*c0 // coefficient de frottement visqueux
    M = [[0, 1];[-k/m, -c/m]];   
    fy = M*y - sign(y(2))*[0;tc/m]; // la force de frottement statique va s'opposer a la force de rappel et est donc de meme signe que le deplacement
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
