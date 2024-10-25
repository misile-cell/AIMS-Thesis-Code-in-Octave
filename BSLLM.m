clear % Clear all variables on workspace
clc   % Clear screen
%close all;

Nt = 5;
Nx = 60;  % Collocation points in \eta and \xi
[D,y] = cheb(Nx); % D is chebyshev diff matrix in \eta
[d,Z] = cheb(Nt); % d is chebyshev diff matrix in \xi

a = 0;
b = 20;
Lx = b - a;

% Time interval
t0 = 0;
tf = 20;
Lt = tf - t0;

% Scaled differentiation matrices
D1 = (2/Lx)*D;
D2 = D1*D1;
D3 = D1^3;
eta = (Lx*y)/2 + (a+b)/2;  % x in [0,Lx]
I  = eye(Nx+1,Nx+1); % Identity matrix

scale_t = 2/(Lt);

% Scaled Differentiation Matrix in t
d1 = scale_t*d;

Xi = (Lt*Z)/2 + (tf + t0)/2;

% Constants
n = 1/2;
w = 1/2;
Pr = 0.7;
Sc = 0.6;

%-------------------------------------------------------------------------
%  Initial guess
%-------------------------------------------------------------------------
[fNt, gNt, hNt] = Initial_Guess_3_System(Nx,Lx,n,Pr,Sc,w);
%-------------------------------------------------------------------------


% Assume that the initial approximation is fNt for all t
F = repmat(fNt,1,Nt+1);
G = repmat(gNt,1,Nt+1);
H = repmat(hNt,1,Nt+1);
%tic


iterations = 20;
tic
%while Tol > Tolerance_level   %for r = 1,2,... iterations - 1
for r = 1:iterations
    
    
    % ----------------------------------------------------------------------
    % SOLVE FOR F
    %----------------------------------------------------------------------
    
    for i = 1:Nt     % for i = 0,1,2,...,Nt
        % For convenience of notation and to ensure that notation in the code
        % Is consistent with that in the write-up we define
        fprev(:,i) = F(:,i);
        
        fr = F(:,i);
        fr1 = D1*F(:,i);
        fr2 = D2*F(:,i);
        fr3 = D3*F(:,i);
        
        gr = G(:,i);
        gr1 = D1*G(:,i);
        gr2 = D2*G(:,i);
        
        hr = H(:,i);
        hr1 = D1*H(:,i);
        hr2 = D2*H(:,i);
        
        
        xi = Xi(i);
        
        sumf = 0;
        sumfprime = 0;
        sumg = 0;
        sumh = 0;
        
        for j = 1:Nt+1
            sumf = sumf + d1(i,j)*F(:,j);
            sumfprime = sumfprime + d1(i,j)*(D1*F(:,j));
            sumg = sumg + d1(i,j)*G(:,j);
            sumh = sumh + d1(i,j)*H(:,j);
        end
        
        
        
        % Define the coefficients
        a0r = 1;
        a1r = (1/4)*(n + 3)*fr + xi + (1/4)*(1 - n)*xi*sumf;
        a2r = -(n + 1)*fr1 - (1/4)*(1 - n)*xi*sumfprime;
        a3r = (1/4)*(n + 3)*fr2;
        a4r = (1/4)*(1 - n)*xi*fr2;
        a5r = -(1/4)*(1 - n)*xi*fr1;
        
        % Define the right hand side
        Ff = fr3 + (1/4)*(n + 3)*fr.*fr2 - (1/2)*(n + 1)*(fr1).^2 + xi*fr2 + (1 - w)*gr + w*hr - (1/4)*(1 - n)*xi*fr1.*sumfprime + (1/4)*(1 - n)*xi*fr2.*sumf;
        Rf = a0r.*fr3 + a1r.*fr2 + a2r.*fr1 + a3r.*fr + a4r.*sumf + a5r.*sumfprime - Ff;
        R1 = Rf - a4r.*(d1(i,Nt+1).*fNt) - a5r.*(d1(i,Nt+1).*(D1*fNt));
        
        % Define matrices
        Ai = diag(a0r)*D3 + diag(a1r)*D2 + diag(a2r)*D1 + diag(a3r);
        Aii = Ai + diag(a4r)*d1(i,i)*I + diag(a5r)*d1(i,i)*D1; % Matrix on main diagonal where i = j
        
        % Impose boundary conditions on the matrix Aii
        
        % f'(\inf) = 0
        Aii(1,:) = D1(1,:);
        R1(1) = 0;
        
        % f(0) = 0
        Aii(Nx+1,:) = 0;
        Aii(Nx+1,Nx+1) = 1;
        R1(Nx+1) = 0;
        
        % f'(0) = 0
        Aii(Nx,:) = D1(Nx+1,:);
        R1(Nx) = 0;
        
        
        % Form the big matrix (call it AA)
        for j = 1:Nt
            % Define matrices off the main diagonal where i NOT j
            Aij = diag(a4r)*d1(i,j)*I + diag(a5r)*d1(i,j)*D1;
            % Impose boundary conditions on each matrix
            Aij(1,:) = 0;
            Aij(Nx,:) = 0;
            Aij(Nx+1,:) = 0;
            % Put the off diagonal matrices on the big matrix AA
            AA((Nx+1)*(i-1)+1:(Nx+1)*i,(Nx+1)*(j-1)+1:(Nx+1)*j) = Aij;
        end
        % Put the Aii matrices on the main diagonal.
        AA((Nx+1)*(i-1)+1:(Nx+1)*i,(Nx+1)*(i-1)+1:(Nx+1)*i) = Aii;
        RRR1(:,i) = R1; % Put R(i) on the columns
        
    end % End of Nt iterations
    
    RR1 = RRR1(:); % convert into columns vector
    Ftemp = AA\RR1;
    
    for i = 1:Nt
        F(:,i) = Ftemp((i-1)*(Nx+1)+1:i*(Nx+1)); % Solution for f at each time level
    end
    
    
    %------------------------------------------------------------------------
    % SOLVE FOR G
    %------------------------------------------------------------------------
    for i = 1:Nt     % for i = 0,1,2,...,Nt
        % For convenience of notation and to ensure that notation in the code
        % Is consistent with that in the write-up we define
        
        gprev(:,i) = G(:,i);
        
        fr = F(:,i);
        fr1 = D1*F(:,i);
        fr2 = D2*F(:,i);
        
        gr = G(:,i);
        gr1 = D1*G(:,i);
        gr2 = D2*G(:,i);
        
        hr = H(:,i);
        hr1 = D1*H(:,i);
        hr2 = D2*H(:,i);
        
        
        xi = Xi(i);
        
        sumg = 0;
        sumf = 0;
        sumh = 0;
        
        for j = 1:Nt+1
            sumf = sumf + d1(i,j)*F(:,j);
            sumg = sumg + d1(i,j)*G(:,j);
            sumh = sumh + d1(i,j)*H(:,j);
        end
        
        % Define the coefficients
        b0r = 1/Pr;
        b1r = (1/4)*(n + 3)*fr + xi + (1/4)*(1 - n)*xi*sumf;
        b2r = 0;
        b3r = -(1/4)*(1 - n)*xi*fr1;
        
        
        % Define the right hand side
        Gg = (gr2/Pr) + (1/4)*(n + 3)*fr.*gr1 + xi*gr1 - (1/4)*(1 - n)*xi*fr1.*sumg + (1/4)*(1 - n)*xi*gr1.*sumf;
        Rg = b0r.*gr2 + b1r.*gr1 + b2r.*gr +  b3r.*sumg - Gg;
        R2 = Rg - b3r.*(d1(i,Nt+1).*gNt);
        
        
        % Define matrices
        Bi = diag(b0r)*D2 + diag(b1r)*D1 + diag(b2r);
        Bii = Bi + diag(b3r)*d1(i,i)*I; % Matrix on main diagonal where i = j
        
        % Impose boundary conditions on the matrix Aii
        
        % g(\inf) = 0
        Bii(1,:) = 0;
        Bii(1,1) = 1;
        R2(1) = 0;
        
        % g(0) = 1
        Bii(Nx+1,:) = 0;
        Bii(Nx+1,Nx+1) = 1;
        R2(Nx+1) = 1;
        
        
        % Form the big matrix (call it AA)
        for j = 1:Nt
            % Define matrices off the main diagonal where i NOT j
            Bij = diag(b3r)*d1(i,j)*I;
            % Impose boundary conditions on each matrix
            Bij(1,:) = 0;
            Bij(Nx+1,:) = 0;
            % Put the off diagonal matrices on the big matrix AA
            BB((Nx+1)*(i-1)+1:(Nx+1)*i,(Nx+1)*(j-1)+1:(Nx+1)*j) = Bij;
        end
        % Put the Aii matrices on the main diagonal.
        BB((Nx+1)*(i-1)+1:(Nx+1)*i,(Nx+1)*(i-1)+1:(Nx+1)*i) = Bii;
        
        RRR2(:,i) = R2; % Put R(i) on the columns
    end % End of Nt iterations
    RR2 = RRR2(:); % convert into columns vector
    Gtemp = BB\RR2;
    
    for i = 1:Nt
        G(:,i) = Gtemp((i-1)*(Nx+1)+1:i*(Nx+1)); % Solution for f at each time level
    end
    
    
    
    
    %------------------------------------------------------------------------
    % SOLVE FOR H
    %------------------------------------------------------------------------
    for i = 1:Nt     % for i = 0,1,2,...,Nt
        % For convenience of notation and to ensure that notation in the code
        % Is consistent with that in the write-up we define
        hprev(:,i) = H(:,i);
        
        fr = F(:,i);
        fr1 = D1*F(:,i);
        fr2 = D2*F(:,i);
        
        gr = G(:,i);
        gr1 = D1*G(:,i);
        gr2 = D2*G(:,i);
        
        hr = H(:,i);
        hr1 = D1*H(:,i);
        hr2 = D2*H(:,i);
        
        
        xi = Xi(i);
        
        sumg = 0;
        sumf = 0;
        sumh = 0;
        
        for j = 1:Nt+1
            sumf = sumf + d1(i,j)*F(:,j);
            sumg = sumg + d1(i,j)*G(:,j);
            sumh = sumh + d1(i,j)*H(:,j);
        end
        
        % Define the coefficients
        c0r = 1/Sc;
        c1r = (1/4)*(n + 3)*fr + xi + (1/4)*(1 - n)*xi*sumf;
        c2r = 0;
        c3r = -(1/4)*(1 - n)*xi*fr1;
        
        Hh = (hr2/Sc) + (1/4)*(n + 3)*fr.*hr1 + xi*hr1 - (1/4)*(1 - n)*xi*fr1.*sumh + (1/4)*(1 - n)*xi*hr1.*sumf;
        Rh = c0r.*hr2 + c1r.*hr1 + c2r.*hr + c3r.*sumh  - Hh;
        R3 = Rh - c3r.*(d1(i,Nt+1).*hNt);
        
        % Define matrices
        Ci = diag(c0r)*D2 + diag(c1r)*D1 + diag(c2r);
        Cii = Ci + diag(c3r)*d1(i,i)*I; % Matrix on main diagonal where i = j
        
        % Impose boundary conditions on the matrix Aii
        
        % h(\inf) = 0
        Cii(1,:) = 0;
        Cii(1,1) = 1;
        R3(1) = 0;
        
        % h(0) = 1
        Cii(Nx+1,:) = 0;
        Cii(Nx+1,Nx+1) = 1;
        R3(Nx+1) = 1;
        
        
        % Form the big matrix (call it AA)
        for j = 1:Nt
            % Define matrices off the main diagonal where i NOT j
            Cij = diag(c3r)*d1(i,j)*I;
            % Impose boundary conditions on each matrix
            Cij(1,:) = 0;
            Cij(Nx+1,:) = 0;
            % Put the off diagonal matrices on the big matrix AA
            CC((Nx+1)*(i-1)+1:(Nx+1)*i,(Nx+1)*(j-1)+1:(Nx+1)*j) = Cij;
        end
        % Put the Aii matrices on the main diagonal.
        CC((Nx+1)*(i-1)+1:(Nx+1)*i,(Nx+1)*(i-1)+1:(Nx+1)*i) = Cii;
        RRR3(:,i) = R3; % Put R(i) on the columns
    end % End of Nt iterations
    
    RR3 = RRR3(:); % convert into columns vector
    Htemp = CC\RR3;
    
    for i = 1:Nt
        H(:,i) = Htemp((i-1)*(Nx+1)+1:i*(Nx+1)); % Solution for f at each time level
    end
    
    for i = 1:Nt
        error_fr(r,i) = norm(F(:,i)-fprev(:,i),inf);
        error_gr(r,i) = norm(G(:,i)-gprev(:,i),inf);
        error_hr(r,i) = norm(H(:,i)-hprev(:,i),inf);
        
    end
    sumfprime = 0;
    sumg = 0;
    sumf = 0;
    sumh = 0;
    
    for j = 1:Nt+1
        sumf = sumf + d1(i,j)*F(:,j);
        sumg = sumg + d1(i,j)*G(:,j);
        sumh = sumh + d1(i,j)*H(:,j);
        sumfprime = sumfprime + d1(i,j)*D1*F(:,j);
    end
    
    for i = 1:Nt
        resf(r,i) = norm(fr3 + (1/4)*(n + 3)*fr.*fr2 - (1/2)*(n + 1)*(fr1).^2 + xi*fr2 + (1 - w)*gr + w*hr - (1/4)*(1 - n)*xi*fr1.*sumfprime + (1/4)*(1 - n)*xi*fr2.*sumf,inf);
        resg(r,i) = norm((gr2/Pr) + (1/4)*(n + 3)*fr.*gr1 + xi*gr1 - (1/4)*(1 - n)*xi*fr1.*sumg + (1/4)*(1 - n)*xi*gr1.*sumf,inf);
        resh(r,i) = norm((hr2/Sc) + (1/4)*(n + 3)*fr.*hr1 + xi*hr1 - (1/4)*(1 - n)*xi*fr1.*sumh + (1/4)*(1 - n)*xi*hr1.*sumf,inf);
    end
    
    
    
end


%-------------------------------------------------------------------------
%           Convergence and error analysis
%-------------------------------------------------------------------------
for r  = 1:iterations
    for i = 1:Nt
        normF(r,i) = max(norm(F(:,i) - fprev(:,i),inf));
        normG(r,i) = max(norm(G(:,i) - gprev(:,i),inf));
        normH(r,i) = max(norm(H(:,i) - hprev(:,i),inf));
        
        xi = Xi(i);
        fr = F(:,i);
        fr1 = D1*F(:,i);
        fr2 = D2*F(:,i);
        fr3 = D3*F(:,i);
        
        gr = G(:,i);
        gr1 = D1*G(:,i);
        gr2 = D2*G(:,i);
        
        hr = H(:,i);
        hr1 = D1*H(:,i);
        hr2 = D2*H(:,i);
        
        
        
        RESF(r,i) = norm(resf(4:end-3),inf);
        RESTS(r,i) = norm(resh(2:end-1),inf);
        RESG(r,i) = norm(resg(2:end-1),inf);
        
    end
    
    % Get norm for all time
    NORMF(r) =  norm(normF(r,:),inf);
    NORMG(r) =  norm(normG(r,:),inf);
    NORMS(r) =  norm(normH(r,:),inf);
    
end

%XI = [0 1 4 9 25];
XI = 0:1/4:tf;
SK = zeros(length(XI),1);
for ii = 1:length(XI);
    xi = XI(ii);
    ZZ = 2*xi/Lt - 1;
    FF1 = 0;
    FF2 = 0;
    FF3 = 0;
    FF  = 0;
    GG1 = 0;
    HH1 = 0;
    for j = 1:Nt+1
        LL = 1;
        for k = 1:Nt+1
            if (k == j), continue; end
            LL = LL .* (ZZ - Z(k))./(Z(j) - Z(k));
        end
        FF = FF  +  LL*F(:,j);
        FF1 = FF1 + LL*D1*F(:,j);
        HH1 = HH1 + LL*D1*H(:,j);
        GG1 = GG1 + LL*D1*G(:,j);
        FF2 = FF2 + LL*D2*F(:,j);
        FF3 = FF3 + LL*D3*F(:,j);
    end
    skin = FF2(Nx+1);
    Nussel_H = HH1(Nx+1);
    Nussel_G = GG1(Nx+1);
    SK(ii) = skin;
    fprintf('%10.2f\t %10.7f\t %10.7f\t %10.7f\n',xi,skin, -Nussel_G,-Nussel_H);
end
toc

% figure(1)
% plot(eta,FF1,'k*','LineWidth',3);
% xlabel('$\eta$','FontSize',20,'InterPreter','Latex')
% ylabel('$f''(\eta)$','FontSize',20,'InterPreter','Latex')
% set(gca,'fontsize',15)
% legend('MD-BSLLM','BSLLM')
% hold on

% figure(4)
% semilogy(1:iterations, RESF(:,1),'LineWidth',3);
% xlabel('iterations','FontSize',20,'InterPreter','Latex')
% ylabel('$E_f$','FontSize',20,'InterPreter','Latex')
% set(gca,'fontsize',20)
% %text(5,10^(-9),'$\xi = 20, 40$','FontSize',20,'InterPreter','Latex')
% hold on

figure(5)
semilogy(1:iterations,normF(:,4),'LineWidth',3);
xlabel('iterations','FontSize',20,'InterPreter','Latex')
ylabel('$E_g$','FontSize',20,'InterPreter','Latex')
set(gca,'fontsize',20)
text(5,10^(-9),'$\xi = 20, 40$','FontSize',20,'InterPreter','Latex')


index = [5 10 15];
VAR = Xi(index); % Get the values of xi at the index collocation points
for i  = 1:3
    
    switch i
        case 1
            str = 'ko-';
        case 2
            str = 'ms--';
        case 3
            str = 'b*:';
    end
    
    figure(1)
    semilogy(1:iterations,normF(:,index(i)),str);
    hold on
    xlabel('iterations','FontSize',18,'InterPreter','Latex')
    ylabel('$Res(f)$','FontSize',18,'InterPreter','Latex')
    %hleg = legend(num2str(VAR(1)),num2str(VAR(2)),num2str(VAR(3)),num2str(VAR(4)));
    hleg = legend(num2str(VAR(1)),num2str(VAR(2)),num2str(VAR(3)));
    htitle = get(hleg,'Title');
    set(htitle,'String','$\xi$','Fontsize',14,'InterPreter','Latex')
    set(gca,'fontsize',14)
end