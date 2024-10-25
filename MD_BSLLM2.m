clear % Clear all variables on workspace
clc   % Clear screen
close all;

Nt = 5;
Nx = 60;  % Collocation points in \eta and \xi
[D,y] = cheb(Nx); % D is chebyshev diff matrix in \eta
[d,Z] = cheb(Nt); % d is chebyshev diff matrix in \xi

a = 0;
b = 20;
Lx = b - a;


% Scaled differentiation matrices
D1 = (2/Lx)*D;
D2 = D1*D1;
D3 = D1^3;
eta = (Lx*y)/2 + (a+b)/2;  % x in [0,Lx]
I  = eye(Nx+1,Nx+1); % Identity matrix

% Time interval
t0 = 0;
tf = 20;
% Number of sub - intervals
p = 40;
Xi = linspace(t0,tf,p+1);  %  Divide the T into p equal intervals
xi_k = zeros(Nt+1,length(Xi)-1);

% Solution at each time interval
Solution_t_F = zeros(Nx+1,Nt+1,length(Xi)-1);
Solution_t_T = zeros(Nx+1,Nt+1,length(Xi)-1);
%Solution_t_H = zeros(Nx+1,Nt+1,length(Xi)-1);


% Constants
n = 1/2;
%w = 1/2;
Pr = 0.7;
%Sc = 0.6;


Tolerance_level = 1e-7;% Convergence tolerance level

%-------------------------------------------------------------------------
%  Initial guess
%-------------------------------------------------------------------------
[fNt, tNt] = Initial_Guess_2_System(Nx,Lx,n,Pr);
%-------------------------------------------------------------------------


% Assume that the initial approximation is fNt for all t
F = repmat(fNt,1,Nt+1);
T = repmat(tNt,1,Nt+1);
%H = repmat(hNt,1,Nt+1);
tic

% Multidomain loop
for k = 1:length(Xi)-1;
    
    scale_t = 2/(Xi(k+1) - Xi(k));
    
    % Scaled Differentiation Matrix in t
    d1 = scale_t*d;
    
    % Store the time intervals like this
    xi_k(:,k) = Z*(Xi(k+1) - Xi(k))/2 + (Xi(k+1) + Xi(k))/2;
    
    
    Tolerance_level = 1e-9;% Convergence tolerance level in
    
    r = 1; Tol = 1;
    iterations = 20;
    
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
            
            tr = T(:,i);
            tr1 = D1*T(:,i);
            tr2 = D2*T(:,i);
            
            xi = xi_k(i,k);
            
            sumf = 0;
            sumfprime = 0;
            sumt = 0;
            
            
            for j = 1:Nt+1
                sumf = sumf + d1(i,j)*F(:,j);
                sumfprime = sumfprime + d1(i,j)*(D1*F(:,j));
                sumt = sumt + d1(i,j)*T(:,j);
                
            end
            
            
            
            % Define the coefficients
            a0r = (1/4)*(n + 7)*fr2;
            a1r = -(n + 1)*fr1 - (1/4)*(1 - n)*xi*sumfprime;
            a2r = (1/4)*(n + 7)*fr + xi + (1/4)*(1 - n)*xi*sumfprime;
            a3r = 1;
            a4r = (1/4)*(1 - n)*xi*fr2;
            a5r = -(1/4)*(1 - n)*xi*fr1;
            
            % Define the right hand side
            Ff = fr3 + (1/4)*(n + 7)*fr.*fr2 - (1/2)*(n + 1)*(fr1).^2 + xi*fr2 + tr - (1/4)*(1 - n)*xi*fr1.*sumfprime + (1/4)*(1 - n)*xi*fr2.*sumf;
            Rf = a3r.*fr3 + a2r.*fr2 + a1r.*fr1 + a0r.*fr + a4r.*sumf + a5r.*sumfprime - Ff;
            R1 = Rf - a4r.*(d1(i,Nt+1).*fNt) - a5r.*(d1(i,Nt+1).*(D1*fNt));
            
            % Define matrices
            Ai = diag(a3r)*D3 + diag(a2r)*D2 + diag(a1r)*D1 + diag(a0r);
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
        % SOLVE FOR T
        %------------------------------------------------------------------------
        for i = 1:Nt     % for i = 0,1,2,...,Nt
            % For convenience of notation and to ensure that notation in the code
            % Is consistent with that in the write-up we define
            
            tprev(:,i) = T(:,i);
            
            fr = F(:,i);
            fr1 = D1*F(:,i);
            fr2 = D2*F(:,i);
            
            tr = T(:,i);
            tr1 = D1*T(:,i);
            tr2 = D2*T(:,i);
            
            xi = xi_k(i,k);
            
            sumt = 0;
            sumf = 0;
            
            for j = 1:Nt+1
                sumf = sumf + d1(i,j)*F(:,j);
                sumt = sumt + d1(i,j)*T(:,j);
                
            end
            
            % Define the coefficients
            b0r = -n*fr1;
            b1r = (1/4)*(n + 7)*fr + xi + (1/4)*(1 - n)*xi*sumf;
            b2r = 1/Pr;
            b3r = -(1/4)*(1 - n)*xi*fr1;
            
            
            % Define the right hand side
            Tt = (tr2/Pr) + (1/4)*(n + 7)*fr.*tr1 -n*fr1.*tr + xi*tr1 - (1/4)*(1 - n)*xi*fr1.*sumt + (1/4)*(1 - n)*xi*tr1.*sumf;
            Rt = b2r.*tr2 + b1r.*tr1 + b0r.*tr +  b3r.*sumt - Tt;
            R2 = Rt - b3r.*(d1(i,Nt+1).*tNt);
            
            
            % Define matrices
            Bi = diag(b2r)*D2 + diag(b1r)*D1 + diag(b0r);
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
        Ttemp = BB\RR2;
        
        for i = 1:Nt
            T(:,i) = Ttemp((i-1)*(Nx+1)+1:i*(Nx+1)); % Solution for f at each time level
        end
        
        %-------------------------------------------------------------------------
        %           Convergence and error analysis
        %-------------------------------------------------------------------------
        
        for i = 1:Nt
            normF(r,i) = max(norm(F(:,i) - fprev(:,i),inf));
            normT(r,i) = max(norm(T(:,i) - tprev(:,i),inf));
            
            xi = xi_k(i,k);
            fr = F(:,i);
            fr1 = D1*F(:,i);
            fr2 = D2*F(:,i);
            fr3 = D3*F(:,i);
            
            tr = T(:,i);
            tr1 = D1*T(:,i);
            tr2 = D2*T(:,i);
            
            
            sumfprime = 0;
            sumt = 0;
            sumf = 0;
            
            for j = 1:Nt+1
                sumf = sumf + d1(i,j)*F(:,j);
                sumt = sumt + d1(i,j)*T(:,j);
                sumfprime = sumfprime + d1(i,j)*D1*F(:,j);
            end
            
            resf = fr3 + (1/4)*(n + 7)*fr.*fr2 - (1/2)*(n + 1)*(fr1).^2 + xi*fr2 + tr - (1/4)*(1 - n)*xi*fr1.*sumfprime + (1/4)*(1 - n)*xi*fr2.*sumf;
            rest = (tr2/Pr) + (1/4)*(n + 7)*fr.*tr1 + xi*tr1 -n*fr1.*tr - (1/4)*(1 - n)*xi*fr1.*sumt + (1/4)*(1 - n)*xi*tr1.*sumf;
            
            RESF(r,i) = norm(resf(4:end-3),inf);
            REST(r,i) = norm(rest(2:end-1),inf);   
        end
        
        % Get norm for all time
        NORMF(r) =  norm(normF(r,:),inf);
        NORMT(r) =  norm(normT(r,:),inf);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save the solution at each interval
    Solution_t_F(:,:,k) = F;
    Solution_t_T(:,:,k) = T;
    
    
    % Use result from previous subinterval as as initial guess for the next interval i.e u_i(t_{n-1})=u_{n-1}{t_{n-1}}
    fNt = F(:,1);
    tNt = T(:,1);
    
    
    % Solution at the end point of the next interval which happens to be initial solution at t_k
    F(:,Nt+1) = fNt;
    T(:,Nt+1) = tNt;
    
    
    % Skin Friction and Nusselt number
    Skinfriction_F = D2*F(:,Nt+1);
    Nusselt_number_T = D1*T(:,Nt+1);
    fprintf('\n %10.2f\t %10.7f\t  %10.7f\t %10.7f\n',Xi(k+1),Skinfriction_F(Nx+1),-Nusselt_number_T(Nx+1));
    
    fprime(:,k) = D1*fNt;
    tprime(:,k) = D1*tNt;
    tp(:,k) = tNt;
    
    conv_error_fr(:,k) = normF(:,1);
    conv_error_tr(:,k) = normT(:,1);
    
end
toc


for k = 1:1:length(Xi)-1
  fprintf('\n iterations number %d\t %10.2f\t %10.7f\t %10.7f\t %10.7f\n',k,Xi(k+1),Skinfriction_F(Nx+1),-Nusselt_number_T(Nx+1));

end

figure(1)
plot(eta,fprime(:,5:5:40),'LineWidth',3);
xlabel('$\eta$','FontSize',20,'InterPreter','Latex')
ylabel('$f''(\eta)$','FontSize',20,'InterPreter','Latex')
set(gca,'fontsize',15)
text(0.5,0.025,'$\zeta = 5, 10, 15, 20, 25, 30, 35, 40$','FontSize',16,'InterPreter','Latex')
print("Figure_1.pdf", "-dpdflatexstandalone")
hold on

figure(2)
plot(eta,tp(:,5:5:40));
xlabel('$\eta$','FontSize',20,'InterPreter','Latex')
ylabel('$\theta(\eta)$','FontSize',20,'InterPreter','Latex')
set(gca,'fontsize',15)
text(0.5,0.6,'$\zeta = 5, 10, 15, 20, 25, 30, 35, 40$','FontSize',16,'InterPreter','Latex')
print ("Figure_2.pdf", "-dpdflatexstandalone")
hold on

figure(3)
semilogy(1:iterations,conv_error_fr(:,5:5:40));
xlabel('iterations','FontSize',20,'InterPreter','Latex')
ylabel('$E_f$','FontSize',20,'InterPreter','Latex')
set(gca,'fontsize',15)
text(3.6,10^(-7),'$\zeta = 5, 10, 15, 20, 25, 30, 35, 40$','FontSize',16,'InterPreter','Latex')
print ("Figure_3.pdf", "-dpdflatexstandalone")
hold on
 
figure(4)
semilogy(1:iterations,conv_error_tr(:,5:5:40));
xlabel('iterations','FontSize',20,'InterPreter','Latex')
ylabel('$E_{\theta}$','FontSize',20,'InterPreter','Latex')
set(gca,'fontsize',15)
text(3.6,10^(-7),'$\zeta = 5, 10, 15, 20, 25, 30, 35, 40$','FontSize',16,'InterPreter','Latex')
print ("Figure_4.pdf", "-dpdflatexstandalone")
hold on

figure(5)
semilogy(1:iterations,RESF)
xlabel('iterations','FontSize',20,'InterPreter','Latex')
ylabel('$||res f||_{\infty}$','FontSize',20,'InterPreter','Latex')
set(gca,'fontsize',15)
text(3.6,10^(-7),'$\zeta = 5, 10, 15, 20, 25, 30, 35, 40$','FontSize',16,'InterPreter','Latex')
print("Figure_5.pdf", "-dpdflatexstandalone")
hold on

figure(6)
semilogy(1:iterations,REST)
xlabel('iterations','FontSize',20,'InterPreter','Latex')
ylabel('$||res f||_{\infty}$','FontSize',20,'InterPreter','Latex')
set(gca,'fontsize',15)
text(3.6,10^(-7),'$\zeta = 5, 10, 15, 20, 25, 30, 35, 40$','FontSize',16,'InterPreter','Latex')
print("Figure_6.pdf", "-dpdflatexstandalone")
hold on




