function [fNt, gNt, hNt] = Initial_Guess_3_System(N,Lx,n,Pr,Sc,w) 
[D,x] = cheb(N);   %D is chebyshev diff matrix in \eta

% Scaled differentiation matrices
D1 = (2/Lx)*D;
D2 = D1*D1;
D3 = D1^3;


eta  = Lx*(x+1)/2;     % \eta in [0,Lx]

fr = (1/2)*(1 + exp(-2*eta)) - exp(-eta);
hr = exp(-eta);
gr = exp(-eta);

fr1 = D1*fr; 
fr2 = D2*fr; 
fr3 = D3*fr;

hr1 = D1*hr; 
hr2 = D2*hr;

gr1 = D1*gr; 
gr2 = D2*gr;


iterations = 5;

for i = 1:iterations
    
    a3r = ones(N+1,1);
    a2r = (1/4)*(n + 3)*fr;
    a1r = -(n + 1)*fr1;
    a0r = (1/4)*(n + 3)*fr2;
    a4r = (1 - w)*ones(N+1,1);
    a5r = w*ones(N+1,1);
    
    b0r = (1/4)*(n + 3)*gr1;
    b1r = (1/Pr)*ones(N+1,1);
    b2r = (1/4)*(n + 3)*fr;
    b3r = zeros(N+1,1);
    
    
    c0r = (1/4)*(n + 3)*hr1;
    c1r = zeros(N+1,1);
    c2r = (1/Sc)*ones(N+1,1);
    c3r = (1/4)*(n + 3)*fr;
    
    FF = fr3 + (1/4)*(n + 3)*fr.*fr2 - (1/2)*(n + 1)*(fr1).^2 + (1 - w)*gr + w*hr;
    GG = (gr2/Pr) + (1/4)*(n + 3)*fr.*gr1;
    HH = (hr2/Sc) + (1/4)*(n + 3)*fr.*hr1;
    
    R1 = a3r.*fr3 + a2r.*fr2 + a1r.*fr1 + a0r.*fr + a4r.*gr + a5r.*hr - FF;
    R2 = b0r.*fr + b1r.*gr2 + b2r.*gr1 - GG;
    R3 = c0r.*fr + c1r.*gr + c2r.*hr2 + c3r.*hr1 - HH;
    
    
    
    A11 = diag(a3r)*D3 + diag(a2r)*D2 + diag(a1r)*D1 + diag(a0r);
    A12 = diag(a4r);
    A13 = diag(a5r);
    
    A21 = diag(b0r);
    A22 = diag(b1r)*D2 + diag(b2r)*D1;
    A23 = zeros(N+1,N+1);
    
    A31 = diag(c0r);
    A32 = diag(c1r);
    A33 = diag(c2r)*D2 + diag(c3r)*D1;
    
    % f'(\inf) = 0
    A11(1,:) = D1(1,:);  
    A12(1,:) = 0; 
    A13(1,:) = 0; 
    R1(1) = 0;
      
    % f'(0) = 0
    A11(N,:) = D1(N+1,:);  
    A12(N,:) = 0; 
    A13(N,:) = 0; 
    R1(N) = 0;
    
    %f(0) = 0
    A11(N+1,:) = 0;  
    A11(N+1,N+1) = 1;  
    A12(N+1,:) = 0;
    A13(N+1,:) = 0;
    R1(N+1) = 0;
    
    %g(infinity) = 0
    A21(1,:) = 0;  
    A22(1,:) = 0;  
    A22(1,1) = 1; 
    A23(1,:) = 0; 
    R2(1) = 0;
    
    % g(0) = 1
    A21(N+1,:) = 0;  
    A22(N+1,:) = 0;  
    A22(N+1,N+1) = 1;
    A23(N+1,:) = 0; 
    R2(N+1) = 1;
    
    
    % h(infinity) = 0
    A31(1,:) = 0;  
    A32(1,:) = 0;  
    A33(1,:) = 0; 
    A33(1,1) = 1; 
    R3(1) = 0;
    
    % h(0) = 1
    A31(N+1,:) = 0;  
    A32(N+1,:) = 0;  
    A33(N+1,:) = 0;
    A33(N+1,N+1) = 1; 
    R3(N+1) = 1;
    
   
    
    AA = [A11 A12 A13;A21 A22 A23;A31 A32 A33]; 
    RR = [R1;R2;R3];
    
    Y = AA\RR;
    
    fr = Y(1:N+1);
    gr = Y(N + 2:2*(N+1));
    hr = Y(2*(N+1)+1:3*(N+1));
    
    fr1 = D1*fr; 
    fr2 = D2*fr; 
    fr3 = D3*fr;
    
    hr1 = D1*hr; 
    hr2 = D2*hr;
    
    gr1 = D1*gr; 
    gr2 = D2*gr;
    
    resf = fr3 + (1/4)*(n + 3)*fr.*fr2 - (1/2)*(n + 1)*(fr1).^2 + (1 - w)*gr + w*hr;
    resg = (1/Pr)*gr2 + (1/4)*(n + 3)*fr.*gr1;
    resh = (hr2/Sc) + (1/4)*(n + 3)*fr.*hr1;
    fprintf('%10.0f\t %10.5e\t %10.5e\t %10.5e\n',i,norm(resf(2:N),inf),norm(resg(2:N),inf),norm(resh(2:N),inf))
    
end
fNt = fr;
gNt = gr;
hNt = hr;
end
