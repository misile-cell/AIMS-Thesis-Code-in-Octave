# Created by Octave 4.0.3, Tue May 08 10:38:06 2018 WAT <misele_kunene@aims_pc40>
function [fNt, tNt] = Initial_Guess_2_System(N,Lx,n,Pr)
[D,x] = cheb(N); %D is a differentiation matrix in \eta

% Scaled differentiation matrices

D1 = (2/Lx)*D;
D2 = D1*D1;
D3 = D1*D1*D1;

eta = Lx*(x + 1)/2; % \eta in [0,Lx]

fr = (1/2)*(1 + exp(-2*eta)) - exp(-eta);
tr = exp(-eta);

fr1 = D1*fr;
fr2 = D2*fr;
fr3 = D3*fr;

tr1 = D1*tr;
tr2 = D2*tr;

iterations = 60;

for i = 1:iterations
  a3r = ones(N + 1,1);
  a2r = (1/4)*(n + 7)*fr2;
  a1r = -(n + 1)*fr1;
  a0r = (1/4)*(n + 7)*fr2;
  a4r = ones(N+1,1);
  a5r = zeros(N+1,1);
  
  b2r = (1/Pr)*ones(N+1,1); 
  b1r = (1/4)*(n + 7)*fr;
  b0r = -n*fr1;
  
  FF = fr3 + (1/4)*fr.*fr2 - (1/2)*(n + 1)*fr1.^2 + tr;
  TT = (1/Pr)*tr2 + (1/4)*(n + 7)*fr.*tr1 - n*fr1.*tr;
  
  R1 = a3r.*fr3 + a2r.*fr2 + a1r.*fr1 + a0r.*fr + a4r.*tr - FF;
  R2 = b2r.*tr2 + b1r.*tr1 + b0r.*tr - TT;
  
  A11 = diag(a3r)*D3 + diag(a2r)*D2 + diag(a1r)*D1 + diag(a0r);
  A12 = diag(a4r);
  A13 = diag(a5r);
  
  A21 = diag(b0r);
  A22 = diag(b1r)*D1 + diag(b2r)*D2;
  A22 = zeros(N+1,N+1);
  
  % f'(\infty) = 0
  A11(1,:) = D1(1,:);
  A12(1,:) = 0;
  R1(1) = 0;
  
  % f'(0) = 0
  A11(N,:) = D1(N+1,:);
  A12(N,:) = 0;
  R1(N) = 0;
  
  %f(0) = 0
  A11(N+1,:) = 0;
  A11(N+1,N+1) = 1;
  A12(N+1,:) = 0;
  R1(N+1) = 0;
  
  %t(\infty) = 0
  A21(1,:) = 0;
  A22(1,:) = 0;
  A22(1,1) = 1;
  R2(1) = 0;
  
  %t(0) = 1
  A21(N+1,:) = 0;
  A22(N+1,:) = 0;
  A22(N+1,N+1) = 1;
  R2(N+1) = 1;
  
  AA = [A11 A12;A21 A22];
  RR = [R1;R2];
  
  Y = pinv(AA)*RR;
  
  fr = Y(1:N+1);
  tr = Y(N+2:2*(N+1));
  
  fr1 = D1*fr;
  fr2 = D2*fr;
  fr3 = D3*fr;
  
  tr1 = D1*tr;
  tr2 = D2*tr;
  
  resf = fr3 + (1/4)*(n + 7)*fr.*fr2 - (1/2)*(n + 1)*fr1.^2 + tr;
  rest = (1/Pr)*tr2 + (1/4)*(fr).*tr1 - n*fr1.*tr; 
  %fprintf('%10.0f\t %10.5e\t %10.5e\t,i,norm(resf(2:N),inf),norm(rest(2:N),inf))
  
end 
fNt = fr;
tNt = tr;
end

