clear all

nw = 32; % number of wheels
velocity = 230*5/18; %velocity of train (m/s)
L_bridge = 25; %length of bridge (m)
L_c =  15; %distance b/w 2 wheelsets (see fig 2) (m)
b = 2.5; %distance b/w 2 wheels on the same bogie 

%%%%%%%%%%%%%%%%%%%%%% TRYYYYYYYYYY! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of observations and dimension of X and Y
[dat]=xlsread('approx_z.xls'); % input measurements from the idkwut.m file
%Y = randn(10,1);
Y = dat(:,1);
T=size(Y,1);  % T is the number of observations in Y
M=1;  % M is the number of series Y

% Number of factors & lags:
p = 2;

% Generate lagged z matrix
ylag = mlag2(Y,p); 
ylag = ylag(p+1:T,:);

K = p*(M^2); % K is the number of elements in the state vector
Z = zeros((T-p)*M,K); % This is the Z matrix given in eqn 11.9 of the paper
for i = 1:T-p
    % No constant
    ztemp = []; %eye(M);
    for j = 1:p        
        xtemp = ylag(i,(j-1)*M+1:j*M);
        xtemp = kron(eye(M),xtemp);
        ztemp = [ztemp xtemp];  
    end
    Z((i-1)*M+1:i*M,:) = ztemp;
end
y = Y(p+1:T,:)';
% Time series observations
t=size(y,2);
%delta_t = 0.001;
delta_t = 0.02;
ntot = 5000;
nburn = 0.4*ntot;

% Finding the matrix A 
[A,xx,~] = findingA(L_bridge, L_c, b, velocity, nw, delta_t );

lambda_A = A(1:t-1, :); % eq(35(a)) of the paper
pp = pinv(lambda_A'*lambda_A)*lambda_A';

%-------- Now set prior means and variances (_prmean / _prvar)
% This is the Kalman filter initial condition for the time-varying parameters B(t)
% B_0 ~ N(B_0_prmean, B_0_prvar))
B_0_prmean = [1.946, -0.9932]';
B_0_prvar = 0.001^2*eye(2);
inv_B_0_prvar = inv(B_0_prvar);

%initial value of Z_{\beta} in eqn 11.10 of paper
Btdrawc = zeros(K,t+1);
Btdraw = zeros(K,t+1);
% covariance of the error term of the eqn z = {Z}{\beta} + {\alpha}A +
% \epsilon; \epsilon ~ N(0, \sigma); \sigma size = 1x1
%sigma_e \sigma_{epsilon}^2 in eqn 11.9 
%Initial parameters required to sample from \sigma_{\epsilon} ~ IW(\sigma_{\epsilon 1}, \nu_{\epsilon 1}, M)
Sigmadraw = 0.001*eye(M); 
Sigmainvdraw = inv(Sigmadraw);
Sigma_prmean = eye(M); % \sigma_{ \epsilon 0 } in eqn 11.11 of the paper
Sigma_prvar = 1; % \nu_{\epsilon 0} in eqn 11.11 of the paper
Ht = zeros(M*t,M);  
for i = 1:t
    Ht((i-1)*M+1:i*M,:) = Sigmadraw;
end

% Q is the covariance of B(t), cov matrix of the state eqn B(i+1) = B(i) +
% v(i); v ~ N(0, Q)
% Q ~ IW(k2_Q*size(subsample)*Var(B_OLS),size(subsample))
Q_prmean = 0.001*eye(p); % [\Sigma{ v0 }] 2x2
Q_prvar = 0.01; %\nu_{\v 0} 1x1
Qdraw = 0.001*eye(K);
Qchol = chol(Q_prmean);
% Storage matrices for posteriors and stuff
Bt_postmean = zeros(K,t+1);
Qmean = zeros(K,K);
Sigmamean = zeros(M,M);

z_alpha = zeros(t-1,0);
alph_A = 0.01*ones(nw,1);
alph_postmean = zeros(32,1);
alph_mid_post = 0;
Sigma_wrt_time = zeros(ntot-nburn,1);
det_Q_wrt_time = zeros(ntot-nburn,1);
tr_Q_wrt_time = zeros(ntot-nburn,1);
el1_Q_wrt_time = zeros(ntot-nburn,1);
alph_A_wrt_time = zeros(nw,ntot-nburn);
for irep = 1:ntot
    %Drawing B_1 
    vbar = zeros(K,K);
    xhy=zeros(K,1);
    i=1;
    zhat1 = Z((i-1)*M+1:i*M,:);
    yhat1 = y(:,i) - zhat1*Btdrawc(:,i) - A(i,:)*alph_A;
    HHat1 = Sigmainvdraw; 
    vbar = vbar + zhat1'*HHat1*zhat1;
    xhy = xhy + zhat1'*HHat1*yhat1;
    for i = 2:t
        zhat1 = Z((i-1)*M+1:i*M,:);
        % include the A_m \alpha in the yhat1 equation
        yhat1 = y(:,i) - zhat1*Btdrawc(:,i) - A(i-1,:)*alph_A;
        HHat1 = Sigmainvdraw;    
        vbar = vbar + zhat1'*HHat1*zhat1;
        xhy = xhy + zhat1'*HHat1*yhat1;
    end
    vbar = inv(vbar + inv_B_0_prvar);
    B0hat = vbar*(inv_B_0_prvar*B_0_prmean + xhy);
    B0draw = mvnrnd(B0hat, vbar)';
    %size(b0draw) = 1x2
    %B0draw = B0hat + chol(vbar)'*randn(K,1); % Draw from the initial condition B(0)
    % Now we're done drawing \bm{ \beta_1 } (here, this is mentioned as B_0
    % for convenience
    
    % Next step is to draw \bm{ \beta_2 },..., \bm{ \beta_M } ie,
    % B_1,...,B_t using the Durbin Koopman Smoother 
    
    ya = zeros(M,t);
    for i = 1:t
        ya(:,i) = y(:,i) - Z((i-1)*M+1:i*M,:)*B0draw;
    end
   
    %Here, the DK smoother doesn't contain the terms A_t and \aplha, so
    %there's no need to modify this code for future purposes
    [Btdrawc,llikkeep] = dk(ya,M,K,t,Qchol,Ht,Qdraw,Z);
    %%%%% Next step is calculating \aplha, step (3) in Secn 2.8 of the
    %%%%% paper
    for i = 2:t
        zhat2 = Z((i-1)*M+1:i*M,:);
        yhat2 = y(:,i) - zhat2*Btdrawc(:,i); % eqn(35(b)) of the paper
        %here, yhat1 is \bm{ Z_{\alpha} }
        z_alpha(i-1,1) = yhat2;
    end
    %%%%%%%%% calculating the value of \bm{ \alpha } at n_g th step according to
    %%%%%% Step (4) using eq (34) of the paper
    alph_A = pp*z_alpha;
   % Add on the initial condition B(0)
    for i = 1:t+1
        Btdraw(:,i) = Btdrawc(:,i) + B0draw;
    end
    %%%%%%%%% Sampling \sigma_{ \epsilon }, step (5) of Secn 2.8 of the
    %%%%%%%%% paper
    
        % Get SSE of the VAR model EQN 11.11 of the paper
    yhat = zeros(M,t);
    for i = 1:t
        yhat(:,i) = y(:,i) - Z((i-1)*M+1:i*M,:)*Btdraw(:,i); %vector of squared errors
    end
    sse_2S = zeros(M,M);
    for i = 1:t
        sse_2S = sse_2S + yhat(:,i)*yhat(:,i)'; %sum of squared errors
    end
    % Sampling \sigma_{\epsilon} ~ IW( \sigma_{\epsilon 1}, \nu_{\epsilon 1}, M )
    Sigmainv = inv(sse_2S + Sigma_prmean);
    Sigmainvdraw = wishrnd(Sigmainv,t+Sigma_prvar);
    Sigmadraw = inv(Sigmainvdraw);  
    Sigmachol = chol(Sigmadraw);
   %%%%%%%%% Sampling \Sigma_{v}, step (6) of Secn 2.8 of the
    %%%%%%%%% paper
    
     % Now get the SSE in the state equation (to estimate the covariance Q)
    Btemp = Btdraw(:,2:t+1)' - Btdraw(:,1:t)'; % EQN 11.13 of the paper
    sse_2 = zeros(K,K);
    for i = 1:t-1
        sse_2 = sse_2 + Btemp(i,:)'*Btemp(i,:);
    end
    % Draw Q, the coavariance matrix of the time-varying parameters B(t)
    %%%%%%%%% Q => \bm{\Sigma_v}
    Qinv = inv(sse_2 + Q_prmean);
    Qinvdraw = wishrnd(Qinv,t+Q_prvar);
    Qdraw = inv(Qinvdraw);      % Qdraw is a draw from the posterior of Q
    Qchol = chol(Qdraw);
    if irep > nburn
    % Save only the means of B(t), Q and SIGMA
    Bt_postmean = Bt_postmean + Btdraw;
    Qmean = Qmean + Qdraw;
    Sigma_wrt_time(irep) = Sigmadraw;
    det_Q_wrt_time(irep) = det(Qdraw);
    tr_Q_wrt_time(irep) = trace(Qdraw);
    el1_Q_wrt_time(irep) = Qdraw(1,1);
    Sigmamean = Sigmamean + Sigmadraw;
    alph_postmean = alph_postmean + alph_A;
    alph_mid_post = alph_mid_post + alph_A(1);
    alph_A_wrt_time(:,irep) = alph_A;
    end
end

Bt_postmean = Bt_postmean./ntot; % Posterior mean of B(t) (VAR regression coeff.)
Qmean = Qmean./ntot;             % Posterior mean of Q (covariance of B(t))
Sigmamean = Sigmamean./ntot;     % Posterior mean of SIGMA (VAR covariance matrix)
alph_postmean = alph_postmean./ntot;
f1 = zeros(t+1,1);
f2 = zeros(t+1,1);
damp1 = zeros(t+1,1);
damp2 = zeros(t+1,1);
for i = 1:t+1
    qq = [ 1 , Bt_postmean(1,i), Bt_postmean(2,i) ];
    s = roots(qq);
    f1(i) = ( 1/(2*pi*delta_t) )*( sqrt( real(log(s(1))).^2 + imag(log(s(1))).^2 ) ) ;
    f2(i) = ( 1/(2*pi*delta_t) )*( sqrt( real(log(s(2))).^2 + imag(log(s(2))).^2 ) ) ;
    damp1(i) = -real(log(s(1)).^2)/( sqrt( real(log(s(1))).^2 + imag(log(s(1))).^2 ) ) ;
    damp2(i) = -real(log(s(2)).^2)/( sqrt( real(log(s(2))).^2 + imag(log(s(2))).^2 ) ) ;
end

subplot(1,2,1) 
plot(delta_t:delta_t:delta_t*(T-1),Bt_postmean(1,:),'LineWidth',2)
hold on
plot(0:delta_t:delta_t*(T-1),dat(:,2),'LineWidth',2)
hold off
title('\beta_1(t)')
legend('\beta_1 -algo', '\beta_1 -actual' )
xlabel('time (s)')

subplot(1,2,2)
plot(delta_t:delta_t:delta_t*(T-1),Bt_postmean(2,:),'LineWidth',2)
hold on
plot(0:delta_t:delta_t*(T-1),dat(:,3),'LineWidth',2)
hold off
title('\beta_2(t)')
legend('\beta_2 -algo', '\beta_2 -actual' )
xlabel('time (s)')

z_calc = zeros(T,1);
for i = 2:T-2
    z_calc(i+1) = Bt_postmean(1,i+1)*z_calc(i) + Bt_postmean(2,i+1)*z_calc(i-1) + A(i+1,:)*alph_A_wrt_time(:,nburn+i+1);
end

figure
plot(0:delta_t:delta_t*(T-1),z_calc,'LineWidth',2) 
hold on
plot(0:delta_t:delta_t*(T-1),dat(:,1),'LineWidth',2)
hold off
title('Displacement  of bridge at mid span')
legend('z -algo', 'z -actual' )
ylabel('z (m)')
xlabel('time (s)')

function [Xlag] = mlag2(X,p)
[Traw,N]=size(X);
Xlag=zeros(Traw,N*p);
    for ii=1:p
        Xlag(p+1:Traw,(N*(ii-1)+1):N*ii)=X(p+1-ii:Traw-ii,1:N);
    end
end


function [atilda,llik0] = dk(y,p,m,t,Qchol,Ht,Qt,Z)

pm = p + m;
wplus = zeros(pm*t,1);
for i = 1:t
    Hchol = chol(Ht((i-1)*p+1:i*p,:));
    wplus((i-1)*pm+1:(i-1)*pm+p,1) = Hchol'*randn(p,1);
    wplus((i-1)*pm+p+1:(i-1)*pm+pm,1) = Qchol'*randn(m,1);
end

[yplus,aplus] = recur(Z,wplus,m,p,t);
[what, ahat,llik0] = kalfilt(y,Z,Ht,Qt,m,p,t);
[whatp, ahatp,llik1] = kalfilt(yplus,Z,Ht,Qt,m,p,t);

atilda = ahat - ahatp + aplus;

end

function [ydraw,alpha] = recur(Z,wdraw,m,p,t)
alpha = zeros(m,t+1);
pm = p + m;
ydraw = zeros(p,t);
for i = 1:t
    ztemp = Z((i-1)*p+1:i*p,:);
    ydraw(:,i) = ztemp*alpha(:,i) + wdraw((i-1)*pm+1:(i-1)*pm+p,1);
    alpha(:,i+1) = alpha(:,i) + wdraw((i-1)*pm+p+1:i*pm,1) ;
end
end

function [what,alph,llik] = kalfilt(y1,Z,Ht,Qt,m,p,t)

%Kalman filter 
Kkeep = zeros(m*t,p);
Lkeep = zeros(m*t,m);
Fkeep = zeros(p*t,p);
a = zeros(m,t+1);
v = zeros(p,t);
Pt = 10e-4*eye(m);
% Pt = zeros(m,m);
llik = 0;
for i = 1:t
    htemp = Ht((i-1)*p+1:i*p,:);
    ztemp = Z((i-1)*p+1:i*p,:);
    v(:,i) = y1(:,i) - ztemp*a(:,i);
    Ft = ztemp*Pt*ztemp' + htemp;
    Ftinv = inv(Ft);
    llik = llik + log(det(Ft)) + v(:,i)'*Ftinv*v(:,i);
    Fkeep ((i-1)*p+1:i*p,:) = Ftinv;
    Kt = Pt*ztemp'*Ftinv ;
    Kkeep((i-1)*m+1:i*m,:) = Kt;
    Ltt = eye(m) - Kt*ztemp;
    Lkeep((i-1)*m+1:i*m,:) = Ltt;
    a(:,i+1) = a(:,i) + Kt*v(:,i);
    Pt = Pt*Ltt' + Qt;
end
llik = -.5*llik;
%Backward recursion to evaluate rt 
rt = zeros(m,t+1);
pm = p+m;
what = zeros(pm*t,1);

for i = t:-1:1
    htemp = Ht((i-1)*p+1:i*p,:);
    ztemp = Z((i-1)*p+1:i*p,:);
    lterm = Lkeep((i-1)*m+1:i*m,:);
    fterm = Fkeep((i-1)*p+1:i*p,:)';
    kterm = Kkeep((i-1)*m+1:i*m,:);
    what((i-1)*pm+1:(i-1)*pm+p,1) = htemp*fterm*v(:,i) - htemp*kterm'*rt(:,i+1);
    what((i-1)*pm+p+1:i*pm,1) = Qt*rt(:,i+1);
    rt(:,i) = ztemp'*fterm*v(:,i) + lterm'*rt(:,i+1);
end

alph = zeros(m,t+1);
for i = 1:t
    alph(:,i+1) = alph(:,i) + Qt*rt(:,i+1);
end
end

function [A, M,arrr] = findingA(L_bridge, L_c, b, velocity, nw, delta_t )
time_of_travel = 8*(L_bridge + L_c + b + b)/velocity; 
% M => total number of time steps ie discretizing time_of_travel into M points
M = time_of_travel/delta_t + 1;
M = ceil(M);
%arrr = 0:delta_t:time_of_travel; %array of time for which vehicle travels
arrr = linspace(0,time_of_travel, M);
tauf = [0, b, b + L_c, b + L_c + b]';%initial position of each wheel before the movement starts
tau2 = tauf + repmat(L_c+2*b,4,1);
tau3 = tau2 + repmat(L_c+2*b,4,1);
tau4 = tau3 + repmat(L_c+2*b,4,1);
tau5 = tau4 + repmat(L_c+2*b,4,1);
tau6 = tau5 + repmat(L_c+2*b,4,1);
tau7 = tau6 + repmat(L_c+2*b,4,1);
tau8 = tau7 + repmat(L_c+2*b,4,1);
tau = [tauf; tau2; tau3; tau4; tau5; tau6; tau7; tau8];
A = zeros(M,nw);
for i = 1:M
    for j = 1:nw
        if (i-1)*delta_t - (tau(j)/velocity)>=0
            h1 = 1;
        else
            h1 = 0;
        end
        if (i-1)*delta_t - ((tau(j) + L_bridge)/velocity)>=0
            h2 = 1;
        else
            h2 = 0;
        end
        A(i,j) = sin( pi*( velocity*(i-1)*delta_t - tau(j))/L_bridge )*( h1 - h2 );
    end
end
end