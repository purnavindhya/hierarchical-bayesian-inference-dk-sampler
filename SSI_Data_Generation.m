clear all
L_bridge = 25;
velocity = 230*5/18;
b = 2.5;

nw = 32;
L_c = 15;
%delta_t = 0.001;
delta_t = 0.02;
[A,M,arrr] = findingA(L_bridge, L_c, b, velocity, nw, delta_t );
vecc = sum(A');


f_actual = zeros(M,1);
damp_actual = zeros(M,1);
for i=1:M
    f_actual(i) = 2.7 - 0.135*vecc(i);
    damp_actual(i) = 0.02 + 0.004*vecc(i);
end

m_b = 22*1000*L_bridge/2;
%m_b = 22*1000/L_bridge;
k_b = (2.7*2*pi)^2*m_b;
c_b = 0.02*sqrt( m_b*k_b );

delta_m_v = 135*1000/9.81;
delta_k_v = (delta_m_v/0.05^2)*(k_b/m_b);
delta_c_v = 0.2*c_b*sqrt(delta_m_v*delta_k_v/(m_b*k_b));

% delta_m_v = 0.1*m_b;
% delta_k_v = 0.1*k_b;
% delta_c_v = 0.1*c_b;
kappa_m = zeros(M,1);
kappa_k = zeros(M,1);
kappa_c = zeros(M,1);
gam = zeros(M,1);
gam(1) = m_b/delta_t^2 + c_b/( 2*delta_t );
beta_1 = zeros(M,1);
beta_2 = zeros(M,1);
beta_1(1) = -( k_b - 2*( 1 )*m_b/delta_t^2 )/gam(1);
beta_2(1) = -( ( 1 )*m_b/delta_t^2 - ( 1 )*c_b/(2*delta_t) )/gam(1);

for i = 1:M
    kappa_m(i) = vecc(i).*vecc(i).*delta_m_v/m_b;
    kappa_k(i) = vecc(i).*vecc(i).*delta_k_v/k_b;
    kappa_c(i) = vecc(i).*vecc(i).*delta_c_v/c_b;
    gam(i+1) = ( 1 + kappa_m(i) )*m_b/( delta_t^2 ) + ( 1 + kappa_c(i) )*c_b/( 2*delta_t );
    beta_1(i+1) = -( k_b - 2*( 1 + kappa_m(i) )*m_b/delta_t^2 )/gam(i+1);
    beta_2(i+1) = -( (( 1 + kappa_m(i) )*m_b/delta_t^2) - (( 1 + kappa_c(i) )*c_b/(2*delta_t)) )/gam(i+1);
end
beta_1 = beta_1(1:M);
beta_2 = beta_2(1:M);

%gam_const = ( m_b/(delta_t^2) + (c_b/(2*delta_t)) );
z = zeros(M,1);
P_S = 135*1000*ones(nw,1); %135 kN
z(2) = beta_1(1)*z(1) + 0 + A(1,:)*P_S/gam(1);
for i = 2:M-1
    z(i+1) = beta_1(i+1)*z(i) + beta_2(i+1)*z(i-1) + A(i+1,:)*P_S/gam(i+1);
end

% for i=2:M-1
%     z(i+1) = z(i)*( -2*(m_b + vecc(i)*delta_m_v )/delta_t^2 + ( k_b + vecc(i)*delta_k_v ) )/( (m_b + vecc(i))/delta_t^2 + (c_b + vecc(i)*delta_c_v)/(2*delta_t) ) ...
%         + z(i-1)*( -(m_b + vecc(i)*delta_m_v )/delta_t^2 + ( c_b + vecc(i)*delta_c_v)/(2*delta_t) )/( (m_b + vecc(i))/delta_t^2 + (c_b + vecc(i)*delta_c_v)/(2*delta_t) ) ...
%         + A(i,:)*P_S/(( (m_b + vecc(i))/delta_t^2 + (c_b + vecc(i)*delta_c_v)/(2*delta_t) ));
% end


figure
subplot(3,2,1);
plot(arrr, vecc)
title('plot of sum(A) over all wheelsets')
xlabel('time (s)')

subplot(3,2,2);
plot(arrr, f_actual)
title('plot of f_{actual}')
xlabel('time (s)')
ylabel('frequency (Hz)')

subplot(3,2,3);
plot(arrr,damp_actual)
title('plot of \eta_{actual}')
xlabel('time (s)')

subplot(3,2,4);
plot(arrr, beta_1)
title('plot of \beta_{1}')
xlabel('time (s) ')

subplot(3,2,5);
plot(arrr, beta_2)
title('plot of \beta_{2}')
xlabel('time (s) ')

subplot(3,2,6);
plot(arrr, z)
title('plot of z_{actual}')
xlabel('time (s)')
ylabel('displacement (m)')

xlswrite('approx_z.xls',[z,beta_1,beta_2,f_actual,damp_actual])

save('idkwut.m', 'z', 'beta_1', 'beta_2', 'f_actual', 'damp_actual', 'M', 'A')

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