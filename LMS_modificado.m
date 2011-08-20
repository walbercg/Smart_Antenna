%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LMS Algorithm example 8.6 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----- Givens -----%
clear all
close all
clc

d = .5;  % element spacing in terms of wavelength d = lambda/2
N = 2;   % number of elements in array
thetaS = input('   What is the desired users AOA (in degrees)?   ');
%thetaI = input('   What is the interferers AOA(in degrees)?   ');

%----- Desired Signal & Interferer -----%

T = 1E-3;
t = (1:100)*T/100;
it = 1:100;
S = cos(2*pi*t/T);
thetaS = thetaS*pi/180;                  % desired user AOA
%I = randn(1,100);  
%thetaI = thetaI*pi/180;                    % interferer AOA
  
%----- Create Array Factors for each user's signal for linear array -----%

vS = []; 
%vI = [];
i = 1:N;
vS = exp(1j*(i-1)*2*pi*d*sin(thetaS)).';
%vI = exp(1j*(i-1)*2*pi*d*sin(thetaI)).';

%----- Solve for Weights using LMS -----%

%snr = 10;       % signal to noise ratio
w = zeros(N,1); %initial weights
X = (vS); %+vI);
Rx = X*X';        % Auto correlation matrix
a = trace(Rx);
mu = 1/(4*real(a)); % step size

wi = zeros(N,max(it));

for n = 1:100
    
    x = S(n)*vS; % + I(n)*vI;
    y = w'*x; %y = w*x.';
    e = conj(S(n)) - y;  %error signal    
    esave(n) = abs(e)^2;
    w = w + mu*conj(e)*x; %LMS solution: w = w +mu*e*conj(x);
    wi(:,n) = w;
    yy(n) = y;
end

w = (w./w(1));    % normalize results to first weight

%----- Plot Results -----%

theta = -pi/2:.01:pi/2;
AF = zeros(1,length(theta));

%--- Determine the array factor for linear array ---%

for i = 1:N
    AF = AF + w(i)'.*exp(1j*(i-1)*2*pi*d*sin(theta));
end

figure
plot(theta*180/pi,abs(AF)/max(abs(AF)),'k')
xlabel('AOA (deg)')
ylabel('|AF_n|')
axis([-90 90 0 1.1])
set(gca,'xtick',[-90 -60 -30 0 30 60 90])
grid on

figure;
plot(it,S,'k',it,yy,'k--')
xlabel('No. of Iterations')
ylabel('Signals')
legend('Desired signal','Array output')

disp('%------------------------------------------------------------------------%')
disp(' ')
disp(['   The weights for the N = ',num2str(N),' ULA are:'])
disp(' ')

for m = 1:length(w)
    disp(['   w',num2str(m),' = ',num2str(w(m))])
end

disp(' ')

figure;
plot(it,abs(wi(1,:)),'kx',it,abs(wi(2,:)),'ko','markersize',2)
%,it,abs(wi(3,:)),'ks',it,abs(wi(4,:)),'k+',it,abs(wi(5,:)),'kd'
xlabel('Iteration no.')
ylabel('|weights|')
figure;plot(it,esave,'k')
xlabel('Iteration no.')
ylabel('Mean square error')