%----------------------------------------------------------
% MATLAB function that implements NLMS Algorithm Implementation
%  
% Original from the paper "Simulation and Performance Analysis of Adaptive
% Filter in Noise Cancellation" by RAJ KUMAR THENUA and Surendra K. Agrawal
% Modified by Petro M. Tshakwanda in Fall 2020.
%----------------------------------------------------------
clear
clc
load 
y = 2*sin(2*pi*0.005*(0:10000-1)');  % .005cycles/sample(one complete cyle in 200 samples)
nvar  = 0.8;  % noise variance
xn = randn(10000,1)*nvar;
nfilt= fir1(20,0.5);
fnoise= filter(nfilt,1,xn);
dn= y+fnoise;  % desired signal(noisy signal)

filter_size=20;
c=0.001;
alpha=0.02;
iterations=10000;

% initialise variables
filter_current = zeros(filter_size,1);
input_vector = zeros(filter_size, 1);

% perform the NLMS algorithm for set number of iterations
for i=1:iterations
   
    % get input sample
    input_vector(1)=xn(i);
    filter_output(i)=dot(filter_current, input_vector);
    error= dn(i)-filter_output(i);
    % calculate step size
    step_size=alpha/(c+(dot(input_vector, input_vector)));
    % update filter tap weight vector
    filter_current = filter_current + step_size*error*input_vector;
    % shift input vector
    for j=filter_size:-1:2
        input_vector(j)=input_vector(j-1);    
    end
    error_signal(i)=error;
end

mse3= (error_signal-y(1:10000)').^2;
msefinal=mean(mse3);

P = (1/10000)*sum(y(1:10000).^2);
Pn = (1/10000)*sum(xn(1:10000).^2);
Pd = (1/10000)*sum(dn(1:10000).^2);
Pe = (1/10000)*sum(error_signal(1:10000).^2);
SNR_pre= 10*log10(Pd/Pn);
SNR_post= 10*log10(Pe/Pn);


error_power= Pe-P;
percentage_removal= ((Pn-error_power)/Pn)*100;

%find moving avarage of db attenuation (averaged to smooth output).
for i=1:iterations-2000
     db(i)=-20*log10((mean(abs(error_signal(i+2000))))/(mean(abs(y(i+2000)))));  
 end
%find total avarage db attenuation
db_avg=mean(db);

subplot(5,1,1),plot(y(1:10000));title('Clean  Signal')
subplot(5,1,3),plot(dn(1:10000));title('Input: Clean  Signal with Random Noise')
subplot(5,1,2),plot(xn(1:10000));title('Input: Random Noise Signal') 
subplot(5,1,5),plot(mse3(1:10000));title('Output:Mean Square Error')
xlabel('No of Iterations(N)','FontSize',10,'FontWeight','bold','Color','r');
subplot(5,1,4),plot(error_signal(1:10000));title('Output:Filtered  Signal')
xlabel('No of Iterations(N)','FontSize',10,'FontWeight','bold','Color','r');
 

snr_pre= 10*log10(mean(y.^2)./mean(fnoise.^2))

snr_post= 10*log10(mean(y.^2)./mean((error_signal-y').^2))
SNR_imprv= snr_post-snr_pre

