%----------------------------------------------------------
% MATLAB function that implements RLS Algorithm Implementation
%  
% Original from the paper "Simulation and Performance Analysis of Adaptive
% Filter in Noise Cancellation" by RAJ KUMAR THENUA and Surendra K. Agrawal
% Modified by Petro M. Tshakwanda in Fall 2020.
%----------------------------------------------------------
clc
clear
filter_size=19;
iterations=8000;
load          % load .mat file in workspace

y = 2*sin(2*pi*0.005*(0:8000-1)');  % .005cycles/sample(one complete cyle in 200 samples)
%  nvar  = 0.8;  % noise variance
%  xn = randn(8000,1)*nvar;
nfilt= fir1(19,.5);
fnoise= filter(nfilt,1,xn);
dn= y+fnoise;  % desired signal(noisy signal)
lambda=1;

% Initialise varaibles
filter_prev = zeros(filter_size,1);
input_vector = zeros(filter_size, 1);
W_inv_prev = eye(filter_size);
Un= zeros(filter_size, 1);
K = zeros(filter_size, 1);

% Perform RLS algorithm for set number of iterations.
for i=1:iterations
 
    input_vector(1)=xn(i);
    Un = W_inv_prev*input_vector;
    K = (1/(lambda+dot(input_vector,Un)))*Un;
    filter_output(i)=dot(filter_prev, input_vector);
    error= dn(i)-filter_output(i);
    filter_prev = filter_prev + K*error;
    W_inv_prev = (1/lambda)*(W_inv_prev - K*((input_vector')*W_inv_prev));
    % Shift along input vector
    for j=filter_size:-1:2
        input_vector(j)=input_vector(j-1);    
    end
    error_signal(i)=error;
    cost(i)=error*error;
    
end
j= error_signal.^2;
mse3= (error_signal-y(1:8000)').^2;
msefinal=mean(mse3);

P = (1/8000)*sum(y(1:8000).^2);
Pn = (1/8000)*sum(xn(1:8000).^2);
Pd = (1/8000)*sum(dn(1:8000).^2);
Pe = (1/8000)*sum(error_signal(1:8000).^2);
% SNR_pre= 10*log10(Pd/Pn)
% SNR_post= 10*log10(Pe/Pn)
error_power= Pe-P;
Noise_removal= 100-((Pe-P)/2)*100;

%find moving avarage of db attenuation (averaged to smooth output).
for i=1:iterations-2000
     db(i)=-20*log10((mean(abs(error_signal(i+2000))))/(mean(abs(y(i+2000)))));  
end
%find total avarage db attenuation
db_avg=mean(db);
subplot(5,1,1),plot(y(1:8000));title('Clean  Signal');
subplot(5,1,3),plot(dn(1:8000));title('Input: Clean  Signal with Random Noise');
subplot(5,1,2),plot(xn(1:8000));title('Input: Random Noise Signal'); 
subplot(5,1,5),plot(mse3(1:8000));title('Output:Mean Square Error');
xlabel('No of Iterations(N)','FontSize',10,'FontWeight','bold','Color','r');
subplot(5,1,4),plot(error_signal(1:8000));title('Output:Filtered  Signal');
xlabel('No of Iterations(N)','FontSize',10,'FontWeight','bold','Color','r')

% SNR calculation
snr_pre= 10*log10(mean(y.^2)./mean(fnoise.^2))

snr_post= 10*log10(mean(y.^2)./mean((error_signal-y').^2))
SNR_imprv= snr_post-snr_pre

