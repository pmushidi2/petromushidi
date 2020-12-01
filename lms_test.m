%----------------------------------------------------------
% MATLAB function that implements LMS Algorithm Implementation 
%  
% Original from the paper "Simulation and Performance Analysis of Adaptive
% Filter in Noise Cancellation" by RAJ KUMAR THENUA and Surendra K. Agrawal
% Modified by Petro M. Tshakwanda in Fall 2020.
%----------------------------------------------------------
clear
clc
%load 
y = 2*sin(2*pi*0.005*(0:10000-1)');  % .005cycles/sample(one complete cyle in 200 samples)
  nvar  = 0.8;  % noise variance
  xn = randn(10000,1)*nvar;
nfilt= fir1(20,0.5);
fnoise= filter(nfilt,1,xn);
dn= y+fnoise;  % desired signal(noisy signal)

filter_size=20;
step_size= 0.0001;
iterations=10000; 
% initialise adaptive filter impulse and input vector to zero vector of length specified at command line
filter_current = zeros(filter_size,1);
input_vector = zeros(filter_size, 1);

% Loop for number of iterations.
for i=1:iterations
   
    input_vector(1)=xn(i);         % insert new sample at beginning of input vector.
    filter_output(i)=dot(filter_current, input_vector);          %Caluclate adaptive filter output
    error= dn(i)-filter_output(i);       % Calculate estimation error
    filter_current = filter_current + 2*step_size*error*input_vector;       % Update filter taps by LMS recursion
    
    % Shfit values of input vector.
    for j=filter_size:-1:2
        input_vector(j)=input_vector(j-1);    
    end
    
    error_signal(i)=error;      % store estimation error  
end
% Find moving average of error squared.
mse= (error_signal-y(1:10000)').^2;
msefinal=mean(mse)

subplot(5,1,1),plot(y(1:10000));title(' Clean  Signal s(n)')
subplot(5,1,2),plot(xn(1:10000));title(' Input: Random Noise Signal x(n)')  
subplot(5,1,3),plot(dn(1:10000));title(' Input: Clean  Signal with Random Noise d(n)')
subplot(5,1,4),plot(error_signal(1:10000));title(' (a) Output:Filtered  Signal')
subplot(5,1,5),plot(mse(1:10000));title('(b) Output:Mean Square Error')
xlabel('No of Iterations(N)','FontSize',10,'FontWeight','bold','Color','r')




