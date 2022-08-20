clear;
clc;
%%%CREATÝNG SÝGNAL%%%
bits=[1 0 0 1 0 1 ]; %Signal Bits
bit_per=1;         
%Converting to digital signal
signal=[]; 
for n=1:1:length(bits)
    if bits(n)==1
       se=ones(1,100);
    else bits(n)=0;
        se=zeros(1,100);
    end
     signal=[signal se];
end
t1=bit_per/100:bit_per/100:100*length(bits)*(bit_per/100);
subplot(5,1,1);
plot(t1,signal,'lineWidth',1);grid on;
axis([ 0 bit_per*length(bits) -.5 1.5]);
xlabel(' time(sec)');
title('Digital Signal');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modulating Signal
A=1;            % Amplitude of carrier signal 
b=2*bits-1; 
T=3;
Eb=T/2;         
fc=3/T;         % Carrier Signal Frequency

t=linspace(0,length(bits),length(bits)*200);
N=length(t);    %Bits size

Nsb=N/length(bits); % Sample Size for each bit
disp(Nsb);
dd=repmat(bits',1,Nsb);
bb=repmat(b',1,Nsb);
dw=dd'; 
dw=dw(:)';
bw=bb';
bw=bw(:)';
w=sqrt(2*Eb/T)*A*cos(2*pi*fc*t+(pi/2));

bpsk_w=bw.*w; % Modulated Signal
subplot(5,1,2)
plot(t,bpsk_w,'lineWidth',1);grid on;hold on;
axis([ 0 length(bits) -A A]);
yticks(-A:A/2:A);
xlabel(' time(sec)');
title('Modulated Signal');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Adding AWGN
reqSNR=2;
signal_energy = norm(bpsk_w(:))^2;                    % energy of the signal
noise_energy = signal_energy/(10^(reqSNR/10));           % energy of noise 
noise_var = noise_energy/(length(bpsk_w(:))-1);     % variance of noise 
noiseStd = sqrt(noise_var);                      % std. deviation of noise 
noise = noiseStd*randn(size(bpsk_w));           % noise
noisySig = bpsk_w+noise;                        % noisy signal

subplot(5,1,3)
plot(t,noise,'r','lineWidth',1);grid on;
xlabel(' time(sec)');
title('AWGN');
subplot(5,1,4)
plot(t,noisySig,t,bpsk_w,'lineWidth',1);grid on;
xlabel(' time(sec)');
title('Modulated Signal with AWGN');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Demodulation

demod=[];
for n=Nsb:Nsb:length(noisySig)
    t=linspace(0,1,200);   
  y=A*cos(2*pi*fc*t+(pi/2));                           
  conv=y.*noisySig((n-(Nsb-1)):n);
  integ=trapz(t,conv)    
  logic=round((2*integ))                                     
  if(logic>0)                                      
                      
    a=1;
  else
    a=0;
  end
  demod=[demod a];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
bit=[];
for n=1:length(demod);
    if demod(n)==1;
       se=ones(1,200);
    else demod(n)==0;
        se=zeros(1,200);
    end
     bit=[bit se];
end
subplot(5,1,5)
t=linspace(0,length(bits),length(bits)*200);
plot(t,bit,'LineWidth',1);grid on;
axis([ 0 n -.5 1.5]);
xlabel(' time(sec)');
title('Recived Information after Demodulation');
