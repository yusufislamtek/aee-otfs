
function [hi,li,ki,taps]=GenerateChannel(N,M,fc,delta_f,T,max_speed)
one_delay_tap = 1/(M*delta_f);
one_doppler_tap = 1/(N*T);

 
tau = [0 30 150 310 370 710 1090 1730 2510]*10^(-9);%EVA
% tau = [0 2.08 4.164 6.246 8.823]*10^(-6);%
taps = length(tau);% number of delay taps

li = round(tau/one_delay_tap);%assuming no fraction for the delay
% li=[0 one_delay_tap];
% li=[0 1 2 3];
% li=[0 sort(randi([0 3],1,taps-1))];
PDP = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9];% EVA PDP


powNorm = 10.^(PDP/10);
% pow_prof = 1/taps; %normalization of power delay profile
powNorm = powNorm/sum(powNorm);%normalization of power delay profile
hi = sqrt(powNorm).*(sqrt(1/2) * (randn(1,taps)+1i*randn(1,taps)));%channel coef. for each path
maxSpeed = max_speed*(1000/3600);
Doppler_vel = (maxSpeed*fc)/physconst('LightSpeed');
max_Doppler_tap = Doppler_vel/one_doppler_tap;
ki = (max_Doppler_tap*cos(2*pi*rand(1,taps)));%Doppler taps using Jake's spectrum
ki = round(ki);
% nu=[0 one_doppler_tap];
% nu=[0 1];
end