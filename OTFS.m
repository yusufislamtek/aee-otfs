% tic
% clc;
% clear all;

function [ERROR, AvgBER, BLERROR, AvgBLER,FERROR,AvgFER] = OTFS(SNR_dB,Q)
%%
% SNR_dB = 10;
nFrame = 10;

%%
N = 32;
M = 32;
% n=3;

% Q = 4;
SignalSet = qammod(0:Q-1,Q,'Gray','UnitAveragePower',1);

SE = log2(Q);
p = SE;

%% Time and frequency resources
fc = 4*10^9; 
delta_f = (15)*10^3; 
T = 1/delta_f; 

SNR = 10.^(SNR_dB/10);
N0 = 1./SNR;
sigma = sqrt(N0/2);

Fn = dftmtx(N);  
Fn = Fn./norm(Fn);  
max_speed = 506;

EYE_MAT = sparse(eye(M*N));
PI_mat = sparse([EYE_MAT(:, 2:end) EYE_MAT(:, 1)]);
Delta_mat = sparse(diag(exp(1j*2*pi.*(0:M*N-1)/M/N)));

%% 
totalError = zeros(1,nFrame);
BER = zeros(1,nFrame);
frameError=0;
for ii = 1:nFrame

    numBits = M*N*SE;
    bits = rand(1,numBits)>0.5;
    bits = reshape(bits,p,length(bits)/p).';
    decoded_bits = zeros(size(bits));
    dec = bi2de(bits,'left-msb')+1;
    X = SignalSet(dec);
    X = reshape(X.',M,N);

    [hi,li,ki,taps]=GenerateChannel(N,M,fc,delta_f,T,max_speed);
    noise = sigma*(randn(M*N,1) + 1i*randn(M*N,1));
    H = zeros(M*N);
    for itao = 1:taps
        H = H + hi(itao)*PI_mat^li(itao)*Delta_mat^ki(itao);
    end  
    Heff=((kron(Fn,eye(M)))*sparse(H)*(kron(Fn',eye(M))));
    H = [];
    noise_DD=kron(Fn,eye(M))*noise;       
    x_vec=reshape(X,N*M,1);
    y_vec=Heff*x_vec+noise_DD;   

    Gmmse = Heff'*(Heff*(Heff')+(N0)*eye(size(Heff)))^(-1);
    Heff = [];
    Gmmse = sparse(Gmmse);
    y_tilda = Gmmse*y_vec;
    Gmmse = [];
%         y_tilda = reshape(y_tilda,M,N);
%         y_tilda = reshape(y_tilda.',1,M*N/1);

    aa = 1;
    blockError=0;
    for jj = 1:1:length(y_tilda)
        metric = zeros(2^p,1);
        for kk = 1:Q
            metric(kk)=norm(y_tilda(jj)-SignalSet(kk)).^2;
        end
        [~,ind] = min(metric);
        decoded_bits(aa,:) = de2bi(ind-1,p,'left-msb');
        aa = aa+1;
    end
    if ~isequal(bits,decoded_bits)
        frameError=frameError+1;
    end
    error = xor(bits,decoded_bits);
    totalError(ii) = sum(error(:));
    BER(ii) = totalError(ii)/numel(bits);
end
FERROR=frameError;
AvgFER=FERROR/nFrame;
ERROR = sum(totalError);
AvgBER = mean(BER);
BLERROR = 0;% sum(totalBLE);
AvgBLER = 0;% mean(BLER);
% toc

