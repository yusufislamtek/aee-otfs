%%
% tic
% clc;
% clear all;

function [ERROR, AvgBER, BLERROR, BLER, FERROR, AvgFER] = AEE_OTFS(SNR_dB, SSSet)
% SSSet = 'SS_1_bps_N4_Q512_5dB_EL64_P51224_L0.0064.mat';
SignalStruct = load(SSSet);
SignalSet = double(SignalStruct.SignalSet);

%%
% SNR_dB = 10;
nFrame = 10;

%%
N = 32;
M = 32;

[Q,n] = size(SignalSet);

G = floor(M*N/n);

p = log2(Q);
SE = p/n;

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

    numBits = G*n*SE;
    bits = rand(1,numBits)>0.5;
    bits = reshape(bits,p,length(bits)/p).';
    decoded_bits = zeros(size(bits));
    dec = bi2de(bits,'left-msb')+1;
    X = SignalSet(dec,:);
%     X = reshape(X.',M,N);
    X = reshape(X.',1,numel(X)); %******
    X = [X,zeros(1,M*N-length(X))]; %****** zero padding for n=3 or n=5
    X = reshape(X.',M,N); %******

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

    Gmmse = Heff'*(Heff*(Heff')+(N0/2)*eye(size(Heff)))^(-1);
    Heff = [];
    Gmmse = sparse(Gmmse);
    y_tilda = Gmmse*y_vec;
    Gmmse = [];
    y_tilda = y_tilda(1:G*n);
    y_tilda = reshape(y_tilda.',n,G).';

    aa = 1;
    for jj = 1:1:G
        cur_y_tilda = y_tilda(jj,:);
        metric = zeros(2^p,1);
        for kk = 1:2^p
            metric(kk)=norm(cur_y_tilda-SignalSet(kk,:)).^2;
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
    totalBLE(ii) = sum(sum(error,2)>0);
    BLER(ii) = totalBLE(ii)/G;
end
FERROR=frameError;
AvgFER=FERROR/nFrame;
ERROR = sum(totalError);
AvgBER = mean(BER);
BLERROR=sum(totalBLE);
BLER = mean(BLER);
% toc





