tic
clear all; clc;

SNR_dB = 10;
repeat = 1;

for ii=1:repeat
      [ERROR(ii), AvgBER(ii), BLERROR(ii), AvgBLER(ii),FERROR(ii),AvgFER(ii)] = OTFS(SNR_dB,2);
      [ERROR1(ii), AvgBER1(ii), BLERROR1(ii), AvgBLER1(ii),FERROR1(ii),AvgFER1(ii)] = AEE_OTFS(SNR_dB, 'SS_1_bps_N4_Q128_5dB_EL64_L0.0065_ME7.6173.mat');
end



Err_sum =       sum(ERROR);
BER_avg =       mean(AvgBER);
Blerr_sum =     sum(BLERROR);
BLER_avg =      mean(AvgBLER);
FER_sum =       sum(FERROR);
FER_avg =       mean(AvgFER);

Err_sum1 =      sum(ERROR1);
BER_avg1 =      mean(AvgBER1);
Blerr_sum1 =    sum(BLERROR1);
BLER_avg1 =     mean(AvgBLER1);
FER_sum1 =      sum(FERROR1);
FER_avg1 =      mean(AvgFER1);



fprintf("0: [BER:[%d %.4e] BLER:[%d %.4e] FER:[%d %.4e]] \n" + ...
        "1: [BER:[%d %.4e] BLER:[%d %.4e] FER:[%d %.4e]] \n", + ...
    Err_sum,BER_avg,Blerr_sum,BLER_avg,FER_sum,FER_avg, ...
    Err_sum1,BER_avg1,Blerr_sum1,BLER_avg1,FER_sum1,FER_avg1);

toc