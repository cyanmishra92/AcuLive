%% Libraries and definitions
close all; clear all; clc
pkg load signal
pkg load communications

total_diff = 0;
pop_wtd_avg = 0;
thres_freq = 1000;
cutoff_freq = 15000; #Decide the cutoff frequrncy of processing

num_files = 65; # Number of files for processing
rec_avg_ar = zeros(1, num_files);
act_avg_ar = zeros(1, num_files);
rec_std_ar = zeros(1, num_files);
act_std_ar = zeros(1, num_files);
################################################################################

%% Reading Files and FFT
for R=1:num_files
  
  filename_rec = strcat('../Dataset/recordedVoice/rec_', num2str(R), '.wav');
  filename_act = strcat('../Dataset/sounds/file', num2str(R), '.wav');
  
  [y_rec, Fs_rec] = audioread(filename_rec); # actual file

  [y_act, Fs_act] = audioread(filename_act); # Recorded file
################################################
## FFT of the recorded voice  
  Y_rec = fft(y_rec);
  Y_rec = abs(Y_rec);
  Y_rec = Y_rec/(max(Y_rec));
  
## Processing the frequency Range
  L_rec = length(Y_rec);
  f_rec = Fs_rec*(0:L_rec-1)/L_rec;
  f_rec = f_rec(1:cutoff_freq);
  P_rec = Y_rec(1:cutoff_freq);
  
################################################
## FFT of the actual voice 
  Y_act = fft(y_act);
  Y_act = abs(Y_act);
  Y_act = Y_act/(max(Y_act));
    
## Processing the frequency Range
  L_act = length(Y_act);
  f_act= Fs_act*(0:L_act-1)/L_act;
  f_act = f_act(1:cutoff_freq);
  P_act = Y_act(1:cutoff_freq);

################################################
## Adding the |f| from 1000 to 90000
  y_rec_hf_avg = (sum(P_rec(thres_freq:cutoff_freq)))/(cutoff_freq - thres_freq +1);
  rec_avg_ar(R) = y_rec_hf_avg;
  rec_std_ar(R) = std(P_rec(thres_freq:cutoff_freq));
  rec_avg_ar =  transpose(rec_avg_ar);
  rec_std_ar =  transpose(rec_std_ar);
  
  y_act_hf_avg = (sum(P_act(thres_freq:cutoff_freq)))/(cutoff_freq - thres_freq +1);
  act_avg_ar(R) = y_act_hf_avg;
  act_std_ar(R) = std(P_act(thres_freq:cutoff_freq));
  act_avg_ar =  transpose(act_avg_ar);
  act_std_ar =  transpose(act_std_ar);
################################################
  #avg_diff = y_rec_hf_avg - y_act_hf_avg
  wtd_avg = (5*y_rec_hf_avg + 5*y_act_hf_avg)/10;
  pop_wtd_avg =  pop_wtd_avg + wtd_avg;
################################################
end
pop_wtd_avg = pop_wtd_avg/R; #For averaging over population
################################################################################

%% Plotting the actual FFT
#length(f_act)
#length(P_act)
plot(f_act, P_act);
grid on;
title('Single-Sided Amplitude Spectrum of Actual Voice')
xlabel('f (Hz)');
ylabel('|(f)|');
figure;
  
#length(f_act)
#length(P_act)
plot(f_rec, P_rec);
grid on;
title('Single-Sided Amplitude Spectrum of Recoeded Voice');
xlabel('f (Hz)');
ylabel('|(f)|');
################################################################################

%% Scoring Policy and Plotting on test data

## Processing Test Data
#test_id = 52
test_id = randint(1, 1, num_files)
#test_id = 1

filename_test = strcat('../Dataset/recordedVoice/rec_', num2str(test_id), '.wav');
filename_test = strcat('../Dataset/sounds/file', num2str(test_id), '.wav');
[y_test, Fs_test] = audioread(filename_test); # test_file
## FFT of the recorded voice  
  Y_test = fft(y_test);
  Y_test = abs(Y_test);
  Y_test = Y_test/(max(Y_test));
  
## Processing the frequency Range
  L_test = length(Y_test);
  f_test = Fs_test*(0:L_test-1)/L_test;
  f_test = f_test(1:cutoff_freq);
  P_test = Y_test(1:cutoff_freq);

  #plot(f_test, )
  
baseline_rec = zeros(1,length(P_test));
#baseline_rec = zeros(1,cutoff_freq);
baseline_rec = baseline_rec + pop_wtd_avg;
baseline_rec = baseline_rec(1:cutoff_freq);
length(P_test);
length(f_test);
length(baseline_rec);
test_score = sum(P_test > baseline_rec);
confidence = ((max(test_score))/(length(test_score)))*100;


figure;
plot(f_test, P_test);
#hold on
#plot(baseline_rec, "linewidth", 5)
grid on
title('Baseline comparision of test data')
xlabel('f (Hz)')
ylabel('|(f)|')

%{
Best-fit values	 
Slope	5.425 ± 0.8111
Y-intercept	0.1050 ± 0.07025
X-intercept	-0.01936
1/Slope	0.1843
 
0.95 Confidence Intervals	 
Slope	3.835 to 7.015
Y-intercept	-0.03266 to 0.2427
X-intercept	-0.06087 to 0.004841
 
Goodness of Fit	 
R square	0.2590
Sy.x	0.4338
 
Is slope significantly non-zero?	 
F	44.74
DFn,DFd	1,128
P Value	< 0.0001
Deviation from horizontal?	Significant
 
Data	 
Number of XY pairs	130
Equation	Y = 5.425*X + 0.1050


####### STD

Best-fit values	 
Slope	2.923 ± 1.258
Y-intercept	0.2674 ± 0.1091
X-intercept	-0.09149
1/Slope	0.3422
 
95% Confidence Intervals	 
Slope	0.4561 to 5.389
Y-intercept	0.05352 to 0.4813
X-intercept	-1.020 to -0.01028
 
Goodness of Fit	 
R square	0.04043
Sy.x	0.4936
 
Is slope significantly non-zero?	 
F	5.394
DFn,DFd	1,128
P Value	0.0218
Deviation from horizontal?	Significant
 
Data	 
Number of XY pairs	130
Equation	Y = 2.923*X + 0.2674

%}
y_test_hf_avg = (sum(P_test(thres_freq:cutoff_freq)))/(cutoff_freq - thres_freq +1);
rgr_avg_confidence = (5.425*y_test_hf_avg + 0.1050)*100
rgr_std_confidence = (2.923*y_test_hf_avg + 0.2674)*100

if rgr_avg_confidence > 50 || rgr_std_confidence > 50 || confidence >50
  display('Recorded voice detected, please authenticate')
end
close all;