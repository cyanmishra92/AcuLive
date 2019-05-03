close all; clear all; clc
pkg load signal

function tt = tri_thresh(data)
  
  % Two endpoints on the curve "data"
  [maxval, maxidx] = max(data);
  x = [maxidx size(data,1)];
  y = [data(maxidx) data(end)];
##  plot(x, y, 'g*');
##  hold on;
  % The slope of the line connecting the two endpoints
  m = ( y(2) - y(1) )/( x(2) - x(1) );
  pm= - 1 / m;
##  plot(x, y);
##  hold on;
  % Point on the curve (xc,yc), point on the line (xl,yl)
  perpDist = zeros(size(data,1)-maxidx,1);
  for i = maxidx:size(data,1)
      xc = i ; yc = data(i);
      yl = ( (m * xc) + (m^2 * yc) - (m * x(1)) + y(1) )/(1+ m^2);
      xl = xc - m*(yl - yc);
      % distance^2
      d2 = (xl - xc)^2 + (yl - yc)^2;
      perpDist(i - maxidx + 1) = d2;
  end
  [val, idx] = max(perpDist);
  
  tt= idx+maxidx;

end

count_rec = 0;
sum_exec_times = 0;
nummm = 3;
tic;

for R=1:nummm
  
  
  fprintf('\nSIGNAL INDEX: %d\n', R);
  
  filename_rec = strcat('../Dataset/recordedVoice/rec_', num2str(R), '.wav');
  filename_act = strcat('../Dataset/sounds/file', num2str(R), '.wav');

  #filename_rec = 'gulsum.wav'
  #filename_act = 'gulsumact.wav';
  
  [y_rec, Fs_rec] = audioread(filename_rec);

  [y_act, Fs_act] = audioread(filename_act);
  
  Y_rec = fft(y_rec);
  Y_rec_abs = abs(Y_rec);
  Y_rec_abs = Y_rec_abs/max(Y_rec_abs);
 
  Y_act = fft(y_act);
  Y_act_abs = abs(Y_act);
  Y_act_abs = Y_act_abs/max(Y_act_abs);
  
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  L_act = length(Y_act_abs);
  T_act = 1/Fs_act;
  t = (0:L_act-1)*T_act;
  P2_act = abs(Y_act_abs);
  P1_act = P2_act(1:L_act/2+1);
  P1_act(2:end-1) = 2*P1_act(2:end-1);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  L_rec = length(Y_rec_abs);
  T_rec = 1/Fs_rec;
  t = (0:L_rec-1)*T_rec;
  P2_rec = abs(Y_rec_abs);
  P1_rec = P2_rec(1:L_rec/2+1);
  P1_rec(2:end-1) = 2*P1_rec(2:end-1);
  
  k = ones(1, 500) / 500;
  mva = conv(P1_rec, k, 'same');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  f_act = Fs_act*(0:(L_act/2))/L_act;
##  figure;
##  plot(f_act,P1_act) ;
##  title('Single-Sided Amplitude Spectrum of Actual Signal')
##  xlabel('f (Hz)')
##  ylabel('|P1(f)|')
##  hold on;
  k2 = ones(1, 500) / 500;
  mva2 = conv(P1_act, k2, 'same');
##  hold off;
##
##  
##  figure;
  ret2 = tri_thresh(mva2);

##  plot(ret2, mva2(ret2), 'r*');
##  hold on;
##  plot(mva2);
##  upperlim = length(mva)+1000;
##  xlim([0 upperlim])
##  title('Moving Average Filtered Actual Signal With Turning Point');
##  ylabel('Magnitude');
##  xlabel('Frequency');
##  hold off;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  f_rec = Fs_rec*(0:(L_rec/2))/L_rec;
##  figure;
##  plot(f_rec,P1_rec);
##  title('Single-Sided Amplitude Spectrum of Recorded Signal')
##  xlabel('f (Hz)')
##  ylabel('|P1(f)|')
##  hold on;
  k = ones(1, 500) / 500;
  mva = conv(P1_rec, k, 'same');
  
##  hold off;
##  figure;
##  #ret = tri_thresh(mva);
##
##  plot(ret2, mva(ret2), 'r*');
##  hold on;
##  plot(mva);
##  xlim([0 upperlim])
##  title('Moving Average Filtered Recorded Signal With Turning Point');
##  ylabel('Magnitude');
##  xlabel('Frequency');
##  hold off;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %MERGED GRAPH
##  figure;
##  
##  plot(mva, 'b');
##  hold on;
##  
##  plot(mva2, 'g');
##  hold on;
##
##  plot(ret2, mva2(ret2), 'r*');
##  u_lim = length(mva);
##  xlim([0 u_lim])
##  title('Moving Average Filtered Signal With Turning Point');
##  ylabel('Magnitude');
##  xlabel('Frequency');
##
##  hold on;
##  plot(ret2, mva(ret2), 'r*');
##
##  hold on;
##  
##  xran = [ret2 ret2];
##  yran = [0 mva(ret2)];
##  plot(xran, yran, 'm');
## 
##  
##  hold off;
##  legend('Recorded Signal', 'Actual Signal');
##  
  ###########################################


  y_act_sum = sum(Y_act_abs(ret2+1:20000));
  y_act_sum2 = sum(Y_act_abs(1:ret2));
  y_act_sum_diff = abs(y_act_sum / sum(Y_act_abs(1 : 20000)) * 100);
    

  y_rec_sum = sum(Y_rec_abs(ret2+1:20000));
  y_rec_sum2 = sum(Y_rec_abs(1:ret2));
  y_rec_sum_diff = abs(y_rec_sum / sum(Y_rec_abs(1 : 20000)) * 100);
    
    
  y = abs( y_act_sum_diff - y_rec_sum_diff)
  
  if(y > 0.18*y_act_sum_diff)
    fprintf('Sound %d is recorded!\n', R);
    count_rec = count_rec + 1;
  else
    fprintf('Sound %d is not recorded!\n', R);
  end
  
  if R == 1
    A = [y_act_sum_diff y_rec_sum_diff y 1.2*y_act_sum_diff];
  else
    A = [A; y_act_sum_diff y_rec_sum_diff y 1.2*y_act_sum_diff];
  end
  

end

  acc = count_rec/nummm*100;
  fprintf('The accuracy is: %f percent.\n', acc);

  sum_exec_times = toc;
  sum_exec_times = sum_exec_times/nummm;

  fprintf('The average execution time for an input is: %f seconds.\n', sum_exec_times);

  
  %figure;
  
  %bar(A);
  %title('Sums and the Threshold Value');
  
