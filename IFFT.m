close all; clear all; clc
pkg load signal

function tt = tri_thresh(data)
  
  % Two endpoints on the curve "data"
  [maxval, maxidx] = max(data);
  x = [maxidx size(data,1)];
  y = [data(maxidx) data(end)];
  plot(x, y, 'g*');
  hold on;
  % The slope of the line connecting the two endpoints
  m = ( y(2) - y(1) )/( x(2) - x(1) );
  pm= - 1 / m;
  plot(x, y);
  hold on;
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


for R=1:3
  
  fprintf('\nSIGNAL INDEX: %d\n', R);
  
  filename_rec = strcat('../Dataset/recordedVoice/rec_', num2str(R), '.wav');
  filename_act = strcat('../Dataset/sounds/file', num2str(R), '.wav');

  #filename_rec = 'gulsum.wav'
  #filename_act = 'gulsumact.wav';
  
  [y_rec, Fs_rec] = audioread(filename_rec);

  [y_act, Fs_act] = audioread(filename_act);
  
  Y_rec = fft(y_rec);
  #Y_rec = Y_rec(1:5000);
  Y_rec_abs = abs(Y_rec);
  Y_rec_abs = Y_rec_abs/max(Y_rec_abs);
  %Y_rec = Y_rec/max(Y_rec_abs);
 
  Y_act = fft(y_act);
  #Y_act = Y_act(1:5000); #For keeping only 5000 points of the frequency line
  Y_act_abs = abs(Y_act);
  %Y_act = Y_act/max(Y_act_abs);
  Y_act_abs = Y_act_abs/max(Y_act_abs);
  
 
  %figure;
  %plot(abs(Y_act),'b')
  %hold on;
  %plot(abs(Y_rec), 'r')

  L_act = length(Y_act);
  T_act = 1/Fs_act;
  t = (0:L_act-1)*T_act;
  #P2_act = abs(Y_act/L_act);
  P2_act = abs(Y_act);
  P1_act = P2_act(1:L_act/2+1);
  P1_act(2:end-1) = 2*P1_act(2:end-1);

  f_act = Fs_act*(0:(L_act/2))/L_act;
  %figure
  %plot(f_act,P1_act) 
  %title('Single-Sided Amplitude Spectrum of Actual')
  %xlabel('f (Hz)')
  %ylabel('|P1(f)|')

  figure;
  k2 = ones(1, 8000) / 8000;
  mva2 = conv(P1_act, k2, 'same');
  ret2 = tri_thresh(mva2);
  plot(ret2, mva2(ret2), 'r*');
  hold on;
  plot(mva2);
  title('Moving Average Filtered Actual Signal With Turning Point');
  ylabel('Magnitude');
  xlabel('Frequency');
  
  y_act_sum = sum(Y_act_abs(ret2+1:20000));
  y_act_sum2 = sum(Y_act_abs(1:ret2));
  y_act_sum_diff = abs(y_act_sum / sum(Y_act_abs(1 : 20000)) * 100);
  
  y_act_avg = y_act_sum;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  
  
  
  L_rec = length(Y_rec);
  T_rec = 1/Fs_rec;
  t = (0:L_rec-1)*T_rec;
  #P2_rec = abs(Y_rec/L_rec);
  P2_rec = abs(Y_rec);
  P1_rec = P2_rec(1:L_rec/2+1);
  P1_rec(2:end-1) = 2*P1_rec(2:end-1);

  f_rec = Fs_rec*(0:(L_rec/2))/L_rec;
  %figure
  %plot(f_rec,P1_rec) 
  %title('Single-Sided Amplitude Spectrum of Recorded Signal')
  %xlabel('f (Hz)')
  %ylabel('|P1(f)|')

  figure;
  k = ones(1, 8000) / 8000;
  mva = conv(P1_rec, k, 'same');
  ret = tri_thresh(mva);
  plot(ret,mva(ret), 'r*');
  hold on;
  plot(mva);
  title('Moving Average Filtered Recorded Signal With Turning Point');
  ylabel('Magnitude');
  xlabel('Frequency');

  y_rec_sum = sum(Y_rec_abs(ret2+1:20000));
  y_rec_sum2 = sum(Y_rec_abs(1:ret2));
  y_rec_sum_diff = abs(y_rec_sum / sum(Y_rec_abs(1 : 20000)) * 100);
  
  y_rec_avg = y_rec_sum;
  
  
  ###########################################


  
  ##
  ##
  #inv_fft = ifft(Y_act);
  #audiowrite ('recovered.wav', inv_fft, Fs)

  diff = y_rec_avg - y_act_avg;
  
    
  y = abs( y_act_sum_diff - y_rec_sum_diff)
  
  if(y > 0.15*y_act_sum_diff)
    fprintf('Sound %d is recorded!\n', R);
  else
    fprintf('Sound %d is not recorded!\n', R);
  end
  
  if R == 1
    %A = [y_act_avg y_rec_avg y 1.2*y_rec_avg];
    A = [y_act_sum_diff y_rec_sum_diff y 1.2*y_act_sum_diff];
  else
    %A = [A; y_act_avg y_rec_avg y 1.2*y_rec_avg];
    A = [A; y_act_sum_diff y_rec_sum_diff y 1.2*y_act_sum_diff];
  end
  

end

  %figure;
  
  %bar(A);
  %title('Sums and the Threshold Value');
  
  
  
  
  

  
  
  
