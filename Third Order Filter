clc;
clear all;
close all;

a = 0.6;
fs = 512;
w= linspace(-pi,pi,fs);


%% Third Order Filter
p1 = 0.6;
p2 = 0.2*pi+1j;
p3 = 0.2*pi-1j;
wp = 0.2*pi;
ws = 0.3*pi;

num1 = [1 -ws];
den1 = [1 -p1];

num2 = [1 0];
den2 = [1 -(0.55+0.4*j)];

num3 = [1 0];
den3 = [1 -(0.55-0.4*j)];

h1 = tf(num1, den1, fs);
h2 = tf(num2, den2, fs);
h3 = tf(num3, den3, fs);
h = h1*h2*h3;

[num den] = tfdata(h,'v');

h = freqz(num,den,w,fs);
figure;
pole_zero(num, den);

figure;
magnitude_response(h, w, wp);

figure;
phase_response(h, w, wp);
 
figure;
group_delay(num, den, fs, w);

figure;
impulse_response(num, den)





%% the zero, pole plot
function pole_zero(num, den)
    zplane (num,den);
end

%% magnitude response
function magnitude_response(h, w, wp1)
    subplot(1, 2, 1)
    maximum = max(20*log10(h));
    plot(w,20*log10(abs((h)))/maximum);
    xlim([-pi pi]);
    title('Magnitude response')
    xlabel('Normalized Frequency (rad/sec)')
    ylabel('Normalized Magnitude Response')
   
    subplot(1,2,2);
    plot(w,20*log10(abs(h))/maximum);
    xlim([-wp1 wp1]);
    title('Magnitude response from -wp to wp')
    xlabel('Normalized Frequency (rad/sec)')
    ylabel('Normalized Magnitude Response')
    

end

%% phase response
function phase_response(h, w, wp1)
    subplot(1, 2, 1)
    phase = angle(h);
    plot(w/pi,phase);
    title('phase response')
    xlabel('Normalized Frequency (rad/sec)')
    ylabel('Phase Response')

    
    subplot(1,2,2);
    plot(w/pi,phase);
    xlim([-wp1 wp1]);
    title('phase response from -wp to wp')
    xlabel('Normalized Frequency (rad/sec)')
    ylabel('Phase Response')

end

%% group delay
function group_delay(num, den, fs, w);
    gd = grpdelay(num,den,fs,w);
    plot(w,gd);
    xlim([-pi pi]);
    title('group delay')
    xlabel('Normalized Frequency (rad/sec)')
    ylabel('Group Delay (sec)')
end

%% impulse response
function impulse_response(num, den)
    [H_z,t] = impz(num,den);
    stem(t,H_z);
    title('impulse response')
    xlabel('Time (sec)')
    ylabel('Impulse Response')
end
