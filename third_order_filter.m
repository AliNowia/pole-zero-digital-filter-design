clc;
close all;

% design params
a = 0.3362; % default = 0.6
r = 0.9355; % by trial and error
r2 = 0.7448;

fs = 1024;
w= linspace(-pi,pi,fs);

%% Third Order Filter

ws = (0.25+0.1/2)*pi;
wp = (0.25-0.1/2)*pi;



z1 = exp(1j*ws);
z2 = exp(-1j*ws);
p1 = a;
p2 = r*exp(1j*wp);
p3 = r*exp(-1j*wp);

num1 = [1 -z1];
den1 = [1 -p1];

num2 = [1 -z2];
den2 = [1 -p2];

num3 = [1 0];
den3 = [1 -p3];

h1 = tf(num1, den1, 1/fs);
h2 = tf(num2, den2, 1/fs);
h3 = tf(num3, den3, 1/fs);
H = h1*h2*h3;

[num, den] = tfdata(H,'v');

h = freqz(num,den,w);


figure;
pole_zero(num, den);

figure;
magnitude_response(h, w, wp, ws);

%%
figure;
phase_response(h, w, wp);
 
figure;
group_delay(num, den, w);

figure;
impulse_response(num, den)

%% Fifth order filter



z3 = exp(1j * (pi+ws)/2);
z4 = exp(-1j * (pi+ws)/2);
p4 = r2*exp(1j * wp/2);
p5 = r2*exp(-1j * wp/2);

num4 = [1 -z3];
den4 = [1 -p4];

num5 = [1 -z4];
den5 = [1 -p5];

h4 = tf(num4, den4, 1/fs);
h5 = tf(num5, den5, 1/fs);

H = H * h4 * h5;

[num, den] = tfdata(H, 'v');

h = freqz(num, den, w);

figure;
pole_zero(num, den);

figure;
magnitude_response(h, w, wp, ws);

figure;
phase_response(h, w, wp);
 
figure;
group_delay(num, den, w);

figure;
impulse_response(num, den)


%% the zero, pole plot
function pole_zero(num, den)
    zplane (num,den);
end

%% magnitude response
function magnitude_response(h, w, wp1, ws1)
    wp1 = abs(wp1);
    ws1 = abs(ws1);
    subplot(1, 2, 1)
    plot(w,20*log10(abs(h) / max(abs(h))));
    xlim([-pi pi]);
    title('Magnitude response')
    xlabel('Normalized Frequency (rad/sec)')
    ylabel('Normalized Magnitude Response')
    xline([wp1, ws1, -wp1, -ws1], 'Color', 'r')
    legend('|H(jw)|', 'transition band')
   
    subplot(1,2,2);
    plot(w,20*log10(abs(h) / max(abs(h))));
    xlim([-wp1 wp1]);
    title('Magnitude response from -wp to wp')
    xlabel('Normalized Frequency (rad/sec)')
    ylabel('Normalized Magnitude Response')
   
end

%% phase response
function phase_response(h, w, wp1)
    wp1 = abs(wp1);
    subplot(1, 2, 1)
    phase = angle(h)*180/pi;
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
function group_delay(num, den, w)
    gd = grpdelay(num,den,w);
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
