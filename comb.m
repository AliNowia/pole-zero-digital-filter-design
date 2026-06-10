clc;
clear;
close all;

%% design params
a = 0; % default = 0.6
r_third = 0.9; % by trial and error
r_fifth = 0.39;

fs = 1024;
w= linspace(-pi,pi,fs);

%% First Order filter

z_First = 0;
p_First = 0.6;
type = "first order LPF";
% filter_plots(z_First, p_First, 1, w, pi, pi, fs, type);

%% Third Order Filter
ws = (0.25+0.1/2)*pi;
wp = (0.25-0.1/2)*pi;

z1 = exp(1j*ws);
z2 = exp(-1j*ws);
p1 = a;
p2 = r_third*exp(1j*wp);
p3 = r_third*exp(-1j*wp);

z_third = [0, z1, z2];
p_third = [a, p2, p3];
n = length(p_third);

type = "third order LPF";
% filter_plots(z_third, p_third, n, w, wp, ws, fs, type);

%% Fifth Order Filter

z3 = exp(1j * (pi+ws)/2);
z4 = exp(-1j * (pi+ws)/2);
p4 = r_fifth*exp(1j * wp/2);
p5 = r_fifth*exp(-1j * wp/2);

z_fifth = [0, z1, z2, z3, z4];
p_fifth = [a, p2, p3, p4, p5];
n = length(p_fifth);

type = "fifth order LPF";
%filter_plots(z_fifth, p_fifth, n, w, wp, ws, fs, type);

%% rotation
%% HP filter (muliply by j^2 = -1)
z_HPF = [0, -z1, -z2, -z3, -z4];
p_HPF = [-a, -p2, -p3, -p4, -p5];
n = length(p_HPF);

type = "fifth order HPF";
% filter_plots(z_HPF, p_HPF, n, w, pi, pi, fs, type);

%% BP filter (multiply by -j)
z_BPF = [0, -1j*z1, -1j*z2, -1j*z3, -1j*z4];
p_BPF = [-1j*a, -1j*p2, -1j*p3, -1j*p4, -1j*p5];
n = length(p_BPF);

type = "fifth order BPF";
 filter_plots(z_BPF, p_BPF, n, w, pi, pi, fs, type);

%% the zero, pole plot
function pole_zero(num, den, type)
    zplane (num,den);
    title(sprintf("Pole-Zero Plot of %s", type), "FontWeight", "bold")
end

%% magnitude response
function magnitude_response(h, w, wp1, ws1, type)
    wp1 = abs(wp1);
    ws1 = abs(ws1);
    subplot(1, 2, 1)
    plot(w,20*log10(abs(h) / max(abs(h))));
    xlim([-pi pi]);
    title('from -\pi to \pi', "FontWeight", "bold")
    xlabel('Normalized Frequency (rad/sec)')
    ylabel('Normalized Magnitude Response')
    xline([0]);
    xline([wp1, ws1, -wp1, -ws1], 'Color', 'r')
    legend('|H(jw)|', 'transition band')
    
    % to get max magnitude ripple
    i = 512;  %index of w = 0
    max_h = max(abs(h));
    maximum = -inf;
    minimum = inf;
    while w(1, i) < wp1 
        if maximum < 20*log10(abs(h(1, i)) / max_h)
            maximum = 20*log10(abs(h(1, i)) / max_h);
        end
        if minimum > 20*log10(abs(h(1, i)) / max_h)
            minimum =  20*log10(abs(h(1, i)) / max_h);
        end
        i = i + 1;
    end
    at_wp = 20*log10(abs(h(1, i)) / max_h);
    if maximum - at_wp > at_wp - minimum
        maximum = maximum - at_wp;
    else
        maximum = at_wp - minimum;
    end

    subplot(1,2,2);
    plot(w,20*log10(abs(h) / max(abs(h))));
    xlim([-wp1 wp1]);
    xline([0]);
    xline([wp1, -wp1], 'Color', 'r')
    title(sprintf('from -wp to wp\nmaximum passband ripple = %.4f', maximum))
    xlabel('Normalized Frequency (rad/sec)')
    ylabel('Normalized Magnitude Response')

    sgtitle(sprintf('Magnitude response of %s', type), "FontWeight", "bold")
end

%% phase response
function phase_response(h, w, wp1, type)
    subplot(1, 2, 1)
    phase = angle(h);
    plot(w/pi,phase);
    title('from -\pi to \pi', "FontWeight", "bold")
    xlabel('Normalized Frequency (rad/sec)')
    ylabel('Phase Response')
    
    
    subplot(1,2,2);
    plot(w/pi,phase);
    xlim([-wp1 wp1]);
    title('from -wp to wp')
    xlabel('Normalized Frequency (rad/sec)')
    ylabel('Phase Response')
    
    sgtitle(sprintf('Phase response of %s', type), "FontWeight", "bold")
end

%% group delay
function group_delay(num, den, w, type)
    gd = grpdelay(num,den,w);
    plot(w,gd);
    xlim([-pi pi]);
    title(sprintf('group delay of %s', type))
    xlabel('Normalized Frequency (rad/sec)')
    ylabel('Group Delay (sec)')
end

%% impulse response
function impulse_response(num, den, type)
    [H_z,t] = impz(num,den);
    stem(t,H_z);
    title(sprintf('impulse response of %s', type))
    xlabel('Time (sec)')
    ylabel('Impulse Response')
end

%% filter function
function filter_plots(z, p, n, w, wp, ws, fs, type)
    num = zeros(n, 2);
    den = zeros(n, 2);
    for i = 1:n
        num(i, :) = [1, -z(1, i)];
        den(i, :) = [1, -p(1, i)];
    end
    h_filter = tf(1, 1, 1/fs);
    for i = 1:n
        h = tf([1, num(i, 2)], [1, den(i, 2)], 1/fs);  
        h_filter = h_filter * h;
    end

    [num_filter, den_filter] = tfdata(h_filter, 'v');

    h_order = freqz(num_filter, den_filter, w);

    figure;
    pole_zero(num_filter, den_filter, type);
    
    figure;
    magnitude_response(h_order, w, wp, ws, type);
    
    figure;
    phase_response(h_order, w, wp, type);
    
    figure;
    group_delay(num_filter, den_filter, w, type);
    
    figure;
    impulse_response(num_filter, den_filter, type);
end
%% Comb
z3 = exp(1j * (pi+ws)/2);
z4 = exp(-1j * (pi+ws)/2);
p4 = r_fifth*exp(1j * wp/2);
p5 = r_fifth*exp(-1j * wp/2);

z_fifth = [0, z1, z2, z3, z4];
p_fifth = [a, p2, p3, p4, p5];
n = length(p_fifth);

type = "fifth order LPF";

num = zeros(n, 2);
    den = zeros(n, 2);
    for i = 1:n
        num(i, :) = [1, -z_fifth(1, i)];
        den(i, :) = [1, -p_fifth(1, i)];
    end
    h_filter = tf(1, 1, 1/fs);
    for i = 1:n
        h = tf([1, num(i, 2)], [1, den(i, 2)], 1/fs);  
        h_filter = h_filter * h;
    end

    [num_filter, den_filter] = tfdata(h_filter, 'v');

    [H_z,t] = impz(num_filter,den_filter);

    %Comb filter using fifth order
    
    L = 8;                          % upsampling factor
    N = length(H_z);
   h_c = zeros(1, length(w)); % allocate output

for m = 1:N
    n = (m-1)*L + 1;           
    h_c(n) = H_z(m);           % place h[m] at position mL
end  
    
    na=6;
    nb=7;

    N = length(h_c);
H_c = fft(h_c, N);              % Frequency response via FFT  

[b, a] = invfreqz(H_c, w, nb, na);

zplane(b, a);