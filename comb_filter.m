clc;
clear all;
close all;

%% design params
a = 0; % default = 0.6
r_third = 0.9; % by trial and error
r_fifth = 0.39;
wp = (0.25-0.1/2)*pi;
ws = (0.25+0.1/2)*pi;

fs = 1024;
w= linspace(-pi,pi,fs);



%% comb Filter

z_comb = [0];
p_comb = [a];
n = 1;

for x = 0:0.25*pi:2*pi

z1 = exp(1j*(ws+x));
z2 = exp(-1j*(ws+x));
z3 = exp(1j * ((pi+ws)/2)+x);
z4 = exp(-1j * ((pi+ws)/2)+x);

p1 = r_third*exp(1j*(wp+x));
p2 = r_third*exp(-1j*(wp+x));
p3 = r_fifth*exp(1j*(wp/2+x));
p4 = r_fifth*exp(-1j*(wp/2+x));



z_comb = [z_comb , z1, z2, z3, z4];
p_comb = [p_comb , p1, p2, p3, p4];


end
n = length(p_comb);

 filter_plots(z_comb, p_comb, n, w, wp, ws, fs);


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
    %xlim([-pi pi]);
    title('Magnitude response')
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
    title(sprintf('Magnitude response from -wp to wp\nmaximum passband ripple = %.4f', maximum))
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
    %xlim([-wp1 wp1]);
    title('phase response from -wp to wp')
    xlabel('Normalized Frequency (rad/sec)')
    ylabel('Phase Response')

end

%% group delay
function group_delay(num, den, w)
    gd = grpdelay(num,den,w);
    plot(w,gd);
    %xlim([-pi pi]);
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

%% filter function
function filter_plots(z, p, n, w, wp, ws, fs)
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
    pole_zero(num_filter, den_filter);
    
    figure;
    magnitude_response(h_order, w, wp, ws);
    
    figure;
    phase_response(h_order, w, wp);
    
    figure;
    group_delay(num_filter, den_filter, w);
    
    figure;
    impulse_response(num_filter, den_filter);
end