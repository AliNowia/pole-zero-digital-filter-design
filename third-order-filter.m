clc;
clear;
close all;

%% design params
a = 0.4; % default = 0.6 | best = 0.4
r_third = 0.86; % 0.885 for 3rd | 0.86 for fifth
r_fifth = 0.2;

fs = 1024;
w= linspace(-pi,pi,fs);

%% First Order filter

z_First = 0;
p_First = a;
type = "first order HPF";
%filter_plots(z_First, p_First, 1, w, pi, pi, fs, type);

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
%filter_plots(z_third, p_third, n, w, wp, ws, fs, type);

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

%% Comb filter testing

L = 8;

% Initialize as empty arrays so we don't keep the original roots
z_comp = []; 
p_comp = [];

for i = 0 : L-1
    % 1. Corrected Magnitude: Added parentheses around (1/L)
    new_mag_z = abs(z_third) .^ (1/L);
    new_mag_p = abs(p_third) .^ (1/L);
    
    % 2. Corrected Phase: Divided the original angle by L
    new_phase_z = (angle(z_third) + 2*pi*i) / L;
    new_phase_p = (angle(p_third) + 2*pi*i) / L;
    
    % 3. Append the new roots for this iteration
    z_comp = [z_comp, new_mag_z .* exp(1j * new_phase_z)];
    p_comp = [p_comp, new_mag_p .* exp(1j * new_phase_p)];
end
n = length(p_comp);

type = "Comb Filter";
filter_plots(z_comp, p_comp, n, w, wp, ws, fs, type);

%% rotation
%% HP filter (muliply by j^2 = -1)
z_HPF = [0, -z1, -z2, -z3, -z4];
p_HPF = [-a, -p2, -p3, -p4, -p5];
n = length(p_HPF);

type = "fifth order HPF";
%filter_plots(z_HPF, p_HPF, n, w, wp+pi/2, pi, fs, type);

%% BP filter (multiply by -j)
z_BPF = [1j*z_fifth, -1j*z_fifth];
p_BPF = [1j*p_fifth, -1j*p_fifth];
n = length(p_BPF);
type = "tenth order BPF";
filter_plots(z_BPF, p_BPF, n, w, wp+pi/2, pi, fs, type);

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
    ylabel('Normalized Magnitude Response (dB)')
    %xline([0]);
    %xline([wp1, ws1, -wp1, -ws1], 'Color', 'r')
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
    xlim([-2.51 -0.63]);
    %xline([0]);
    %xline([wp1, -wp1], 'Color', 'r')
    title(sprintf('from -wp to wp'))
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
    xlim([-2.51 -0.63]/pi);
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


function [a1, r11, r22] = find_params(w, fs, wp, ws)
    a = linspace(0.3, 0.9, 30);
    r1 = linspace(0.3, 0.9, 30);
    r2 = linspace(0.3, 0.9, 30);
    
    a1 = 0; r11 = 0; r22 = 0;
    
    z1 = exp(1j*ws);
    z2 = exp(-1j*ws);
    z3 = exp(1j * (pi+ws)/2);
    z4 = exp(-1j * (pi+ws)/2);
    
    z_fifth = [0, z1, z2, z3, z4]; 
    
    num_coeffs = poly(z_fifth); 
    
    total_iterations = length(a) * length(r1) * length(r2);

    avos = zeros(1, total_iterations);
    avms = zeros(1, total_iterations);
    c = 1;
        
    fprintf("Beginning iteration...\n");
    
    for i = 1:length(a)
        p1 = a(i);
        for j = 1:length(r1)
            p2 = r1(j)*exp(1j*wp);
            p3 = r1(j)*exp(-1j*wp);
            for k = 1:length(r2)
                p4 = r2(k)*exp(1j * wp/2);
                p5 = r2(k)*exp(-1j * wp/2);
                
                p_fifth = [p1, p2, p3, p4, p5];
                
                den_coeffs = poly(p_fifth); 
                
                h_o = freqz(num_coeffs, den_coeffs, w(415:620));
                
                h_o_db = 20*log10(abs(h_o) / max(abs(h_o)));
                
                avo = max(h_o_db); % This will always be 0
                avm = min(h_o_db); % This will be the deepest negative dip
                
                avos(c) = avo;
                avms(c) = avm;
                
                fprintf("Iteration: %d | Avo = %.2f | Avm = %.2f\n", c, avo, avm);
                c = c + 1;
                 
                if abs(avo - avm) <= 2
                    a1 = a(i);
                    r11 = r1(j);
                    r22 = r2(k);
                    fprintf("Match found! Your values: [%.2f, %.2f, %.2f]\n", a1, r11, r22);
                    return;
                end
            end
        end
    end
    
    fprintf("Sorry, couldn't find satisfying values for the ripple constraint.\n");
    
    figure;
    plot(1:total_iterations, abs(avos - avms));
    title('Ripple vs Iteration');
    xlabel('Iteration');
    ylabel('Ripple (dB)');
end


[a_1, r__1, r__2] = find_params(w, fs, wp, ws);