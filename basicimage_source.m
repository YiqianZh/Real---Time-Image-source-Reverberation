% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% basic image source model
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all; close all;

%% Parameters %%

Fs = 48000;                             % Initialising Sample rate
Cair = 343;                                % Speed of sound in air (m/s)

alpha = 0.1;                             % Acoustic Absorption
R = sqrt(1 - alpha);                     % Reflection coefficient
%% Room dimensions (in meters) %%

% Lengths along x,y and z directions
Lx = 5;                                
Ly = 4;
Lz = 3;

V = Lx*Ly*Lz;                            % Volume of room

A1 = Lx*Ly ; 
A2 = Lx*Lz ; 
A3 = Ly*Lz;    % Area of room walls


%% Calculation of T60 of room %%

alpha_bar = 2*alpha*(A1 + A2 + A3);            % Formula for alpha_bar
T60 = (24*log(10)*V)/(Cair*alpha_bar);     % Formula for T60

%N = (Cair*T60)/min([Lx,Ly,Lz]);            % Number of reflections 

% Number of reflections Nx,Ny,Nz 
Nx = Cair*T60/Lx;
Ny = Cair*T60/Ly;
Nz = Cair*T60/Lz;


%Initialising vector for impulse respone output
IR = zeros(ceil(T60*Fs),1);
IR_x = zeros(ceil(T60*Fs),1);
IR_y = zeros(ceil(T60*Fs),1);
IR_z = zeros(ceil(T60*Fs),1);

IR_v = zeros(ceil(T60*Fs),1);
IR_t = zeros(ceil(T60*Fs),1);
IR_r = zeros(ceil(T60*Fs),1);
IR_s = zeros(ceil(T60*Fs),1);
IR_u = zeros(ceil(T60*Fs),1);

IR_q = zeros(ceil(T60*Fs),1);
IR_o = zeros(ceil(T60*Fs),1);
IR_m = zeros(ceil(T60*Fs),1);
IR_k = zeros(ceil(T60*Fs),1);
IR_l = zeros(ceil(T60*Fs),1);
IR_n = zeros(ceil(T60*Fs),1);
IR_p = zeros(ceil(T60*Fs),1);

nv = zeros(ceil(T60*Fs),3);

%% Source and Receiver Positions %%

% Position of source in cartesian coordinates(p,q,r)
p = 2;
q = 2;
r = 2;

% Position of receiver in cartesian coordinates(a,b,c)
a = Lx/sqrt(2);
b = Ly/sqrt(2);
c = Lz/sqrt(2); 

%% Computation %%

for d = -Nx:Nx
    if rem(d,2)~=0                          % Odd case
        Ad = (d+1)*Lx-p-a;
    else                                    % Even case
        Ad = d*Lx+p-a;
    end
    
    
    for e =  -Ny:Ny
        if rem(e,2)~=0                      % Odd case
            Be = (e+1)*Ly-q-b;
        else                                % Even case
            Be = e*Lx+q-b;
        end
    
        
        for f = -Nz:Nz
            if rem(f,2)~=0                  % Odd case
                Cf = (f+1)*Lz-r-c;
            else                            % Even case
                Cf = f*Lz+r-c;
            end
            
            %% Calculation of impulse response %%
            
            % Calculation of distances
            l_def = sqrt(Ad*Ad + Be*Be + Cf*Cf); 
            
            nv = [Ad,Be,Cf]./l_def;
            
            % 0-order
            w = 1;
            
            % 1-order
            x = Ad/l_def;
            y = Be/l_def;
            z = Cf/l_def;
         
            % 2-order
            v_2 = sqrt(3)*x*y;
            t_2 = sqrt(3)*y*z;
            r_2 = (1/2)*(3* z^2 - 1);
            s_2 = sqrt(3)*x*z;
            u_2 = sqrt(3/4)*(x^2 - y^2);
            
            % 3-order
            q_3 = sqrt(5/8)*y*(3*x^2-y^2);
            o_3 = sqrt(15)*x*y*z;
            m_3 = sqrt(3/8)*y*(5*z^2-1);
            k_3 = (1/2)*z*(5*z^2-3);
            l_3 = sqrt(3/8)*x*(5*z^2-1);
            n_3 = sqrt(15/4)*z*(x^2-y^2);
            p_3 = sqrt(5/8)*x*(x^2-3*y^2);
            
            % Number of reflections
            omiga = abs(d) + abs(e) + abs(f);
            
            % Calculation of reflected impulse magnitude
            g = power(R,omiga)/l_def;
            
            % Calculation of arrival times
            t = l_def/Cair;
                                                         
            bin = round(t*Fs);           % Calculation of sample bin 
            if bin <= length(IR)
                IR(bin) = IR(bin) + g*w;   % Impulse response and its magnitude
            end               
            
            if bin <= length(IR_x)
                IR_x(bin) = IR_x(bin) + g*x;   % Impulse response and its magnitude
            end   
            if bin <= length(IR_y)
                IR_y(bin) = IR_y(bin) + g*y;   % Impulse response and its magnitude
            end   
            if bin <= length(IR_z)
                IR_z(bin) = IR_z(bin) + g*z;   % Impulse response and its magnitude
                IR_k(bin) = IR(bin) + g*k_3;
                IR_r(bin) = IR(bin) + g*r_2;
            end   
            
            if bin <= length(IR_x) && bin <= length(IR_y)
                IR_v(bin) = IR(bin) + g*v_2;
                IR_u(bin) = IR(bin) + g*u_2;
                IR_q(bin) = IR(bin) + g*q_3;
                IR_p(bin) = IR(bin) + g*p_3;
            end
            
            if bin <= length(IR_x) && bin <= length(IR_z)
                IR_s(bin) = IR(bin) + g*s_2;
                IR_l(bin) = IR(bin) + g*l_3;
            end
            
            if bin <= length(IR_y) && bin <= length(IR_z)
                IR_t(bin) = IR(bin) + g*t_2;
                IR_m(bin) = IR(bin) + g*m_3;
            end
            
            if bin <= length(IR_x) && bin <= length(IR_y) && bin <= length(IR_z)
                IR_o(bin) = IR(bin) + g*o_3;
                IR_n(bin) = IR(bin) + g*n_3;
            end
        end
    end
end


%% build in the ambisonics outputs in order to get 4 channels of audio
amb_audio = zeros(ceil(T60*Fs),16);

output_audio(:,1) = IR;   % w
output_audio(:,2) = IR_x; % x
output_audio(:,3) = IR_y; % y
output_audio(:,4) = IR_z; % z


output_audio(:,5) = IR_v; % x
output_audio(:,6) = IR_t; % y
output_audio(:,7) = IR_r; % z
output_audio(:,8) = IR_s; % x
output_audio(:,9) = IR_u; % y

output_audio(:,10) = IR_q; % z
output_audio(:,11) = IR_o; % x
output_audio(:,12) = IR_m; % y
output_audio(:,13) = IR_k; % z
output_audio(:,14) = IR_l; % x
output_audio(:,15) = IR_n; % y
output_audio(:,16) = IR_p; % z

audiowrite("IR.wav",output_audio,Fs); %write audio   

%%  Convolution with input audio %%




% 
% %Input Audio
% [input_audio,Fs] = audioread('guzheng.wav');
% % In case of stereo,convert to mono
% input_audio = 0.5*sum(input_audio,2);   
% 
% % Output Audio
% output_audio = conv(input_audio,IR);
%             
% % Play Output Audio
% soundsc(output_audio, Fs);
% 
% % Plot impulse response
% plot([0:length(IR)-1]/Fs, IR);
% title('Plot of Impluse Response vs. time');
% xlabel('Time(s)');ylabel('Magnitude');








