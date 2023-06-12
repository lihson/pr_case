clear all; 

n=30;
use_sh_lambda = false;
senario = 0;% 0 for LOS; 1 for NLOS

% PREALLOCATING
Pr_case_lambda = zeros(3,n);
for i = 1:n
    % DENSITY %
    lambda_l = 10/pi; % Line density 10/pi per km
    lambda_v = 15; % Number of receiving vehicular nodes per km
    lambda_1 = 5; % Number of tier 1 nodes per km
    lambda_2 = i*0.5; % Number of tier 2 nodes per km
    
    % TRANSMISSION POWER %
    p1_dbm = 46; % Transmit power of tier 1 node in dBm
    p1 = 10^(.1*(p1_dbm-30));
    p2_dbm = 20; % Transmit power of tier 2 node in dBm
    p2 = 10^(.1*(p2_dbm-30));
    pv_dbm = 20;
    pv = 10^(.1*(pv_dbm-30));
    
    % GAIN %
    mlg_2 = 1; % Mainlobe gain TIER 2 node
    slg_2 = .01; % Sidelobe gain TIER 2 node
    mlg_1 = 1; % Mainlobe gain TIER 1 node
    %slg_1 = .01; % Sidelobe gain TIER 1 node, not applicable
    v_m = 0.01;
    v_s0 = 1;
    v_s1 = 0.01;
    
    % BIAS %
    b1_db = 0; % Selection bias of TIER 1 node in dB
    b1 = 10^(.1*b1_db);
    b2_db = 0; % Selection bias of TIER 2 node in dB
    b2 = 10^(.1*b2_db);

    A_MS = p1*b1*mlg_1/(p2*b2*mlg_2); %dl
    B_MS = pv*v_m/(pv*v_s0);          %ul
    
    if senario == 0
        alpha1 = 2;
        alpha2 = 2;
    else
        alpha1 = 4;
        alpha2 = 4;
    end

    % NAKAGAMI-M FADING PARAMETERS %
    m1 = 1; % Tier 1
    m20 = 2; % Tier 2 typical line
    m21 = 1; % Tier 2 other lines
    
    % SHADOWING PARAMETERS %
    ln_mu1 = 0;% Mean of log-normal shadowing gain in dB for TIER 1
    ln_sig1 = 4; % Std deviation of log-normal shadowing gain in dB for TIER 1
    ln_mu20 = 0; % Mean of log-normal shadowing gain in dB for TIER 2 TYP LINE
    ln_sig20 = 2; % Std deviation of log-normal shadowing gain in dB for TIER 2 TYP LINE
    ln_mu21 = 0;% Mean of log-normal shadowing gain in dB for TIER 2 OTHER LINES
    ln_sig21 = 4;% Std deviation of log-normal shadowing gain in dB for TIER 2 OTHER LINES
    
    sf_c = exp(2/alpha1*log(10)/10*ln_mu1 + (2/alpha1*log(10)/10)^2*ln_sig1^2/2);   % tier 1
    sf_u0 = exp(1/alpha2*log(10)/10*ln_mu20 + (1/alpha2*log(10)/10)^2*ln_sig20^2/2);% tier 2 line 
    sf_a = exp(2/alpha2*log(10)/10*ln_mu21 + (2/alpha2*log(10)/10)^2*ln_sig21^2/2); % tier 2 other

    sh_lambda_1 = sf_c*lambda_1;  
    sh_lambda_2 = sf_u0*lambda_2; 
    
    if use_sh_lambda
        lambda_1_int = sh_lambda_1;
        lambda_2_int = sh_lambda_2;
    else
        lambda_1_int = lambda_1;
        lambda_2_int = lambda_2;
    end
    f1 = @(x) 2*lambda_2_int*exp(-lambda_1_int*pi*B_MS^(2/alpha1)*x.^(2*alpha2/alpha1) - 2*lambda_2_int*x);
    f2 = @(x) 2*lambda_2_int*exp(-lambda_1_int*pi*A_MS^(2/alpha1)*x.^(2*alpha2/alpha1) - 2*lambda_2_int*x);
    f1_int = integral(f1, 0, Inf);
    f2_int = integral(f2, 0, Inf);
    
    pr_int_case1 = 1 - f1_int;
    pr_int_case2 = f1_int - f2_int;
    pr_int_case4 = f2_int;

    Pr_case_lambda(1,i) = pr_int_case1;
    Pr_case_lambda(2,i) = pr_int_case2;
    Pr_case_lambda(3,i) = pr_int_case4;
    x(i) = i*0.5/5;
    i = i + 1;
end
 
figure
plot(x,Pr_case_lambda(1,:));
hold on;
plot(x,Pr_case_lambda(2,:));
hold on;
plot(x,Pr_case_lambda(3,:));
grid on;
axis([0 2.5 0 1]);

legend('Case 1','Case 2','Case 4');
xlabel("λs/λm")
ylabel("Joint Probability of Association")
if senario == 0
    if use_sh_lambda
        title('LOS-sh');
    else
        title('LOS');
    end
else
    ylim([0 0.7]);
    if use_sh_lambda
        title('NLOS-sh');
    else
        title('NLOS');
    end
end
