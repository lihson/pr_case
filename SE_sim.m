clear;

r=5; % radius
n=65; % simulation of sbs density upper bound 
senario = 0; % 0 for LOS; 1 for NLOS
num_simulations = 10; % number of simulation for each λs/λm

% PREALLOCATING
pr_NLOS_decouple = zeros(3,n);
sir_NLOS_decouple = zeros(6,n);

for i = 1:4:n
    fprintf('i = %d\n', i);
    
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
    
    if senario == 0
        alpha1 = 4;
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
    
    
    % NUMBER OF TIER 1 CELL AND ROADS
    numc = poissrnd(lambda_1*pi*(r^2),1,num_simulations); % Number of TIER 1 nodes in observation window for all iterations   服从lambade分布的  
    numl = poissrnd(lambda_l*(2*r*pi),1,num_simulations); % Number of LINES in observation window for all iterations
    
    % PREALLOCATING
    pr_case = zeros(3, num_simulations);
    sir_case = zeros(6, num_simulations);

    for iS=1:num_simulations
        Nl = numl(iS); % Number of lines
        Nc = numc(iS); % Number of TIER 1 nodes

        % GENERATING THE RANDOM LINES
        rho1 = [0 -r+2*r*rand(1,Nl) ];
        [~,ind2] = sort(abs(rho1)); % Sort rho
        rho = rho1(ind2) - 2*r*(rho1(ind2)>r) + 2*r*(rho1(ind2)<-r);
        theta =pi*rand(1,length(rho));
        len_chord = 2*(r^2-(rho.^2)).^.5; % Length of chord of each line
        
        % GENERATING THE VEHICULAR
        numv = poissrnd(lambda_v*len_chord);
        numv = numv + (numv==0);
        total_v = sum(numv);
        
        clens = cumsum(numv);% Cumlative number of vehiculars of each road
        idx=zeros(1,clens(end));idx2=idx; idx3=idx;
        idx([1 clens(1:end-1)+1]) = diff([0 len_chord]); %%mei条路之间长度差值
        len_c = cumsum(idx);% Road length for each vehiculars
        idx2([1 clens(1:end-1)+1]) = diff([0 rho]);
        rho_vec = cumsum(idx2);% Road rho for each vehiculars
        idx3([1 clens(1:end-1)+1]) = diff([0 theta]);
        theta_vec = cumsum(idx3);% Road theta for each vehiculars
        
        x1 = (-len_c/2 + len_c.*rand(1,total_v));% Relative position on the road each vehiculars located
        x = zeros(size(x1));
        for j=1:length(rho) % Sort vehiculars by relative position on each road 
            x( sum(numv(1:j-1))+1: sum(numv(1:j))) = sort(x1(sum(numv(1:j-1))+1: sum(numv(1:j))));%
        end
        beta = atan(x./rho_vec);
        gamma = (theta_vec+beta);
        u = (rho_vec.^2+x.^2).^.5;
        signm = sign(rho_vec) + 1*(rho_vec==0);
        pts_v = ((u.*cos(gamma).*signm) + 1j*(u.*sin(gamma).*signm)).';
        
        % GENERATING TIER 2 NODES
        numu = poissrnd(lambda_2*len_chord);
        numu = numu + (numu==0);
        total_u = sum(numu);
        
        clens = cumsum(numu);
        idx=zeros(1,clens(end));idx2=idx; idx3=idx;
        idx([1 clens(1:end-1)+1]) = diff([0 len_chord]);
        len_c = cumsum(idx);
        idx2([1 clens(1:end-1)+1]) = diff([0 rho]);
        rho_vec = cumsum(idx2);
        idx3([1 clens(1:end-1)+1]) = diff([0 theta]);
        theta_vec = cumsum(idx3);
        
        x1 =  -len_c/2 + len_c.*rand(1,total_u);
        x = zeros(size(x1));
        for j=1:length(rho)
            x( sum(numu(1:j-1))+1: sum(numu(1:j))) = sort(x1(sum(numu(1:j-1))+1: sum(numu(1:j))));
        end
        beta = atan(x./rho_vec);
        gamma = (theta_vec+beta);
        u = (rho_vec.^2+x.^2).^.5;
        signm = sign(rho_vec) + 1*(rho_vec==0); 
        pts_u = (( u.*cos(gamma).*signm) + 1j*( u.*sin(gamma).*signm)).';
        dist_u = abs(pts_u);  
        
        % GENERATING TIER 1 NODES
        rad_c = sqrt(rand(Nc,1))*r; 
        phi_c = rand(Nc,1)*(2*pi); 
        xc =  rad_c.*cos(phi_c);
        yc =  rad_c.*sin(phi_c);
        pts_c = xc+1j*yc; 
        dist_c = abs(pts_c); 
        
        % PREALLOCATING 
        pr_decouple = zeros(2, total_v);
        sir_decouple = zeros(2, total_v);

        for iL = 1: numl(iS)+1
            nowu_pre = sum(numu(1:iL-1));
            nowu_aft = sum(numu(iL+1:end));
            nowv_pre = sum(numv(1:iL-1));
            nowv_aft = sum(numv(iL+1:end));
            nowu = numu(iL);
            nowv = numv(iL);
            for iV = 1:nowv

                % DL
                % SBS
                shad_u0 = 10.^(.1*(ln_mu20 + ln_sig20*randn(nowu,1))); % typical line     
                shad_u0u = [10.^(.1*(ln_mu21 + ln_sig21*randn(nowu_pre,1)));...
                          shad_u0;...
                          10.^(.1*(ln_mu21 + ln_sig21*randn(nowu_aft,1)))];
                % MBS
                shad_c = 10.^(.1*(ln_mu1 + ln_sig1*randn(length(dist_c), 1))); % 
                
                % GAINS
                bf_gains = [slg_2*ones(nowu_pre,1);
                           mlg_2*ones(nowu,1);
                           slg_2*ones(nowu_aft,1)];
                      
                % CHANNEL PARAMETER
                chnl_u = [gamrnd(m21, 1/m21, nowu_pre,1)
                          gamrnd(m20, 1/m20, nowu, 1); 
                          gamrnd(m21, 1/m21, nowu_aft,1) ];
                chnl_c = gamrnd(m1, 1/m1, length(dist_c), 1); 
                      
               
                [DLval,DLind] = max([p2*mlg_2*b2*shad_u0.*(abs(pts_v(nowv_pre+iV)-pts_u(nowu_pre+1:nowu_pre+nowu))).^(-alpha2);...
                    p1*mlg_1*b1*shad_c.*(abs(pts_v(nowv_pre+iV)-pts_c)).^(-alpha1)]);
                if DLind <= nowu  
                    pr_decouple(1,nowv_pre+iV)= 1;% index 1 for DL, 2 for UL; value 1 for sbs, 2 for mbs
                    sig_pwr = p2*bf_gains(DLind+nowu_pre)*shad_u0u(DLind+nowu_pre)*chnl_u(DLind+nowu_pre).*...
                              abs(pts_v(nowv_pre+iV)-pts_u(nowu_pre+DLind)).^(-alpha2);
                    total_pwr = [p2*bf_gains.*shad_u0u.*chnl_u.*(abs(pts_v(nowv_pre+iV)-pts_u)).^(-alpha2);...
                                     p1*mlg_1.*chnl_c.*shad_c.*(abs(pts_v(nowv_pre+iV)-pts_c)).^(-alpha1)] ;
                    sir = sig_pwr/sum( [total_pwr(1:find(total_pwr==sig_pwr)-1);total_pwr(find(total_pwr==sig_pwr)+1:end)]);
                    sir2 = sig_pwr/(sum(total_pwr) - sig_pwr);
                    sir_decouple(1,nowv_pre+iV)= sir;
                else
                    pr_decouple(1,nowv_pre+iV)= 2;
                    sig_pwr = p1*mlg_1*shad_c(DLind-nowu).*chnl_c(DLind-nowu).*(abs(pts_v(nowv_pre+iV)-pts_c(DLind-nowu))).^(-alpha1) ;   
                    total_pwr = [p2*bf_gains.*shad_u0u.*chnl_u.*(abs(pts_v(nowv_pre+iV)-pts_u)).^(-alpha2);...
                                     p1*mlg_1.*chnl_c.*shad_c.*abs(pts_v(nowv_pre+iV)-pts_c).^(-alpha1)];
                    sir = sig_pwr/(sum(total_pwr) - sig_pwr);
                    sir_decouple(1,nowv_pre+iV)= sir;  
                end
                

                % UL
                UL_shad_u0u = [10.^(.1*(ln_mu21 + ln_sig21*randn(nowv_pre,1)));...
                          10.^(.1*(ln_mu20 + ln_sig20*randn(nowv,1)));...
                          10.^(.1*(ln_mu21 + ln_sig21*randn(nowv_aft,1)))];
                UL_shad_c = 10.^(.1*(ln_mu1 + ln_sig1*randn(length(pts_v), 1)));
                
                % GAINS
                bf_gains = [v_s1*ones(nowv_pre,1);
                           v_s0*ones(nowv,1);
                           v_s1*ones(nowv_aft,1)];

                % CHANNEL PARAMETERS 
                chnl_u = [gamrnd(m21, 1/m21, nowv_pre,1)
                          gamrnd(m20, 1/m20, nowv, 1); 
                          gamrnd(m21, 1/m21, nowv_aft,1) ];
                chnl_c = gamrnd(m1, 1/m1, length(pts_v), 1);
                      
                [ULval,ULind] = max([pv*v_s0*shad_u0.*(abs(pts_v(nowv_pre+iV)-pts_u(nowu_pre+1:nowu_pre+nowu))).^(-alpha2);...
                    pv*v_m*shad_c.*(abs(pts_v(nowv_pre+iV)-pts_c)).^(-alpha1)]);
                if ULind <= nowu  
                    pr_decouple(2,nowv_pre+iV)= 1;
                    sig_pwr = pv*bf_gains(iV+nowv_pre).*UL_shad_u0u(iV+nowv_pre).*chnl_u(iV+nowv_pre).*...
                                 abs(pts_v(nowv_pre+iV)-pts_u(nowu_pre+ULind)).^(-alpha2);
                    total_pwr = pv*bf_gains.*UL_shad_u0u.*chnl_u.*(abs(pts_v-pts_u(nowu_pre+ULind))).^(-alpha2);
                    sir = sig_pwr/sum([total_pwr(1:find(total_pwr==sig_pwr)-1);total_pwr(find(total_pwr==sig_pwr)+1:end)]);
                    sir_decouple(2,nowv_pre+iV)= sir;
                else
                    pr_decouple(2,nowv_pre+iV)= 2;
                    sig_pwr = pv*v_m*UL_shad_c(nowv_pre+iV)*chnl_c(nowv_pre+iV)*...  
                             (abs(  pts_v(nowv_pre+iV) - pts_c(ULind-nowu))  ).^(-alpha1);
                    total_pwr = pv*v_m.*UL_shad_c.*chnl_c.*(abs(pts_v-pts_c(ULind-nowu))).^(-alpha1);         
                    sir = sig_pwr/(sum(total_pwr) - sig_pwr);
                    sir_decouple(2,nowv_pre+iV)= sir;   
                end
            end  
        end  
        
        [case1, posti_case1] = find_case(pr_decouple , 2, 2); % d m   u-M
        [case2, posti_case2] = find_case(pr_decouple , 2, 1); % d-M   u-S
    %    [case3 posti_case3] = find_case(pr_decouple , 1, 2); % d-S   u-M
        [case4, posti_case4] = find_case(pr_decouple , 1, 1); % d-S   u-S

        pr_case(1:3,iS) = [case1;case2;case4]./ (length(posti_case1)+length(posti_case2)+length(posti_case4));
        sir_case(1:6,iS) = [sum(log(1+sir_decouple(1,posti_case1)))/length(posti_case1);...
                            sum(log(1+sir_decouple(2,posti_case1)))/length(posti_case1);...
                            sum(log(1+sir_decouple(1,posti_case2)))/length(posti_case2);...
                            sum(log(1+sir_decouple(2,posti_case2)))/length(posti_case2);...
                            sum(log(1+sir_decouple(1,posti_case4)))/length(posti_case4);
                            sum(log(1+sir_decouple(2,posti_case4)))/length(posti_case4)];
    end
    
    % Means
    pr_NLOS_decouple(1:3,i) = mean(pr_case,2); 
    sir_NLOS_decouple(1:6,i) = mean(sir_case,2);
end

figure
se_decouple_NLOS = ((sir_NLOS_decouple(1,1:4:n)+sir_NLOS_decouple(2,1:4:n)).*pr_NLOS_decouple(1,1:4:n)...
                +(sir_NLOS_decouple(3,1:4:n)+sir_NLOS_decouple(4,1:4:n)).*pr_NLOS_decouple(2,1:4:n)+...
                (sir_NLOS_decouple(5,1:4:n)+sir_NLOS_decouple(6,1:4:n)).*pr_NLOS_decouple(3,1:4:n));
plot(0.1:0.4:n/10,se_decouple_NLOS);

xlim([0 6.5]);
xlabel("λs/λm")
ylabel("Spectral Efficiency (nats/HZ)")
if senario == 0
    ylim([0 4]);
    title("LOS")
else
    ylim([1 9]);
    title("NLOS")
end

function [num ,ind] = find_case(A,a,b)
  ii = 1;
  ind = [];
  for i=1:length(A)
      if A(1,i) == a && A(2,i)==b
         ind(ii)=i;
         ii= ii+1;
      end
  end
  num =  ii -1;
end
