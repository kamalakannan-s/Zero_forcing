clear all
close all
clc
%% set OFDM parameters
N_user = 2;
OFDM1 = data_pkt(20e6,64); %SCS = 512; BW = 20MHz data for UE1
OFDM2 = data_pkt(20e6,64); %data for UE2
OFDM1.A = Tx_gain(OFDM1);
OFDM2.A = Tx_gain(OFDM2);
%% set the location of AP and UE
AP_xloc = 0;
AP_yloc = 0;
dist = 5:5:30;
ang = 0:1:360;
for i = 1:length(dist)
    for j = 1:length(ang)
UE1_xloc = AP_xloc+dist(i)*sind(ang(j));
UE1_yloc = AP_yloc+dist(i)*cosd(ang(j));


UE2_xloc = 0;
UE2_yloc = 5;
G = 1;
G_ss_dB = 4;
G_ss = 10^(G_ss_dB/10);

%% create class and obj for AP and UE
AP = element([1,2],[AP_xloc,AP_yloc,0],0,G);
UE1 = element([1,1],[UE1_xloc,UE1_yloc,0],0,G);
UE2 = element([1,1],[UE2_xloc,UE2_yloc,0],0,G);
AP.array_pos = array_creator(AP);

% AP.array_phase = getphase(AP);
UE1.array_pos = array_creator(UE1);
% UE1.array_phase = getphase(UE1);
UE2.array_pos = array_creator(UE2);
% UE2.array_phase = getphase(UE2);
angle = UE1.angle_between_UE(UE2);

%% getting channel
H1 = H_response(AP,UE1,1,OFDM1);
H2 = H_response(AP,UE2,1,OFDM2);

%% estimate channel
% OFDM1.d = initial_delay(AP,UE1);
% OFDM2.d = initial_delay(AP,UE2);

H1est = channel_est(OFDM1,H1,AP,UE1);
H2est = channel_est(OFDM2,H2,AP,UE2);

%% getting precoder
W = ZF_precoder(OFDM1,N_user,H1est,H2est);

%% input data
% for i = 1:N_user
    X = [OFDM1.A*OFDM1.data;OFDM2.A*OFDM2.data];
    H = [H1;H2];
    for k = 1:OFDM1.N
        tx_data(:,:,k) = W(:,:,k)*X(:,k);
        rx_data(:,:,k) = H(:,:,k)*tx_data(:,:,k);
    end

rx_data1(:) = rx_data(1,1,:);
rx_data2(:) = rx_data(2,1,:);

error1 = rx_data1-OFDM1.data;
error2 = rx_data2-OFDM2.data;
evm1 = sqrt(sum(error1.*conj(error1))/(sum(OFDM1.data.*conj(OFDM1.data))));
evm2 = sqrt(sum(error2.*conj(error2))/(sum(OFDM2.data.*conj(OFDM2.data))));
% SNDR1 = 10*log10(1/evm1^2);
EVM(i,j) = 20*log10(evm1);
SNDR1(i,j) = 10*log10(1/evm1^2);
SNDR2(i,j) = 10*log10(1/evm2^2);
    end
end

[D,Theta] = meshgrid(dist,ang);
X = D.*cosd(Theta);
Y = D.*sind(Theta);
x = reshape(X,1,[]);
y = reshape(Y,1,[]);

figure(6)
surf(X,Y,SNDR1');
shading interp
colorbar
view(2);
hold on;
plot(UE2_yloc,UE2_xloc,'r*');

figure(7)
surf(X,Y,SNDR2');
shading interp
colorbar
view(2);
hold on;
plot(UE2_yloc,UE2_xloc,'r*');
% hold on;
% plot(ss_x2,ss_y2,'r*');
% hold on;
% plot(ue1_x,ue1_y,'k-o');
% hold off


% plot(rx_data2,'.', 'markersize', 8)
% H = [H1;H2];
% H_plus = (H'*H)\H';

% angle = acosd((UE1_xloc*UE2_xloc+UE1_yloc*UE2_yloc)/(norm([UE1_xloc UE1_yloc])*norm([UE2_xloc UE2_yloc])));

