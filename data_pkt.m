classdef data_pkt
    properties
        BW; %Bandwidth
        N; %subcarriers
        data; %data bits in BPSK
        del_f; %sub-carrier spacing
        data_bins; %data bin for OFDM
        d; %LOS delay to be compensated for
        variance;
        noise_power;
        NF = 5.5;
        smoothening = 1;
        smooth_filt = [1 1 1]/3;
        ref_pow = 30;
        A;
    end
    methods
        function obj = data_pkt(varargin)
            obj.BW = varargin{1};
            obj.N = varargin{2};
            obj.del_f = obj.BW/obj.N;
            obj.data = 2*randi([0 1],1,obj.N)-1;
            obj.noise_power = -174+10*log10(obj.BW)+obj.NF;
            obj.variance = 50e-3*10^(obj.noise_power/10);
%             obj.variance = 0;
%             obj.data_bins = [-obj.N/2+1:obj.N/2];
        end
        function d = initial_delay(Tx,Rx)
            d = norm(Tx.coord-Rx.coord);
        end
        function H_est = channel_est(varargin)
            OFDM = varargin{1};
            H = varargin{2};
            Tx = varargin{3};
            Rx = varargin{4};
            
            %% create ofdm data plus noise
            [row,col,scs] = size(H);
            H_est = zeros(row,col,OFDM.N);
            for i = 1:row
                for j = 1:col
                    if scs == 1
                        temp(:) = H(i,j)*ones(1,OFDM.N);
                    else
                    temp(:) = H(i,j,:);
                    end
                    rx_data = temp.*OFDM.data;
                    data_time = OFDM.A*ifft(rx_data)+sqrt(OFDM.variance/2)*(randn(size(rx_data))+1i*randn(size(rx_data)));
                    H_est(i,j,:) = fft(data_time).*OFDM.data;
                    if OFDM.smoothening
                        temp(:) = H_est(i,j,:);
                        smooth_data = conv(temp,OFDM.smooth_filt);
                        H_est(i,j,:) = [temp(1:2) smooth_data(3:end-2)];
                    end
                        
                end
            end
        end
        function W = ZF_precoder(varargin)
            OFDM = varargin{1};
            N_users = varargin{2};
            H_est = [];
            for i = 1:N_users
                H_est = [H_est;varargin{i+2}];
            end
%             H_est = varargin{1};
            for k = 1:OFDM.N
                W(:,:,k) = (H_est(:,:,k)'*H_est(:,:,k))\H_est(:,:,k)';
            end
        end
        
        function A = Tx_gain(obj)
            rms_volt = rms(ifft(obj.data));
            pow_db = 10*log10(rms_volt^2/(50e-3));
            A = 10^((obj.ref_pow-pow_db)/20);
        end
        
%         function H = channel_estimate
        
            
    end
end