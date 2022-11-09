classdef element
    properties
        freq = 5e9;
        lambda 
        dist 
        array_size
        coord
        array_pos
        surface
        gain
        array_phase
        rcs = 1;
    end
    methods
        function obj = element(varargin) %constructor to initialize class
            %varargin - [array_size,coord,surface,gain,rcs]
            obj.array_size = varargin{1};
            obj.coord = varargin{2};
            obj.surface = varargin{3};
            obj.lambda = 3e8/obj.freq;
            obj.dist = obj.lambda/2;
            obj.gain = varargin{4};
            obj.array_phase = ones(obj.array_size);
            if nargin>4
                obj.rcs = varargin{5};
            end
        end
        function array_pos = array_creator(obj)
            %% creating the array from first coordinates
            M = obj.array_size(1);
            N = obj.array_size(2);
            array_pos = zeros(M,N,3); %xyz coordinates for rows and columns
            if obj.surface
                x=3;y=2;z=1;  
            else
                x=1;y=2;z=3;  
            end
            array_pos(:,:,z) = obj.coord(z)*ones(M,N); %Z coordinate is same across the array
            for i = 1:N
                array_pos(:,i,x) = obj.coord(x)+(-2*obj.surface+1)*(i-1)*obj.dist; %keeping reference at left of the start
                for j = 1:M
                    array_pos(j,:,y) = obj.coord(y)+(j-1)*obj.dist;
                end
            end
        end
        function H_val = H_response(varargin)
            obj1 = varargin{1};
            obj2 = varargin{2};
            intermediate = varargin{3};
            is_WB = 0;
            if nargin>3
                is_WB = 1;
%                 ofdm_data = varargin{4};
                OFDM = varargin{4}; %obj for ofdm
            end
        %H_rx_ant = zeros(size(obj2.array_phase));
        [tx_row, tx_col] = size(obj1.array_phase);
        [rx_row, rx_col] = size(obj2.array_phase);
%         h = zeros(size(obj2.array_phase));
        row_channel = 0;
        for i = 1:rx_row
            for j = 1:rx_col
%                 for k = 1:tx_row
%                     for l = 1:tx_col
                        row_channel = row_channel+1;
                        d = obj1.finding_seperation(obj2.array_pos(i,j,:));
                        phase_multiplier = reshape(obj1.array_phase.',1,numel(d));
                        d_mat = reshape(d',1,numel(d));
%                         if (obj2.is_multipath == 1)
%                          h(row_channel,:) = phase_multiplier.*(obj1.gain*obj1.lambda*exp(-1i*2*pi*d_mat/obj1.lambda)./(4*pi*d_mat));   
                        if is_WB
                            k = -OFDM.N/2+1:OFDM.N/2;
                            for scs = 1:OFDM.N
                                h(row_channel,:,scs) = exp(-1i*2*pi*d_mat*OFDM.del_f*k(scs)/3e8).*phase_multiplier.*(obj1.rcs*obj1.gain*obj2.gain*obj1.lambda*exp(-1i*2*pi*d_mat/obj1.lambda)./(4*pi*d_mat)); 
                            end
                        else
                            h(row_channel,:) = phase_multiplier.*(obj1.rcs*obj1.gain*obj2.gain*obj1.lambda*exp(-1i*2*pi*d_mat/obj1.lambda)./(4*pi*d_mat));
                        end
%                 end
            end
        end
        if intermediate
            H_val = h;
        else
            H_val = sum(sum(h.*obj2.array_phase));
        end
        end
        function d = finding_seperation(obj,rx)
            d = sqrt((obj.array_pos(:,:,1)-rx(:,:,1)).^2+(obj.array_pos(:,:,2)-rx(:,:,2)).^2+(obj.array_pos(:,:,3)-rx(:,:,3)).^2);
        end
        function arr_phase = getphase(obj)
            arr_phase = ones(obj.array_size);
        end
        function beam_phase = compute_phase(SS,AP,UE,bypass_AP)
            if isequal(SS,AP)
                [N,M,~] = size(SS.array_pos);
            else
            [M, N, ~] = size(SS.array_pos);
            end
            if bypass_AP
                dist_travel = SS.finding_seperation(UE.array_pos(1,1,:));
            else
                dist_travel = SS.finding_seperation(AP.array_pos(1,1,:))+SS.finding_seperation(UE.array_pos(1,1,:));
            end
            if M > 1
                if isequal(SS,AP)
                    extra_path = dist_travel(1,:)-dist_travel(1,1);
                else
%                 extra_path = mean(dist_travel(2:end,1)-dist_travel(1:end-1,1));
                  extra_path = dist_travel(:,1)-dist_travel(1,1);
                end
            else
                extra_path = 0;
            end
            alpha = asind(extra_path/SS.dist);
            if isequal(SS,AP)
                beam_phase = zeros(N,M);
            else
            beam_phase = zeros(M,N);
            end
            for i = 1:M
%           phase_shift_array(:,i) = exp(1i*2*pi*(i-1)*dist*cosd(theta)/lambda);
%                 beam_phase(i,:) = exp(1i*2*pi*(M-i)*SS.dist*cosd(alpha+90)/SS.lambda);
                if isequal(SS,AP)
                    beam_phase(:,i) = exp(1i*2*pi*extra_path(i)/SS.lambda);
                else
                beam_phase(i,:) = exp(1i*2*pi*extra_path(i)/SS.lambda);
                end
            end
        end 
        function angle_meas = angle_between_UE(UE1,UE2)
            angle_meas = acosd(dot(UE1.coord,UE2.coord)/(norm(UE1.coord)*norm(UE2.coord)));
        end
      end
 end  