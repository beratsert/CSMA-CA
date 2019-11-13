%-----------A.Berat Sert 2015401069--------------%
%-------------2015401069 Boun EE-----------------%
%------------Senior Design Project---------------%

clc;
clear;
close all;

%---CONFIGS---%
Time=0;                     % Process Time
CW_min = 16;                % Minumum Contention Window value
CW = 16;                    % Default value of Contention Window
CW_max = 1024;              % Maximum Contention Window value
CW_sim_max = 0;            % The maximum value for Contention Window in the simulation
Total_data = 0;             % A variable that includes the total data transfer.
Simulation_count=0;         % A variable that counts event number.
Time_Array=[];              % Time vector for data transfer events
Total_Data_Transfer=[];     % Throughput vector for data transfer events
data_duration=0;            % The calculated duration for data in each event.
Collision_count=0;          % The number of total collision
Collision = [];             % Collision matrix which includes when n-th collision happened.

% unit : sec
SlotTime = 9*(10^-6);       % Minimum Time Interval that can be sensed
SIFS = 16*(10^-6);          % Short Interframe Space
DIFS = 34*(10^-6);          % DIFS = SIFS + (2 * Slot Time)
RTS = 100*(10^-6);           % Request to Send
CTS = 34*(10^-6);           % Clear to Send
ACK = 28*(10^-6);           % Acknowledgement Packet

DataRate = 866.7*(10^6)/8; % 160Mhz 256-QAM Channels' max speed is 866.7 Mbps 
MaximumDataPacket = 5500*8; % Maximum Data Size for a data packet in Bytes
Clients = {};
BackOff = [];
plot_tmp=0;
plot_time=0;
plot_data=0;
Collision_data = 0;
data_1st_user = 0;
data_tmp=0;
collision_tmp=0;
%Reads credentials from input_network.txt
fid = fopen('input_network.txt');
while ~feof(fid)                                         % Writes to the Clients which is a cell array.      
      credentials = fgetl(fid);
      Clients = [Clients {strsplit(credentials, ' ')}];  % Parses the data 
end
fclose(fid);


NumberOfUser = size(Clients,2);                          % Calculates the total number of user
BackOff = zeros(1, NumberOfUser);                        % includes the backoff numbers of each user with their index number.
Transfer_Counter = zeros(1,NumberOfUser);
Contention_Window(1:NumberOfUser)= CW;                   % Contention Window for each users

while 1
    data_duration = 0;  
    for i = 1:NumberOfUser
        
        if (Time > str2double(Clients{1,i}{1,1})) && (BackOff(1,i) == 0) && (0 ~= str2double(Clients{1,i}{1,2}))
        BackOff(1,i) = randi([1,Contention_Window(i)]); 
        end
        
    end
    
    tmp=0;
    
    for i = 1:NumberOfUser
        if (Time < str2double(Clients{1,i}{1,1}))
            BackOff(1, i) = 2*CW_max;
        end
        if (str2double(Clients{1,i}{1,2})== 0 )
            BackOff(1, i) = 2*CW_max;
            tmp=tmp+1;
        end
    end
    if tmp==NumberOfUser
        fprintf("\n%d User completed its transmission!\n" ,NumberOfUser);
        break;
    end
    
    Time = Time + DIFS;
    [Min_BackOff, IndexOfMin] = min(BackOff);
    Time = Time + Min_BackOff*SlotTime+DIFS;
    count=0;
    Max_Index= zeros(1,NumberOfUser);
    
    for i = 1: size(BackOff,2)
        
        if Min_BackOff == BackOff(1,i)
            count = count +1 ;
            Max_Index(count) = i;
        end
    end
    
    if count>1 && Min_BackOff ~= 2*CW_max %if there is a collision, it increases the CW.
        for i= 1:count 
            if CW_max/2+1 > Contention_Window(Max_Index(i))
                Contention_Window(Max_Index(i))=2*Contention_Window(Max_Index(i)); 
            
            if CW_sim_max< Contention_Window(i)
               CW_sim_max = Contention_Window(i);
            end
            
            end
        end
        
        Time = Time + RTS + DIFS + CTS;
        Collision_count = Collision_count + 1;
        Collision( Collision_count,1) = Collision_count;      
        Collision( Collision_count,2) = Time;
        
    elseif count == NumberOfUser && Min_BackOff == 2*CW_max
        
        Time = Time + RTS + DIFS + CTS;

    else 
        Contention_Window(IndexOfMin) = CW_min;
%         if Transfer_Counter(IndexOfMin) == 0
%             Transfer_Counter(IndexOfMin) = Transfer_Counter(IndexOfMin) + 1;
%             Clients_Transfer{1,IndexOfMin}{Transfer_Counter(IndexOfMin),1} = 0;
%             Clients_Transfer{1,IndexOfMin}{Transfer_Counter(IndexOfMin),2} = 0;
%         end
%         
        if MaximumDataPacket > (10^6)*str2double(Clients{1,IndexOfMin}(2))
            
            data_duration = str2double(Clients{1,IndexOfMin}(2))/DataRate;
            Clients{1,IndexOfMin}(2) = {'0'}; 
            fprintf("%d. User has completed its transfer at %f\n",IndexOfMin, Time);
%             Transfer_Counter(IndexOfMin) = Transfer_Counter(IndexOfMin) + 1;
%             Clients_Transfer{1,IndexOfMin}{Transfer_Counter(IndexOfMin),1} = Time; 
%             Clients_Transfer{1,IndexOfMin}{Transfer_Counter(IndexOfMin),2} = Clients_Transfer{1,IndexOfMin}{Transfer_Counter(IndexOfMin)-1,2}  + data_duration*DataRate;
        
        else
            
            data_duration = MaximumDataPacket/DataRate;
            Clients{1,IndexOfMin}{1,2} = num2str( str2double(Clients{1,IndexOfMin}(2)) - MaximumDataPacket*1e-6);
%             Transfer_Counter(IndexOfMin) = Transfer_Counter(IndexOfMin) + 1;
%             Clients_Transfer{1,IndexOfMin}{Transfer_Counter(IndexOfMin),1} = Time; 
%             Clients_Transfer{1,IndexOfMin}{Transfer_Counter(IndexOfMin),2} = Clients_Transfer{1,IndexOfMin}{Transfer_Counter(IndexOfMin)-1,2}  + data_duration*DataRate;
            if IndexOfMin == 1
                
                data_1st_user = data_1st_user + MaximumDataPacket;
                
            end
        end
        
        Time = Time + RTS + SIFS + CTS + SIFS + data_duration + SIFS + ACK;
        
    end
    
    BackOff(1,1:NumberOfUser) = BackOff(1,1:NumberOfUser) - Min_BackOff;
    Simulation_count = Simulation_count +1;
    Total_data=Total_data + data_duration*DataRate;
    Time_Array(Simulation_count)=Time;
    Total_Data_Transfer(Simulation_count)=Total_data;
    
    if  plot_time < Time 
        plot_tmp = plot_tmp+1;
        if plot_time == 0
            throughput2(plot_tmp,2) = 0;
            throughput2(plot_tmp,1) = 0;
            collision2(plot_tmp,2) = 0;
            collision2(plot_tmp,1) = 0;
            transfer_1st_user(plot_tmp,1)=0;
            transfer_1st_user(plot_tmp,2)=0;
        else
            throughput2(plot_tmp,2) = 8*(Total_data - plot_data)/(Time-Time1);
            throughput2(plot_tmp,1) = Time;
            collision2(plot_tmp,1) = Time;
            collision2(plot_tmp,2) = (Collision_count-collision_tmp)/(Time-Time1);
            transfer_1st_user(plot_tmp,1) = Time;
            transfer_1st_user(plot_tmp,2) = 8*(data_1st_user-data_tmp)/(Time-Time1);
        end
        plot_data = Total_data;
        plot_time = plot_time + 1;
        Time1 = Time;
        data_tmp = data_1st_user;
        collision_tmp = Collision_count;
    end
    if (Time > 300)
        break
    end
end
% subplot(2,1,1);
% plot(1:17,throughput2/1e6);

%plot(transfer_1st_user(1:plot_tmp,1),transfer_1st_user(1:plot_tmp,2)/1e6,'*-b');
%xlabel('The number of clients')
%ylabel('The throughput of the first user in Mbps')


% plot(Time_Array,Throughput_Array./1e6,'.');
% xlim([0 Time]);
% xlabel('Time')
% ylabel('Data in MB')
% title('Total Transferred Data vs Time')



%plot(Collision(:,2),Collision(:,1));
%xlim([0 Time]);
%xlabel('Time')
%ylabel('the Number of collution')
%title('the Number of Collution vs Time')

% 
% Throughput_of_1st_User(1,1) = 0;
% Throughput_of_1st_User(1,2) = 0;
% for i = 2:Transfer_Counter(1)
%     Throughput_of_1st_User(i,1) = Clients_Transfer{1,1}{i,1};
%     Throughput_of_1st_User(i,2) = (Clients_Transfer{1,1}{i,2} - Clients_Transfer{1,1}{i-1,2})/(Clients_Transfer{1,1}{i,1}-Clients_Transfer{1,1}{i-1,1}); 
%     
% end
% plot(Throughput_of_1st_User(:,1),Throughput_of_1st_User(:,2)*8/1e6);
% xlim([0 Time]);
% ylim([0 1000]);
% xlabel('Time')plot(throughput2(1:plot_tmp,1),throughput2(1:plot_tmp,2)/1e6,'*-b');
%xlim([0 throughput2(plot_tmp,1)]);
%ylim([300 500]);
%xlabel('Time(s)')
%ylabel('Throughput in Mbps')
% ylabel('Throughput in Mbps')
% title('the Throughput of 1st User vs Time')

% Throughput(1,1) = 0;
% Throughput(1,2) = 0;
% for i = 2:fix(Simulation_count)
%     j=i;
%     Throughput(i,1) = Time_Array(j);
%     Throughput(i,2) = (Total_Data_Transfer(j) -Total_Data_Transfer(j-1))/(Time_Array(j)-Time_Array(j-1)); 
%      
% end
% plot(Throughput(:,1),Throughput(:,2)*8/1e6,'o--r');
% subplot(2,1,1);
plot(throughput2(1:plot_tmp,1),throughput2(1:plot_tmp,2)/1e6,'*-b');
xlim([0 throughput2(plot_tmp,1)]);
xticks([0;5;10;15;20;25;30;35;40;45;50;55;60;65;70;75;80;85;90;95;100;105;110;115;120;125;130;135;140;145;150;155;160;165;170;175;180;185;190;195;200;205;210;215;220;225;230;235;240;245;250;255;260;265;270;275;280;285;290;295;])
xticklabels({  '1'   '2'   '3'   '4'   '5'   '6'   '7'   '8'   '9'  '10'  '11'  '12'  '13'  '14'  '15'  '16'  '17'  '18'  '19' '20'   '21'   '22'   '23'   '24'   '25'   '26'   '27'   '28'   '29'  '30' '31'   '32'   '33'   '34'   '35'   '36'   '37'   '38'   '39'  '40'  '41'  '42'  '43'  '44'  '45'  '46'  '47'  '48'  '49' '50'   '51'   '52'   '53'   '54'   '55'   '56'   '57'   '58'   '59'  '60'
})
ylim([300 500]);
xlabel('Time(s)')
ylabel('Throughput in Mbps')
% 
% subplot(2,1,2);
% plot(collision2(1:plot_tmp,1),collision2(1:plot_tmp,2),'*-b');
% xlim([0 throughput2(plot_tmp,1)]);
% xticks([0;5;10;15;20;25;30;35;40;45;50;55;60;65;70;75;80;85;90;95;100;105;110;115;120;125;130;135;140;145])
% xticklabels({  '1'   '2'   '3'   '4'   '5'   '6'   '7'   '8'   '9'  '10'  '11'  '12'  '13'  '14'  '15'  '16'  '17'  '18'  '19' '20'   '21'   '22'   '23'   '24'   '25'   '26'   '27'   '28'   '29'  '30'
% })
% xlabel('Number of Clients')
% ylabel('Number of Collision per sec')
% xlim([0 throughput2(plot_tmp,1)]);
% xlabel('Time(s)');
% ylabel('Total Collision');