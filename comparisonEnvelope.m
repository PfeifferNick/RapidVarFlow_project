clear;clc;close all

%% plot different envelopes for different pipe material


H_min_env_GRP=cell2mat(struct2cell(load('H_min_env_GRP')));
H_min_env_Concrete=cell2mat(struct2cell(load('H_min_env_Concrete')));
H_min_env_Iron=cell2mat(struct2cell(load('H_min_env_Iron')));
H_min_env_steel=cell2mat(struct2cell(load('H_min_env_steel')));
H_max_env_GRP=cell2mat(struct2cell(load('H_max_env_GRP')));
H_max_env_Concrete=cell2mat(struct2cell(load('H_max_env_Concrete')));
H_max_env_Iron=cell2mat(struct2cell(load('H_max_env_Iron')));
H_max_env_steel=cell2mat(struct2cell(load('H_max_env_steel')));

xaxisH = 0:55:1100;

figure(1)
         plot(xaxisH,H_max_env_GRP,'r',xaxisH,H_min_env_GRP,'r'),xlabel ('Pipe Length [m]'),ylabel('Pressue Head H [m]')
        hold on 
          plot(xaxisH,H_max_env_Concrete,'b',xaxisH,H_min_env_Concrete,'b')
          hold on 
          plot(xaxisH,H_max_env_Iron,'g',xaxisH,H_min_env_Iron,'g')
         hold on 
         plot(xaxisH,H_max_env_steel,'m',xaxisH,H_min_env_steel,'m')

          legend('max env. GRP','min env. GRP','max env. Concrete','min env. Concrete','max env. Iron','min env. Iron','max env. Steel','min env. Steel' ,'Location', 'best')%,'Location','southoutside','orientation','horizontal')
         grid on;
       
         print ('-f1','comparison_envelope','-depsc');