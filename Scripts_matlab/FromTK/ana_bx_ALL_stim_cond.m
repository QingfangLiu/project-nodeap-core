clear all
close all 
root = '/Users/daria/Desktop/NODEAP_Behavioral_Data_reupload';

subjects = {'NODEAP_06', 'NODEAP_07', 'NODEAP_08', 'NODEAP_09', 'NODEAP_10', 'NODEAP_12', 'NODEAP_13', 'NODEAP_15' ...
'NODEAP_16', 'NODEAP_17', 'NODEAP_24', 'NODEAP_26', 'NODEAP_30', 'NODEAP_34', 'NODEAP_35', 'NODEAP_36' ...
'NODEAP_39', 'NODEAP_41', 'NODEAP_43', 'NODEAP_44', 'NODEAP_45', 'NODEAP_46', 'NODEAP_48', 'NODEAP_51' ...
'NODEAP_53', 'NODEAP_55', 'NODEAP_56', 'NODEAP_59', 'NODEAP_60', 'NODEAP_63', 'NODEAP_64', 'NODEAP_65' ...
'NODEAP_66', 'NODEAP_68', 'NODEAP_70', 'NODEAP_73', 'NODEAP_74', 'NODEAP_76', 'NODEAP_77', 'NODEAP_78' ...
'NODEAP_79', 'NODEAP_80', 'NODEAP_83', 'NODEAP_84', 'NODEAP_85', 'NODEAP_86','NODEAP_87', 'NODEAP_88'}


nsubs = length(subjects);
nsess = 3;
nruns = 5; % conditioning


% 1 = stim sham, 2= sham stim, 3 = sham sham
stim_cond = [3,2,1; 3,1,2; 1,3,2; 3,2,1; 3,2,1; 1,3,2; 2,3,1,; 2,3,1; ...
    2,3,1; 2,1,3; 3,2,1; 2,1,3; 2,1,3; 1,3,2; 1,2,3; 1,2,3; ... 
    3,2,1; 2,1,3; 1,2,3; 1,2,3; 3,1,2; 2,1,3; 1,2,3; 3,2,1; ... 
    2,3,1; 3,1,2; 3,2,1; 2,1,3; 1,3,2; 3,1,2; 1,3,2; 1,2,3; ...
    2,1,3; 1,3,2; 3,1,2; 1,2,3; 3,1,2; 2,1,3; 1,3,2; 2,3,1; ...
    2,3,1; 2,3,1; 3,1,2; 3,1,2; 1,3,2; 1,2,3; 3,2,1; 2,3,1];


cols = [ones(1,3)'*[1:nsubs]/nsubs]'-0.15;
%cols = [[1 0 1]; [0 1 0]; [0 1 1]; [0 0 0]];
% sw_deval = [2,1];
% sa_deval = [1,2];

deval = [2,1; 1,2; 2,1; 1,2; 1,2; 1,2; 2,1; 2,1; ...
    1,2; 1,2; 2,1; 2,1; 2,1; 2,1; 1,2; 1,2; ...
    1,2; 1,2; 1,2; 2,1; 2,1; 1,2; 2,1; 2,1; ...
    1,2; 2,1; 1,2; 2,1; 1,2; 2,1; 1,2; 2,1; ...
    2,1; 2,1; 1,2; 2,1; 2,1; 1,2; 1,2; 2,1; ...
    1,2; 2,1; 1,2; 2,1; 2,1; 1,2; 2,1; 2,1];
%deval = [1,2; 2,1; 2,1; 1,2]; % 1=non-devalued, 2=devalued; (1=SW 2=SA)
%deval = [2,1,2;2,1,2;1,2,1;1,2,1]; %1 = SW
%deval_mid = [2,1; 1,2; 1,2; 2,1];
%deval = [1,2,1; 2,1,2; 2,1,2; 1,2,1] % [1,2,1] = SA devalued for sess 1 & 3
deval_mid = [1,2; 2,1; 1,2; 2,1; 2,1; 2,1; 1,2; 1,2; ...
    2,1; 2,1; 1,2; 1,2; 1,2; 1,2; 2,1; 2,1; ...
    2,1; 2,1; 2,1; 1,2; 1,2; 2,1; 1,2; 1,2; ...
    2,1; 1,2; 2,1; 1,2; 2,1; 1,2; 2,1; 1,2; ... 
    1,2; 1,2; 2,1; 1,2; 1,2; 2,1; 2,1; 1,2; ...
    2,1; 1,2; 2,1; 2,1; 1,2; 2,1; 1,2; 1,2];



for s=1:nsubs
    
    sn = subjects{s};
  
    for sess=1:nsess
        
        % load data
        load(fullfile(root, sn, sprintf('complete_DAY2_POST_%s_SESS_%02d.mat',sn,sess)));
        
        % Pleasantenss ratings SCREEN
        d = res.initial_pleas_ratings;
        % 1     selected odor (1=sweet, 2=savory)
        % 2     onset of cue onset from t0
        % 3     onset of odor trigger from t0
        % 4     onset of sniff cue from t0
        % 5     pleasantness rating
        % 6     rating RT (ms)
        % 7     time rating submitted from t0
        % 8     rating scale flipped? (0=no, 1=yes);
        
        for o = 1:8
            pleasRat_initial(s,o) = nanmean(d(d(:,1)==o,5));
        end
        
        disp([res.odornames{res.odors(1),2} ' and ' res.odornames{res.odors(2),2}])

        
        
        % Pleasantenss ratings SCREEN
        d = res.pleasantness_ratings_SCREEN;
        % 1     selected odor (1=sweet, 2=savory)
        % 2     onset of cue onset from t0
        % 3     onset of odor trigger from t0
        % 4     onset of sniff cue from t0
        % 5     pleasantness rating
        % 6     rating RT (ms)
        % 7     time rating submitted from t0
        % 8     rating scale flipped? (0=no, 1=yes)
            
        for o = 1:2
            pleasRat_screen(s,deval(s,o)) = nanmean(d(d(:,1)==o,5));
        end
        
        d = res.intensity_ratings;
        for o = 1:2
          intensRat_screen(s,deval(s,o)) = nanmean(d(d(:,1)==o,5));
        end
        
        % conditioning data
%         for r=1:nruns
%             d = res.conditioning{r,sess};
%             % 1     cue pair (1:16; 1:4 = SW/0 A, 5:8 = SW/0 set B, 9:12 = SA/0 set A, 13:16 = SA/0 set B)
%             % 2     odor identity (1=sweet; 2=savory)
%             % 3     odor set (1=A, 2=B (B is not probed on DAY2 at PRE))
%             % 4     left/right (1=odor cue on left, i.e. not flipped, 2=flipped l/r)
%             % 5     choice (1=left, 2=right)
%             % 6     odor chosen? (1=yes, 0=no)
%             % 7     choice stimuli onset from t0
%             % 8     response time (from onset of choice)
%             % 9     odor trigger onset from t0
%             % 10    sniff cue onset from t0
%             % 11    cue # on the left
%             % 12    cue # on the right
%             % 13    selected cue #
%             % 14    selected odor (0=none, 1=sweet, 2, savory))
%          if floor(sess/2)==sess/2            
%            for o=1:2
%                 cond(s,r,deval_mid(s,o),sess) = nanmean(d(d(:,2)==o,6));
%            end
%          else
%             for o=1:2
%                 cond(s,r,deval(s,o),sess) = nanmean(d(d(:,2)==o,6));
%             end
%          end
%         end
%         
        
                % Pleasantenss ratings PRE and POST
                d = res.pleasantness_ratings_PRE{sess};
            if floor(sess/2)==sess/2
               
                for o = 1:2
                    pleasRat_pre(s,o,sess) = nanmean(d(d(:,1)==deval_mid(s,o),5));
                end
                
                d = res.pleasantness_ratings_POST{sess};
                for o = 1:2
                    pleasRat_post(s,o,sess) = nanmean(d(d(:,1)==deval_mid(s,o),5));
                end
            else
                for o = 1:2
                    pleasRat_pre(s,o,sess) = nanmean(d(d(:,1)==deval(s,o),5));
                end
                
                d = res.pleasantness_ratings_POST{sess};
                for o = 1:2
                    pleasRat_post(s,o,sess) = nanmean(d(d(:,1)==deval(s,o),5));
                end
            end
            
                
                        
        
                % choice task PRE and POST
                % 1     cue pair (1:16; 1:4 = SW/0 set A, 5:8 = SA/0 set A, 9:16 = SW/SA set A)
                % 2     odor identity (1=sweet; 2=savory; 3=both S+/+)
                % 3     odor set (1=A, 2=B (B is not shown on DAY2 at PRE))
                % 4     left/right (1=odor cue on left, i.e. not flipped, 2=flipped l/r)
                % 5     choice (1=left, 2=right)
                % 6     odor chosen? (1=yes, 0=no for S+/-; 1=sweet,2=savory for S+/+)
                % 7     choice stimuli onset from t0
                % 8     response time (from onset of choice)
                % 9     odor trigger onset from t0
                % 10    sniff cue onset from 
                % 11    cue # on the left
                % 12    cue # on the right
                % 13    selected cue #
                % 14    selected odor (0=none, 1=sweet, 2, savory))
                
                d = res.choice_PRE{sess};
                            
            if sess == 2 %%% mod 
                if deval_mid(s,1) == 2
                    d(d(:,2)==3,6) = -(d(d(:,2)==3,6)-1)+2; % recode 1,2=SW,SA to 1,2=NON-DEVAL,DEVAL
                end
                choice_pre_post(s,1,sess) = nanmean(d(d(:,2)==3,6))-1;
                
                choice_pre_ab(s,1,sess) = nanmean(d(d(:,2)==3 & d(:,3)==1,6)); %Set A
                choice_pre_ab(s,2,sess) = nanmean(d(d(:,2)==3 & d(:,3)==2,6)); %Set B
        
                
                d = res.choice_POST{sess};
                if deval_mid(s,1) == 2
                    d(d(:,2)==3,6) = -(d(d(:,2)==3,6)-1)+2; % recode 1,2=SW,SA to 1,2=NON-DEVAL,DEVAL
                end
                
            else
               
                if deval(s,1) == 2
                    d(d(:,2)==3,6) = -(d(d(:,2)==3,6)-1)+2; % recode 1,2=SW,SA to 1,2=NON-DEVAL,DEVAL
                end
                choice_pre_post(s,1,sess) = nanmean(d(d(:,2)==3,6))-1;
                
                choice_pre_ab(s,1,sess) = nanmean(d(d(:,2)==3 & d(:,3)==1,6)); %Set A
                choice_pre_ab(s,2,sess) = nanmean(d(d(:,2)==3 & d(:,3)==2,6)); %Set B
        
                
                d = res.choice_POST{sess};
                if deval(s,1) == 2
                    d(d(:,2)==3,6) = -(d(d(:,2)==3,6)-1)+2; % recode 1,2=SW,SA to 1,2=NON-DEVAL,DEVAL
                end
            end
            
                choice_pre_post(s,2,sess) = nanmean(d(d(:,2)==3,6))-1;

                choice_post_ab(s,1,sess) = nanmean(d(d(:,2)==3 & d(:,3)==1,6)); %Set A
                choice_post_ab(s,2,sess) = nanmean(d(d(:,2)==3 & d(:,3)==2,6)); %Set B  
        
    end
end
% pleasantness data all odors
figure, ah = axes;
hold on
bh = bar(mean(pleasRat_initial)); set(bh, 'BarWidth', 0.4, 'FaceColor', [.75 .75 .75])
eh = errorbar(mean(pleasRat_initial), std(pleasRat_initial)/sqrt(nsubs)); set(eh, 'LineStyle',  'none', 'LineWidth',2)
set(ah,'xLim', [0.5,8.5], 'XTick', [1:8], 'XTickLabel', {res.odornames{1:8,1}}); xlabel('odors'); ylabel('pleasantness');
title('ALL Initial pleasantness ratings SCREEN')


% selected odors
figure, ah = subplot(1,2,1); hold on
% for s=1:nsubs
%     ph = plot(pleasRat_screen(s,:));
%     set(ph, 'color', cols(s,:), 'LineWidth',2)
% end
bh = bar(mean(pleasRat_screen)); set(bh, 'BarWidth', 0.4, 'FaceColor', [.75 .75 .75])
eh = errorbar(mean(pleasRat_screen), std(pleasRat_screen)/sqrt(nsubs)); set (eh, 'Linestyle',  'none', 'LineWidth',2)
set(ah,'xLim', [0.5,2.5], 'XTick', [1,2], 'XTickLabel', {'NON-DEVAL','DEVAL'}); 
xlabel('odors'); ylabel('pleasantness');
%legend({'sub 1','sub 2', 'sub 3', 'sub 4'},'Location', 'southeast')
title('ALL Pleasantness ratings selected odors SCREEN')

ah = subplot(1,2,2); hold on
% for s=1:nsubs
%     ph = plot(intensRat_screen(s,:));
%     set(ph, 'color', cols(s,:), 'LineWidth',2);
% end
bh = bar(mean(intensRat_screen)); set(bh, 'BarWidth', 0.4, 'FaceColor', [.75 .75 .75])
eh = errorbar(mean(intensRat_screen), std(intensRat_screen)/sqrt(nsubs)); set (eh, 'Linestyle',  'none', 'LineWidth',2)
set(ah,'xLim', [0.5,2.5], 'XTick', [1,2], 'XTickLabel', {'NON-DEVAL','DEVAL'}); 
xlabel('odors'); ylabel('intensity');
%legend({'sub 1','sub 2', 'sub 3', 'sub 4'},'Location', 'southeast')
title('ALL Intensity ratings selected odors SCREEN')

% % Conditioning
% figure, ah = axes; hold on
% cond_av = mean(cond,5);
% % for s=1:nsubs
% %     ph = plot(mean(cond_av(s,:,:),3));
% %     set(ph, 'color', cols(s,:),  'LineWidth', 1);
% % end
% eh1 = errorbar(mean(cond_av(:,:,1)), std(cond_av(:,:,1))./sqrt(nsubs));set(eh1, 'linewidth',2, 'color', 'b');
% eh2 = errorbar(mean(cond_av(:,:,2)), std(cond_av(:,:,2))./sqrt(nsubs));set(eh2, 'linewidth',2, 'color', 'r');
% legend([eh1,eh2],{'NON-DEVAL','DEVAL'}, 'Location', 'southeast'); 
% set(ah,'XTick', 0:6); xlabel('runs'); ylabel('percent odor chosen'); xlim([0 6]);
% % title('ALL Conditioning average')
% 
% figure,
% cond_sess = squeeze(mean(cond,3));
% for sess=1:nsess
%     ah = subplot(1,3,sess); hold on
% %     for s=1:nsubs
% %         ph = plot(cond_sess(s,:,sess));
% %         set(ph, 'color', cols(s,:), 'LineWidth',1.5);
% %     end
%     eh = errorbar(mean(cond_sess(:,:,sess)), std(cond_sess(:,:,sess))./sqrt(nsubs));
%     set(eh, 'linewidth',2, 'color', [.75, .75, .75]);
%     set(ah,'XTick', 0:6); xlabel('runs'); ylabel('percent odor chosen');
%     title(sprintf('ALL Conditioning session %d',sess))
%     %legend({'sub 1','sub 2', 'sub 3', 'sub 4', 'avg'}, 'Location', 'southeast')
% end

% Pleasantness ratings PRE vs POST
pleasRat_pre_av = mean(pleasRat_pre,3);
pleasRat_post_av = mean(pleasRat_post,3);
pleasRat_pre_post_av = [pleasRat_pre_av(:,1), pleasRat_post_av(:,1), pleasRat_pre_av(:,2),pleasRat_post_av(:,2)];
figure, ah = subplot(1,2,1); hold on
% for s=1:nsubs
%        ph = plot(1:2, pleasRat_pre_post_av(s,1:2));
%        set(ph, 'color', cols(s,:), 'LineWidth',2);
%        ph = plot(4:5, pleasRat_pre_post_av(s,3:4));
%        set(ph, 'color', cols(s,:), 'LineWidth',2);
% end
bh = bar([1:2,4:5],mean(pleasRat_pre_post_av)); set(bh, 'BarWidth', 0.4, 'FaceColor', [.75 .75 .75])
eh = errorbar([1:2,4:5],mean(pleasRat_pre_post_av), std(pleasRat_pre_post_av)/sqrt(nsubs)); set (eh, 'Linestyle',  'none')
set(ah,'xLim', [0.5,5.5], 'XTick', [1:2,4:5], 'XTickLabel', {'NON pre', 'NON post', 'DEV pre', 'DEV post'}); 
xlabel('odors'); ylabel('pleasantness');
title('ALL Pleasantness ratings PRE vs POST average')


pleasRat_diff_av = pleasRat_post_av-pleasRat_pre_av;
ah = subplot(1,2,2); hold on
% for s=1:nsubs
%     ph = plot([1,2], pleasRat_diff_av(s,:));
%     set(ph, 'color', cols(s,:), 'LineWidth', 2);
% end
bh = bar([1:2],mean(pleasRat_diff_av)); set(bh, 'BarWidth', 0.4, 'FaceColor', [.75 .75 .75])
eh = errorbar([1:2],mean(pleasRat_diff_av), std(pleasRat_diff_av)/sqrt(nsubs)); set (eh, 'Linestyle',  'none', 'LineWidth',2)
set(ah,'xLim', [0.5,2.5], 'XTick', [1:2], 'XTickLabel', {'NON', 'DEV'}); 
%legend({'sub 1','sub 2', 'sub 3', 'sub 4'})
xlabel('odors'); ylabel('pleasantness POST-PRE');
title('ALL Pleasantness ratings PRE vs POST average')

% Pleasantness ratinsg PRE vs POST per SESS
pleasRat_diff = squeeze(pleasRat_post-pleasRat_pre);
figure, ah = axes; hold on
xx=[1:2; 4:5; 7:8];
% for s=1:nsubs
%     for ii=1:nsess
%         ph(s) = plot(xx(ii,:), pleasRat_diff(s,:,ii));
%         set(ph(s), 'color', cols(s,:), 'LineWidth', 2);
%     end
% end
%legend(ph,{'sub 1','sub 2', 'sub 3', 'sub 4'}, 'Location', 'southeast')
bh = bar([1:2, 4:5, 7:8],mean([pleasRat_diff(:,:,1), pleasRat_diff(:,:,2), pleasRat_diff(:,:,3)])); set(bh, 'BarWidth', 0.4, 'FaceColor', [.75 .75 .75])
eh = errorbar([1:2, 4:5, 7:8],mean([pleasRat_diff(:,:,1), pleasRat_diff(:,:,2), pleasRat_diff(:,:,3)]), std([pleasRat_diff(:,:,1), pleasRat_diff(:,:,2), pleasRat_diff(:,:,3)])/sqrt(nsubs)); set (eh, 'Linestyle',  'none', 'LineWidth',2)
set(ah,'xLim', [0.5,8.5], 'XTick', [1:2, 4:5, 7:8], 'XTickLabel', {'NON_1', 'DEV_1', 'NON_2', 'DEV_2', 'NON_3', 'DEV_3'}); 
xlabel('odors'); ylabel('pleasantness POST-PRE');
title('ALL Pleasantness ratings PRE vs POST per session')


% choices PRE vs POST
figure, ah = subplot(1,2,1); hold on;
choice_pre_post_av = mean(choice_pre_post,3);
% for s=1:nsubs
%     ph = plot(1:2, choice_pre_post_av(s,:));
%     set(ph, 'color', cols(s,:), 'LineWidth', 2);
% end
bh = bar(1:2,mean(choice_pre_post_av));
set(bh, 'BarWidth', 0.4, 'FaceColor', [.75 .75 .75]); hold on; plot([0.5,2.5],[0.5,0.5], '--')
eh = errorbar([1:2],mean(choice_pre_post_av), std(choice_pre_post_av)/sqrt(nsubs)); 
% legend({'sub 1','sub 2', 'sub 3', 'sub 4'}, 'Location','southeast')
set(eh, 'Linestyle',  'none', 'LineWidth',2)
set(ah,'xLim', [0.5,2.5], 'XTick', 1:2, 'XTickLabel', {'PRE', 'POST'}); ylabel('Chocie DEVAL > NON-DEVAL');
title('ALL Choice S+/S+ PRE vs POST average')

ah = subplot(1,2,2); hold on;
xx=[1:2; 4:5; 7:8];
% for s=1:nsubs
%     for ii=1:nsess
%         ph(s) = plot(xx(ii,:), choice_pre_post(s,:,ii));
%         set(ph(s), 'color', cols(s,:), 'LineWidth',2);
%     end
% end
xx=[1:2, 4:5, 7:8];
choice_pre_post = [choice_pre_post(:,:,1),choice_pre_post(:,:,2),choice_pre_post(:,:,3)];
bh = bar(xx,mean(choice_pre_post)); 
set(bh, 'BarWidth', 0.4, 'FaceColor', [.75 .75 .75]); hold on; plot([0.5,8.5],[0.5,0.5], '--')
eh = errorbar(xx,mean(choice_pre_post), std(choice_pre_post)/sqrt(nsubs)); 
set(eh, 'Linestyle',  'none', 'LineWidth',2)
%legend(ph,{'sub 1','sub 2','sub 3', 'sub 4'}, 'Location', 'southeast')
set(ah,'xLim', [0.5,8.5], 'XTick', xx, 'XTickLabel', {'PRE_1', 'POST_1','PRE_2', 'POST_2','PRE_3', 'POST_3'}); ylabel('Chocie DEVAL > NON-DEVAL');
title('ALL Choice S+/S+ PRE vs POST per session')

% choices per SET
figure, ah = axes; hold on
choice_pre_post_ab_av = [mean(choice_pre_ab(:,1,:),3), mean(choice_post_ab,3)]-1;
xx= [1:3];
% for s=1:nsubs
%     ph = plot(xx, choice_pre_post_ab_av(s,:));
%     set(ph, 'color', cols(s,:), 'LineWidth',1.5);
% end
bh = bar(xx,mean(choice_pre_post_ab_av)); set(bh, 'BarWidth', 0.4, 'FaceColor', [.75 .75 .75]); hold on; plot([0.5,3.5],[0.5,0.5], '--')
eh = errorbar(xx,mean(choice_pre_post_ab_av), std(choice_pre_post_ab_av)/sqrt(nsubs)); 
set(eh, 'Linestyle',  'none', 'LineWidth',2)
set(ah,'xLim', [0.5,3.5], 'XTick', xx, 'XTickLabel', {'SET A pre', 'SET A post', 'SET B post'}); 
xlabel('odors'); ylabel('Chocie DEVAL > NON-DEVAL');
%legend({'sub 1','sub 2', 'sub 3', 'sub 4'})
title('ALL Choice S+/S- PRE vs POST average per SET')

figure, ah = axes; hold on;
choice_pre_post_ab = [choice_pre_ab(:,1,1), choice_post_ab(:,:,1), choice_pre_ab(:,1,2), choice_post_ab(:,:,2),choice_pre_ab(:,1,3), choice_post_ab(:,:,3)]-1;
xx=[1:3; 5:7; 9:11];
% for s=1:nsubs
%     for ii=1:nsess
%         ph = plot(xx(ii,:), [choice_pre_ab(s,1,ii),choice_post_ab(s,:,ii)]-1);
%         set(ph, 'color', cols(s,:), 'LineWidth',2);
%     end
% end
xx=[1:3, 5:7, 9:11];
bh = bar(xx,mean(choice_pre_post_ab)); 
set(bh, 'BarWidth', 0.4, 'FaceColor', [.75 .75 .75]); hold on; plot([0.5,11.5],[0.5,0.5], '--')
eh = errorbar(xx,mean(choice_pre_post_ab), std(choice_pre_post_ab)/sqrt(nsubs)); 
set(eh, 'Linestyle',  'none', 'LineWidth',2)
set(ah,'xLim', [0.5,11.5], 'XTick', xx, 'XTickLabel', {'A pre_1', 'A post_1', 'B post_1','A pre_2', 'A post_2', 'B post_2','A pre_3', 'A post_3', 'B post 3'}); 
%legend({'sub 1','sub 2','sub 3','sub 4'});
ylabel('Chocie DEVAL > NON-DEVAL');
title('ALL Choice S+/S+ PRE vs POST per SET per session')
