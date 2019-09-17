function [p_wtn, c_wtn, p_btwn, c_btwn]=within_group_box_plot_face_color_exp_nov(groups, boxes, inputs, plot_title, test, ctype)

    x = [];
    g = []; iii=1;num=1;
    for ii=1:size(inputs,1)
        for jj=1:size(inputs,2)            
            if isrow(inputs{ii,jj})
                inputs{ii,jj} = inputs{ii,jj}'; 
            end
            x = vertcat(x, inputs{ii,jj}); 
            for kk = 1:length(inputs{ii,jj})
                g(iii) = num; 
                iii=iii+1;
            end
            if isempty(inputs{ii,jj})
                g(iii) = num; 
                x = vertcat(x, NaN);
                iii=iii+1;
            end
            num = num+1;
        end 
    end
    
    len_box = length(boxes); 
    positions = []; 
    if isempty(boxes) 
        positions = [1, 1.25];%1:1:length(groups); 
        tick_loc = positions;
        data = groups;
        sig_between = 0; 
        draw_leg = 0; 
    else
        for i=1:length(groups)
           positions = [positions, linspace(i-0.25, i+0.25, len_box)]; 
        end
        tick_loc = []; j = 1; 
        for i=1:length(groups)
           tick_loc = [tick_loc, mean(positions(j:j+len_box-1))]; 
           j = j+len_box; 
        end 
        data = boxes; 
        sig_between = 1; 
        draw_leg = 1; 
    end
    
    figure()
    hold on
    
    boxplot(x,g, 'positions', positions);
    set(gca,'xtick',tick_loc)
    set(gca,'xticklabel',groups,'FontSize',14,'FontWeight','bold')
    xt = get(gca, 'XTick');
    set(gca, 'XTickLabel','')                                           % Turn Off X-Labels
    xts = groups;
    fntsz = 20;
    text(xt, [-.05, -.05], xts, 'FontSize', fntsz, 'FontWeight', 'bold', 'HorizontalAlignment','center', 'VerticalAlignment','top')

    color = ['g'];
    color = repmat(color, 1, length(groups)); 
    
    h = findobj(gca,'Tag','Box');
    
    for i=1:length(h)
       patch(get(h(i),'XData'),get(h(i),'YData'),color(i),'FaceAlpha',.6);
    end

    leg_c = get(gca, 'Children');

    [p_wtn, c_wtn] = significance_within(inputs, data, positions, test, ctype);
    p_btwn= NaN; c_btwn = NaN; 
    if sig_between
        [p_btwn, c_btwn] = significance_between(inputs, groups, test, ctype);
    end
        
%     title(plot_title)
    if draw_leg
        hleg1 = legend(leg_c(1:3), boxes, 'FontSize',13, 'Orientation', 'horizontal');
    end
    hold off;
    
end
