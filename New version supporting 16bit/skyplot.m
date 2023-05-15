function [] = skyplot(ids,az,el)
%SKYPLOT
SAT_SYS_IDX = ['G','R','E','C','J','I','S'];
syscolor = {'.r' '.g' '.b' '.y' '.m' '.k' '.c'};%a color for each sys
title_satnr = '';
figure;
handle = axes;
hseries_skyplot_sats = [];
htext_dop = [];
htext_skyplot = cell(size(SAT_SYS_IDX,2),1);
plot_polar(handle);
if ~isempty(az)
    a=90-az;%make azimuth point 0 to North (up)
    r=0.01*abs(el-90);%map elevation to length (0 is 90º)
    svx = r.*cos(a*pi/180);
    svy = r.*sin(a*pi/180);
    lblidx = 0;
    totalsats = 0;
    for i=1:length(SAT_SYS_IDX)
        csys = ids(:,1)==SAT_SYS_IDX(i);
        satnr = sum(csys);
        if satnr~=0
            lblidx = lblidx + 1;
            csvx = svx(csys);
            csvy = svy(csys);
            sats = ids(csys,:);
            for j=1:satnr
                idx = [];
                if ~isempty(hseries_skyplot_sats)
                    idx = find(sum(hseries_skyplot_sats==sats(j,:),2)==3);
                end
                if isempty(idx)
                    %add new sat to  list
                    idx = size(hseries_skyplot_sats,1)+1;
                    hseries_skyplot_sats(idx,:) = sats(j,:);
                    hseries_skyplot(idx) = cell(1);
                end
                satsvisible(idx) = true;
                if ~isempty(hseries_skyplot{idx})
                    hseries_skyplot{idx}.XData = [hseries_skyplot{idx}.XData csvx(j)]; hseries_skyplot{idx}.YData = [hseries_skyplot{idx}.YData csvy(j)];
                else
                    hp = plot(handle,csvx(j),csvy(j),syscolor{i},'MarkerSize',14);
                    hseries_skyplot{idx} = hp;
                end
            end

            if ~isempty(htext_skyplot{i})
                delete(htext_skyplot{i});
            end
            htext_skyplot{i} = text(handle,csvx+0.01,csvy+0.01,ids(csys,:));
            title_satnr = strcat(title_satnr,SAT_SYS_IDX(i),'#',num2str(satnr),', ');
            totalsats = totalsats + satnr;
        else
            if ~isempty(htext_skyplot{i})
                delete(htext_skyplot{i});
            end
        end
    end
    if any(~satsvisible)
        points = hseries_skyplot(~satsvisible);
        hseries_skyplot(~satsvisible) = [];
        hseries_skyplot_sats(~satsvisible,:) = [];
        ptnr = size(points,1);
        for i=1:ptnr
            delete(points{i});
        end
    end
    if ~isempty(htext_dop)
        delete(htext_dop);
    end
    title_satnr = strcat(title_satnr,'sat#',num2str(totalsats));
    htext_dop = text(handle,-0.5,-1.1,title_satnr);
end

end

