function data = gcfd
%  GCFD get current figure data
% INRA\Olivier Vitrac 28/10/02 - rev. 09/11/09

% REVISION HISTORY
% 29/10/02 release candidate
% 09/11/09 add text

propAXESlist  = {'sXlabel','sYlabel','nXscale','nYscale','nXlim','nYlim'};
propLINElist  = {'Color','LineStyle','LineWidth','Marker','MarkerSize','MarkerEdgeColor','MarkerFaceColor','Xdata','Ydata','Zdata'};
propIMAGElist = {'Xdata','Ydata','CData','AlphaData','CDataMapping','AlphaDataMapping','tag'};
propTEXTlist  = {'BackgroundColor','Color','EdgeColor','FontAngle','FontName','FontSize','FontUnits','FontWeight','HorizontalAlignment','Interpreter',...
    'LineStyle', 'LineWidth','Margin','Position','Rotation','String','Units','VerticalAlignment'};
typeOBJECTlist = {'line','image','text'};
excludedTAG = {'legend','Colorbar'};
data = [];
set(gcf,'units','normalized','position',[0 0 .5 .5])

hc = get(gcf,'children');
indisaxes = find(strcmp(get(hc,'type'),'axes'));
if any(indisaxes)
    n = length(indisaxes);
    disp(sprintf('...%d axes objects found in %d children',n,length(hc)))
else
    disp('no axes objects'), return
end
i = 0;
% For each axes
for hi = flipud(hc(indisaxes))'
    tagname = get(hi,'tag');
    if ismember(tagname,excludedTAG)
        disp(sprintf('AXES [%d] = ''%s'' is excluded',i+1,tagname))
    else
        i = i+1; disp(sprintf('AXES [%d]',i))
        hp = flipud(get(hi,'children'));
        % Axes properties
        for p = propAXESlist
            if any(get(hi,p{1}(2:end)))
                if p{1}(1)=='s', data = setfield(data,{i},p{1}(2:end),get(get(hi,p{1}(2:end)),'string'));
                elseif p{1}(1)=='n', data = setfield(data,{i},p{1}(2:end),get(hi,p{1}(2:end)));
                end
            else
                data = setfield(data,{i},p{1}(2:end),[]);
            end
        end
        % Object properties
        for typeOBJECT = typeOBJECTlist
            switch lower(typeOBJECT{1})
            case 'line' , propOBJECTlist = propLINElist;
            case 'image', propOBJECTlist = propIMAGElist;
            case 'text' , propOBJECTlist = propTEXTlist;
            end
            indisobject = find(strcmp(get(hp,'type'),typeOBJECT));
            m = length(indisobject);
            if any(indisobject)
                disp(sprintf('...%d ''%s'' objects found in %d children',m,typeOBJECT{1},length(hp)));
                j = 0;
                for hl=hp(indisobject)'
                    j= j+1;
                    for p = propOBJECTlist
                        value = get(hl,p{1});
                         if ~isempty(value)
                            data = setfield(data,{i},typeOBJECT{1},{j},p{1},value);
                        else
                             data = setfield(data,{i},typeOBJECT{1},{j},p{1},[]);
                         end
                    end
                end
            end
        end
    end
end
