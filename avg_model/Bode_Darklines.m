%Author: Jerry Francis
%Date: October 5 2008
%
%Makes lines on the bode plots darker.
%Bode_Darklines( width, hFig )
% width - width of line
% hFig - handle of figure of bode plot to darken
%
%Example: 
% g1 = rss(4,1,1); g2 = rss(5,1,1);  %random state space models
% figure;
% bode(g1,g2);
% legend('g1','g2','Location','NorthWest');
% grid on;
% bode_darklines(2);
%
function Bode_Darklines( width, hFig, varargin )

if( nargin <2 )
    hFig = gcf;
end
if( nargin <1 )
    width = 3;
end

%Child Window Array
hFig_Ch = get(hFig, 'Children');
nhFig_Ch = length(hFig_Ch);  %Number of children

for n = 1:nhFig_Ch
    hhFigCh = hFig_Ch(n);
    hhTyp = get(hhFigCh,'Type');
    disp(hhTyp);
    %Test to find axes
    if( strcmp(hhTyp, 'axes') ) 
        disp(hhFigCh);
        hhFigChCh = get(hhFigCh,'Children');
        disp(hhFigChCh);
        
        %Test to find group
        for nn = 1:length(hhFigChCh)
            hhTyp = get(hhFigChCh(nn),'Type');
            disp(hhTyp);
            if( strcmp(hhTyp,'hggroup'))
                hhFigChChCh = get(hhFigChCh(nn),'Children');
                disp(hhFigChChCh);
                
                %Test to find lines
                for( nnn = 1:length(hhFigChChCh))
                    hhTyp = get(hhFigChChCh(nnn),'Type');
                    disp(hhTyp);
                    if( strcmp(hhTyp, 'line'))
                        %Format line properties
                        disp('Found Bode Curve');
                        set(hhFigChChCh(nnn),'LineWidth',width,varargin{:});
                    end
                end
            end
        end
    end
end

return;
