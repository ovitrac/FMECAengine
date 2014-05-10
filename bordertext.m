% Usage:
%       htext = bordertext( 'Position', 'YourText');
%       htext = bordertext( 'Position', 'YourText', '[FormatString]' );
%
%       run 'bordertext' without arguments to see an example.
%   
% Position:
%                'topleft'           'top'          'topright'
%       'lefttop'|-------------------------------------------|'righttop'
%                |                                           |
%          'left'|                  'center'                 |'right'
%                |                                           |
%    'leftbottom'|-------------------------------------------|'rightbottom'
%                'bottomleft'       'bottom'     'bottomright'
% 
%       This 'Position'-string can be prefixed with 'inner' or 'figure'.
%       e.g: 'innertopleft' or 'figurebottomright'
%
%       'Position' can also be an [ X Y ]-position wrt the current axis.
%
% Text:
%       Your own string:    'Hello world'
%       or with more lines: 'Hello World\nNice weather\nisn''t it?'
%       or with more lines: {'Hello World','Nice weather','isn''t it?'}
%
% FormatString:
%       FontSize:           [%f] (any positive floatingpoint value)
%       FontAttributes:     [B]  [I]  [L]  [T] (Bold, Italic, Latex, Tex)
%       FontColors:          [w]  [k]  [r]  [g]  [b]  [c]  [m]  [y]
%       BackgroundColors:   [bw] [bk] [br] [bg] [bb] [bc] [bm] [by]
%       EdgeColors:         [ew] [ek] [er] [eg] [eb] [ec] [em] [ey]
%       TextRotation:       [R%f] (any floatingpoint value)    
%       HorizontalAlignment:[Hl] [Hc] [Hr]  (normally you don't need these)
%       VerticalAlignment:  [Vb] [Vm] [Vt]  (normally you don't need these)
%
% To remove a bordertext-string from the figure:
%       bordertext( 'Position', '' );
%
% Examples:
%       bordertext( 'innertopleft', '1^r^s^t line\n2^n^d line.', '[B][I][T][r][bg][ek][14][R10]' );
%   or  bordertext( 'innertopleft', '1^r^s^t line\n2^n^d line.', 'B I T r bg ek 14 R10' );
%       (Bold, Italic, Tex-interpreter, red chars, green background, black edge, fontsize=14, rotation=10degr)
% 
%       to remove: bordertext( 'innertopleft', '' );
%
% Last release:
% V 2.05  16-04-2013  Make '[ X Y ]'-text invisible outside axis.
%
% Author:   Jan G. de Wilde     email: jan.dewilde@nl.thalesgroup.com
%           Thales Nederland B.V.

% V 1.00  22-09-2006  First release.
% V 2.00  01-04-2013  Text stays on axis, also while zooming, panning or resizing.
% V 2.01  04-04-2013  Fix for reverse X- and Y-axes.
% V 2.02  05-04-2013  Some bugfixes.
% V 2.03  10-04-2013  More flexible in interpreting 'FormatString'.
% V 2.04  12-04-2013  Now able to handle linked-axes.
% V 2.05  16-04-2013  Make '[ X Y ]'-text invisible outside axis.

function h = bordertext( pos, txt, FormatStr )

if nargin == 0, % Help, Test and Debug mode
    doc bordertext
    
    figure(999);	clf('reset');
    plot( randn(1,100), '+' );	axis([ 0 100 -3 3 ]);   pan on;
%     set(gca,'XDir','reverse');   set(gca,'YDir','reverse');
%     imagesc(rand(25),[-1 1]); colormap(gray); axis auto tight square
    title('title');   xlabel('xlabel');   ylabel('ylabel');

    bordertext('topleft','topleft');
    bordertext('top','top');
    bordertext('topright','topright');
    bordertext('righttop','righttop');
    bordertext('right','right');
    bordertext('rightbottom','rightbottom');
    bordertext('bottomright','bottomright');
    bordertext('bottom','bottom');
    bordertext('bottomleft','bottomleft');
    bordertext('leftbottom','leftbottom');
    bordertext('left','left');
    bordertext('lefttop','lefttop');
    bordertext('center','center');
    
    bordertext('innertopleft','innertopleft');
    bordertext('innertop','innertop');
    bordertext('innertopright','innertopright');
    bordertext('innerrighttop','innerrighttop');
    bordertext('innerright','innerright');
    bordertext('innerrightbottom','innerrightbottom');
    bordertext('innerbottomright','innerbottomright');
    bordertext('innerbottom','innerbottom');
    bordertext('innerbottomleft','innerbottomleft');
    bordertext('innerleftbottom','innerleftbottom');
    bordertext('innerleft','innerleft');
    bordertext('innerlefttop','innerlefttop');
    
    bordertext('figuretopleft','figuretopleft');
    bordertext('figuretop','figuretop');
    bordertext('figuretopright','figuretopright');
    bordertext('figurerighttop','figurerighttop');
    bordertext('figureright','figureright');
    bordertext('figurerightbottom','figurerightbottom');
    bordertext('figurebottomright','figurebottomright');
    bordertext('figurebottom','figurebottom');
    bordertext('figurebottomleft','figurebottomleft');
    bordertext('figureleftbottom','figureleftbottom');
    bordertext('figureleft','figureleft');
    bordertext('figurelefttop','figurelefttop');
    bordertext('figurecenter','figurecenter');
    
    bordertext([50 2],'h = bordertext( ''Position'', ''Text'' [,''FormatString'' ])','[r][14][Hc]');

    return
end

ax      = axis;
curr_ax = gca;
fnr     = gcf;
Xrev    = strcmp( get( curr_ax, 'XDir' ), 'reverse' );   % V 2.01
Yrev    = strcmp( get( curr_ax, 'YDir' ), 'reverse' );

if nargin == 1, % This is the CallBack-function while zooming, panning or resizing
    if strcmp( pos, 'callback' ),
        
        % from GetLinkedAxes() ...
        tmp = getappdata( curr_ax );    % V 2.04
        if isfield( tmp, 'graphics_linkaxes' ),
            Targets = get( tmp.graphics_linkaxes, 'Targets' );
            
            linked_axes = NaN( 1, length( Targets ) );
            for iTargets = 1 : length( Targets ),
                tmpChildren             = get( Targets( iTargets ), 'Children' );
                linked_axes( iTargets ) = get( tmpChildren(1), 'Parent' );
            end
        else
            linked_axes = curr_ax;
        end

        for curr_ax = linked_axes,    % V 2.04
            UserData    = get( curr_ax, 'UserData' );
            if ~isfield( UserData, 'BordertextHandles' ),   return;   end
            
            bt_handles  = UserData.BordertextHandles;
            if isfield( bt_handles, 'userdefined' ),    % V 2.05
                for htxt = UserData.BordertextHandles.userdefined,
                    if ishandle( htxt )
                        pos = get( htxt, 'Position' );
                        if pos(1) < ax(1) || pos(1) > ax(2) || ...
                                pos(2) < ax(3) || pos(2) > ax(4),
                            set( htxt, 'Visible', 'off' );
                        else
                            set( htxt, 'Visible', 'on' );
                        end
                    end
                end
            end
            
            AllFields = fieldnames( bt_handles );
            for iflds = 1 : length( AllFields ),
                FieldName = AllFields( iflds ); FieldName = FieldName{1};
                th = UserData.BordertextHandles.( FieldName );
                if ~ishandle( th ),   continue;   end
                
                if strncmpi( FieldName, 'inner', 5 ),
                    FieldName = FieldName( 6 : end );
                end
                
                Xnew = NaN;   Ynew = NaN;   Znew = 0;
                switch FieldName,
                    case 'topleft',     Xnew = ax(1+Xrev);          Ynew = ax(4-Yrev);
                    case 'top',         Xnew = (ax(1)+ax(2))/2;     Ynew = ax(4-Yrev);
                    case 'topright',    Xnew = ax(2-Xrev);          Ynew = ax(4-Yrev);
                    case 'righttop',    Xnew = ax(2-Xrev);          Ynew = ax(4-Yrev);
                    case 'right',       Xnew = ax(2-Xrev);          Ynew = (ax(3)+ax(4))/2;
                    case 'rightbottom', Xnew = ax(2-Xrev);          Ynew = ax(3+Yrev);
                    case 'bottomright', Xnew = ax(2-Xrev);          Ynew = ax(3+Yrev);
                    case 'bottom',      Xnew = (ax(1)+ax(2))/2;     Ynew = ax(3+Yrev);
                    case 'bottomleft',  Xnew = ax(1+Xrev);          Ynew = ax(3+Yrev);
                    case 'leftbottom',  Xnew = ax(1+Xrev);          Ynew = ax(3+Yrev);
                    case 'left',        Xnew = ax(1+Xrev);          Ynew = (ax(3)+ax(4))/2;
                    case 'lefttop',     Xnew = ax(1+Xrev);          Ynew = ax(4-Yrev);
                    case 'center',      Xnew = (ax(1)+ax(2))/2;     Ynew = (ax(3)+ax(4))/2;
                    otherwise,          % skip
                end
                if ~isnan( Xnew ),   set( th, 'Position', [ Xnew Ynew Znew ] );   end
            end
        end
    else
        fprintf('Usage: htext = bordertext( ''Position'', ''YourText'', ''[FormatString]'' );\n');
        error('Too few arguments.');
    end
    
    return
end

if nargin == 2,   FormatStr = [];   end

% Make once in each figure an invisisible axes for text on 'figure'-level.
UserData = get( fnr, 'UserData' );
if ~isfield( UserData, 'BordertextHandles' ) || ~ishandle( UserData.BordertextHandles ),
    UserData.BordertextHandles.h_axes = axes( 'Unit', 'Normalized', ...
        'Position', [0 0 1 1], 'Visible', 'off', 'Nextplot', 'add', 'Hittest','off');
    set( fnr, 'UserData', UserData, 'CurrentAxes', curr_ax );
end
h_axes = UserData.BordertextHandles.h_axes;

% Some defaults
ALIGN_H     = '';
ALIGN_V     = '';
BOLD        = false;
ITALIC      = false;
INTERPRETER = 'none';
BCOL        = 'none';
ECOL        = 'none';
FCOL        = [ 0 0 0 ];
ROTATION    = '';

if ~ischar( pos ),
   UserPos = pos;
   pos     = 'userdefined';
end

if strncmpi( pos, 'figure', 6 ),    FontSize = 12;
else                                FontSize =  8;
end

%% Interpret/decode the FormatString
if ~isempty( FormatStr ),
    FormatStr( FormatStr == '[' ) = ' ';
    FormatStr( FormatStr == ']' ) = ' ';
    
    SpacePos = find( FormatStr == ' ' );
    SpacePos = [ 0 SpacePos length(FormatStr)+1 ];
    
    for i = 1 : length( SpacePos ) - 1,        
        cmd = FormatStr( SpacePos(i)+1 : SpacePos(i+1)-1 ) ;
        cmd( isspace( cmd ) ) = '';
        if isempty( cmd ), 	 continue;   end

        switch cmd,            
            case 'B',   BOLD     = true;        % Bold
            case 'I',   ITALIC   = true;        % Italic
            case 'L',   INTERPRETER = 'Latex';  % Latex interpreter
            case 'T',   INTERPRETER = 'Tex';    % Tex interpreter
                
            case 'w',   FCOL = [ 1 1 1 ];       % FontColors
            case 'k',   FCOL = [ 0 0 0 ];
            case 'r',   FCOL = [ 1 0 0 ];
            case 'g',   FCOL = [ 0 1 0 ];
            case 'b',   FCOL = [ 0 0 1 ];
            case 'c',   FCOL = [ 0 1 1 ];
            case 'm',   FCOL = [ 1 0 1 ];
            case 'y',   FCOL = [ 1 1 0 ];
                
            case 'bw',  BCOL = [ 1 1 1 ];       % BackgroundColors
            case 'bk',  BCOL = [ 0 0 0 ];
            case 'br',  BCOL = [ 1 0 0 ];
            case 'bg',  BCOL = [ 0 1 0 ];
            case 'bb',  BCOL = [ 0 0 1 ];
            case 'bc',  BCOL = [ 0 1 1 ];
            case 'bm',  BCOL = [ 1 0 1 ];
            case 'by',  BCOL = [ 1 1 0 ];
                
            case 'ew',  ECOL = [ 1 1 1 ];       % EdgeColors
            case 'ek',  ECOL = [ 0 0 0 ];
            case 'er',  ECOL = [ 1 0 0 ];
            case 'eg',  ECOL = [ 0 1 0 ];
            case 'eb',  ECOL = [ 0 0 1 ];
            case 'ec',  ECOL = [ 0 1 1 ];
            case 'em',  ECOL = [ 1 0 1 ];
            case 'ey',  ECOL = [ 1 1 0 ];
                
            case 'Hl',	ALIGN_H = 'left';       % Alignments
            case 'Hc',  ALIGN_H = 'center';
            case 'Hr',  ALIGN_H = 'right';
                
            case 'Vb',  ALIGN_V = 'bottom';
            case 'Vm',  ALIGN_V = 'middle';
            case 'Vt',  ALIGN_V = 'top';
                
            otherwise,
                if cmd(1) == 'R',               % Rotation
                    ROTATION =   str2double( cmd( 2 : end ) );
                    continue;
                end
                
                tmp = abs( str2double( cmd ) ); % FontSize
                if isnan( tmp ),
                    fprintf('%s :: Unrecognized substring: ''%s''\n', ...
                        mfilename, cmd );
                else
                    FontSize = tmp;
                end
        end
    end
end

% Coordinates in HiddenAxis
if strncmpi( pos, 'figure', 6 ),
    set( fnr, 'CurrentAxes', h_axes );
    xpos_left = 0.0;    xpos_middle = 0.5;    xpos_right  = 1.0;
    ypos_top  = 1.0;    ypos_middle = 0.5;    ypos_bottom = 0.0;
end

% Handle multi-line
if ~iscell( txt )
    ind = strfind( txt, '\n' );
    if ~isempty( ind ),
        str = [];
        ind = [ -1 ind length(txt)+1 ];
        for i = 2 : length( ind ),
            str{i-1} = txt(ind(i-1)+2:ind(i)-1); %#ok<AGROW>
        end
        txt = str;   clear str;
    end
end

% Translate 'pos' to [X Y]-coordinates
switch lower( pos ),
    case 'topleft',
        h = text( ax(1+Xrev), ax(4-Yrev), txt );
        set( h, 'HorizontalAlignment', 'left');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 0 );
        
    case 'top',
        h = text( (ax(1)+ax(2))/2, ax(4-Yrev), txt );
        set( h, 'HorizontalAlignment', 'center');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 0 );
        
    case 'topright',
        h = text( ax(2-Xrev), ax(4-Yrev), txt );
        set( h, 'HorizontalAlignment', 'right');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 0 );
        
    case 'righttop',
        h = text( ax(2-Xrev), ax(4-Yrev), txt );
        set( h, 'HorizontalAlignment', 'right');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 90 );
        
    case 'right',
        h = text( ax(2-Xrev), (ax(3)+ax(4))/2, txt );
        set( h, 'HorizontalAlignment', 'center');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 90 );
        
    case 'rightbottom',
        h = text( ax(2-Xrev), ax(3+Yrev), txt );
        set( h, 'HorizontalAlignment', 'left');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 90 );
        
    case 'bottomright',
        h = text( ax(2-Xrev), ax(3+Yrev), txt );
        set( h, 'HorizontalAlignment', 'right');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 0 );
        
    case 'bottom',
        h	= text( (ax(1)+ax(2))/2, ax(3+Yrev), txt );
        set( h, 'HorizontalAlignment', 'center');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 0 );
        
    case 'bottomleft',
        h = text( ax(1+Xrev), ax(3+Yrev), txt );
        set( h, 'HorizontalAlignment', 'left');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 0 );
        
    case 'leftbottom',
        h = text( ax(1+Xrev), ax(3+Yrev), txt );
        set( h, 'HorizontalAlignment', 'left');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 90 );
        
    case 'left',
        h = text( ax(1+Xrev), (ax(3)+ax(4))/2, txt );
        set( h, 'HorizontalAlignment', 'center');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 90 );
        
    case 'lefttop',
        h = text( ax(1+Xrev), ax(4-Yrev), txt );
        set( h, 'HorizontalAlignment', 'right');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 90 );
        
    case 'center',
        h = text( (ax(1)+ax(2))/2, (ax(3)+ax(4))/2, txt );
        set( h, 'HorizontalAlignment', 'center');
        set( h, 'VerticalAlignment'  , 'middle');
        set( h, 'Rotation'           , 0 );
        
    case 'innertopleft',
        h = text( ax(1+Xrev), ax(4-Yrev), txt );
        set( h, 'HorizontalAlignment', 'left');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 0 );
        
    case 'innertop',
        h = text( (ax(1)+ax(2))/2, ax(4-Yrev), txt );
        set( h, 'HorizontalAlignment', 'center');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 0 );
        
    case 'innertopright',
        h = text( ax(2-Xrev), ax(4-Yrev), txt );
        set( h, 'HorizontalAlignment', 'right');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 0 );
        
    case 'innerrighttop',
        h = text( ax(2-Xrev), ax(4-Yrev), txt );
        set( h, 'HorizontalAlignment', 'right');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 90 );
        
    case 'innerright',
        h = text( ax(2-Xrev), (ax(3)+ax(4))/2, txt );
        set( h, 'HorizontalAlignment', 'center');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 90 );
        
    case 'innerrightbottom',
        h = text( ax(2-Xrev), ax(3+Yrev), txt );
        set( h, 'HorizontalAlignment', 'left');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 90 );
        
    case 'innerbottomright',
        h = text( ax(2-Xrev), ax(3+Yrev), txt );
        set( h, 'HorizontalAlignment', 'right');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 0 );
        
    case 'innerbottom',
        h = text( (ax(1)+ax(2))/2, ax(3+Yrev), txt );
        set( h, 'HorizontalAlignment', 'center');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 0 );
        
    case 'innerbottomleft',
        h = text( ax(1+Xrev), ax(3+Yrev), txt );
        set( h, 'HorizontalAlignment', 'left');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 0 );
        
    case 'innerleftbottom',
        h = text( ax(1+Xrev), ax(3+Yrev), txt );
        set( h, 'HorizontalAlignment', 'left');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 90 );
        
    case 'innerleft',
        h = text( ax(1+Xrev), (ax(3)+ax(4))/2, txt );
        set( h, 'HorizontalAlignment', 'center');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 90 );
        
    case 'innerlefttop',
        h = text( ax(1+Xrev), ax(4-Yrev), txt );
        set( h, 'HorizontalAlignment', 'right');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 90 );

    case 'figuretopleft',
        h = text( xpos_left, ypos_top, txt );
        set( h, 'HorizontalAlignment', 'left');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 0 );
        
    case 'figuretop',
        h = text( xpos_middle, ypos_top, ax(4), txt );
        set( h, 'HorizontalAlignment', 'center');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 0 );
        
    case 'figuretopright',
        h = text( xpos_right, ypos_top, txt );
        set( h, 'HorizontalAlignment', 'right');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 0 );
        
    case 'figurerighttop',
        h = text( xpos_right, ypos_top, [ txt ' ' ] );
        set( h, 'HorizontalAlignment', 'right');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 90 );
        
    case 'figureright',
        h = text( xpos_right, ypos_middle, txt );
        set( h, 'HorizontalAlignment', 'center');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 90 );
        
    case 'figurerightbottom',
        h = text( xpos_right, ypos_bottom, txt );
        set( h, 'HorizontalAlignment', 'left');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 90 );
        
    case 'figurebottomright',
        h = text( xpos_right, ypos_bottom, txt );
        set( h, 'HorizontalAlignment', 'right');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 0 );
        
    case 'figurebottom',
        h = text( xpos_middle, ypos_bottom, txt );
        set( h, 'HorizontalAlignment', 'center');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 0 );
        
    case 'figurebottomleft',
        h = text( xpos_left, ypos_bottom, txt );
        set( h, 'HorizontalAlignment', 'left');
        set( h, 'VerticalAlignment'  , 'bottom');
        set( h, 'Rotation'           , 0 );
        
    case 'figureleftbottom',
        h = text( xpos_left, ypos_bottom, txt );
        set( h, 'HorizontalAlignment', 'left');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 90 );
        
    case 'figureleft',
        h = text( xpos_left, ypos_middle, txt );
        set( h, 'HorizontalAlignment', 'center');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 90 );
        
    case 'figurelefttop',
        h = text( xpos_left, ypos_top, [ txt ' ' ] );
        set( h, 'HorizontalAlignment', 'right');
        set( h, 'VerticalAlignment'  , 'top');
        set( h, 'Rotation'           , 90 );
        
    case 'figurecenter',
        h = text( xpos_middle, ypos_middle, txt );
        set( h, 'HorizontalAlignment', 'center');
        set( h, 'VerticalAlignment'  , 'middle');
        set( h, 'Rotation'           , 0 );
        
    case 'userdefined',
        h = text( UserPos(1), UserPos(end), txt );
        
    otherwise,
        fprintf('%s :: Unrecognized ''Position'': ''%s''\n', ...
                        mfilename, pos );
        clear h;   return;
end

% Set Text-Properties
if BOLD,     set( h, 'FontWeight', 'Bold' );   end

if ITALIC,   set( h, 'FontAngle', 'Italic' );   end

set( h, 'Interpreter', INTERPRETER );

set( h, 'Color',           FCOL );
set( h, 'BackGroundColor', BCOL );
set( h, 'EdgeColor',       ECOL );

if ~isempty( ALIGN_H  ),   set( h, 'HorizontalAlignment', ALIGN_H  );   end
if ~isempty( ALIGN_V  ),   set( h, 'VerticalAlignment',   ALIGN_V  );   end

if ~isempty( ROTATION ),   set( h, 'Rotation',            ROTATION );   end

if ~isnan( FontSize ),   set( h, 'FontSize', FontSize );   end

set( fnr, 'CurrentAxes', curr_ax ); % Set scope to current axis

% Store handles in the figure of in the subplot(axis)
if strcmpi( pos, 'userdefined' ),
    UserData = get( curr_ax, 'UserData' );    % V 2.05
    if isfield( UserData, 'BordertextHandles' ),
        if ~isfield( UserData.BordertextHandles, 'userdefined' ),
            UserData.BordertextHandles.userdefined = h;
        else
            UserData.BordertextHandles.userdefined( length( UserData.BordertextHandles.userdefined ) + 1 ) = h;
        end
    else
        UserData.BordertextHandles.userdefined = h;
    end
    set( curr_ax, 'UserData', UserData );
else
    if strncmpi( pos, 'figure', 6 ),    h_curr = fnr;
    else                                h_curr = curr_ax;
    end
    UserData = get( h_curr, 'UserData' );
    if isfield( UserData, 'BordertextHandles' ),
        if isfield( UserData.BordertextHandles, pos ),
            if ishandle( UserData.BordertextHandles.(pos) ),
                delete( UserData.BordertextHandles.(pos) );
            end
            UserData.BordertextHandles = rmfield( UserData.BordertextHandles, pos );
        end
    end
    UserData.BordertextHandles.(pos) = h;
    set( h_curr, 'UserData', UserData );
end

% Set Callback functions for zooming and panning
set( zoom( fnr ), 'actionpostcallback',     'bordertext(''callback'')');
set( pan(  fnr ), 'actionpostcallback',     'bordertext(''callback'')');
set(       fnr, 'WindowButtonMotionFcn',    'bordertext(''callback'')');
set(       fnr, 'ResizeFcn',                'bordertext(''callback'')');

if nargout == 0;   clear h;   end

return