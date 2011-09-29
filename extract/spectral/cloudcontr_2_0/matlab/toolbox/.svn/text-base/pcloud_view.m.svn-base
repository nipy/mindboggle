% h = POINTLAB(P) point cloud visualizer.
%
% P: a struct where P.points contains 3D point coordinates
%
function h_out = pointlab( P )

% setup
h = gca;
hold on;
 
% modify the position of the axis in the current figure
hg = hgtransform;
setappdata( h, 'hg', hg );
set( h,'Position',[0 0.05 1 1-0.05 ] );
set( gcf, 'color', 'white');

% do not give output if not requested
if nargout == 0
    h_out = [];
else
    h_out = h;
end

%%%%%%%%%%%%%%%%%%%%%%%% menu setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%
menu_view = uimenu('Label','View');
    uimenu( menu_view,'Label','View Samples','Callback',@view_samples_callback);
    uimenu( menu_view,'Label','View Normals','Callback',@view_normals_callback);
        
menu_tools = uimenu('Label','Tools');
    uimenu( menu_tools,'Label','Pick sample', 'Callback', @point_pick_callback);
    uimenu( menu_tools,'Label','Snapshot', 'Callback', @snapshot_callback);
    uimenu( menu_tools,'Label','Record Fly Around', 'Callback', @record_callback);
% save status variables
h_samples = draw_samples();  setappdata(h,'h_samples',h_samples);

% set drawing limits and vis setup
xlim([-1,1]);
ylim([-1,1]);
axis equal
axis off

function view_samples_callback( IGNORE1, IGNORE2 ) %#ok<INUSD>
    if strcmp( get( h_samples, 'Visible' ), 'on')
        set( h_samples, 'Visible', 'off' )
    else
        set( h_samples, 'Visible', 'on' )
    end
end
% this is slighgly different...
% normals are not drawn by default (they would take too long)
% I just draw them if specifically requested.
function view_normals_callback( IGNORE1, IGNORE2 )  %#ok<INUSD>
    for i=1:size(P.points,1)
        myline( P.points(i,:), P.normals(i,:), .1, 'parent', hg );
    end
end

function point_pick_callback( IGNORE1, IGNORE2 ) %#ok<INUSD>
    % convert selection point in homogeneous and apply transform
	disp( find_closest( h, P.points ) );
end

function snapshot_callback( IGNORE1, IGNORE2 ) %#ok<INUSD>
    outfilename = sprintf('results/snapshot_%s_%s', datestr(now,'ddmmyy'), datestr(now,'HHMMss') );
    print('-dpng', outfilename);
end
function record_callback( IGNORE1, IGNORE2 ) %#ok<INUSD>
    disp('- modify the basic view direction');
    filename = input('- insert the name of the output: ', 's');
    
    prev_matrix = get( hg, 'matrix' );
    for alpha = linspace(0,2*pi,100);
        rot_matrix = makehgtform('yrotate',alpha);
        set(hg,'Matrix',rot_matrix*prev_matrix);
        drawnow;
        gif_add_frame(gca,filename,2);
    end
end

function h_samples = draw_samples()
    h_samples = plot3( P.points(:,1), P.points(:,2), P.points(:,3), '.b', 'parent', hg, 'MarkerSize', .25 );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ARCBALL STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% register callbacks ad 
set( gcf, 'WindowButtonMotionFcn', @motion_callback );
set( gcf, 'WindowButtonDownFcn', @buttondown_callback );
set( gcf, 'WindowButtonUpFcn', @buttonup_callback );
mousestatus = 'buttonup';
START = [0,0,0];
M_previous = get( hg, 'Matrix' );

%%%%%%%%%%%% CALLBACKS %%%%%%%%%%%%%%
% motion callback, event "every" mouse movement
function motion_callback(src, event)
    % retrieve the current point
    currp = get(gcf,'CurrentPoint');
    % retrieve window geometry
    HGEOM = get( src, 'Position');
    % transform in sphere coordinates (3,4) = (WIDTH, HEIGHT)
    currp = point_on_sphere( currp, HGEOM(3), HGEOM(4) );
    
    % workaround condition (see point_on_sphere)
    if isnan(currp)
       return; 
    end
    
    %%%%% ARCBALL COMPUTATION %%%%%
    if strcmp(mousestatus, 'buttondown')
        % compute angle and rotation axis
        rot_dir = cross( START, currp ); rot_dir = rot_dir / norm( rot_dir );
        rot_ang = acos( dot( currp, START ) );
       
        % convert direction in model coordinate system
        M_tr = inv( M_previous );
        rot_dir = M_tr*[rot_dir,0]';
        rot_dir = rot_dir(1:3);
        rot_dir = rot_dir / norm( rot_dir ); % renormalize
        % construct matrix
        R_matrix = makehgtform('axisrotate',rot_dir,rot_ang);
        % set hgt matrix
        set(hg,'Matrix',M_previous*R_matrix);
        % refresh drawing
        drawnow;
    end
end

% only 1 event on click
function buttondown_callback( src, evnt )
    % change status
    mousestatus = 'buttondown';
    % retrieve the current point
    currp = get(gcf,'CurrentPoint');
    % retrieve window geometry
    HGEOM = get( src, 'Position');
    % SET START POSITION
    START = point_on_sphere( currp, HGEOM(3), HGEOM(4) );
    % SET START MATRIX
    M_previous = get( hg, 'Matrix' );    
end
function buttonup_callback( src, evnt )
    % change status
    mousestatus = 'buttonup';
    % reset the start position
    START = [0,0,0];
end

%%%%%%%%%%%% UTILITY FUNCTION %%%%%%%%%%%%%
function currp = point_on_sphere( currp, width, height )
    currp(3) = 0;
    
    % determine radius of the sphere
    R = min(width, height)/2;
    
    % TRANSFORM the point in window coordinate into 
    % the coordinate of a sphere centered in middle window
    ORIGIN = [width/2, height/2, 0];
    currp = currp - ORIGIN;
    
    % normalize position to [-1:1] WRT unit sphere 
    % centered at the origin of the window
    currp = currp / R;
    
    % if position is out of sphere, normalize it to  
    % unit length
    L = sqrt( currp*currp' );
    if L > 1
       % currp = nan; % workaround to stop evaluation
       % disp('out of sphere');
       
       currp = currp / L; 
       currp(3) = 0;
    else
       % add the Z coordinate to the scheme
       currp(3) = sqrt( 1 - currp(1)^2 - currp(2)^2 );
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END ARCBALL STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % END OF CURVELAB
