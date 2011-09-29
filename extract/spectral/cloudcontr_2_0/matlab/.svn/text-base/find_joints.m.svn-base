function [joints, segments] = find_joints(pts, spls, corresp, A, show_joints)
joints = zeros(0,2);% [index, diameter]
segments = zeros(size(spls,1),1);%- for joints, NaN for branches, 0 for NaN.
segments(:) = NaN;
jidx = 0;
for i=1:length(segments)
    if ~isnan(segments(i)), continue, end;    
    
    if isnan(spls(i,1))
        segments(i) = 0;
        continue;
    end
        
    links = find( A(i,:)==1 );
    if length(links) > 3 % joint
        jidx = jidx - 1;
        segments(i) = jidx;
        [tmp, diam] = GS.compute_bbox(pts(corresp==i,:));
        joints(end+1,:) = [i,diam];
    end
end

if show_joints
    figure('Name','Find joints','NumberTitle','off');set(gcf,'color','white');movegui('southwest');
    plot_skeleton(spls, A);hold on;
    scatter3(spls(joints(:,1),1),spls(joints(:,1),2),spls(joints(:,1),3),200,'b', 'filled');
    axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(-90,0);view3d rot;
end