function edge_num = plot_connectivity(pts, A, sizee,colore)
    A( logical(eye(size(A)))) = 0;
    tmax = max(max(A))*1.0;
    edge_num = 0;

    if nargin<3            
        for i=1:(size(A,1)-1)
            for j=(i+1):size(A,2)
                if( A(i,j)>0 )
                    idx = [i;j];
                    line( pts(idx,1),pts(idx,2),pts(idx,3), 'LineWidth', 150*A(i,j)/tmax, 'Color', [A(i,j)/tmax,.0,.0]);
%                       text(pts(i,1),pts(i,2),pts(i,3),int2str(i)) ;
                    edge_num = edge_num +1;
                end
            end
        end
    else
        for i=1:(size(A,1)-1)
            for j=(i+1):size(A,2)
                if( A(i,j)>0 )
                    idx = [i;j];
                    line( pts(idx,1),pts(idx,2),pts(idx,3), 'LineWidth', sizee, 'Color', colore);
%                       text(pts(i,1),pts(i,2),pts(i,3),int2str(i)) ;
                    edge_num = edge_num +1;
                end
            end
        end
    end
end