function [ centers ] = convert_edges_2_centers( edges  )
%Converts bin-edges to bin-centers 
centers = conv(edges, [0.5 0.5], 'valid');
end

