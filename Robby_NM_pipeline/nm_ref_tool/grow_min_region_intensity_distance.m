function mask = grow_min_region_intensity_distance(data, idx, max_size, alpha)
    % GROW_MIN_REGION grows a 3D region from a seed index, balancing intensity and centroid distance.
    % INPUTS:
    %   data     - 3D matrix of values (e.g., intensities)
    %   idx      - starting voxel index [x, y, z]
    %   max_size - maximum number of voxels in region
    %   alpha    - weight for intensity (0 to 1); distance weight is (1 - alpha)

    sz = size(data);
    mask = false(sz);
    mask_idx = sub2ind(sz, idx(1), idx(2), idx(3));
    mask(mask_idx) = true;
    current_points = idx;

    [dx, dy, dz] = ndgrid(-1:1, -1:1, -1:1);
    neighbors = [dx(:), dy(:), dz(:)];
    neighbors(all(neighbors == 0, 2), :) = [];

    while nnz(mask) < max_size
        disp(nnz(mask));
        % Compute centroid of current region
        [x_list, y_list, z_list] = ind2sub(sz, find(mask));
        centroid = [mean(x_list), mean(y_list), mean(z_list)];

        candidate_indices = [];
        candidate_values = [];
        candidate_distances = [];

        for i = 1:size(current_points, 1)
            p = current_points(i, :);
            for j = 1:size(neighbors, 1)
                n = p + neighbors(j, :);

                if all(n >= 1) && all(n <= sz)
                    lin_idx = sub2ind(sz, n(1), n(2), n(3));
                    if ~mask(lin_idx)
                        candidate_indices(end+1, :) = n; %#ok<AGROW>
                        val = data(n(1), n(2), n(3));
                        dist = norm(n - centroid);
                        candidate_values(end+1) = val; %#ok<AGROW>
                        candidate_distances(end+1) = dist; %#ok<AGROW>
                    end
                end
            end
        end

        if isempty(candidate_indices)
            break;
        end

        % Normalize values and distances
        norm_vals = (candidate_values - min(candidate_values)) / (max(candidate_values) - min(candidate_values) + eps);
        norm_dists = (candidate_distances - min(candidate_distances)) / (max(candidate_distances) - min(candidate_distances) + eps);

        % Combined score
        scores = alpha * norm_vals + (1 - alpha) * norm_dists;

        % Choose the candidate with the lowest score
        [~, min_idx] = min(scores);
        next_point = candidate_indices(min_idx, :);
        lin_idx = sub2ind(sz, next_point(1), next_point(2), next_point(3));

        mask(lin_idx) = true;
        current_points(end+1, :) = next_point;
    end
end
