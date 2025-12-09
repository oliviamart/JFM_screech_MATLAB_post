function [eigvals_tracked] = sort_eigenvalues(eigvals, eigvecs, freqs, U_inf)

[N, F] = size(eigvals);
mac_threshold = 0.80;
tracks = struct('freq', {}, 'eigval', {}, 'eigvec', {}, 'active', {});
tracks_new = struct('freq', {}, 'eigval', {}, 'eigvec', {}, 'active', {});
next_track_index = 1;

% get rid of eigenvalues outside range of interest
% this means if 
% (1) acoustic
% (2) very high negative wavenumber
% (3) evanescent
% (4) NaN
for i = 1:N
    for j = 1:F
        freq = -1*freqs(j)*2*pi/(1-U_inf);
        if (real(eigvals(i,j)) >= freq) | (abs(imag(eigvals(i,j))) > 10^-4 ...
                | real(eigvals(i,j)) < -10 | isnan(eigvals(i,j)))
            eigvals(i,j) = NaN;
            eigvecs{j}(:,i) = NaN;
        end
    end
end

% Initialize with first frequency
for i = 1:N
    tracks(next_track_index).freq = freqs(1);
    tracks(next_track_index).eigval = eigvals(i,1);
    tracks(next_track_index).eigvec = eigvecs{1}(:,i);
    tracks(next_track_index).active = true;
    next_track_index = next_track_index + 1;
end

% Track across frequencies
for f = 2:F
    curr_vals = eigvals(:,f);
    curr_vecs = eigvecs{f};
    assigned = false(N,1);

    for t = 1:length(tracks)
        if isnan(tracks(t).eigval(end))
            tracks(t).active = false;
            continue;
        end

        if ~tracks(t).active
           continue;
        end

        best_mac = 0;
        best_mac_next = 0;
        best_j = 0;

        for j = 1:N
            if assigned(j), continue; end
            %if isnan(eigvals(j,f)), continue; end

            v1 = tracks(t).eigvec;
            v2 = curr_vecs(:,j);
            v1 = v1 / norm(v1);
            v2 = v2 / norm(v2);
            phase = angle(v1' * v2);
            v2_aligned = v2 * exp(-1i * phase);
            mac = abs(v1' * v2_aligned)^2 / (norm(v1)^2 * norm(v2_aligned)^2);
            
            if mac > best_mac
                best_mac_next = best_mac;
                best_mac = mac;
                best_j = j;
            end
        end

        if best_mac > mac_threshold
            j = best_j;
            tracks(t).freq(end+1) = freqs(f);
            tracks(t).eigval(end+1) = curr_vals(j);
            tracks(t).eigvec = curr_vecs(:,j);
            assigned(j) = true;
        else
            tracks(t).active = false;
        end
    end

    % New tracks for unassigned eigenvectors
    for j = 1:N
        if ~assigned(j) & ~isnan(curr_vals(j))
            tracks(next_track_index).freq = freqs(f);
            tracks(next_track_index).eigval = curr_vals(j);
            tracks(next_track_index).eigvec = curr_vecs(:,j);
            tracks(next_track_index).active = true;
            next_track_index = next_track_index + 1;
        end
    end
end

tracks_new = tracks;

% Convert to padded matrix
eigvals_tracked = NaN(length(tracks_new), F);

for t = 1:length(tracks_new)
    for k = 1:length(tracks_new(t).freq)
        % Find the closest frequency index (handles floating point issues)
        [~, idx] = min(abs(freqs - tracks_new(t).freq(k)));
        eigvals_tracked(t, idx) = tracks_new(t).eigval(k);
    end
end


end

