function [alleles, which_seqs, weight_contribution] = call_alleles_coarse_grain(aligned_seqs, aligned_seq_weights, dominant_only)

    if (nargin < 2)
        aligned_seq_weights = ones(size(aligned_seqs));
    else        
        assert(size(aligned_seqs,1) == size(aligned_seq_weights,1));
    end
    
    valid_mask = find(~cellfun(@isempty, aligned_seqs));
    aligned_seqs = aligned_seqs(valid_mask);
    aligned_seq_weights = aligned_seq_weights(valid_mask);    
    
    alleles = [];
    which_seqs = [];
    weight_contribution = [];
    
    if (isempty(aligned_seqs))
        return;
    end

    if (nargin < 3)
        dominant_only = true;
    end

    events = cellfun(@(x) x.get_event_structure(), aligned_seqs, 'un', false);
    events = vertcat(events{:});
    [events, ~, which_event] = unique_by_freq(events, aligned_seq_weights);            
    event_weight = accumarray(which_event, aligned_seq_weights);
    assert(issorted(event_weight, 'descend'));
    
    dominant_frac=0.5;
    if (dominant_only)
        if (event_weight(1) / sum(event_weight) > dominant_frac)
            N = 1;
        else
            return;
        end
    else
        N = size(events,1);
    end
    
    alleles = cell(N,1);
    which_seqs = cell(N,1);
    weight_contribution = cell(N,1);
        
    for i = 1:N

        event_mask = find(which_event == i);
        seqs_for_event = cellfun(@(x) x.get_seq(), aligned_seqs(event_mask), 'un', false);
        seq_weights_for_event = aligned_seq_weights(event_mask);
        seq_lengths = cellfun(@length, seqs_for_event);

        % Empirically we find almost all reads grouped by UMI/CB are of the same length
        % so insertions/deletions due to sequencing and RT are negligible. Although we 
        % could do a more sophisticated length-dependent multialign, in practice it's
        % not worth the trouble just to include a few more reads.
        
        if (length(unique(seq_lengths)) == 1)
            consensus = weighted_consensus(seqs_for_event, seq_weights_for_event);
            [~, ref_ind] = max(seq_weights_for_event);
            alleles{i} = CARLIN_def.desemble_sequence(consensus, aligned_seqs{event_mask(ref_ind)}.get_ref());
            which_seqs{i} = valid_mask(event_mask);
            weight_contribution{i} = seq_weights_for_event;
        else
            mode_length = weighted_mode_scalar(seq_lengths, seq_weights_for_event);
            length_mask = find(seq_lengths == mode_length);
            if (dominant_only && (N==1) && (sum(seq_weights_for_event(length_mask)) / sum(event_weight) <= dominant_frac))
                break;
            end
            consensus = weighted_consensus(seqs_for_event(length_mask), seq_weights_for_event(length_mask));
            [~, ref_ind] = max(seq_weights_for_event(length_mask));
            ref_ind = length_mask(ref_ind);
            alleles{i} = CARLIN_def.desemble_sequence(consensus, aligned_seqs{event_mask(ref_ind)}.get_ref());
            which_seqs{i} = valid_mask(event_mask(length_mask));
            weight_contribution{i} = seq_weights_for_event(length_mask);
        end
    end
    
    if (N == 1)
        alleles = alleles{1};
        which_seqs = which_seqs{1};
        weight_contribution = weight_contribution{1};
    end
end

function consensus = weighted_consensus(seqs, weights)
    % Compute per-column consensus from weighted sequences without expansion.
    % Instead of repelem+vertcat+mode (O(total_weight * seq_len)),
    % uses weighted character frequencies (O(num_unique * seq_len)).
    seq_mat = vertcat(seqs{:});
    seq_len = size(seq_mat, 2);
    consensus = blanks(seq_len);
    for col = 1:seq_len
        chars = seq_mat(:, col);
        [uc, ~, ic] = unique(chars);
        char_weights = accumarray(ic, weights);
        [~, max_idx] = max(char_weights);
        consensus(col) = uc(max_idx);
    end
end

function m = weighted_mode_scalar(values, weights)
    % Compute weighted mode of a numeric vector without expansion.
    [uv, ~, ic] = unique(values);
    val_weights = accumarray(ic, weights);
    [~, max_idx] = max(val_weights);
    m = uv(max_idx);
end