function [SEQ, read_SEQ, QC, Nreads] = parse_preprocessed_h5(h5_file, cfg)
    % Parse preprocessed FASTQ data from HDF5 file (bulk sequencing)
    %
    % This static method loads HDF5 files created by fastq_memory_preprocessor.py
    % for bulk sequencing data (no cell barcodes).
    %
    % Input:
    %   h5_file: Path to HDF5 file (string)
    %   cfg: Configuration structure (unused for HDF5, included for compatibility)
    %
    % Output:
    %   SEQ:      Cell array of unique read sequences (sorted by frequency desc)
    %   read_SEQ: Index array mapping reads to SEQ
    %   QC:       Cell array of quality scores (one per read)
    %   Nreads:   Total number of reads

    fprintf('Parsing preprocessed bulk FASTQ from HDF5: %s\n', h5_file);

    % Check file exists
    assert(isfile(h5_file), sprintf('HDF5 file not found: %s', h5_file));

    % Read arrays from HDF5 (stored as uint8)
    SEQ_raw = h5read(h5_file, '/sequences');         % [N_SEQ x max_seq_len] uint8
    QC_raw = h5read(h5_file, '/quality_scores');     % [N_reads or N_SEQ x max_len] uint8

    % MATLAB may load HDF5 datasets with transposed dimensions; fix if needed
    if size(SEQ_raw, 1) < size(SEQ_raw, 2)
        SEQ_raw = SEQ_raw';
    end
    if size(QC_raw, 1) < size(QC_raw, 2)
        QC_raw = QC_raw';
    end

    % Read index array (uint32, 0-indexed from Python)
    read_SEQ = h5read(h5_file, '/read_to_seq');      % [N_reads] uint32

    % Convert from 0-based to 1-based indexing for MATLAB
    read_SEQ = read_SEQ + 1;

    % Get total number of reads
    Nreads = length(read_SEQ);

    % Convert uint8 arrays back to strings (char arrays)
    % Remove trailing null bytes (uint8(0))

    % Sequences
    SEQ = cell(size(SEQ_raw, 1), 1);
    for i = 1:size(SEQ_raw, 1)
        seq_bytes = SEQ_raw(i, :);
        % Find last non-zero byte
        last_nonzero = find(seq_bytes ~= 0, 1, 'last');
        if isempty(last_nonzero)
            SEQ{i} = '';
        else
            SEQ{i} = char(seq_bytes(1:last_nonzero));
        end
    end

    % Quality scores - check if per-read or per-sequence
    if size(QC_raw, 1) == Nreads
        % QC is per-read (matches number of reads)
        fprintf('  Using per-read quality scores\n');
        QC = cell(size(QC_raw, 1), 1);
        for i = 1:size(QC_raw, 1)
            qc_bytes = QC_raw(i, :);
            last_nonzero = find(qc_bytes ~= 0, 1, 'last');
            if isempty(last_nonzero)
                QC{i} = '';
            else
                QC{i} = char(qc_bytes(1:last_nonzero));
            end
        end
    else
        % QC is per-sequence; expand to per-read using read_SEQ indices
        fprintf('  Expanding per-sequence quality scores to per-read\n');
        QC_per_seq = cell(size(QC_raw, 1), 1);
        for i = 1:size(QC_raw, 1)
            qc_bytes = QC_raw(i, :);
            last_nonzero = find(qc_bytes ~= 0, 1, 'last');
            if isempty(last_nonzero)
                QC_per_seq{i} = '';
            else
                QC_per_seq{i} = char(qc_bytes(1:last_nonzero));
            end
        end
        % Map to per-read using indices
        QC = QC_per_seq(read_SEQ);
    end

    % Log metadata from HDF5 attributes if available
    try
        h5_info = h5info(h5_file);

        % Extract attributes
        attrs = h5_info.Attributes;
        metadata = struct();
        for i = 1:length(attrs)
            attr_name = attrs(i).Name;
            attr_value = h5readatt(h5_file, '/', attr_name);
            metadata.(attr_name) = attr_value;
        end

        % Print statistics
        fprintf('  Total reads: %d\n', metadata.total_reads);
        fprintf('  Valid reads (indexed): %d\n', Nreads);
        fprintf('  Filtered reads: %d\n', metadata.filtered_reads);
        fprintf('  Unique sequences: %d\n', metadata.unique_seqs);
        if isfield(metadata, 'platform')
            fprintf('  Platform: %s\n', metadata.platform);
        end

    catch
        % Metadata not available, continue silently
        fprintf('  Loaded %d reads with %d unique sequences\n', ...
            Nreads, length(SEQ));
    end

end
