function [CB, read_CB, UMI, read_UMI, SEQ, read_SEQ, QC, Nreads] = parse_preprocessed_h5(h5_file, cfg)
    % Parse preprocessed FASTQ data from HDF5 file
    %
    % This static method loads HDF5 files created by fastq_memory_preprocessor.py
    % and returns data in the same format as other parse_*_fastq.m functions.
    %
    % Input:
    %   h5_file: Path to HDF5 file (string)
    %   cfg: Configuration structure (unused for HDF5, included for compatibility)
    %
    % Output:
    %   CB:        Cell array of unique cell barcodes (sorted by frequency desc)
    %   read_CB:   Index array mapping reads to CB
    %   UMI:       Cell array of unique UMI sequences (sorted by frequency desc)
    %   read_UMI:  Index array mapping reads to UMI
    %   SEQ:       Cell array of unique read sequences (sorted by frequency desc)
    %   read_SEQ:  Index array mapping reads to SEQ
    %   QC:        Cell array of quality scores (one per sequence in SEQ)
    %   Nreads:    Total number of reads
    
    fprintf('Parsing preprocessed FASTQ from HDF5: %s\n', h5_file);
    
    % Check file exists
    assert(isfile(h5_file), sprintf('HDF5 file not found: %s', h5_file));
    
    % Read arrays from HDF5 (stored as uint8)
    CB_raw = h5read(h5_file, '/cell_barcodes');      % [N_CB x CB_length] uint8
    UMI_raw = h5read(h5_file, '/umis');              % [N_UMI x UMI_length] uint8
    SEQ_raw = h5read(h5_file, '/sequences');         % [N_SEQ x max_seq_len] uint8
    QC_raw = h5read(h5_file, '/quality_scores');     % [N_SEQ x max_seq_len] uint8

    % MATLAB may load HDF5 datasets with transposed dimensions; fix if needed
    if size(CB_raw, 1) < size(CB_raw, 2)
        CB_raw = CB_raw';
    end
    if size(UMI_raw, 1) < size(UMI_raw, 2)
        UMI_raw = UMI_raw';
    end
    if size(SEQ_raw, 1) < size(SEQ_raw, 2)
        SEQ_raw = SEQ_raw';
    end
    if size(QC_raw, 1) < size(QC_raw, 2)
        QC_raw = QC_raw';
    end
    
    % Read index arrays (uint32)
    read_CB = h5read(h5_file, '/read_to_cb');        % [N_reads] uint32 (0-indexed)
    read_UMI = h5read(h5_file, '/read_to_umi');      % [N_reads] uint32 (0-indexed)
    read_SEQ = h5read(h5_file, '/read_to_seq');      % [N_reads] uint32 (0-indexed)
    
    % Convert from 0-based to 1-based indexing for MATLAB
    read_CB = read_CB + 1;
    read_UMI = read_UMI + 1;
    read_SEQ = read_SEQ + 1;
    
    % Get total number of reads
    Nreads = length(read_CB);
    
    % Convert uint8 arrays back to strings (char arrays)
    % Remove trailing null bytes (uint8(0))
    
    % Cell barcodes
    CB = cell(size(CB_raw, 1), 1);
    for i = 1:size(CB_raw, 1)
        cb_bytes = CB_raw(i, :);
        % Find last non-zero byte
        last_nonzero = find(cb_bytes ~= 0, 1, 'last');
        if isempty(last_nonzero)
            CB{i} = '';
        else
            CB{i} = char(cb_bytes(1:last_nonzero));
        end
    end
    
    % UMIs
    UMI = cell(size(UMI_raw, 1), 1);
    for i = 1:size(UMI_raw, 1)
        umi_bytes = UMI_raw(i, :);
        last_nonzero = find(umi_bytes ~= 0, 1, 'last');
        if isempty(last_nonzero)
            UMI{i} = '';
        else
            UMI{i} = char(umi_bytes(1:last_nonzero));
        end
    end
    
    % Sequences
    SEQ = cell(size(SEQ_raw, 1), 1);
    for i = 1:size(SEQ_raw, 1)
        seq_bytes = SEQ_raw(i, :);
        last_nonzero = find(seq_bytes ~= 0, 1, 'last');
        if isempty(last_nonzero)
            SEQ{i} = '';
        else
            SEQ{i} = char(seq_bytes(1:last_nonzero));
        end
    end
    
    % Quality scores
    if size(QC_raw, 1) == length(read_CB)
        % QC already per-read
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
        % QC per-sequence; expand to per-read
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
        QC = QC(read_SEQ);
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
        fprintf('  Unique CBs: %d\n', metadata.unique_cbs);
        fprintf('  Unique UMIs: %d\n', metadata.unique_umis);
        fprintf('  Unique sequences: %d\n', metadata.unique_seqs);
        fprintf('  Platform: %s\n', metadata.platform);
        
    catch
        % Metadata not available, continue silently
        fprintf('  Loaded %d reads with %d unique CBs, %d UMIs, %d sequences\n', ...
            Nreads, length(CB), length(UMI), length(SEQ));
    end
    
end
