function [CB, read_CB, UMI, read_UMI, SEQ, read_SEQ, QC, Nreads] = parse_splitpipe_fastq(fastq_file, cfg)
% PARSE_SPLITPIPE_FASTQ Parse output from split-pipe --mode comb
%
% This function handles Parse Biosciences data that has been preprocessed
% with split-pipe in "comb" mode. The split-pipe output format:
%   - Single FASTQ file with cell barcode + UMI in read header
%   - Read header format: @readname_CB:XXXXXXXX_UMI:YYYYYYYY
%   - Sequence contains the actual insert (CARLIN amplicon)
%
% Input:
%   fastq_file - Cell array with path to split-pipe combined FASTQ file(s)
%                Can be Nx1 (single files) or Nx2 where column 2 is ignored
%   cfg        - Configuration struct with CB and UMI length parameters
%
% Output:
%   CB, read_CB   - Unique cell barcodes and their read indices
%   UMI, read_UMI - Unique UMIs and their read indices  
%   SEQ, read_SEQ - Unique sequences and their read indices
%   QC            - Quality scores for each read
%   Nreads        - Total number of reads processed

    % Handle input format - split-pipe outputs single combined file
    if size(fastq_file, 2) == 2
        % If given as pairs, use only first column (combined file)
        fastq_file = fastq_file(:, 1);
    end
    
    Nfastqs = size(fastq_file, 1);
    Nreads = cell(Nfastqs, 1);
    
    CB = cell(Nfastqs, 1);
    UMI = cell(Nfastqs, 1);
    QC = cell(Nfastqs, 1);
    SEQ = cell(Nfastqs, 1);
    
    % Get expected lengths from config
    cb_length = cfg.CB{1}.length;
    umi_length = cfg.UMI.length;
    
    for i = 1:Nfastqs
        
        assert(exist(fastq_file{i}, 'file') == 2, ...
            sprintf('Missing FASTQ file: %s', fastq_file{i}));
        
        fprintf('Parsing split-pipe FASTQ: %s\n', fastq_file{i});
        
        [fastq_file(i), ext] = FastQData.maybe_unzip(fastq_file(i));
        
        info = fastqinfo(fastq_file{i});
        Nreads{i} = info.NumberOfEntries;
        
        fprintf('Reads in FASTQ: %d\n', Nreads{i});
        
        % Pre-allocate
        CB{i} = cell(Nreads{i}, 1);
        UMI{i} = cell(Nreads{i}, 1);
        QC{i} = cell(Nreads{i}, 1);
        SEQ{i} = cell(Nreads{i}, 1);
        
        % Read FASTQ
        [headers, sequences, qualities] = fastqread(fastq_file{i});
        
        % Parse cell barcodes and UMIs from headers
        [CB{i}, UMI{i}] = parse_splitpipe_headers(headers, cb_length, umi_length);
        
        SEQ{i} = sequences(:);
        QC{i} = qualities(:);
        
        FastQData.maybe_clear_unzipped(fastq_file(i), ext);
        
    end
    
    Nreads = cell2mat(Nreads);
    Nreads = sum(Nreads);
    
    if (Nfastqs > 1)
        fprintf('Total reads in FASTQs: %d\n', Nreads);
    end
    
    QC = vertcat(QC{:});
    CB = vertcat(CB{:});
    UMI = vertcat(UMI{:});
    SEQ = vertcat(SEQ{:});
    
    [CB, ~, read_CB] = unique_by_freq(CB);
    [UMI, ~, read_UMI] = unique_by_freq(UMI);
    [SEQ, ~, read_SEQ] = unique_by_freq(SEQ);

end


function [CB, UMI] = parse_splitpipe_headers(headers, cb_length, umi_length)
% PARSE_SPLITPIPE_HEADERS Extract CB and UMI from split-pipe header format
%
% split-pipe --mode comb header formats supported:
%   Format 1: @readname_CB:XXXXXXXXXXXXXXXXXXXXXXXX_UMI:YYYYYYYYYY
%   Format 2: @readname CB:XXXXXXXXXXXXXXXXXXXXXXXX UMI:YYYYYYYYYY
%   Format 3: @CB:XXXXXXXXXXXXXXXXXXXXXXXX_UMI:YYYYYYYYYY_readname
%
% The function auto-detects the format from the first header.

    n_reads = length(headers);
    CB = cell(n_reads, 1);
    UMI = cell(n_reads, 1);
    
    % Detect format from first header
    first_header = headers{1};
    
    % Try different regex patterns for split-pipe output formats
    patterns = {
        'CB:([ACGTN]+).*UMI:([ACGTN]+)',     % Standard format
        'CB[_:]([ACGTN]+).*UMI[_:]([ACGTN]+)', % Variant with underscore
        '_([ACGTN]{%d})_([ACGTN]{%d})$'       % Suffix format (CB_UMI at end)
    };
    
    % Replace length placeholders in pattern 3
    patterns{3} = sprintf(patterns{3}, cb_length, umi_length);
    
    matched_pattern = 0;
    for p = 1:length(patterns)
        tokens = regexp(first_header, patterns{p}, 'tokens', 'once');
        if ~isempty(tokens)
            matched_pattern = p;
            break;
        end
    end
    
    if matched_pattern == 0
        error(['Cannot parse split-pipe header format. First header: %s\n' ...
               'Expected formats:\n' ...
               '  @readname_CB:XXXX_UMI:YYYY\n' ...
               '  @readname CB:XXXX UMI:YYYY\n' ...
               '  @readname_XXXX_YYYY (CB_UMI suffix)'], first_header);
    end
    
    fprintf('Detected split-pipe header format: pattern %d\n', matched_pattern);
    
    % Parse all headers
    pattern = patterns{matched_pattern};
    
    % Vectorized regexp on entire cell array
    all_tokens = regexp(headers, pattern, 'tokens', 'once');
    valid = ~cellfun(@isempty, all_tokens);

    % Extract CB and UMI from valid tokens
    CB(valid) = cellfun(@(t) t{1}, all_tokens(valid), 'un', false);
    UMI(valid) = cellfun(@(t) t{2}, all_tokens(valid), 'un', false);

    % Fill invalid entries with N-pads
    CB(~valid) = {repmat('N', 1, cb_length)};
    UMI(~valid) = {repmat('N', 1, umi_length)};

    % Validate lengths - replace wrong-length entries
    bad_cb = cellfun(@length, CB) ~= cb_length;
    bad_umi = cellfun(@length, UMI) ~= umi_length;
    CB(bad_cb) = {repmat('N', 1, cb_length)};
    UMI(bad_umi) = {repmat('N', 1, umi_length)};

end
