function [so1, so2, ro] = pairwise_align(a1, a2)

    assert(isa(a1, 'AlignedSEQ') && isa(a2, 'AlignedSEQ'));
    
    si1 = a1.get_seq();
    si2 = a2.get_seq();

    ri1 = a1.get_ref();
    ri2 = a2.get_ref();

    assert(isequal(ri1(ri1~='-'), ri2(ri2~='-')));

    L1 = length(si1);
    assert(L1 == length(ri1));
    
    L2 = length(si2);
    assert(L2 == length(ri2));
    
    if (isequal(ri1, ri2) && isequal(si1, si2))
        so1 = si1;
        so2 = si2;
        ro = ri1;
        return;
    end 
    
    p1 = 1;
    p2 = 1;

    % Pre-allocate output arrays to maximum possible length (L1+L2)
    max_len = L1 + L2;
    so1 = blanks(max_len);
    so2 = blanks(max_len);
    ro  = blanks(max_len);
    k = 0;  % write pointer (number of characters written so far)

    while (p1 <= L1 || p2 <= L2)
        if (p1 > L1)
            left = L2-p2+1;
            assert(all(ri2(p2:end) == '-'));
            ro(k+1:k+left)  = '-';
            so1(k+1:k+left) = '-';
            so2(k+1:k+left) = si2(p2:L2);
            k = k + left;
            break;
        end
        if (p2 > L2)
            left = L1-p1+1;
            assert(all(ri1(p1:end) == '-'));
            ro(k+1:k+left)  = '-';
            so1(k+1:k+left) = si1(p1:L1);
            so2(k+1:k+left) = '-';
            k = k + left;
            break;
        end
        k = k + 1;
        if (ri1(p1) == '-')
            if (ri2(p2) == '-')
                ro(k)  = '-';
                so1(k) = si1(p1);
                so2(k) = si2(p2);
                p1 = p1+1;
                p2 = p2+1;
            else
                ro(k)  = '-';
                so1(k) = si1(p1);
                so2(k) = '-';
                p1 = p1+1;
            end
        else
            if (ri2(p2) == '-')
                ro(k)  = '-';
                so1(k) = '-';
                so2(k) = si2(p2);
                p2 = p2+1;
            else
                assert(ri1(p1) == ri2(p2))
                ro(k)  = ri1(p1);
                so1(k) = si1(p1);
                so2(k) = si2(p2);
                p1 = p1+1;
                p2 = p2+1;
            end
        end
    end

    % Truncate to actual length
    so1 = so1(1:k);
    so2 = so2(1:k);
    ro  = ro(1:k);

    L = k;
    assert(length(so1) == L);
    assert(length(so2) == L);
    
    assert(all(ro(ro~='-') == ri1(ri1~='-')));
    assert(all(ro(ro~='-') == ri2(ri2~='-')));
    assert(all(so1(so1~='-') == si1(si1~='-')));
    assert(all(so2(so2~='-') == si2(si2~='-')));
end