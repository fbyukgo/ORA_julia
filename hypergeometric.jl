using BigCombinatorics

function hypergeometric(k, N, K, n)
    pval = (Binomial(K, k) * Binomial(N-K, n-k)) / Binomial(N, n)
    return pval
end


