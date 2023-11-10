using BigCombinatorics

function hypergeometric(k, N, K, n)
    (Binomial(K, k) * Binomial(N-K, n-k)) / Binomial(N, n)
end

@time hypergeometric(17, 1000, 500, 20)
