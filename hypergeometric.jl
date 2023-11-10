
function hypergeometric(k, N, K, n)
    (binomial(K, k) * binomial(N-K, n-k)) / binomial(N, n)
end


: