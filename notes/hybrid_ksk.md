# Hybrid key switching

Let $p1 , p2\in R_{Q}$. We need to multiply $p1$ with encryption (or if you prefer, masked version) of $p2$. Let $ct = p2 + as + e$. If we multiply $p1$ with $ct$, we will get:
$$p1p2 + p1\cdot as + p1 \cdot e$$
But now the noise has blown up to $p1 \cdot e$, of which each coefficient is as big as $Q$.

So how can we eliminate noise blow up? One way is to decompose $p1$ into smaller values and somehow recompose it while multiplying it with (potentially multiple) ciphertexts of $p2$.

Since decomposition/recomposition of polynomial in $R_{Q}$ equals decomposition and recomposition of each coefficient of the polynomial, we abuse the notation to denote a coefficient of polynomial $p \in R_{Q}$ with $p$.

Since coefficient $p1 \in Q$ is in RNS basis, naturally we can decompose $p1$ into smaller moduli $Q_j$s such that:
$$\prod Q_j = Q$$
where $Q_j$ consist of at-most dnum primes.

For example, if $Q = \prod q_0 \times \cdots \times q_4$ and dnum = 3, then $Q_0 = q_0 \times q_1 \times q_3$ and $Q_1 = q_4$

Let
$$decomp(p1) = [[p1]_{Q_0}, ..., [p1]_{Q_{\alpha-1}}]$$
and the corresponding recomposition factors be
$$[[p2[(\frac{Q}{Q_0})^{-1}]_{Q_0}\frac{Q}{Q_0}]_Q, ..., p2[(\frac{Q}{Q_{alpha-1}})^{-1}]_{Q_{alpha-1}}\frac{Q}{Q_{alpha-1}}]_Q]$$
Then we can recompose as:
$$[p1p2]_Q = \sum_{i=0}^{alpha-1} [p1]_{Ql_i} \cdot [p2(\frac{Ql}{Ql_i})^{-1}]_{Ql_i}\frac{Ql}{Ql_i} \mod{Q}$$
If we were to encrypt each element of recomposition factor as a ciphertext, then the noise of the resulting ciphertext encrypting $p1p2$ after summation will be $max(\log{Q_j})$.

To make the noise growth of key switching negligible, we will have need some way to eliminate $max(\log{Q_j})$. One way to do it is by dividing the ciphertext by a value that roughly has $max(\log{Q_j})$ bits.

Define $P_s$ (called specialP) as,
$$P_s = \prod p_j$$
such that $\log{P_s} \approx \log{Q_j}$

Again let
$$decomp(p1) = [[p1]_{Q_0}, ..., [p1]_{Q_{\alpha-1}}]$$
and the corresponding recomposition factors, this time, be
$$[[P_s \cdot p2[(\frac{Q}{Q_0})^{-1}]_{Q_0}\frac{Q}{Q_0}]_Q, ..., P_s\cdot p2[(\frac{Q}{Q_{alpha-1}})^{-1}]_{Q_{alpha-1}}\frac{Q}{Q_{alpha-1}}]_Q]$$

Then we can recompose as:
$$[P_sp1p2]_{QP_s} = \sum_{i=0}^{alpha-1} [p1]_{Q_i} \cdot (P_sp2[(\frac{Ql}{Ql_0})^{-1}]_{Ql_i}\frac{Ql}{Ql_i}) \mod QP_s$$
Again, if recomposition factors were encrypted, then the error of output ciphertext will roughly equal $max(\log{Q_j})$. But now we scale the output polynomials by $1/P_s$ and switch the modulus to $Q$. If $\log{P_s} \approx max(\log{Q_j})$, the output ciphertext has negligible error.

Note that error in output ciphertext depends on:

1. Difference between $\log{P_s}$, $max(\log{Q_j})$
2. and error in recomposition factors ciphertexts

### Parameters

The main parameter is dnum (called decomposition number).

Let $Q$ be product of k primes $\prod_{i=0} q_i$. $Q$ is decomposed into $\alpha$ smaller modulis $Q_j$s each consisting of disjoint subset of primes in $Q$ of count at-most dnum. Hence, $\alpha = ceil(k/dnum)$. For convenience we set $Q_j = q_i \times ... \times q_{min(i+dnum-1, k-1)}$.

$P_s$ is usually chosen such that $\log{P_s} \approx max(\log{Q_j})$.

There's a trade off between security and key switching time complexity when choosing dnum:

1. High dnum, and hence lower NTTs in key switching, requires $P_s$ to be high (to reduce key switching noise) as well. Therefore, to maintain same security, Q, i.e. depth of the circuit, must be reduced.
2. On other hand low dnum increases number of NTTs in key switching, hence time complexity.

### Key generation

For parameters $dnum, Q, P_s$ such that
$$Q = \prod^{k-1}_{i=0} q_i$$
$$P_s = \prod{p_j}$$
find $\alpha$ smaller moduli chains as $Q_0, ..., Q_{\alpha-1}$ where $\alpha = ceil(L/dnum)$
$$Q_i = q_i \times ... \times q_{min((i + dnum - 1), k-1)}$$

Precompute recomposition factor corresponding factor $\gamma_j$ corresponding to each $Q_j$ as:
$$\gamma_j = P_s \cdot [(\frac{Q}{Q_j})^{-1}]_{Q_j} \cdot \frac{Q}{Q_j})^{-1}]_{QP_s}$$

To generate key switching key for input polynomial $p \in R_{Q}$,

1. Extend $p \in R_{Q}$ to $R_{QPs}$.
2. For each $gamma_j$ encrypt its product with $p$ and set it is key switching key ciphertext $j$ as:
   $$ksk_j = (a_js + e_j + \gamma_j p, -a_j) \in R^2_{QP_s}$$
   Output key switching key as collection of ciphertexts $[ksk_0, ..., ksk_{alpha-1}] \in R^{\alpha}_{QP_s}$

Implementation notes:

1. $\gamma_i \mod P_s = 0$. Hence when calculating $\gamma_jp$ in RNS basis, for any $p_j \in P_s \prod{p_j}$, $\gamma_jp \mod{p_j} = 0$. We can skip calculating (only!) the product for part $P_s$. Therefore, we don't really need to find representation of $p$, that is extend $p$, in $P_s$.
2. We sample $-a$ pseudorandomly and convert it to $a$ in key generation. This approach much cleaner if ciphertexts are seeded since pseudo random part of ciphertext, $-a$, can be sampled directly without any additional operations.

## Key switch

To key switch polynomial $p1 \in R_{Q}$,

1. Decompose $p1$ into $R_{Q_j}$s. This requires finding the representation of each coefficient in $R_Q$ in $R_{Q_j}$, hence requires no additional operation.
2. For each $Q_j$, extend $[p1]_{Q_j} \in R_{Q_j}$ to $R_{QP_s}.$ This is done by switching CRT basis of each coefficient from $Q_j$ to $QP_s$ (using SwitchCRTBasis, but since $Q_j$ is part of $Q$ we can skip calculating representation for primes that are in $Q_j$).

    **Switch $[p1]_{Q_j} \in Q_j$ to $QP_s$:**
    We will use the same algorithm as SwitchCRTBasis but with a few modifications. Observe that $Q_j$ consists of some primes in $Q$. For example, if $Q_j = q_0 \times q_1 \times q_2$, then all $q_0, q_1, q_2$ exist in $QP_s$. In SwichCRTBasis we come across 2 cases.

    1. Case 1: when $q_k$ in $QP_s$ is not in $Q_j$. In this case, we need to find the representation of $x \in Q_j$ for $q_k$.
    2. Case 2: when $q_k$ in $QP_s$ is in $Q_j$ and $q_k == q_{j,i}$. In this case the representation of $x$ remains unchanged and equals $x_{j,i}$. This is because, in summation for $q_k$, the case when $q_k \neq q_{j,i}$, $q_k$ is in $Q_j/q_{j,i}$ and hence the value vanishes (i.e. becomes 0). And for the case when $q_k = q_{j,i}$, $(Q_j/q_{j,i})^{-1} \times Q_j/q_{j,i}$ equals 1. Hence, the summation outputs $x_{j,i}$.

3. Denote representation of $[p1]_{Q_j}$ in $QP_s$ as $p1_j$. Calculate inner product:
   $$(c_0, c_1) = \sum_{j=0}^{\alpha-1} p1_j (ksk_j)$$
4. Output ciphertext $ct' = (c_0, c_1)$ $\in R^2_{QP_s}$. $ct'$ is such that:
   $$c_0 + c_1s = P_sp1p2 + E$$
   For some error $E$.

    We now scale $ct'$ by $1/P_s$ and switch the modulus to $Q$:
    $$ct = [\lceil \frac{1}{P}ct' \rfloor]_Q$$

    Now $ct \in R_Q$ is such that
    $$c_0 + c_1s = p1p2 + e$$
    where e is small.

Implementation notes:

1. Key switching key ciphertexts are used only in inner sum. Hence, they should stay in Evaluation representation to reduce NTTs at key switching runtime.
2. To calculate $[\lceil \frac{1}{P_s}c_0 \rfloor]_Q$ in RNS basis where $c_0 \in R_{QP_s}$ we use algorithm 2 of 2018/931 (we refer to it as ApproximateModDown)
