# BFV

### Secret key encryption

Given encoded message $m(X)$, we first calculate $\Delta m$ as
$$\Delta m = t^{-1}[Qm]_t \mod{Q}$$
where $Q = \prod q_i$

Then we encrypt $\Delta m$ as ciphertext ct:
$$ct = (a s + e + \Delta m, -a)$$

Notes:

1. We precompute $t^{-1} \mod q_i$ for each $q_i$ in $Q$ and scale $[[Qm]_t]_{q_i}$ by $t^{-1}_{q_i}$ to output $[t^{-1}[Qm]_t]_{q_i}$

### Public key Encryption

TODO

### BFV decryption

TODO

### BFV multiplication

We implement Algorithm 2 of 2021/204. Following is a brief description of the algorithm:

1. Given $ct, ct' \in R^2_{Q}$ as input:

    1. Switch the CRT basis of $ct$ from Q to P to extend the basis of $ct$ to $R_{QP}$.
        1. We use SwitchCRTBasis function to switch basis Q to P.
        2. Basis Q and basis P collectively represent $ct \in {R_QP}$
    2. Calculate:
       $$ct'' = [\lceil \frac{P}{Q}ct'\rfloor]_{P}$$
       We use FastConvertPOverQ function to calculate $ct''$. Note that $ct''$ is in basis P.

        Then we switch the basis of $ct''$ from P to Q, using SwitchCRTBasis, to extend the basis of $ct''$ to $R_{QP}$.

        After (1) and (2) we will have:
        $$ct, ct'' \in R_{QP}$$

2. Calculate degree 2 ciphertext $ct^o = ct \otimes ct'' \in R_{QP}$ , where $\otimes$ is tensor operation.
    1. $ct^{o}$ is a degree 2 ciphertext that satisfies:
       $$ct^o = \frac{QP}{t^2} (mm') + E = c^o_0 + c^o_1s + c^o_2s^2$$
3. We scale down $ct^o$ by $t/P$ and round it as $[\lceil \frac{t}{P} \cdot ct^o \rfloor]_{Q}$ to finally output ciphertext product $\in R^2_Q$.
    1. We use ScaleAndRound to calculate $[\lceil \frac{t}{P} \cdot ct^o \rfloor]_{Q}$
    2. Note that the output ciphertext is still a degree 2 ciphertext.

### Relinerization

Relinearization reduces a degree 2 ciphertext $ct \in R^3_{Q}$ to a degree 1 ciphertext $ct_{r} \in R^2_{Q}$. It multiples $c_2 \in R_{Q}$ with encryption of $s^2$. Since multiplying a polynomial in $R_{Q}$ with a ciphertext will blow up the noise in ciphertext, we cannot directly multiply the polynomial with the ciphertext.

Recall that using hybrid key switching procedure one can multiply polynomial $p1 \in R_{Q}$ at level $l$ with encryption of another polynomial $p2 \in R_{Q}$ to output a ciphertext $ct_{o}$ s.t
$$c^o_{0} + c^o_{1}s   = p1p2$$
with negligible increase in noise.

We use hybrid key switching to multiply $c_2$ with encryption of $s^2$. For relinerization at level $l$, hybrid key switching key is generated for polynomial $s^2$ . We refer to the key as relinearization key as $rlk_{l}$. In general, $rlk_{l} = hybridKsk(s^2, l)$.

Let $ct$ be a degree 2 ciphertext, such that
$$ct = c_0 + c1s + c2s^2. $$

Given $rlk_{l}$ for level $l$ we key switch $c_2$ (i.e. multiply $c_2$ with $s^2$ homomorphically) to output:
$$c'_0, c'_1 = \text{keyswitch}(c_2)$$
such that,
$$c'_0 + c'_1s = c_2s^2$$
Output relinearized ciphertext $ct_{r}$ as:
$$ct_{r} = (c_0 + c'_0, c_1 + c'_1)$$

# Polynomial operations

Let
$$Q = q_0 \cdot q_i \cdots q_n$$
$$P = p_0 \cdot p_i \cdots p_m$$

**SwitchCRTBasis**

Given $[x]_Q = [x_0, ..., x_n]$ we need to calculate $[x_Q] \mod{P}$. To do so, For each $p_j$ we calculate
$$[[x]_Q + uQ]_{p_j} = \sum [x_i (\frac{Q}{q_i})^{-1}]_{q_i} \frac{Q}{q_i} \mod{p_j}$$
where $u < l/2$

We eliminate the overflow caused $uQ$ by estimating $u$:
$$u = \sum [x_i (\frac{Q}{q_i})^{-1}]_{q_i} \frac{1}{q_i}$$

and then subtracting $uQ \mod{p_j}$ from $[[x]_Q + uQ]_{p_j}$ to get $[[x]_Q]_{p_j}$

**Where does the equation for u comes from?**

TODO

**Why not use equation 3 of HPS?**

One may wonder why not use equation 3 of HPS. However, calculating $v'$ is not as straigthward as it is in equation 2. Although one may benefit from grouping $[\frac{Q}{q_i}^{-1} \cdot \frac{Q}{q_i}]_{p_j}$, but to estimate $v'$ we need to calculate $\frac{Q}{q_i} \cdot x_i \cdot \frac{1}{q_i}$ without any modular reduction. This means it requires, recalling $\frac{Q}{q_i}^{-1}$ is a big integer, big integer arithmetic. Hence, we prefer equation 2 over equation 3 because in equation 2 estimating $v$ only requires simple 64 bit arithmetic.

**SimpleScaleAndRound**

Simple scale and round is used in BFV decryption

Given a polynomial $x \in R_Q$ it outputs
$$[\lceil\frac{t}{Q}x\rfloor]_t$$
We implement simple scaling procedure using digit decomposition technique outlined in 2021/204.

Digit decomposition is helpful to contain error induced by $x_i \cdot (\theta_i + \epsilon_i)$ where $x_i$ is a value in $q_i$ and $\epsilon_i$ is error induced due to precision loss when calculating $\theta_i$ (note that $\theta_i$ is a fraction).

Value of $x_i \epsilon_i$ increases with value of $x_i$, which is in range $(0, q_i]$. In general, when using double floats to store $\theta_i$ and $q_i >= 51$, the error induced gets big enough to impact the resulting polynomial $\in R_t$.

The idea of digit decomposition is to decompose $x_i$ into $x_i^0$ and $x_i^1$ with some base $\beta$ such that $x_i^0 + \beta x_i^1 = x_i$. So instead of calculating $\sum x_i v_i$ we instead calculate
$$\sum \sum_{k=0}^1{} x_i^k \cdot \beta^k v_i = \sum v_i \sum x_i^k \beta^k = \sum v_i x_i$$

Now breaking $\beta^kv_i$ into its rational and fraction parts as $\omega_i^k$ and $\theta_i^k + \epsilon_i^k$ such that:
$$\sum x_i^k \beta^kv_i = \sum x_i^k \cdot \omega_i^k + \sum x_i^k \cdot (\theta_i^k + \epsilon_i^k)$$
we are able to limit the error accumulated due to precision loss, $\sum x^k_i \epsilon_i^k$ by limiting the value of $x_i^k$ to base $\beta$.

In practice, since $q_i$ usually is $<= 60$ it suffices to set $k = 2$ and $\beta = 2^{\frac{\log_{q_i}}{2}}$.

**FastConvertPOverQ**

Refer to appendix E of 2021/204.

**ScaleAndRound**

Given $x \in R_{QP}$ we need to calculate:
$$[\lceil\frac{t}{P}x\rfloor]_Q$$

We implement complex scaling of HPS 2018/117.

Recall that
$$x = \sum x_i (\frac{QP}{q_i})^{-1}_{q_i} \frac{QP}{q_i} + \sum x_k (\frac{QP}{p_k})^{-1}_{p_k} \frac{QP}{p_k} - vQP$$
therefore,
$$\frac{t}{P} x = \sum x_i (\frac{QP}{q_i})^{-1}_{q_i} \frac{tQ}{q_i} + \sum x_k (\frac{QP}{p_k})^{-1}_{p_k} \frac{tQ}{p_k} - vtQ$$
To calculate $[\lceil\frac{t}{P}x\rfloor]_Q$, let's separate first and second part. We will ignore $-vtQ$ since it vanishes over $Q$.

Since $q_i$ divides $Q$, $x_i (\frac{QP}{q_i})^{-1} \frac{tQ}{q_i}$ will not have a fractional part. Thus there's no rounding operation in the first part. Now notice that for any $q_j$ the first part equals,
$$[\sum x_i (\frac{QP}{q_i})^{-1}_{q_i} \frac{tQ}{q_i}]_{q_j} = [x_j (\frac{QP}{q_j})^{-1}_{q_j} \frac{tQ}{q_j}]_{q_j}$$
This is because $\forall x_i (\frac{QP}{q_i})^{-1} \frac{tQ}{q_i} | q_i \neq q_j$ $q_j$ exists in numerator thus the term vanishes over $q_j$.

To calculate the second part, ie. $[\lceil \sum x_k (\frac{QP}{p_k})^{-1}_{p_k} \frac{tQ}{p_k} \rfloor]_{q_j}$ we will have to separate calculation of fractional part and rational part. For a single term
$$x_k (\frac{QP}{p_k})^{-1}_{p_k} \frac{tQ}{p_k} $$
let $\omega_k = \lfloor \frac{(\frac{QP}{p_k})^{-1}_{p_k} tQ}{p_k} \rfloor$ be the rational part and $\theta_k = \frac{(\frac{QP}{p_k})^{-1}_{p_k} tQ \mod {p_k}}{p_k}$  be the rational part. Hence, the second part can be re-written as: 
$$[\sum x_k \cdot \omega_k + \sum x_k \theta_k]_{q_j}$$
Finally, 
 $$[\lceil\frac{t}{P}x\rfloor]_{q*j} = [x_j (\frac{QP}{q_j})^{-1}*{q*j} \frac{tQ}{q_j}]*{q*j} + [\sum x_k \cdot \omega_k + \sum x_k \theta_k]*{q_j}$$
**LastModDown**

This uses the same algorithm as ApproximateModDown (i.e. algorithm 2 of 2018/931) but scales by $q_l$ where $q_l$ is last prime in moduli chain of $Q$.

On input $x \in R_Q$, output $x' \in R_{Q'}$ where $Q'$ equals $Q$ after dropping last prime $q_l$ such that
$$x'_i = [\lceil \frac{1}{q_l} x_i \rfloor]_{Q'}$$

-   Input can either be in Evaluation or Coefficient form. Try to accommodate it into a single function. If not, create two separate functions as LastModDownCoeff and LastModDownEval

**ApproximateModDown**

Given input $x \in R_{QP}$ outputs $x' \in R_Q$ such that
$$x'_i = [\lceil \frac{1}{P} x_i \rfloor]_Q$$

-   Implement algorithm 2 of 2018/931
-   This is only used in hybrid key switching where polynomials are in evaluation representation (after inner product). Hence, the implementation must accept inputs in evaluation form.
