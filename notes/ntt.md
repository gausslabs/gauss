## Numbertheoretic transform and negacylic wrapped convolution

### Positive wrapped convolution $(X^N-1)$ - DFT

Let $R$ denote ring of polynomials $\in Q/X^N-1$. One can multiply two polynomials $f(x) \cdot g(x)$ to output $h(x) \in R$ using an FFT algorithm.

To understand how, let $\omega$ be primitive $N^{th}$ root of unity in $Q$, that is $\omega^{N} = 1$ and $\omega^{k} \neq 1 \forall k < N$.

Now note that $X^N - 1$ splits as product of linear factors with powers of $\omega = [1, \omega, ..., \omega^{N-1}]$. Therefore we can write:

$$X^N - 1 = \prod_{i=0}^{N-1}{X - \omega^i} $$
Since $X - 1, X - \omega, ..., X - \omega^{N-1}$ are co-prime (none of them can be factorised any further), we can use CRT to define the following isomorphism
$$X^N - 1 \cong \prod X - \omega^1 $$
Let $R^S = \prod X - \omega^i$ where in by each element we refer to any constant $\in X - \omega^i$ (polynomial $\in R$ reduced by either of the linear factor is a constant).

Now notice that complexity of multiplication $\in R$ is $N^2$ but in $R^S$ it is $N$. This is because each element in $R^S$ is a constant polynomial. So if one can find an efficient way to reduce $f(X) \in R$ by linear factors (i.e. find an efficient map) in $R^S$, along with an inverse map, then they can map polynomials from $R$ to $R^S$ multiply their representation in $R^S$ in linear complexity and then map product $\in R^S$ to $R$.

**Mapping $R$ to $R^S$ and back**

Notice that given f$(X) \in R$, then one can write a linear polynomial $X - \alpha$:
$$f(X) = q(X)(X-\alpha) + r(X)$$
where $r(X)$ is unique representation of $f(X) \in \mod{X-\alpha}$. Additionally note that,
$$f(\alpha) -  q(\alpha)(\alpha-\alpha)= f(\alpha) =   r(\alpha)$$
Give $r(X)$ is a constant and equals $f(\alpha)$ at $X = \alpha$, $r(X) = f(\alpha)$.

Thus one can easily find unique representation of $f(X) \in R$ in linear polynomials of $R^S$ by evaluating $f(X)$ at powers of $\omega: [1, \omega, \omega^{N-1}]$.

Now the question is does there exist an efficient way to evaluate $f(X)$ at multiple powers of $\omega$ with less than naive time complexity of $N^2$. This is where FFT comes into picture. Using FFT one can evaluate $f(X)$ at different powers of $\omega$ with time complexity $N\log{N}$.

But how does one map back from $R^S$ to $R$? Given that $R \rightarrow R^R$, inverse map from constants in $R^S$ to coefficients of polynomial in $R$ is a simple interpolation. Again, IFFT, that is inverse FFT, can be used to interpolate in $N\log{N}$.

Note: Viewing FFTs via this lens makes it quite easy to recognise that FFT corresponds to matrix vector product where the matrix contains powers of $\omega$ and vector consists of coefficients of polynomial. In IFFT, it is again a matrix vector product with matrix equal to inverse of the one used in FFT and vector containing evaluations at different powers of $\omega$.

With FFT and IFFT we can multiply two polynomials $f(X)$ and $g(X)$ in $R$ as:
$$h(X) = IFFT(FFT(f(X))\otimes FFT(g(X)))$$
where $\otimes$ represents Hadamard product between two vectors.

### Negative wrapped convolution $(X^N+1)$ using DFT

Let $R'$ be the ring polynomial $Q/X^N+1$.

Notice that unlike $R, R'$ does not splits into linear factors with powers $N^{th}$ root of unity $\omega$. One way to tackle the problem is to first observe that $X^{2N} - 1 = (X^N-1) \cdot (X^N+1)$, which implies each polynomial in $R'$ has unique representation in $R_{2N}$. So to multiply in $R'$, one can pad polynomials in $R'$ with $N$ $0$s to make them of degree $2N$, multiply them in $R_{2N}$ with $N\log{N}$ time complexity, and then reduce the result by $X^N+1$ (which translates to replacing $X^N$ with -1). But using this naive method, we pay double for multiplication than one should pay for degree $N$.

Let $\rho$ be $2N^{th}$ root of unity. Therefore, $\omega = \rho^2$. This is because $\sqrt \omega^{2N} = \omega^N = 1$. Also note that $\rho^N = -1$ (to understand why, view $\rho$ as dividing unit circle in complex & real plane in $2N$ equal parts. Raising $\rho$ to $N$ only gets halfway through the circle, that is -1).

Notice that if,
$$h(X) \in R' = \sum_{i=0}^{N-1} a_i X^{i}$$
then,
$$h(\rho X) = \sum_{i=0}^{N-1} a_i \rho^i X^{i}$$

Now observe that if $h(X) = X^N + 1$ then $h(\rho X) = -(X^N-1)$.

Now let $h(X)$ be some random polynomial and let its representation in $R'$ be $r(X)$, then we can write the following equation:
$$h(X) = q(X)(X^N+1) + r(X)$$
Now apply the map $X \rightarrow \rho X$ as
$$h(\rho X) = -q(\rho X)(X^N-1) + r(\rho X)$$
where $r(\rho X)$ is representation of $h(\rho X)$ in $R$.

Since $r(\rho X)$ and $r(X)$ have same degree, mapping $r(\rho X)$ to $r(X)$ is a simple linear transformation as:
$$[1, \rho^{-1}, ..., \rho^{N-1}] \otimes r(\rho X) = r(X)$$

Now substituting $h(X)$ with $f(X)g(X)$,
$$f(X)g(X) = q(X)(X^N+1) + r(X)$$
and
$$f(\rho X)g(\rho X) = -q(\rho X)(X^N-1) + r(\rho X)$$

This implies that to multiply $f(X)\cdot g(X) \mod{X^N+1}$, we should

1. first find $f(\rho X)$ and $g(\rho X)$
2. then multiply them $\mod{X^N-1}$, to find $r(\rho X) = f(\rho X) \cdot g(\rho X) \mod{X^N-1}$
3. Map $r(\rho X)$ to $r(X)$ where $r(X) = f(X) \cdot g(X) \mod{X^N+1}$

In practice this translates to,

1. Mapping $f(X),g(X)$ to $f(\rho X),g(\rho X)$ is simple element wise product between coefficients of polynomials and vector $[1, \rho, ..., \rho^{N-1}]$
2. Calculating $r(\rho X) = f(\rho X) \cdot g(\rho X) \mod X^N-1  = IFFT(FFT(f(\rho X)) \otimes FFT(g(\rho X)))$.
3. Mapping $r(\rho X)$ to $r(X)$ is simple element wise multiplication $r(X) = [1, \rho^{-1}, ..., \rho^{N-1}] \otimes r(\rho X)$.

### Number theoretic transform (NTTs)

NTTs are generalisation of DFTs over coefficient ring $Z_q$ for some prime $q$ (i.e $Z_q/X^N+1$,$Z_q/X^N-1$). NTTs like DFTs use FFTs as a fast subroutine for forward and inverse mapping with the difference being that all operations are performed in ring $Z_q$ instead over rationals. For $N^{th}$ or $2N^{th}$ primitive root of unity to exists, we are restricted to certan choice $q$ for given $N$. In particular for,

1. Positive wrapped convolution, ie $Z_q/X^N-1$: There must exist $N^{th}$ primitive root of unity, that is $\omega_N^{N} = 1$ mod q.
2. Negative wrapped convolution, ie $Z_q/X^N+1$: There must exist $2N^{th}$ primitive root of unity, that is $\omega_{2N}^{2N} = 1$ mod q.

**NTTs Implementation**

In RLWE we are only interested in polynomial rings $R_q = Z_q/X^N+1$. Thus for this section we will only bother ourselves with NTT implementations for ring $R_q$.

As mentioned before, we can perform multiplication $f(X) \cdot g(X)$ where $f(X), g(X) \in R_q$ using "Positive wrapped convolution" (i.e. $\mod X^N-1$). To do so we:

1. Map $f(X)$, $g(X)$ to $f(\rho X), g(\rho X)$ where $\rho^{2N} = 1 \mod q$.
2. Calculate $r(\rho X) = f(\rho X) \cdot g(\rho X) \mod{Z_q/X^N-1}$ as $INTT(NTT(f(\rho X)g(\rho X)))$, where INTT/NTT are FFT algorithms that work over $Z_q$.
3. Map $r(\rho X)$ to $r(X)$.

Let $\psi = [1, \rho, \rho^2, ..., \rho^{N-1}]$ and let [f] represent coefficient vector of polynomial f(X). Then step (1) requires two element-wise products as $\psi \otimes [f]$ and $\psi \otimes [g]$ and step (3) requires a single element wise product: $\psi^{-1} \otimes [r(\rho X)]$. Additionally, for $\rho$ to exist $q$ must satisfy $q = 1 \mod 2N$.

Implementations of NTT and INTT are based on the approach outlined above, but apply optimisations such as

1. Merging el-wise multiplication with $\psi/\psi^{-1}$ in NTT/INTT and pipelining SO -> BO -> SO:
    1. Instead of multiplying $\psi/\psi^{-1}$ vector with coefficient vector of polynomials before & after NTT/INTT operations, it is possible to merge $\psi/\psi^{-1}$ multiplication NTT/INTT operations.
    2. In particular butterfly operation of Cooley-Tukey (CT) can be modified to multiply $\psi$ during NTT operation. Although, in practice we use vector of $\psi =[1, \rho, ..., \rho^{N-1}]$ stored in bit-reversed order, section 3 of [5] gives intuition of why merging works.

        Algorithm 1 of [5] is based on CT FFT algorithm that does not uses pre-computed roots instead computes roots on the fly, starting with $\omega = 1$. Since, $\rho = \omega_{N/2} = \sqrt{\omega_N}$ one can evaluate butterfly operation of CT with following sequence of root powers $\{\omega_{N/2}^{2j+1}\ | \forall j \in [0,N-1]\} = [\sqrt{\omega}, \sqrt{\omega}\omega, \sqrt{\omega}\omega^2, ..., \sqrt{\omega} \omega^{N-1}]$ to produce evaluations of $[f(\rho), f(\rho \omega), ... f(\rho \omega^{N-1})]$. This can be done by setting $\omega = \omega_{N/2}$ in algorithm 1.

    3. However. CT's butterfly cannot be adjusted for multiplication with $\psi^{-1}$ in INTT, but Gentleman-Sande (GS) DIF algorithm can be.
    4. CT, GS are tranformed to accept inputs in SO, BO and produce outputs in BO, SO respectively. Due to this, one can evaluate the entire process of $f = INTT_{GS}(NTT_{CT}(f))$ without any bit-reversal at run-time.
    5. For more information on optimisation related to merging and eliminating bit reversals check section 3 of [3].
2. Use of redundant representation of pre-computes (for ex, $\psi/\psi^{-1}$ stored in bit reversed order) for faster modular multiplication, such as barrett, shoup and montogomery multiplication, to speed modular arithmetic inside the butterfly routines of CT and GS. For more information refer to [2] and [4].

Note: For comprehensive view on all optimisations and their original sources I will direct you to section 2.2 of [2].

**Constant time NTTs**

Only required during secret key related operations. Have left this as a later task. But note that one only need to make the butterfly operations constant time.

### Finding $n^{th}$ primitive root of unity

Let $\alpha$ be primitive $n^{th}$ root of unity in $F_q$ for some prime q. Then $[1, \alpha, ..., \alpha^{n-1}]$ forms a subgroup in $Z_q$. Recall Langrange's group theorem that order of subgroup divides order of main group, then condition for $n^{th}$ primitive root of unity to exist in $F_q$ is:
$$n | (q - 1)$$

Notice that, due to Fermat's little theorem: for any $x \in F_q$, $x^{q-1} = 1 \mod{q}$
$$x^{\frac{q-1}{n}}$$
is $n^{th}$ root of unity for any $x \in F_q$

But we are interested in finding primitive $n^{th}$ root of unity. That is find $x$ such that $x^n = 1$ and for any $0 < k < n$, $x^k \neq 1$. We can further restrict $k$ to values that divide $n$ because if $x^n$ is not primitive root then $x^k$ can be primitive root only if $k$ divides $n$ since
$$x^n = (x^k)^{\frac{n}{k}} = 1$$

To confirm that $x$ such that $x^n = 1$ is $n^{th}$ primitive root of unity we must rule out that $x^{k} \neq 1 \forall k | n$

Assuming that $n$ is a power of 2, instead of checking $x^k \neq 1$ for each $k$ that divides $n$, we can just check $x^{n/2} \neq 1$. This is because, since only power of 2 can divide $n$, $k$ can only be a power of 2. For any power of two $k/2$ if $x^{k/2} = 1$ then its square is $(x^{k/2})^2 = 1$ as well. Therefore, if $x^{n/2} \neq 1$ then we can conclude that $x^{n/4} \neq 1$ as well. If it were the case that $x^{n/4} = 1$ and $x^{n/2} \neq 1$ will be contradiction because $(x^{n/4})^2 = x^{n/2} = 1^2 = 1$. We can recursively apply this argument to rule out that for any power of 2 $k < n$, $x^k \neq 1$, if $x^{n/2} \neq 1$. Thus $n$ is the smallest value for which $x^n = 1$.

We use the following procedure to find $n^{th}$ primitive root of unity [5]:

1. Randomly select $x \in F_q$.
2. Calculate $g = x^{\frac{q-1}{n}}$
3. If $g^{n/2} \neq 1$, then return $g$ otherwise back to (1)

### References

1. https://pages.cpsc.ucalgary.ca/~eberly/Courses/CPSC599/Exercises/exercise_06.pdf
2. Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based Cryptography - https://eprint.iacr.org/2016/504.pdf
3. https://www.cs.bham.ac.uk/~sinharos/papers/High-PerformanceIdealLattice-BasedCryptographyon8-bitAVR.pdf
4. FASTER ARITHMETIC FOR NUMBER-THEORETIC TRANSFORMS - https://arxiv.org/pdf/1205.2926.pdf
5. https://crypto.stackexchange.com/a/63616
