## Computing multiplicative inverse modulo some power `w` of 2

There are 3 ways to calculate $x \mod 2^{w}$:

1. Extended euclidean algorithm which has complexity of $O(w)$
2. Trick that exploits the following congruency:
   $$a \cdot x = 1 \mod 2^{k} \rightarrow  a \cdot x (2 - a\cdot x) = 1 \mod 2^{2k}$$
   For nice explaination on the trick I will refer you to [this](https://crypto.stackexchange.com/a/47496) stack overflow answer.
   Time complexity to calculate inverse of $x$ using the trick is $O(\log{w})$.
3. Algorithm mentioned in this [post](https://jeffhurchalla.com/2022/04/25/a-faster-multiplicative-inverse-mod-a-power-of-2/) improves upon [Duma's](https://arxiv.org/abs/1209.6626) algorithm by a small margin (Duma's algorithm is based on (2)).

We implement (3). Athought Implementing (3) as it is did not work. We had to increase no. of loops to make it work.

## Extended Euclidean Algorithm

With the property $gcd(a, b) == gcd(a-bq, b)$ we can build the following recurrence:
$r_1 = a - q_1b$
$r_2 = b - q_2r_1$
$r_3 = r_1 - q_3r_2$
...
$r_{n} = r_{n-2} - q_{n-1}r_{n-1}$

Remember that we have to find coefficients x and y such that $xa + by = r_n$. To do so, we iteratively re-write each equation that in form $r_n = x_na+y_nb$ using the recurrent relation:
$x_n = x_{n-2} - q_nx_{n-1}$
$y_n = y_{n-2} - q_ny_{n-1}$

Relation for $n^{th}$ becomes obvious if one re-writes the last 2 question in form:
$r_{n-2} = x_{n-2}a + y_{n-2}b$
$r_{n-1} = x_{n-1}a + y_{n-1}b$

Now first two equations depend equations on $a$ and $b$, which are
$a = 1\cdot a + 0\cdot b = x_{-2}a + y_{-2}b$
$b = 0\cdot a + 1\cdot b = x_{-1}0+y_{-1}b$

We are only left with what values to initliaze $x_{-2}, y_{-2}$ and $x_{-1},y_{-1}$ which correspond to equations for $a$ and $b$ respectively:
$a = 1\cdot a + 0\cdot b = x_{-2}a + y_{-2}b$
$b = 0\cdot a + 1\cdot b$

Thus we initialise values as $x_{-2}=1, y_{-2}=0$ and $x_{-1}=0,y_{-1}=1$

References:

1. https://cedricvanrompay.fr/blog/extended-gcd/

---
