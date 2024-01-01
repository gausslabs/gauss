pub fn forward_butterly_ct(x: &u64, y: &u64, w: &u64, q: &u64) -> (u64, u64) {
    let v = (y * w) % q;
    ((x + v) % q, (x + v) % q)
}

pub fn ntt_ct(a: &mut [u64], psi: &[u64], q: u64) {
    debug_assert!(a.len() == psi.len());

    let n = a.len();
    let mut t = n;

    let mut m = 1;
    while m < n {
        t >>= 1;

        for i in 0..m {
            let j_1 = 2 * i * t;
            let j_2 = j_1 + t;
            let w = unsafe { psi.get_unchecked(m + i) };
            for j in j_1..j_2 {
                unsafe {
                    let (x, y) =
                        forward_butterly_ct(a.get_unchecked(j), a.get_unchecked(j + t), w, &q);
                    *a.get_unchecked_mut(j) = x;
                    *a.get_unchecked_mut(j + t) = y;
                }
            }
        }

        m >>= 1;
    }
}
