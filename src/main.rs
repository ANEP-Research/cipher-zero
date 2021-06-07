mod polynomials;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Graph {
    pub n: usize,
    pub data: Vec<(usize, usize)>,
}

pub fn substitute(x: i64) -> i64 {
    2*x+1
}

pub fn inv_substitute(x: i64) -> i64 {
    assert_eq!(x % 2, 1);
    (x-1)/2
}

const EBF: i64 = 1 << 9;
const BF: i64 = 1 << 8;

use polynomials::Polynomial;

impl Graph {
    pub fn new(n: usize) -> Self {
        Self {
            n,
            data: vec![],
        }
    }

    pub fn add_edge(&mut self, u: usize, v: usize) {
        self.data.push((u, v));
    }

    pub fn expand_sequence(&self) -> Vec<i64> {
        let mut a = self.data.clone();
        a.sort();
        let mut res = vec![];
        for (u, v) in a {
            res.push(1 % BF);
            res.push((u as i64) % BF);
            res.push((v as i64) % BF);
        }
        res
    }

    pub fn generate_key(&self) -> Polynomial {
        let mut res = Polynomial::new(vec![1], EBF);
        for (u, v) in self.data.clone() {
            res = res.mult(&Polynomial::new(vec![1, substitute(u as i64)% EBF, substitute(v as i64) % EBF], EBF));
        }
        res
    }
}

fn encrypt(p: Polynomial, g: Graph) -> (Polynomial, Vec<i64>) {
    let k = g.generate_key();
    let s = g.expand_sequence();
    let mut p_t = p.clone();
    p_t.extend(EBF);
    for i in 0..=p_t.deg() {
        p_t.data[i] ^= s[i % s.len()];
        p_t.data[i] = substitute(p_t.data[i]);
    }
    let (mut c, mut q) = p_t.rem(&k);
    for i in 0..=c.deg() {
        c.data[i] ^= s[i % s.len()];
    }
    for i in 0..q.len() {
        q[i] ^= s[i % s.len()];
    }
    (c, q)
}

fn decrypt(r: Polynomial, q: Vec<i64>, g: Graph) -> Polynomial {
    let k = g.generate_key();
    let s = g.expand_sequence();
    let mut res = r.data.clone();
    assert_eq!(res.len() - k.data.len() <= 1, true);
    let mut q = q.clone();
    for i in 0..=r.deg() {
        res[i] ^= s[i % s.len()];
    }
    for i in 0..q.len() {
        q[i] ^= s[i % s.len()];
    }
    for i in 0..=(r.deg() - k.deg()) {
        for j in 0..=k.deg() {
            res[i+j] += q[i]*k.data[j];
            res[i+j] %= EBF;
        }
    }
    for i in 0..res.len() {
        res[i] = inv_substitute(res[i]);
        res[i] %= BF;
        res[i] ^= s[i % s.len()];
    }
    Polynomial::new(res, BF)
}

pub fn string2poly(s: String) -> Polynomial {
    let mut data = vec![];
    for c in s.chars() {
        data.push(c as i64);
    }
    Polynomial::new(data, BF)
}

pub fn poly2string(s: Polynomial) -> String {
    assert_eq!(s.p, BF);
    let mut res = String::new();
    for c in s.data.clone() {
        res.push((c as u8) as char);
    }
    res
}

pub fn poly2hex(s: Polynomial) -> String {
    let mut res = String::new();
    for c in s.data.clone() {
        res.push_str(&format!("{:03x}", c));
    }
    res
}

pub fn hex2poly(s: String, field: i64) -> Polynomial {
    assert_eq!(s.len() % 3, 0);
    let mut data = vec![];
    for i in 0..(s.len()/3) {
        let x = i64::from_str_radix(s.get(3*i..3*(i+1)).unwrap(), 16).unwrap();
        data.push(x);
    }
    Polynomial::new(data, field)
}

fn main() {
    // TODO: Remove below test
    let mut g = Graph::new(5);
    g.add_edge(1, 2);
    g.add_edge(3, 2);
    g.add_edge(2, 4);
    let plaintext = String::from("ABCDEFG");
    let (c, q) = encrypt(string2poly(plaintext.clone()), g.clone());
    println!("Cipher text: {}", poly2hex(c.clone()));
    println!("Quotient: {}", q[0]);
    assert_eq!(c.clone(), hex2poly(poly2hex(c.clone()), c.p));
    assert_eq!(plaintext, poly2string(decrypt(c, q, g)));
}