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
const BLOCK_SIZE: usize = 16;
const ROUNDS: usize = 8;

use polynomials::Polynomial;

impl Graph {
    pub fn new(n: usize) -> Self {
        Self {
            n,
            data: vec![],
        }
    }

    pub fn add_edge(&mut self, u: usize, v: usize) {
        assert_eq!(u <= self.n, true);
        assert_eq!(v <= self.n, true);
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
        res.pop();
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

fn get_round(p: Polynomial, g: Graph, num: usize) -> Polynomial {
    let seq = g.expand_sequence();
    let k = g.generate_key();
    assert_eq!(p.p, EBF); // Must be EBF
    let mut before = p.clone();
    for i in 0..before.data.len() {
        before.data[i] ^= seq[i % seq.len()];
    }
    let (c, mut q) = before.rem(&k);
    let mut new_data = c.data.get(0..k.deg()).unwrap().to_vec();
    new_data.append(&mut q);
    if num % 2 == 1 {
        new_data.reverse();
    }
    for i in 0..new_data.len() {
        new_data[i] ^= seq[i % seq.len()];
    }
    Polynomial::new(new_data, EBF)
}

fn get_round_inv(p: Polynomial, g: Graph, num: usize) -> Polynomial {
    let seq = g.expand_sequence();
    let k = g.generate_key();
    assert_eq!(p.p, EBF); // Must be EBF
    let mut data = p.data.clone();
    for i in 0..data.len() {
        data[i] ^= seq[i % seq.len()];
    }
    if num % 2 == 1 {
        data.reverse();
    }
    let c_len = k.deg();
    let c: Vec<i64> = data.get(0..c_len).unwrap().to_vec();
    let q: Vec<i64> = data.get(c_len..).unwrap().to_vec();
    let mut res = vec![0;data.len()];
    for i in 0..c.len() {
        res[i] += c[i];
        res[i] %= EBF;
    }
    for i in 0..=(p.deg() - k.deg()) {
        for j in 0..=k.deg() {
            res[i+j] += q[i]*k.data[j];
            res[i+j] %= EBF;
        }
    }
    for i in 0..res.len() {
        res[i] ^= seq[i % seq.len()];
    }
    Polynomial::new(res, EBF)
}

fn encrypt(p: Polynomial, g: Graph) -> Polynomial {
    let s = g.expand_sequence();
    let mut p_t = p.clone();
    p_t.extend(EBF);
    for i in 0..=p_t.deg() {
        p_t.data[i] ^= s[i % s.len()];
        p_t.data[i] = substitute(p_t.data[i]);
    }
    for num in 0..ROUNDS {
        p_t = get_round(p_t.clone(), g.clone(), num);
    }
    p_t
}

fn decrypt(p: Polynomial, g: Graph) -> Polynomial {
    let k = g.generate_key();
    let s = g.expand_sequence();
    let mut res = p.clone();
    let mut round = ROUNDS;
    while round > 0 {
        res = get_round_inv(res.clone(), g.clone(), round-1);
        round -= 1;
    }
    let mut r = res.data;
    for i in 0..r.len() {
        r[i] = inv_substitute(r[i]);
        r[i] %= BF;
        r[i] ^= s[i % s.len()];
    }
    Polynomial::new(r, BF)
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

const TEST_VECTOR: &'static str = "crack me plz";

#[test]
fn test_simple_encryption() {
    let mut g = Graph::new(5);
    g.add_edge(1, 2);
    g.add_edge(3, 2);
    g.add_edge(2, 4);
    g.add_edge(3, 5);
    g.add_edge(1, 3);
    let plaintext = String::from(TEST_VECTOR);
    let c = encrypt(string2poly(plaintext.clone()), g.clone());
    println!("Ciphertext: {}", poly2hex(c.clone()));
    assert_eq!(poly2hex(c.clone()), "1b61ae1be0d50120c404409f04218e1191dd");
    assert_eq!(c.clone(), hex2poly(poly2hex(c.clone()), c.p));
    assert_eq!(plaintext, poly2string(decrypt(c, g)));
}

#[test]
fn test_wrong_graph_decryption() {
    let mut g1 = Graph::new(5);
    g1.add_edge(1, 2);
    g1.add_edge(3, 2);
    g1.add_edge(2, 4);
    g1.add_edge(3, 5);
    g1.add_edge(1, 3);
    let mut g2 = Graph::new(5);
    g2.add_edge(1, 2);
    g2.add_edge(3, 2);
    g2.add_edge(2, 4);
    g2.add_edge(3, 5);
    g2.add_edge(1, 4);
    let plaintext = String::from(TEST_VECTOR);
    let c = encrypt(string2poly(plaintext.clone()), g1.clone());
    println!("Ciphertext: {}", poly2hex(c.clone()));
    assert_eq!(poly2hex(c.clone()), "1b61ae1be0d50120c404409f04218e1191dd");
    assert_eq!(c.clone(), hex2poly(poly2hex(c.clone()), c.p));
    assert_ne!(plaintext, poly2string(decrypt(c.clone(), g2.clone())));
    println!("Decrypted1(Real): {}", poly2hex(decrypt(c.clone(), g1.clone())));
    println!("Decrypted2(Wrong): {}", poly2hex(decrypt(c.clone(), g2.clone())));
}

fn main() {
    
}