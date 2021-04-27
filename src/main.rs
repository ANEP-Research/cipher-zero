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
    let mut p_t = p.clone();
    p_t.extend(EBF);
    for i in 0..=p_t.deg() { 
        p_t.data[i] = substitute(p_t.data[i]); 
    }
    p_t.rem(&k)
}

fn decrypt(r: Polynomial, q: Vec<i64>, g: Graph) -> Polynomial {
    let k = g.generate_key();
    let mut res = r.data.clone();
    assert_eq!(res.len() - k.data.len() <= 1, true);
    for i in 0..=(r.deg() - k.deg()) {
        for j in 0..=k.deg() {
            res[i+j] += q[i]*k.data[j];
            res[i+j] %= EBF;
        }
    }
    for i in 0..res.len() {
        res[i] = inv_substitute(res[i]);
        res[i] %= BF;
    }
    Polynomial::new(res, BF)
}

fn main() {
    // TODO: Remove below test
    let mut g = Graph::new(5);
    g.add_edge(1, 2);
    g.add_edge(3, 2);
    g.add_edge(2, 4);
    let (c, q) = encrypt(Polynomial::new(vec!['A' as i64,'B' as i64,'C' as i64,'D' as i64,'E' as i64,'F' as i64,'G' as i64,'H' as i64], BF), g.clone());
    dbg!(c);
}