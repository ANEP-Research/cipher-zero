#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Polynomial {
    pub data: Vec<i64>,
    pub p: i64,
}

const EBF: i64 = 1 << 9;

fn egcd(mut x: i64, mut y: i64) -> (i64, i64, i64) {
    let (mut a0, mut a1, mut b0, mut b1) = (1, 0, 0, 1);
    while y != 0 {
        let (q,r)  = (x / y, x % y);
        let (c, d) = ( a0 - q * a1, b0 - q * b1  ); 
        x  = y;
        y  = r;
        a0 = a1;
        a1 = c;
        b0 = b1; 
        b1 = d;
    }
    (x, a0, b0)
}

fn mod_inv(a: i64, m: i64) -> i64 {
    let (gcd, a, b) = egcd(a, m);
    ((a % m) + m)%m
}

impl Polynomial {
    pub fn new(data: Vec<i64>, p: i64) -> Self {
        Self {
            data,
            p,
        }
    }

    pub fn deg(&self) -> usize {
        self.data.len()-1
    }

    pub fn mult(&self, other: &Polynomial) -> Polynomial {
        assert_eq!(self.p, other.p);
        let mut res = vec![0;self.deg()+other.deg()+1];
        for i in 0..self.data.len() {
            for j in 0..other.data.len() {
                res[i+j] += (self.data[i]*other.data[j])%self.p;
                res[i+j] %= self.p;
            }
        }
        Polynomial::new(res, self.p)
    }

    pub fn rem(&self, other: &Polynomial) -> (Polynomial, Vec<i64>) {
        let mut res = self.data.clone();
        assert_eq!(self.p, other.p);
        let mut quot = vec![];
        for i in (other.deg()..=self.deg()).rev() {
            let d = (mod_inv(other.data[other.deg()], EBF)*res[i])%EBF;
            quot.push(d);
            for j in 0..=other.deg() {
                res[i-j] -= d*other.data[other.deg()-j];
                res[i-j] %= self.p;
                res[i-j] += self.p;
                res[i-j] %= self.p;
            }
        }
        quot.reverse();
        (Polynomial::new(res, self.p), quot)
    }

    pub fn extend(&mut self, p: i64) {
        assert_eq!(p % self.p, 0);
        self.p = p;
    }
}

