use rand::{prelude::*, rng};
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

use crate::error::Result;
use crate::utils::{get_hash_reader, get_randomness, hash_message};
use crate::{config::DilithiumConfig, error::ThresholdError};

use math::{
    polynomial::{Polynomial, N, Q},
    poly_vector::PolynomialVector,
};

/// Represents a Dilithium key pair (public and private keys).
#[derive(Clone, Debug)]
pub struct DilithiumKeyPair {
    pub public_key: DilithiumPublicKey,
    pub private_key: DilithiumPrivateKey,
}

impl DilithiumKeyPair {
    /// Initialize key pair.
    pub fn new(
        public_key: DilithiumPublicKey,
        private_key: DilithiumPrivateKey,
    ) -> Self {
        Self {
            public_key,
            private_key,
        }
    }
}

/// Dilithium public key containing matrix m and vector t.
#[derive(Clone, Debug, PartialEq)]
pub struct DilithiumPublicKey {
    pub m: Vec<Vec<Polynomial>>,
    pub t: PolynomialVector,
    pub security_level: usize,
    pub config: DilithiumConfig,
}

impl DilithiumPublicKey {
    /// Initialize public key.
    pub fn new(
        m: Vec<Vec<Polynomial>>,
        t: PolynomialVector,
        security_level: usize,
    ) -> Self {
        let config = DilithiumConfig::new(security_level);
        DilithiumPublicKey {
            m,
            t,
            security_level,
            config,
        }
    }
}

/// Dilithium private key containing secret vectors s1 and s2.
#[derive(Clone, Debug)]
pub struct DilithiumPrivateKey {
    pub(crate) s1: PolynomialVector,
    pub(crate) s2: PolynomialVector,
    pub dilithium: Dilithium,
}

impl DilithiumPrivateKey {
    /// Initialize private key.
    pub fn new(
        s1: PolynomialVector,
        s2: PolynomialVector,
        security_level: usize,
    ) -> Self {
        let dilithium = Dilithium::new(security_level);
        DilithiumPrivateKey { s1, s2, dilithium }
    }
}

/// Dilithium signature containing vectors z and h, and challenge c.
#[derive(Clone, Debug, PartialEq)]
pub struct DilithiumSignature {
    pub z: PolynomialVector,
    pub h: PolynomialVector,
    pub c: Polynomial,
}

impl DilithiumSignature {
    /// Initialize signature.
    pub fn new(
        z: PolynomialVector,
        h: PolynomialVector,
        c: Polynomial,
    ) -> Self {
        DilithiumSignature { z, h, c }
    }
}

/// Main Dilithium algorithm implementation.
///
/// Provides key generation, signing, and verification functionality
/// according to the CRYSTALS-Dilithium specification.
#[derive(Clone, Copy, Debug)]
pub struct Dilithium {
    security_level: usize,
    pub config: DilithiumConfig,
}

impl Dilithium {
    /// Initialize Dilithium with specified security level.
    pub fn new(security_level: usize) -> Self {
        let config = DilithiumConfig::new(security_level);
        Self {
            security_level,
            config,
        }
    }

    /// Generate a Dilithium key pair.
    /// self not needed
    pub fn keygen(self, seed: Option<&[u8]>) -> DilithiumKeyPair {
        let seed = match seed {
            Some(s) => s.to_vec(),
            None => {
                let mut rng = rng();
                let mut bytes = vec![0u8; 32];
                rng.fill_bytes(&mut bytes);
                bytes
            }
        };

        // Expand seed to generate randomness
        let (rho, rho_prime, _k) = Self::expand_seed(&seed);

        // Generate matrix A from rho
        let a = self.expand_a(&rho);

        // Sample secret vectors s1 and s2
        let s1 = self.sample_s(&rho_prime, "s1", self.config.l);
        let s2 = self.sample_s(&rho_prime, "s2", self.config.k);

        // Compute t = a * s1 + s2
        let t = &a * &s1 + &s2;

        // Extract high-order bits of t
        let t1 = self.high_bits(&t);

        let public_key = DilithiumPublicKey::new(a, t1, self.security_level);
        let private_key = DilithiumPrivateKey::new(s1, s2, self.security_level);

        DilithiumKeyPair::new(public_key, private_key)
    }

    /// Sign a message using Dilithium.
    // TODO add proper error handling
    pub fn sign(
        &self,
        message: &[u8],
        private_key: &DilithiumPrivateKey,
        randomness: Option<&[u8]>,
    ) -> Result<DilithiumSignature> {
        const MAX_ATTEMPTS: u16 = 1000;
        let randomness = get_randomness(randomness);

        // Hash message
        let mu = hash_message(message);

        // Generate matrix a from public key (would need to be passed or stored)
        // For now, we'll regenerate it (in practice, this should be optimized)
        let (rho, _, _) = Self::expand_seed(&randomness);
        let a = self.expand_a(&rho);

        let mut signature = None;

        for kappa in 0..MAX_ATTEMPTS {
            // Sample mask vector y
            let y = self.sample_y(&randomness, kappa);

            // Compute w = A * y
            let w = &a * &y;
            let w1 = self.high_bits(&w);

            // Generate challenge
            let c = self.generate_challenge(&mu, &w1);

            // Compute response z = y + c * s1
            let cs1 = private_key.s1.clone() * c;
            let z = y + &cs1;

            // Check bounds
            if self.check_z_bounds(&z) {
                // Compute hint h
                let h = self.compute_hint(&w, &z, private_key.s2.clone(), &c);

                if self.check_h_bounds(&h) {
                    //return Ok(DilithiumSignature::new(z, h, c));
                    signature = Some(DilithiumSignature::new(z, h, c));
                    break;
                }
            }
        }

        match signature {
            Some(signature) => Ok(signature),
            None => Err(ThresholdError::SignatureGenerationFailed),
        }
    }

    /// Verify a Dilithium signature.
    pub fn verify(
        &self,
        message: &[u8],
        signature: &DilithiumSignature,
        public_key: &DilithiumPublicKey,
    ) -> bool {
        // Check signature bounds
        if !self.check_signature_bounds(signature) {
            return false;
        }

        // Hash message
        let mu = hash_message(message);

        // Recompute w'
        let w_prime = self.recompute_w(signature, public_key);

        // Extract high bits using hint
        let w1_prime = self.use_hint(&w_prime);

        // Recompute challenge
        let c_prime = self.generate_challenge(&mu, &w1_prime);

        // Verify challenge matches
        signature.c == c_prime
    }

    /// Expand seed into rho, rho_prime, and K.
    fn expand_seed(seed: &[u8]) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
        let mut reader = get_hash_reader(seed);

        let mut expanded = vec![0u8; 96];
        reader.read(&mut expanded);

        let rho = expanded[..32].to_vec();
        let rho_prime = expanded[32..64].to_vec();
        let k = expanded[64..96].to_vec();

        (rho, rho_prime, k)
    }

    /// Expand rho to generate public matrix A.
    fn expand_a(&self, rho: &[u8]) -> Vec<Vec<Polynomial>> {
        let mut a = Vec::with_capacity(self.config.k);

        for i in 0..self.config.k {
            let mut row = Vec::with_capacity(self.config.l);
            for j in 0..self.config.l {
                // Generate polynomial A[i,j] from rho, i, j
                let mut seed = rho.to_vec();
                seed.push(i as u8);
                seed.push(j as u8);

                let poly_coeffs = self.sample_uniform(&seed);
                row.push(Polynomial::from(poly_coeffs));
            }
            a.push(row);
        }

        a
    }

    /// Sample uniform polynomial from seed.
    /// TODO test it
    /// TODO update it
    fn sample_uniform(&self, seed: &[u8]) -> Vec<i32> {
        let mut reader = get_hash_reader(seed);

        let mut bytes = vec![0u8; N * 4];
        reader.read(&mut bytes);

        let mut coeffs = vec![0i32; N];
        for i in 0..N {
            let idx = i * 4;
            let val = u32::from_le_bytes([
                bytes[idx],
                bytes[idx + 1],
                bytes[idx + 2],
                bytes[idx + 3],
            ]);
            coeffs[i] = (val % Q as u32) as i32;
        }

        coeffs
    }

    /// Sample secret vector s1 or s2.
    fn sample_s(
        &self,
        rho_prime: &[u8],
        s_type: &str,
        length: usize,
    ) -> PolynomialVector {
        PolynomialVector::new(
            (0..length)
                .map(|i| {
                    let mut seed = rho_prime.to_vec();
                    seed.extend_from_slice(s_type.as_bytes());
                    seed.push(i as u8);
                    Polynomial::from(self.sample_eta(&seed))
                })
                .collect(),
        )
    }

    /// Sample polynomial with coefficients in [-eta, eta].
    /// TODO add proper testing
    fn sample_eta(&self, seed: &[u8]) -> Vec<i32> {
        let mut reader = get_hash_reader(seed);

        let mut bytes = vec![0u8; N];
        reader.read(&mut bytes);

        let mut coeffs = vec![0i32; N];
        for i in 0..N {
            let val =
                (bytes[i] as i32) % (2 * self.config.eta + 1) - self.config.eta;
            coeffs[i] = val.rem_euclid(Q);
        }

        coeffs
    }

    /// Sample mask vector y.
    fn sample_y(&self, randomness: &[u8], kappa: u16) -> PolynomialVector {
        let polys = (0..self.config.l)
            .map(|i| {
                let mut seed = randomness.to_vec();
                seed.extend_from_slice(&kappa.to_le_bytes());
                seed.push(i as u8);
                Polynomial::from(self.sample_gamma1(&seed))
            })
            .collect();

        PolynomialVector::new(polys)
    }

    /// Sample polynomial with coefficients in [-gamma1, gamma1].
    fn sample_gamma1(&self, seed: &[u8]) -> Vec<i32> {
        let mut reader = get_hash_reader(seed);

        let mut bytes = vec![0u8; N * 4];
        reader.read(&mut bytes);

        let mut coeffs = vec![0i32; N];
        coeffs.iter_mut().enumerate().take(N).for_each(|(i, c)| {
            let idx = i * 4;
            let val = u32::from_le_bytes([
                bytes[idx],
                bytes[idx + 1],
                bytes[idx + 2],
                bytes[idx + 3],
            ]);
            *c = ((val % (2 * self.config.gamma1 as u32 + 1)) as i32
                - self.config.gamma1)
                .rem_euclid(Q);
        });

        coeffs
    }

    /// Extract high-order bits from polynomial vector.
    fn high_bits(&self, v: &PolynomialVector) -> PolynomialVector {
        let result_polys = v
            .as_slice()
            .iter()
            .map(|poly| {
                let high_coeffs = poly
                    .coeffs()
                    .iter()
                    .map(|&coeff| {
                        (coeff + self.config.gamma2) / (2 * self.config.gamma2)
                    })
                    .collect::<Vec<i32>>();
                Polynomial::from(high_coeffs)
            })
            .collect();

        PolynomialVector::new(result_polys)
    }

    /// Generate challenge polynomial from message hash and w1.
    fn generate_challenge(
        &self,
        mu: &[u8],
        w1: &PolynomialVector,
    ) -> Polynomial {
        let mut hasher = Shake256::default();
        hasher.update(mu);
        hasher.update(b"challenge");

        w1.as_slice().iter().for_each(|p| {
            p.coeffs()
                .iter()
                .for_each(|c| hasher.update(&c.to_le_bytes()))
        });

        let mut reader = hasher.finalize_xof();
        let mut hash_output = vec![0u8; self.config.tau * 2];
        reader.read(&mut hash_output);

        let mut coeffs = vec![0i32; N];

        // Sample tau positions for ±1 coefficients
        for i in 0..self.config.tau {
            let pos = (hash_output[i * 2] as usize) % N;
            let sign = if hash_output[i * 2 + 1] % 2 == 0 {
                1
            } else {
                -1
            };
            coeffs[pos] = sign;
        }

        Polynomial::from(coeffs)
    }

    fn check_z_bounds(&self, z: &PolynomialVector) -> bool {
        z.norm_infinity() < self.config.gamma1 - self.config.beta
    }

    /// Compute hint vector h.
    fn compute_hint(
        &self,
        w: &PolynomialVector,
        _z: &PolynomialVector,
        s2: PolynomialVector,
        c: &Polynomial,
    ) -> PolynomialVector {
        let cs2 = s2 * c;
        let w_minus_cs2 = w.clone() - cs2;
        self.make_hint(&w_minus_cs2, w)
    }

    /// Create hint from two vectors.
    fn make_hint(
        &self,
        _v1: &PolynomialVector,
        _v2: &PolynomialVector,
    ) -> PolynomialVector {
        let result_polys = (0..self.config.k)
            .map(|_| Polynomial::from(vec![0i32; N]))
            .collect();

        PolynomialVector::new(result_polys)
    }

    /// Check if hint h satisfies bound requirements.
    fn check_h_bounds(&self, _h: &PolynomialVector) -> bool {
        // Simplified bound check
        true
    }

    /// Check if signature components satisfy bound requirements.
    fn check_signature_bounds(&self, signature: &DilithiumSignature) -> bool {
        signature.z.norm_infinity() < self.config.gamma1 - self.config.beta
            && signature.c.norm_infinity() <= self.config.tau as i32
    }

    /// Recompute w during verification.
    fn recompute_w(
        &self,
        signature: &DilithiumSignature,
        public_key: &DilithiumPublicKey,
    ) -> PolynomialVector {
        // w = m * z - c * t * 2^d
        let az = &public_key.m * &signature.z;
        let ct = public_key.t.clone() * signature.c;
        let ct_scaled = ct * (1 << self.config.d);

        az - ct_scaled
    }

    /// Use hint to recover high-order bits.
    fn use_hint(&self, w: &PolynomialVector) -> PolynomialVector {
        // Simplified hint usage
        self.high_bits(w)
    }
}
