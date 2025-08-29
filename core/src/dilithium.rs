use rand::prelude::*;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

use crate::error::Result;
use crate::utils::{get_hash_reader, get_randomness, hash_message};
use crate::{error::ThresholdError, params::DilithiumConfig};
use math::traits::FiniteField;

use math::prelude::*;

/// Represents a Dilithium key pair (public and private keys).
#[derive(Clone, Debug)]
pub struct DilithiumKeyPair<'a, FF: FiniteField> {
    pub public_key: DilithiumPublicKey<'a, FF>,
    pub private_key: DilithiumPrivateKey<'a, FF>,
}

impl<'a, FF: FiniteField> DilithiumKeyPair<'a, FF> {
    /// Initialize key pair.
    pub fn new(
        public_key: DilithiumPublicKey<'a, FF>,
        private_key: DilithiumPrivateKey<'a, FF>,
    ) -> Self {
        Self {
            public_key,
            private_key,
        }
    }
}

/// Dilithium public key containing matrix m and vector t.
#[derive(Clone, Debug, PartialEq)]
pub struct DilithiumPublicKey<'a, FF: FiniteField> {
    // TODO define a standalone type
    pub m: Vec<Vec<Polynomial<'a, FF>>>,
    pub t: PolynomialVector<'a, FF>,
    pub security_level: usize,
    pub config: DilithiumConfig,
}

impl<'a, FF: FiniteField> DilithiumPublicKey<'a, FF> {
    /// Initialize public key.
    pub fn new(
        m: Vec<Vec<Polynomial<'a, FF>>>,
        t: PolynomialVector<'a, FF>,
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
pub struct DilithiumPrivateKey<'a, FF: FiniteField> {
    pub(crate) s1: PolynomialVector<'a, FF>,
    pub(crate) s2: PolynomialVector<'a, FF>,
    pub dilithium: Dilithium,
}

impl<'a, FF: FiniteField> DilithiumPrivateKey<'a, FF> {
    /// Initialize private key.
    pub fn new(
        s1: PolynomialVector<'a, FF>,
        s2: PolynomialVector<'a, FF>,
        security_level: usize,
    ) -> Self {
        let dilithium = Dilithium::new(security_level);
        DilithiumPrivateKey { s1, s2, dilithium }
    }
}

/// Dilithium signature containing vectors z and h, and challenge c.
#[derive(Clone, Debug, PartialEq)]
pub struct DilithiumSignature<'a, FF: FiniteField> {
    pub z: PolynomialVector<'a, FF>,
    pub h: PolynomialVector<'a, FF>,
    pub c: Polynomial<'a, FF>,
}

impl<'a, FF: FiniteField> DilithiumSignature<'a, FF> {
    /// Initialize signature.
    pub fn new(
        z: PolynomialVector<'a, FF>,
        h: PolynomialVector<'a, FF>,
        c: Polynomial<'a, FF>,
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
    pub fn keygen<FF: FiniteField>(
        self,
        seed: Option<&[u8]>,
    ) -> DilithiumKeyPair<'static, FF> {
        let seed = match seed {
            Some(s) => s.to_vec(),
            None => {
                let mut rng = rand::thread_rng();
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
    pub fn sign<FF: FiniteField>(
        &self,
        message: &[u8],
        private_key: &DilithiumPrivateKey<'static, FF>,
        randomness: Option<&[u8]>,
    ) -> Result<DilithiumSignature<'static, FF>> {
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

            let s1 = private_key.s1.clone();
            // Compute response z = y + c * s1
            let cs1 = s1 * &c;
            let z = y + &cs1;

            // Check bounds
            if self.check_z_bounds(&z) {
                // Compute hint h
                let s2 = private_key.s2.clone();
                let h = self.compute_hint(&w, s2, &c);

                if self.check_h_bounds() {
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
    pub fn verify<FF: FiniteField>(
        &self,
        message: &[u8],
        signature: &DilithiumSignature<'static, FF>,
        public_key: &DilithiumPublicKey<'static, FF>,
    ) -> bool {
        // Check signature bounds
        if !self.check_signature_bounds(&signature.z, &signature.c) {
            return false;
        }

        // Hash message
        let mu = hash_message(message);

        // Recompute w'
        let w_prime = self.recompute_w(signature.clone(), public_key);

        // Extract high bits using hint
        let w1_prime = self.high_bits(&w_prime);

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
    fn expand_a<FF: FiniteField>(
        &self,
        rho: &[u8],
    ) -> Vec<Vec<Polynomial<'static, FF>>> {
        let mut a = Vec::with_capacity(self.config.k);

        for i in 0..self.config.k {
            let mut row = Vec::with_capacity(self.config.l);
            for j in 0..self.config.l {
                // Generate polynomial A[i,j] from rho, i, j
                let mut seed = rho.to_vec();
                seed.push(i as u8);
                seed.push(j as u8);

                let poly_coeffs = self.sample_uniform(&seed);
                row.push(poly![poly_coeffs]);
            }
            a.push(row);
        }

        a
    }

    /// Sample uniform polynomial from seed.
    /// TODO test it
    /// TODO update it
    /// TODO refactor
    fn sample_uniform(&self, seed: &[u8]) -> Vec<i32> {
        let mut reader = get_hash_reader(seed);

        let mut bytes = vec![0u8; N * 4];
        reader.read(&mut bytes);

        let mut coeffs = vec![0i32; N];
        for (i, c) in coeffs.iter_mut().enumerate().take(N) {
            let idx = i * 4;
            let val = u32::from_le_bytes([
                bytes[idx],
                bytes[idx + 1],
                bytes[idx + 2],
                bytes[idx + 3],
            ]);
            *c = (val % Q as u32) as i32;
        }

        coeffs
    }

    /// Sample secret vector s1 or s2.
    fn sample_s<FF: FiniteField>(
        &self,
        rho_prime: &[u8],
        s_type: &str,
        length: usize,
    ) -> PolynomialVector<'static, FF> {
        let polys = (0..length)
            .map(|i| {
                let mut seed = rho_prime.to_vec();
                seed.extend_from_slice(s_type.as_bytes());
                seed.push(i as u8);
                Polynomial::from(self.sample_eta(&seed))
            })
            .collect();
        poly_vec!(polys)
    }

    /// Sample polynomial with coefficients in [-eta, eta].
    /// TODO add proper testing
    /// TODO refactor
    fn sample_eta(&self, seed: &[u8]) -> Vec<i32> {
        let mut reader = get_hash_reader(seed);

        let mut bytes = vec![0u8; N];
        reader.read(&mut bytes);

        let mut coeffs = vec![0i32; N];
        for (i, c) in coeffs.iter_mut().enumerate().take(N) {
            let val = (bytes[i] as i32) % (2 * self.config.eta as i32 + 1)
                - self.config.eta as i32;
            *c = val.rem_euclid(Q);
        }

        coeffs
    }

    /// Sample mask vector y.
    fn sample_y<FF: FiniteField>(
        &self,
        randomness: &[u8],
        kappa: u16,
    ) -> PolynomialVector<'static, FF> {
        let polys = (0..self.config.l)
            .map(|i| {
                let mut seed = randomness.to_vec();
                seed.extend_from_slice(&kappa.to_le_bytes());
                seed.push(i as u8);
                poly!(self.sample_gamma1(&seed))
            })
            .collect();

        poly_vec!(polys)
    }

    /// Sample polynomial with coefficients in [-gamma1, gamma1].
    pub(crate) fn sample_gamma1<FF: FiniteField>(
        &self,
        seed: &[u8],
    ) -> Polynomial<'static, FF> {
        let mut reader = get_hash_reader(seed);

        let mut bytes = vec![0u8; N * 4];
        reader.read(&mut bytes);

        //let mut coeffs = vec![0i32; N];
        let coeffs: Vec<FF> = (0..N)
            .map(|i| {
                let idx = i * 4;
                let val = u32::from_le_bytes([
                    bytes[idx],
                    bytes[idx + 1],
                    bytes[idx + 2],
                    bytes[idx + 3],
                ]);
                let sample = val % (2 * self.config.gamma1 + 1);
                let coeff = sample as i32 - self.config.gamma1 as i32;
                <i32 as Into<FF>>::into(coeff)
            })
            .collect();

        poly!(coeffs)
    }

    /// Extract high-order bits from polynomial vector.
    fn high_bits<FF: FiniteField>(
        &self,
        v: &PolynomialVector<'static, FF>,
    ) -> PolynomialVector<'static, FF> {
        let gamma2: FF = self.config.gamma2.into();
        let double_gamma2: FF = (2 * self.config.gamma2).into();

        let result_polys = v
            .as_slice()
            .iter()
            .map(|poly| {
                let high_coeffs: Vec<FF> = poly
                    .coefficients()
                    .iter()
                    .map(|&coeff| (coeff + gamma2) / double_gamma2)
                    .collect();
                poly!(high_coeffs)
            })
            .collect();

        poly_vec!(result_polys)
    }

    /// Generate challenge polynomial from message hash and w1.
    fn generate_challenge<FF: FiniteField>(
        &self,
        mu: &[u8],
        w1: &PolynomialVector<'static, FF>,
    ) -> Polynomial<'static, FF> {
        let mut hasher = Shake256::default();
        hasher.update(mu);
        hasher.update(b"challenge");

        w1.as_slice().iter().for_each(|p| {
            p.coefficients()
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

        poly![coeffs]
    }

    fn check_z_bounds<FF: FiniteField>(
        &self,
        z: &PolynomialVector<'static, FF>,
    ) -> bool {
        z.norm_infinity() < self.config.gamma1 - self.config.beta
    }

    /// Compute hint vector h.
    fn compute_hint<FF: FiniteField>(
        &self,
        w: &PolynomialVector<'static, FF>,
        s2: PolynomialVector<'static, FF>,
        c: &Polynomial<'static, FF>,
    ) -> PolynomialVector<'static, FF> {
        let cs2 = s2 * c;
        let w_minus_cs2 = w.clone() - cs2;
        self.make_hint(&w_minus_cs2, w)
    }

    /// Create hint from two vectors.
    fn make_hint<'a, FF: FiniteField>(
        &self,
        _v1: &PolynomialVector<'a, FF>,
        _v2: &PolynomialVector<'a, FF>,
    ) -> PolynomialVector<'static, FF> {
        let result_polys = (0..self.config.k).map(|_| poly![0; N]).collect();
        poly_vec!(result_polys)
    }

    /// Check if hint h satisfies bound requirements.
    /// TODO remove it
    fn check_h_bounds(&self) -> bool {
        // TODO review it
        // Simplified bound check
        true
    }

    /// Check if signature components satisfy bound requirements.
    fn check_signature_bounds<FF: FiniteField>(
        &self,
        z: &PolynomialVector<'static, FF>,
        c: &Polynomial<'static, FF>,
    ) -> bool {
        z.norm_infinity() < self.config.gamma1 - self.config.beta
            && c.norm_infinity() <= self.config.tau as u32
    }

    /// Recompute w during verification.
    fn recompute_w<FF: FiniteField>(
        &self,
        signature: DilithiumSignature<'static, FF>,
        public_key: &DilithiumPublicKey<'static, FF>,
    ) -> PolynomialVector<'static, FF> {
        // w = m * z - c * t * 2^d
        let az = &public_key.m * &signature.z;
        let ct = public_key.t.clone() * signature.c;
        let scalar = 1u64 << self.config.d;
        let ct_scaled = ct * scalar;

        az - ct_scaled
    }
}
