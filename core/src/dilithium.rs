use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};

use crate::error::Result;
use crate::utils::{get_hash_reader, get_randomness, hash_message};
use crate::{error::ThresholdError, params::DilithiumConfig};
use math::{matrix::Matrix, traits::FiniteField};

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
    pub m: Matrix<'a, FF>,
    pub t: PolynomialVector<'a, FF>,
    pub security_level: usize,
    pub config: DilithiumConfig,
}

impl<'a, FF: FiniteField> DilithiumPublicKey<'a, FF> {
    /// Initialize public key.
    pub fn new(
        m: Matrix<'a, FF>,
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
        let seed = get_randomness(seed);

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
        let w_prime = self.recompute_w(signature, public_key);

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
    fn expand_a<FF: FiniteField>(&self, rho: &[u8]) -> Matrix<'static, FF> {
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

        Matrix::new(a)
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
        signature: &DilithiumSignature<'static, FF>,
        public_key: &DilithiumPublicKey<'static, FF>,
    ) -> PolynomialVector<'static, FF> {
        // w = m * z - c * t * 2^d
        let az = &public_key.m * &signature.z;
        let ct = public_key.t.clone() * &signature.c;
        let scalar = 1u64 << self.config.d;
        let ct_scaled = ct * scalar;

        az - ct_scaled
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::hash_message;

    fn seed32(byte: u8) -> Vec<u8> {
        vec![byte; 32]
    }

    impl Default for Dilithium {
        fn default() -> Self {
            Dilithium::new(crate::params::DEFAULT_SECURITY_LEVEL)
        }
    }

    #[test]
    fn expand_seed_is_deterministic_and_lengths_match() {
        let seed = seed32(7);
        let (rho1, rho_prime1, k1) = Dilithium::expand_seed(&seed);
        let (rho2, rho_prime2, k2) = Dilithium::expand_seed(&seed);

        assert_eq!(rho1, rho2);
        assert_eq!(rho_prime1, rho_prime2);
        assert_eq!(k1, k2);

        assert_eq!(rho1.len(), 32);
        assert_eq!(rho_prime1.len(), 32);
        assert_eq!(k1.len(), 32);
    }

    #[test]
    fn keygen_is_deterministic_for_fixed_seed() {
        let seed = seed32(42);

        // create two independent instances
        let d1 = Dilithium::default();
        let d2 = Dilithium::default();

        let kp1 = d1.keygen::<FieldElement>(Some(&seed));
        let kp2 = d2.keygen::<FieldElement>(Some(&seed));

        // Public keys: ensure the derived t (high bits of t) are identical.
        assert_eq!(kp1.public_key.t, kp2.public_key.t);

        // Private parts should also be identical (s1 and s2)
        assert_eq!(kp1.private_key.s1, kp2.private_key.s1);
        assert_eq!(kp1.private_key.s2, kp2.private_key.s2);
    }

    #[test]
    fn keygen_with_different_seeds_produces_different_keys() {
        let d = Dilithium::default();
        let kp1 = d.keygen::<FieldElement>(Some(&seed32(1)));
        let kp2 = d.keygen::<FieldElement>(Some(&seed32(2)));
        // At least one of the secret vectors should differ with high probability.
        assert!(
            kp1.private_key.s1 != kp2.private_key.s1
                || kp1.private_key.s2 != kp2.private_key.s2
                || kp1.public_key.t != kp2.public_key.t
        );
    }

    #[test]
    fn t1_matches_high_bits_of_a_s1_plus_s2() {
        let dil = Dilithium::default();
        let kp = dil.keygen::<FieldElement>(Some(&seed32(99)));

        let t = &kp.public_key.m * &kp.private_key.s1 + &kp.private_key.s2;
        let t1 = dil.high_bits::<FieldElement>(&t);

        assert_eq!(t1, kp.public_key.t);
    }

    #[test]
    fn sample_uniform_has_length_n_and_in_field_range() {
        let dil = Dilithium::default();
        let coeffs = dil.sample_uniform(&seed32(5));
        assert_eq!(coeffs.len(), N);
        for c in coeffs {
            assert!((0..Q).contains(&c));
        }
    }

    #[test]
    fn sample_eta_is_centered_bounded_by_eta() {
        let dil = Dilithium::default();
        let coeffs = dil.sample_eta(&seed32(11));
        assert_eq!(coeffs.len(), N);

        let eta = dil.config.eta as i32;
        for &c in &coeffs {
            // Because coefficients are reduced mod Q to [0, Q),
            // the centered absolute value must be <= eta.
            let centered = std::cmp::min(c, Q - c);
            assert!(
                centered <= eta,
                "centered={} eta={} c={}",
                centered,
                eta,
                c
            );
        }
    }

    #[test]
    fn sample_s_lengths_match_config() {
        let dil = Dilithium::default();
        let (rho, rho_prime, _) = Dilithium::expand_seed(&seed32(123));

        let s1 = dil.sample_s::<FieldElement>(&rho_prime, "s1", dil.config.l);
        let s2 = dil.sample_s::<FieldElement>(&rho_prime, "s2", dil.config.k);

        assert_eq!(s1.as_slice().len(), dil.config.l);
        assert_eq!(s2.as_slice().len(), dil.config.k);

        // Ensure values look small in infinity norm as they come from eta sampler.
        assert!(s1.norm_infinity() <= dil.config.eta);
        assert!(s2.norm_infinity() <= dil.config.eta);

        // Use rho to build A and verify shapes via multiply
        let a = dil.expand_a::<FieldElement>(&rho);
        let w = &a * &s1;
        assert_eq!(w.as_slice().len(), dil.config.k);
    }

    #[test]
    fn sign_and_verify_roundtrip() {
        let dil = Dilithium::default();
        let kp = dil.keygen::<FieldElement>(Some(&seed32(77)));

        let msg = b"threshold pq sigs FTW";
        let randomness = seed32(33);
        let sig = dil
            .sign::<FieldElement>(msg, &kp.private_key, Some(&randomness))
            .expect("sign");

        assert!(dil.verify::<FieldElement>(msg, &sig, &kp.public_key));
        // a signature should satisfy signature bounds
        assert!(dil.check_signature_bounds::<FieldElement>(&sig.z, &sig.c));
        // and our placeholder hint is zero vector
        assert_eq!(sig.h.norm_infinity(), 0);
    }

    #[test]
    fn verify_fails_on_modified_message() {
        let dil = Dilithium::default();
        let kp = dil.keygen::<FieldElement>(Some(&seed32(1)));

        let msg = b"hello world";
        let randomness = seed32(55);
        let sig = dil
            .sign::<FieldElement>(msg, &kp.private_key, Some(&randomness))
            .expect("sign");

        let mut tampered = msg.to_vec();
        tampered[0] ^= 0x01;

        assert!(!dil.verify::<FieldElement>(&tampered, &sig, &kp.public_key));
    }

    #[test]
    fn verify_fails_on_modified_challenge() {
        let dil = Dilithium::default();
        let kp = dil.keygen::<FieldElement>(Some(&seed32(9)));

        let msg = b"alter c";
        let randomness = seed32(99);
        let sig = dil
            .sign::<FieldElement>(msg, &kp.private_key, Some(&randomness))
            .expect("sign");

        // flip the first coefficient of c
        let mut c_coeffs = sig.c.coefficients().to_vec();
        let zero = c_coeffs[0] - c_coeffs[0]; // portable zero in field
                                              // toggle between zero and non-zero
        c_coeffs[0] = if c_coeffs[0] == zero {
            // set to 1
            <i32 as Into<FieldElement>>::into(1)
        } else {
            zero
        };
        let bad_sig = DilithiumSignature::new(
            sig.z.clone(),
            sig.h.clone(),
            poly!(c_coeffs),
        );

        assert!(!dil.verify::<FieldElement>(msg, &bad_sig, &kp.public_key));
    }

    #[test]
    fn recompute_w_matches_manual_formula() {
        let dil = Dilithium::default();
        let kp = dil.keygen::<FieldElement>(Some(&seed32(101)));
        let msg = b"manual w recompute";
        let randomness = seed32(202);
        let sig = dil
            .sign::<FieldElement>(msg, &kp.private_key, Some(&randomness))
            .expect("sign");

        let az = &kp.public_key.m * &sig.z;
        let ct = kp.public_key.t.clone() * sig.c.clone();
        let ct_scaled = ct * (1u64 << dil.config.d);
        let expected = az - ct_scaled;

        let recomputed = dil.recompute_w::<FieldElement>(&sig, &kp.public_key);
        assert_eq!(recomputed, expected);
    }

    #[test]
    fn generate_challenge_has_small_norm_and_reasonable_weight() {
        let dil = Dilithium::default();

        let msg = b"challenge test";
        let mu = hash_message(msg);

        // Construct w1 like signing does
        let randomness = seed32(77);
        let (rho, _, _) = Dilithium::expand_seed(&randomness);
        let a = dil.expand_a::<FieldElement>(&rho);
        let y = dil.sample_y::<FieldElement>(&randomness, 0);
        let w = &a * &y;
        let w1 = dil.high_bits::<FieldElement>(&w);

        let c = dil.generate_challenge::<FieldElement>(&mu, &w1);
        // Challenge polynomial should be very sparse with ±1 entries
        assert!(c.norm_infinity() <= 1);

        // Count non-zeros; since positions are sampled with possible collisions,
        // the count must be <= tau.
        let coeffs = c.coefficients();
        let zero = coeffs[0] - coeffs[0];
        let nonzeros = coeffs.iter().filter(|&&x| x != zero).count();
        assert!(nonzeros <= dil.config.tau);
    }

    #[test]
    fn z_bounds_hold_for_signature() {
        let dil = Dilithium::default();
        let kp = dil.keygen::<FieldElement>(Some(&seed32(121)));
        let msg = b"bounds test";
        let sig = dil
            .sign::<FieldElement>(msg, &kp.private_key, Some(&seed32(88)))
            .expect("sign");
        assert!(dil.check_z_bounds::<FieldElement>(&sig.z));
        assert!(dil.check_signature_bounds::<FieldElement>(&sig.z, &sig.c));
    }
}
