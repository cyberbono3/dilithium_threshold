#[cfg(test)]
mod integration_tests {
    use core::{
        dilithium::Dilithium,
        error::Result,
        threshold::{ThresholdKeyShare, ThresholdSignature},
    };

    // TODO debug it properly and make it pass
    #[ignore]
    #[test]
    fn test_complete_workflow() -> Result<()> {
        // Set up test fixtures
        let threshold = 3;
        let participants = 5;
        let security_level = 3;
        let message = b"Test message for threshold signature";

        let ts = ThresholdSignature::new(
            threshold,
            participants,
            Some(security_level),
        )
        .expect("ThresholdSignatureFailed");
        //println!("threshold signature: {:?}", ts);

        // 1. Distributed key generation
        let key_shares = ts
            .distributed_keygen(None)
            .expect("Distributed keygen failed");
       // println!("key_shares: {:?}", key_shares);
        println!("key_shares_len: {:?}", key_shares.len());
        assert_eq!(key_shares.len(), participants);

        // Verify all shares have the same public key
        let public_key = &key_shares[0].public_key;
        for share in &key_shares {
            assert_eq!(share.public_key.m, public_key.m);
            assert_eq!(share.public_key.t, public_key.t);
        }

        // 2. Partial signing
        // let signing_participants: Vec<&ThresholdKeyShare> = key_shares.iter().take(threshold).collect();
        let signing_participants = &key_shares[..threshold];
        println!("signing_participants_len: {:?}", signing_participants.len());
        let mut partial_signatures = Vec::with_capacity(threshold);

        for share in signing_participants {
            let partial_sig = ts
                .partial_sign(message, share, None)
                .expect("Partial sign failed");
            let is_valid =
                ts.verify_partial_signature(message, &partial_sig, share);
            assert!(
                is_valid,
                "Partial signature from participant {} is invalid",
                share.participant_id
            );
            partial_signatures.push(partial_sig);
        }

        // 3. Signature combination
        let combined_signature = ts
            .combine_signatures(&partial_signatures, public_key)
            .expect("Combine signatures failed");

        // 4. Verification using standard Dilithium
        let dilithium = Dilithium::new(security_level);
        let is_valid =
            dilithium.verify(message, &combined_signature, public_key);

        assert!(is_valid, "Combined signature verification failed");

        Ok(())
    }
}
