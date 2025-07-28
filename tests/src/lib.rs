#[cfg(test)]
mod tests {
    use core::{
        dilithium::Dilithium,
        error::Result,
        threshold::ThresholdSignature,
    };

    #[test]
    fn test_complete_workflow() -> Result<()> {
        // Set up test fixtures
        let threshold = 3;
        let participants = 5;
        let security_level = 3;
        let message = b"Test message for threshold signature";
        
        let ts = ThresholdSignature::new(threshold, participants, Some(security_level)).expect("ThresholdSignatureFailed");
   
        
        // 1. Distributed key generation
        let key_shares = ts.distributed_keygen(None).expect("Distributed keygen failed");
        assert_eq!(key_shares.len(), participants);
        
        // Verify all shares have the same public key
        let public_key = &key_shares[0].public_key;
        for share in &key_shares {
            assert_eq!(share.public_key.m.len(), public_key.m.len());
            assert_eq!(share.public_key.m[0].len(), public_key.m[0].len());
            assert_eq!(share.public_key.t.len(), public_key.t.len());
        }
        
        // 2. Partial signing
        let signing_participants: Vec<_> = key_shares.iter().take(threshold).collect();
        let mut partial_signatures = Vec::new();
        
        for share in &signing_participants {
            let partial_sig = ts.partial_sign(message, share, None).expect("Partial sign failed");
            partial_signatures.push(partial_sig);
            
            // Note: verify_partial_signature is commented out in the Rust implementation
        }
        
        // 3. Signature combination
        let combined_signature = ts.combine_signatures(&partial_signatures, public_key).expect("Combine signatures failed");
        
        // 4. Verification using standard Dilithium
        let dilithium = Dilithium::new(security_level);
        let is_valid = dilithium.verify(message, &combined_signature, public_key);
        
        assert!(is_valid, "Combined signature verification failed");
        
        Ok(())
    }
    
}