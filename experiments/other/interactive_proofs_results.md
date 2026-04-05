# Interactive Proofs / Short Certificates for nth Prime: Results

**Script:** interactive_proofs.py

## What Was Tested
Whether verifying "p is the nth prime" is faster than computing it. Five sub-investigations: Pratt certificates + counting certificates, SNARG-style verification of pi(x) = n, proof compression analysis, Arthur-Merlin / hint-based protocols, and verifiable computation certificates.

## Key Findings
- Pratt certificates prove primality in O(log^2(p)) but certifying the INDEX (that exactly n-1 primes exist below p) requires a pi(x) certificate
- A pi(x) certificate would need to encode O(sqrt(x)) bits of zeta-zero information
- SNARG-style verification could work in theory but the circuit for pi(x) has size O(x^{2/3}), making the proof O(polylog) to verify but O(x^{2/3}) to generate
- Arthur-Merlin protocols with hints: the hint must contain ~170 bits of information (the random part), which is itself the hard part
- Verification is strictly easier than computation, but the certificate size is still O(sqrt(x)) in the worst case

## Verdict
**CLOSED** -- Failure Mode: E (Equivalence)

## One-Line Summary
Verifying p(n) requires a certificate encoding ~170 bits of zeta-zero information; verification is easier than computation but not polylog.
