# %%
# Example usage for all scripts
from src.dvq.statistical import (
    average_kmer_jaccard_similarity,
    similarity_wen,
    moment_of_inertia,
    persistence_homology,
    compare_persistence_homology,
    denq_entropy_generalised,
    calculate_deng_entropies_multiprocess,
    deng_KL_divergence,
    kl_divergence,
    js_divergence,
    wasserstein_distance,
    wasserstein_distance_matrix
)

from src.dvq.functional import (
    find_genomic_repeats,
    power_spectrum,
    spectral_distance,
    power_spectra_batch,
    wavelet_energy,
    wavelet_distance
)


# %%
def main():

    # Test sequences
    seqs = [
        "ACGTACGTACGTACGTACGT",
        "AAAAGGGGCCCC",
        "ATCGATCGATCGATCGATCG",
        "AGCTTTCGAAGCTTTCGAAGCTTTCGAA",
        "CGTAGCTAGCTAGCTAGCTA",
        "TTGACGTCGTAGCTAGCTAAGT",
    ]

    repeats = find_genomic_repeats(seqs)
    print(repeats)

    # tests for Statistical module 
    # Example sequences
    seq_1 = "ACGT" * 1000
    seq_2 = "GCTA" * 1000
    seqs = [seq_1, seq_2, "TTAA" * 1000]

      
    # Test Deng Entropy
    print("\n--- Testing Deng Entropy ---")
    entropy = denq_entropy_generalised(seq_1)
    print(f"Deng entropy for sequence: {entropy}")

    entropies = calculate_deng_entropies_multiprocess(seqs)
    print(f"Entropies for sequences: {entropies}")

    # Test Deng KL Divergence
    print("\n--- Testing Deng KL Divergence ---")
    kl_div = deng_KL_divergence(seq_1, seq_2)
    print(f"Deng KL divergence: {kl_div}")
    
    # Test Standard KL Divergence
    print("\n--- Testing Standard KL Divergence ---")
    kl_div_standard = kl_divergence(seq_1, seq_2)
    print(f"Standard KL divergence: {kl_div_standard}")

    # Test Wen's Method
    print("\n--- Testing Wen's method ---")
    similarity = similarity_wen(seq_1, seq_2)
    print(f"Similarity using Wen's method: {similarity}")

    inertia = moment_of_inertia(seq_1)
    print(f"Moment of inertia for sequence: {inertia}")

    # Test Persistent Homology
    print("\n--- Testing Persistent Homology ---")
    persistence = persistence_homology(seq_1, plot=True)
    print(f"Persistent homology for sequence: {persistence}")

    homology_distance = compare_persistence_homology(seq_1, seq_2)
    print(f"Persistent homology distance: {homology_distance}")

    # Test k-mer functions
    print("\n--- Testing k-mer functions ---")
    overlap_group = average_kmer_jaccard_similarity(seqs, seqs)
    print(f"Self similarity: {overlap_group}")

    overlap_comparison = average_kmer_jaccard_similarity(seq_1, seq_2)
    print(f"Kmer similarity: {overlap_comparison}")
    
    # Test JS Divergence
    print("\n--- Testing JS Divergence ---")
    js = js_divergence(seq_1, seq_2)
    print(f"JS divergence: {js}")

    # Test Wasserstein Distance
    print("\n--- Testing Wasserstein Distance ---")
    w_dist = wasserstein_distance(seq_1, seq_2)
    print(f"Wasserstein distance: {w_dist}")

    w_matrix = wasserstein_distance_matrix(seqs)
    print(f"Wasserstein distance matrix:\n{w_matrix}")

    # Test Fourier
    print("\n--- Testing Fourier ---")
    spectrum = power_spectrum(seq_1)
    print(f"Power spectrum shape: {spectrum.shape}")

    spectra = power_spectra_batch(seqs)
    print(f"Power spectra batch shape: {spectra.shape}")

    spec_dist = spectral_distance(seq_1, seq_2)
    print(f"Spectral distance (euclidean): {spec_dist}")

    spec_dist_corr = spectral_distance(seq_1, seq_2, metric="correlation")
    print(f"Spectral distance (correlation): {spec_dist_corr}")

    # Test Wavelet
    print("\n--- Testing Wavelet ---")
    energy = wavelet_energy(seq_1)
    print(f"Wavelet energy vector: {energy}")

    wav_dist = wavelet_distance(seq_1, seq_2)
    print(f"Wavelet distance (euclidean): {wav_dist}")

    wav_dist_corr = wavelet_distance(seq_1, seq_2, metric="correlation")
    print(f"Wavelet distance (correlation): {wav_dist_corr}")

    return

if __name__ == "__main__":
    main()


## TODO: 
# 
# Add features:
# JS divergence <- done
# Wasserstein distance <- done
# Fourier <- done
# Wavelet <- done

# %%
