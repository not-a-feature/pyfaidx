import os
import pytest

from pyfaidx import FastaVariant


path = os.path.dirname(__file__)
os.chdir(path)


@pytest.fixture(scope="module")
def ensure_tiny_vcf_tbi():
    """Compress and tabix-index tests/data/tiny.vcf -> tiny.vcf.gz + .tbi."""
    try:
        import pysam
    except ImportError:
        pytest.skip("pysam not installed.")

    vcf_txt = os.path.join("data", "tiny.vcf")
    vcf_gz = vcf_txt + ".gz"
    tbi = vcf_gz + ".tbi"

    if os.path.exists(vcf_gz) and os.path.exists(tbi):
        # Do not recreate if already present
        return vcf_gz
    
    # (Re)create compressed and indexed VCF to be safe
    pysam.tabix_compress(vcf_txt, vcf_gz, force=True)
    pysam.tabix_index(vcf_gz, preset="vcf", force=True)

    assert os.path.exists(vcf_gz)
    assert os.path.exists(tbi)
    return vcf_gz


def _make_fasta_variant(vcf_gz, **kwargs):
    return FastaVariant(os.path.join("data", "tiny.fa"), vcf_gz, as_raw=True, hom=True, het=True, **kwargs)


def test_combined_indels_window(ensure_tiny_vcf_tbi):
    """Window 1..48 (0-based [0:48]) includes all variants on chr1 and chr2."""
    try:
        fasta = _make_fasta_variant(ensure_tiny_vcf_tbi)
        s1 = fasta['chr1'][0:48]
        assert s1 == 'AGAACCCCGGTTGGTTTAAACCCCGGGGTTTTAAAACCCCGGGGTTTT'

        s2 = fasta['chr2'][0:48]
        assert s2 == "ACGTACGTANACGTACGTACGTACGTACGTACGTACGTACGT"

    except (ImportError, IOError):
        pytest.skip("pysam or PyVCF not installed or VCF indexing unavailable.")


def test_chr1_deletion_only_window(ensure_tiny_vcf_tbi):
    """Window 12..19 (0-based [11:19]) includes only the deletion at pos15 (TTA->T)."""
    try:
        fasta = _make_fasta_variant(ensure_tiny_vcf_tbi)
        s = fasta['chr1'][11:19]
        assert s == 'GTTTAA'  # 8 bp -> 6 bp (net -2)

    except (ImportError, IOError):
        pytest.skip("pysam or PyVCF not installed or VCF indexing unavailable.")


def test_chr1_insertion_only_window(ensure_tiny_vcf_tbi):
    """Window 9..12 (0-based [8:12]) includes only the insertion at pos10 (G->GTT).
    Length remains 4 due to truncation to request size; content reflects insertion."""
    try:
        fasta = _make_fasta_variant(ensure_tiny_vcf_tbi)

        s = fasta['chr1'][8:12]
        assert s == 'GGTT'  # was 'GGGG'

    except (ImportError, IOError):
        pytest.skip("pysam or PyVCF not installed or VCF indexing unavailable.")
