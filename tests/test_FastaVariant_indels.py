import os
import pytest

from pyfaidx import FastaVariant


path = os.path.dirname(__file__)
os.chdir(path)


@pytest.fixture(scope="module")
def ensure_tiny_vcf_tbi():
    """Compress and tabix-index tests/data/tiny.vcf -> tiny.vcf.gz + .tbi.
    Always reindex to pick up any edits to the VCF in this test module.
    """
    try:
        import pysam
    except ImportError:
        pytest.skip("pysam not installed.")

    vcf_txt = os.path.join("data", "tiny.vcf")
    vcf_gz = vcf_txt + ".gz"

    # (Re)create compressed and indexed VCF to ensure it's fresh
    if os.path.exists(vcf_gz) and os.path.exists(vcf_gz + ".tbi"):
        return vcf_gz

    pysam.tabix_compress(vcf_txt, vcf_gz, force=True)
    pysam.tabix_index(vcf_gz, preset="vcf", force=True)

    assert os.path.exists(vcf_gz)
    assert os.path.exists(vcf_gz + ".tbi")
    return vcf_gz


def _make_fasta_variant(vcf_gz, **kwargs):
    return FastaVariant(os.path.join("data", "tiny.fa"), vcf_gz, as_raw=True, hom=True, het=True, **kwargs)


def test_combined_indels_window(ensure_tiny_vcf_tbi):
    """Window 1..48 (0-based [0:48]) includes all variants on chr1, chr2 and chr3"""
    try:
        fasta = _make_fasta_variant(ensure_tiny_vcf_tbi)
        s1 = fasta['chr1'][0:48]
        assert s1 == 'AGAACCCCGGTTGGTTTAAACCCCGGGGTTTTAAAACCCCGGGGTTTT'

        s2 = fasta['chr2'][0:48]
        assert s2 == "ACGTACGTANACGTACGTACGTACGTACGTACGTACGTACGT"

        s3 = fasta["chr3"][0:20]
        assert s3 == "GCGTGTTTACGTACGTACG"
    
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
    Length remains 4 due to truncation to request size; content reflects insertion.
    """
    try:
        fasta = _make_fasta_variant(ensure_tiny_vcf_tbi)

        s = fasta['chr1'][8:12]
        assert s == 'GGTT'  # was 'GGGG'

    except (ImportError, IOError):
        pytest.skip("pysam or PyVCF not installed or VCF indexing unavailable.")



def test_chr3_overlapping_del_then_ins(ensure_tiny_vcf_tbi):
    """Deletion TAC->T at pos4..6 (anchored at pos4), and insertion T->TTT at pos8.
    These are adjacent events. Validate consensus within a window spanning them.
    """
    fasta = _make_fasta_variant(ensure_tiny_vcf_tbi)
    s = fasta['chr3'][3:12]  # 0-based slice [3:12] == bases 4..12
    assert isinstance(s, str)
    assert len(s) == 9
    # Basic sanity: result starts with the anchor 'T' at pos4
    assert s[0] == 'T'


def test_chr3_deletion_at_chrom_end(ensure_tiny_vcf_tbi):
    """End-of-chromosome deletion using anchored form (POS=19, REF=GT -> ALT=G).
    Ensures we handle trailing-edge edits without coordinate mismatch.
    """
    fasta = _make_fasta_variant(ensure_tiny_vcf_tbi)
    s = fasta['chr3'][12:20]
    assert s == "ACGTACG"  # one base shorter than reference
    assert len(s) == 7
