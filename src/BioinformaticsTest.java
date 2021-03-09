
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class BioinformaticsTest {

	@Test
	public void testSlowRNAScore() {
		assertEquals(1,Bioinformatics.slowRNAScore("ACCCCCU"));
		assertEquals(2,Bioinformatics.slowRNAScore("ACCCCGU"));
		assertEquals(3,Bioinformatics.slowRNAScore("ACUGAGCCCU"));
		assertEquals(4,Bioinformatics.slowRNAScore("AAUUGCGC"));
		
		// Artificial test cases with scores of 5, 6, and 7
		// that complete execution even with the slowRNAScore method.
		assertEquals(5, Bioinformatics.slowRNAScore("GACCUGGGCAUUA"));
		assertEquals(6, Bioinformatics.slowRNAScore("UACUGAGCUACGAU"));
		assertEquals(7, Bioinformatics.slowRNAScore("CGUAUGACGCUAGCCG"));
		
		// Slow, but will finish
		assertEquals(8,Bioinformatics.slowRNAScore("ACUGAGCCCUGUUAGCUAA"));

	}	

	@Test
	public void testFastRNAScore() {
		assertEquals(1,Bioinformatics.fastRNAScore("ACCCCCU"));
		assertEquals(2,Bioinformatics.fastRNAScore("ACCCCGU"));
		assertEquals(3,Bioinformatics.fastRNAScore("ACUGAGCCCU"));
		assertEquals(4,Bioinformatics.fastRNAScore("AAUUGCGC"));
		
		// Artificial test cases with scores of 5, 6, and 7
		// Test the same extra cases that were tested the slowRNAScore
		assertEquals(5, Bioinformatics.fastRNAScore("GACCUGGGCAUUA"));
		assertEquals(6, Bioinformatics.fastRNAScore("UACUGAGCUACGAU"));
		assertEquals(7, Bioinformatics.fastRNAScore("CGUAUGACGCUAGCCG"));
		
		assertEquals(8,Bioinformatics.fastRNAScore("ACUGAGCCCUGUUAGCUAA"));
		assertEquals(52,Bioinformatics.fastRNAScore("GGAUACGGCCAUACUGCGCAGAAAGCACCGCUUCCCAUCCGAACAGCGAAGUUAAGCUGCGCCAGGCGGUGUUAGUACUGGGGUGGGCGACCACCCGGGAAUCCACCGUGCCGUAUCCU"));
		assertEquals(68,Bioinformatics.fastRNAScore("AAAGAUCGGGUGAGAUAGUAGAGAUAGUAUGUGUCUCUCAUCUACUAUCGGGUAGAUUUCAUCUACUAUCGGGUAUAUCGGGUAAAAUCGGGUAAGAUUCUCUCUCAUCUACUGUGUCUCUCAUCUACUAUCGGGUAUAUCGGGUAAAAUCGGGUAA"));

	}

	@Test
	public void testSlowDNAScore() {
		assertEquals(4,Bioinformatics.slowDNAScore("ATCTGAT","TGCATA"));
		assertEquals(5,Bioinformatics.slowDNAScore("ATGTTAT","ATCGTAC"));
		assertEquals(7,Bioinformatics.slowDNAScore("ATATATAT","TATATATA"));
		
		// Additional tests that have a longest subsequence of 8, 9, 10 respectively
		assertEquals(8, Bioinformatics.slowDNAScore("GACACGAGAAA", "GTTTTCAGCAGAA"));
		assertEquals(9, Bioinformatics.slowDNAScore("CGAGAATGCAGA", "CAGCAGTTATGTCA"));
		assertEquals(10, Bioinformatics.slowDNAScore("TGAATCAGTAAGTG", "GTTAATGTAACTGCAT"));

	}

	@Test
	public void testFastDNAScore() {
		assertEquals(4,Bioinformatics.fastDNAScore("ATCTGAT","TGCATA"));
		assertEquals(5,Bioinformatics.fastDNAScore("ATGTTAT","ATCGTAC"));
		assertEquals(7,Bioinformatics.fastDNAScore("ATATATAT","TATATATA"));
		assertEquals(23,Bioinformatics.fastDNAScore("TCCCAGTTATGTCAGGGGACACGAGAATGCAGAGAC","AATTGCCGCCGTCGTTTTCAGCAGTTATGTCAGATC"));
		assertEquals(75,Bioinformatics.fastDNAScore("GCGCGTGCGCGGAAGGAGCCAAGGTGAAGTTGTAGCAGTGTGTCAGAAGAGGTGCGTGGCACCATGCTGTCCCCCGAGGCGGAGCGGGTGCTGCGGTACCTGGTCGAAGTAGAGGAGTTG","GACTTGTGGAACCTACTTCCTGAAAATAACCTTCTGTCCTCCGAGCTCTCCGCACCCGTGGATGACCTGCTCCCGTACACAGATGTTGCCACCTGGCTGGATGAATGTCCGAATGAAGCG"));

		// Tests from slowDNAScore
		assertEquals(8, Bioinformatics.fastDNAScore("GACACGAGAAA", "GTTTTCAGCAGAA"));
		assertEquals(9, Bioinformatics.fastDNAScore("CGAGAATGCAGA", "CAGCAGTTATGTCA"));
		assertEquals(10, Bioinformatics.fastDNAScore("TGAATCAGTAAGTG", "GTTAATGTAACTGCAT"));
		// Three additional tests that fastDNAScore can handle
		assertEquals(21, Bioinformatics.fastDNAScore("GCGTGCGCGGAAGGAGCCAAGGTGAAGTGTTG", "CTTGTGGAACCTACTTCCTGAAAATAACCTTCTGTCCTCCGAGCTCTCCGCACCAGCG"));
		assertEquals(31, Bioinformatics.fastDNAScore("ACTTGCCACATGCGTGCGCGGAAGGAGCCAAGGTGAAGTGTTATCGGA", "ACCGCTTGTGGAACCTACTTCCTGAAAATAACCTTCTGTCCTCCGAGCTCTCCGCCAAACCAGCGATT"));
		assertEquals(59, Bioinformatics.fastDNAScore("AATACAGTATCGTTTGCCACATGCGTGAGAGTCGCTTCAAGGTGAGAAGATTATCATGATTAAGTGGAGCCAAGGTGAATAGTGTTATCGAGA", "AATCGGTCCGCTTGATGGATAACCTACTTCACTAGAATGCAATAACCTTCTGTCCTCCGAGCTCTCCGACCAAACCAGCGAAAGAATACAGGT"));
	}
}
