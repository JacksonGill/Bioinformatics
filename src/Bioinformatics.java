import java.util.Arrays;

public class Bioinformatics {

	/**
	 * Determines the length of the longest common subsequence
	 * between two DNA strands.
	 * 
	 * @param dna1 A DNA string containing only the letters A, T, C, and G.
	 * @param dna2 Another DNA string containing only the letters A, T, C, and G.
	 * @return Length of the longest common subsequence in a global alignment
	 *         of the two DNA strands.
	 */
	public static int slowDNAScore(String dna1, String dna2) {
		// The kick-off specifies that last indexes of each string so that both complete strings are checked
		return slowDNAScore(dna1, dna2, dna1.length() - 1, dna2.length() - 1);
	}

	/**
	 * Determines the length of the longest common subsequence between two DNA strands
	 * using recursion. The base case checks to make sure that i and j are -1. If so, 0 is
	 * returned because you are comparing to an empty string. Next, the result can be 
	 * defined by the maximum of case1, case2, and case3. The first case ignores the last 
	 * character in dna1. The second case ignores the last character in dna2. If the these
	 * last characters are the same, the third case is invoked. Since these characters are the
	 * same, there is a common subsequence and the length of the subsequence is incremented by 1.
	 * Next, the last character of dna1 and dna2 are ignored. This process continues recursively and the ending
	 * result is the longest common subsequence.
	 * @param dna1
	 * @param dna2
	 * @param i
	 * @param j
	 * @return int result (longest common subsequence)
	 */
	private static int slowDNAScore(String dna1, String dna2, int i, int j) {
		if(i == -1 || j == -1) return 0; //  Base case
		else {
			// Case 1: Ignore the last character of seq1. 
			//	This length is represented by S(seq1,seq2,i-1,j)
			int case1 = slowDNAScore(dna1, dna2, i-1, j);
			// Case 2: Ignore the last character of seq2. 
			//	This length is represented by S(seq1,seq2,i,j-1)
			int case2 = slowDNAScore(dna1, dna2, i, j-1);
			
			int result = Math.max(case1, case2);
			
			// Case 3: If the last characters of each string are the same, 
			// then they can be part of a common subsequence. In this case, 
			// the length of your subsequence increases by one, and you move 
			// on to the next character in each string. This quantity is 
			// represented by S(seq1,seq2,i-1,j-1)+1
			if(dna1.charAt(i) == dna2.charAt(j)) {
				int case3 = slowDNAScore(dna1,dna2,i-1,j-1)+1;
				result = Math.max(result, case3);
			}
			
			return result;
		}
	}	
	
	/**
	 * This method should return the exact same answers as your
	 * slowDNAScore method. However, it must use memoization to
	 * perform faster. This means that previously calculated
	 * results will be saved for reuse later. You should save
	 * your previously calculated results in a 2D array that
	 * is passed between all recursive method calls. Otherwise,
	 * the code is very similar to the slowDNAScore method.
	 * 
	 * Note: This method is simply the kick-off
	 *       method for your actual recursive algorithm. 
	 *       Also, do NOT use the substring method.
	 * 
	 * @param dna1 A DNA string containing only the letters A, T, C, and G.
	 * @param dna2 Another DNA string containing only the letters A, T, C, and G.
	 * @return Length of the longest common subsequence in a global alignment
	 *         of the two DNA strands.
	 */
	public static int fastDNAScore(String dna1, String dna2) {
		// known[i][j] will store results of fastDNAScore(known, dna1, dna2, i, j)
		int[][] known = new int[dna1.length()][dna2.length()];
		for(int i = 0; i < known.length; i++) {
			Arrays.fill(known[i], -1);
		}
		return fastDNAScore(known, dna1, dna2, dna1.length() - 1, dna2.length() - 1);
	}

	/**
	 * This method uses the same functionality as the slowDNAScore. However, this method
	 * uses memoization by storing all results in a 2D array representing the previously calculated
	 * results. This 2D array is passed between each recursive call and is checked to see if
	 * a result is already stored in the array. If so, the result is immediately returned and no
	 * more recursive calls are called. This saves time because the method does not have to calculate
	 * results that have already been previously calculated.
	 * @param known
	 * @param dna1
	 * @param dna2
	 * @param i
	 * @param j
	 * @return int result (longest common subsequence)
	 */
	private static int fastDNAScore(int[][] known, String dna1, String dna2, int i, int j) {
		if(i == -1 || j == -1) return 0; //  Base case
		else if(known[i][j] != -1) { // Base case: I already know the answer
			return known[i][j];
		} else {
			// Case 1: Ignore the last character of seq1. 
			int case1 = fastDNAScore(known, dna1, dna2, i-1, j);
			// Case 2: Ignore the last character of seq2. 
			int case2 = fastDNAScore(known, dna1, dna2, i, j-1);
			int result = Math.max(case1, case2);
			
			// Case 3: If the last characters of each string are the same, 
			// then they can be part of a common subsequence. In this case, 
			// the length of your subsequence increases by one, and you move 
			// on to the next character in each string. This quantity is 
			// represented by S(seq1,seq2,i-1,j-1)+1
			if(dna1.charAt(i) == dna2.charAt(j)) {
				int case3 = fastDNAScore(known, dna1,dna2,i-1,j-1)+1;
				result = Math.max(result, case3);
			}
			known[i][j] = result; // store calculated result in known for memoization
			return result;
		}
	}	
	
	/**
	 * Method determines the maximum number of base pair matches 
	 * in a folded RNA strand with no pseudo-knots (bases that
	 * are part of different loops cannot pair up). This slow
	 * version of the method simply uses the recursive score
	 * calculation, which makes it very inefficient.
	 * 
	 * Note: This method is simply the kick-off
	 *       method for your actual recursive algorithm. 
	 *       Also, do NOT use the substring method.
	 * 
	 * @param rna String representing RNA strand. Can only contain
	 *            the letters A, U, C, and G.
	 * @return Maximum number of base pairings in folded strand
	 *         with no pseudo-knots.
	 */
	public static int slowRNAScore(String rna) {
		// Initial start and end indexes comprise the entire String
		return slowRNAScore(rna, 0, rna.length() - 1);
	}	

	/**
	 * This methods determines the largest number of base pair matches given
	 * and RNA strand recursively. The base case checks to see if both indices
	 * representing the first and last characters don't overlap and that they are not equal
	 * to each other. If so, 0 is returned. Next, four cases are considered. The first ignores the
	 * first character, and the second ignores the last character. If both characters at their
	 * respective indices can be a base pair (using the helper function), the third case is invoked
	 * and increments the number of matching base pairs and ignores the first and last characters in the rna
	 * string. Next, the fourth case considers every possible way of splitting the rna string into two substrings.
	 * This is accomplished by using a for loop which tracks the value of k iterating between the indices of i and j.
	 * This process continues recursively and the ending result is the maximum number of base pairs.
	 * @param rna
	 * @param i (starting at first character)
	 * @param j (starting at last character)
	 * @return int rnaResult (maximum number of base pairs)
	 */
	private static int slowRNAScore(String rna, int i, int j) {
		if (i <= j) { // Base case, makes sure that indexes don't overlap and are not equal
			// Case 1: Ignore the first character in rna
			int case1 = slowRNAScore(rna, i+1, j);
			// Case 2: Ignore the last character in rna
			int case2 = slowRNAScore(rna, i, j-1);
			int rnaResult = Math.max(case1, case2);
			// Case 3: If the characters can be a pair using the helper function,
			// the first character and last character are ignored. The number
			// of matching base pairs is incremented
			if (basesCanPair(rna.charAt(i), rna.charAt(j))) {
				int case3 = 1 + slowRNAScore(rna, i+1, j-1);
				rnaResult = Math.max(rnaResult, case3);
			}
			// Case 4: Consider ways that rna can be split into two substrings.
			// Find the maximum value of each successive substring represented by k
			for (int k=i; k+1 < j; k++) {
				int case4 = slowRNAScore(rna, i, k) + slowRNAScore(rna, k+1, j);
				rnaResult = Math.max(rnaResult, case4);
			}
			return rnaResult;
		}
		return 0;
	}
					
	
	/**
	 * Helper method that returns true if the bases can pair. A pairs with U and
	 * G pairs with C.
	 * 
	 * @param b1 A RNA base. Must be A, U, G, or C
	 * @param b2 A RNA base. Must be A, U, G, or C
	 * @return true if b1 can pair with b2 and false otherwise
	 */
	private static boolean basesCanPair(char b1, char b2) {
		if (b1 == 'A' && b2 == 'U') {
			return true;
		} else if (b1 == 'G' && b2 == 'C') {
			return true;
		} else if (b1 == 'U' && b2 == 'A') {
			return true;
		} else if (b1 == 'C' && b2 == 'G') {
			return true;
		}
		return false; // bases don't pair return false
	}
	
	/**
	 * This method should return the exact same answers as your
	 * slowRNAScore method. However, it must use memoization to
	 * perform faster. This means that previously calculated
	 * results will be saved for reuse later. You should save
	 * your previously calculated results in a 2D array that
	 * is passed between all recursive method calls. Otherwise,
	 * the code is very similar to the slowRNAScore method.
	 * 
	 * Note: This method is simply the kick-off
	 *       method for your actual recursive algorithm. 
	 *       Also, do NOT use the substring method.
	 * 
	 * @param rna String representing RNA strand. Can only contain
	 *            the letters A, U, C, and G.
	 * @return Maximum number of base pairings in folded strand
	 *         with no pseudo-knots.
	 */
	public static int fastRNAScore(String rna) {
		// known[i][j] will store results of fastRNAScore(known, rna, i, j)
		int[][] known = new int[rna.length()][rna.length()];
		for(int i = 0; i < known.length; i++) {
			Arrays.fill(known[i], -1);
		}
		return fastRNAScore(known, rna.toUpperCase(), 0, rna.length() - 1);
	}

	/**
	 * This method uses the same functionality as the slowRNAScore. However, this method
	 * uses memoization by storing all results in a 2D array representing the previously calculated
	 * results. This 2D array is passed between each recursive call and is checked to see if
	 * a result is already stored in the array. If so, the result is immediately returned and no
	 * more recursive calls are called. This saves time because the method does not have to calculate
	 * results that have already been previously calculated.
	 * @param known
	 * @param rna
	 * @param i (starting at first character)
	 * @param j (starting at last character)
	 * @return int rnaResult (maximum number of base pairs)
	 */
	private static int fastRNAScore(int[][] known, String rna, int i, int j) {
		if (i <= j) { // Base case, makes sure that indexes don't overlap and are not equal
			if (known[i][j] != -1) { // if the calculation is stored in the 2D array, simply return the result
				return known[i][j];
			}
			// Case 1: Ignore the first character in rna
			int case1 = fastRNAScore(known, rna, i+1, j);
			// Case 2: Ignore the last character in rna
			int case2 = fastRNAScore(known, rna, i, j-1);
			int rnaResult = Math.max(case1, case2);
			// Case 3: If the characters can be a pair using the helper function,
			// the first character and last character are ignored. The number
			// of matching base pairs is incremented
			if (basesCanPair(rna.charAt(i), rna.charAt(j))) {
				int case3 = 1 + fastRNAScore(known, rna, i+1, j-1);
				rnaResult = Math.max(rnaResult, case3);
			}
			// Case 4: Consider ways that rna can be split into two substrings.
			// Find the maximum value of each successive substring represented by k
			for (int k=i; k+1 < j; k++) {
				int case4 = fastRNAScore(known, rna, i, k) + fastRNAScore(known, rna, k+1, j);
				rnaResult = Math.max(rnaResult, case4);
			}
			known[i][j] = rnaResult; // store calculated result in known for memoization
			return rnaResult;
		}
		return 0;
	}
}
