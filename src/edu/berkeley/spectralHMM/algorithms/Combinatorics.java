 /*
    This file is part of spectralHMM.

    spectralHMM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    spectralHMM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with spectralHMM.  If not, see <http://www.gnu.org/licenses/>.
  */

package edu.berkeley.spectralHMM.algorithms;

import java.math.BigInteger;
import java.util.HashMap;

import edu.berkeley.spectralHMM.datatypes.Pair;

// Contains some combinatoric functions
public class Combinatorics {
	private static HashMap<Pair<Integer, Integer>, BigInteger> chooseCache = new HashMap<Pair<Integer, Integer>, BigInteger>();
	private static HashMap<Integer, BigInteger> factorialCache = new HashMap<Integer, BigInteger>();
	
	public static BigInteger choose(int n, int k){	//returns n choose k
		if (n < 0 || k < 0 || n < k)	return BigInteger.ZERO;
		if (n == 0 || k == 0) return BigInteger.ONE;
		Pair<Integer, Integer> pr = new Pair<Integer, Integer>(n, k);
		if (chooseCache.containsKey(pr)) {
			return chooseCache.get(pr);
		}
		
		BigInteger ret = choose(n - 1, k - 1).add(choose (n - 1, k));
		chooseCache.put(pr, ret);
		return ret;
	}
	
	public static BigInteger factorial(int n){
		if (n < 0)	return BigInteger.ZERO;
		if (n == 0 || n == 1)	return BigInteger.ONE;
		if (factorialCache.containsKey(n))	return factorialCache.get(n);
		BigInteger ret = new BigInteger("" + n);
		ret = ret.multiply(factorial(n-1));
		factorialCache.put(n, ret);
		return ret;
	}

	public static double logChoose(int n, int k) {
		if (k > n - k)	k = n - k;
		double ret = 0.;
		int j = k;
		for (int i = n; i >= n - k + 1; i--, j--)	{
			ret = ret + Math.log(i) - Math.log(j);
		}
		return ret;
	}
}
