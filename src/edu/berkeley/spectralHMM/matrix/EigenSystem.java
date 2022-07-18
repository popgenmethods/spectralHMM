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

package edu.berkeley.spectralHMM.matrix;

import java.math.BigDecimal;
import java.util.Arrays;

import edu.berkeley.spectralHMM.datatypes.Pair;

public class EigenSystem {
	// variables
	public final BigDecimal[] eigenValues;
	public final BigDecimal[][] eigenVectors;
	
	public EigenSystem(BigDecimal[][] eigenVectors, BigDecimal[] eigenValues)	{
		int dim = eigenValues.length;
		
		@SuppressWarnings("unchecked")
		Pair<BigDecimal, Integer>[] prs = new Pair [dim];
		
		//sort the eigenvalues in increasing order
		for (int i = 0; i < dim; i++)	{
			prs[i] = new Pair<BigDecimal, Integer>(eigenValues[i], i);
		}
		// sort them
		Arrays.sort(prs);
		
		this.eigenVectors = new BigDecimal[dim][];
		this.eigenValues = new BigDecimal[dim];
		
		for (int i = 0; i < dim; i++)	{
			this.eigenVectors[i] = eigenVectors[prs[i].second()];
			this.eigenValues[i] = eigenValues[prs[i].second()];
		}
	}
}