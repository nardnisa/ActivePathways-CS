package org.jurilab.apw.internal.demo5;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.distribution.*;
import org.apache.commons.math3.linear.RealMatrix;

public class StatCalculator {
	private void Hypergeomatric(Integer[] counts){
		int m=counts[0]+counts[1];
		int n=counts[2]+counts[3];
		int k=counts[0]+counts[2];
		int x=counts[0];
		int N = m+n+k+x;
		int populationSize=N;
		int numberOfSuccesses=x-1;
		int sampleSize=k;
		
		if(counts.length<0){
			System.out.print("counts contains negative values. Something went very wrong.");	
		}else{
			//System.out.println(m+" "+n+" "+k+" "+x);
		}
	}
	
	public static strictfp double getSDpopulation(double[] values) {
		double mean = getMean(values);
		double n = values.length;
		double dv = 0;
		for (double d : values) {
			double dm = d - mean;
			dv += dm * dm;
		}
		return Math.sqrt(dv / n);
	}
	
	
	public static strictfp double getMean(double[] values) {
		return getSum(values) / values.length;
	}
	
	public static strictfp double getSum(double[] values) {
		if (values == null || values.length == 0) {
			throw new IllegalArgumentException("The data array either is null or does not contain any data.");
		}else{
			double sum = 0;
			for (int i = 0; i < values.length; i++) {
				sum += values[i];
			}
			return sum;
		}
	}
	
	public static double getMax(double[] values) {
		double maximum = values[0];
		try {
			for (double d: values) {
		        maximum = Math.max(maximum, d);
		    }	
		} catch (Exception e) {
			System.out.println("Remove NA");
		}
	    return maximum;
	}
	 
	public static double getMin(double[] values) {
		double minimum = values[0];
		try {
		    for (double d: values) {
		    	minimum = Math.min(minimum, d);
		    }	
		} catch (Exception e) {
			System.out.println("Remove NA");
		}
		
	    return minimum;
	} 
	
	public static double getMedian(double[] values){
		Arrays.sort(values);
		double median;
		if (values.length % 2 == 0)
		    median = ((double)values[values.length/2] + (double)values[values.length/2 - 1])/2;
		else
		    median = (double) values[values.length/2];
		return median;
	}
	
	public static double[][] makeMatrixFromFile(String fin, Integer startrow, Integer startcol){
		int nrow=countRow(fin, startrow);
		int ncol=countCol(fin, startcol)-startcol;
		int ln =0;
		double[][] matrix = new double[nrow][ncol];
		int x=0, y=0;
		try {
            FileInputStream fi 	= new FileInputStream(fin);
  			BufferedReader 	br 	= new BufferedReader(new InputStreamReader(fi));  			
  			String sl = new String();
  			while ((sl=br.readLine())!=null){
  				if(ln>=startrow){
  					String sp[]=sl.split("\t");
  					for (int i = startcol; i < sp.length; i++) {
  						double dat= Double.parseDouble(sp[i]);
  						matrix[x][y]=dat;
  						y++;
					}
  					x++;
  					y=0;
  					ln++;
  				}else{
  					ln++;
  				}
  			}
  			br.close();
  			fi.close();
		}catch (Exception e) {
			System.out.println("Error"+e);
		}
		return matrix;
	}
	
	
	public static Integer countRow(String fin,Integer startrow){
		int nrow=0;
		int ln=0;
		try {
            FileInputStream fi 	= new FileInputStream(fin);
  			BufferedReader 	br 	= new BufferedReader(new InputStreamReader(fi));  			
  			String sl = new String();
  			while ((sl=br.readLine())!=null){
  				if(ln>=startrow){
  					nrow++;
  				}else{
  					ln++;
  				}
  			}
  			br.close();
  			fi.close();
		}catch (Exception e) {
			System.out.println("Error"+e);
		}
		return nrow;
	}
	
	public static Integer countCol(String fin,Integer startcol){
		int ncol=0;
		
		try {
            FileInputStream fi 	= new FileInputStream(fin);
  			BufferedReader 	br 	= new BufferedReader(new InputStreamReader(fi));  			
  			String sl = new String();
  			while ((sl=br.readLine())!=null){
  				ncol=sl.split("\t").length-startcol+1;
  				break;
  			}
  			br.close();
  			fi.close();
		}catch (Exception e) {
			System.out.println("Error"+e);
		}
		return ncol;
	}
	
	public double[][] transpose(double[][]matrix) {
		int M = matrix.length;
		int N = matrix[0].length;
		double[][] tmatrix  = new double[N][M];
		for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			tmatrix[j][i] = matrix[i][j];
		return tmatrix;
	 }
	
	public static double[] removeZeros(double[] dat){
		
		int targetIndex = 0;
		for(int sourceIndex=0;sourceIndex<dat.length;sourceIndex++ )
		{
		    if( dat[sourceIndex]!= 0 )
		        dat[targetIndex++] = dat[sourceIndex];
		}
		double[] newdat = new double[targetIndex];
		System.arraycopy(dat, 0, newdat, 0, targetIndex );
		return newdat;
	}
	
	public static double[] removeOnes(double[] dat){
		int targetIndex = 0;
		for(int sourceIndex=0;sourceIndex<dat.length;sourceIndex++ )
		{
		    if( dat[sourceIndex]!=1.0)
		        dat[targetIndex++]=dat[sourceIndex];
		}
		double[] newdat = new double[targetIndex];
		System.arraycopy(dat,0,newdat,0,targetIndex);
		return newdat;
	}
	
	public Hashtable<Double,Double> EmpiricalDistribution(double observes[]){
		Hashtable<Double,Integer> freqmap = new Hashtable<Double,Integer>();
		Hashtable<Double,Double> cdfmap = new Hashtable<Double,Double>();
	
		for (int i = 0; i < observes.length; i++) {
			if(freqmap.size()>0){
				if(freqmap.containsKey(observes[i])){
					int countn=freqmap.get(observes[i])+1;
					freqmap.put(observes[i],countn);	
				}else{
					freqmap.put(observes[i],1);	
				}
			}else{
				freqmap.put(observes[i],1);	
			}
		}
		double cdf[] = new double[freqmap.size()];
		List<Double> key = new ArrayList<Double>(freqmap.keySet());
	    Collections.sort(key);
	    double ntotal = observes.length;
	    double cumdist=0.0;
	    int count=0;
	    for (Double freq:key) {
	    	double prob=freqmap.get(freq)/ntotal;
	    	cumdist=cumdist+prob;
    		cdfmap.put(freq, cumdist);
    		cdf[count]=cumdist;
    		count++;
	    }
	    if (Math.abs(cdf[cdf.length - 1]-1.0) > 1E-7) {
            throw new IllegalArgumentException("The sum of probabilities is not 1.");
        }
		return cdfmap;
	}
	
	public double[] getColumn(double[][] array, int colindex){
		double[] column = new double[array.length]; //  assuming a rectangular 2D array!
	    for(int i=0; i<column.length; i++){
	    	column[i] = array[i][colindex];
	    }
	    return column;
	}
	
	public double[] getlowertri(RealMatrix matrix, boolean diag){ 
		int nrow = matrix.getRowDimension();
		int ncol= matrix.getColumnDimension();
		ArrayList<Double> entrylist = new ArrayList<Double>();
		int i, j; 
		for (i = 0; i < nrow; i++){ 
			for (j = 0; j < ncol; j++){ 
			  if(diag){
				  if (i < j){ 
			          //System.out.print("0" + " ");
			      }else{
			    	  //System.out.print(matrix[i][j] + " ");
			    	  entrylist.add(matrix.getEntry(i,j));
			      }
			  }else{
			      if (i < j){ 
			         // System.out.print("0" + " ");
			      }else if(i==j){
			    	 // System.out.print("0" + " "); 
			      }else{
			    	  //System.out.print(matrix.getEntry(i,j)+ " ");
			    	  entrylist.add(matrix.getEntry(i,j));
			      }
			  }
		  }
		}
		Double entry[] = new Double[entrylist.size()];
		double entryfinal[] = ArrayUtils.toPrimitive(entrylist.toArray(entry));
		return entryfinal;	
	} 

	
	// Method to form upper 
	// triangular matrix 
	public static void uppertri(int matrix[][],int row, int col,boolean diag){ 
		int i, j; 
		for (i = 0; i < row; i++){ 
		  for (j = 0; j < col; j++){
			  if(diag){
				  if (i > j){ 
			          //System.out.print("0" + " "); 
			      }else{
			    	  //System.out.print(matrix[i][j] + " ");
			      }
			  }else{
				  if (i > j){ 
			          //System.out.print("0" + " "); 
			      }else if(i==j){
			    	  //System.out.print("0" + " "); 
			      }else{
			    	  //System.out.print(matrix[i][j] + " ");
			      }  
			  }
		  } 
		} 
	} 
	
	public double chiSquare(double[] expected, long[] observed)throws IllegalArgumentException {
	/**
     * {@inheritDoc}
     * <p><strong>Note: </strong>This implementation rescales the 
     * <code>expected</code> array if necessary to ensure that the sum of the
     * expected and observed counts are equal.</p>
     * 
     * @param observed array of observed frequency counts
     * @param expected array of expected frequency counts
     * @return chi-square test statistic
     * @throws IllegalArgumentException if preconditions are not met
     * or length is less than 2
     */
	    if ((expected.length < 2) || (expected.length != observed.length)) {
            throw new IllegalArgumentException("observed, expected array lengths incorrect");
        }
        if (!isPositive(expected) || !isNonNegative(observed)) {
            throw new IllegalArgumentException("observed counts must be non-negative and expected counts must be postive");
        }
        double sumExpected = 0d;
        double sumObserved = 0d;
        for (int i = 0; i < observed.length; i++) {
            sumExpected += expected[i];
            sumObserved += observed[i];
        }
        double ratio = 1.0d;
        boolean rescale = false;
        if (Math.abs(sumExpected - sumObserved) > 10E-6) {
            ratio = sumObserved / sumExpected;
            rescale = true;
        }
        double sumSq = 0.0d;
        double dev = 0.0d;
        for (int i = 0; i < observed.length; i++) {
            if (rescale) {
                dev = ((double) observed[i] - ratio * expected[i]);
                sumSq += dev * dev / (ratio * expected[i]);
            } else {
                dev = ((double) observed[i] - expected[i]);
                sumSq += dev * dev / expected[i];
            }
        }
        return sumSq;
    }	
	
	
	 public double chiSquareTest(double[] expected, long[] observed, ChiSquaredDistribution distribution)throws IllegalArgumentException {
	 /**
	  * {@inheritDoc}
	  * <p><strong>Note: </strong>This implementation rescales the 
	  * <code>expected</code> array if necessary to ensure that the sum of the
	  * expected and observed counts are equal.</p>
	  * 
	  * @param observed array of observed frequency counts
	  * @param expected array of expected frequency counts
	  * @return p-value
	  * @throws IllegalArgumentException if preconditions are not met
	  * @throws MathException if an error occurs computing the p-value
	  */
	  distribution= new ChiSquaredDistribution(expected.length - 1.0);
	  return 1.0 - distribution.cumulativeProbability(chiSquare(expected, observed));
    }
	 
	 private boolean isNonNegative(long[] in) {
		 for (int i = 0; i < in.length; i ++) {
			 if (in[i] < 0) {
				 return false;
			 }
		 }
		 return true;
	 }
	 
	 private boolean isPositive(double[] in) {
		 for (int i = 0; i < in.length; i ++) {
			 if (in[i] <= 0) {
				 return false;
			 }
		 }
		 return true;
	 }
}
