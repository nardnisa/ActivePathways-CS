package org.jurilab.apw.internal.demo5;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.Map;
 
public class PValueCorrection {
	
	public static double[] pAdjust(ArrayList<Double> pvaluelist, String correction_method) {
		double[] pvalues = new double[pvaluelist.size()];
		for (int i = 0; i < pvaluelist.size(); i++) {
			pvalues[i]=pvaluelist.get(i);
		}
		int size = pvalues.length;
		if (size < 1) throw new IllegalArgumentException("pAdjust requires at least one element");
		int type;
		switch (correction_method.toLowerCase()) {
			case "bh":
		    case "fdr": type=0; break;
		    case "by": type=1; break;
		    case "bonferroni": type=2; break;
		    case "hochberg": type=3; break;
		    case "holm": type=4; break;
		    case "hommel": type=5; break;
		    default: throw new IllegalArgumentException(correction_method + " doesn't match any accepted FDR types");
		 }
	    
		if (type==2) {  // Bonferroni method
			double[] result = new double[size];
			for (int i = 0; i < size; i++) {
				double b = pvalues[i]*size;
				if (b >= 1) {
					result[i] = 1;
				}else if (0 <= b && b < 1) {
					result[i] = b;
				}else{
					throw new RuntimeException(b+" is outside [0, 1)");
				}
			}
	        return result;
	        
		}else if (type==4) {  // Holm method
			int[] o = Order(pvalues, false);
			double[] o2Double = intToDouble(o);
			double[] cummaxInput = new double[size];
			for (int i = 0; i < size; i++){
				cummaxInput[i] = (size-i)*pvalues[o[i]];
			}
			int[] ro = Order(o2Double, false);
			double[] cummaxOutput = cummax(cummaxInput);
	        double[] pmin = pminx(cummaxOutput,1.0);
	        double[] result = new double[size];
	        for (int i = 0; i < size; i++) {
	        	result[i] = pmin[ro[i]];
	        }
	        return result;
		}else if (type==5) {
			int[] indices=seqLen(size, size);
			int[] o = Order(pvalues, false);
			double[] p = new double[size];
			for (int i = 0; i < size; i++) {
				p[i] = pvalues[o[i]];
			}
			double[] o2Double = intToDouble(o);
			int[] ro = Order(o2Double, false);
			double[] q = new double[size];
			double[] pa = new double[size];
			double[] npi = new double[size];
			for (int i = 0; i < size; i++) {
				npi[i] = p[i]*size/indices[i];
			}
			double min = doubleArrayMin(npi);
			Arrays.fill(q, min);
			Arrays.fill(pa, min);
			for (int j = size; j>=2; --j) {
				int[] ij = seqLen(1, size-j+1);
				for (int i = 0; i < size-j+1; i++) {
                    ij[i]--;
				}
	            int i2Length = j-1;
	            int[] i2 = new int[i2Length];
	            for (int i = 0; i < i2Length; i++) {
	            	i2[i] = size-j+2+i-1;
	            }
	            double q1 = j * p[i2[0]]/2.0;
	            for (int i = 1; i < i2Length; i++) {
	            	double temp_q1 = p[i2[i]]*j/(2.0+i);
	            	if (temp_q1 < q1) q1=temp_q1;
	            }
	            for (int i = 0; i < size-j+1; i++) {
	            	q[ij[i]] = Math.min(p[ij[i]]*j,q1);
	            }
	            for (int i = 0; i < i2Length; i++) {
	            	q[i2[i]] = q[size-j];
	            }
	            for (int i = 0; i < size; i++) {
	            	if (pa[i] < q[i]) {
	            		pa[i] = q[i];
	            	}
	            }
			}
			for (int i = 0; i < size; i++) {
				q[i] = pa[ro[i]];
			}
	        return q;
		}
		
		double[] ni = new double[size];
		int[] o = Order(pvalues, true);
		double[] oDouble = intToDouble(o);
		for (int i = 0; i < size; i++) {
			if (pvalues[i]<0||pvalues[i]>1){
				//System.out.println("Warning: array[" + i + "] = "+pvalues[i]+" is outside [0, 1]");
			}
			ni[i] = (double)size/(size-i);
		}
		int[] ro = Order(oDouble, false);
		double[] cumminInput = new double[size];
		if (type==0) {// BH method
			for (int i = 0; i < size; i++) {
				cumminInput[i] = ni[i] * pvalues[o[i]];
			}
		}else if(type==1){  // BY method
			double q = 0;
			for (int i = 1; i < size+1; i++) {
				q += 1.0/i;
			}
	        for (int i = 0; i < size; i++) {
	        	cumminInput[i] = q*ni[i]*pvalues[o[i]];
	        }
		}else if (type==3) { // Hochberg method
			for (int i = 0; i < size; i++) {
				cumminInput[i] = (i+1)*pvalues[o[i]];
			}
		}
		double[] cumminArray = cummin(cumminInput);
	    double[] pmin = pminx(cumminArray,1.0);
	    double[] result = new double[size];
	    for (int i = 0; i < size; i++) {
	    	result[i] = pmin[ro[i]];
	    }
	    return result;
	}
	 
	private static int[] seqLen(int start, int end) {
        int[] result;
        if (start == end) {
            result = new int[end+1];
            for (int i = 0; i < result.length; i++) {
                result[i] = i+1;
            }
        } else if (start < end) {
            result = new int[end-start+1];
            for (int i = 0; i < result.length; i++) {
                result[i] = start+i;
            }
        } else {
            result = new int[start-end+1];
            for (int i = 0; i < result.length; i++) {
                result[i] = start-i;
            }
        }
        return result;
    }
 
    private static int[] Order(double[] array, boolean decreasing) {
    	ArrayList<Map.Entry<Integer, Double>> sortedarray = new ArrayList<Map.Entry<Integer,Double>>();
		Hashtable<Integer,Double> hsorted = new Hashtable<Integer, Double>();
    	int[] idx = new int[array.length];
    	for (int i = 0; i < array.length; i++) {
        	hsorted.put(i,array[i]);
        }
    	sortedarray=sortHashMapValue(hsorted,decreasing);
    	for (int i = 0; i < sortedarray.size(); i++) {
    		idx[i]=sortedarray.get(i).getKey();
    	}
    	return idx;
    }
    
    private static double[] cummin(double[] array) {
        if (array.length < 1) throw new IllegalArgumentException("cumulative min. value requires at least one element");
        double[] output = new double[array.length];
        double cumulativeMin = array[0];
        for (int i = 0; i < array.length; i++) {
            if (array[i] < cumulativeMin) cumulativeMin = array[i];
            output[i] = cumulativeMin;
        }
        return output;
    }
 
    private static double[] cummax(double[] array) {
        if (array.length < 1) throw new IllegalArgumentException("cumulative max. value requires at least one element");
        double[] output = new double[array.length];
        double cumulativeMax = array[0];
        for (int i = 0; i < array.length; i++) {
            if (array[i] > cumulativeMax) cumulativeMax = array[i];
            output[i] = cumulativeMax;
        }
        return output;
    }
 
    private static double[] pminx(double[] array, double x) {
        if (array.length < 1) throw new IllegalArgumentException("pmin requires at least one element");
        double[] result = new double[array.length];
        for (int i = 0; i < array.length; i++) {
            if (array[i] < x) {
                result[i] = array[i];
            } else {
                result[i] = x;
            }
        }
        return result;
    }
    
    private static double[] intToDouble(int[] array) {
        double[] result = new double[array.length];
        for (int i = 0; i < array.length; i++) {
            result[i] = array[i];
        }
        return result;
    }
 
    private static double doubleArrayMin(double[] array) {
        if (array.length < 1) throw new IllegalArgumentException("pAdjust requires at least one element");
        return Arrays.stream(array).min().orElse(Double.NaN);
    }
    
    private static ArrayList<Map.Entry<Integer, Double>>  sortHashMapValue(Hashtable<Integer,Double> array, final boolean decreasing){
    	ArrayList<Map.Entry<Integer, Double>> l = new ArrayList(array.entrySet());
    	Collections.sort(l, new Comparator<Map.Entry<Integer, Double>>(){
    		public int compare(Map.Entry<Integer, Double> o1, Map.Entry<Integer, Double> o2) {
    			if(decreasing){
    				return o2.getValue().compareTo(o1.getValue()); //max2min  
    			}else{
    	    	   return o1.getValue().compareTo(o2.getValue()); //min2max
    			}
    		}
    	});
    	return l;
	}
} 
