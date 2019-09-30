package org.jurilab.apw.internal.demo5;
import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Random;

import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.stat.correlation.Covariance;


public class RunActivePathways  {
	private IOhandling io = new IOhandling();
	private static StatCalculator stat = new StatCalculator();
	private PValueCorrection pcorrect = new PValueCorrection();
	private Hashtable<String, String> hgmt_names= new Hashtable<String, String>();
	private Hashtable<String, String> hgmt_genelist= new Hashtable<String, String>();
	private ArrayList<String> genelistUpdate = new ArrayList<String>();
	public int TIME_VISIBLE = 30000;  //30 seconds
	private int ndigitLimit=20;
	private int progress = 0;
	
	private int repaintBar(JPanel progressPane, JProgressBar progressBar,int progress, int pmax) {
		try {
			progressPane.revalidate();
			progressPane.repaint();
			Random random = new Random();
			int diff=Math.abs(pmax-progress);
			if(diff==0){
				diff++;
			}
			progress+= random.nextInt(diff);
			if(progress==pmax) {
				progress=pmax--;
			}
			progress=Math.min(progress, pmax);
			progressBar.setValue(progress);	
		} catch (Exception e) {
		}
		return progress;
	}
	public Hashtable<String, String> mainmodule(String fscores, String fgmt, Hashtable<String, String> hgmt_names,Hashtable<String, String> hgmt_genelist, double gsetmin, double gsetmax, String mergemethod, String correctmethod,String returnall, String cytofiles,double cutoff,double significant,ArrayList<Integer> genesetfilter,JPanel progressPane, JProgressBar progressBar){
		if(progress>=100) {
			progress=0;
		}
		progress=repaintBar(progressPane,progressBar,progress,5);
		
		Hashtable<String,Double> mergedscores= new Hashtable<String,Double>();
		Hashtable<String,String> enrichresults= new Hashtable<String,String>();
		Hashtable<String,String> enrichresults_sig = new Hashtable<String,String>();
		Hashtable<String,String> finalenrichresults= new Hashtable<String,String>();
		Hashtable<String,String> hgmt_names_sig = new Hashtable<String,String>();
		Hashtable<String,String> hgmt_genelist_sig = new Hashtable<String,String>();
		ArrayList<Map.Entry<String, Double>> sortedmergedscores = new ArrayList<Map.Entry<String,Double>>();
		ArrayList<String> orderedgenelist = new ArrayList<String>();
		ArrayList<Double> sigpbrown = new ArrayList<Double>();
		ArrayList<Double> enrpval = new ArrayList<Double>();
		ArrayList<String> GOIDlist= new ArrayList<String>();
		Hashtable<String, Boolean> sigindx = new Hashtable<String, Boolean>();
		Hashtable<String, String> sigcols = new Hashtable<String, String>();
		ArrayList<String> background = new ArrayList<String>();
		boolean checkValidFormat = validate_scores(fscores);
		boolean checkValidSigIndex=true;
		boolean checkVector = isVector(fscores);
		
		progress=repaintBar(progressPane,progressBar,progress,10);
		if(checkValidFormat){
			double[][] matrix=stat.makeMatrixFromFile(fscores,1,1);
			RealMatrix datamatrix = MatrixUtils.createRealMatrix(matrix);
			ArrayList<String> genelist = getGenelist(fscores,1);
			readGMT(fgmt,hgmt_names,hgmt_genelist,genesetfilter);
			
			progress=repaintBar(progressPane,progressBar,progress,15);
			background = makeBackground(fgmt, background);
			
			progress=repaintBar(progressPane,progressBar,progress,20);
			double[][] newmatrix=remove_gene_notinbackground(datamatrix, genelist,background,genelistUpdate);
			
			progress=repaintBar(progressPane,progressBar,progress,25);
			mergedscores= merge_pval(newmatrix,genelistUpdate, mergemethod,checkVector);
			
			progress=repaintBar(progressPane,progressBar,progress,30);
			if(mergedscores.size()==0){
				JOptionPane.showMessageDialog(null,"No genes made the cutoff","Warning",JOptionPane.ERROR_MESSAGE);
				checkValidSigIndex=false;
			}else{
				// Sort genes by p-value
				progress=repaintBar(progressPane,progressBar,progress,35);
				sortedmergedscores =sortHashMapValue(mergedscores,false); //min to max (most to least sig)
				for (int i = 0; i < sortedmergedscores.size(); i++) {
					if(sortedmergedscores.get(i).getValue()<=cutoff){
						orderedgenelist.add(sortedmergedscores.get(i).getKey());
						sigpbrown.add(sortedmergedscores.get(i).getValue());
					}
				}
			}
			progress=repaintBar(progressPane,progressBar,progress,40);
			if(checkValidSigIndex){
				enrichresults=enrichmentAnalysis(orderedgenelist,fgmt,"",hgmt_names,hgmt_genelist,background);
				Hashtable<String, Double> hadjpval = new Hashtable<String, Double>(); 
				ArrayList<String> newInfo=new ArrayList<String>();
				Enumeration e = enrichresults.keys();
				while (e.hasMoreElements()) {
					Object GOID = e.nextElement();
					GOIDlist.add(GOID.toString());
					String info[]=enrichresults.get(GOID.toString()).split("\t");
					double pval=Double.parseDouble(info[1]);
					enrpval.add(pval);
					StringBuffer retext = new StringBuffer();
					retext.append(info[0]);
					for (int i = 1; i < info.length; i++) {
						if(i==1) {
							retext.append("\ttmppvalue");
						}else{
							retext.append("\t"+info[i]);
						}
					}
					newInfo.add(retext.toString());
				}
				
				progress=repaintBar(progressPane,progressBar,progress,45);
				double adjpval[] = pcorrect.pAdjust(enrpval,correctmethod);
				for (int i = 0; i < adjpval.length; i++) {
					enrichresults.put(GOIDlist.get(i),newInfo.get(i).replace("tmppvalue", Double.toString(adjpval[i])));
				}
				
				/******************************/
				// extract only sig adjpval
				/******************************/
				progress=repaintBar(progressPane,progressBar,progress,55);
				for (int i = 0; i < adjpval.length; i++) {
					if(adjpval[i]<=significant){
						sigindx.put(GOIDlist.get(i),true);
					}
				}
				
				if (sigindx.size() == 0) {
					checkValidSigIndex=false;
					if(returnall.toLowerCase().equals("no")){
						JOptionPane.showMessageDialog(null,"No significant terms were found\r\nCytoscape files were not written","Warning",JOptionPane.ERROR_MESSAGE);
					}
					if(cytofiles.isEmpty()){
						JOptionPane.showMessageDialog(null,"Cytoscape files were not found and written","Warning",JOptionPane.ERROR_MESSAGE);
					}
				}
				
				progress=repaintBar(progressPane,progressBar,progress,65);
				if(checkValidSigIndex) {
					/******************************/
					// Column contribution
					/******************************/
					boolean contribution = check_contribution(fscores);
					boolean checksigcols=true;
					if(contribution){
						progress=repaintBar(progressPane,progressBar,progress,75);
						sigcols = columnSignificance(fscores,fgmt,datamatrix,hgmt_names,hgmt_genelist,background,cutoff,significant,correctmethod,adjpval);
						Enumeration e2 = enrichresults.keys();
						while (e2.hasMoreElements()) {
							Object GOID = e2.nextElement();
							String resinfo = enrichresults.get(GOID.toString());
							String sigcolinfo=sigcols.get(GOID.toString());
							enrichresults.put(GOID.toString(),resinfo+"\t"+sigcolinfo);
						}
					}else{
						checksigcols=false;
					}
					/******************************/
					// Prepare Cytoscape file for only Significant enrichment pvalues 
					/******************************/
					// Get onyl sig.pvalue
					progress=repaintBar(progressPane,progressBar,progress,80);
					Enumeration e3 = enrichresults.keys();
					while (e3.hasMoreElements()) {
						Object GOID = e3.nextElement();
						if(sigindx.containsKey(GOID.toString())){
							if(sigindx.get(GOID.toString())==true){
								String info[]= enrichresults.get(GOID).toString().split("\t");
								String termname = info[0];
								String adjp = info[1];
								enrichresults_sig.put(GOID.toString(),termname+"\t"+adjp);
								hgmt_names_sig.put(GOID.toString(),hgmt_names.get(GOID.toString()));
								hgmt_genelist_sig.put(GOID.toString(),hgmt_genelist.get(GOID.toString()));
							}	
						}
						
					}
					progress=repaintBar(progressPane,progressBar,progress,85);
					if(returnall.toLowerCase().contains("yes")){
						finalenrichresults=enrichresults;
						System.out.println("Final Enrichment size: "+enrichresults.size());
					}else{
						System.out.println("Return only sig. enrichment analysis results");
						finalenrichresults=enrichresults_sig;
						System.out.println("Final EnrichmentSig size: "+enrichresults_sig.size());
					}
					progress=repaintBar(progressPane,progressBar,progress,95);
					if(!cytofiles.isEmpty()&&sigindx.size()>0){
						prepareCytoscape(finalenrichresults,hgmt_names_sig,hgmt_genelist_sig,sigcols,cytofiles);
					}
				}
			}
		}
		progress=repaintBar(progressPane,progressBar,progress,100);
		return finalenrichresults;
	}	
	
	private boolean validate_scores(String fscores){
		/* check:
		 * is.matrix
		 * is.numeric
		 * is.na
		 * is.scores<0 || is.scores>1
		 * is.anyduplicate
		*/
		boolean checkAllValid=true;
		boolean ismatrix=check_isMatrix(fscores);
		boolean isnumeric=check_isNumeric(fscores);
		boolean isNA=check_isNA(fscores);
		boolean isvalidscores=check_isValidScores(fscores);
		boolean isduplicated=check_isDuplicatedGene(fscores);
		
		StringBuffer validateMsg = new StringBuffer();
		validateMsg.append("Please check your input files:\r\n");
		if(!ismatrix){ System.out.println("Scores must be a matrix"); validateMsg.append("- Scores must be a matrix\r\n");}
		if(isNA){System.out.println("Scores may not contain missing values");validateMsg.append("- Scores may not contain missing values\r\n");}
		if(!isnumeric){ System.out.println("Scores must be numeric");validateMsg.append("- Scores must be numeric\r\n"); }
		if(!isvalidscores){ System.out.println("All values in scores must be in [0,1]");validateMsg.append("- All values in scores must be in [0,1]\r\n"); }
		if(isduplicated){System.out.println("Genes cannot be duplicated");validateMsg.append("- Genes cannot be duplicated\r\n"); }
		
		if(!ismatrix ||isNA|| !isnumeric || !isvalidscores || isduplicated) {
			JOptionPane.showMessageDialog(null,validateMsg.toString(),"Warning",JOptionPane.ERROR_MESSAGE);
			checkAllValid=false;
		}
		return checkAllValid;
	}
	
	private Boolean check_isMatrix(String fscore){
		boolean check_equalncol=true;
		int firstncol = 0;
		int firstln=0;
		try {
            FileInputStream fi 	= new FileInputStream(fscore);
  			BufferedReader 	br 	= new BufferedReader(new InputStreamReader(fi));  			
  			String sl = new String();
  			while ((sl=br.readLine())!=null){
  				int ncol=sl.split("\t").length;
  				if(firstln>0){
  					if(ncol!=firstncol){
  						check_equalncol=false;
  						break;
  					}  					
  				}else{
  					firstncol=ncol;
  					firstln=1;
  				}
  			}
  			br.close();
  			fi.close();
		}catch (Exception e) {
			System.out.println("Error"+e);
		}
		return check_equalncol;
	}
	private Boolean check_isNumeric(String fscore){
		int firstln=0;
		boolean check_isnumeric=true;
		try {
            FileInputStream fi 	= new FileInputStream(fscore);
  			BufferedReader 	br 	= new BufferedReader(new InputStreamReader(fi));  			
  			String sl = new String();
  			while ((sl=br.readLine())!=null){
  				if(firstln>0){
  					String sp[]=sl.split("\t");
  					for (int i = 1; i < sp.length; i++) {
  						if(!isNumeric(sp[i])){
  		  	 				check_isnumeric=false;
  		  	  				break;
  		  	  			}	
					}
  				}else{
  					firstln=1;
  				}
  			}
  			br.close();
  			fi.close();
		}catch (Exception e) {
			System.out.println("Error"+e);
		}
		
		return check_isnumeric;
	}
	private Boolean check_isNA(String fscore){
		int firstln=0;
		boolean check_isNA=false;
		try {
            FileInputStream fi 	= new FileInputStream(fscore);
  			BufferedReader 	br 	= new BufferedReader(new InputStreamReader(fi));  			
  			String sl = new String();
  			while ((sl=br.readLine())!=null){
  				if(firstln>0){
  					String gname=sl.split("\t")[0];
  					String values=sl.replace(gname+"\t","");
  					if(values.toUpperCase().contains("\tNA")){
  	  					check_isNA=true;
  	  					break;
  	  				}	
  				}else{
  					firstln=1;
  				}
  			}
  			br.close();
  			fi.close();
		}catch (Exception e) {
			System.out.println("Error"+e);
		}
		return check_isNA;
	}
	private Boolean check_isValidScores(String fscore){
		int firstln=0;
		boolean check_validscores=true;
		try {
            FileInputStream fi 	= new FileInputStream(fscore);
  			BufferedReader 	br 	= new BufferedReader(new InputStreamReader(fi));  			
  			String sl = new String();
  			while ((sl=br.readLine())!=null && check_validscores){
  				if(firstln>0){
  					String sp[]=sl.split("\t");
  	  				for (int i = 1; i < sp.length; i++) {
  	  					if(isNumeric(sp[i])){
  	  						try {
	  	  						double df = Double.parseDouble(sp[i]);
	  	  						if(df<0||df>1){
	  	  							check_validscores=false;
	  	  							break;
	  	  						}	
							} catch (Exception e) {
								check_validscores=false;
							}
  	  					}else{
	  	  					check_validscores=false;
  	  					}
  					}	
  				}else{
  					firstln=1;
  				}
  			}
  			br.close();
  			fi.close();
		}catch (Exception e) {
			System.out.println("Error"+e);
		}
		
		return check_validscores;
	}
	private Boolean check_isDuplicatedGene(String fscore){
		ArrayList<String> glist = new ArrayList<String>();
		int firstln=0;
		boolean check_duplicated=false;
		try {
            FileInputStream fi 	= new FileInputStream(fscore);
  			BufferedReader 	br 	= new BufferedReader(new InputStreamReader(fi));  			
  			String sl = new String();
  			while ((sl=br.readLine())!=null){
  				if(firstln>0){
  					String genename = sl.split("\t")[0];
  					if(glist.size()>0){
  						if(glist.contains(genename)){
  							System.out.println(genename);
  							check_duplicated=true;
  							break;
  						}else{
  							glist.add(genename);	
  						}
  					}else{
  						glist.add(genename);
  					}
  				}else{
  					firstln=1;
  				}
  			}
  			br.close();
  			fi.close();
		}catch (Exception e) {
			System.out.println("Error"+e);
		}
		return check_duplicated;
	}
	private void validate_thresholdvalues(String fscore){
		
	}
	private Boolean check_thresholds(String cutoff, String significant){
		boolean checkcutoff=true;
		boolean checksig=true;
		boolean checkthreshold=false;
		
		try {
			double cutoff_value=Double.parseDouble(cutoff);
			if(cutoff_value<0 || cutoff_value>1){
				System.out.println("cutoff must be a value in [0,1]");
				checkcutoff=false;
			}
		} catch (Exception e) {
			System.out.println(e);
		}
		try {
			double sig_value=Double.parseDouble(significant);
			if(sig_value<0 || sig_value>1){
				System.out.println("significant must be a value in [0,1]");
				checksig=false;
			}
		} catch (Exception e) {
			System.out.println(e);
		}
		
		if(checkcutoff && checksig){
			checkthreshold=true;
		}
		return checkthreshold;
	}
	private void readGMT(String fgmt, Hashtable<String, String> hgmt_names,Hashtable<String, String> hgmt_genelist, ArrayList<Integer>genesetfilter){
		int count=0;
		try {
            FileInputStream fi 	= new FileInputStream(fgmt);
  			BufferedReader 	br 	= new BufferedReader(new InputStreamReader(fi));
  			String sl = new String();
  			while ((sl=br.readLine())!=null){
  				if(sl.split("\t").length>=3){ // id name genelist
  					String sp[]=sl.split("\t");
  					String id=sp[0];
  					String name=sp[1];
  					String glist=sl.replace(id+"\t"+name+"\t", "");
  					int ngene=glist.split("\t").length;
  					if(ngene>=genesetfilter.get(0) && ngene<=genesetfilter.get(1)) {
  						hgmt_names.put(id,name);
  	  					hgmt_genelist.put(id,glist);
  	  					count++;
  					}
  				}else{
  					System.out.println("gmt is not a valid GMT object");
  					System.out.println("No pathways in gmt made the geneset.filter");
  					JOptionPane.showMessageDialog(null,"No pathways in gmt made the geneset.filter","Warning",JOptionPane.ERROR_MESSAGE);
  					break;
  				}
  			}
  			br.close();
  			fi.close();
		}catch (Exception e) {
			System.out.println("Error"+e);
		}
	}
	private void writeGMT(String fout, Hashtable<String, String> hgmt_names,Hashtable<String, String> hgmt_genelist, ArrayList<String>sortedEntry){
		io.check_valid_fout(fout);
		try {
			for (int i = 0; i < sortedEntry.size(); i++) {
				String name=hgmt_names.get(sortedEntry.get(i));
				String glist=hgmt_genelist.get(sortedEntry.get(i));
				io.writeFile(fout,sortedEntry.get(i)+"\t"+name+"\t"+glist+"\r\n");
			}	
		} catch (Exception e) {
			// TODO: handle exception
		}
	}
	private ArrayList<String> makeBackground(String fgmt, ArrayList<String>background) {
		background = new ArrayList<String>();
		Hashtable<String, String> hgmt_background= new Hashtable<String, String>();
		try {
            FileInputStream fi 	= new FileInputStream(fgmt);
  			BufferedReader 	br 	= new BufferedReader(new InputStreamReader(fi));
  			String sl = new String();
  			while ((sl=br.readLine())!=null){
  				String sp [] = sl.split("\t");
  				for (int i = 2; i < sp.length; i++) {
  					hgmt_background.put(sp[i], "");
				}
  			}
  			br.close();
  			fi.close();
		}catch (Exception e) {
			System.out.println("Error"+e);
		}
	    List<String> bkg = new ArrayList<String>(hgmt_background.keySet());
	    Collections.sort(bkg);
	    for (String bkglist:bkg) {
	    	background.add(bkglist);
	    }
	    return background;
	 }
	private Boolean check_genesetfilter(String genesetLower, String genesetUpper){
		boolean check_genefilter=false;
		if(!isNumeric(genesetLower) || !isNumeric(genesetUpper)){
			System.out.println("geneset.filter must be numeric");
		}else{
			if((!isInteger(genesetLower) || !isInteger(genesetUpper)) && (genesetLower.startsWith("-") || genesetUpper.startsWith("-"))){
				System.out.println("geneset.filter must be integer and positive");
			}else{
				if(!isInteger(genesetLower) || !isInteger(genesetUpper)){
					System.out.println("geneset.filter must be integer");
				}else if(genesetLower.startsWith("-") || genesetUpper.startsWith("-")){
					System.out.println("geneset.filter must be positive");	
				}else{
					check_genefilter=true;	
				}
			}
		}
		return check_genefilter;
	}
	private Boolean check_contribution(String fscore){
		boolean checkcontribution = false;
		int firstln=0;
		try {
            FileInputStream fi 	= new FileInputStream(fscore);
  			BufferedReader 	br 	= new BufferedReader(new InputStreamReader(fi));  			
  			String sl = new String();
  			while ((sl=br.readLine())!=null){
  				int ncol=sl.split("\t").length;
  				if(firstln>0){
  					if(ncol>2){
  						checkcontribution=true;
  					}else{
  						System.out.println("scores contains only one column. Column contributions will not be calculated");
  						break;
  					}
  				}else{
  					firstln=1;
  				}
  			}
  			br.close();
  			fi.close();
		}catch (Exception e) {
			System.out.println("Error"+e);
		}
		return checkcontribution;
	} 
	private String check_cytoscapefilenames(String fcyto, String fscore){
		String newfcyto=fcyto;
		if(fcyto!=null){
			if(check_contribution(fscore)&& fcyto.split(",").length!=3){
				System.out.println("Must supply 3 file names to cytoscape.filenames");
			}
		}
		if(!check_contribution(fscore)){
			if(fcyto.split(",").length<2){
				System.out.println("Must supply 2 file names to cytoscape.filenames");
			}
			if(fcyto.split(",").length==3){
				System.out.println("Column contributions will not be evaluated\r\nso the contribution matrix is not being written.\r\ncytoscape.filenames[2] will be ignored");
				String sp[]=fcyto.split(",");
				newfcyto=sp[0]+","+sp[2];
			}
		}
		return newfcyto;
	}
	private Hashtable<String, String> filter_gmt(String genesetLower, String genesetUpper, Hashtable<String, String> hgmt_genelist,Hashtable<String, String> hgmtfiltered_genelist){
		if(genesetLower!=null && genesetUpper!=null){
			int lower=Integer.parseInt(genesetLower);
			int upper=Integer.parseInt(genesetUpper);
			Enumeration e = hgmt_genelist.keys();
			while (e.hasMoreElements()) {
				Object id = e.nextElement();
				String glist=hgmt_genelist.get(id.toString());
				long gsize= glist.split("\t").length;
				if(gsize>=lower && gsize<=upper){
					hgmtfiltered_genelist.put(id.toString(),glist);
				}				
			}
			if(hgmtfiltered_genelist.size()==0){
				System.out.println("No pathways in gmt made the geneset.filter");
			}
			if(hgmtfiltered_genelist.size()<hgmt_genelist.size()){
				System.out.println((hgmt_genelist.size()-hgmtfiltered_genelist.size())+" terms were removed from gmt because they did not make the geneset.filter");
			}
		}
		return hgmtfiltered_genelist;
	}

	private double[][] remove_gene_notinbackground(RealMatrix datamatrix,ArrayList<String>genelist,ArrayList<String>background, ArrayList<String>genelistUpdate){
		ArrayList<Integer> validInd = new ArrayList<Integer>(); 
		Collections.sort(background);
		for (int i = 0; i < genelist.size(); i++) {
			if(Collections.binarySearch(background,genelist.get(i))>=0){
				validInd.add(i);
				genelistUpdate.add(genelist.get(i));
			}
		}
		double newdata[][]= new double[validInd.size()][datamatrix.getColumnDimension()];
		for (int i = 0; i < newdata.length; i++) {
			newdata[i]=datamatrix.getRow(validInd.get(i));
		}
		return newdata;
	}
	
	private Hashtable<String, Double> merge_pval(double[][] matrix,ArrayList<String>genelist,String method, boolean checkVector){
		int firstln=0;
		double newscores [] = new double[1];
		double pbrown[] = new double[matrix.length];
		Hashtable<String,Double> mergescores= new Hashtable<String, Double>();
		if(method.toLowerCase().equals("fisher")){method="sumlog";}  //sumlog = combine p-values by the sum of logs (Fisher's)method
		if(checkVector){
			if(method.toLowerCase().equals("brown")){
				JOptionPane.showMessageDialog(null,"Brown's method cannot be used with a single list of p-values","Warning",JOptionPane.ERROR_MESSAGE);
			}
		}else{
			//# scores is a matrix
			if(method.toLowerCase().equals("brown")){
				double[][] tmatrix=stat.transpose(matrix);
				RealMatrix covmatrix = calculateCovariances(tmatrix);
				RealMatrix datamatrix = MatrixUtils.createRealMatrix(matrix);
				for (int i = 0; i < datamatrix.getRowDimension(); i++) {
					double eachrow[] = datamatrix.getRow(i);
					pbrown[i]=brownsMethod(eachrow, covmatrix);
					mergescores.put(genelist.get(i),pbrown[i]);
				}
			}
		}
		return mergescores;
	}
	public Boolean isInteger(String str){
		if (str == null) {
	        return false;
	    }
	    int length = str.length();
	    if (length == 0) {
	        return false;
	    }
	    int i = 0;
	    if (str.charAt(0) == '-') {
	        if (length == 1) {
	            return false;
	        }
	        i = 1;
	    }
	    for (; i < length; i++) {
	        char c = str.charAt(i);
	        if (c < '0' || c > '9') {
	            return false;
	        }
	    }
	    return true;
	}
	public Boolean isNumeric(String str){
		if (str == null) {
	        return false;
	    }
	    int length = str.length();
	    if (length == 0) {
	        return false;
	    }
	    int i = 0;
	    if (str.charAt(0) == '-') {
	        if (length == 1) {
	            return false;
	        }
	        i = 1;
	    }
	    if(!str.matches("^[0-9eE_.-]+$")){
	        return false;
	    }
	    
		return true;
	}
	public Integer countLine(String fin){
		int firstln=0;
		int ln=0;
		try {
            FileInputStream fi 	= new FileInputStream(fin);
  			BufferedReader 	br 	= new BufferedReader(new InputStreamReader(fi));  			
  			String sl = new String();
  			while ((sl=br.readLine())!=null){
  				if(firstln==0){
  					ln++;
  				}else{
  					firstln=1;
  				}
  			}
  			br.close();
  			fi.close();
		}catch (Exception e) {
			System.out.println("Error"+e);
		}
		return ln;
	}
	public Boolean isVector(String fin){
		//count line ==1
		boolean checkvector=false;
		int lines = 0;
		try {
			BufferedReader reader = new BufferedReader(new FileReader("fin"));
			while (reader.readLine() != null){
				lines++;
				if(lines>2){
					checkvector=true; 
					break;
				}
			}
			reader.close();	
		} catch (Exception e) {
			// TODO: handle exception
		}
		return checkvector;
	}
	public Double [] format_pval01(String fin){
		// Some metap functions don't like p-values that are 0 or 1 so make them (0, 1) to avoid errors
		int firstln=0;
		Hashtable<String,ArrayList<Double>> hnewscores=new Hashtable<String, ArrayList<Double>>();
		ArrayList<Double> newscores = new ArrayList<Double>();
		try {
            FileInputStream fi 	= new FileInputStream(fin);
  			BufferedReader 	br 	= new BufferedReader(new InputStreamReader(fi));  			
  			String sl = new String();
  			while ((sl=br.readLine())!=null){
  				if(firstln>0){
  					ArrayList<Double> pvallist= new ArrayList<Double>();
  					String sp [] = sl.split("\t");
  					for (int i =1; i < sp.length; i++) {
  						try {
  							double newpval=0.0;
  							double pval=Double.parseDouble(sp[i]);
  							if(pval==0){
  								newpval=1e-16;
  							}else if(pval==1){
  								newpval=1-(1e-16);
  							}else{
  								newpval=pval;
  							}
  							pvallist.add(newpval);
						} catch (Exception e) {
							// TODO: handle exception
						}
					}
  					hnewscores.put(sp[0], pvallist);
  				}else{
  					firstln=1;
  				}
  			}
  			br.close();
  			fi.close();
		}catch (Exception e) {
			System.out.println("Error"+e);
		}
		Enumeration e = hnewscores.keys();
		while (e.hasMoreElements()) {
			Object gene = e.nextElement();
			ArrayList<Double> values=hnewscores.get(gene.toString());
			for (int i = 0; i < values.size(); i++) {
				newscores.add(values.get(i));	
			}
		}
		
	    Double scores_array[] = new Double[newscores.size()];
	    scores_array = newscores.toArray(scores_array);
		
	    ArrayList<String> al = new ArrayList<String>();
	    String ar[] = new String[al.size()];
	    ar = al.toArray(ar);
	    
	    
		return scores_array;
	}
	public RealMatrix calculateCovariances(double tmatrix[][]){
		double tdatamatrix[][]=new double[tmatrix[0].length][tmatrix.length];
		for (int i = 0; i < tmatrix.length; i++) {
			double tdata[]=transformData(tmatrix[i]);
			for (int j = 0; j < tdata.length; j++) {
				tdatamatrix[j][i]=tdata[j];
			}
		}
		RealMatrix matrix = MatrixUtils.createRealMatrix(tdatamatrix);
		RealMatrix covmatrix = new Covariance(matrix).getCovarianceMatrix();
		return covmatrix;
	}
	public double[] transformData(double[] dat){
		/* If all values in dat are the same (equal to y), return dat. 
		 * The covariance matrix will be the zero matrix, 
		 * and brown's method gives the p-value as y
		 * Otherwise (dat - dvm) / dvsd is NaN and ecdf throws and error
		 */
		Hashtable<Double,Double> cdfmap = new Hashtable<Double, Double>();
		double transformeddata[] = new double[dat.length];
		double s[] = new double[dat.length];
		double dvm;
		double dvsd;
		double max=stat.getMax(dat);
		double min=stat.getMin(dat);
		
		if((max-min)<1.5E-8){
			transformeddata=dat;
			System.out.println("near equality");
		}else{
			System.out.println("need transformation");
			dvm=stat.getMean(dat);
			dvsd=stat.getSDpopulation(dat);
			for (int i = 0; i < dat.length; i++) {
				s[i]=(dat[i]-dvm)/dvsd;
			}
			
			cdfmap=stat.EmpiricalDistribution(s);
			for (int i = 0; i < s.length; i++) {
				if(cdfmap.containsKey(s[i])){
					double newvalue=-2*Math.log(cdfmap.get(s[i]));
					if(newvalue>1E-7){
						transformeddata[i]=newvalue;
					}else{
						transformeddata[i]=0;
					}
				}
			}
		}
		return transformeddata;
	}
	
	public double brownsMethod(double[] each_rowscores, RealMatrix covmatrix){
		double logpval[]=new double[each_rowscores.length];
		double covsum= 0.0;
		double var = 0.0;
		double sf = 0.0;
		double df =0.0;
		double x=0.0;
		double q=0.0;
		double cdf = 0.0;
		int N = covmatrix.getColumnDimension();
		int expected = 2*N;
		
		if(each_rowscores.length==0 && !covmatrix.isSquare()){
			System.out.println("Either data.matrix or cov.matrix must be supplied");
		}
		if(each_rowscores.length>0 && covmatrix.isSquare()){
			System.out.println("Both data.matrix and cov.matrix were supplied. Ignoring data.matrix");
		}
		
		double entry[]=stat.getlowertri(covmatrix, false);
		covsum = 2*stat.getSum(entry);
		var=(4*N)+covsum;
		sf=var/(2*expected);
		df=(2*Math.pow(expected,2))/var;
		if(df>(2*N)){
			df=2*N;
			sf=1;
		}
		for (int i = 0; i < each_rowscores.length; i++) {
			try {
				logpval[i]=-Math.log(each_rowscores[i]);	
			} catch (Exception e) {
				System.out.println("Error (pval): "+e);
			}
		}
		x=2*stat.getSum(logpval);
		q=x/sf;
		cdf = Gamma.regularizedGammaQ((df/2),(q/2));
		return cdf;
	}
	
	public Hashtable<String,String> enrichmentAnalysis(ArrayList<String>genelist,String fgmt,String fout, Hashtable<String, String> hgmt_names, Hashtable<String, String> hgmt_genelist, ArrayList<String> background){
		System.out.println("Running Enrichment Analysis ...");
		io.check_valid_fout(fout);
		ArrayList<String> annotations = new ArrayList<String>();
		ArrayList<Double> orderpval = new ArrayList<Double>();
		Hashtable<String, String> results = new Hashtable<String, String>();
		Hashtable<Integer,String> initialBkg = new Hashtable<Integer, String>();
		
		int bkgsize_begin = background.size();
		for (int i = 0; i < background.size(); i++) {
			initialBkg.put(i, background.get(i));
		}		
		
		Collections.sort(background);
		Enumeration e = hgmt_names.keys();
		while (e.hasMoreElements()) {
			Object term= e.nextElement();
			String GOID = term.toString();
			annotations = new ArrayList<>(Arrays.asList(hgmt_genelist.get(GOID).split("\t")));
			Collections.sort(annotations);
			if(background.size()!=bkgsize_begin){
				background=new ArrayList<String>(initialBkg.values());
				Collections.sort(background);
			}
			orderpval=orderedHypergeometric(genelist,background,annotations);
			double pval = orderpval.get(0);
			double index = orderpval.get(1);
			int ind = (int)index;
			StringBuffer overlap=new StringBuffer();
			for (int i = 0; i < ind+1; i++){
				if(Collections.binarySearch(annotations,genelist.get(i))>=0){
					if(i==0){
						overlap.append(genelist.get(0));
					}else{
						overlap.append(","+genelist.get(i));	
					}	
				}
			}
			
			results.put(GOID,hgmt_names.get(GOID)+"\t"+Double.toString(pval)+"\t"+annotations.size()+"\t"+annotations);
		}
		return results;
	}
	
	public double hypergeometric(RealMatrix counts){
		/* x, q = vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls.
		 * m = the number of white balls in the urn.
		 * n = the number of black balls in the urn.
		 * k = the number of balls drawn from the urn.
		 * phyper(x-1, m, n, k, lower.tail=FALSE)
		 */
		int m=(int) (counts.getEntry(0,0)+counts.getEntry(1,0));
		int n=(int) (counts.getEntry(0,1)+counts.getEntry(1,1));
		int k=(int) (counts.getEntry(0,0)+counts.getEntry(0,1));
		int x=(int) (counts.getEntry(0,0));
		int populationSize=m+n;
		int numberOfSuccesses=m;
		int sampleSize=k;
		HypergeometricDistribution hyper=new HypergeometricDistribution(populationSize,numberOfSuccesses,sampleSize);
		double pval=hyper.upperCumulativeProbability(x);
		return pval;
	 }
 
	public ArrayList<Double> orderedHypergeometric(ArrayList<String> genelist,ArrayList<String> background,ArrayList<String> annotations){
		 ArrayList<Integer> whichin = new ArrayList<Integer>();
		 ArrayList<String> cl = new ArrayList<String>();
		 ArrayList<String> annotvalid_cl = new ArrayList<String>();
		 ArrayList<Double> scorelist = new ArrayList<Double>();
		 ArrayList<Double> results = new ArrayList<Double>();
		 RealMatrix counts = null;
		 double pval=1.0;
		 int ind = -1;
		 double scores=-1.0;
		 
		 Collections.sort(background);
		 Collections.sort(annotations);
		 for (int i = 0; i < genelist.size(); i++) {
			 int index=Collections.binarySearch(annotations,genelist.get(i));
			 if(index>=0){
				 whichin.add(i);
			 }
		}
		
		if(whichin.size()==0){
			pval=1.0; ind=0;	
		}else{
			Collections.sort(whichin);
			int gl_length=whichin.get(0)+1;
			try {
				for (int i = 0; i < gl_length; i++) {
					String gl=genelist.get(i);
					int bkgvalidIND = Collections.binarySearch(background,gl);
					if(bkgvalidIND>=0){
						background.remove(bkgvalidIND);	
					}
				}
				cl =background;
				annotvalid_cl=new ArrayList<String>();
				for (int j = 0; j < cl.size(); j++) {
					int annotvalidIND = Collections.binarySearch(annotations,cl.get(j));
					if(annotvalidIND>=0){
						annotvalid_cl.add(annotations.get(annotvalidIND));
					}
				}
				int genelist0=gl_length-1;
				int complement1=annotvalid_cl.size();
				int complement0=cl.size()-complement1;
				counts = new Array2DRowRealMatrix(new double[][] {{1,complement1},{genelist0,complement0}},false);
				scores = hypergeometric(counts);
				scorelist.add(scores);
			} catch (Exception e) {
			}
			if(whichin.size()==1){
				pval=scores; ind=whichin.get(0);
				scorelist.add(scores);
			}else{
				/*
				# Update counts and recalculate score for the rest of the indeces in which.in
			    # The genes in genelist[which.in[i]:which.in[i-1]] are added to the genes
			    # being tested and removed from the complement. Of these, 1 will always be
			    # in annotations and the rest will not. Therefore we can just modify the
			    # contingency table rather than recounting which genes are in annotations
				*/
				for (int i = 1; i < whichin.size(); i++){
					int diff =whichin.get(i)-whichin.get(i-1);
					counts = new Array2DRowRealMatrix(new double[][] {{i+1,(counts.getEntry(0, 1)-1)},{(counts.getEntry(1, 0)+diff-1),(counts.getEntry(1, 1)-diff+1)}},false);
					scores=hypergeometric(counts);
					scorelist.add(scores);
				}
				double minscore=scorelist.get(0);
				for (int i =1; i < scorelist.size(); i++) {
					minscore=Math.min(minscore, scorelist.get(i));
				}
				
				int index=0;
				for (int i = 0; i < scorelist.size(); i++) {
					if(scorelist.get(i)==minscore){
						index=i;
					}
				}
				ind=whichin.get(index);
				pval=minscore;
			}
		}
		//# Return the lowest p-value and the associated index
		results.add(pval);
		results.add((double)ind);
		return results;
	 }
	
	public ArrayList<String> getGenelist(String fin, Integer startrow){
			int ln =0;
			 ArrayList<String> genelist = new ArrayList<String>();
			try {
	            FileInputStream fi 	= new FileInputStream(fin);
	  			BufferedReader 	br 	= new BufferedReader(new InputStreamReader(fi));  			
	  			String sl = new String();
	  			while ((sl=br.readLine())!=null){
	  				if(ln>=startrow){
	  					genelist.add(sl.split("\t")[0].toUpperCase());
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
			
			return genelist;
	}
	
	public static ArrayList<Map.Entry<String, Double>>  sortHashMapValue(Hashtable<String,Double> mergedscores,final boolean decreasing){
		       ArrayList<Map.Entry<String, Double>> l = new ArrayList(mergedscores.entrySet());
		       Collections.sort(l, new Comparator<Map.Entry<String, Double>>(){
	    	       public int compare(Map.Entry<String, Double> o1, Map.Entry<String, Double> o2) {
	    	    	   if(decreasing){
	    	    		   return o2.getValue().compareTo(o1.getValue()); //max2min  
	    	    	   }else{
	    	    		   return o1.getValue().compareTo(o2.getValue()); //min2max
	    	    	   }
	    	       }});
		       return l;
		}
	
	public Hashtable<String, String> columnSignificance(String fscores, String fgmt, RealMatrix datamatrix, Hashtable<String, String> hgmt_names,Hashtable<String, String> hgmt_genelist, ArrayList<String> background, double cutoff,double significant, String correctionmethod, double[] pvals){
			background = makeBackground(fgmt, background);
			ArrayList<String> genelist = getGenelist(fscores,1);
			ArrayList<String> colnames=getColnames(fscores);
			ArrayList<String> GOIDlist = new ArrayList<String>();
			Hashtable<String,Hashtable<Integer, String>>pval_eachcol= new Hashtable<String, Hashtable<Integer,String>>();
			Hashtable<String, String> res = new Hashtable<String, String>();
			Hashtable<Integer, String> newres = new Hashtable<Integer, String>();
			Hashtable<String, String> results = new Hashtable<String, String>();
			
			int[] nrow = new int[(datamatrix.getRowDimension())]; 
			int[] ncol = new int[(datamatrix.getColumnDimension())];
			for (int i = 0; i < nrow.length; i++) {nrow[i]=i;}
			for (int i = 0; i < ncol.length; i++) {ncol[i]=i;}

			for (int i = 0; i < datamatrix.getColumnDimension(); i++) {
				res = new Hashtable<String, String>();
				newres = new Hashtable<Integer, String>();
				Hashtable<String, Double> hcolscores= new Hashtable<String, Double>();
				ArrayList<Map.Entry<String, Double>> sorted_colscores= new ArrayList<Map.Entry<String,Double>>();
				RealMatrix colscores = datamatrix.getSubMatrix(nrow,new int[]{i});
				for (int j = 0; j < colscores.getRowDimension(); j++) {
					if(colscores.getEntry(j,0)<=cutoff){
						hcolscores.put(genelist.get(j), colscores.getEntry(j,0));
					}
				}
				ArrayList<String> glist = new ArrayList<String>();
				sorted_colscores=sortHashMapValue(hcolscores,false);
				for (int j = 0; j < sorted_colscores.size(); j++) {
					glist.add(sorted_colscores.get(j).getKey());
				}
				res=enrichmentAnalysis(glist, fgmt,"",hgmt_names,hgmt_genelist,background);
				ArrayList<Double> enrpval = new ArrayList<Double>();
				Enumeration e = res.keys();
				while (e.hasMoreElements()) {
					Object GOID = e.nextElement();
					String info[]=res.get(GOID.toString()).split("\t");
					double pval=Double.parseDouble(info[1]);
					enrpval.add(pval);
				}
				int n=0;
				double adjpval[] = pcorrect.pAdjust(enrpval, correctionmethod);
				Enumeration e1 = res.keys();
				while (e1.hasMoreElements()) {
					Object GOID = e1.nextElement();
					String info[]=res.get(GOID.toString()).split("\t");
					double pval = adjpval[n];
					info[1]=Double.toString(pval);

					if(pval>cutoff){
						info[info.length-1]="NA"; //overlap
					}
					StringBuffer newinfo= new StringBuffer();
					newinfo.append(info[0]);
					for (int j = 1; j < info.length; j++) {
						newinfo.append("\t"+info[j]);
					}
					res.put(GOID.toString(),newinfo.toString());
					newres.put(n, GOID.toString()+"\t"+newinfo.toString());
					if(i==0){
						GOIDlist.add(GOID.toString());	
					}
					n++;
				}
				pval_eachcol.put(colnames.get(i), newres);
			}
			/*******************************/
			//    Get Evidence
			/*******************************/
			for (int line = 0; line < res.size(); line++) {
				String eachrow[] = new String[colnames.size()];
				String GOID =GOIDlist.get(line); 
				int count_overlap=0;
				for (int j = 0; j < pval_eachcol.size(); j++) {
					Hashtable<Integer, String> eachcolres = pval_eachcol.get(colnames.get(j));
					String info[]=eachcolres.get(line).split("\t");
					String value = info[info.length-1];
					if(value.toUpperCase().equals("NA")){
						eachrow[j]="NA";
					}else{
						eachrow[j]=colnames.get(j);
						count_overlap++;
					}
				}	
				String ev="";
				StringBuffer evbuf =new StringBuffer();
				if(count_overlap==0){
					if(pvals[line]<=significant){
						ev="combined";
					}else{
						ev="none";
					}
				}else{
					for (int l = 0; l <colnames.size(); l++) {
						if(!eachrow[l].equals("NA") && !eachrow[l].isEmpty()){
							if(l==0){
								evbuf.append(eachrow[l]);	
							}else{
								evbuf.append(","+eachrow[l]);
							}
						}
					}
					ev=evbuf.toString();
				}
				results.put(GOID,res.get(GOID)+"\t"+ev);
			}
			return results;
		}
	
	private ArrayList<String> getColnames(String fin){
		ArrayList<String> colnames= new ArrayList<String>();
		try {
            FileInputStream fi 	= new FileInputStream(fin);
  			BufferedReader 	br 	= new BufferedReader(new InputStreamReader(fi));  			
  			String sl[] = br.readLine().split("\t");
  			for (int i = 1; i < sl.length; i++) {
				colnames.add(sl[i]);
			}
  			br.close();
  			fi.close();
		}catch (Exception e) {
			System.out.println("Error"+e);
		}
		return colnames;
	}
	
	private void prepareCytoscape(Hashtable<String, String> enrichresults,Hashtable<String, String> hgmt_names_sig,Hashtable<String, String> hgmt_genelist_sig, Hashtable<String, String> sigcols, String fcyto_list){
		String fcyto[] = fcyto_list.split(",");
		String fcyto1=fcyto[0]; io.check_valid_fout(fcyto1);
		String fcyto2=fcyto[1]; io.check_valid_fout(fcyto2);
		String fcyto3=fcyto[2]; io.check_valid_fout(fcyto3);
		
		Hashtable<String, String> colsig =new Hashtable<String, String>();
		Hashtable<String, Double> printout_cyto1 =new Hashtable<String, Double>();
		Hashtable<String, Double> printout_cyto3 =new Hashtable<String, Double>();
		ArrayList<Map.Entry<String, Double>> sorted_cyto1 = new ArrayList<Map.Entry<String,Double>>();
		ArrayList<Map.Entry<String, Double>> sorted_cyto3 = new ArrayList<Map.Entry<String,Double>>();
		Hashtable<String, Integer> valid_ev = new Hashtable<String, Integer>();
		StringBuffer ev_all= new StringBuffer();
		
		if(sigcols.size()>0){
	    	Enumeration e = sigcols.keys();
	    	while (e.hasMoreElements()) {
				Object GOID = e.nextElement();
				String info[]=sigcols.get(GOID.toString()).split("\t");
				String ev[] = info[info.length-1].split(",");
				for (int i = 0; i < ev.length; i++) {
					if(!ev[i].equals("none")&&!ev[i].isEmpty()){
						valid_ev.put(ev[i],0);	
					}	
				}
			}
	    	
	    	int ind=0; 
	    	List<String> ev = new ArrayList<String>(valid_ev.keySet());
		    Collections.sort(ev);
		    for (String each: ev) {
		    	valid_ev.put(each,ind);
		    	if(each.contains("combined")) {
		    		if(ind==0){
			    		ev_all.append(each);
			    	}else{
			    		ev_all.append(","+each);
			    	}	
		    	}else {
		    		if(ind==0){
			    		ev_all.append("Genes_"+each);
			    	}else{
			    		ev_all.append(",Genes_"+each);
			    	}	
		    	}
		    	ind++;
		    }
		        		
		    // Use pichart
	    	String colcode=genRGBColorCode(valid_ev.size());
	    	String colnames=ev_all.toString().replace("Genes_", "");
	    	String instruct = "piechart: attributelist=\""+colnames+"\" colorlist=\""+colcode+"\" showlabels=FALSE";
	    	int line=0;
	    	try {
	    		io.writeFile(fcyto1, "term.id\tterm.name\tadjusted.p.val\r\n");
	    		io.writeFile(fcyto3, "term.id\t"+colnames.replace(",", "\t")+"\tinstruct\r\n");
	    		
	    		String colindx []=colnames.split(",");
	    		Enumeration e1 = sigcols.keys();
		    	while (e1.hasMoreElements()) {
		    		Object GOID = e1.nextElement();
		    		if(enrichresults.containsKey(GOID.toString())){
		    			String info[]=sigcols.get(GOID.toString()).split("\t");
			    		String ovlp=info[3].replace("[","").replace("]", "").replace(" ","");
			    		String evch[]=info[info.length-1].split(",");
			    		StringBuffer value = new StringBuffer();
			    		Arrays.sort(evch); 
						for (int i = 0; i < colindx.length; i++) {
							if(Arrays.binarySearch(evch,colindx[i])>=0){
								if(i==0) {
									value.append("1");
								}else{
									value.append("\t1");
								}
							}else {
								if(i==0) {
									value.append("0");
								}else{
									value.append("\t0");
								}
							}
						}
						colsig.put(GOID.toString(),value.toString()+"\t"+instruct);
						String sepInfo[]=enrichresults.get(GOID.toString()).split("\t");
						double adjpval=Double.parseDouble(sepInfo[1]);
						printout_cyto1.put(GOID.toString()+"\t"+sepInfo[0],adjpval);
						printout_cyto3.put(GOID.toString()+"\t"+value.toString()+"\t"+instruct,adjpval);
						line++;	
		    		}else{
		    			//System.out.println("No significant GOID found from Enrichment Analysis. Cytoscape files will not be written.");
		    		}
		    	}
		    	DecimalFormat decimalformatter = new DecimalFormat("#0.0000");
		    	decimalformatter.setMaximumFractionDigits(ndigitLimit);
				sorted_cyto1 =sortHashMapValue(printout_cyto1,false); //min to max (most to least sig)
		    	sorted_cyto3 =sortHashMapValue(printout_cyto3,false); 
		    	ArrayList<String> sortedEntry =new ArrayList<String>();
		    	for (int i = 0; i < sorted_cyto1.size(); i++) {
		    		io.writeFile(fcyto1, sorted_cyto1.get(i).getKey()+"\t"+decimalformatter.format(sorted_cyto1.get(i).getValue())+"\r\n");
		    		io.writeFile(fcyto3, sorted_cyto3.get(i).getKey()+"\r\n");
		    		sortedEntry.add(sorted_cyto1.get(i).getKey().split("\t")[0]);
				}
		    	writeGMT(fcyto2, hgmt_names_sig, hgmt_genelist_sig,sortedEntry);
		    	
			} catch (Exception e2) {
				// TODO: handle exception
			}
	    }
	}
	
	private String genRGBColorCode(int ncolor){
		Hashtable<Integer,Color> RGBpalette = new Hashtable<Integer, Color>();
		Color color = new Color(0);
		StringBuffer hexlist = new StringBuffer();
		int paletteCode=0;
		if(ncolor<=10){
			paletteCode=1;
		}else if(ncolor<=20){
			paletteCode=2;
		}else if(ncolor>20){
			paletteCode=3;
		}
		 
		RGBpalette=getRGBpalette(paletteCode, ncolor);
		for (int i = 0; i < ncolor; i++) {
			color=RGBpalette.get(i);
			String hex = String.format("#%02X%02X%02X", color.getRed(),color.getGreen(),color.getBlue());  
			if(i==0) {
				hexlist.append(hex);		
			}else {
				hexlist.append(","+hex);	
			}
		}
		return hexlist.toString();
	}
	
	private Hashtable<Integer,Color> getRGBpalette(int hexSeries, int ncolor){
		//25 color codes
		Hashtable<Integer,Color> hexTable = new Hashtable<Integer,Color>();
		//Debug: 
		/*hexTable.put(0,"#01B8AA");
		hexTable.put(1,"#374649");
		hexTable.put(2,"#FD625E");
		hexTable.put(3,"#F2C80F");
		hexTable.put(4,"#5F6B6D");
		hexTable.put(5,"#8AD4EB");
		hexTable.put(6,"#FE9666");
		hexTable.put(7,"#A66999");
		hexTable.put(8,"#3599B8");
		hexTable.put(9,"#DFBFBF");
		*/

		if(hexSeries==1) {
			//Default 10 groups
			hexTable=gen_Hex10(hexTable);
		}else if(hexSeries==2 ||hexSeries==3){
			//Default 20 groups
			hexTable=gen_Hex20(hexTable);
		}else if(hexSeries==3){
			// if ncolor >20
			// randomly generate RGB codes
			hexTable=gen_Hex20(hexTable);
			hexTable=gen_Hex21up(hexTable,ncolor);
		}
		return hexTable;
	}
	
	private Hashtable<Integer,Color> gen_Hex10(Hashtable<Integer,Color> hexTable){
		hexTable.put(0,new Color(1,184,170));
		hexTable.put(1,new Color(55,70,73));
		hexTable.put(2,new Color(253,98,94));
		hexTable.put(3,new Color(242,200,15));
		hexTable.put(4,new Color(95,107,109));
		hexTable.put(5,new Color(138,212,235));
		hexTable.put(6,new Color(254,150,102));
		hexTable.put(7,new Color(166,105,153));
		hexTable.put(8,new Color(53,153,184));
		hexTable.put(9,new Color(223,191,191));
		return hexTable;
	}
	
	private Hashtable<Integer,Color> gen_Hex20(Hashtable<Integer,Color> hexTable){
		hexTable.put(0,new Color(31,119,180));
		hexTable.put(1,new Color(174,199,232));
		hexTable.put(2,new Color(255,127,14));
		hexTable.put(3,new Color(255,187,120));
		hexTable.put(4,new Color(44,160,44));
		hexTable.put(5,new Color(152,223,138));
		hexTable.put(6,new Color(214,39,40));
		hexTable.put(7,new Color(255,152,150));
		hexTable.put(8,new Color(148,103,180));
		hexTable.put(9,new Color(197,176,213));
		hexTable.put(10,new Color(140,86,75));
		hexTable.put(11,new Color(196,156,148));
		hexTable.put(12,new Color(227,119,194));
		hexTable.put(13,new Color(247,182,210));
		hexTable.put(14,new Color(127,127,127));
		hexTable.put(15,new Color(199,199,199));
		hexTable.put(16,new Color(188,189,34));
		hexTable.put(17,new Color(219,219,141));
		hexTable.put(18,new Color(23,190,207));
		hexTable.put(19,new Color(158,218,229));
		return hexTable;
	}
	
	private Hashtable<Integer,Color> gen_Hex21up(Hashtable<Integer,Color> hexTable,int ncolor){
		int nround=0;
		while (nround<ncolor) {
			int R=(int)(Math.random()*256+0);
			int G=(int)(Math.random()*256+0);
			int B=(int)(Math.random()*256+0);
			if(!hexTable.containsValue(new Color(R,G,B))) {
				hexTable.put(nround,new Color(R,G,B));
				nround++;
			}
		}
		return hexTable;
	}
}


