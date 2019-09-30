package org.jurilab.apw.internal.demo5;

import javax.swing.*;
import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Map;
import java.util.Random;
import java.util.Scanner;
import java.util.regex.Pattern;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.text.DefaultEditorKit;
import javax.swing.text.JTextComponent;
import javax.swing.text.TextAction;

import org.cytoscape.application.CyApplicationManager;
import org.cytoscape.command.CommandExecutorTaskFactory;
import org.cytoscape.event.CyEventHelper;
import org.cytoscape.service.util.CyServiceRegistrar;
import org.cytoscape.view.model.CyNetworkView;
import org.cytoscape.view.model.CyNetworkViewManager;
import org.cytoscape.view.vizmap.VisualMappingManager;
import org.cytoscape.view.vizmap.VisualStyle;
import org.cytoscape.work.SynchronousTaskManager;
import org.cytoscape.work.FinishStatus;
import org.cytoscape.work.ObservableTask;
import org.cytoscape.work.TaskIterator;
import org.cytoscape.work.TaskObserver;
import org.jdesktop.swingx.JXFrame;
import org.jdesktop.swingx.JXTitledPanel;


public class APWManager implements TaskObserver{
	final CyServiceRegistrar serviceRegistrar;
	final CyEventHelper eventHelper;
	public CommandExecutorTaskFactory commandTaskFactory = null;
	public SynchronousTaskManager<?> taskManager = null;
	public int TIME_VISIBLE = 2000;  
	private String modelName = null;
	private int modelNumber = -1;
	private int threadState = -1;
	private IOhandling io = new IOhandling();
	private DecimalFormat decimalformatter = new DecimalFormat("#0.0000");
	private DecimalFormat integerformatter=new DecimalFormat("###");
	private RunActivePathways APW = new RunActivePathways();
	private Hashtable<String, String> hgmt_names= new Hashtable<String, String>();
	private Hashtable<String, String> hgmt_genelist= new Hashtable<String, String>();
	private String currentPath = Paths.get("").toAbsolutePath().toString();
	private String mergelist [] = new String [] {"Brown", "Fisher", "logitp", "meanp","sump", "sumz", "sumlog"};
	private String correctlist [] = new String [] {"fdr", "holm", "hochberg", "hommel","bonferroni", "BH", "BY", "none"};
	private String yesno [] = new String[]{"no","yes"};
	private String datasetedgelist[] = new String []{"Automatic", "Separate edge for each data set (denser)", "Combine edges across data sets (sparser)"};
	private String metriclist[] = new String []{"Jaccard + Overlap Combined","Jaccard","Overlap"};
	private String gsetmin = new String("5");
	private String gsetmax = new String("1000");
	private String cutoffapw= new String("0.1");
	private String sigapw= new String("0.05");
	private String nodeFDR= new String("0.05");
	private String nodePval= new String("0.05");
	private String edgeCutoff= new String("0.66");
	private JXFrame mainFrame = new JXFrame();
	private JFrame progressFrame = new JFrame();
	private JPanel progressPane = new JPanel();
	private JProgressBar progressBar = new JProgressBar(0,100);
	private JTextField fin;
	private JTextField gmt;
	private JTextField fcytoedit1;
	private JTextField fcytoedit2;
	private JTextField fcytoedit3;
	private JTextField networkname;
	private JFormattedTextField gset_min;
	private JFormattedTextField gset_max;
	private JFormattedTextField cutoff_apw;
	private JFormattedTextField sig_apw;
	private JFormattedTextField fdr_enr_node;
	private JFormattedTextField pval_enr_node;
	private JFormattedTextField cutoff_edge;
	private JComboBox<String> mergemethod;
	private JComboBox<String> correctmethod;
	private JComboBox<String> returnall;
	private JComboBox<String> datasetedge;
	private JComboBox<String> metric;
	private JCheckBox fcytotik1;
	private JCheckBox fcytotik2;
	private JCheckBox fcytotik3;
	private JButton jbfin;
	private JButton jbgmt;
	private JButton fcytosave1;
	private JButton fcytosave2;
	private JButton fcytosave3;
	private JButton advanceOptionButton;
	private JButton simpleOptionButton;
	private JButton buildButton;
	private JButton resetButton;
	private JButton cancelButton;
	private JSlider metricadjust;
	private JLabel adjustState = new JLabel("Jaccard (50%) + Overlap (50%)");
	private StringBuffer summary;
	private ArrayList<Integer> genesetfilter;
	private String fin_tip="A tab delimited file format that consists of a numerical matrix of p-values where each row is a gene and each column is a test. Rownames should be the genes and colnames the names of the tests. All values must be in [0,1] with missing values removed or converted to 1.";
	private String gmt_tip="A tab delimited file format that describes pathways and corresponding gene sets where rows represent pathway IDs and columns correspond to gene names of individual pathways.";
	private String gset_tip ="A numeric vector of length two giving the lower and upper limits for the size of the annotated geneset to pathways in gmt."; 
	private String gsetmin_tip="Minimum size of gene set found in pathways.";
	private String gsetmax_tip="Maximum size of gene set found in pathways.";
	private String gsetcutoff_tip="A maximum p-value for a gene to be used for enrichment analysis. Any genes with q.val > significant will be discarded before testing.";
	private String gsetsig_tip="A number in [0,1] denoting the maximum p-value for a pathway to be considered significantly enriched."; 
	private String mergemethod_tip="Method to merge p-values.";
	private String correctmethod_tip="Method to correct p-values.";
	private String runAll_tip="Whether to return results for all terms or only significant terms.";
	private String cytofile1_tip="A list of significant terms and the associated p-value. Only terms with q.val less than or equal to significant are written to this file.";
	private String cytofile2_tip="A shortened version of the supplied gmt file, containing only the terms in terms.txt";
	private String cytofile3_tip="A matrix indicating whether the significant pathways are found to be significant when considering only one column from input matrix. One indicates that that term is significant using only that column to test for enrichment analysis.";
	private String networkName_tip="The name of the EnrichmentMap network. If not provided then EnrichmentMap will automatically generate a name for the network based on the name of the first data set."; 
	private String edgestrategy_tip = "AUTOMATIC: EnrichmentMap decides which of the previous options to use.<br>"
			+ "<br>SEPARATE: Create separate edges for each data set when appropriate. A separate similarity score will be calculated for each data set.<br>"
			+ "<br>COMBINE: Gene sets with the same name are combined (set union) and then the similarity score is calculated.";
	private String edgecutoff_tip ="A similarity cutoff in [0,1]. Edges with a similarity score lower than the one entered will not be included in the network.";
	private String edgemetriclabel_tip = "A formula for calculating the similarity score.";
	private String edgemetricoption_tip = "To determine what percentage to use for "
			+ "JACCARD and OVERLAP when combining their value [0,1]<br>"
			+ "0 means 100% JACCARD and 0% OVERLAP.<br>"
			+ "1 means 0% JACCARD and 100% OVERLAP.";
	private String node_qval_tip="Value is between [0.1]. Gene set nodes with a q-value lower than the one entered will not be included in the network.";
	private String gif = "/images/progress.gif";
	private String visualStyle="/resources/APW_styles.xml";
	private ImageIcon progressImage;
	private JLabel progressImageLabel;
	private Thread mainThread;
	private Thread progressbarThread;
	private Task task;
	private boolean checkActiveProgressBar=true;
	public long runtime;
	public int progressCount=1;
	public long eachMin=60000;
	
	public APWManager(final CyServiceRegistrar cyRegistrar) {
		this.serviceRegistrar = cyRegistrar;
		this.eventHelper = serviceRegistrar.getService(CyEventHelper.class);
	}
	
	public void initComponents_simple_initial() {
	        /*****************************************************************/
		/*       Create GUI + Setup Default EnrichmentMap parameters
		/****************************************************************/
		mainFrame = new JXFrame("JFrame");
		mainFrame.setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
		mainFrame.setTitle("ActivePathways Demo V.0.5");
		JPanel contentPane = new JPanel(new GridBagLayout());
		JSplitPane moduleSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
		moduleSplitPane.setBorder(BorderFactory.createEmptyBorder());
		moduleSplitPane.setTopComponent(create_APW_tab_simple());
		moduleSplitPane.setDividerLocation(300);
		contentPane.add(moduleSplitPane, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(11, 11, 12, 12), 0, 0));
		advanceOptionButton = new JButton("Advance");
		buildButton = new JButton(" Build ");
		resetButton = new JButton("Reset");
		cancelButton = new JButton("Cancel");
		advanceOptionButton.setPreferredSize(new Dimension(30,30));
		buildButton.setPreferredSize(new Dimension(30,30));
		resetButton.setPreferredSize(new Dimension(30,30));
		cancelButton.setPreferredSize(new Dimension(30,30));
		JPanel buttonPanel =  new JPanel(new GridBagLayout());
		JPanel buttonPanelRight =  new JPanel(new GridBagLayout());
		JPanel buttonPanelLeft =  new JPanel(new GridBagLayout());
		buttonPanelRight.add(buildButton, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 12, 5, 5), 0, 0));
		buttonPanelRight.add(cancelButton,new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 5, 0), 0, 0));
		buttonPanelLeft.add(resetButton,new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 5, 5), 0, 0));
		buttonPanelLeft.add(advanceOptionButton,new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 5, 5), 0, 0));
		buttonPanel.add(buttonPanelLeft,new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 5, 200), 0, 0));
		buttonPanel.add(buttonPanelRight,new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0, GridBagConstraints.EAST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 5, 0), 0, 0));
		contentPane.add(buttonPanel, new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 5, 12), 0, 0));
		mainFrame.setContentPane(contentPane);
		mainFrame.setSize(400,610);
		mainFrame.setVisible(true);
		mainFrame.setLocationRelativeTo(null);
		mainFrame.pack();
        

        
		/*****************************************************************/
		/*  		Add ActionListener to all components
		/****************************************************************/
	        mainFrame.addWindowListener(new WindowAdapter(){
        	public void windowClosing(WindowEvent we){
	           	if(mainThread.isAlive()) {
	           		int option = JOptionPane.showConfirmDialog((Component) null, "Stop running APW ?","Alert", JOptionPane.YES_NO_CANCEL_OPTION);
	           		if(option==JOptionPane.YES_OPTION){
	           			mainFrame.dispose();
	           			progressFrame.dispose();
	           			mainThread.stop();
	           			progressbarThread.stop();
	           		}
	           	}else{
	           		int option = JOptionPane.showConfirmDialog((Component) null, "Exit?","Alert", JOptionPane.YES_NO_CANCEL_OPTION);
	           		if(option==JOptionPane.YES_OPTION){
	           			mainFrame.dispose();
	           		}
	           	}
        	}
        });

        JPopupMenu textEditMenu = new JPopupMenu();
        Action cut = new DefaultEditorKit.CutAction();
        cut.putValue(Action.NAME, "Cut");
        cut.putValue(Action.ACCELERATOR_KEY, KeyStroke.getKeyStroke("control X"));
        textEditMenu.add( cut );

        Action copy = new DefaultEditorKit.CopyAction();
        copy.putValue(Action.NAME, "Copy");
        copy.putValue(Action.ACCELERATOR_KEY, KeyStroke.getKeyStroke("control C"));
        textEditMenu.add( copy );

        Action paste = new DefaultEditorKit.PasteAction();
        paste.putValue(Action.NAME, "Paste");
        paste.putValue(Action.ACCELERATOR_KEY, KeyStroke.getKeyStroke("control V"));
        textEditMenu.add( paste );

        Action selectAll = new SelectAll();
        textEditMenu.add( selectAll );
        
        fin.setComponentPopupMenu(textEditMenu);
        fin.grabFocus();
        mainFrame.getRootPane().getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0), "Cancel");
        mainFrame.getRootPane().getActionMap().put("Cancel",new AbstractAction(){ 
        	public void actionPerformed(ActionEvent e){
        		mainFrame.dispose();
        	}
        });
        
        fin.setFocusTraversalKeysEnabled(false);
        fin.addKeyListener(new KeyAdapter() {
        	public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				if(keyCode==KeyEvent.VK_V){
					String filename=fin.getText();
					if(!filename.equals("")){
						boolean valid = isValidFile(filename,"",true);
						if(valid){
							Path path = Paths.get(filename).toAbsolutePath();
							String finpath = path.getParent().toString();
				        	String fcyto1=finpath+File.separator+"terms.txt";
				        	String fcyto2=finpath+File.separator+"abridged.gmt";
				        	String fcyto3=finpath+File.separator+"groups.txt";
			    	      	jbfin.grabFocus();
						}
					}
				}
			}
		
        	public void keyPressed(KeyEvent e){
        		int keyCode= e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String filename=fin.getText();
        			fin.setText(filename);
        	        	if(!filename.equals("")){
            			boolean valid = isValidFile(filename,"",true);
            			if(valid){
            				File path = new File(filename);
            	        	String finpath=path.getParent();
            	        	String fcyto1=finpath+File.separator+"terms.txt";
            	        	String fcyto2=finpath+File.separator+"abridged.gmt";
            	        	String fcyto3=finpath+File.separator+"groups.txt";
                	      	jbfin.grabFocus();	
            			}
            		}else{
            			JOptionPane.showMessageDialog(null,"Please select a Valid Input File","Warning",JOptionPane.ERROR_MESSAGE);
            		}
        		}
			}
		});
        jbfin.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
	        	JFileChooser fileopen = new JFileChooser();
	    	    FileFilter filter = new FileNameExtensionFilter("csv files", "csv");
	    	    fileopen.addChoosableFileFilter(filter);
	    	    int ret = fileopen.showDialog(null, "Import");
	    	    if (ret == JFileChooser.APPROVE_OPTION) {
	    	    	File file = fileopen.getSelectedFile();
	    	      	String filename = file.getPath();
	    	      	fin.setText(filename);
	    	    }
	        }
	    });
        
        gmt.setComponentPopupMenu(textEditMenu);
        gmt.setFocusTraversalKeysEnabled(false);
		gmt.addKeyListener(new KeyAdapter() {
        	public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				switch(keyCode){
				case KeyEvent.VK_V:
					if(!gmt.getText().equals("")){
						String filename=gmt.getText();
						isValidFile(filename,"",true);
					}
				}
			}
			
        	public void keyPressed(KeyEvent e) {
        		int keyCode= e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String filename=gmt.getText();
        			gmt.setText(filename);
            		if(!filename.equals("")){
            			boolean valid = isValidFile(filename,"",true);
            			if(valid){
            				jbgmt.grabFocus();	
            			}
            		}else{
            			JOptionPane.showMessageDialog(null,"Please select a Valid GMT File","Warning",JOptionPane.ERROR_MESSAGE);
            		}
        		}
			}
		});
	
		jbgmt.addActionListener(new ActionListener(){
	        public void actionPerformed(ActionEvent e) {
	        	JFileChooser fileopen = new JFileChooser();
	    	    FileFilter filter = new FileNameExtensionFilter("gmt files", "gmt");
	    	    fileopen.addChoosableFileFilter(filter);
	    	    int ret = fileopen.showDialog(null, "Import");
	    	    if (ret == JFileChooser.APPROVE_OPTION) {
	    	    	File file = fileopen.getSelectedFile();
	    	      	String filename = file.getPath();
	    	      	gmt.setText(filename);
	    	    }  	
	        }
	    });
		
		gset_min.setComponentPopupMenu(textEditMenu);
		gset_min.setFocusTraversalKeysEnabled(false);
		gset_min.addKeyListener(new KeyAdapter() {
			public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				if(keyCode==KeyEvent.VK_V){
					String value=gset_min.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							if(!Pattern.compile("[.E]").matcher(value.toUpperCase()).find()){
								int val=Integer.parseInt(value);
								if(val>0){
									gset_min = new JFormattedTextField(integerformatter);
									gset_min.setValue(val);
									gset_max.grabFocus();
									valid=false;
									gset_min.setFocusTraversalKeysEnabled(true);
								}else{
									JOptionPane.showMessageDialog(null,"Gene Set Number must be > 0","Warning",JOptionPane.ERROR_MESSAGE);
									gset_min.grabFocus();
								}	
							}else{
								JOptionPane.showMessageDialog(null,"Please put Integer Number","Warning",JOptionPane.ERROR_MESSAGE);
								gset_min.grabFocus();
							}
						}else{
							gset_min.setText(gsetmin);
						}
					}
				}
			}
			
        	public void keyPressed(KeyEvent e) {
        		int keyCode = e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String value=gset_min.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							if(!Pattern.compile("[.E]").matcher(value.toUpperCase()).find()){
								int val=Integer.parseInt(value);
								if(val>0){
									gset_min = new JFormattedTextField(integerformatter);
									gset_min.setValue(val);
									gset_max.grabFocus();
									gset_min.setFocusTraversalKeysEnabled(true);
									valid=false;	
								}else{
									JOptionPane.showMessageDialog(null,"Gene Set Number must be > 0","Warning",JOptionPane.ERROR_MESSAGE);
									gset_min.grabFocus();
								}	
							}else{
								JOptionPane.showMessageDialog(null,"Please put Integer Number","Warning",JOptionPane.ERROR_MESSAGE);
								gset_min.grabFocus();
							}
						}else{
							gset_min.setText(gsetmin);
						}
					}
        		}
			}
		});
		
		gset_max.setComponentPopupMenu(textEditMenu);
		gset_max.setFocusTraversalKeysEnabled(false);
		gset_max.addKeyListener(new KeyAdapter() {
			public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				if(keyCode==KeyEvent.VK_V){
					String value=gset_max.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							if(!Pattern.compile("[.E]").matcher(value.toUpperCase()).find()){
								int val=Integer.parseInt(value);
								if(val>0){
									if(val<Integer.parseInt(gset_min.getText())){
										JOptionPane.showMessageDialog(null,"Gene Set Maximum value must be >= Mininum value","Warning",JOptionPane.ERROR_MESSAGE);
										gset_max.grabFocus();
									}else{
										gset_max = new JFormattedTextField(integerformatter);
										gset_max.setValue(val);
										cutoff_apw.grabFocus();
										gset_max.setFocusTraversalKeysEnabled(true);
										valid=false;	
									}
								}else{
									JOptionPane.showMessageDialog(null,"Gene Set Number must be > 0","Warning",JOptionPane.ERROR_MESSAGE);
									gset_max.grabFocus();
								}	
							}else{
								JOptionPane.showMessageDialog(null,"Please put Integer Number","Warning",JOptionPane.ERROR_MESSAGE);
								gset_max.grabFocus();
							}
						}else{
							gset_max.setText(gsetmax);
						}
					}
				}
			}
			
        	public void keyPressed(KeyEvent e) {
        		int keyCode = e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String value=gset_max.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							if(!Pattern.compile("[.E]").matcher(value.toUpperCase()).find()){
								int val=Integer.parseInt(value);
								if(val>0) {
									if(val<Integer.parseInt(gset_min.getText())){
										JOptionPane.showMessageDialog(null,"Gene Set Maximum value must be >= Mininum value","Warning",JOptionPane.ERROR_MESSAGE);
										gset_max.grabFocus();
									}else{
										gset_max = new JFormattedTextField(integerformatter);
										gset_max.setValue(val);
										cutoff_apw.grabFocus();
										gset_max.setFocusTraversalKeysEnabled(true);
										valid=false;
									}
								}else{
									JOptionPane.showMessageDialog(null,"Gene Set Number must be > 0","Warning",JOptionPane.ERROR_MESSAGE);
									gset_max.grabFocus();
								}	
							}else{
								JOptionPane.showMessageDialog(null,"Please put Integer Number","Warning",JOptionPane.ERROR_MESSAGE);
								gset_max.grabFocus();
							}
						}else{
							gset_max.setText(gsetmax);
						}
					}
        		}
			}
		});
		
		cutoff_apw.setComponentPopupMenu(textEditMenu);
		cutoff_apw.setFocusTraversalKeysEnabled(false);
		cutoff_apw.addKeyListener(new KeyAdapter() {
			public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				if(keyCode==KeyEvent.VK_V){
					String value=cutoff_apw.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							double val = Double.parseDouble(value);
							if(val>=0 && val<=1.00){
								cutoff_apw = new JFormattedTextField(decimalformatter);
								cutoff_apw.setValue(val);
								sig_apw.grabFocus();
								cutoff_apw.setFocusTraversalKeysEnabled(true);
								valid=false;
							}else{
								JOptionPane.showMessageDialog(null,"Cut-off value must be >= 0 and <=1 ","Warning",JOptionPane.ERROR_MESSAGE);
								cutoff_apw.grabFocus();
							}
						}else{
							cutoff_apw.setText(cutoffapw);
						}
					}
				}
			}
			
        	public void keyPressed(KeyEvent e) {
        		int keyCode = e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String value=cutoff_apw.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							double val = Double.parseDouble(value);
							if(val>=0 && val<=1.00){
								cutoff_apw = new JFormattedTextField(decimalformatter);
								cutoff_apw.setValue(val);
								sig_apw.grabFocus();
								cutoff_apw.setFocusTraversalKeysEnabled(true);
								valid=false;
							}else{
								JOptionPane.showMessageDialog(null,"Cut-off value must be >= 0 and <=1 ","Warning",JOptionPane.ERROR_MESSAGE);
								cutoff_apw.grabFocus();
							}
						}else{
							cutoff_apw.setText(cutoffapw);
						}
					}
        		}
			}
		});
		
		sig_apw.setComponentPopupMenu(textEditMenu);
		sig_apw.setFocusTraversalKeysEnabled(false);
		sig_apw.addKeyListener(new KeyAdapter() {
			public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				if(keyCode==KeyEvent.VK_V){
					String value=sig_apw.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							double val = Double.parseDouble(value);
							if(val>=0 && val<=1.00){
								if(val>Double.parseDouble(cutoff_apw.getText())){
									JOptionPane.showMessageDialog(null,"Significant value must be <= cutoff","Warning",JOptionPane.ERROR_MESSAGE);
									sig_apw.grabFocus();
								}else{
									sig_apw = new JFormattedTextField(decimalformatter);
									sig_apw.setValue(val);
									mergemethod.grabFocus();
									sig_apw.setFocusTraversalKeysEnabled(true);	
								}
							}else{
								JOptionPane.showMessageDialog(null,"Significant value must be >= 0 and <=1 ","Warning",JOptionPane.ERROR_MESSAGE);
								sig_apw.grabFocus();
							}
						}else{
							sig_apw.setText(sigapw);
						}
					}
				}
			}
			
        	public void keyPressed(KeyEvent e) {
        		int keyCode = e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String value=sig_apw.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							double val = Double.parseDouble(value);
							if(val>=0 && val<=1.00){
								if(val>Double.parseDouble(cutoff_apw.getText())){
									JOptionPane.showMessageDialog(null,"Significant value must be <= cutoff","Warning",JOptionPane.ERROR_MESSAGE);
									sig_apw.grabFocus();
								}else {
									sig_apw = new JFormattedTextField(decimalformatter);
									sig_apw.setValue(val);
									mergemethod.grabFocus();
									sig_apw.setFocusTraversalKeysEnabled(true);
								}
							}else{
								JOptionPane.showMessageDialog(null,"Significant value must be >= 0 and <=1 ","Warning",JOptionPane.ERROR_MESSAGE);
								sig_apw.grabFocus();
							}
						}else{
							sig_apw.setText(sigapw);
						}
					}
        		}
			}
		});
		
		buildButton.addActionListener(new ActionListener(){
	        public void actionPerformed(ActionEvent e){
	        	String inputParams[]= get_userAgrs().split(",");
	        	boolean checkInputs = checkValidInputs(inputParams);
	        	buildButton.setEnabled(false);
    	        buildButton.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
    			mainThread = new Thread() {
	        		public void run() {
	        	    	long t0 = System.currentTimeMillis();
	    	    		if(checkInputs){
	    		        	String cytofiles = inputParams[9]+","+inputParams[10]+","+inputParams[11];
	    		        	genesetfilter = new ArrayList<Integer>();
	    	    			genesetfilter.add(Integer.parseInt(inputParams[2]));
	    	    			genesetfilter.add(Integer.parseInt(inputParams[3]));
	    	    			APW.mainmodule(inputParams[0], inputParams[1], hgmt_names, hgmt_genelist, Double.parseDouble(inputParams[2]), Double.parseDouble(inputParams[3]), inputParams[6], inputParams[7],inputParams[8], cytofiles,Double.parseDouble(inputParams[4]),Double.parseDouble(inputParams[5]),genesetfilter, progressPane, progressBar);
	    		        	if(isValidFile(inputParams[9],"Terms.txt",true)){  //double check
	    		        		runEMCommand(inputParams[12],inputParams[9], inputParams[10],inputParams[11], inputParams[4],inputParams[5],inputParams[13],inputParams[14],Double.parseDouble(inputParams[15]),inputParams[16]);
	    		        	}else{
	    		        		JOptionPane.showMessageDialog(null,"No significant terms were found.\r\nCytoscape files were not written","Warning",JOptionPane.ERROR_MESSAGE);
	    		        	}
	    	        	}
	    	    		repaintBar(progressPane, progressBar, 100);   
	    	    		try {
	    	    			Thread.sleep(TIME_VISIBLE);
	    	    			progressFrame.setVisible(false);
	    	    			mainFrame.setVisible(false);
	    	    			
						} catch (Exception e2) {
							// TODO: handle exception
						}
	    	    		
	    	    		long t1 = System.currentTimeMillis();
	    	    		float diff = (t1 - t0)/1000.0F;
	    	    		JOptionPane pane = new JOptionPane("End of process ("+diff+" seconds)",JOptionPane.INFORMATION_MESSAGE);
	    	    		final JDialog dialog = pane.createDialog(null, "Message");
	    	    		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
	    	    		dialog.setModal(false);
	    	    		dialog.setVisible(true);
	    	    		try {
	    	    			Thread.sleep(TIME_VISIBLE);
	    	    			dialog.setVisible(false);
	    	    			progressbarThread.stop();
    	    				mainThread.stop();
						} catch (Exception e3) {
							// TODO: handle exception
						}
	        		}
	        	};
	        	
	        	progressbarThread = new Thread() {
	        		public void run() {
	        			if(checkInputs){
	        				progressFrame=new JFrame();
		    	        	progressFrame.setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
		    	        	progressFrame.setTitle("Please wait...");
		    	        	progressFrame.setSize(new Dimension(50,20));
		    	            progressBar.setValue(0);
		    	            progressBar.setStringPainted(true);
		    	            progressBar.setPreferredSize(new Dimension(250,25));
		    	            progressPane.add(progressBar);
		    	            
		    	            progressFrame.addWindowListener(new WindowAdapter(){
		    	            	public void windowClosing(WindowEvent we){
		    	            	if(mainThread.isAlive()) {
		    	         		   int option = JOptionPane.showConfirmDialog((Component) null, "Hide progress bar ?","Alert", JOptionPane.YES_NO_CANCEL_OPTION);
		    	        		   if(option==JOptionPane.YES_OPTION){
		    	        			   progressFrame.setVisible(false);
		    	        		   }else {
		    	        			   progressFrame.setVisible(true);
		    	        		   }
		    	   	        	}
		    	   	         }
		    	            });
		    	            
		    	            task = new Task();
		    	            task.addPropertyChangeListener(new PropertyChangeListener() {
		    	            	@Override
		    	        		public void propertyChange(PropertyChangeEvent evt) {
		    	           	        if ("progress" == evt.getPropertyName()) {
		    	           	            int progress = (Integer) evt.getNewValue();
		    	           	            progressBar.setValue(progress);
		    	           	            if(progressBar.getValue()==100){
		    	           	            	progressFrame.setVisible(false);
		    	           	            }
		    	           	        }
		    	        		}
		    	            });
		    	            task.execute();
		    	            progressFrame.getContentPane().add(progressPane);
		    	    		progressFrame.setLocationRelativeTo(null);
		    	    		progressFrame.setVisible(true);
		    	    		progressFrame.pack();
	    	    		}
	        		}
	        	};
	        	
	        	mainThread.start();
	        	progressbarThread.start();
	        	
	        	if(!checkInputs){
    				mainFrame.setVisible(true);	
    				buildButton.setEnabled(true);
    			}
	        }
		});
	
		advanceOptionButton.addActionListener(new ActionListener() { 
        	public void actionPerformed(ActionEvent evt) {
        		mainFrame.setVisible(false);
        		String APW_fin = fin.getText();
	        	String GMT_fin = gmt.getText();
	        	String gsetmin = gset_min.getText();
	        	String gsetmax = gset_max.getText();
	        	String cutoffapw = cutoff_apw.getText();
	        	String sigapw = sig_apw.getText();
        		initComponents_advance(mainFrame,APW_fin,GMT_fin,gsetmin,gsetmax,cutoffapw,sigapw);
        	}
		});
		
		
        resetButton.addActionListener(new ActionListener() { 
        	public void actionPerformed(ActionEvent evt) {
        		fin.setText("");
        		gmt.setText("");
        		gset_min.setText("5");
        		gset_max.setText("1000");
        		cutoff_apw.setText("0.10");
        		sig_apw.setText("0.05");
        		mergemethod.setSelectedIndex(0); 
        		correctmethod.setSelectedIndex(0);
        		returnall.setSelectedIndex(0);
        		fcytotik1.setSelected(true);
        		fcytotik2.setSelected(true);
        		fcytotik3.setSelected(true);
        		fcytoedit1.setText("terms.txt");
        		fcytoedit2.setText("abridged.gmt");
        		fcytoedit3.setText("groups.txt");
        		networkname.setText("Networgk 1");
        		datasetedge.setSelectedIndex(0);
        		cutoff_edge.setValue("0.66");
        		metric.setSelectedIndex(0);
        		metricadjust.setValue(50);
            }
        });
       
        cancelButton.addActionListener(new ActionListener() {
           public void actionPerformed(ActionEvent evt) {
        	   if(mainThread.isAlive()) {
        		   int option = JOptionPane.showConfirmDialog((Component) null, "Stop running APW ?","Alert", JOptionPane.YES_NO_CANCEL_OPTION);
        		   if(option==JOptionPane.YES_OPTION){
        			   mainFrame.dispose();
        			   progressFrame.dispose();
        			   mainThread.stop();
        			   progressbarThread.stop();
        		   }
        	   }else {
        		   int option = JOptionPane.showConfirmDialog((Component) null, "Exit?","Alert", JOptionPane.YES_NO_CANCEL_OPTION);
        		   if(option==JOptionPane.YES_OPTION){
        			   mainFrame.dispose();
        		   }
        	   }
    	   }
       });
       
       fin.grabFocus();
	}
	
	public void initComponents_simple(JXFrame frame,String APW_fin, String GMT_fin, String gsetmin, String gmax, String cutoffapw, String sigapw){
        /*****************************************************************/
	/*  	Repaint the GUI for simple option
	/****************************************************************/
		frame.validate();
		frame.repaint();
		frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		frame.setTitle("ActivePathways Demo V.0.5");
		JPanel contentPane = new JPanel(new GridBagLayout());
		JSplitPane moduleSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
        moduleSplitPane.setBorder(BorderFactory.createEmptyBorder());
        moduleSplitPane.setTopComponent(create_APW_tab_simple());
        moduleSplitPane.setDividerLocation(300);
        contentPane.add(moduleSplitPane, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(11, 11, 12, 12), 0, 0));
        advanceOptionButton = new JButton("Advance");
        buildButton = new JButton(" Build ");
        resetButton = new JButton("Reset");
        cancelButton = new JButton("Cancel");
        advanceOptionButton.setPreferredSize(new Dimension(30,30));
        buildButton.setPreferredSize(new Dimension(30,30));
        resetButton.setPreferredSize(new Dimension(30,30));
        cancelButton.setPreferredSize(new Dimension(30,30));
        JPanel buttonPanel =  new JPanel(new GridBagLayout());
        JPanel buttonPanelRight =  new JPanel(new GridBagLayout());
        JPanel buttonPanelLeft =  new JPanel(new GridBagLayout());
        buttonPanelRight.add(buildButton, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 12, 5, 5), 0, 0));
        buttonPanelRight.add(cancelButton,new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 5, 0), 0, 0));
        buttonPanelLeft.add(resetButton,new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 5, 5), 0, 0));
        buttonPanelLeft.add(advanceOptionButton,new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 5, 5), 0, 0));
        buttonPanel.add(buttonPanelLeft,new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 5, 200), 0, 0));
        buttonPanel.add(buttonPanelRight,new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0, GridBagConstraints.EAST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 5, 0), 0, 0));
        contentPane.add(buttonPanel, new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 5, 12), 0, 0));
        frame.setContentPane(contentPane);
        frame.setSize(400,610);
        frame.setVisible(true);
        frame.setLocationRelativeTo(null);
        frame.pack();
        
        if(!APW_fin.isEmpty()){
        	fin.setText(APW_fin);
        	gmt.setText(GMT_fin);
        	gset_min.setText(gsetmin);
        	gset_max.setText(gsetmax);
        	cutoff_apw.setText(cutoffapw);
        	sig_apw.setText(sigapw);
        	Path path = Paths.get(APW_fin).toAbsolutePath();
			String finpath = path.getParent().toString();
        	String fcyto1=finpath+File.separator+"terms.txt";
        	String fcyto2=finpath+File.separator+"abridged.gmt";
        	String fcyto3=finpath+File.separator+"groups.txt";
        	fcytoedit1.setText(fcyto1);
	      	fcytoedit2.setText(fcyto2);
	      	fcytoedit3.setText(fcyto3);
	      	fcytoedit1.setEnabled(true);
	      	fcytoedit2.setEnabled(true);
	      	fcytoedit3.setEnabled(true);
	      	buildButton.grabFocus();
        }else {
        	buildButton.grabFocus();
        }
        
        fin.grabFocus();
        frame.getRootPane().getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0), "Cancel");
        frame.getRootPane().getActionMap().put("Cancel",new AbstractAction(){ 
        	public void actionPerformed(ActionEvent e){
        		frame.dispose();
        	}
        });
        
        fin.setFocusTraversalKeysEnabled(false);
        fin.addKeyListener(new KeyAdapter() {
        	public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				if(keyCode==KeyEvent.VK_V){
					String filename=fin.getText();
					if(!filename.equals("")){
						boolean valid = isValidFile(filename,"",true);
						if(valid){
							Path path = Paths.get(filename).toAbsolutePath();
							String finpath = path.getParent().toString();
				        	String fcyto1=finpath+File.separator+"terms.txt";
				        	String fcyto2=finpath+File.separator+"abridged.gmt";
				        	String fcyto3=finpath+File.separator+"groups.txt";
			    	      	jbfin.grabFocus();
						}
					}
				}
			}
		
        	public void keyPressed(KeyEvent e){
        		int keyCode= e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String filename=fin.getText();
        			fin.setText(filename);
        			if(!filename.equals("")){
            			boolean valid = isValidFile(filename,"",true);
            			if(valid){
            				File path = new File(filename);
            	        	String finpath=path.getParent();
            	        	String fcyto1=finpath+File.separator+"terms.txt";
            	        	String fcyto2=finpath+File.separator+"abridged.gmt";
            	        	String fcyto3=finpath+File.separator+"groups.txt";
                	      	jbfin.grabFocus();	
            			}
            		}else{
            			JOptionPane.showMessageDialog(null,"Please select a Valid Input File","Warning",JOptionPane.ERROR_MESSAGE);
            		}
        		}
			}
		});
        jbfin.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
	        	JFileChooser fileopen = new JFileChooser();
	    	    FileFilter filter = new FileNameExtensionFilter("csv files", "csv");
	    	    fileopen.addChoosableFileFilter(filter);
	    	    int ret = fileopen.showDialog(null, "Import");
	    	    if (ret == JFileChooser.APPROVE_OPTION) {
	    	    	File file = fileopen.getSelectedFile();
	    	      	String filename = file.getPath();
	    	      	fin.setText(filename);
	    	    }else{
	        }  	
	        }
	    });
        
        gmt.setFocusTraversalKeysEnabled(false);
		gmt.addKeyListener(new KeyAdapter() {
        	public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				switch(keyCode){
				case KeyEvent.VK_V:
					if(!gmt.getText().equals("")){
						String filename=gmt.getText();
						isValidFile(filename,"",true);
					}
				}
			}
			
        	public void keyPressed(KeyEvent e) {
        		int keyCode= e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String filename=gmt.getText();
        			gmt.setText(filename);
            		if(!filename.equals("")){
            			boolean valid = isValidFile(filename,"",true);
            			if(valid){
            				jbgmt.grabFocus();	
            			}
            		}else{
            			JOptionPane.showMessageDialog(null,"Please select a Valid GMT File","Warning",JOptionPane.ERROR_MESSAGE);
            		}
        		}
			}
		});
	
		jbgmt.addActionListener(new ActionListener(){
	        public void actionPerformed(ActionEvent e) {
	        	JFileChooser fileopen = new JFileChooser();
	    	    FileFilter filter = new FileNameExtensionFilter("gmt files", "gmt");
	    	    fileopen.addChoosableFileFilter(filter);
	    	    int ret = fileopen.showDialog(null, "Import");
	    	    if (ret == JFileChooser.APPROVE_OPTION) {
	    	    	File file = fileopen.getSelectedFile();
	    	      	String filename = file.getPath();
	    	      	gmt.setText(filename);
	        }  	
	        }
	    });
		
		gset_min.setFocusTraversalKeysEnabled(false);
		gset_min.addKeyListener(new KeyAdapter() {
			public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				if(keyCode==KeyEvent.VK_V){
					String value=gset_min.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							if(!Pattern.compile("[.E]").matcher(value.toUpperCase()).find()){
								int val=Integer.parseInt(value);
								if(val>0){
									gset_min = new JFormattedTextField(integerformatter);
									gset_min.setValue(val);
									gset_max.grabFocus();
									valid=false;
									gset_min.setFocusTraversalKeysEnabled(true);
								}else{
									JOptionPane.showMessageDialog(null,"Gene Set Number must be > 0","Warning",JOptionPane.ERROR_MESSAGE);
									gset_min.grabFocus();
								}	
							}else{
								JOptionPane.showMessageDialog(null,"Please put Integer Number","Warning",JOptionPane.ERROR_MESSAGE);
								gset_min.grabFocus();
							}
						}else{
							gset_min.setText(gsetmin);
						}
					}
				}
			}
			
        	public void keyPressed(KeyEvent e) {
        		int keyCode = e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String value=gset_min.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							if(!Pattern.compile("[.E]").matcher(value.toUpperCase()).find()){
								int val=Integer.parseInt(value);
								if(val>0){
									gset_min = new JFormattedTextField(integerformatter);
									gset_min.setValue(val);
									gset_max.grabFocus();
									gset_min.setFocusTraversalKeysEnabled(true);
									valid=false;	
								}else{
									JOptionPane.showMessageDialog(null,"Gene Set Number must be > 0","Warning",JOptionPane.ERROR_MESSAGE);
									gset_min.grabFocus();
								}	
							}else{
								JOptionPane.showMessageDialog(null,"Please put Integer Number","Warning",JOptionPane.ERROR_MESSAGE);
								gset_min.grabFocus();
							}
						}else{
							gset_min.setText(gsetmin);
						}
					}
        		}
			}
		});
		
		gset_max.setFocusTraversalKeysEnabled(false);
		gset_max.addKeyListener(new KeyAdapter() {
			public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				if(keyCode==KeyEvent.VK_V){
					String value=gset_max.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							if(!Pattern.compile("[.E]").matcher(value.toUpperCase()).find()){
								int val=Integer.parseInt(value);
								if(val>0){
									if(val<Integer.parseInt(gset_min.getText())){
										JOptionPane.showMessageDialog(null,"Gene Set Maximum value must be >= Mininum value","Warning",JOptionPane.ERROR_MESSAGE);
										gset_max.grabFocus();
									}else{
										gset_max = new JFormattedTextField(integerformatter);
										gset_max.setValue(val);
										cutoff_apw.grabFocus();
										gset_max.setFocusTraversalKeysEnabled(true);
										valid=false;	
									}
								}else{
									JOptionPane.showMessageDialog(null,"Gene Set Number must be > 0","Warning",JOptionPane.ERROR_MESSAGE);
									gset_max.grabFocus();
								}	
							}else{
								JOptionPane.showMessageDialog(null,"Please put Integer Number","Warning",JOptionPane.ERROR_MESSAGE);
								gset_max.grabFocus();
							}
						}else{
							gset_max.setText(gsetmax);
						}
					}
				}
			}
			
        	public void keyPressed(KeyEvent e) {
        		int keyCode = e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String value=gset_max.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							if(!Pattern.compile("[.E]").matcher(value.toUpperCase()).find()){
								int val=Integer.parseInt(value);
								if(val>0) {
									if(val<Integer.parseInt(gset_min.getText())){
										JOptionPane.showMessageDialog(null,"Gene Set Maximum value must be >= Mininum value","Warning",JOptionPane.ERROR_MESSAGE);
										gset_max.grabFocus();
									}else{
										gset_max = new JFormattedTextField(integerformatter);
										gset_max.setValue(val);
										cutoff_apw.grabFocus();
										gset_max.setFocusTraversalKeysEnabled(true);
										valid=false;
									}
								}else{
									JOptionPane.showMessageDialog(null,"Gene Set Number must be > 0","Warning",JOptionPane.ERROR_MESSAGE);
									gset_max.grabFocus();
								}	
							}else{
								JOptionPane.showMessageDialog(null,"Please put Integer Number","Warning",JOptionPane.ERROR_MESSAGE);
								gset_max.grabFocus();
							}
						}else{
							gset_max.setText(gsetmax);
						}
					}
        		}
			}
		});
		
		cutoff_apw.setFocusTraversalKeysEnabled(false);
		cutoff_apw.addKeyListener(new KeyAdapter() {
			public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				if(keyCode==KeyEvent.VK_V){
					String value=cutoff_apw.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							double val = Double.parseDouble(value);
							if(val>=0 && val<=1.00){
								cutoff_apw = new JFormattedTextField(decimalformatter);
								cutoff_apw.setValue(val);
								sig_apw.grabFocus();
								cutoff_apw.setFocusTraversalKeysEnabled(true);
								valid=false;
							}else{
								JOptionPane.showMessageDialog(null,"Cut-off value must be >= 0 and <=1 ","Warning",JOptionPane.ERROR_MESSAGE);
								cutoff_apw.grabFocus();
							}
						}else{
							cutoff_apw.setText(cutoffapw);
						}
					}
				}
			}
			
        	public void keyPressed(KeyEvent e) {
        		int keyCode = e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String value=cutoff_apw.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							double val = Double.parseDouble(value);
							if(val>=0 && val<=1.00){
								cutoff_apw = new JFormattedTextField(decimalformatter);
								cutoff_apw.setValue(val);
								sig_apw.grabFocus();
								cutoff_apw.setFocusTraversalKeysEnabled(true);
								valid=false;
							}else{
								JOptionPane.showMessageDialog(null,"Cut-off value must be >= 0 and <=1 ","Warning",JOptionPane.ERROR_MESSAGE);
								cutoff_apw.grabFocus();
							}
						}else{
							cutoff_apw.setText(cutoffapw);
						}
					}
        		}
			}
		});
		
		sig_apw.setFocusTraversalKeysEnabled(false);
		sig_apw.addKeyListener(new KeyAdapter() {
			public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				if(keyCode==KeyEvent.VK_V){
					String value=sig_apw.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							double val = Double.parseDouble(value);
							if(val>=0 && val<=1.00){
								if(val>Double.parseDouble(cutoff_apw.getText())){
									JOptionPane.showMessageDialog(null,"Significant value must be <= cutoff","Warning",JOptionPane.ERROR_MESSAGE);
									sig_apw.grabFocus();
								}else{
									sig_apw = new JFormattedTextField(decimalformatter);
									sig_apw.setValue(val);
									mergemethod.grabFocus();
									sig_apw.setFocusTraversalKeysEnabled(true);	
								}
							}else{
								JOptionPane.showMessageDialog(null,"Significant value must be >= 0 and <=1 ","Warning",JOptionPane.ERROR_MESSAGE);
								sig_apw.grabFocus();
							}
						}else{
							sig_apw.setText(sigapw);
						}
					}
				}
			}
			
        	public void keyPressed(KeyEvent e) {
        		int keyCode = e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String value=sig_apw.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							double val = Double.parseDouble(value);
							if(val>=0 && val<=1.00){
								if(val>Double.parseDouble(cutoff_apw.getText())){
									JOptionPane.showMessageDialog(null,"Significant value must be <= cutoff","Warning",JOptionPane.ERROR_MESSAGE);
									sig_apw.grabFocus();
								}else {
									sig_apw = new JFormattedTextField(decimalformatter);
									sig_apw.setValue(val);
									mergemethod.grabFocus();
									sig_apw.setFocusTraversalKeysEnabled(true);
								}
							}else{
								JOptionPane.showMessageDialog(null,"Significant value must be >= 0 and <=1 ","Warning",JOptionPane.ERROR_MESSAGE);
								sig_apw.grabFocus();
							}
						}else{
							sig_apw.setText(sigapw);
						}
					}
        		}
			}
		});
		
		buildButton.addActionListener(new ActionListener(){
	        public void actionPerformed(ActionEvent e){
	        	mainThread = new Thread() {
	        		public void run() {
	    	        	String APW_fin = fin.getText();
	    	        	String GMT_fin = gmt.getText();
	    	        	String gsetmin = gset_min.getText();
	    	        	String gsetmax = gset_max.getText();
	    	        	String cutoffapw = cutoff_apw.getText();
	    	        	String sigapw = sig_apw.getText();
	    	        	String merge="Brown";
	    	        	String correct="fdr";
	    	        	String returnAPW="no";
	    	        	
	    	        	Path path = Paths.get(APW_fin).toAbsolutePath();
						String finpath = path.getParent().toString();
			        	String fcyto1=finpath+File.separator+"terms.txt";
			        	String fcyto2=finpath+File.separator+"abridged.gmt";
			        	String fcyto3=finpath+File.separator+"groups.txt";
		    	      			
	    	        	String mapName=APW_fin.split(File.separator)[APW_fin.split(File.separator).length-1].split("\\.")[0];
	    	        	String edgeDataset= "AUTOMATIC";
	    	        	String nodeFDR = "0.05";
	    	        	String nodePval = "0.05";
	    	        	String edgeCutoff = "0.66";
	    	        	String edgeMetric= "COMBINED";
	    	        	Double combinedConstant=0.50;
	    	        	String edgeStrategy = "AUTOMATIC";
	    	        	String inputParams[]= {APW_fin,GMT_fin,gsetmin,gsetmax,cutoffapw,sigapw,merge,correct,returnAPW,fcyto1,fcyto2,fcyto3,mapName,edgeDataset,nodeFDR,nodePval,edgeCutoff,edgeCutoff};
	    	        	long t0 = System.currentTimeMillis();
	    	    		if(checkValidInputs(inputParams)){
	    		        	String cytofiles = fcyto1+","+fcyto2+","+fcyto3;
	    		        	genesetfilter = new ArrayList<Integer>();
	    	    			genesetfilter.add(Integer.parseInt(gsetmin));
	    	    			genesetfilter.add(Integer.parseInt(gsetmax));
	    	    			APW.mainmodule(inputParams[0], inputParams[1], hgmt_names, hgmt_genelist, Double.parseDouble(inputParams[2]), Double.parseDouble(inputParams[3]), inputParams[6], inputParams[7],inputParams[8], cytofiles,Double.parseDouble(inputParams[4]),Double.parseDouble(inputParams[5]),genesetfilter, progressPane, progressBar);
	    	    			if(isValidFile(fcyto1,"Terms.txt",true)){  //double check
	    		        		runEMCommand(mapName,fcyto1, fcyto2,fcyto3, nodePval,nodeFDR,edgeCutoff,edgeMetric,combinedConstant,edgeStrategy);
	    		        	}else{
	    		        		JOptionPane.showMessageDialog(null,"No significant terms were found.\r\nCytoscape files were not written","Warning",JOptionPane.ERROR_MESSAGE);
	    		        	}
	    	        	}
	    	    		
	    	    		long t1 = System.currentTimeMillis();
	    	    		float diff = (t1 - t0)/1000.0F;
	    	    		JOptionPane pane = new JOptionPane("End of process ("+diff+" seconds)",JOptionPane.INFORMATION_MESSAGE);
	    	    		final JDialog dialog = pane.createDialog(null, "Message");
	    	    		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		    	        dialog.setModal(false);
	    	    		dialog.setVisible(true);
	    	    		new Timer(TIME_VISIBLE, new ActionListener() {
	    	    			@Override
	    	    			public void actionPerformed(ActionEvent e) {
	    	    				dialog.setVisible(false);
	    	    			}
	    	    		}).start();
	    	    		}
	        	};
	        	
	        	progressbarThread = new Thread() {
	        		public void run() {
	        			JFrame progressFrame = new JFrame();
	    	        	progressFrame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
	    	        	progressFrame.setTitle("Please wait...");
	    	        	progressFrame.setSize(new Dimension(50,20));
	    	        	JPanel progressPane = new JPanel();
	    	    		progressPane.setLayout(new BorderLayout());
	    	    		progressImageLabel=new JLabel();
	    	    		progressImage = new ImageIcon(getClass().getResource(gif));
	    	    		progressImageLabel.setIcon(progressImage);
	    	    		progressPane.add(progressImageLabel,BorderLayout.CENTER);
	    	    		progressFrame.getContentPane().add(progressPane);
	    	    		progressFrame.setLocationRelativeTo(null);
	    	    		progressFrame.setVisible(true);
	    	    		progressFrame.pack();
	        		}
	        	};
	        	mainThread.start();
	        }
		});
		
		
		advanceOptionButton.addActionListener(new ActionListener() { 
        	public void actionPerformed(ActionEvent evt) {
        		frame.setVisible(false);
        		String APW_fin = fin.getText();
	        	String GMT_fin = gmt.getText();
	        	String gsetmin = gset_min.getText();
	        	String gsetmax = gset_max.getText();
	        	String cutoffapw = cutoff_apw.getText();
	        	String sigapw = sig_apw.getText();
        		initComponents_advance(frame,APW_fin,GMT_fin,gsetmin,gsetmax,cutoffapw,sigapw);
        	}
		});
		
		
        resetButton.addActionListener(new ActionListener() { 
        	public void actionPerformed(ActionEvent evt) {
        		fin.setText("");
        		gmt.setText("");
        		gset_min.setText("5");
        		gset_max.setText("1000");
        		cutoff_apw.setText("0.10");
        		sig_apw.setText("0.05");
        		mergemethod.setSelectedIndex(0); 
        		correctmethod.setSelectedIndex(0);
        		returnall.setSelectedIndex(0);
        		fcytotik1.setSelected(true);
        		fcytotik2.setSelected(true);
        		fcytotik3.setSelected(true);
        		fcytoedit1.setText("terms.txt");
        		fcytoedit2.setText("abridged.gmt");
        		fcytoedit3.setText("groups.txt");
        		networkname.setText("Networgk 1");
        		datasetedge.setSelectedIndex(0);
        		cutoff_edge.setValue("0.66");
        		metric.setSelectedIndex(0);
        		metricadjust.setValue(50);
            }
        });
       
        cancelButton.addActionListener(new ActionListener() {
           public void actionPerformed(ActionEvent evt) {
    		   int option = JOptionPane.showConfirmDialog((Component) null, "Exit?","Alert", JOptionPane.YES_NO_CANCEL_OPTION);
    		   if(option==JOptionPane.YES_OPTION){
    			   frame.dispose();
    		   }
    	   }
       });
       fin.grabFocus();
	}
	
	private void initComponents_advance(JXFrame frame,String APW_fin, String GMT_fin, String gsetmin, String gmax, String cutoffapw, String sigapw) {
        /*****************************************************************/
	/*  	Repaint the GUI for advance option
	/****************************************************************/
		frame.getContentPane().validate();
		frame.getContentPane().repaint();
		frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		frame.setTitle("ActivePathways Demo V.0.5");
		JPanel contentPane = new JPanel(new GridBagLayout());
        JSplitPane moduleSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
        moduleSplitPane.setBorder(BorderFactory.createEmptyBorder());
        moduleSplitPane.setTopComponent(create_APW_tab());
        moduleSplitPane.setDividerLocation(300);
        contentPane.add(moduleSplitPane, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(11, 11, 12, 12), 0, 0));
        moduleSplitPane.setBottomComponent(create_ENR_tab());
        simpleOptionButton = new JButton("Simple");
        buildButton = new JButton(" Build ");
        resetButton = new JButton("Reset");
        cancelButton = new JButton("Cancel");
        simpleOptionButton.setPreferredSize(new Dimension(30,30));
        buildButton.setPreferredSize(new Dimension(30,30));
        resetButton.setPreferredSize(new Dimension(30,30));
        cancelButton.setPreferredSize(new Dimension(30,30));
        JPanel buttonPanel =  new JPanel(new GridBagLayout());
        JPanel buttonPanelRight =  new JPanel(new GridBagLayout());
        JPanel buttonPanelLeft =  new JPanel(new GridBagLayout());
        buttonPanelRight.add(buildButton, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 12, 5, 5), 0, 0));
        buttonPanelRight.add(cancelButton,new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 5, 0), 0, 0));
        buttonPanelLeft.add(resetButton,new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 5, 5), 0, 0));
        buttonPanelLeft.add(simpleOptionButton,new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 5, 5), 0, 0));
        buttonPanel.add(buttonPanelLeft,new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 5, 200), 0, 0));
        buttonPanel.add(buttonPanelRight,new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0, GridBagConstraints.EAST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 5, 0), 0, 0));
        contentPane.add(buttonPanel, new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 5, 12), 0, 0));
        frame.setContentPane(contentPane);
        frame.setSize(700,610);
        frame.setVisible(true);
        frame.setLocationRelativeTo(null);
        frame.pack();
        
        /*****************************************************************/
	/*  		Add ActionListener to all components
	/****************************************************************/
        if(!APW_fin.isEmpty()){
        	fin.setText(APW_fin);
        	gmt.setText(GMT_fin);
        	gset_min.setText(gsetmin);
        	gset_max.setText(gsetmax);
        	cutoff_apw.setText(cutoffapw);
        	sig_apw.setText(sigapw);
        	Path path = Paths.get(APW_fin).toAbsolutePath();
			String finpath = path.getParent().toString();
        	String fcyto1=finpath+File.separator+"terms.txt";
        	String fcyto2=finpath+File.separator+"abridged.gmt";
        	String fcyto3=finpath+File.separator+"groups.txt";
        	fcytoedit1.setText(fcyto1);
	      	fcytoedit2.setText(fcyto2);
	      	fcytoedit3.setText(fcyto3);
	      	fcytoedit1.setEnabled(true);
	      	fcytoedit2.setEnabled(true);
	      	fcytoedit3.setEnabled(true);
        	networkname.grabFocus();
        }else {
        	fin.grabFocus();
        }
        
        frame.getRootPane().getInputMap(JComponent.WHEN_IN_FOCUSED_WINDOW).put(KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0), "Cancel");
        frame.getRootPane().getActionMap().put("Cancel",new AbstractAction(){ 
        	public void actionPerformed(ActionEvent e){
        		frame.dispose();
        	}
        });
        
        fin.setFocusTraversalKeysEnabled(false);
        fin.addKeyListener(new KeyAdapter() {
        	public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				if(keyCode==KeyEvent.VK_V){
					String filename=fin.getText();
					String mapName=filename.split(File.separator)[filename.split(File.separator).length-1].split("\\.")[0];
		        	networkname.setText(mapName);
					if(!filename.equals("")){
						boolean valid = isValidFile(filename,"",true);
						if(valid){
							Path path = Paths.get(filename).toAbsolutePath();
							String finpath = path.getParent().toString();
				        	String fcyto1=finpath+File.separator+"terms.txt";
				        	String fcyto2=finpath+File.separator+"abridged.gmt";
				        	String fcyto3=finpath+File.separator+"groups.txt";
			    	      	fcytoedit1.setText(fcyto1);
			    	      	fcytoedit2.setText(fcyto2);
			    	      	fcytoedit3.setText(fcyto3);
			    	      	fcytoedit1.setEnabled(true);
			    	      	fcytoedit2.setEnabled(true);
			    	      	fcytoedit3.setEnabled(true);
			    	      	jbfin.grabFocus();
						}
					}
				}
			}
		
        	public void keyPressed(KeyEvent e){
        		int keyCode= e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String filename=fin.getText();
        			fin.setText(filename);
        			String mapName=filename.split(File.separator)[filename.split(File.separator).length-1].split("\\.")[0];
		        	networkname.setText(mapName);
					
            		if(!filename.equals("")){
            			boolean valid = isValidFile(filename,"",true);
            			if(valid){
            				File path = new File(filename);
            	        	String finpath=path.getParent();
            	        	String fcyto1=finpath+File.separator+"terms.txt";
            	        	String fcyto2=finpath+File.separator+"abridged.gmt";
            	        	String fcyto3=finpath+File.separator+"groups.txt";
                	      	fcytoedit1.setText(fcyto1);
                	      	fcytoedit2.setText(fcyto2);
                	      	fcytoedit3.setText(fcyto3);
                	      	fcytoedit1.setEnabled(true);
			    	      	fcytoedit2.setEnabled(true);
			    	      	fcytoedit3.setEnabled(true);
        					jbfin.grabFocus();	
            			}
            		}else{
            			JOptionPane.showMessageDialog(null,"Please select a Valid Input File","Warning",JOptionPane.ERROR_MESSAGE);
            		}
        		}
			}
		});
        jbfin.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
	        	JFileChooser fileopen = new JFileChooser();
	    	    FileFilter filter = new FileNameExtensionFilter("csv files", "csv");
	    	    fileopen.addChoosableFileFilter(filter);
	    	    int ret = fileopen.showDialog(null, "Import");
	    	    if (ret == JFileChooser.APPROVE_OPTION) {
	    	    	File file = fileopen.getSelectedFile();
	    	      	String filename = file.getPath();
	    	      	fin.setText(filename);
	    	      	String mapName=filename.split(File.separator)[filename.split(File.separator).length-1].split("\\.")[0];
		        	networkname.setText(mapName);
		        }  	
	        }
	    });
        
        gmt.setFocusTraversalKeysEnabled(false);
		gmt.addKeyListener(new KeyAdapter() {
        	public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				switch(keyCode){
				case KeyEvent.VK_V:
					if(!gmt.getText().equals("")){
						String filename=gmt.getText();
						isValidFile(filename,"",true);
					}
				}
			}
			
        	public void keyPressed(KeyEvent e) {
        		int keyCode= e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String filename=gmt.getText();
        			gmt.setText(filename);
            		if(!filename.equals("")){
            			boolean valid = isValidFile(filename,"",true);
            			if(valid){
            				jbgmt.grabFocus();	
            			}
            		}else{
            			JOptionPane.showMessageDialog(null,"Please select a Valid GMT File","Warning",JOptionPane.ERROR_MESSAGE);
            		}
        		}
			}
		});
	
		jbgmt.addActionListener(new ActionListener(){
	        public void actionPerformed(ActionEvent e) {
	        	JFileChooser fileopen = new JFileChooser();
	    	    FileFilter filter = new FileNameExtensionFilter("gmt files", "gmt");
	    	    fileopen.addChoosableFileFilter(filter);
	    	    int ret = fileopen.showDialog(null, "Import");
	    	    if (ret == JFileChooser.APPROVE_OPTION) {
	    	    	File file = fileopen.getSelectedFile();
	    	      	String filename = file.getPath();
	    	      	gmt.setText(filename);
		        }  	
	        }
	    });
		
		gset_min.setFocusTraversalKeysEnabled(false);
		gset_min.addKeyListener(new KeyAdapter() {
			public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				if(keyCode==KeyEvent.VK_V){
					String value=gset_min.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							if(!Pattern.compile("[.E]").matcher(value.toUpperCase()).find()){
								int val=Integer.parseInt(value);
								if(val>0){
									gset_min = new JFormattedTextField(integerformatter);
									gset_min.setValue(val);
									gset_max.grabFocus();
									valid=false;
									gset_min.setFocusTraversalKeysEnabled(true);
								}else{
									JOptionPane.showMessageDialog(null,"Gene Set Number must be > 0","Warning",JOptionPane.ERROR_MESSAGE);
									gset_min.grabFocus();
								}	
							}else{
								JOptionPane.showMessageDialog(null,"Please put Integer Number","Warning",JOptionPane.ERROR_MESSAGE);
								gset_min.grabFocus();
							}
						}else{
							gset_min.setText(gsetmin);
						}
					}
				}
			}
			
        	public void keyPressed(KeyEvent e) {
        		int keyCode = e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String value=gset_min.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							if(!Pattern.compile("[.E]").matcher(value.toUpperCase()).find()){
								int val=Integer.parseInt(value);
								if(val>0){
									gset_min = new JFormattedTextField(integerformatter);
									gset_min.setValue(val);
									gset_max.grabFocus();
									gset_min.setFocusTraversalKeysEnabled(true);
									valid=false;	
								}else{
									JOptionPane.showMessageDialog(null,"Gene Set Number must be > 0","Warning",JOptionPane.ERROR_MESSAGE);
									gset_min.grabFocus();
								}	
							}else{
								JOptionPane.showMessageDialog(null,"Please put Integer Number","Warning",JOptionPane.ERROR_MESSAGE);
								gset_min.grabFocus();
							}
						}else{
							gset_min.setText(gsetmin);
						}
					}
        		}
			}
		});
		
		gset_max.setFocusTraversalKeysEnabled(false);
		gset_max.addKeyListener(new KeyAdapter() {
			public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				if(keyCode==KeyEvent.VK_V){
					String value=gset_max.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							if(!Pattern.compile("[.E]").matcher(value.toUpperCase()).find()){
								int val=Integer.parseInt(value);
								if(val>0){
									if(val<Integer.parseInt(gset_min.getText())){
										JOptionPane.showMessageDialog(null,"Gene Set Maximum value must be >= Mininum value","Warning",JOptionPane.ERROR_MESSAGE);
										gset_max.grabFocus();
									}else{
										gset_max = new JFormattedTextField(integerformatter);
										gset_max.setValue(val);
										cutoff_apw.grabFocus();
										gset_max.setFocusTraversalKeysEnabled(true);
										valid=false;	
									}
								}else{
									JOptionPane.showMessageDialog(null,"Gene Set Number must be > 0","Warning",JOptionPane.ERROR_MESSAGE);
									gset_max.grabFocus();
								}	
							}else{
								JOptionPane.showMessageDialog(null,"Please put Integer Number","Warning",JOptionPane.ERROR_MESSAGE);
								gset_max.grabFocus();
							}
						}else{
							gset_max.setText(gsetmax);
						}
					}
				}
			}
			
        	public void keyPressed(KeyEvent e) {
        		int keyCode = e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String value=gset_max.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							if(!Pattern.compile("[.E]").matcher(value.toUpperCase()).find()){
								int val=Integer.parseInt(value);
								if(val>0) {
									if(val<Integer.parseInt(gset_min.getText())){
										JOptionPane.showMessageDialog(null,"Gene Set Maximum value must be >= Mininum value","Warning",JOptionPane.ERROR_MESSAGE);
										gset_max.grabFocus();
									}else{
										gset_max = new JFormattedTextField(integerformatter);
										gset_max.setValue(val);
										cutoff_apw.grabFocus();
										gset_max.setFocusTraversalKeysEnabled(true);
										valid=false;
									}
								}else{
									JOptionPane.showMessageDialog(null,"Gene Set Number must be > 0","Warning",JOptionPane.ERROR_MESSAGE);
									gset_max.grabFocus();
								}	
							}else{
								JOptionPane.showMessageDialog(null,"Please put Integer Number","Warning",JOptionPane.ERROR_MESSAGE);
								gset_max.grabFocus();
							}
						}else{
							gset_max.setText(gsetmax);
						}
					}
        		}
			}
		});
		
		cutoff_apw.setFocusTraversalKeysEnabled(false);
		cutoff_apw.addKeyListener(new KeyAdapter() {
			public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				if(keyCode==KeyEvent.VK_V){
					String value=cutoff_apw.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							double val = Double.parseDouble(value);
							if(val>=0 && val<=1.00){
								cutoff_apw = new JFormattedTextField(decimalformatter);
								cutoff_apw.setValue(val);
								sig_apw.grabFocus();
								cutoff_apw.setFocusTraversalKeysEnabled(true);
								valid=false;
							}else{
								JOptionPane.showMessageDialog(null,"Cut-off value must be >= 0 and <=1 ","Warning",JOptionPane.ERROR_MESSAGE);
								cutoff_apw.grabFocus();
							}
						}else{
							cutoff_apw.setText(cutoffapw);
						}
					}
				}
			}
			
        	public void keyPressed(KeyEvent e) {
        		int keyCode = e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String value=cutoff_apw.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							double val = Double.parseDouble(value);
							if(val>=0 && val<=1.00){
								cutoff_apw = new JFormattedTextField(decimalformatter);
								cutoff_apw.setValue(val);
								sig_apw.grabFocus();
								cutoff_apw.setFocusTraversalKeysEnabled(true);
								valid=false;
							}else{
								JOptionPane.showMessageDialog(null,"Cut-off value must be >= 0 and <=1 ","Warning",JOptionPane.ERROR_MESSAGE);
								cutoff_apw.grabFocus();
							}
						}else{
							cutoff_apw.setText(cutoffapw);
						}
					}
        		}
			}
		});
		
		sig_apw.setFocusTraversalKeysEnabled(false);
		sig_apw.addKeyListener(new KeyAdapter() {
			public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				if(keyCode==KeyEvent.VK_V){
					String value=sig_apw.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							double val = Double.parseDouble(value);
							if(val>=0 && val<=1.00){
								if(val>Double.parseDouble(cutoff_apw.getText())){
									JOptionPane.showMessageDialog(null,"Significant value must be <= cutoff","Warning",JOptionPane.ERROR_MESSAGE);
									sig_apw.grabFocus();
								}else{
									sig_apw = new JFormattedTextField(decimalformatter);
									sig_apw.setValue(val);
									mergemethod.grabFocus();
									sig_apw.setFocusTraversalKeysEnabled(true);	
								}
							}else{
								JOptionPane.showMessageDialog(null,"Significant value must be >= 0 and <=1 ","Warning",JOptionPane.ERROR_MESSAGE);
								sig_apw.grabFocus();
							}
						}else{
							sig_apw.setText(sigapw);
						}
					}
				}
			}
			
        	public void keyPressed(KeyEvent e) {
        		int keyCode = e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String value=sig_apw.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							double val = Double.parseDouble(value);
							if(val>=0 && val<=1.00){
								if(val>Double.parseDouble(cutoff_apw.getText())){
									JOptionPane.showMessageDialog(null,"Significant value must be <= cutoff","Warning",JOptionPane.ERROR_MESSAGE);
									sig_apw.grabFocus();
								}else {
									sig_apw = new JFormattedTextField(decimalformatter);
									sig_apw.setValue(val);
									mergemethod.grabFocus();
									sig_apw.setFocusTraversalKeysEnabled(true);
								}
							}else{
								JOptionPane.showMessageDialog(null,"Significant value must be >= 0 and <=1 ","Warning",JOptionPane.ERROR_MESSAGE);
								sig_apw.grabFocus();
							}
						}else{
							sig_apw.setText(sigapw);
						}
					}
        		}
			}
		});
		
		
		fcytotik1.addActionListener(new ActionListener(){
			 public void actionPerformed(ActionEvent e) {
		       if(fcytotik1.isSelected()){
		    	   fcytoedit1.setEditable(true);
		    	  if(fcytoedit1.getText().equals("")){
		    		  fcytoedit1.setText("terms.txt");
				  }else{
					  String f=fcytoedit1.getText();	
				  }
		       }else{
		    	   fcytoedit1.setText("");
		       }
			 }
		});
		
		fcytotik2.addActionListener(new ActionListener(){
			 public void actionPerformed(ActionEvent e) {
		       if(fcytotik2.isSelected()){
		    	   fcytoedit2.setEditable(true);
		    	   if(fcytoedit2.getText().equals("")){
		    		   	fcytoedit2.setText("abridged.gmt");
					  }else{
						  String f=fcytoedit2.getText();	
					  }
		       }else{
		    	   fcytoedit2.setText("");
		       }
			 }
		});
		
		fcytotik3.addActionListener(new ActionListener(){
			 public void actionPerformed(ActionEvent e) {
		       if(fcytotik3.isSelected()){
		    	   fcytoedit3.setEditable(true);
		    	   if(fcytoedit3.getText().equals("")){
		    		   fcytoedit3.setText("groups.txt");
		    	   }else{
				      String f=fcytoedit3.getText();	
		    	   }
		       }else{
		    	   fcytoedit3.setText("");
		       }
			 }
		});
		
		fcytosave1.addActionListener(new ActionListener(){
	        public void actionPerformed(ActionEvent e) {
	        	JFileChooser fileopen = new JFileChooser();
	    	    FileFilter filter = new FileNameExtensionFilter("txt files", "txt");
	    	    fileopen.addChoosableFileFilter(filter);
	    	    int ret = fileopen.showDialog(null, "Select");
	    	    if (ret == JFileChooser.APPROVE_OPTION) {
	    	    	File file = fileopen.getSelectedFile();
	    	      	String filename = file.getPath();
	    	      	fcytoedit1.setText(filename);
	        }  	
	        }
	    });
		
		fcytosave2.addActionListener(new ActionListener(){
	        public void actionPerformed(ActionEvent e) {
	        	JFileChooser fileopen = new JFileChooser();
	    	    FileFilter filter = new FileNameExtensionFilter("txt files", "txt");
	    	    fileopen.addChoosableFileFilter(filter);
	    	    int ret = fileopen.showDialog(null, "Select");
	    	    if (ret == JFileChooser.APPROVE_OPTION) {
	    	    	File file = fileopen.getSelectedFile();
	    	      	String filename = file.getPath();
	    	      	fcytoedit2.setText(filename);
	        }  	
	        }
	    });
		
		fcytosave3.addActionListener(new ActionListener(){
	        public void actionPerformed(ActionEvent e) {
	        	JFileChooser fileopen = new JFileChooser();
	    	    FileFilter filter = new FileNameExtensionFilter("txt files", "txt");
	    	    fileopen.addChoosableFileFilter(filter);
	    	    int ret = fileopen.showDialog(null, "Select");
	    	    if (ret == JFileChooser.APPROVE_OPTION) {
	    	    	File file = fileopen.getSelectedFile();
	    	      	String filename = file.getPath();
	    	      	fcytoedit3.setText(filename);
	        }  	
	        }
	    });
		
		networkname.addActionListener(new ActionListener(){
        	public void actionPerformed(ActionEvent e) {
        		networkname.setText(networkname.getText());
        	}
		});
		

		cutoff_edge.grabFocus();
		cutoff_edge.setFocusTraversalKeysEnabled(false);
		cutoff_edge.addKeyListener(new KeyAdapter() {
			public void keyReleased(KeyEvent e) {
				int keyCode = e.getKeyCode();
				if(keyCode==KeyEvent.VK_V){
					String value=cutoff_edge.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							double val = Double.parseDouble(value);
							if(val>=0 && val<=1.00){
								cutoff_edge = new JFormattedTextField(decimalformatter);
								cutoff_edge.setValue(val);
								metric.grabFocus();
								cutoff_edge.setFocusTraversalKeysEnabled(true);
								valid=false;
							}else{
								JOptionPane.showMessageDialog(null,"Cut-off value must be >= 0 and <=1 ","Warning",JOptionPane.ERROR_MESSAGE);
								cutoff_edge.grabFocus();
							}
						}else{
							cutoff_edge.setText(edgeCutoff);
						}
					}
				}
			}
			
        	public void keyPressed(KeyEvent e) {
        		int keyCode = e.getKeyCode();
        		if(keyCode==KeyEvent.VK_ENTER||keyCode==KeyEvent.VK_TAB){
        			String value=cutoff_edge.getText();
					if(!value.equals("")){
						boolean valid = isNumberFormat(value,true);
						if(valid){
							double val = Double.parseDouble(value);
							if(val>=0 && val<=1.00){
								cutoff_edge = new JFormattedTextField(decimalformatter);
								cutoff_edge.setValue(val);
								metric.grabFocus();
								cutoff_edge.setFocusTraversalKeysEnabled(true);
								valid=false;
							}else{
								JOptionPane.showMessageDialog(null,"Cut-off value must be >= 0 and <=1 ","Warning",JOptionPane.ERROR_MESSAGE);
								cutoff_edge.grabFocus();
							}
						}else{
							cutoff_edge.setText(edgeCutoff);
						}
					}
        		}
			}
		});
		
		metricadjust.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				String adjState=Integer.toString(metricadjust.getValue());
				adjustState.setText("Jaccard ("+(100-metricadjust.getValue())+"%) + Overlap ("+adjState+"%)");
			}
		});

		
		buildButton.addActionListener(new ActionListener(){
	        public void actionPerformed(ActionEvent e){
	        	mainThread = new Thread() {
	        		public void run() {
	    	        	String APW_fin = fin.getText();
	    	        	String GMT_fin = gmt.getText();
	    	        	String gsetmin = gset_min.getText();
	    	        	String gsetmax = gset_max.getText();
	    	        	String cutoffapw = cutoff_apw.getText();
	    	        	String sigapw = sig_apw.getText();
	    	        	String merge=mergemethod.getSelectedItem().toString();
	    	        	String correct=correctmethod.getSelectedItem().toString();
	    	        	String returnAPW=returnall.getSelectedItem().toString();
	    	        	String fcyto1=currentPath+File.separator+"terms.txt"; //defualt
	    	        	String fcyto2=currentPath+File.separator+"abridged.gmt"; //defualt
	    	        	String fcyto3=currentPath+File.separator+"groups.txt";//defualt
	    	        			
	    	        	if(fcytotik1.isSelected()){
	    	        		fcyto1=fcytoedit1.getText();
	    	        	}
	    	        	if(fcytotik2.isSelected()){
	    	        		fcyto2=fcytoedit2.getText();
	    	        	}
	    	        	if(fcytotik3.isSelected()){
	    	        		fcyto3=fcytoedit3.getText();
	    	        	}
	    	        	
	    	        	String mapName=networkname.getText();
	    	        	String edgeDataset= datasetedge.getSelectedItem().toString();
	    	        	String edgeCutoff = cutoff_edge.getText();
	    	        	String edgeMetric= metric.getSelectedItem().toString();
	    	        	if(edgeMetric.equals(metriclist[0])){
	    	        		edgeMetric="COMBINED";
	    	        	}else if(edgeMetric.equals(metriclist[1])){
	    	        		edgeMetric="JACCARD";
	    	        	}else if(edgeMetric.equals(metriclist[2])){
	    	        		edgeMetric="OVERLAP";
	    	        	}
	    	        	Double combinedConstant=Double.parseDouble(Integer.toString(metricadjust.getValue()))/100;
	    	        	String edgeStrategy = datasetedge.getSelectedItem().toString();
	    	        	if(edgeStrategy.equals(datasetedgelist[0])){
	    	        		edgeStrategy="AUTOMATIC";
	    	        	}else if(edgeStrategy.equals(datasetedgelist[1])){
	    	        		edgeStrategy="DISTINCT";
	    	        	}else if(edgeStrategy.equals(datasetedgelist[2])){
	    	        		edgeStrategy="COMPOUND";
	    	        	}
	    	        	
	    	        	String inputParams[]= {APW_fin,GMT_fin,gsetmin,gsetmax,cutoffapw,sigapw,merge,correct,returnAPW,fcyto1,fcyto2,fcyto3,mapName,edgeDataset,nodeFDR,nodePval,edgeCutoff,edgeCutoff};
	    	        	
	    	        	long t0 = System.currentTimeMillis();
	    	    		if(checkValidInputs(inputParams)){
	    		        	String cytofiles = fcyto1+","+fcyto2+","+fcyto3;
	    		        	genesetfilter = new ArrayList<Integer>();
	    	    			genesetfilter.add(Integer.parseInt(gsetmin));
	    	    			genesetfilter.add(Integer.parseInt(gsetmax));
	    	    			APW.mainmodule(inputParams[0], inputParams[1], hgmt_names, hgmt_genelist, Double.parseDouble(inputParams[2]), Double.parseDouble(inputParams[3]), inputParams[6], inputParams[7],inputParams[8], cytofiles,Double.parseDouble(inputParams[4]),Double.parseDouble(inputParams[5]),genesetfilter, progressPane, progressBar);
	    		        	if(isValidFile(fcyto1,"Terms.txt",true)){  //double check
	    		        		runEMCommand(mapName,fcyto1, fcyto2,fcyto3, nodePval,nodeFDR,edgeCutoff,edgeMetric,combinedConstant,edgeStrategy);
	    		        	}else{
	    		        		JOptionPane.showMessageDialog(null,"No significant terms were found.\r\nCytoscape files were not written","Warning",JOptionPane.ERROR_MESSAGE);
	    		        	}
	    	        	}
	    	    		long t1 = System.currentTimeMillis();
	    	    		float diff = (t1 - t0)/1000.0F;
	    	    		JOptionPane pane = new JOptionPane("End of process ("+diff+" seconds)",JOptionPane.INFORMATION_MESSAGE);
	    	    		final JDialog dialog = pane.createDialog(null, "Message");
	    	    		dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
	    	            dialog.setModal(false);
	    	    		dialog.setVisible(true);
	    	    		new Timer(TIME_VISIBLE, new ActionListener() {
	    	    			@Override
	    	    			public void actionPerformed(ActionEvent e) {
	    	    				dialog.setVisible(false);
	    	    			}
	    	    		}).start();
	    	    		dialog.dispose();
	        		}
	        		
	        	};
	        	
	        	progressbarThread = new Thread() {
	        		public void run() {
	        			JFrame progressFrame = new JFrame();
	    	        	progressFrame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
	    	        	progressFrame.setTitle("Please wait...");
	    	        	progressFrame.setSize(new Dimension(50,20));
	    	        	JPanel progressPane = new JPanel();
	    	    		progressPane.setLayout(new BorderLayout());
	    	    		progressImageLabel=new JLabel();
	    	    		progressImage = new ImageIcon(getClass().getResource(gif));
	    	    		progressImageLabel.setIcon(progressImage);
	    	    		progressPane.add(progressImageLabel,BorderLayout.CENTER);
	    	    		progressFrame.getContentPane().add(progressPane);
	    	    		progressFrame.setLocationRelativeTo(null);
	    	    		progressFrame.setVisible(true);
	    	    		progressFrame.pack();
	        		}
	        	};
	        	
	        	mainThread.start();
	        }
		});
		
		simpleOptionButton.addActionListener(new ActionListener() { 
        	public void actionPerformed(ActionEvent evt) {
        		frame.setVisible(false);
        		String APW_fin = fin.getText();
	        	String GMT_fin = gmt.getText();
	        	String gsetmin = gset_min.getText();
	        	String gsetmax = gset_max.getText();
	        	String cutoffapw = cutoff_apw.getText();
	        	String sigapw = sig_apw.getText();
        		initComponents_simple(frame,APW_fin,GMT_fin,gsetmin,gsetmax,cutoffapw,sigapw);
        	}
		});
		
        resetButton.addActionListener(new ActionListener() { 
        	public void actionPerformed(ActionEvent evt) {
        		fin.setText("");
        		gmt.setText("");
        		gset_min.setText("5");
        		gset_max.setText("1000");
        		cutoff_apw.setText("0.10");
        		sig_apw.setText("0.05");
        		mergemethod.setSelectedIndex(0); 
        		correctmethod.setSelectedIndex(0);
        		returnall.setSelectedIndex(0);
        		fcytotik1.setSelected(true);
        		fcytotik2.setSelected(true);
        		fcytotik3.setSelected(true);
        		fcytoedit1.setText("terms.txt");
        		fcytoedit2.setText("abridged.gmt");
        		fcytoedit3.setText("groups.txt");
        		networkname.setText("Networgk 1");
        		datasetedge.setSelectedIndex(0);
        		cutoff_edge.setValue("0.66");
        		metric.setSelectedIndex(0);
        		metricadjust.setValue(50);
            }
        });
       
        cancelButton.addActionListener(new ActionListener() {
           public void actionPerformed(ActionEvent evt) {
        	   int option = JOptionPane.showConfirmDialog((Component) null, "Exit?","Alert", JOptionPane.YES_NO_CANCEL_OPTION);
    		   if(option==JOptionPane.YES_OPTION){
    			   frame.dispose();
    		   }
    	   }
       });
        
       if(APW_fin.isEmpty()) {
       	 fin.grabFocus();
       }else {
       	 networkname.grabFocus();
       }
	}
      
	private JXTitledPanel create_APW_tab(){
		/*************************************************/
		/* Part1: ActivePathways Options
		/*************************************************/
		/* 0. Jpanel: ALL 
		 * 1. JTextField + Jbutton: Input file
		 * 2. JTextField + Jbutton: GMT file
		 * 3. JFormattedTextField: 2 geneset.filter Min[5], Max[1000]
		 * 4. JFormattedTextField: cutoff 
		 * 5. JFormattedTextField: Significant
		 * 6. JCombobox: merge.method (dropdown menu: "Brown", "Fisher", "logitp", "meanp","sump", "sumz", "sumlog")
		 * 7. JCombobox: correction.method (dropdown menu: "fdr","holm","hochberg", "hommel","bonferroni", "BH", "BY", "none")
		 * 8. JCombobox: return.all (yes/no)
		 * 9. JCheckbox + JTextField + Jbutton : 3 cytoscape.filenames [terms.txt, gmt.txt, groups.txt]
		*//**********************************************/
		
		JPanel detailsPanel = new JPanel(new GridBagLayout());
		JPanel jpfin = new JPanel();
		fin = new JTextField();
		jbfin= new JButton("...");
		fin.setEditable(true);
		fin.setPreferredSize(new Dimension(500,25));

		fin.setToolTipText("<html><p width=\"400\">" +fin_tip+"</p></html>");
		jbfin.setPreferredSize(new Dimension(32,25));
		jpfin.add(fin);
		jpfin.add(jbfin);
		detailsPanel.add(new JLabel("Input file:"), new GridBagConstraints(0, 0, 1, 1, 1.0,0.0, GridBagConstraints.NORTHEAST, GridBagConstraints.NONE, new Insets(7, 0, 0, 0), 0, 0));
		detailsPanel.add(jpfin, new GridBagConstraints(1, 0, 1, 1, 1.0,0.0, GridBagConstraints.NORTHEAST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
		
		JPanel jpgmt = new JPanel();
		gmt = new JTextField();
		jbgmt = new JButton("...");
		gmt.setEditable(true);
		gmt.setPreferredSize(new Dimension(500,25));
		gmt.setToolTipText("<html><p width=\"400\">" +gmt_tip+"</p></html>");
		jbgmt.setPreferredSize(new Dimension(32,25));
		jpgmt.add(gmt);
		jpgmt.add(jbgmt);
		detailsPanel.add(new JLabel("GMT file:"), new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0, GridBagConstraints.NORTHEAST, GridBagConstraints.NONE, new Insets(10, 0, 0, 0), 0, 0));
		detailsPanel.add(jpgmt, new GridBagConstraints(1,1,1,1,0.0,0.0, GridBagConstraints.NORTHEAST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
		gset_min = new JFormattedTextField();
		gset_min.setText("5");
		gset_min.setEditable(true);
		gset_min.setPreferredSize(new Dimension(90,25));
		gset_min.setToolTipText("<html><p width=\"300\">" +gsetmin_tip+"</p></html>");
		gset_max = new JFormattedTextField();
		gset_max.setText("1000");
		gset_max.setEditable(true);
		gset_max.setPreferredSize(new Dimension(90,25));
		gset_max.setToolTipText("<html><p width=\"300\">" +gsetmax_tip+"</p></html>");
		cutoff_apw = new JFormattedTextField();
		cutoff_apw.setText("0.1");
		cutoff_apw.setEditable(true);
		cutoff_apw.setPreferredSize(new Dimension(90,25));
		cutoff_apw.setToolTipText("<html><p width=\"400\">" +gsetcutoff_tip+"</p></html>");
		sig_apw = new JFormattedTextField();
		sig_apw.setText("0.05");
		sig_apw.setEditable(true);
		sig_apw.setPreferredSize(new Dimension(90,25));
		sig_apw.setToolTipText("<html><p width=\"400\">" +gsetsig_tip+"</p></html>");
	
        JPanel jpgset = new JPanel();
		jpgset.add(new JLabel("min"));
		jpgset.add(gset_min);
		jpgset.add(new JLabel("max"));
		jpgset.add(gset_max);
		jpgset.add(new JLabel("cutoff"));
		jpgset.add(cutoff_apw);
		jpgset.add(new JLabel("significant"));
		jpgset.add(sig_apw);
		JLabel gsetfilter = new JLabel("Geneset filter:");
		gsetfilter.setToolTipText("<html><p width=\"400\">" +gset_tip+"</p></html>");
        detailsPanel.add(gsetfilter, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0, GridBagConstraints.NORTHEAST, GridBagConstraints.NONE, new Insets(10, 0, 0, 0), 0, 0));
        detailsPanel.add(jpgset, new GridBagConstraints(1,2,1,1,0.0,0.0, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));

        JPanel jpmethod  = new JPanel();
		mergemethod=new JComboBox<String>(mergelist); 
		mergemethod.setToolTipText("<html><p width=\"200\">" +mergemethod_tip+"</p></html>");
		
		correctmethod =new JComboBox<String>(correctlist);
		correctmethod.setToolTipText("<html><p width=\"200\">" +correctmethod_tip+"</p></html>");
		
		
		jpmethod.add(mergemethod);
		
		JPanel jpreturn = new JPanel();
		returnall = new JComboBox<String>(yesno);
		returnall.setToolTipText("<html><p width=\"400\">" +runAll_tip+"</p></html>");
		jpreturn.add(new JLabel("Return all?"));
		jpreturn.add(returnall);

        detailsPanel.add(new JLabel("Merge method:"), new GridBagConstraints(0, 3, 1, 1, 1.0, 2.0, GridBagConstraints.NORTHEAST, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 0, 0));
        detailsPanel.add(mergemethod, new GridBagConstraints(1,3,1,1,1.0,2.0, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 2, 5), 0, 0));
        detailsPanel.add(new JLabel("Correction method:"), new GridBagConstraints(0, 4, 1, 1, 1.0, 2.0, GridBagConstraints.NORTHEAST, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 0, 0));
        detailsPanel.add(correctmethod, new GridBagConstraints(1,4,1,1,1.0,2.0, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 2, 5), 0, 0));
        detailsPanel.add(new JLabel("Return All:"), new GridBagConstraints(0, 5, 1, 1, 1.0, 2.0, GridBagConstraints.NORTHEAST, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 0, 0));
        detailsPanel.add(returnall, new GridBagConstraints(1,5,1,1,1.0,2.0, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5,2,5), 0, 0));
        
		JPanel jpcyto1=new JPanel();
		JPanel jpcyto2=new JPanel();
		JPanel jpcyto3=new JPanel();
		fcytotik1 = new JCheckBox();
		fcytotik2 = new JCheckBox();
		fcytotik3 = new JCheckBox();
		fcytotik1.setSelected(true);
		fcytotik2.setSelected(true);
		fcytotik3.setSelected(true);
		fcytoedit1 = new JTextField();
		fcytoedit2 = new JTextField();
		fcytoedit3 = new JTextField();
		
		fcytoedit1.setText("terms.txt");
		fcytoedit1.setEditable(true);
		fcytoedit1.setPreferredSize(new Dimension(380,25));
		fcytoedit1.setToolTipText("<html><p width=\"400\">" +cytofile1_tip+"</p></html>");
		
		fcytoedit2.setText("abridged.gmt");
		fcytoedit2.setEditable(true);
		fcytoedit2.setPreferredSize(new Dimension(380,25));
		fcytoedit2.setToolTipText("<html><p width=\"400\">" +cytofile2_tip+"</p></html>");
				
		fcytoedit3.setText("groups.txt");
		fcytoedit3.setEditable(true);
		fcytoedit3.setPreferredSize(new Dimension(380,25));
		fcytoedit3.setToolTipText("<html><p width=\"400\">" +cytofile3_tip+"</p></html>");
		
		fcytosave1 = new JButton("...");
		fcytosave1.setPreferredSize(new Dimension(10,25));
		fcytosave2 = new JButton("...");
		fcytosave2.setPreferredSize(new Dimension(10,25));
		fcytosave3 = new JButton("...");
		fcytosave3.setPreferredSize(new Dimension(10,25));

		jpcyto1.add(fcytotik1);
		jpcyto1.add(new JLabel("Terms"));
		jpcyto1.add(fcytoedit1);
		jpcyto1.add(fcytosave1);

		jpcyto2.add(fcytotik2);
		jpcyto2.add(new JLabel("GMT"));
		jpcyto2.add(fcytoedit2);
		jpcyto2.add(fcytosave2);
		
		jpcyto3.add(fcytotik3);
		jpcyto3.add(new JLabel("Grouping"));
		jpcyto3.add(fcytoedit3);
		jpcyto3.add(fcytosave3);
		
		JPanel	detailsPanel_cytofiles= new JPanel(new GridBagLayout());
		detailsPanel_cytofiles.add(fcytotik1, new GridBagConstraints(0,1,1,1,1.0,1.0, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(2, 0, 0, 0), 0, 0));
		detailsPanel_cytofiles.add(new JLabel("Terms"), new GridBagConstraints(1,1,1,1,1.0,0.0, GridBagConstraints.NORTHWEST, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 0, 0));
		detailsPanel_cytofiles.add(fcytoedit1, new GridBagConstraints(2,1,1,1,1.0,1.0, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(2, 0, 0, 2), 0, 0));
		detailsPanel_cytofiles.add(fcytosave1, new GridBagConstraints(3,1,1,1,1.0,1.0, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(2, 0, 2, 0), 0, 0));
		detailsPanel_cytofiles.add(fcytotik2, new GridBagConstraints(0,2,1,1,1.0,1.0, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
		detailsPanel_cytofiles.add(new JLabel("GMT"), new GridBagConstraints(1,2,1,1,1.0,0.0, GridBagConstraints.NORTHWEST, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 0, 0));
		detailsPanel_cytofiles.add(fcytoedit2, new GridBagConstraints(2,2,1,1,1.0,1.0, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 2), 0, 0));
		detailsPanel_cytofiles.add(fcytosave2, new GridBagConstraints(3,2,1,1,1.0,1.0, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 2, 0), 0, 0));
		detailsPanel_cytofiles.add(fcytotik3, new GridBagConstraints(0,3,1,1,1.0,1.0, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
		detailsPanel_cytofiles.add(new JLabel("Grouping"), new GridBagConstraints(1,3,1,1,1.0,0.0, GridBagConstraints.NORTHWEST, GridBagConstraints.NONE, new Insets(5, 0, 0, 0), 0, 0));
		detailsPanel_cytofiles.add(fcytoedit3, new GridBagConstraints(2,3,1,1,1.0,1.0, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 2), 0, 0));
		detailsPanel_cytofiles.add(fcytosave3, new GridBagConstraints(3,3,1,1,1.0,1.0, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 2, 0), 0, 0));
        detailsPanel.add(new JLabel("Cytoscape files:"), new GridBagConstraints(0, 6, 1, 1, 0.0, 0.0, GridBagConstraints.NORTHEAST, GridBagConstraints.NONE, new Insets(2, 0, 0, 0), 0, 0));
        detailsPanel.add(detailsPanel_cytofiles, new GridBagConstraints(1,6,1,1,1.0,1.0, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 4), 0, 0));
        
		JXTitledPanel detailPanel = new JXTitledPanel("Run ActivePathways");
        detailPanel.setContentContainer(detailsPanel);
        return detailPanel;
	}
	
	private JXTitledPanel create_APW_tab_simple(){
		/*************************************************/
		/* Part1: ActivePathways Options
		/*************************************************/
		/* 0. Jpanel: ALL 
		 * 1. JTextField + Jbutton: Input file
		 * 2. JTextField + Jbutton: GMT file
		 * 3. JFormattedTextField: 2 geneset.filter Min[5], Max[1000]
		 * 4. JFormattedTextField: cutoff 
		 * 5. JFormattedTextField: Significant
		 * 6. JCombobox: merge.method (dropdown menu: "Brown", "Fisher", "logitp", "meanp","sump", "sumz", "sumlog")
		 * 7. JCombobox: correction.method (dropdown menu: "fdr","holm","hochberg", "hommel","bonferroni", "BH", "BY", "none")
		 * 8. JCombobox: return.all (yes/no)
		 * 9. JCheckbox + JTextField + Jbutton : 3 cytoscape.filenames [terms.txt, gmt.txt, groups.txt]
		*//**********************************************/
		
		JPanel detailsPanel = new JPanel(new GridBagLayout());
		JPanel jpfin = new JPanel();
		fin = new JTextField();
		jbfin= new JButton("...");
		fin.setEditable(true);
		fin.setPreferredSize(new Dimension(500,25));

		fin.setToolTipText("<html><p width=\"400\">" +fin_tip+"</p></html>");
		jbfin.setPreferredSize(new Dimension(32,25));
		jpfin.add(fin);
		jpfin.add(jbfin);
		detailsPanel.add(new JLabel("Input file:"), new GridBagConstraints(0, 0, 1, 1, 1.0,0.0, GridBagConstraints.NORTHEAST, GridBagConstraints.NONE, new Insets(7, 0, 0, 0), 0, 0));
		detailsPanel.add(jpfin, new GridBagConstraints(1, 0, 1, 1, 1.0,0.0, GridBagConstraints.NORTHEAST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
		
		JPanel jpgmt = new JPanel();
		gmt = new JTextField();
		jbgmt = new JButton("...");
		gmt.setEditable(true);
		gmt.setPreferredSize(new Dimension(500,25));
		gmt.setToolTipText("<html><p width=\"400\">" +gmt_tip+"</p></html>");
		jbgmt.setPreferredSize(new Dimension(32,25));
		jpgmt.add(gmt);
		jpgmt.add(jbgmt);
		detailsPanel.add(new JLabel("GMT file:"), new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0, GridBagConstraints.NORTHEAST, GridBagConstraints.NONE, new Insets(10, 0, 0, 0), 0, 0));
		detailsPanel.add(jpgmt, new GridBagConstraints(1,1,1,1,0.0,0.0, GridBagConstraints.NORTHEAST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
		gset_min = new JFormattedTextField();
		gset_min.setText("5");
		gset_min.setEditable(true);
		gset_min.setPreferredSize(new Dimension(90,25));
		gset_min.setToolTipText("<html><p width=\"300\">" +gsetmin_tip+"</p></html>");
		gset_max = new JFormattedTextField();
		gset_max.setText("1000");
		gset_max.setEditable(true);
		gset_max.setPreferredSize(new Dimension(90,25));
		gset_max.setToolTipText("<html><p width=\"300\">" +gsetmax_tip+"</p></html>");
		cutoff_apw = new JFormattedTextField();
		cutoff_apw.setText("0.1");
		cutoff_apw.setEditable(true);
		cutoff_apw.setPreferredSize(new Dimension(90,25));
		cutoff_apw.setToolTipText("<html><p width=\"400\">" +gsetcutoff_tip+"</p></html>");
		sig_apw = new JFormattedTextField();
		sig_apw.setText("0.05");
		sig_apw.setEditable(true);
		sig_apw.setPreferredSize(new Dimension(90,25));
		sig_apw.setToolTipText("<html><p width=\"400\">" +gsetsig_tip+"</p></html>");
	
        JPanel jpgset = new JPanel();
		jpgset.add(new JLabel("min"));
		jpgset.add(gset_min);
		jpgset.add(new JLabel("max"));
		jpgset.add(gset_max);
		jpgset.add(new JLabel("cutoff"));
		jpgset.add(cutoff_apw);
		jpgset.add(new JLabel("significant"));
		jpgset.add(sig_apw);
		JLabel gsetfilter = new JLabel("Geneset filter:");
		gsetfilter.setToolTipText("<html><p width=\"400\">" +gset_tip+"</p></html>");
        detailsPanel.add(gsetfilter, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0, GridBagConstraints.NORTHEAST, GridBagConstraints.NONE, new Insets(10, 0, 0, 0), 0, 0));
        detailsPanel.add(jpgset, new GridBagConstraints(1,2,1,1,0.0,0.0, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
		JXTitledPanel detailPanel = new JXTitledPanel("Run ActivePathways");
        detailPanel.setContentContainer(detailsPanel);
        return detailPanel;
	}
	
	private JXTitledPanel create_ENR_tab(){
		/*************************************************/
		/* Part2: EnrichmentMap Options (Analysis type: Generic/gProfiler analysis)
		/*************************************************/
		/* 2.1 Number of Nodes (gene-set filtering)
		 * **************
		 * 0. JPanel: ALL
		 * 1. JTextField: Network Name
		 * 3. JFormattedTextField: FDR q-value cutoff [0.1]
		 * 4. JFormattedTextField: p-value cutoff [1.0]
		 * 
		 *************** 
		 * 2.2 JPanel: Number of Edges (gene-set similarity filtering)
		 ***************
		 * 0. JPanel: ALL 
		 * 1. JCombobox: Data Set Edges [Automatic, Separate edge for each data set (denser, Combine edges across data sets (sparser)]
		 * 2. JFormat: Cutoff (allow both decimal and scientific notation)
		 * 3. JCombobox: Metic[Jaccard, Overlap, Jaccard + Overlap Combined]
		 * 4. JSlider: Jaccard (50%) + Overlap (50%)
		 *//************************************************/
		
		JPanel detailsPanel = new JPanel(new GridBagLayout());
		JPanel jpname = new JPanel();
		networkname = new JTextField("Network 1");
		networkname.setEditable(true);
		networkname.setPreferredSize(new Dimension(550,25));
		networkname.setToolTipText("<html><p width=\"400\">" +networkName_tip+"</p></html>");
		jpname.add(new JLabel("Network Name:"));
		jpname.add(networkname);
		
	JXTitledPanel EdgeFilter = new JXTitledPanel();
        EdgeFilter= create_ENR_EdgeFilter();
        detailsPanel.add(jpname, new GridBagConstraints(1, 0, 1, 1, 1.0, 0.05, GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(5, 0, 5, 1), 0, 0));
        detailsPanel.add(EdgeFilter, new GridBagConstraints(1, 1, 1, 1, 1.0, 1.0, GridBagConstraints.NORTHEAST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0,0), 0, 0));
        	JXTitledPanel detailPanel = new JXTitledPanel("Create Enrichment Map");
		detailPanel.setContentContainer(detailsPanel);
		return detailPanel;
	}
	
	private JXTitledPanel create_ENR_NodeFilter(){
		JPanel detailsPanel = new JPanel(new GridBagLayout());
		fdr_enr_node = new JFormattedTextField();
		pval_enr_node = new JFormattedTextField();
		
		try {
			fdr_enr_node.setValue("0.1");
			fdr_enr_node.setEditable(true);
			fdr_enr_node.setPreferredSize(new Dimension(60,25));
			fdr_enr_node.setToolTipText("<html><p width=\"400\">" +node_qval_tip+"</p></html>");
			
			pval_enr_node.setValue("1.0");
			pval_enr_node.setEditable(true);
			pval_enr_node.setPreferredSize(new Dimension(60,25));
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		detailsPanel.add(new JLabel("FDR q-value cutoff:"), new GridBagConstraints(1, 0, 1, 1, 1.0, 0.02, GridBagConstraints.NORTHEAST, GridBagConstraints.NONE, new Insets(10, 0, 0, 5), 0, 0));
		detailsPanel.add(fdr_enr_node, new GridBagConstraints(2, 0, 1, 1, 1.0, 0.02, GridBagConstraints.NORTHEAST, GridBagConstraints.HORIZONTAL, new Insets(7, 0, 0, 5), 0, 0));
		
		JXTitledPanel detailPanel = new JXTitledPanel("Number of Nodes (geneset filtering)");
	    detailPanel.setContentContainer(detailsPanel);
	    return detailPanel;
	}
		
	private JXTitledPanel create_ENR_EdgeFilter(){
		JPanel detailsPanel = new JPanel(new GridBagLayout());
		
		datasetedge = new JComboBox<String>(datasetedgelist);
		datasetedge.setToolTipText("<html><p width=\"400\">" +edgestrategy_tip+"</p></html>");
		cutoff_edge= new JFormattedTextField();
		JLabel metriclabel= new JLabel("Metric:"); 
		metriclabel.setToolTipText("<html><p width=\"300\">" +edgemetriclabel_tip+"</p></html>");
		try {
			cutoff_edge = new JFormattedTextField();
			cutoff_edge.setText("0.66");
			cutoff_edge.setEditable(true);
			cutoff_edge.setPreferredSize(new Dimension(60,25));
			cutoff_edge.setToolTipText("<html><p width=\"400\">" +edgecutoff_tip+"</p></html>");
		} catch (Exception e) {
			e.printStackTrace();
		}
		metric = new JComboBox<String>(metriclist);
		metricadjust =new JSlider(JSlider.HORIZONTAL,0,100,50);
		metric.setToolTipText("<html><p width=\"400\">" +edgemetricoption_tip+"</p></html>");
		detailsPanel.add(new JLabel("Data Set Edges:"), new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0, GridBagConstraints.NORTHEAST, GridBagConstraints.NONE, new Insets(7, 5, 0, 5), 0, 0));
		detailsPanel.add(datasetedge, new GridBagConstraints(2, 0, 1, 1, 1.0, 1.0, GridBagConstraints.NORTH, GridBagConstraints.HORIZONTAL, new Insets(5, 0, 5, 5), 0, 0));
		
		detailsPanel.add(new JLabel("Cutoff:"), new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0, GridBagConstraints.NORTHEAST, GridBagConstraints.NONE, new Insets(2, 0, 0, 5), 0, 0));
		detailsPanel.add(cutoff_edge, new GridBagConstraints(2, 1, 1, 1, 1.0, 1.0, GridBagConstraints.NORTH, GridBagConstraints.HORIZONTAL, new Insets(0, 0,0, 5), 0, 0));
		detailsPanel.add(metriclabel, new GridBagConstraints(1, 3, 1, 1, 0.0, 0.0, GridBagConstraints.NORTHEAST, GridBagConstraints.NONE, new Insets(7, 0, 0, 5), 0, 0));
		detailsPanel.add(metric, new GridBagConstraints(2, 3, 1, 1, 1.0, 1.0, GridBagConstraints.NORTH, GridBagConstraints.HORIZONTAL, new Insets(5, 0,0, 5), 0, 0));
		detailsPanel.add(adjustState, new GridBagConstraints(2, 4, 1, 1, 0.0, 0.0, GridBagConstraints.NORTH, GridBagConstraints.NONE, new Insets(5, 0,0, 5), 0, 0));
		detailsPanel.add(metricadjust, new GridBagConstraints(2, 5, 1, 1, 1.0, 0.0, GridBagConstraints.NORTH, GridBagConstraints.HORIZONTAL, new Insets(5, 0, 5, 5), 0, 0));
		JXTitledPanel detailPanel = new JXTitledPanel("Number of Edges (geneset similarity filtering)");
	    detailPanel.setContentContainer(detailsPanel);
	    return detailPanel;
	}
	
	private boolean isValidFile(String filename, String prefix, boolean popupOption){
		boolean result=false;
		File path = new File(filename);
		if(path.isFile() || path.isDirectory()){
			if(path.exists() && path.getParent()!=null){
				result=true;
			}
		}
		
		if(popupOption){
			if(result){
				if(!path.isFile()&&path.isDirectory()){
					JOptionPane.showMessageDialog(null,"Please select Valid File","Warning",JOptionPane.ERROR_MESSAGE);	
				}else if(!path.isFile()&&!path.isDirectory()){
					JOptionPane.showMessageDialog(null,"Path not found","Warning",JOptionPane.ERROR_MESSAGE);
				}
			}else{
				JOptionPane.showMessageDialog(null,prefix+"File not found","Warning",JOptionPane.ERROR_MESSAGE);
			}
		}
		
		return result;
	}
	
	private boolean isNumberFormat(String keyIn, boolean popupOption){
		boolean	result = true;
		try {
	    	new BigDecimal(keyIn);
	    } catch (NumberFormatException e) {
	    	if(popupOption){
	    		JOptionPane.showMessageDialog(null,"Please put a Valid Number","Warning",JOptionPane.ERROR_MESSAGE);	
	    	}	    	
	        result=false;
	    }
		return result;
	}
	
	private boolean checkValidInputs(String inputParams []){
    	boolean resulttrue = false;
    	boolean resultfalse = false;
    	summary = new StringBuffer();
    	summary.append("Arguments are not valid!\r\n");
    	Integer fnames_ind[] = new Integer[]{0,1};
		String fnames[] = new String[]{"Input ","GMT ","Terms ","abridged GMT ","Grouping "};
		Integer numberFormat_ind[]= new Integer[]{2,3,4,5,13};
		for (int i = 0; i < fnames_ind.length; i++) {
			String prefix=fnames[i];
			if (isValidFile(inputParams[fnames_ind[i]],prefix,false)){
				resulttrue=true;
			}else{
				summary.append("-"+prefix+" file"+"\r\n");
				resultfalse=true;
			}
		}
		if(resulttrue){
			for (int i = 0; i < numberFormat_ind.length; i++){
				if(isNumberFormat(inputParams[numberFormat_ind[i]], false)){
					if(numberFormat_ind[i]==2 || numberFormat_ind[i]==3) {
						int val=Integer.parseInt(inputParams[numberFormat_ind[i]]);
						if(val>0) {
							resulttrue=true;
						}else{
							if(numberFormat_ind[i]==2) {
								summary.append("-geneset min: "+inputParams[numberFormat_ind[i]]+"\r\n");	
							}else {
								summary.append("-geneset max: "+inputParams[numberFormat_ind[i]]+"\r\n");
							}
							resultfalse=true;
						}
					}
					if(numberFormat_ind[i]==4||numberFormat_ind[i]==5||numberFormat_ind[i]==13) {
						double val=Double.parseDouble(inputParams[numberFormat_ind[i]]);
						if(val>=0 && val<=1.0) {
							resulttrue=true;
						}else {
							if(numberFormat_ind[i]==4) {
								summary.append("-cutoff: "+inputParams[numberFormat_ind[i]]+"\r\n");	
							}else if(numberFormat_ind[i]==5) {
								summary.append("-significant: "+inputParams[numberFormat_ind[i]]+"\r\n");
							}else {
								summary.append("-edge's cut-off: "+inputParams[numberFormat_ind[i]]+"\r\n");
							}
							resultfalse=true;
						}
					}
				}else{
					if(numberFormat_ind[i]==2) {
						summary.append("-geneset min: "+inputParams[numberFormat_ind[i]]+"\r\n");	
					}else if(numberFormat_ind[i]==3){
						summary.append("-geneset max: "+inputParams[numberFormat_ind[i]]+"\r\n");
					}if(numberFormat_ind[i]==4) {
						summary.append("-cutoff: "+inputParams[numberFormat_ind[i]]+"\r\n");	
					}else if(numberFormat_ind[i]==5) {
						summary.append("-significant: "+inputParams[numberFormat_ind[i]]+"\r\n");
					}else {
						summary.append("-edge's cut-off: "+inputParams[numberFormat_ind[i]]+"\r\n");
					}
					resultfalse=true;
				} 
			}	
		}
		if(resultfalse) {
			JOptionPane.showMessageDialog(null,summary.toString(),"Warning",JOptionPane.ERROR_MESSAGE);
			return false;
		}else {
			return true;
		}
	}
	
	public void runEMCommand(String mapName, String fcyto, String fgmt, String finstruct, String pval,String qval,String similaritycutoff,String coeffecients,double combinedConstant, String edgeStrategy){
		Map<String, Object>build_args = new HashMap<>();
		Map<String, Object>chart_args = new HashMap<>();
		Map<String, Object>visualStyle_args = new HashMap<>();
		Map<String, Object>visualApply = new HashMap<>();
		
		build_args.put("networkName", mapName);
		build_args.put("analysisType", "generic");
		build_args.put("enrichmentsData1",fcyto);
		build_args.put("gmtFile",fgmt);
		build_args.put("pvalue",Double.parseDouble(pval));
		build_args.put("qvalue",Double.parseDouble(qval));
		build_args.put("similaritycutoff",Double.parseDouble(similaritycutoff));
		build_args.put("coeffecients",coeffecients.toUpperCase());
		build_args.put("edgeStrategy",edgeStrategy.toUpperCase());
		build_args.put("combinedConstant",combinedConstant);
		executeCommand("enrichmentmap", "build", build_args, null);
		
		chart_args.put("file", finstruct);
		chart_args.put("keyColumnIndex",1);
		chart_args.put("startLoadRow",1);
		chart_args.put("firstRowAsColumnNames",true);
		chart_args.put("DataTypeTargetForNetworkCollection","Node Table Columns");
		chart_args.put("KeyColumnForMapping","shared name");
		chart_args.put("WhereImportTable","To selected networks only");
		chart_args.put("targetNetworkList",mapName);
		executeCommand("table","import file", chart_args, null);

		//Set default visual style
		String visURL = "";
		try {
			visURL = readStringFromURL(visualStyle);
			visualStyle_args.put("file",visURL);
			executeCommand("vizmap", "load file", visualStyle_args, null);
		} catch (IOException e) {
			e.printStackTrace();
		}

		visualApply.put("styles","APW_Visual_Style");
		executeCommand("vizmap","apply",visualApply,null);
	}

	public void executeCommand(String namespace, String command, Map<String, Object> args, TaskObserver observer) {
		if (commandTaskFactory == null)
			commandTaskFactory = getService(CommandExecutorTaskFactory.class);
		if (taskManager == null)
			taskManager = getService(SynchronousTaskManager.class);
		TaskIterator ti = commandTaskFactory.createTaskIterator(namespace, command, args, observer);
		taskManager.execute(ti);
	}
	
	public <S> S getService(Class<S> serviceClass) {
		return serviceRegistrar.getService(serviceClass);
	}

	public <S> S getService(Class<S> serviceClass, String filter) {
		return serviceRegistrar.getService(serviceClass, filter);
	}
	
	public void allFinished(FinishStatus finishStatus) {
	}

	public void taskFinished(ObservableTask task) {
		String models = task.getResults(String.class);
		int offset = models.indexOf(' ');
		if (offset >= 0) {
			String model = models.substring(1, offset);
			modelName = new String(models.substring(offset+1, models.length()-1));
			try {
				modelNumber = Integer.parseInt(model);
 			} catch (Exception e) {}
		} else {
			modelNumber = -2;
		}
	}
	
	/***************************************
	 *	create progress bar 
	/***************************************/
	class Task extends SwingWorker<Void, Void> {
	/* Main task. Executed in background thread.*/
		@Override
		public Void doInBackground() {
			Random random = new Random();
			int progress = 0;
			//Initialize progress property.
			setProgress(0);
			while (progress < 100 && checkActiveProgressBar) {
				//Sleep for up to one second.
				try {
					Thread.sleep(random.nextInt(10000));
				} catch (InterruptedException ignore) {}
	                //Make random progress.
					progress += random.nextInt(1);
	                setProgress(Math.min(progress, 100));
	            }
	            return null;
		}
		
		/*Executed in event dispatching thread */
		@Override
		public void done() {
			Toolkit.getDefaultToolkit().beep();
			//buildButton.setEnabled(true);
			buildButton.setCursor(null); //turn off the wait cursor
		}
	}
	
	static class SelectAll extends TextAction{
        public SelectAll(){
            super("Select All");
            putValue(Action.ACCELERATOR_KEY, KeyStroke.getKeyStroke("control S"));
        }
        public void actionPerformed(ActionEvent e){
            JTextComponent component = getFocusedComponent();
            component.selectAll();
            component.requestFocusInWindow();
        }
    }
	
	private void repaintBar(JPanel progressPane, JProgressBar progressBar,int percentage) {
		progressPane.revalidate();
		progressPane.repaint();
		progressBar.setValue(percentage);
	}

	
	public String get_userAgrs() {
		String APW_fin = fin.getText();
    	String GMT_fin = gmt.getText();
    	String gsetmin = gset_min.getText();
    	String gsetmax = gset_max.getText();
    	String cutoffapw = cutoff_apw.getText();
    	String sigapw = sig_apw.getText();
    	String merge="Brown";
    	String correct="fdr";
    	String returnAPW="no";
    	Path path = Paths.get(APW_fin).toAbsolutePath();
		String finpath = path.getParent().toString();
    	String fcyto1=finpath+File.separator+"terms.txt";
    	String fcyto2=finpath+File.separator+"abridged.gmt";
    	String fcyto3=finpath+File.separator+"groups.txt";
      	String mapName=APW_fin.split(File.separator)[APW_fin.split(File.separator).length-1].split("\\.")[0];
    	String edgeCutoff = "0.66";
    	String edgeMetric= "COMBINED";
    	Double combinedConstant=0.50;
    	String edgeStrategy = "AUTOMATIC";
    	String inputParams= APW_fin+","+GMT_fin+","+gsetmin+","+gsetmax+","+cutoffapw+","+sigapw+","+merge+","+correct+","+returnAPW+","+fcyto1+","+fcyto2+","+fcyto3+","+mapName+","+edgeCutoff+","+edgeMetric+","+combinedConstant+","+edgeStrategy;
    	return inputParams;
	}
	
	public static String readStringFromURL(String requestURL) throws IOException{
	    try (Scanner scanner = new Scanner(new URL(requestURL).openStream(),StandardCharsets.UTF_8.toString())){
	        scanner.useDelimiter("\\A");
	        return scanner.hasNext() ? scanner.next() : "";
	    }
	}
	

}
