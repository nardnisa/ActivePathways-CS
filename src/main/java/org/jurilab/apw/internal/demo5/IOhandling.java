package org.jurilab.apw.internal.demo5;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import javax.swing.JOptionPane;

public class IOhandling {
	public String loadFile(String filename)	{
    	StringBuffer content = new StringBuffer();
        try {
            FileInputStream 	fi 	= new FileInputStream(filename);
  				BufferedReader 	br 	= new BufferedReader(new InputStreamReader(fi));  			
  				String st = new String();
  				while ((st=br.readLine())!=null)	{
  					content.append(st.trim());
  					content.append("\r\n");
  				}
  				br.close();
  				fi.close();
        } 
        catch (IOException e) {
        	JOptionPane.showMessageDialog(null, "Loading file error!","Error",JOptionPane.ERROR_MESSAGE);
        }
        return content.toString();
    }
	
	public String check_valid_fout(String fout){
		File file = new File(fout);
		if(file.exists()){
			file.delete();
		}
		return fout;
	}
	
	public String get_fullPath(String filename){
		String fullpath= new String();
		try {
			File pathName = new File(filename);
			fullpath=pathName.getCanonicalPath();
		} catch (Exception e) {
			JOptionPane.showMessageDialog(null, "File not found","Error",JOptionPane.ERROR_MESSAGE);
		}
		return fullpath;
	}
	
	public static void writeFile(String filename, String content) throws Exception {
	    File file = new File(filename);
	    long fileLength = file.length();
	    RandomAccessFile raf = new RandomAccessFile(file, "rw");
	    raf.seek(fileLength);
	    raf.writeBytes(content);
	    raf.close();
	  }
}
