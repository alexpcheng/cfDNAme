package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map.Entry;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import shared.Parser;
import shared.PreParser;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.MultiCros2;
import stream.Read;
import structures.ListNum;


/**
 * This class is designed to handle very large numbers of output files.
 * @author Brian Bushnell
 * @date May 1, 2019
 *
 */
public class DemuxByName2 {

	public static void main(String[] args){
		
		final int oldCap=Shared.numBuffers(), oldZipThreads=ReadWrite.MAX_ZIP_THREADS, oldZl=ReadWrite.ZIPLEVEL;
		final boolean oldPigz=ReadWrite.USE_PIGZ, oldUnpigz=ReadWrite.USE_UNPIGZ;
		
		Timer t=new Timer();
		DemuxByName2 x=new DemuxByName2(args);
		x.process(t);
		
		Shared.setBuffers(oldCap);
		ReadWrite.ZIPLEVEL=oldZl;
		ReadWrite.USE_PIGZ=oldPigz;
		ReadWrite.USE_UNPIGZ=oldUnpigz;
		ReadWrite.MAX_ZIP_THREADS=oldZipThreads;
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public DemuxByName2(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		boolean setInterleaved=false; //Whether it was explicitly set.
		
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=false;
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		ReadWrite.ZIPLEVEL=2;
		
		
		Parser parser=new Parser();
		parser.overwrite=overwrite;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
//				align2.FastaReadInputStream2.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("names") || a.equals("name") || a.equals("affixes")){
				if(b!=null){
					String[] x=b.split(",");
					for(String s : x){
						names.put(s, s);
					}
				}
			}else if(a.equalsIgnoreCase("length") || a.equalsIgnoreCase("len") || a.equalsIgnoreCase("affixlength") || a.equalsIgnoreCase("affixlen")){
				fixedAffixLength=Integer.parseInt(b);
			}else if(a.equalsIgnoreCase("readsperbuffer") || a.equalsIgnoreCase("rpb")){
				readsPerBuffer=Tools.parseIntKMG(b);
			}else if(a.equalsIgnoreCase("bytesPerBuffer") || a.equalsIgnoreCase("bpb")){
				bytesPerBuffer=Tools.parseIntKMG(b);
			}else if(a.equals("prefixmode") || a.equals("prefix") || a.equals("pm")){
				prefixMode=Tools.parseBoolean(b);
			}else if(a.equals("suffixmode") || a.equals("suffix") || a.equals("sm")){
				prefixMode=!Tools.parseBoolean(b);
			}else if(a.equals("column")){
				column=Integer.parseInt(b);
				assert(column>0 || column==-1) : "Column is 1-based; must be 1+ or else -1 to disable.";
				column--;
			}else if(a.equals("hdist") || a.equals("hamming") || a.equals("hammingdistance")){
				hdist=Integer.parseInt(b);
			}else if(a.equals("substringmode") || a.equals("substring")){
				substringMode=Tools.parseBoolean(b);
			}else if(a.equals("outu") || a.equals("outu1")){
				outu1=b;
			}else if(a.equals("outu2")){
				outu2=b;
			}else if(a.equals("pattern")){
				parser.out1=b;
			}else if(a.equals("perheader") || a.equals("persequence") || a.equals("pername")){
				perheader=Tools.parseBoolean(b);
			}else if(a.equals("delimiter")){
				if(b==null){delimiter=null;}
				
				//Convenience characters
				else if(b.equalsIgnoreCase("space")){
					delimiter=" ";
				}else if(b.equalsIgnoreCase("tab")){
					delimiter="\t";
				}else if(b.equalsIgnoreCase("whitespace")){
					delimiter="\\s+";
				}else if(b.equalsIgnoreCase("pound")){
					delimiter="#";
				}else if(b.equalsIgnoreCase("greaterthan")){
					delimiter=">";
				}else if(b.equalsIgnoreCase("lessthan")){
					delimiter="<";
				}else if(b.equalsIgnoreCase("equals")){
					delimiter="=";
				}else if(b.equalsIgnoreCase("colon")){
					delimiter=":";
				}else if(b.equalsIgnoreCase("semicolon")){
					delimiter=";";
				}else if(b.equalsIgnoreCase("bang")){
					delimiter="!";
				}else if(b.equalsIgnoreCase("and") || b.equalsIgnoreCase("ampersand")){
					delimiter="&";
				}else if(b.equalsIgnoreCase("quote") || b.equalsIgnoreCase("doublequote")){
					delimiter="\"";
				}else if(b.equalsIgnoreCase("singlequote") || b.equalsIgnoreCase("apostrophe")){
					delimiter="'";
				}
				
				//Java meta characters
				else if(b.equalsIgnoreCase("backslash")){
					delimiter="\\\\";
				}else if(b.equalsIgnoreCase("hat") || b.equalsIgnoreCase("caret")){
					delimiter="\\^";
				}else if(b.equalsIgnoreCase("dollar")){
					delimiter="\\$";
				}else if(b.equalsIgnoreCase("dot")){
					delimiter="\\.";
				}else if(b.equalsIgnoreCase("pipe") || b.equalsIgnoreCase("or")){
					delimiter="\\|";
				}else if(b.equalsIgnoreCase("questionmark")){
					delimiter="\\?";
				}else if(b.equalsIgnoreCase("star") || b.equalsIgnoreCase("asterisk")){
					delimiter="\\*";
				}else if(b.equalsIgnoreCase("plus")){
					delimiter="\\+";
				}else if(b.equalsIgnoreCase("openparen")){
					delimiter="\\(";
				}else if(b.equalsIgnoreCase("closeparen")){
					delimiter="\\)";
				}else if(b.equalsIgnoreCase("opensquare")){
					delimiter="\\[";
				}else if(b.equalsIgnoreCase("opencurly")){
					delimiter="\\{";
				}
				
				else{
					delimiter=b;
				}
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		if(perheader){
			ReadWrite.USE_PIGZ=false;
		}else{
			String[] x=names.keySet().toArray(new String[names.size()]);
			names.clear();
			for(String s : x){
				File f=new File(s);
				if(f.exists() && f.isFile()){
					TextFile tf=new TextFile(s);
					String[] lines=tf.toStringLines();
					for(String s2 : lines){
						names.put(s2, s2);
					}
				}else{
					names.put(s, s);
				}
			}

//			for(String s : names){
//				assert(affixLen<0 || affixLen==s.length()) : "All names must have the same length.";
//				affixLen=s.length();
//			}
//			assert(affixLen>0) : "Must include at least one non-zero-length affix (name).";
			
			{
				BitSet bs=new BitSet();
				if(fixedAffixLength>0){
					bs.set(fixedAffixLength);
				}
				for(String s : names.keySet()){
					bs.set(s.length());
				}
				affixLengths=new int[bs.cardinality()];
				for(int i=0, bit=-1; i<affixLengths.length; i++){
					bit=bs.nextSetBit(bit+1);
					affixLengths[i]=bit;
				}
				Arrays.sort(affixLengths);
				Tools.reverseInPlace(affixLengths);
			}
			
			assert((affixLengths.length>0 && affixLengths[0]>0) || delimiter!=null || perheader) : 
				"Must include at least one non-zero-length affix (name), or a delimiter.";
			ReadWrite.MAX_ZIP_THREADS=Tools.max(1, Tools.min(ReadWrite.MAX_ZIP_THREADS, (Shared.threads()*2-1)/Tools.max(1, names.size())));
			if(names.size()>8){ReadWrite.USE_PIGZ=false;}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;
			qfin1=parser.qfin1;
			qfin2=parser.qfin2;

			out1=parser.out1;
			out2=parser.out2;
			qfout1=parser.qfout1;
			qfout2=parser.qfout2;
			
			extin=parser.extin;
			extout=parser.extout;
		}

		assert(out1==null || out1.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(out2==null || out2.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(qfout1==null || qfout1.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		assert(qfout2==null || qfout2.contains("%")) : "Output filename must contain '%' symbol, which will be replaced by affix.";
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(out1!=null && out2==null && out1.indexOf('#')>-1){
			out2=out1.replace("#", "2");
			out1=out1.replace("#", "1");
		}
		if(outu1!=null && outu2==null && outu1.indexOf('#')>-1){
			outu2=outu1.replace("#", "2");
			outu1=outu1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(out1==null && out2!=null){throw new RuntimeException("Error - cannot define out2 without defining out1.");}
		
		if(!setInterleaved){
			assert(in1!=null && (out1!=null || out2==null)) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\nout2="+out2+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}else{ //There is one input stream.
				if(out2!=null){
					FASTQ.FORCE_INTERLEAVED=true;
					FASTQ.TEST_INTERLEAVED=false;
					outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
				}
			}
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		if(out2!=null && out2.equalsIgnoreCase("null")){out2=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1, out2, outu1, outu2)){
			outstream.println((out1==null)+", "+(out2==null)+", "+out1+", "+out2+", "+outu1+", "+outu2);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+out2+", "+outu1+", "+outu2+"\n");
		}

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
		
		if(ffin1!=null && out1!=null && ffin1.samOrBam()){
			String ext=ReadWrite.rawExtension(out1);
			useSharedHeader=FileFormat.isSamOrBam(ext);
		}
		
		if(perheader){
			fixedAffixLength=-1;
			substringMode=false;
			affixLengths=null;
			delimiter=null;
		}
		
		if(hdist>0 && names.size()>0){
			int mutants=mutate(names, hdist);
			System.err.println("Added "+mutants+" mutants.");
		}
		
		if(affixLengths.length>1) {fixedAffixLength=-1;}
		for(Entry<String, String> e : names.entrySet()){
			nameList.add(e.getKey());
		}
	}
	
	private static int mutate(HashMap<String, String> names, int hdist){
		ArrayList<String> list=new ArrayList<String>(names.size());
		list.addAll(names.values());
		int mutants=0;
		for(String name : list){
			mutants+=addMutants(name.getBytes(), name, hdist, names);
		}
		return mutants;
	}
	
	private static int addMutants(final byte[] key, final String value, final int hdist, final HashMap<String, String> names){
		assert(hdist>0);
		int added=0;
		for(int i=0; i<key.length; i++){
			final byte old=key[i];
			if(AminoAcid.isACGTN(old)){
				for(byte b : symbols){
					if(b!=old){
						key[i]=b;
						{
							String s=new String(key);
							if(!names.containsKey(s)){
								names.put(s, value);
								added++;
							}else{
								assert(value.equals(names.get(s))) : "Collision between "+value+" and "+names.get(s)+" for mutant "+s;
							}
						}
						if(hdist>1){
							added+=addMutants(key, value, hdist-1, names);
						}
					}
				}
				key[i]=old;
			}
		}
		return added;
	}
	
	//Returns the value
	private String getName(Read r){
		String name=getName0(r);
		if(name==null){return null;}
		return nameList.isEmpty() ? name : names.get(name);
	}
	
	//Returns the key
	private String getName0(Read r){
		final String id=r.id;
		final int idlen=id.length();
		
		if(nameList.size()>0){
			if(substringMode){
				for(String s : nameList){
					if(id.contains(s)){return s;}
				}
				return null;
			}else if(affixLengths.length>0){
				for(int affixLen : affixLengths){
					final String sub=idlen>=affixLen ? prefixMode ? id.substring(0, affixLen) : id.substring(idlen-affixLen) : id;
					if(names.containsKey(sub)) {return sub;}
				}
				return null;
			}
		}
		
		final String name;
		if(fixedAffixLength>0){
			name=(id.length()<=fixedAffixLength ? id : prefixMode ? id.substring(0, fixedAffixLength) : id.substring(idlen-fixedAffixLength));
		}else if(delimiter!=null){
			String[] split=id.split(delimiter);
			assert(split.length>1) : "Delimiter '"+delimiter+"' was not found in name '"+id+"'";
			if(column>-1){
				int col=Tools.min(column, split.length-1);
				name=split[col];
				if(col!=column && !warned){
					System.err.println("*** WARNING! ***\n"
							+ "Only "+(col+1)+" columns for record "+id+"\n"
							+ "Further warnings will be suppressed.\n");
					warned=true;
					assert(errorState=true); //Change error state to true if assertions are enabled.
				}
			}else{
				name=split[prefixMode ? 0 : split.length-1];
			}
		}else{
			assert(perheader);
			name=id;
		}
		
		return name;
	}
	
	void process(Timer t){
		
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, useSharedHeader, ffin1, ffin2, qfin1, qfin2);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		
		final MultiCros2 mcros;
		if(out1!=null){
			final int buff=1;
			
			mcros=new MultiCros2(out1, out2, overwrite, append, true, useSharedHeader, FileFormat.FASTQ, buff);
			mcros.readsPerBuffer=readsPerBuffer;
			mcros.bytesPerBuffer=bytesPerBuffer;
			
			if(paired && out2==null && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			assert(out2==null || (!out2.equalsIgnoreCase(in1) && !out2.equalsIgnoreCase(in2))) : "out1 and out2 have same name.";
		}else{
			mcros=null;
		}
		
		final ConcurrentReadOutputStream rosu;
		if(outu1!=null){
			final int buff=4;
			
			FileFormat ffout1=FileFormat.testOutput(outu1, FileFormat.FASTQ, extout, true, overwrite, append, false);
			FileFormat ffout2=(outu2==null ? null : FileFormat.testOutput(outu2, FileFormat.FASTQ, extout, true, overwrite, append, false));
			rosu=ConcurrentReadOutputStream.getStream(ffout1, ffout2, null, null, buff, null, true);
			rosu.start();
		}else{
			rosu=null;
		}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		long readsOut=0;
		long basesOut=0;
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			
			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				
				ArrayList<Read> unmatched=new ArrayList<Read>();
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					String name=getName(r1);
					if(name!=null){
						mcros.add(r1, name);
						
						readsOut+=r1.pairCount();
						basesOut+=r1.pairLength();
					}else{
						unmatched.add(r1);
					}
					
					readsProcessed+=r1.pairCount();
					basesProcessed+=r1.pairLength();
				}
				if(rosu!=null){rosu.add(unmatched, ln.id);}
				
				cris.returnList(ln);
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadStats.writeAll();
		
		errorState|=ReadWrite.closeStream(cris);
		
		if(mcros!=null){
			mcros.dumpAll();
			errorState|=mcros.errorState();
		}
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);
		
		outstream.println("Time:               "+t);
		outstream.println("Reads Processed:    "+readsProcessed+" \t"+String.format(Locale.ROOT, "%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+basesProcessed+" \t"+String.format(Locale.ROOT, "%.2fm bases/sec", bpnano*1000));
		outstream.println("Reads Out:          "+readsOut);
		outstream.println("Bases Out:          "+basesOut);
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String in2=null;
	
	private String qfin1=null;
	private String qfin2=null;

	private String out1=null;
	private String out2=null;

	private String outu1=null;
	private String outu2=null;

	private String qfout1=null;
	private String qfout2=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
//	private boolean exclude=true;

	private String delimiter=null;
	private boolean prefixMode=true;
	private boolean substringMode=false;
	private boolean perheader=false;
	private int column=-1;
	private boolean warned=false;
//	private int affixLen=-1;

	private int readsPerBuffer=1000;
	private int bytesPerBuffer=3000000;
	private int hdist=0;
	
	private int fixedAffixLength=-1;
	
	private int[] affixLengths;

	private HashMap<String, String> names=new HashMap<String, String>();
	private ArrayList<String> nameList=new ArrayList<String>();
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;
	
	private static final byte[] symbols={'A', 'C', 'G', 'T', 'N'};
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	private boolean useSharedHeader=false;
	
}
