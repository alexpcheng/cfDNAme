package stream;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Shared;
import shared.Tools;
import structures.ListNum;

/**
 * Allows output of reads to multiple different output streams.
 * @author Brian Bushnell
 * @date Apr 12, 2015
 *
 */
public class MultiCros2 {
	
	public static void main(String[] args){
		String in=args[0];
		String pattern=args[1];
		ArrayList<String> names=new ArrayList<String>();
		for(int i=2; i<args.length; i++){
			names.add(args[i]);
		}
		final int buff=Tools.max(16, 2*Shared.threads());
		MultiCros2 mcros=new MultiCros2(pattern, null, false, false, false, false, FileFormat.FASTQ, buff);
		
		ConcurrentReadInputStream cris=ConcurrentReadInputStream.getReadInputStream(-1, true, false, in);
		cris.start();
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning

			for(Read r1 : reads){
				mcros.add(r1, r1.barcode(true));
			}
//			System.err.println("x");
			cris.returnList(ln);
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln);
		ReadWrite.closeStreams(cris);
		mcros.dumpAll();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	public MultiCros2(String pattern1_, String pattern2_,
			boolean overwrite_, boolean append_, boolean allowSubprocess_, boolean useSharedHeader_, int defaultFormat_, int maxSize_){
		assert(pattern1_!=null && pattern1_.indexOf('%')>=0);
		assert(pattern2_==null || pattern1_.indexOf('%')>=0);
		if(pattern2_==null && pattern1_.indexOf('#')>=0){
			pattern1=pattern1_.replaceFirst("#", "1");
			pattern2=pattern1_.replaceFirst("#", "2");
		}else{
			pattern1=pattern1_;
			pattern2=pattern2_;
		}
		
		overwrite=overwrite_;
		append=append_;
		allowSubprocess=allowSubprocess_;
		useSharedHeader=useSharedHeader_;
		
		defaultFormat=defaultFormat_;
		maxSize=maxSize_;
		
		bufferMap=new LinkedHashMap<String, Buffer>();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	public String fname(){return pattern1;}
	
	/** Return true if this stream has detected an error */
	public boolean errorState(){
		return errorState;
	}

	public boolean finishedSuccessfully(){
		return true; //TODO
	}
	
	public void add(Read r, String name){
		Buffer b=bufferMap.get(name);
		if(b==null){
			b=new Buffer(name);
			bufferMap.put(name, b);
		}
		b.add(r);
//		System.err.println("Added "+name);
	}
	
	public long dumpAll(){
		long dumped=0;
		for(Entry<String, Buffer> e : bufferMap.entrySet()){
			dumped+=e.getValue().dump();
		}
		return dumped;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------        Inner Classes         ----------------*/
	/*--------------------------------------------------------------*/
	
	private class Buffer {
		
		Buffer(String name_){
			name=name_;
			String s1=pattern1.replaceFirst("%", name);
			String s2=pattern2==null ? null : pattern2.replaceFirst("%", name);
			ff1=FileFormat.testOutput(s1, defaultFormat, null, allowSubprocess, false, true, false);
			ff2=FileFormat.testOutput(s2, defaultFormat, null, allowSubprocess, false, true, false);
			
			list=new ArrayList<Read>(readsPerBuffer);
		}
		
		void add(Read r){
			list.add(r);
			long size=r.countBytes()+r.countMateBytes();
			currentBytes+=size;
			bytesInFlight+=size;
			readsInFlight+=r.pairCount();
			if(list.size()>=readsPerBuffer || currentBytes>=bytesPerBuffer){dump();}
		}
		
		long dump(){
//			System.err.println("Dumping "+name);
			if(list.isEmpty()){return 0;}
			if(counter==0 && overwrite){
				delete(ff1);
				delete(ff2);
			}
//			assert(false) : counter+", "+overwrite;
			
			final long size0=list.size();
			
			ConcurrentReadOutputStream ros=ConcurrentReadOutputStream.getStream(ff1, ff2, maxSize, null, useSharedHeader && counter==0);
			ros.start();
			ros.add(list, 0);
			errorState=ReadWrite.closeStream(ros) | errorState;
//			System.err.println("Closed stream "+name);
			
			bytesInFlight-=currentBytes;
			readsInFlight-=size0;
			readsWritten+=size0;
			currentBytes=0;
			counter++;
			list=new ArrayList<Read>(readsPerBuffer);
			return size0;
		}
		
		private void delete(FileFormat ff){
			if(ff==null){return;}
			File f=new File(ff.name());
			if(f.exists()){f.delete();}
		}
		
		ArrayList<Read> list;
		final String name;
		final FileFormat ff1;
		final FileFormat ff2;
		long readsWritten=0;
		long currentBytes=0;
		long counter=0;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/
	
	public final String pattern1, pattern2;
	public final LinkedHashMap<String, Buffer> bufferMap;
	
	boolean errorState=false;
	boolean started=false;
	final boolean overwrite;
	final boolean append;
	final boolean allowSubprocess;
	final int defaultFormat;
	final int maxSize;
	final boolean useSharedHeader;

	public int readsPerBuffer=1000;
	public int bytesPerBuffer=3000000;
	long maxBytesInFlight=4000000000L;
	
	long readsInFlight=0;
	long bytesInFlight=0;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static boolean verbose=false;

}
