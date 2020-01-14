package coolname;

import java.util.*;
import java.util.regex.*;
import java.io.*;

public class GRN
{
  public int nGenes;
  public int nRootGenes; //number of unregulated genes.
  public String[] geneNames; //index to name,
  public Map<String,Integer> index; //and name to index
  /* NB : The index order must respect two things :
     1) The root genes must be first in the list, i.e. have the smaller indices. (This can be
        changed as long as the "isRegulated" function changes accordingly)
     2) The index of any gene must be greater than those of its regulators (i.e. indices follow
        a topological ordering of the network). In a bipartite network where all regulators are
        unregulated, this trivially follows from the first point.
   */
  public int[][] inhibitors, activators; // target-indexed regulator lists
  public StringBuilder buildLog=new StringBuilder(); // A log of what happened during the creation of this GRN.
  
  /** Read a bipartite GRN from a file
  */
  public GRN(File networkFile) throws IOException
  {
		index=new HashMap<String,Integer>();
		HashMap<String,int[]> activatorMap=new HashMap<>();
		HashMap<String,int[]> inhibitorMap=new HashMap<>();
		String line;
    int id=0;
		
    // fill Map "index" with (name,id) pairs of regulators then store ids in array "regulators"
    // fill Maps 'activatorMap' and 'inhibitorMap' with the regulator sets of target genes
    BufferedReader br=new BufferedReader(new FileReader(networkFile));
    Pattern linePattern=Pattern.compile("([^;]*);([^;]*);([^;]*)");
    Matcher m=null;
		br.readLine(); // discard header line
		while((line=br.readLine()) != null) 
     if((m=linePattern.matcher(line)).matches())
     {
      //extract info from line
      String targetName=m.group(1).trim();
			String[] activatorNames=m.group(2).trim().split("\\s+");
			String[] inhibitorNames=m.group(3).trim().split("\\s+");
      if(activatorNames[0].equals("")) activatorNames=new String[0]; // no word really means NO word, not 'a single, empty, word'
      if(inhibitorNames[0].equals("")) inhibitorNames=new String[0];

      if(activatorMap.containsKey(targetName))
      { buildLog.append("Ignoring the following line (redefinition of "+targetName+"\'s regulators):\n"+line+"\n");
        continue;
      } 
			activatorMap.put(targetName,new int[activatorNames.length]);
			inhibitorMap.put(targetName,new int[inhibitorNames.length]);
      //trim regulator names and store them if not already known
      for(int i=0;i<activatorNames.length; i++)
			{ if (!index.containsKey(activatorNames[i]=activatorNames[i].trim()))
          index.put(activatorNames[i],id++);
        activatorMap.get(targetName)[i]=index.get(activatorNames[i]);
      }
			for(int i=0;i<inhibitorNames.length; i++)
			{ if (!index.containsKey(inhibitorNames[i]=inhibitorNames[i].trim()))
          index.put(inhibitorNames[i],id++);
			  inhibitorMap.get(targetName)[i]=index.get(inhibitorNames[i]);
      } 
     }
		br.close();

    nRootGenes=id;

    // Now 'index' contains all regulators, add the targets to it.
    for(String name : activatorMap.keySet())
			if (index.containsKey(name))
        buildLog.append("Ignoring regulators of "+name+" : a regulator can not be regulated.\n");
      else
        index.put(name, id++);

    // Now that 'index' is full, make the "names" array from it
    nGenes=id;
    geneNames=new String[nGenes];  
    for(Map.Entry<String,Integer> entry : index.entrySet())
      geneNames[entry.getValue()]=entry.getKey();

    // And now, to the Network, i.e. the list of activators and inhibitors of each target gene
		inhibitors=new int[nGenes][];  
    activators=new int[nGenes][];
    for(int i=0;i<nGenes;i++)
     if(isRegulated(i))// (i.e. if i>=nRootGenes)
     {activators[i]=activatorMap.get(geneNames[i]);
      inhibitors[i]=inhibitorMap.get(geneNames[i]);
     }
     else
      activators[i]=inhibitors[i]=new int[0];
  }
  
  /** Read a bipartite GRN from three files describing it
  */
  public GRN(File targetsFile, File regulatorsFile, File networkFile) throws IOException
  {
		index=new HashMap<String,Integer>();
		
		String line;
    int id=0;
		
    // fill Map "index" with (name,id) pairs from regulators file, then store ids in array "regulators"
    BufferedReader br=new BufferedReader(new FileReader(regulatorsFile));
		while((line=br.readLine()) != null) 
			if (index.containsKey(line))
        System.out.printf("Ignoring duplicate gene %s.\n",line);
      else
        index.put(line, id++);
		br.close();

    nRootGenes=id;

    // fill Map "index" with (name,id) pairs from targets file, then store ids in array "targets"
    br=new BufferedReader(new FileReader(targetsFile));
		while((line=br.readLine()) != null)
			if (index.containsKey(line))
        System.out.printf("ignoring duplicate gene %s.\n",line);
      else
        index.put(line, id++);
		br.close();

    // now that we know the total number of genes, make the "names" array
    nGenes=id;
    geneNames=new String[nGenes];  
    for(Map.Entry<String,Integer> entry : index.entrySet())
      geneNames[entry.getValue()]=entry.getKey();

    // And now, to the Network, i.e. the list of activators and inhibitors of each target gene
		inhibitors=new int[nGenes][];  
    activators=new int[nGenes][];

    br=new BufferedReader(new FileReader(networkFile));
		br.readLine(); // discard header line

		while((line=br.readLine()) != null)
		{
			Pattern p=Pattern.compile("([^;]*);([^;]*);([^;]*)"); 
      // The pattern that lines follow is :
      // target;activators, space-separated);(inhibitors, space separated)
			Matcher m=p.matcher(line);
			if (!m.matches()) 
			{System.out.printf("In network file %s, ignoring illegal line : \n %s \n",networkFile.getName(),line);
			 continue;
			}
      String targetName=m.group(1).trim();
      if( !index.containsKey(targetName))
        System.out.println("Gene "+targetName+" unknown");
			int target=index.get(targetName);
			String[] activatorNames=m.group(2).trim().split(" +");
      if(activatorNames[0].equals("")) activatorNames=new String[0];
			String[] inhibitorNames=m.group(3).trim().split(" +");
      if(inhibitorNames[0].equals("")) inhibitorNames=new String[0];
			activators[target]=new int[activatorNames.length]; // these will be overwritten if the same...
			inhibitors[target]=new int[inhibitorNames.length]; // ...target appears twice in the newtork file
			
      for(int i=0;i<activatorNames.length; i++)
			  if (index.containsKey(activatorNames[i]))
          activators[target][i]=index.get(activatorNames[i]);
        else
				  System.out.println("Regulator "+activatorNames[i]+" not found (line :\n"+line+"\n)");
			for(int i=0;i<inhibitorNames.length; i++)
			  if (index.containsKey(inhibitorNames[i]))
			    inhibitors[target][i]=index.get(inhibitorNames[i]);
        else
				  System.out.println("Regulator "+inhibitorNames[i]+" not found (line :\n"+line+"\n)");
		}
		br.close();
  }

  /** Generate a bipartite GRN with nGenes genes, of which nRegulators Regulators, and mean activator number=mean inhibitor number=d
  */
  public GRN(int nGenes, int nRegulators, double d)
  {
    if(nRegulators < 2) throw new IllegalArgumentException("Can't make a network with less than 2 regulators : each gene needs 1 activator and 1 inhibitor at least.");
    this.nGenes=nGenes;
    this.nRootGenes=nRegulators;
		index=new HashMap<String,Integer>();
    geneNames=new String[nGenes]; 
    int regulatorNumberDigits=(nRegulators+"").length(); 
    int targetNumberDigits=(nGenes-nRegulators+"").length(); 
    for(int i=0;i<nRegulators;i++)
    {
      String name=String.format("R%0"+regulatorNumberDigits+"d",i+1);
      geneNames[i]=name;
      index.put(name,i);
    }
    for(int i=nRegulators;i<nGenes;i++)
    {
      String name=String.format("T%0"+targetNumberDigits+"d",i-nRegulators+1);
      geneNames[i]=name;
      index.put(name,i);
    }
		
    // And now, to the Network : list of activators and inhibitors of each target gene
		inhibitors=new int[nGenes][];  
    activators=new int[nGenes][];

    double[] cumulative=cumulativePoisson(d-1,20);
    for(int i=nRegulators;i<nGenes;i++)
		{
			// find him some regulators
      Random random=new Random();
      int na,ni;
      do
      {
        na=draw(cumulative,random)+1;
        ni=draw(cumulative,random)+1;
      }while(na+ni>nRegulators);
      int[] regulators=randomPart(nRegulators, na+ni, random);
      inhibitors[i]=Arrays.copyOfRange(regulators,0,ni);
      activators[i]=Arrays.copyOfRange(regulators,ni,ni+na);
		}
  }

  /** Tells whether gene with index g is regulated.
  */
  public boolean isRegulated(int g)
  {
     return (g>=nRootGenes);
  }
  
  /** First values of the cumulative Poisson distribution with mean d.
  */
  public static double[] cumulativePoisson(double d, int length)
  {
    double[] t=new double[length];
    double p=Math.exp(-d);
    double s=t[0]=p;
    for(int i=1;i<length;i++)
    {
      p*=d/i;
      s+=p;
      t[i]=s;
    }
    return t;
  }
  /** draw a number according to this cumulative distribution
  */
  public static int draw(double[] cumulative, Random random)
  {
    double x=random.nextDouble();
    int i=0;
    while(i<cumulative.length && cumulative[i]<=x)
      i++;
    return i;
  }

  /** return p distinct integers between 0 and n-1
      (more efficient than fisher-yates when p*p < n)
      Precondition : p <= n
  */
  public static int[] randomPart(int n, int p, Random random)
  {
    int[] raw=new int[p];
    int[] result=new int[p];
    for(int i=0;i<p;i++)
    {
      int x=raw[i]=random.nextInt(p-i);
      for(int j=i-1;j>=0;j--)
        if(x>=raw[j])
          x++;
      result[i]=x;
    }
    return result;
  }
}
