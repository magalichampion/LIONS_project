package coolname;

import java.io.*;
import java.util.*;
import joptsimple.*;

public class EM {
	
	final static int nStates=3;

	//inferred parameters
	public double[] mus,sigmas,alphas;
  public double epsilon;
  public boolean fixedEpsilon;
	//data 
	NetworkGraph network;
	int nGenes;
	int nIndiv=0;
	String[] indivNames={};
	double[][] expressions; // expression[indiv][gene]
  //Belief propagation
  public int BPIterations;
	double[][][] beliefs;   // beliefs[indiv][gene][state]
	double[][] deregulationBeliefs;   // deregulationBeliefs[indiv][gene]
	
  /** Initialize fields which depend only on the network, not on the dataset.
  */
  public EM(NetworkGraph network, OptionSet options)
  {
    this.network=network;
    nGenes=network.grn.nGenes;
    setOptions(options);
  }

  /** Initialize fields that are commandline options.
  */
  public void setOptions(OptionSet options)
  {
    alphas=new double[nStates];
    List temp=options.valuesOf("ialpha");
    for(int s=0;s<nStates;s++)
      alphas[s]=(Double)temp.get(s);
    normalize(this.alphas);
    mus=new double[nStates];
    temp=options.valuesOf("imu");
    for(int s=0;s<nStates;s++)
      mus[s]=(Double)temp.get(s);
    sigmas=new double[nStates];
    temp=options.valuesOf("isigma");
    for(int s=0;s<nStates;s++)
      sigmas[s]=(Double)temp.get(s);
    epsilon=(Double)options.valueOf("iepsilon");
    BPIterations=(Integer)options.valueOf("bp");
    fixedEpsilon=options.has("fixedepsilon");
  }

  /** Get ready to work on a specific dataset.
      Fill the datastuctures expressions, indivNames, nIndiv, beliefs and
      deregulationBeliefs.
  */
  public void init(double[][] expressions, String[] indivNames)
  {
    if(indivNames.length != expressions.length)
      throw new IllegalArgumentException("IndivNames and expression table have different lengths.");
    if(expressions.length>0 && expressions[0].length != nGenes)
      throw new IllegalArgumentException("Expression table has wrong number of genes.");
    this.expressions=expressions;
    this.indivNames=indivNames;
    nIndiv=indivNames.length;
    beliefs=new double[nIndiv][nGenes][];
    deregulationBeliefs=new double[nIndiv][nGenes];
  }

  /** Get ready to work on a specific dataset.
      Fill the datastuctures expressions, indivNames, nIndiv, beliefs and
      deregulationBeliefs.
  */
  public void init(double[][] expressions)
  {
    String[] indivNames=new String[expressions.length];
    init(expressions,indivNames);
  }

  /** Get ready to work on a specific dataset.
      Fill the datastuctures expressions, indivNames, nIndiv, beliefs and
      deregulationBeliefs.
  */
	public void init(File expressionFile) throws FileNotFoundException
	{
		BufferedReader br=new BufferedReader(new FileReader(expressionFile));
		try
		{
			//the first line tells which column goes with which gene.
      int[] rank=new int[nGenes];
      Arrays.fill(rank,-1);
			String line=br.readLine();
			String[] words=line.split("\\s+");
      int length=words.length;
			for(int g=0;g<length;g++)
				if(network.grn.index.containsKey(words[g]))
          rank[network.grn.index.get(words[g])]=g;
      // we must check all genes are present
      for(int g=0;g<nGenes; g++)
        if(rank[g]==-1)
          throw new IllegalArgumentException("Expression of gene "+network.grn.geneNames[g]+"not found.");
		  // the following lines are expression profiles
			nIndiv=0; // first count them ;-)
			while((line=br.readLine()) != null)
				nIndiv++;
			br.close();
			beliefs=new double[nIndiv][nGenes][];
			deregulationBeliefs=new double[nIndiv][nGenes];
			// now read and save them.
			indivNames=new String[nIndiv];
			expressions=new double[nIndiv][nGenes];
			br=new BufferedReader(new FileReader(expressionFile));
			br.readLine();
			for(int ind=0;ind<nIndiv;ind++)
			{
				words=br.readLine().split("\\s+");
				indivNames[ind]=words[0];
				if(words.length!=length+1) 
          throw new IllegalArgumentException("Line "+ind+" of expression file has length "+words.length+" instead of "+(length+1)+".");
				for(int g=0;g<nGenes;g++)
					expressions[ind][g]=Double.parseDouble(words[rank[g]+1]);
			}
		}
		catch(IOException e)
		{throw new RuntimeException("IOException while parsing expression file :\n"+e.getStackTrace());}
	}
	
  /** Holding the parameters fixed, assign the beliefs to the posterior marginals, or an
      approximation of them.
  */
	public void stepE()
	{
    // 0) set the prior deregulation proba equal to epsilon
		for(int g=0;g<nGenes;g++)
      if( network.grn.isRegulated(g) )
          network.deregulationBits[g].setInput(epsilon);

		for(int ind=0;ind<nIndiv;ind++)
		{
			//1) set inputs (computed from parameters) into Network
			for(int g=0;g<nGenes;g++)
				network.geneTrits[g].setInput(computeInput(expressions[ind][g], network.grn.isRegulated(g)));
			//2) run Belief Propagation
      network.getSolver().setNumIterations(BPIterations);
			network.solve();
			//3) retrieve beliefs (dimple returns them in normalized form)
			for (int g=0;g<nGenes;g++)
				beliefs[ind][g]=network.geneTrits[g].getBelief();
			for (int g=0;g<nGenes;g++)
        if(network.grn.isRegulated(g))
				  deregulationBeliefs[ind][g]=network.deregulationBits[g].getP1();
		}
	}
	
  /** Keeping the beliefs fixed, choose parameters that maximize the expected log-likelihood.
      See the paper for how these maximizing values are obtained.
  */
	public void stepM()
	{
                    double oldL,newL;//these indented lines are here for debugging (detect when logL decreases during step M, and which part of it)
    // find new alphas, using beliefs of non-regulated genes
                    oldL=logL1();
    for(int s=0;s<nStates;s++)
    {
      alphas[s]=0;
      for(int g=0; g<nGenes; g++)
        if( ! network.grn.isRegulated(g) )
          for(int ind=0;ind<nIndiv;ind++)
            alphas[s]+=beliefs[ind][g][s];
    }
    normalize(alphas);
                    newL=logL1();
                    if(newL<oldL) System.out.println("LogL1 decreased by "+(oldL-newL)+", ");

    if(!fixedEpsilon)
    {
                    oldL=logL3();
     // find new epsilon, using deregulation beliefs
     epsilon=0;
     for(int g=0; g<nGenes; g++)
       if( network.grn.isRegulated(g))
         for(int ind=0;ind<nIndiv;ind++)
           epsilon+=deregulationBeliefs[ind][g];
     epsilon/=nIndiv*(nGenes-network.grn.nRootGenes);
                     newL=logL3();
                     if(newL<oldL) System.out.println("LogL3 decreased by "+(oldL-newL)+", ");
    }

    // for mus and sigmas, all genes are used
                    oldL=logL2();
		for(int s=0;s<nStates;s++)
		{
			double num=0, den=0;
			for(int g=0;g<nGenes;g++)
			for(int ind=0;ind<nIndiv;ind++)
			{
			  num+=beliefs[ind][g][s]*expressions[ind][g];
			  den+=beliefs[ind][g][s];
			}
	    //if s is fully absent, mus[s] and sigmas[s] don't matter in the expected log-likelihood. Stop here to avoid a division by zero.
		  if(den==0) continue; 
 			mus[s]=num/den;
			num=0;
			for(int g=0;g<nGenes;g++)
			for(int ind=0;ind<nIndiv;ind++)
			{
				double delta=expressions[ind][g]-mus[s];
				num+=beliefs[ind][g][s]*delta*delta;
			}
			sigmas[s]=Math.sqrt(num/den);
		}
                    newL=logL2();
                    if(newL<oldL) System.out.println("LogL2 decreased by "+(oldL-newL)+", ");
	}

  //returns false if the vector's sum is 0 (in which case it is not changed)
	public static boolean normalize(double[] b)
	{
		double s=0;
		for(double d: b)
			s+=d;
		if(s==0)
		  return false;
		for(int i=0;i<b.length;i++)
			b[i]/=s;
		return true;
	}

// The following two functions are not used.
  /** Rearranges mus, sigmas and alphas so that mus are in increasing order.
  */
  public String sortParameters()
  {
    String log="";
    boolean changed=false;
    if(mus[0]>mus[1]){ swapParameters(0,1); changed=true; log+="swapped -1  with 0 ; ";}
    if(mus[1]>mus[2]){ swapParameters(1,2); changed=true; log+="swapped +1  with 0 ; ";}
    if(mus[0]>mus[1]){ swapParameters(0,1); changed=true; log+="swapped -1  with 0 ; ";}
    if( changed ) return log;
    else return null;
  }
  private void swapParameters(int i, int j)
  {
    double temp=mus[i];
    mus[i]=mus[j];
    mus[j]=temp;
    temp=sigmas[i];
    sigmas[i]=sigmas[j];
    sigmas[j]=temp;
    temp=alphas[i];
    alphas[i]=alphas[j];
    alphas[j]=temp;
  }

	/**
	 * Use Bayes' Law to compute the input factor associated with an expression value
	 * @param expr The expression value
	 * @return an array proportionnal to alphas[s]^(regulated)*exp(-z^2/2)/sigmas[s], where z=(mus[s]-expr)/sigmas[s]
	 */
	public double[] computeInput(double expr, boolean regulated)
	{
	  double[] result= new double[nStates];
	  double maxlog=Double.NEGATIVE_INFINITY;
	  for(int s=0;s<nStates;s++)
	  {
      if (sigmas[s]==0)
       result[s]= (expr==mus[s]) ? Double.MAX_VALUE : Double.MIN_VALUE;
      else
       result[s]=logNormal(mus[s],sigmas[s],expr) +  (regulated ? 0 : Math.log(alphas[s]));
       
      if(result[s]>maxlog)
        maxlog=result[s];
	  }
	  for(int s=0;s<nStates;s++)
	  {
		result[s]-=maxlog;
		result[s]=Math.exp(result[s]);
	  }
	  normalize(result);
	  return result;
	}

  public double logNormal(double mu, double sigma, double x)
  {
    double z=(x-mu)/sigma;
    return (-z*z/2 - Math.log(sigma));
  }


  // The following is useful only for debugging

  public double logL()
  {
    return logL1()+logL2()+logL3();
  }

  public double logL1() // Warning : always make sure alphas sum to one when calling this method !
  {
    double result=0;
    //Part one : hidden states of regulator genes
    for(int s=0;s<nStates;s++)
    {
      double coef=0;
      for(int ind=0;ind<nIndiv;ind++)
        for(int g=0; g<nGenes; g++)
          if(! network.grn.isRegulated(g))
            coef+=beliefs[ind][g][s];
      result+=coef*Math.log(alphas[s]);
    }
    return result;
  }

  public double logL2()
  {
    double result=0;
    //Part two : expression values
		for(int s=0;s<nStates;s++)
			for(int ind=0;ind<nIndiv;ind++)
			for(int g=0;g<nGenes;g++)
			{
			  result+=beliefs[ind][g][s]*logNormal(mus[s],sigmas[s],expressions[ind][g]);
			}
    return result;
  }

  public double logL3()
  {
    
    // part three : the deregulation bits priors, and random choice of deregulated genes' expression.
    double deregulation=0;
    int total=0;
		for(int ind=0;ind<nIndiv;ind++)
		for(int g=0;g<nGenes;g++)
      if(network.grn.isRegulated(g))
		  {
        deregulation+=deregulationBeliefs[ind][g];
        total++;
		  }
    return deregulation*(Math.log(epsilon)+Math.log(0.5))+(total-deregulation)*Math.log(1-epsilon);
  }

  /** for an experiment : see the contents of schedule */
  public void printScheduleComposition()
  {
    Map<String, Integer> classes=new HashMap<>();
    for (Object entry: network.getSchedule())
    { String name=((com.analog.lyric.dimple.schedulers.scheduleEntry.NodeScheduleEntry) entry).getNode().getClass().getName();
      Integer previous=classes.get(name);
      classes.put(name, previous==null ? 1 : (previous+1));
    }
    System.out.println(classes);
  }
}
