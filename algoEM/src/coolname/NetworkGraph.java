package coolname;

import com.analog.lyric.dimple.model.core.FactorGraph;
import com.analog.lyric.dimple.model.variables.Bit;
import com.analog.lyric.dimple.model.factors.Factor;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.*;

public class NetworkGraph extends FactorGraph{

	final static int[][] andIndices={{0,0,0},{0,1,1},{0,2,1},{1,0,1},{1,1,1},{1,2,1},{2,0,1},{2,1,1},{2,2,2}};
	final static double[] andValues={1,1,1,1,1,1,1,1,1};
	final static int[][] licornIndices={{0,0,1},{0,1,0},{0,2,0},{1,0,2},{1,1,1},{1,2,0},{2,0,2},{2,1,2},{2,2,0}};
	final static double[] licornValues={1,1,1,1,1,1,1,1,1};
	final static int[][] muxIndices={{0,0,0},{1,1,0},{2,2,0}, {0,1,1},{0,2,1},{1,0,1},{1,2,1},{2,0,1},{2,1,1}};
	private static double[] muxValues={1,1,1,.5,.5,.5,.5,.5,.5};

	public Trit[] geneTrits;
	public Bit[] deregulationBits;
  public GRN grn;
	
	
	public NetworkGraph(GRN grn)
	{
		super();
    setOption(com.analog.lyric.dimple.options.BPOptions.updateApproach, com.analog.lyric.dimple.solvers.optimizedupdate.UpdateApproach.NORMAL);//use this in Dimple Release 0.07 to avoid using the "optimized" solve, which on our graphs is not optimized at all.
    this.grn=grn;
		geneTrits=new Trit[grn.nGenes];
		deregulationBits=new Bit[grn.nGenes];

		for(int i=0;i<grn.nGenes;i++) // here we need the topological ordering of genes in grn : if a gene comes before one of its regulators, a NullPointerException is expected soon.
    {
      geneTrits[i]=new Trit();
      if( grn.isRegulated(i) )
      {
         Trit[] activators=new Trit[grn.activators[i].length];
         Trit[] inhibitors=new Trit[grn.inhibitors[i].length];
         for(int j=0;j<activators.length;j++)
           if((activators[j]=geneTrits[grn.activators[i][j]])==null) System.out.println("gene "+i+" regulated by gene "+grn.activators[i][j]);
         for(int j=0;j<inhibitors.length;j++)
           if((inhibitors[j]=geneTrits[grn.inhibitors[i][j]])==null) System.out.println("gene "+i+" regulated by gene "+grn.inhibitors[i][j]);
         Trit licornResult=new Trit();
         addLicorn(activators, inhibitors, licornResult);
         deregulationBits[i]=new Bit();
         addMux(licornResult, geneTrits[i], deregulationBits[i]);
      }
    }
	}

  public HashMap<Integer,Integer> durations=new HashMap<>();
  public void solve()
  {
    long start=System.currentTimeMillis();
    super.solve();
    int d=(int)(System.currentTimeMillis()-start);
    Integer occ=durations.get(d);
    durations.put(d,occ==null? 1 : occ+1);
  }

	public Trit tritAnd(Trit a, Trit b)
	{
		Trit c=new Trit();
		addFactor(andIndices,andValues,a,b,c);
		return c;
	}
	
	public Trit tritsAnd(Trit[] inputs)
	{
		switch(inputs.length)
		{
			case 0:
				return new Trit();
			case 1:
				return inputs[0];
			default:
				Trit result=inputs[0];
				for(int i=1;i<inputs.length;i++)
					result=tritAnd(result,inputs[i]);
        return result;
		}
	}
	
	public Factor addLicorn(Trit [] activators, Trit[] inhibitors, Trit target)
	{
		Trit coA=tritsAnd(activators);
		Trit coI=tritsAnd(inhibitors);
		return addFactor(licornIndices,licornValues,coA,coI,target);
	}

  public Factor addMux(Trit licornResult,Trit target, Bit deregulation)
  {
     return addFactor(muxIndices, muxValues, licornResult, target, deregulation);
  }

}
