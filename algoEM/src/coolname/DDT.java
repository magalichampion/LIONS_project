package coolname;

import joptsimple.*;
import java.io.*;
import java.util.*;

public class DDT
{
    static  GRN grn;
    static  NetworkGraph network;
		static  EM em;
    static PrintStream log,plot;

	public static void main(String[] args)
	{
    OptionParser parser=new OptionParser();
    parser.accepts("expression").withRequiredArg(); // expression profiles
    parser.accepts("network").withRequiredArg(); // network file
    parser.accepts("iterations").withRequiredArg().ofType(Integer.class).defaultsTo(1000); // max number of EM iterations
    parser.accepts("convergence").withRequiredArg().ofType(Double.class).defaultsTo(0.00001); // alternatively, small variation in epsilon to declare convergence
    parser.accepts("bp").withRequiredArg().ofType(Integer.class).defaultsTo(10); // number of BP passes in each E step
    parser.accepts("ialpha").withRequiredArg().ofType(Double.class).withValuesSeparatedBy(':').defaultsTo(0.2,0.6,0.2); // for inference
    parser.accepts("imu").withRequiredArg().ofType(Double.class).withValuesSeparatedBy(':').defaultsTo(-1.,0.,1.); // for inference
    parser.accepts("isigma").withRequiredArg().ofType(Double.class).withValuesSeparatedBy(':').defaultsTo(0.7,0.7,0.7); // for inference
    parser.accepts("iepsilon").withRequiredArg().ofType(Double.class).defaultsTo(0.1); // for inference
    parser.accepts("plot").withRequiredArg();//output file for plot
    parser.accepts("scores").withRequiredArg();//output file for scores
    parser.accepts("log").withRequiredArg();//output file for log
    OptionSet options=parser.parse(args);

    //initialize log,plot, network
    //initialize em with expressions bp, mu, sigma, alpha, epsilon
    //run it with iterations, convergence
    //output scores.

    log=System.out;
    if(options.has("log"))
      try { log=new PrintStream(new File((String)options.valueOf("log"))); }
      catch(IOException e)
      {System.err.println("IOException opening output file \""+options.valueOf("log")+"\". \n"); e.printStackTrace();
      }

    plot=new PrintStream(new NullOutputStream());
    if(options.has("plot"))
      try { plot=new PrintStream(new File((String)options.valueOf("plot"))); }
      catch(IOException e)
      {System.err.println("IOException opening output file \""+options.valueOf("plot")+"\". \n"); e.printStackTrace();
      }

	  if(options.has("network")) // use a file-specified network
    {
		  try
      { log.print("Parsing network file...");
        grn=new GRN(new File((String)options.valueOf("network")));
        log.println("done.");
      }
      catch(IOException e)
      {
        e.printStackTrace(); 
        System.exit(0);
      }
    }
    else
     { printUsage("Missing network");
       System.exit(0);
     }
    
    network=new NetworkGraph(grn);
		em=new EM(network,options);

    if(options.has("expression"))
    {
     String path=(String) options.valueOf("expression");
     log.print("Initializing EM with expression file "+path+"...");
     try
     { em.init(new File(path)); 
       log.println("done.");
     }
     catch(FileNotFoundException e)
     {System.err.println("Expression file \""+path+" not found."); System.exit(0); }
    }
    else
    { printUsage("Missing expression file");
      System.exit(0);
    }
    
    plot.println(Arrays.toString(args)); 
    int iterations=(Integer)options.valueOf("iterations");
    double convergence=(Double)options.valueOf("convergence");

    System.out.println("Starting EM iterations");
    printStart();
		for(int i=0;i<iterations;i++)
		{
      double oldepsilon=em.epsilon;
			em.stepE();
			em.stepM();
      printAfterM();
      log.println("Iteration "+(i+1)+" done.");

      if(Math.abs(oldepsilon-em.epsilon)<convergence)
      {
        log.println("Reached cconvergence.");
        break;
      }
		} 
    if(options.has("scores"))
      outputScoresTable(new File((String) options.valueOf("scores")));
	}

  static void printAfterM()
  {
   StringBuilder line=new StringBuilder();
   for (double d: em.alphas)// 1-3 : alphas
     line.append(d+" ");
   for (double d: em.mus)// 4-6 : mus
     line.append(d+" ");
   for (double d: em.sigmas)// 7-9: sigmas
     line.append(d+" ");
   line.append(em.epsilon+" ");  //10 : epsilon
   line.append(em.logL1()+" "); //11-14 logL (reg, expr, dereg, all)
   line.append(em.logL2()+" ");
   line.append(em.logL3()+" ");
   line.append(em.logL()+" ");
   plot.println(line);
  }

  static void printStart()
  {
   StringBuilder line=new StringBuilder();
   for (double d: em.alphas)// 1-3 : alphas
     line.append(d+" ");
   for (double d: em.mus)// 4-6 : mus
     line.append(d+" ");
   for (double d: em.sigmas)// 7-9: sigmas
     line.append(d+" ");
   line.append(em.epsilon+" ");  //10 : epsilon
   plot.println(line);
  }
 
  public static void outputScoresTable(File file)
  {
    PrintWriter printer=null;
    try
    {printer=new PrintWriter(new FileWriter(file));}
    catch(IOException e)
    {e.printStackTrace(); System.exit(0);}

    StringBuilder sb=new StringBuilder();
    for(int ind=0;ind<em.nIndiv;ind++)
      sb.append(", "+em.indivNames[ind]);
    sb.deleteCharAt(0);
    printer.println(sb);

    for(int g=0;g<grn.nGenes;g++)
    {
      sb=new StringBuilder(grn.geneNames[g]);
      for(int ind=0;ind<em.nIndiv;ind++)
        sb.append(","+em.deregulationBeliefs[ind][g]);
      printer.println(sb);
    }
    printer.close();
  }

  public static void printUsage(String message)
  {
    System.err.println(message);
    System.err.println("Required parameters :");
    System.err.println("-network <path>  : the file must have one line per regulated gene, e.g. to say that gene SLIT2 has two activators SETBP1 and VGLL1 and one inhibitor GATA3 :");
    System.err.println("                SLIT2 ; SETBP1 VGLL1 ; GATA3");
    System.err.println("-expression <path> : the file must have gene names (space-separated) on the first line, then every line starts with a sample's name and contains its expression values, space-separated.");
    System.err.println("-scores <path> : where the posterior deregulation probabilities will be written");
    System.err.println("");
    System.err.println("Optional parameters :");
    System.err.println("-convergence <amount> : default is 0.00001. Iteration stops when epsilon varies by less than this amount during one EM iteration.");
    System.err.println("-iterations <n> : even if the above convergence criterion is not reached, stop after <n> iterations. Default is 1000.");
    System.err.println("-bp <n> : number of belief propagation passes done in step E. Default is 10.");
    System.err.println("-plot <path> : where to output the evolution of the parameters with each iteration.");
    System.err.println("-log <path> : where to output log messages about the program execution instead of standard output.");
  }
}
