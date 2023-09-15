import java.io.*;
import java.util.Objects;

/* args hold the command line arguments

    Example:-
    sim bimodal 6  gcc_trace.txt
    args[0] = "bimodal"
    args[1] = "6"
    ... and so on
*/
public class sim
{
    public static void main(String[] args) {

        bp_params params = new bp_params();             // check bp_params.java class for the definition of class bp_params 
        String trace_file = "";                         // Variable that holds trace file name
        char outcome;                                   // Variable holds read/write type read from input file
        long addr;                                      // Variable holds the address read from input file

        if (!(args.length == 3 || args.length == 4 || args.length == 6))    // Checks if correct number of inputs have been given. Throw error and exit if wrong
        {
            System.out.println("Error: Wrong number of inputs:" + args.length);
            System.exit(0);
        }

        params.bp_name = args[0];
        
        if(Objects.equals(params.bp_name, "bimodal"))           // Bimodal
        {
            if(args.length != 3)
            {
                System.out.printf("Error: %s wrong number of inputs:%d%n", params.bp_name, args.length);
                System.exit(0);
            }
            params.M2       = Long.parseLong(args[1]);
            trace_file      = args[2];
            System.out.printf("COMMAND%n./sim %s %d %s%n", params.bp_name, params.M2, trace_file);
        }
        else if(Objects.equals(params.bp_name, "gshare"))       // Gshare
        {
            if(args.length != 4)
            {
                System.out.printf("Error: %s wrong number of inputs:%d%n", params.bp_name, args.length);
                System.exit(0);
            }
            params.M1       = Long.parseLong(args[1]);
            params.N        = Long.parseLong(args[2]);
            trace_file      = args[3];
            System.out.printf("COMMAND%n./sim %s %d %d %s%n", params.bp_name, params.M1, params.N, trace_file);

        }
        else if(Objects.equals(params.bp_name, "hybrid"))       // Hybrid
        {
            if(args.length != 6)
            {
                System.out.printf("Error: %s wrong number of inputs:%d%n", params.bp_name, args.length);
                System.exit(0);
            }
            params.K        = Long.parseLong(args[1]);
            params.M1       = Long.parseLong(args[2]);
            params.N        = Long.parseLong(args[3]);
            params.M2       = Long.parseLong(args[4]);
            trace_file      = args[5];
            System.out.printf("COMMAND%n./sim %s %d %d %d %d %s%n", params.bp_name, params.K, params.M1, params.N, params.M2, trace_file);

        }
        else
        {
            System.out.printf("Error: Wrong branch predictor name:%s%n", params.bp_name);
            System.exit(0);
        }

        // Read file line by line
        try (BufferedReader br = new BufferedReader(new FileReader(trace_file)))
        {
            String line;
            while ((line = br.readLine()) != null) {
                String [] split = line.split("\\s+");               // split line at whitespace
                outcome = split[1].charAt(0);                       // gets t/n char from String split
                addr = Long.parseLong(split[0], 16);                // gets address from String split and converts to long. parseLong uses base 16
                if (outcome == 'n')
                    System.out.printf("%x %s%n", addr, outcome);    // Print and test if file is read correctly 
                else if (outcome == 't')
                    System.out.printf("%x %s%n", addr, outcome);    // Print and test if file is read correctly 
                /* ************************************
                  Add branch predictor code here
                **************************************/
            }
        }
        catch (IOException x)                                       // Throw error if file I/O fails
        {
            System.err.format("IOException: %s%n", x);
        }
    }
}
